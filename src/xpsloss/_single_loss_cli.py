"""CLI schema and reproducibility records for the two M4a commands."""

import csv
from dataclasses import fields
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
import io
import json
import math
from pathlib import Path
import platform
import time

import numpy as np
import scipy

from . import __version__
from . import xps_single_loss as api
from ._loss_panels import validate_settings
from ._single_loss_types import COLUMNS, COMPONENTS, QUANTITIES, SingleLossInputError, vector
from .config import _reject_constant, _unique_object, keys, parse_config
from .dielectric import DIELECTRIC_MODEL_ID
from .model import BOHR_RADIUS_A, HBAR2_OVER_2ME_EV_A2, finite_number
from .xps_depth_averaged import DEPTH_AVERAGED_MODEL_ID
from .xps_fixed_depth import EPSILON_MOMENTA, FIXED_DEPTH_MODEL_ID

COMMANDS = ("xps-loss-window-area", "xps-single-loss")


def add_commands(sub):
    for name in COMMANDS:
        child = sub.add_parser(name)
        child.add_argument("config")
        child.add_argument("--out", type=Path, required=True)
        child.add_argument("--loss-window-eV", type=float, nargs=2, required=True)
        child.add_argument("--weight-length-A", type=float, required=True)
        child.add_argument("--weight-length-source")
        child.add_argument("--epsilon-momentum", choices=EPSILON_MOMENTA, required=True)
        child.add_argument("--theta-deg", type=float, default=0.)
        child.add_argument("--hole-charge", type=float, default=1.)
        defaults = api.AREA_DEFAULTS
        if name == COMMANDS[1]:
            child.add_argument("--spectrum", required=True)
            defaults = api.SPECTRUM_DEFAULTS
        for opt, default in defaults.items():
            child.add_argument("--"+opt.replace("_", "-"), type=type(default))


def _read(path):
    raw = Path(path).read_bytes()
    obj = json.loads(raw.decode("utf-8-sig"), object_pairs_hook=_unique_object, parse_constant=_reject_constant)
    return obj, hashlib.sha256(raw).hexdigest()


def _edges(spec, max_cells):
    if isinstance(spec, dict):
        keys(spec, ("min", "max", "step"), label="output_edges_eV")
        lo, hi, step = (finite_number(spec[k], k) for k in ("min", "max", "step"))
        if not 0 <= lo < hi or step <= 0:
            raise ValueError("invalid output edge range")
        # The compact JSON form describes a decimal grid. Check its exact
        # interval count before rounding the expanded endpoints to binary64.
        intervals = (Fraction(str(hi))-Fraction(str(lo)))/Fraction(str(step))
        if intervals.denominator != 1:
            raise ValueError("nonintegral output edge count")
        if intervals+2 > max_cells:
            raise SingleLossInputError("cell_limit_exceeded",
                                       message_details=dict(cells_at_least=int(intervals)+2, max_cells=max_cells))
        spec = np.linspace(lo, hi, int(intervals)+1)
    return vector(spec, "invalid_output_edges", increasing=True)


def _spectrum_input(path, max_cells):
    try:
        spec, digest = _read(path)
        keys(spec, ("schema", "axis", "primary", "output_edges_eV"), label="spectrum")
        if spec["schema"] != "xpsloss-single-loss-spectrum-v0" or spec["axis"] != "kinetic_energy_eV":
            raise ValueError("unknown spectrum schema or axis")
        p = spec["primary"]
        cls = {"binned": api.BinnedPrimary, "delta": api.DeltaPrimary}.get(p.get("kind"))
        if cls is None:
            raise ValueError("invalid primary kind")
        pos = "edges_eV" if cls is api.BinnedPrimary else "energies_eV"
        keys(p, ("kind", pos, "intensities"), ("intensity_unit", "description", "included_broadening",
                                               "truncation", "intrinsic_satellites"), label="primary")
        for k in (pos, "intensities"):
            for v in p[k]:
                finite_number(v, k)
        primary = cls(**{k: v for k, v in p.items() if k != "kind"})
        edge_spec = spec["output_edges_eV"]
        if isinstance(edge_spec, list):
            for v in edge_spec:
                finite_number(v, "output_edges_eV")
        edges = _edges(edge_spec, max_cells)
    except SingleLossInputError:
        raise
    except (ValueError, TypeError, KeyError, AttributeError, OverflowError) as exc:
        raise SingleLossInputError("invalid_schema", detail=str(exc)) from exc
    return primary, edges, dict(file=str(path), sha256=digest, output_edges_sha256_float64_le=
                                hashlib.sha256(edges.astype("<f8", copy=False).tobytes()).hexdigest())


def json_value(value):
    if isinstance(value, dict):
        return {str(k): json_value(v) for k, v in value.items()}
    if isinstance(value, np.ndarray):
        return json_value(value.tolist())
    if isinstance(value, (list, tuple)):
        return [json_value(v) for v in value]
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return str(value)
    return value


def _csv(rows):
    buf = io.StringIO(newline="")
    writer = csv.writer(buf, lineterminator="\n")
    for row in rows:
        writer.writerow(format(v, ".17g") if isinstance(v, (float, np.floating)) else v for v in row)
    return buf.getvalue()


def _areas(result):
    yield ("column", "window_area_per_A", "window_area_error_per_A", "window_area_budget_per_A",
           "cutoff_sensitivity", "cutoff_sensitivity_error")
    for c in COLUMNS:
        yield (c, result.window_area[c], result.window_area_error[c], result.window_area_budget[c],
               result.cutoff_sensitivity[c], result.cutoff_sensitivity_error[c])


def _regions(result):
    names = [*QUANTITIES, *(f"cutoff_{s}_{c}" for s in ("plus", "minus") for c in COMPONENTS)]
    yield ["region", "edge_low_eV", "edge_high_eV", *(k+err for k in names for err in ("", "_error"))]
    edges = result.output_edges_eV
    for region in ["low", *range(len(edges)-1), "high"]:
        if isinstance(region, int):
            row = [region, edges[region], edges[region+1]]
            row.extend(getattr(result, k+e)[region] for k in QUANTITIES for e in ("", "_error"))
            row.extend(getattr(result, "cutoff_sensitivity_bins"+e)[s][c][region]
                       for s in ("plus", "minus") for c in COMPONENTS for e in ("", "_error"))
        else:
            outside = result.outside[region]
            row = ([region, None, edges[0]] if region == "low" else [region, edges[-1], None])
            row.extend(outside[k+e] for k in QUANTITIES for e in ("", "_error"))
            row.extend(outside["cutoff_sensitivity"+e][s][c]
                       for s in ("plus", "minus") for c in COMPONENTS for e in ("", "_error"))
        yield row


def run(args):
    if args.out.exists() or args.out.is_symlink():
        raise FileExistsError(f"output path already exists: {args.out}")
    start = time.perf_counter()
    raw, digest = _read(args.config)
    config = parse_config(raw, digest)
    spectrum = args.command == COMMANDS[1]
    settings = (api.SPECTRUM_DEFAULTS if spectrum else api.AREA_DEFAULTS).copy()
    sources = dict.fromkeys(settings, "default")
    for k in settings:
        if getattr(args, k) is not None:
            settings[k], sources[k] = getattr(args, k), "cli"
        elif k.startswith("kernel_"):
            field = k.removeprefix("kernel_")
            settings[k] = getattr(config, field)
            if field in raw.get("integration", {}):
                sources[k] = "config"
    validate_settings(settings)
    definition = dict(loss_window_eV=args.loss_window_eV, epsilon_momentum=args.epsilon_momentum,
                      theta_deg=args.theta_deg, hole_charge=args.hole_charge)
    kw = dict(**definition, weight_length_A=args.weight_length_A, **settings)
    spectrum_input, primary_record, edges = None, None, None
    if spectrum:
        primary, edges, spectrum_input = _spectrum_input(args.spectrum, settings["max_cells"])
        result = api.single_loss_spectrum(config.material, config.electron_energy_eV, primary, edges, **kw)
        primary_record = {f.name: getattr(primary, f.name) for f in fields(primary)}
        primary_record.update(kind=api._primary_values(primary)["kind"], primary_integral=result.primary_integral,
                              primary_abs_integral=result.primary_abs_integral, support=result.diagnostics["primary_support"])
    else:
        result = api.loss_window_area(config.material, config.electron_energy_eV, **kw)
    files = {"window_area.csv": _csv(_areas(result))}
    if spectrum:
        files["single_loss_spectrum.csv"] = _csv(_regions(result))
    # prereg 4.4: results carries every field of 4.1, including the input copies.
    not_results = {"diagnostics", "warnings", "primary"}
    units = dict(loss_window_eV="eV", kernel_energy_eV="eV", theta_deg="degree from surface normal",
                 weight_length_A="angstrom", window_area="angstrom^-1", inverse_window_area_A="angstrom",
                 cutoff_sensitivity="dimensionless", single_loss_weight="dimensionless")
    if spectrum:
        units.update(bin_integral=primary.intensity_unit, density=primary.intensity_unit+"/eV")
    manifest = dict(schema_version=1, status="complete", command=args.command, created_utc=datetime.now(timezone.utc).isoformat(),
                    elapsed_seconds=time.perf_counter()-start, input_sha256=digest, input_file=args.config,
                    package_version=__version__, model_id=api.SINGLE_LOSS_MODEL_ID,
                    dependent_model_ids=dict(depth_averaged=DEPTH_AVERAGED_MODEL_ID, fixed_depth=FIXED_DEPTH_MODEL_ID,
                                             dielectric=DIELECTRIC_MODEL_ID),
                    runtime=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__),
                    config=config.public_dict(), units=units, constants=dict(bohr_radius_A=BOHR_RADIUS_A,
                    hbar2_over_2me_eV_A2=HBAR2_OVER_2ME_EV_A2), files=list(files), warnings=result.warnings,
                    spectrum_input=spectrum_input,
                    spectrum_definition=dict(**definition, weight=dict(length_A=result.weight_length_A,
                        length_kind="external_constant", length_source=args.weight_length_source),
                        kernel_energy_eV=config.electron_energy_eV, axis="kinetic_energy_eV", config_loss_grid_used=False),
                    effective_settings=settings, effective_settings_sources=sources, primary=primary_record,
                    output_edges_eV=edges, results={k: v for k, v in vars(result).items() if k not in not_results},
                    diagnostics=result.diagnostics, comparison=result.diagnostics["comparison"])
    files["run.json"] = json.dumps(json_value(manifest), ensure_ascii=False, allow_nan=False, indent=2)+"\n"
    # Serialize everything before claiming the directory, then write the
    # completion marker last. Roll back only files this invocation created.
    args.out.mkdir(parents=True, exist_ok=False)
    created = []
    try:
        for filename, content in files.items():
            path = args.out/filename
            with path.open("x", encoding="utf-8", newline="") as stream:
                created.append(path)
                stream.write(content)
    except BaseException:
        for path in created:
            path.unlink()
        args.out.rmdir()
        raise
    return dict(status="complete", model_id=api.SINGLE_LOSS_MODEL_ID, output_dir=str(args.out.resolve()),
                files=list(files), warning_count=len(result.warnings))

