"""CLI schema and reproducibility records for the two M4a commands."""

import csv
from dataclasses import asdict, fields
from datetime import datetime, timezone
from decimal import Decimal
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
from . import xps_multiple_loss as multiple
from . import _spectrum_schema as schema
from ._order_types import ORDER_DEFAULTS, MultipleLossInputError
from ._loss_panels import validate_settings
from ._single_loss_types import COLUMNS, COMPONENTS, QUANTITIES, SingleLossInputError
from .config import _reject_constant, _unique_object, keys, parse_config
from .dielectric import DIELECTRIC_MODEL_ID
from .model import BOHR_RADIUS_A, HBAR2_OVER_2ME_EV_A2, finite_number
from .xps_depth_averaged import DEPTH_AVERAGED_MODEL_ID
from .xps_fixed_depth import EPSILON_MOMENTA, FIXED_DEPTH_MODEL_ID

MULTIPLE = "xps-multiple-loss"
COMMANDS = ("xps-loss-window-area", "xps-single-loss", MULTIPLE)


def add_commands(sub):
    for name in COMMANDS:
        child = sub.add_parser(name)
        child.add_argument("config")
        child.add_argument("--out", type=Path, required=True)
        child.add_argument("--loss-window-eV", type=float, nargs=2, required=True)
        weight = child.add_mutually_exclusive_group(required=True)
        weight.add_argument("--weight-length-A", type=float)
        weight.add_argument("--weight-length-inverse-13a-window-area", action="store_true")
        child.add_argument("--weight-length-source")
        child.add_argument("--extra-inverse-length-per-A", type=float)
        child.add_argument("--extra-inverse-length-source")
        child.add_argument("--weight-length-seed-A", type=float)
        multiplier = child.add_mutually_exclusive_group()
        multiplier.add_argument("--multiplier-A", type=float)
        multiplier.add_argument("--multiplier-inverse-window-area", action="store_true")
        child.add_argument("--multiplier-source")
        child.add_argument("--epsilon-momentum", choices=EPSILON_MOMENTA, required=True)
        child.add_argument("--theta-deg", type=float, default=0.)
        child.add_argument("--hole-charge", type=float, default=1.)
        defaults = api.AREA_DEFAULTS
        if name != COMMANDS[0]:
            child.add_argument("--spectrum", required=True)
            defaults = api.SPECTRUM_DEFAULTS
        for opt, default in defaults.items():
            child.add_argument("--"+opt.replace("_", "-"), type=type(default))
        if name == MULTIPLE:
            child.add_argument("--model", required=True)
            child.add_argument("--max-order", required=True)
            for opt, default in ORDER_DEFAULTS.items():
                child.add_argument("--"+opt.replace("_", "-"), type=float if default is None else type(default))


def _read(path, *, measured_errors=False):
    raw = Path(path).read_bytes()
    try:
        obj = json.loads(raw.decode("utf-8-sig"), object_pairs_hook=_unique_object, parse_constant=_reject_constant)
    except ValueError as exc:
        # A malformed file cannot be dispatched from a parsed JSON object.
        # Preserve legacy errors unless it identifies the measured schema.
        if measured_errors and schema.declares_measured_schema(raw):
            raise schema.schema_error(exc) from exc
        raise
    return obj, hashlib.sha256(raw).hexdigest()


def _edge_hash(edges):
    return hashlib.sha256(edges.astype("<f8", copy=False).tobytes()).hexdigest()


def _spectrum_input(path, max_cells, higher):
    try:
        spec, digest = _read(path, measured_errors=higher)
        if higher and isinstance(spec, dict) and (spec.get("schema") == schema.SCHEMA or
                (spec.get("schema") != "xpsloss-single-loss-spectrum-v0" and "measurement" in spec)):
            primary, edges, axis, measurement = schema.parse_measured(spec, max_cells)
            record = dict(file=str(path), sha256=digest, schema=schema.SCHEMA,
                          binding_output_edges_sha256_float64_le=_edge_hash(edges))
            return primary, edges, record, axis, measurement
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
        edges = schema.edges(edge_spec, max_cells)
    except (SingleLossInputError, MultipleLossInputError):
        raise
    except (ValueError, TypeError, KeyError, AttributeError, OverflowError) as exc:
        raise SingleLossInputError("invalid_schema", detail=str(exc)) from exc
    return primary, edges, dict(file=str(path), sha256=digest, output_edges_sha256_float64_le=_edge_hash(edges)), None, None


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
    derived = args.weight_length_inverse_13a_window_area
    if derived and args.weight_length_source is not None:
        raise ValueError("--weight-length-source requires --weight-length-A")
    if not derived and any(v is not None for v in (args.extra_inverse_length_per_A, args.extra_inverse_length_source,
                                                   args.weight_length_seed_A)):
        raise ValueError("extra attenuation and seed options require --weight-length-inverse-13a-window-area")
    if args.extra_inverse_length_source is not None and args.extra_inverse_length_per_A is None:
        raise ValueError("--extra-inverse-length-source requires --extra-inverse-length-per-A")
    if args.multiplier_source is not None and args.multiplier_A is None:
        raise ValueError("--multiplier-source requires --multiplier-A")
    if args.out.exists() or args.out.is_symlink():
        raise FileExistsError(f"output path already exists: {args.out}")
    start = time.perf_counter()
    raw, digest = _read(args.config)
    config = parse_config(raw, digest)
    higher = args.command == MULTIPLE
    spectrum = args.command != COMMANDS[0]
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
    convention = "external" if args.multiplier_A is not None else (
        "inverse_window_area" if args.multiplier_inverse_window_area else "weight_length")
    kw = dict(**definition, weight_length_A=args.weight_length_A, **settings,
              weight_length_kind="inverse_13a_window_area" if derived else "external_constant",
              extra_inverse_length_per_A=args.extra_inverse_length_per_A or 0.,
              weight_length_seed_A=args.weight_length_seed_A, multiplier=convention, multiplier_A=args.multiplier_A)
    spectrum_input, primary_record, edges = None, None, None
    if spectrum:
        try:
            primary, edges, spectrum_input, axis, measurement = _spectrum_input(args.spectrum, settings["max_cells"], higher)
        except SingleLossInputError as exc:
            if higher and exc.reason == "invalid_schema":
                raise MultipleLossInputError("invalid_schema", detail=getattr(exc,"detail",str(exc))) from exc
            raise
        if higher:
            try:
                requested = int(args.max_order)
            except ValueError:
                requested = args.max_order
            for k, default in ORDER_DEFAULTS.items():
                value = getattr(args, k)
                settings[k] = default if value is None else value
                sources[k] = "default" if value is None else "cli"
            expansion = multiple.multiple_loss_spectrum(config.material, config.electron_energy_eV, primary,
                edges, **kw, model=args.model, max_order=requested, **{k:settings[k] for k in ORDER_DEFAULTS},
                **({"energy_axis": axis} if axis is not None else {}))
            result = expansion.single_loss
        else:
            result = api.single_loss_spectrum(config.material, config.electron_energy_eV, primary, edges, **kw)
        primary_record = {f.name: getattr(primary, f.name) for f in fields(primary)}
        if axis is None:
            primary_record.update(kind=api._primary_values(primary)["kind"], primary_integral=result.primary_integral,
                                  primary_abs_integral=result.primary_abs_integral, support=result.diagnostics["primary_support"])
        else:
            mask = primary.intensities > 0
            binned = isinstance(primary, multiple.BindingBinnedPrimary)
            points = (primary.binding_edges_eV[np.r_[mask, False] | np.r_[False, mask]]
                      if binned else primary.binding_energies_eV[mask])
            primary_record.update(axis=schema.AXIS, kind="binned" if binned else "delta",
                                  primary_integral=result.primary_integral, primary_abs_integral=result.primary_abs_integral,
                                  support=dict(min_eV=points.min() if points.size else None,
                                               max_eV=points.max() if points.size else None))
            spectrum_input["kinetic_output_edges_sha256_float64_le"] = _edge_hash(result.output_edges_eV)
    else:
        result = api.loss_window_area(config.material, config.electron_energy_eV, **kw)
    files = ({"multiple_loss_spectrum.csv": _csv(_order_regions(expansion)),
              "order_areas.csv": _csv(_order_areas(expansion))}
             if higher else {"window_area.csv": _csv(_areas(result))})
    if spectrum and not higher:
        files["single_loss_spectrum.csv"] = _csv(_regions(result))
    # prereg 4.4: results carries every field of 4.1, including the input copies.
    not_results = {"diagnostics", "warnings", "primary"}
    units = dict(loss_window_eV="eV", kernel_energy_eV="eV", theta_deg="degree from surface normal",
                 weight_length_A="angstrom", window_area="angstrom^-1", inverse_window_area_A="angstrom",
                 cutoff_sensitivity="dimensionless", single_loss_weight="dimensionless")
    if spectrum:
        units.update(bin_integral=primary.intensity_unit, density=primary.intensity_unit+"/eV")
    units.update(multiplier_A="angstrom", multiplier_interval_A="angstrom", weight_length_interval_A="angstrom",
                 seed_A="angstrom", area13="angstrom^-1", error13="angstrom^-1",
                 extra_inverse_length_per_A="angstrom^-1")
    weight = dict(length_A=result.weight_length_A, length_kind=result.weight_length_kind,
                  length_source=args.weight_length_source)
    if derived:
        weight.update(result.diagnostics["weight_length"], extra_inverse_length_source=args.extra_inverse_length_source)
    manifest = dict(schema_version=1, status="complete", command=args.command, created_utc=datetime.now(timezone.utc).isoformat(),
                    elapsed_seconds=time.perf_counter()-start, input_sha256=digest, input_file=args.config,
                    package_version=__version__, model_id=result.model_id,
                    dependent_model_ids=dict(depth_averaged=DEPTH_AVERAGED_MODEL_ID, fixed_depth=FIXED_DEPTH_MODEL_ID,
                                             dielectric=DIELECTRIC_MODEL_ID),
                    runtime=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__),
                    config=config.public_dict(), units=units, constants=dict(bohr_radius_A=BOHR_RADIUS_A,
                    hbar2_over_2me_eV_A2=HBAR2_OVER_2ME_EV_A2), files=list(files), warnings=result.warnings,
                    spectrum_input=spectrum_input,
                    spectrum_definition=dict(**definition, weight=weight, multiplier=dict(
                        convention=convention, value_A=result.multiplier_A, interval_A=result.multiplier_interval_A,
                        source=args.multiplier_source),
                        kernel_energy_eV=config.electron_energy_eV, axis="kinetic_energy_eV", config_loss_grid_used=False),
                    effective_settings=settings, effective_settings_sources=sources, primary=primary_record,
                    output_edges_eV=edges, results={k: v for k, v in vars(result).items() if k not in not_results},
                    diagnostics=result.diagnostics, comparison=result.diagnostics["comparison"])
    if higher:
        nested = asdict(expansion)
        if axis is None:
            del nested["energy_axis"]
        else:
            manifest["axis"] = manifest["spectrum_definition"]["axis"] = schema.AXIS
            manifest["measurement"] = measurement
            manifest["units"]["binding_energy"] = "eV"
            axis_record = nested["energy_axis"]
            axis_record["input_primary"]["kind"] = primary_record["kind"]
            manifest["energy_axis"] = axis_record
        nested["single_loss"]["comparison"] = result.diagnostics["comparison"]
        manifest["dependent_model_ids"]["single_loss"] = result.model_id
        manifest["spectrum_definition"]["max_order"] = requested
        manifest["units"].update(area="dimensionless", order_area_ratio="dimensionless")
        manifest.update({key:getattr(expansion,key) for key in ("model_id","depends_on_model_id","repeated_kernel",
            "hole_treatment","kernel_energy_policy","truncation","order_area_ratio","convention_ratio",
            "diagnostics","warnings")})
        manifest.update(results=nested, comparison=dict.fromkeys(("orders","areas","spectrum"),"not_matched"),
            physical_validity="not_established", orders_definition=[
                _order_definition(o.order,result) for o in expansion.orders])
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
    return dict(status="complete", model_id=manifest["model_id"], output_dir=str(args.out.resolve()),
                files=list(files), warning_count=len(manifest["warnings"]))


def _order_areas(result):
    yield ("order","column","area","area_error","cutoff_sensitivity","cutoff_sensitivity_error")
    for o in result.orders:
        for c in COLUMNS:
            yield (o.order,c,o.area[c],o.area_error[c],o.cutoff_sensitivity[c],o.cutoff_sensitivity_error[c])


def _order_definition(order, single):
    power = Fraction(single.multiplier_A)**order
    # Decimal(int) also converts exact coefficients beyond the interpreter's
    # decimal digit limit for str(int), without changing process-wide settings.
    return dict(order=order, multiplier_power=f"{Decimal(power.numerator)}/{Decimal(power.denominator)}",
                multiplier_A=single.multiplier_A, multiplier_interval_A=single.multiplier_interval_A)


def _order_regions(result):
    axis = result.energy_axis
    fields = [("zero_loss", result.single_loss if axis is None else result, "zero_loss")]
    fields.extend((f"order{o.order}_{c}",o,c) for o in result.orders for c in COMPONENTS)
    fields.append(("spectrum",result,"spectrum"))
    edges = result.single_loss.output_edges_eV
    if axis is None:
        columns = ["edge_low_eV", "edge_high_eV"]
        ends = {"low": (None, edges[0]), "high": (edges[-1], None)}
        first, last = "low", "high"
    else:
        columns = ["binding_edge_low_eV", "binding_edge_high_eV", "kinetic_edge_low_eV", "kinetic_edge_high_eV"]
        be = axis["input_binding_edges_eV"]
        first, last = "low_binding", "high_binding"
        ends = {first: (None, min(be[0], be[-1]), edges[-1], None),
                last: (max(be[0], be[-1]), None, None, edges[0])}
        if be[0] > be[-1]:
            first, last = last, first
    yield ["region", *columns, *(label+s for label,_,_ in fields for s in ("","_error"))]
    for i in (first, *range(len(edges)-1), last):
        if isinstance(i, int):
            if axis is None:
                limits = (edges[i], edges[i+1])
            else:
                j = len(edges)-2-i if be[0] < be[-1] else i
                limits = (*sorted((be[i], be[i+1])), edges[j], edges[j+1])
            yield [i,*limits,*(getattr(obj,k+s)[i] for _,obj,k in fields for s in ("","_error"))]
        else:
            yield [i,*ends[i],*(obj.outside[i][k+s] for _,obj,k in fields for s in ("","_error"))]

