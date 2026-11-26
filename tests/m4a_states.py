"""Shared harness of the M4a tests of R-A (prereg 4.7, 5, 7): rows, kernels, judgements.

The implementation is reached only through the surface of prereg section 4:
the public names of ``xpsloss`` and the test ports ``_call_kernel``,
``integrate_window``, ``fold``, ``Panel.from_lobatto``, ``_single_loss_spectrum``
and ``_load_comparison_table``.

Conventions fixed here where the registration is silent (reported in ra_report.md):

* ``side`` of ``kernel(T, side)`` names the direction of the one-sided limit:
  ``"left"`` is T-0, ``"right"`` is T+0.  ``"interior"`` at a declared jump of an
  analytic kernel is a contract violation and raises ``KernelContractError``.
* ``features`` always carries the five keys; undeclared ones are ``()`` or ``None``.
"""
from contextlib import redirect_stderr, redirect_stdout
from fractions import Fraction as Fr
import hashlib
import io
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time
import traceback
import zipfile

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
TABLE_PATH = ROOT/"benchmarks/numerics/m4a/expected_states_m4a.json"
TABLE_SHA256 = "250b45f863ccbd1c4c94474f23769976380c5dad5c2671328712ef4acdc7f408"
RA_PATH = ROOT/"benchmarks/numerics/m4a/references_ra_m4a.json"
EPS = 2.0**-52
COLUMNS = ("term_13a", "term_13b", "term_14", "term_15_T1", "term_15_T2", "term_15_T3",
           "term_16_e", "term_16_h", "no_hole", "hole", "total")
NO_HOLE = ("term_13a", "term_13b", "term_15_T1", "term_15_T3", "term_16_e")
HOLE = ("term_14", "term_15_T2", "term_16_h")
QUANTITIES = ("zero_loss", "single_loss_total", "single_loss_no_hole", "single_loss_hole", "spectrum")
OWN_KINDS = ("state", "identity", "exact_zero_or_null", "replay", "diagnostic", "record")
AREA_SETTINGS = ("epsrel", "area_epsabs", "panel_limit", "kernel_epsabs", "kernel_epsrel", "kernel_limit")
SPECTRUM_SETTINGS = ("epsrel", "area_epsabs", "spectrum_epsabs_fraction", "panel_limit", "max_cells",
                     "kernel_epsabs", "kernel_epsrel", "kernel_limit")
WINDOW_CSV = "column,window_area_per_A,window_area_error_per_A,window_area_budget_per_A,cutoff_sensitivity,cutoff_sensitivity_error"
_CACHE = {}


class KernelContractError(AssertionError):
    """The analytic kernel was asked for an undefined value (see the module docstring)."""


def table():
    if "table" not in _CACHE:
        raw = TABLE_PATH.read_bytes()
        _CACHE["table_sha256"] = hashlib.sha256(raw).hexdigest()
        _CACHE["table"] = json.loads(raw)
    return _CACHE["table"]


def table_sha256():
    table()
    return _CACHE["table_sha256"]


def rows():
    return {r["id"]: r for r in table()["rows"]}


def full_run():
    return os.environ.get("XPSLOSS_M4A_FULL") == "1"


def selected_ids():
    """Rows of this run in the registered execution order."""
    only = os.environ.get("XPSLOSS_M4A_ROWS")
    out = []
    for rid in table()["execution_order"]:
        row = rows()[rid]
        if only:
            if rid in only.split(","):
                out.append(rid)
        elif row["ci"] or full_run():
            out.append(rid)
    return out


def plan(row):
    return [tuple(p.split("|")) for p in row["check_plan"]]


def fval(obj):
    if isinstance(obj, dict):
        if "invalid" in obj:
            return float(obj["invalid"])
        value = float.fromhex(obj["hex"])
        if value != float(obj["value"]):
            raise AssertionError(obj)
        return value
    if isinstance(obj, str):
        return float.fromhex(obj) if obj.startswith("0x") else float(obj)
    return float(obj)


def G(obj, name):
    return obj[name] if isinstance(obj, dict) else getattr(obj, name)


def settings_of(row):
    return {**table()["defaults"]["settings"], **row["settings"]}


def conditions(row):
    i, d = row["inputs"], table()["defaults"]
    return dict(window=tuple(fval(x) for x in i["loss_window_eV"]),
                L=fval(i.get("weight_length_A", d["weight_length_A"])),
                charge=fval(i.get("hole_charge", d["hole_charge"])),
                E0=fval(i.get("electron_energy_eV", d["electron_energy_eV"])),
                momentum=i.get("epsilon_momentum", d["epsilon_momentum"]))


# ---------------------------------------------------------------- the implementation (lazy)

def impl():
    import xpsloss
    from xpsloss import xps_single_loss
    return xpsloss, xps_single_loss


def material(name):
    from xpsloss import Material, Oscillator
    m = table()["materials"][name]
    return Material(name, float.fromhex(m["gap_hex"]),
                    tuple(Oscillator(*(float.fromhex(h) for h in o)) for o in m["oscillators_hex"]))


def material_dict(name):
    m = table()["materials"][name]
    return dict(name=name, gap_eV=m["gap_eV"], oscillators=[
        dict(zip(("energy_eV", "strength_eV2", "width_eV", "dispersion"), o)) for o in m["oscillators"]])


def primary_arrays(name):
    p = table()["primaries"][name]
    key = "energies_eV" if p["kind"] == "delta" else "edges_eV"
    return p["kind"], [float.fromhex(h) for h in p[key]["hex"]], [float.fromhex(h) for h in p["intensities"]["hex"]]


def make_primary(spec):
    """A primary of the table by name, or the raw (possibly invalid) one of an N row."""
    xs = impl()[1]
    if isinstance(spec, str):
        kind, x, s = primary_arrays(spec)
    elif spec["kind"] == "both":
        return (make_primary(spec["delta"]), make_primary(spec["binned"]))
    else:
        kind = spec["kind"]
        x = [fval(v) for v in spec["energies_eV" if kind == "delta" else "edges_eV"]]
        s = [fval(v) for v in spec["intensities"]]
    return (xs.DeltaPrimary if kind == "delta" else xs.BinnedPrimary)(np.array(x, float), np.array(s, float))


def make_edges(spec):
    if isinstance(spec, str):
        return np.array([float.fromhex(h) for h in table()["output_edges"][spec]["hex"]])
    return np.array([fval(v) for v in spec["values"]])


# ---------------------------------------------------------------- analytic kernels (layer rows)

_SIGN = {"left": -1, "right": 1, "interior": 0}


class AnalyticKernel:
    """kernel(T, side) -> (values[11, n], errors[11, n]); counts its calls."""
    def __init__(self, row):
        self.k = row["inputs"]["kernel"]
        self.charge = fval(row["inputs"].get("hole_charge", 1.0))
        self.calls, self.points = 0, 0

    def _inside(self, t, sign, lo, hi):
        if np.any(((t == lo) | (t == hi)) & (sign == 0)):
            raise KernelContractError(f"side 'interior' at a jump of the kernel ({lo}, {hi})")
        return ((t > lo) & (t < hi)) | ((t == lo) & (sign > 0)) | ((t == hi) & (sign < 0))

    def terms(self, t, sign):
        k, form = self.k, self.k["form"]
        if form == "rectangle":
            lo, hi = map(float, k["support_eV"])
            mask = self._inside(t, sign, lo, hi)
            terms = k["terms"] if "terms" in k else {"term_13a": k["term_13a"]}
            return {name: np.where(mask, float(v), 0.0) for name, v in terms.items()}
        if form == "piecewise_linear":
            mask = self._inside(t, sign, 2.0, 5.0)
            return {"term_13a": 0.1*np.maximum(0.0, np.minimum(t - 2.0, (5.0 - t)/2.0)),
                    "term_14": np.where(mask, 0.02*(t - 3.0), 0.0)}
        if form == "c_over_T":
            return {"term_13a": float(k["c"])/t}
        if form == "lorentzian":
            if "k_hex" in k:
                c, t0, g = float.fromhex(k["k_hex"]), 2.4, 0.5
            elif "c" in k:
                c, t0, g = float(k["c"]), fval(k["t0_eV"]), fval(k["gamma_eV"])
            else:
                c, t0, g = 1.0, 2.25, 0.3
            return {k["column"]: c/(1.0 + ((t - t0)/g)**2)}
        if form == "chebyshev":
            coefficients = np.zeros(33)
            coefficients[32] = 1.0
            return {k["column"]: fval(k["kappa"])*(1.0 - np.polynomial.chebyshev.chebval(2.0*t - 5.0, coefficients))}
        if form == "constant":
            return {k["column"]: np.full(t.shape, fval(k["value"]))}
        if form == "bump":
            c, h = fval(k["c_eV"]), fval(k["h_eV"])
            u = (t - c)/h
            return {k["column"]: np.where((u >= 0) & (u <= 1), 4.0*float(k["K0"])*u*(1.0 - u), 0.0)}
        raise AssertionError(form)

    def __call__(self, T, side):
        self.calls += 1
        t = np.asarray(T, float)
        self.points += t.size
        sides = [side]*t.size if isinstance(side, str) else list(side)
        sign = np.array([_SIGN[s] for s in sides])
        values = np.zeros((11, t.size))
        for name, v in self.terms(t, sign).items():
            values[COLUMNS.index(name)] = v
        values[8] = sum(values[COLUMNS.index(n)] for n in NO_HOLE)
        values[9] = sum(values[COLUMNS.index(n)] for n in HOLE)
        values[10] = values[8] + self.charge*values[9]
        return values, np.zeros_like(values)


def features_of(row):
    f = row["inputs"]["declared_features"]
    return dict(jumps_eV=tuple(fval(x) for x in f.get("jumps_eV", ())),
                kinks_eV=tuple(fval(x) for x in f.get("kinks_eV", ())),
                singularities=tuple((fval(s["center_eV"]), fval(s["half_width_eV"])) for s in f.get("singularities", ())),
                min_half_width_eV=fval(f["min_half_width_eV"]) if "min_half_width_eV" in f else None,
                min_variation_scale_eV=fval(f["min_variation_scale_eV"]) if "min_variation_scale_eV" in f else None)


def injection_of(row):
    spec = row["inputs"].get("injection")
    if spec is None:
        return None
    fired = {"interpolation": 0, "verification": 0}
    if row["id"] == "RA14":
        def inject(role, T, values, errors):
            fired[role] += 1
            values = np.array(values, float)
            return (np.zeros_like(values) if role == "interpolation" else values), values.copy()
    elif row["id"] in ("N26", "N27"):
        def inject(role, T, values, errors):
            fired[role] += 1
            values = np.array(values, float)
            return values, np.full(values.shape, 1e-4 if role == "interpolation" else 0.0)
    elif row["id"] == "N34":
        def inject(role, T, values, errors):
            fired[role] += 1
            values = np.array(values, float)
            return values, np.broadcast_to(np.exp(-((np.asarray(T, float) - 2.5)/0.001)**2), values.shape).copy()
    else:
        raise AssertionError(row["id"])
    inject.fired = fired
    return inject


# ---------------------------------------------------------------- outcome of a row

class Outcome:
    def __init__(self, row_id):
        self.row_id, self.status, self.exc = row_id, None, None
        self.result = self.panels = self.folded = self.kernel = self.cli = None
        self.calls, self.extras, self.seconds, self.harness_error = [], {}, 0.0, None

    def summary(self):
        out = dict(row=self.row_id, status=self.status, seconds=round(self.seconds, 4))
        if self.exc is not None:
            out.update(exception=type(self.exc).__name__, reason=getattr(self.exc, "reason", None),
                       message=str(self.exc)[:400])
        if self.harness_error:
            out["harness_error"] = self.harness_error[-600:]
        return out


def classify(exc):
    from xpsloss import NumericalConvergenceError
    return "not_converged" if isinstance(exc, NumericalConvergenceError) else "error"


def _recorder(out, original, inject=None):
    def wrapper(material_, electron_energy_eV, loss_eV, **kw):
        loss = np.array(loss_eV, float, ndmin=1)
        out.calls.append(dict(loss_eV=loss.copy(), kwargs=dict(kw), electron_energy_eV=float(electron_energy_eV)))
        if inject is not None:
            inject(loss, kw)
        return original(material_, electron_energy_eV, loss_eV, **kw)
    return wrapper


class patched:
    """Wrap _call_kernel (recording) and optionally replace _load_comparison_table."""
    def __init__(self, out, comparison_table=None, inject=None):
        self.out, self.table, self.inject = out, comparison_table, inject

    def __enter__(self):
        xs = impl()[1]
        self.xs, self.saved = xs, (xs._call_kernel, xs._load_comparison_table)
        xs._call_kernel = _recorder(self.out, xs._call_kernel, self.inject)
        if self.table is not None:
            xs._load_comparison_table = lambda: json.loads(json.dumps(self.table))
        return self

    def __exit__(self, *exc):
        self.xs._call_kernel, self.xs._load_comparison_table = self.saved
        return False


def run_layer(row, out):
    from xpsloss._loss_panels import integrate_window
    c, s, i = conditions(row), settings_of(row), row["inputs"]
    out.kernel = AnalyticKernel(row)
    inject = injection_of(row)
    out.extras["inject"] = inject
    kw = dict(features=features_of(row), weight_length_A=c["L"], hole_charge=c["charge"],
              epsrel=s["epsrel"], area_epsabs=s["area_epsabs"],
              spectrum_epsabs_fraction=s["spectrum_epsabs_fraction"], panel_limit=s["panel_limit"],
              max_cells=s["max_cells"])
    if "primary" in i:
        kw.update(primary=make_primary(i["primary"]), output_edges_eV=make_edges(i["output_edges_eV"]))
    if "initial_edges_eV" in i:
        kw["initial_edges"] = [float(x) for x in i["initial_edges_eV"]]
    if inject is not None:
        kw["inject"] = inject
    out.result, out.panels = integrate_window(out.kernel, c["window"], **kw)
    if row["id"] == "RA1":
        from xpsloss._fold import fold
        out.extras["fold3"] = fold(out.panels, kw["primary"], kw["output_edges_eV"], multiplier=3.0,
                                   hole_charge=c["charge"])


def run_layer_b(row, out):
    from xpsloss._loss_panels import Panel
    from xpsloss._fold import fold
    c, i = conditions(row), row["inputs"]
    left, right = map(float, i["declared_features"]["panels_eV"])
    values, errors = np.zeros((11, 17)), np.zeros((11, 17))
    k = i["kernel"]
    if row["id"] == "RA10":
        for name in ("term_13a", "term_14", "no_hole", "hole", "total"):
            values[COLUMNS.index(name)] = float(k[name])
    else:
        values[[0, 8, 10]] = float(k["term_13a"])
        errors[:, 2] = fval(i["injection"]["point_error"])
    panel = Panel.from_lobatto(left, right, values, errors)
    out.folded = fold([panel], make_primary(i["primary"]), make_edges(i["output_edges_eV"]),
                      multiplier=c["L"], hole_charge=c["charge"])
    G(out.folded, "outputs")


def api_call(row, out, kernel_batch_size=None, settings=None):
    xpsloss, xs = impl()
    c, i = conditions(row), row["inputs"]
    kw = dict(loss_window_eV=c["window"], weight_length_A=c["L"], epsilon_momentum=c["momentum"])
    if "hole_charge" in i:
        kw["hole_charge"] = c["charge"]
    kw.update(row["settings"] if settings is None else settings)
    mat = material(i["material"])
    if row["entry"] == "api_area":
        return xpsloss.loss_window_area(mat, c["E0"], **kw)
    primary, edges = make_primary(i["primary"]), make_edges(i["output_edges_eV"])
    if kernel_batch_size is not None:
        return xs._single_loss_spectrum(mat, c["E0"], primary, edges, kernel_batch_size=kernel_batch_size, **kw)
    return xpsloss.single_loss_spectrum(mat, c["E0"], primary, edges, **kw)


# ---------------------------------------------------------------- CLI rows

def comparison_table(name):
    """The test versions of the completed comparison table (prereg 4.6)."""
    if name is None:
        return None
    t = table()
    regions = ["outside.low"] + [f"bin[{i}]" for i in range(7)] + ["outside.high"]
    status = {"areas": {k: {"status": "pass", "reason": None} for k in t["key_groups"]["areas"]},
              "spectrum": {f"{r}.{q}": {"status": "pass", "reason": None}
                           for r in regions for q in t["key_groups"]["spectrum_per_region"]}}
    key = "outside.low.single_loss_total"
    if name == "test_table_C08_missing_key":
        del status["spectrum"][key]
    elif name == "test_table_C08_no_reference":
        status["spectrum"][key] = {"status": "no_reference", "reason": "test table: no reference for this key"}
    elif name != "test_table_C02":
        raise AssertionError(name)
    return dict(version="test", registration_sha256=TABLE_SHA256, reference_sha256=None,
                records=[dict(entry="single_loss_spectrum", conditions=t["conditions_example_C02"], key_status=status)])


def cli_inputs(row, folder):
    i = row["inputs"]
    cfg = folder/"input.json"
    cfg.write_bytes(json.dumps(dict(material=material_dict(i["config_material"]),
                                    electron_energy_eV=i["config_electron_energy_eV"],
                                    loss_grid_eV=i["config_loss_grid_eV"],
                                    integration=i["config_integration"])).encode("ascii"))
    out = folder/"output"
    argv = ["xps-loss-window-area" if row["entry"] == "cli_area" else "xps-single-loss", str(cfg)]
    spectrum = None
    if "spectrum" in i:
        s = i["spectrum"]
        kind, x, strength = primary_arrays(s["primary"])
        primary = {"kind": kind, "energies_eV" if kind == "delta" else "edges_eV": x, "intensities": strength,
                   "intensity_unit": "arb. unit", "description": s.get("description"),
                   "included_broadening": s.get("included_broadening"), "truncation": None,
                   "intrinsic_satellites": None}
        edges = table()["output_edges"][s["output_edges_eV"]]
        obj = {"schema": "xpsloss-single-loss-spectrum-v0",
               "axis": "binding_energy_eV" if "defect" in s else "kinetic_energy_eV", "primary": primary,
               "output_edges_eV": edges["grid"] if s["edges_form"] == "min_max_step" else edges["values"]}
        if "json_whitespace_and_key_order" in s:
            obj["primary"] = dict(reversed(list(primary.items())))
            text = json.dumps(dict(reversed(list(obj.items()))), indent=3, separators=(" ,", " : "))
        else:
            text = json.dumps(obj)
        spectrum = folder/"spectrum.json"
        spectrum.write_bytes(text.encode("ascii"))
        argv += ["--spectrum", str(spectrum)]
    for name, value in i["argv"].items():
        argv.append(name)
        argv += [str(v) for v in value] if isinstance(value, list) else [str(value)]
    return argv + ["--out", str(out)], cfg, spectrum, out


def run_cli(row, out):
    from xpsloss import cli
    inject = None
    if row["id"] == "N29":
        from xpsloss import DepthAverageConvergenceError
        a, b = (float(x) for x in row["inputs"]["argv"]["--loss-window-eV"])

        def inject(loss, kw):
            if "injected" in out.extras:
                return
            index = int(np.flatnonzero((loss > a) & (loss < b))[0])
            exc = DepthAverageConvergenceError("budget_not_met", loss[index], index, column="total",
                                               estimate=1.25e-3, budget=5e-4, subdivisions=3, panels=4)
            out.extras["injected"] = dict(vars(exc))
            raise exc
    with tempfile.TemporaryDirectory() as tmp:
        argv, cfg, spectrum, folder = cli_inputs(row, Path(tmp))
        stdout, stderr = io.StringIO(), io.StringIO()
        with patched(out, comparison_table(row["inputs"].get("comparison_table")), inject):
            with redirect_stdout(stdout), redirect_stderr(stderr):
                try:
                    code = cli.main(argv)
                except SystemExit as exc:
                    code = exc.code
            record = dict(exit_code=code, stdout=stdout.getvalue(), stderr=stderr.getvalue(),
                          output_dir_after="created" if folder.exists() else "absent", argv=argv,
                          config_sha256=hashlib.sha256(cfg.read_bytes()).hexdigest(),
                          spectrum_sha256=hashlib.sha256(spectrum.read_bytes()).hexdigest() if spectrum else None,
                          spectrum_name=spectrum.name if spectrum else None, files={})
            if folder.exists():
                record["files"] = {p.name: p.read_bytes() for p in folder.iterdir()}
            out.cli = record
            record["calls_in_cli"] = len(out.calls)
            if code == 0 and "run.json" in record["files"]:
                record["manifest"] = json.loads(record["files"]["run.json"])
                try:
                    if row["id"] in ("C01", "C02"):
                        out.extras["api_replay"] = cli_replay(row, record["manifest"])
                    if row["id"] == "C08":
                        out.extras["second_display"] = c08_second(row)
                except Exception:
                    out.extras["in_row_error"] = traceback.format_exc()
    try:
        record["stderr_json"] = json.loads(record["stderr"]) if record["stderr"].strip() else None
    except ValueError:
        record["stderr_json"] = "not JSON"
    if code == 0:
        out.status = "complete"
    else:
        payload = record["stderr_json"] if isinstance(record["stderr_json"], dict) else {}
        out.status = payload.get("status", "error")


def cli_replay(row, manifest):
    """The API call rebuilt from run.json alone (prereg 4.3: effective_settings)."""
    xpsloss, _ = impl()
    d = manifest["spectrum_definition"]
    mat = manifest["config"]["material"]
    from xpsloss import Material, Oscillator
    m = Material(mat["name"], mat["gap_eV"], tuple(Oscillator(**o) for o in mat["oscillators"]))
    kw = dict(loss_window_eV=tuple(d["loss_window_eV"]), weight_length_A=d["weight"]["length_A"],
              epsilon_momentum=d["epsilon_momentum"], hole_charge=d["hole_charge"], theta_deg=d["theta_deg"],
              **manifest["effective_settings"])
    if row["entry"] == "cli_area":
        return xpsloss.loss_window_area(m, d["kernel_energy_eV"], **kw)
    p = manifest["primary"]
    primary = (xpsloss.DeltaPrimary(np.array(p["energies_eV"]), np.array(p["intensities"])) if "energies_eV" in p
               else xpsloss.BinnedPrimary(np.array(p["edges_eV"]), np.array(p["intensities"])))
    return xpsloss.single_loss_spectrum(m, d["kernel_energy_eV"], primary, np.array(manifest["output_edges_eV"]), **kw)


def c08_second(row):
    """C08 once more with the key present as no_reference: the same display must come out."""
    from xpsloss import cli
    dummy = Outcome("C08b")
    with tempfile.TemporaryDirectory() as tmp:
        argv, _, _, folder = cli_inputs(row, Path(tmp))
        with patched(dummy, comparison_table("test_table_C08_no_reference")):
            with redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
                code = cli.main(argv)
        if code:
            return dict(exit_code=code)
        manifest = json.loads((folder/"run.json").read_bytes())
        return dict(exit_code=0, comparison=manifest["comparison"])


# ---------------------------------------------------------------- V rows: tag v0.2.0 in one child process

def _canonical(x):
    if isinstance(x, np.ndarray):
        return [x.dtype.str, list(x.shape), hashlib.sha256(x.tobytes()).hexdigest()]
    if isinstance(x, (float, np.floating)):
        return float(x).hex()
    if isinstance(x, (bool, np.bool_)):
        return bool(x)
    if isinstance(x, np.integer):
        return int(x)
    if isinstance(x, dict):
        return {str(k): _canonical(v) for k, v in x.items()}
    if isinstance(x, (tuple, list)):
        return [_canonical(v) for v in x]
    if isinstance(x, BaseException):
        return dict(type=type(x).__name__, message=str(x))
    return x


def v_snapshot(profile, folder):
    """Runs in the child process with PYTHONPATH on the version under test (uses no M4a name)."""
    from dataclasses import fields
    import xpsloss
    from xpsloss import cli, effective_cross_section_depth_averaged, effective_cross_section_fixed_depth
    from xpsloss.config import load_config
    t = table()
    profiles = t["v_profiles"][profile]
    ga_full = 0.2*np.arange(1, 500)
    ga = ga_full if profile == "full" else ga_full[::4]
    if hashlib.sha256(ga.astype("<f8").tobytes()).hexdigest() != profiles["GA"]["sha256_float64_le"] \
            or len(ga) != profiles["GA"]["length"]:
        raise AssertionError("GA differs from the registered profile")
    precision = dict(epsabs=1e-12, epsrel=1e-9, limit=200)
    records = dict(module_file=str(Path(xpsloss.__file__).resolve()), version=xpsloss.__version__,
                   profile=profile, points=dict(GA=len(ga)))

    def capture(function):
        try:
            r = function()
        except Exception as exc:  # the state is compared afterwards
            return dict(type=type(exc).__name__, message=str(exc), attrs=_canonical(vars(exc)),
                        reason=getattr(exc, "reason", None))
        return {f.name: _canonical(getattr(r, f.name)) for f in fields(r)}

    def fixed(name, loss):
        return capture(lambda: effective_cross_section_fixed_depth(
            material(name), 1000., loss, 10., epsilon_momentum="k_total", theta_deg=0., hole_charge=1., **precision))

    def averaged(name, loss):
        return capture(lambda: effective_cross_section_depth_averaged(
            material(name), 1000., loss, weight_length_A=20., epsilon_momentum="k_total", theta_deg=0.,
            hole_charge=1., **precision))
    records["V01"], records["V02"] = fixed("Al", ga), fixed("Si", ga)
    records["V03"], records["V04"] = averaged("Al", ga), averaged("Si", ga)
    records["V07"] = averaged("P1", np.array([1.999]))
    records["aux_i"] = fixed("P1", np.array([1.999]))
    records["aux_ii"] = averaged("Si", np.r_[ga_full[:10], 1.12])

    def main(argv):
        stdout, stderr = io.StringIO(), io.StringIO()
        with redirect_stdout(stdout), redirect_stderr(stderr):
            try:
                code = cli.main(argv)
            except SystemExit as exc:
                code = exc.code
        return code, stdout.getvalue(), stderr.getvalue()
    for rid, command, option, csv in (("V05", "xps-fixed-depth", ["--depth-A", "10"], "xps_fixed_depth.csv"),
                                      ("V06", "xps-depth-averaged", ["--weight-length-A", "20"], "xps_depth_averaged.csv")):
        cfg, out = folder/f"{rid}.json", folder/f"{rid}_out"
        cfg.write_bytes(json.dumps(dict(material=material_dict("Al"), electron_energy_eV=1000.0,
                                        loss_grid_eV=profiles["GC"]["config_loss_grid_eV"],
                                        integration=precision)).encode("ascii"))
        grid = load_config(cfg).grid()
        if hashlib.sha256(np.asarray(grid, "<f8").tobytes()).hexdigest() != profiles["GC"]["sha256_float64_le"]:
            raise AssertionError("GC differs from the registered profile")
        records["points"]["GC"] = len(grid)
        code, stdout, stderr = main([command, str(cfg), *option, "--epsilon-momentum", "k_total", "--out", str(out)])
        record = dict(exit_code=code, stderr=stderr, output_dir_after="created" if out.exists() else "absent")
        if code == 0:
            payload = json.loads(stdout)
            record["output_dir_matches"] = payload.pop("output_dir") == str(out.resolve())
            manifest = json.loads((out/"run.json").read_bytes())
            record["package_version_matches"] = manifest.pop("package_version") == xpsloss.__version__
            for key in ("created_utc", "elapsed_seconds", "input_file"):
                manifest.pop(key)
            record.update(stdout=payload, manifest=manifest, files=sorted(p.name for p in out.iterdir()),
                          csv_sha256=hashlib.sha256((out/csv).read_bytes()).hexdigest())
        records[rid] = record
    code, stdout, stderr = main(["capabilities"])
    records["V08"] = dict(exit_code=code, stdout=json.loads(stdout) if code == 0 else None, stderr=stderr)
    return records


def v_available():
    """None if the tag can be used, otherwise ('skip' | 'fail', reason) as in prereg 5.2."""
    tag = table()["tag"]
    if not (ROOT/".git").exists():
        return "skip", "the V rows need the Git checkout; unavailable in an unpacked sdist"
    if shutil.which("git"):
        got = subprocess.run(["git", "rev-parse", tag["name"] + "^{commit}"], cwd=ROOT, capture_output=True, text=True)
        if got.returncode == 0 and got.stdout.strip() == tag["commit"]:
            return None
        reason = f"tag {tag['name']} does not resolve to {tag['commit']}: {got.stdout.strip() or got.stderr.strip()}"
    else:
        reason = "git is not available"
    return "fail", reason


def v_records():
    """{'tag': records, 'current': records, 'seconds': {...}} or {'unavailable': (...)}; computed once."""
    if "v" in _CACHE:
        return _CACHE["v"]
    state = v_available()
    if state is not None:
        _CACHE["v"] = dict(unavailable=state)
        return _CACHE["v"]
    profile = "full" if full_run() else "ci"
    result = dict(profile=profile, seconds={})
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        raw = subprocess.check_output(["git", "archive", "--format=zip", table()["tag"]["commit"], "src/xpsloss"], cwd=ROOT)
        with zipfile.ZipFile(io.BytesIO(raw)) as archive:
            archive.extractall(folder/"tag")
        for label, source in (("tag", folder/"tag"/"src"), ("current", ROOT/"src")):
            work = folder/f"work_{label}"
            work.mkdir()
            dest = work/"snapshot.json"
            start = time.perf_counter()
            env = {k: v for k, v in os.environ.items() if k not in ("PYTHONUTF8", "PYTHONIOENCODING")}
            proc = subprocess.run([sys.executable, str(Path(__file__).resolve()), "--v-snapshot", profile, str(dest)],
                                  cwd=ROOT, env=env | {"PYTHONPATH": str(source)}, capture_output=True, text=True)
            result["seconds"][label] = time.perf_counter() - start
            if proc.returncode or not dest.exists():
                result[label] = dict(child_failed=(proc.stderr or proc.stdout)[-2000:])
                continue
            records = json.loads(dest.read_bytes())
            records["module_path_ok"] = Path(records["module_file"]).resolve().is_relative_to(Path(source).resolve())
            result[label] = records
    _CACHE["v"] = result
    return result


# ---------------------------------------------------------------- running the rows (cached, in order)

def run_row(rid):
    if rid in _CACHE.setdefault("outcomes", {}):
        return _CACHE["outcomes"][rid]
    row = rows()[rid]
    out = Outcome(rid)
    start = time.perf_counter()
    try:
        entry = row["entry"]
        if entry.startswith("existing"):
            v = v_records()
            out.extras["v"] = v
            if "unavailable" in v:
                out.status = "not_run"
            else:
                current = v.get("current", {})
                record = current.get(rid)
                if record is None:
                    out.status, out.harness_error = "harness_error", str(current.get("child_failed"))
                elif entry == "existing_api":
                    out.status = "success" if "type" not in record else (
                        "not_converged" if "Convergence" in record["type"] else "error")
                else:
                    out.status = "complete" if record["exit_code"] == 0 else "error"
        elif entry in ("cli_area", "cli_spectrum"):
            run_cli(row, out)
        else:
            try:
                if entry in ("layer_ab", "layer_a"):
                    run_layer(row, out)
                elif entry == "layer_b":
                    run_layer_b(row, out)
                else:
                    with patched(out):
                        if entry == "api_spectrum_pointwise":
                            c02 = run_row("C02")
                            out.result = api_call(row, out, kernel_batch_size=1,
                                                  settings=c02.cli["manifest"]["effective_settings"])
                        else:
                            out.result = api_call(row, out)
                out.status = "success"
            except (ImportError, AttributeError, TypeError, KeyError, KernelContractError, NameError):
                raise
            except Exception as exc:  # the registered states are exceptions of the implementation
                out.exc, out.status = exc, classify(exc)
                out.extras["traceback"] = traceback.format_exc()
    except Exception:
        out.status, out.harness_error = "harness_error", traceback.format_exc()
    out.seconds = time.perf_counter() - start
    _CACHE["outcomes"][rid] = out
    return out


# ---------------------------------------------------------------- reading keys

def split_region(key):
    """('outside.low' | 'bin[i]' | 'outside.high', rest)."""
    if key.startswith("outside."):
        cut = key.index(".", 8)
    else:
        cut = key.index(".")
    return key[:cut], key[cut+1:]


def read_key(result, key):
    """(value, error or None) of a key of prereg 7.1 from a result or from run.json results."""
    if key.startswith("window_area."):
        c = key[12:]
        return G(result, "window_area")[c], G(result, "window_area_error")[c]
    if key == "single_loss_weight":
        return G(result, "single_loss_weight"), G(result, "single_loss_weight_error")
    if key.startswith("inverse_window_area_interval_A."):
        interval = G(result, "inverse_window_area_interval_A")
        return (None if interval is None else interval[0 if key.endswith("lower") else 1]), None
    if key.startswith("cutoff_sensitivity."):
        c = key[19:]
        return G(result, "cutoff_sensitivity")[c], G(result, "cutoff_sensitivity_error")[c]
    region, rest = split_region(key)
    side = region[8:] if region.startswith("outside.") else None
    index = None if side else int(region[4:-1])
    if rest.startswith("cutoff_sensitivity."):
        _, s, c = rest.split(".")
        if side:
            box = G(result, "outside")[side]
            return box["cutoff_sensitivity"][s][c], box["cutoff_sensitivity_error"][s][c]
        return (G(result, "cutoff_sensitivity_bins")[s][c][index],
                G(result, "cutoff_sensitivity_bins_error")[s][c][index])
    if side:
        box = G(result, "outside")[side]
        return box[rest], box[rest + "_error"]
    return G(result, rest)[index], G(result, rest + "_error")[index]


def bits(x):
    return None if x is None else float(x).hex()


def region_names(result_or_n):
    n = result_or_n if isinstance(result_or_n, int) else len(G(result_or_n, "output_edges_eV")) - 1
    return ["outside.low"] + [f"bin[{i}]" for i in range(n)] + ["outside.high"]


def warning_codes(diagnostics):
    return [w["code"] for w in diagnostics["warnings"] if w.get("source") == "xps_single_loss"]


def diagnostics_of(out):
    """The diagnostics of the result, of the manifest, or of the exception snapshot."""
    if out.result is not None:
        return G(out.result, "diagnostics")
    if out.cli and "manifest" in out.cli:
        return out.cli["manifest"]["diagnostics"]
    if out.exc is not None:
        return getattr(out.exc, "diagnostics", None)
    return None


def zero_loss_expected(primary_name, edges_name, region):
    """G_D(0) of a delta primary (exact; used by Z1)."""
    kind, x, s = primary_arrays(primary_name)
    if kind != "delta":
        raise AssertionError(primary_name)
    edges = [float.fromhex(h) for h in table()["output_edges"][edges_name]["hex"]]
    total = 0.0
    for e, strength in zip(x, s):
        if region == "outside.low":
            hit = e < edges[0]
        elif region == "outside.high":
            hit = e > edges[-1]
        else:
            i = int(region[4:-1])
            hit = edges[i] <= e < edges[i+1] or (i == len(edges) - 2 and e == edges[-1])
        total += strength if hit else 0.0
    return total


# ---------------------------------------------------------------- judgements of the state test

class Judgement:
    """Collects (row, key, kind) -> status; a key that raises is a failure, never a pass."""
    def __init__(self):
        self.items = {}

    def run(self, rid, key, kind, check):
        try:
            verdict = check()
        except Exception as exc:  # a missing field is a failure (never a skip)
            verdict = f"{type(exc).__name__}: {str(exc)[:300]}"
        if verdict is True or verdict is None:
            self.items[(rid, key, kind)] = ("pass", None)
        elif isinstance(verdict, tuple):
            self.items[(rid, key, kind)] = verdict
        else:
            self.items[(rid, key, kind)] = ("fail", verdict)

    def failures(self, rid=None):
        return {k: v for k, v in self.items.items() if v[0] not in ("pass", "recorded") and (rid is None or k[0] == rid)}

    def counts(self):
        out = {}
        for (_, _, kind), (status, _) in self.items.items():
            out.setdefault(kind, {}).setdefault(status, 0)
            out[kind][status] += 1
        return out


def need(condition, message):
    if not condition:
        raise AssertionError(message)
    return True


def check_state(row, out):
    """The terminal state of prereg 5.2 (type, reason, attributes, zero evaluations)."""
    e, rid = row["expect"], row["id"]
    need(out.status != "harness_error", f"unregistered error (harness or implementation): {out.harness_error}")
    if e["status"] == "record_only":
        return ("recorded", out.summary())
    if row["entry"].startswith("existing"):
        return check_v(row, out)
    need(out.status == e["status"], f"status {out.status!r}, expected {e['status']!r}: {out.summary()}")
    if "exception" in e and row["entry"] not in ("cli_area", "cli_spectrum"):
        xpsloss, xs = impl()
        need(type(out.exc).__name__ == e["exception"], f"exception {type(out.exc).__name__}: {out.exc}")
        if e["exception"] == "SingleLossInputError":
            need(isinstance(out.exc, ValueError) and type(out.exc) is xs.SingleLossInputError, "not the public type")
        if e["exception"] == "SingleLossConvergenceError":
            need(isinstance(out.exc, xpsloss.NumericalConvergenceError) and type(out.exc) is xs.SingleLossConvergenceError,
                 "not the public type")
            need(isinstance(out.exc.diagnostics, dict), "no diagnostics snapshot on the exception")
        if e.get("reason") is not None:
            need(out.exc.reason == e["reason"], f"reason {out.exc.reason!r}")
        for text in e.get("message_contains", ()):
            need(text in str(out.exc), f"message lacks {text!r}: {out.exc}")
        if e.get("raised_in") == "first M3a call":
            need(len(out.calls) == 1, f"{len(out.calls)} calls of _call_kernel")
            need(upstream_equal(row, out), "type, message or public attributes differ from the direct M3a call")
    if "kernel_evaluations" in e:
        need(evaluations(row, out) == e["kernel_evaluations"], f"{evaluations(row, out)} kernel evaluations")
    if row["entry"] in ("cli_area", "cli_spectrum"):
        check_cli_state(row, out)
    if e.get("success_definition") == "layer_b_returned_on_fixed_panels":
        need(isinstance(G(out.folded, "outputs"), dict) and G(out.folded, "outputs"), "fold returned no outputs")
    for name, value in e.get("attributes", {}).items():
        check_attribute(row, out, name, value)
    if "warnings_include" in e or "warnings_exclude" in e:
        codes = warning_codes(G(out.result, "diagnostics"))
        need(len(G(out.result, "warnings")) >= len(codes), "result.warnings is shorter than the list of codes")
        for code in ("window_only", "weight_not_probability", "kernel_fixed_at_E0", *e.get("warnings_include", ())):
            need(code in codes, f"warning {code} missing: {codes}")
        for code in e.get("warnings_exclude", ()):
            need(code not in codes, f"warning {code} present")
    return True


def evaluations(row, out):
    return out.kernel.calls if row["entry"].startswith("layer") else len(out.calls)


def upstream_equal(row, out):
    """N08, N23: the exception equals the one of the direct M3a call with the recorded arguments."""
    from xpsloss import effective_cross_section_depth_averaged
    call = out.calls[0]
    rename = {"kernel_epsabs": "epsabs", "kernel_epsrel": "epsrel", "kernel_limit": "limit"}
    try:
        effective_cross_section_depth_averaged(material(row["inputs"]["material"]), call["electron_energy_eV"],
                                               call["loss_eV"], **{rename.get(k, k): v for k, v in call["kwargs"].items()})
    except Exception as direct:
        return (type(direct) is type(out.exc) and str(direct) == str(out.exc)
                and _canonical(vars(direct)) == _canonical(vars(out.exc)))
    return False


def check_attribute(row, out, name, value):
    rid = row["id"]
    d = diagnostics_of(out)
    if name == "required_panels":
        got = out.exc.required_panels if out.exc is not None else d["required_panels"]
        need(got == value and isinstance(got, (int, np.integer)), f"required_panels {got!r}")
    elif name == "resolution_basis":
        need(d["resolution_basis"] == value, f"resolution_basis {d['resolution_basis']!r}")
    elif name == "refinement_stages_min":
        need(d["stages"] >= value, f"stages {d['stages']!r}")
    elif name == "dominant_part":
        need(out.exc.dominant_part == value, f"dominant_part {out.exc.dominant_part!r}")
    elif name in ("index", "loss_eV", "node_kind", "stage", "inner_exception", "inner_reason"):
        source = out.cli["stderr_json"] if out.cli else out.exc
        need(G(source, name) == value, f"{name} {G(source, name)!r}")
    elif name == "diagnostics.error_budget.outputs":
        need(G(out.cli["stderr_json"] if out.cli else out.exc, "diagnostics")["error_budget"]["outputs"] == {}, "outputs not empty")
    elif name == "diagnostics.error_budget.max_ratio":
        need(G(out.cli["stderr_json"] if out.cli else out.exc, "diagnostics")["error_budget"]["max_ratio"] is None, "max_ratio not None")
    elif name == "index_loss_column_and_inner_attributes":
        check_n29(out, "all")
    elif name == "panel_width_ulp":
        need(out.exc.panel_width_ulp == value, f"panel_width_ulp {out.exc.panel_width_ulp!r}")
    elif name == "panel_eV":
        check_panel(row, out)
    elif name == "inverse_reason":
        need(G(out.result, "inverse_reason") == value, "inverse_reason")
    elif name == "inverse_window_area_A":
        need(G(out.result, "inverse_window_area_A") is None and G(out.result, "inverse_window_area_interval_A") is None,
             "inverse not None")
    elif name in ("kernel_derived_exact_zero", "spectrum_equals_zero_loss", "all_spectrum_outputs_exact_zero",
                  "area_keys_judged_as_usual", "budget_flag", "budget_ratio_of_zero_keys", "primary_support_fields"):
        pass  # judged key by key through the check plan of Z1 and Z2
    else:
        raise AssertionError(f"unknown registered attribute {name}")


def check_panel(row, out):
    a, b = conditions(row)["window"]
    left, right = out.exc.panel_eV
    need(a <= left < right <= b and right == math.nextafter(left, math.inf), f"panel_eV {out.exc.panel_eV!r}")
    return True


def check_n29(out, which):
    injected, got = out.extras["injected"], out.cli["stderr_json"]
    if which in ("all", "index"):
        need(got["index"] == injected["index"], "index")
    if which in ("all", "loss_eV"):
        need(got["loss_eV"] == injected["loss_eV"], "loss_eV")
    if which in ("all", "column"):
        need(got["column"] == injected["column"], "column")
    if which in ("all", "inner_attributes"):
        attrs = got["inner_attributes"]
        for name in ("reason", "loss_eV", "index", "column", "estimate", "budget"):
            need(attrs[name] == injected[name], f"inner_attributes[{name}] {attrs.get(name)!r}")
        for name, value in injected.items():
            if name in attrs and value is not None:
                need(attrs[name] == value, f"inner_attributes[{name}]")
    return True


def check_cli_state(row, out):
    e, c = row["expect"], out.cli
    need(c["exit_code"] == e["exit_code"], f"exit code {c['exit_code']}: {c['stderr'][:300]}")
    need(c["output_dir_after"] == e["output_dir_after"], f"output directory {c['output_dir_after']}")
    need((c["stdout"] + c["stderr"]).isascii(), "stdout or stderr is not ASCII")
    if e["exit_code"] == 0:
        need(sorted(c["files"]) == sorted(e.get("files", c["files"])), f"files {sorted(c['files'])}")
        m = c["manifest"]
        need(m["status"] == "complete", "run.json status")
        for path, value in e.get("run_json", {}).items():
            node = m
            for part in path.split("."):
                node = node[part]
            need(node == value and type(node) is type(value), f"run.json {path} = {node!r}")
        for group, value in e.get("comparison", {}).items():
            if group in ("areas", "spectrum"):
                need(m["comparison"][group] == value and m["diagnostics"]["comparison"][group] == value,
                     f"comparison.{group} = {m['comparison'].get(group)!r}")
    else:
        s = c["stderr_json"]
        need(isinstance(s, dict), f"stderr is not a JSON object: {c['stderr'][:300]}")
        need(not c["stdout"].strip(), "stdout not empty on failure")
        expect = dict(e.get("stderr_json", {}))
        if "exception" in e:
            expect.update(error_type=e["exception"], reason=e["reason"], status=e["status"])
        for name, value in expect.items():
            if value == "absent":
                need(name not in s, f"stderr has {name}")
            else:
                need(s.get(name) == value, f"stderr {name} = {s.get(name)!r}")
    return True


def v_unavailable_verdict(state):
    """Prereg 5.2: only an unpacked sdist (no Git checkout) skips; a checkout without the tag fails."""
    kind, reason = state
    if kind == "skip":
        return "skipped", reason
    return "fail", reason


def check_v(row, out):
    rid, e, v = row["id"], row["expect"], out.extras["v"]
    if "unavailable" in v:
        return v_unavailable_verdict(v["unavailable"])
    for label in ("tag", "current"):
        need("child_failed" not in v[label], f"{label} child failed: {v[label].get('child_failed')}")
        need(v[label]["module_path_ok"], f"{label}: wrong module path {v[label]['module_file']}")
    need(v["tag"]["version"] == "0.2.0", f"tag version {v['tag']['version']}")
    record = v["current"][rid]
    need(out.status == e["status"], f"status {out.status}")
    if "exception" in e:
        need(record["type"] == e["exception"] and record["reason"] == e["reason"], f"{record['type']}, {record['reason']}")
    if "exit_code" in e:
        need(record["exit_code"] == e["exit_code"], f"exit code {record['exit_code']}: {record['stderr'][:300]}")
    if e.get("output_dir_after"):
        need(record["output_dir_after"] == e["output_dir_after"], "output directory")
    if rid == "V08":
        need(record["stdout"]["status"] == e["stdout_json_status"], "capabilities status")
    return True


def check_v_replay(row, out):
    """All fields equal to the run of the tag (prereg 5.2)."""
    rid, v = row["id"], out.extras["v"]
    if "unavailable" in v:
        return v_unavailable_verdict(v["unavailable"])
    old, new = v["tag"][rid], v["current"][rid]
    if rid in ("V05", "V06"):
        for record in (old, new):
            need(record["output_dir_matches"], "stdout output_dir is not the real output directory")
            need(record["package_version_matches"], "run.json package_version is not __version__")
    if rid == "V08":
        o, n = dict(old["stdout"]), dict(new["stdout"])
        need(o.pop("version") == v["tag"]["version"] and n.pop("version") == v["current"]["version"], "version")
        so, sn = o.pop("supported"), n.pop("supported")
        need(sn[:len(so)] == so, f"supported prefix changed: {sn}")
        need(sn[len(so):] == row["expect"]["appended_supported"], f"appended names {sn[len(so):]}")
        for key, value in o.items():
            need(key in n and n[key] == value, f"existing block {key} changed")
        return True
    if old != new:
        if isinstance(old, dict) and isinstance(new, dict):
            different = sorted(k for k in set(old) | set(new) if old.get(k) != new.get(k))
        else:
            different = "records"
        raise AssertionError(f"differs from the tag in {different}")
    return True


# ---------------------------------------------------------------- key handlers by check kind

def tolerance_sum(*terms):
    return sum(Fr(float(t)) for t in terms)


def judge_row(row, harness, judgement):
    """Judge every key of the check plan that the state test owns."""
    rid = row["id"]
    out = harness(rid)
    for key, kind, oracle in plan(row):
        if kind not in OWN_KINDS:
            continue
        if out.status == "harness_error":
            judgement.items[(rid, key, kind)] = ("fail", "unregistered error (harness or implementation): " + str(out.harness_error).strip().splitlines()[-1][:400])
            continue
        judgement.run(rid, key, kind, lambda: dispatch(row, out, key, kind, oracle, harness))
    if rid == "RB1" and full_run():
        # in-row check of the full run: the display against the packaged completed comparison table
        try:
            shown = G(out.result, "diagnostics")["comparison"]
            out.extras["matched_display"] = ("pass" if shown["areas"] == shown["spectrum"] == "matched" else "fail", dict(shown))
        except Exception as exc:
            out.extras["matched_display"] = ("fail", f"{type(exc).__name__}: {exc}")
    return out


def dispatch(row, out, key, kind, oracle, harness):
    rid = row["id"]
    if kind == "state":
        return check_state(row, out)
    if kind == "record":
        return ("recorded", record_value(out, key))
    if kind == "replay":
        return check_replay(row, out, key, harness)
    if kind == "identity":
        return check_identity(row, out, key, harness)
    if kind == "exact_zero_or_null":
        return check_zero(row, out, key)
    return check_diagnostic(row, out, key, oracle, harness)


def record_value(out, key):
    if key.startswith("terminal"):
        return out.summary()
    if out.result is None:
        d = diagnostics_of(out)
        return None if d is None else _canonical(d.get("error_budget"))
    if key.startswith("cutoff_sensitivity"):
        return {c: [bits(G(out.result, "cutoff_sensitivity")[c]), bits(G(out.result, "cutoff_sensitivity_error")[c])]
                for c in COLUMNS}
    return _canonical(G(out.result, "diagnostics")["error_budget"])


def same_bits(a, b, key):
    (va, ea), (vb, eb) = read_key(a, key), read_key(b, key)
    need(bits(va) == bits(vb), f"value {va!r} != {vb!r}")
    need(bits(ea) == bits(eb), f"error {ea!r} != {eb!r}")
    return True


def success_result(harness, rid):
    other = harness(rid)
    need(other.status in ("success", "complete"), f"{rid} did not succeed: {other.summary()}")
    return other


def check_replay(row, out, key, harness):
    rid = row["id"]
    if row["entry"].startswith("existing"):
        return check_v_replay(row, out)
    if rid in ("C01", "C02"):
        need("api_replay" in out.extras, "no API replay: " + str(out.extras.get("in_row_error", "the CLI run did not complete"))[-400:])
        return same_bits(out.cli["manifest"]["results"], out.extras["api_replay"], key)
    need(out.result is not None, f"no result: {out.summary()}")
    if rid == "I3":
        return same_bits(out.result, success_result(harness, "RB1").result, key)
    if rid == "I8":
        return same_bits(out.result, success_result(harness, "C02").extras["api_replay"], key)
    if rid == "I4b":
        for other in ("RB1", "I4a"):
            need(bits(read_key(out.result, key)[0]) == bits(read_key(success_result(harness, other).result, key)[0]),
                 f"zero loss differs from {other}")
        return True
    raise AssertionError(f"no replay rule for {rid}")


def check_identity(row, out, key, harness):
    rid = row["id"]
    if rid == "RA1":
        name = key.split("@")[0]
        j3 = Fr(float(G(out.extras["fold3"], "outputs")[name]["value"]))
        j2 = Fr(float(read_key(out.result, name)[0]))
        need(abs(j3 - Fr(3, 2)*j2) <= 4*Fr(EPS)*abs(Fr(3, 2)*j2), f"J3 {float(j3)!r}, J2 {float(j2)!r}")
        return True
    if rid == "RA4b":
        ra = json.loads(RA_PATH.read_bytes())["auxiliary"]["RA4b.single_loss_weight@difference"]
        wa = Fr(float(G(success_result(harness, "RA4a").result, "single_loss_weight")))
        wb = Fr(float(G(out.result, "single_loss_weight")))
        need(abs(wa - wb - Fr(ra["value"])) <= Fr(ra["width"]), f"difference {float(wa - wb)!r}")
        return True
    need(out.result is not None, f"no result: {out.summary()}")
    if rid == "I1b":
        a, b = success_result(harness, "I1a").result, success_result(harness, "RB1").result
        (y, e), (ya, ea), (yb, eb) = read_key(out.result, key), read_key(a, key), read_key(b, key)
        y, ya, yb = Fr(float(y)), Fr(float(ya)), Fr(float(yb))
        bound = tolerance_sum(e, ea, 2*float(eb)) + 8*Fr(EPS)*(abs(ya) + 2*abs(yb))
        need(abs(y - (ya + 2*yb)) <= bound, f"{float(y)!r} vs {float(ya + 2*yb)!r}")
        return True
    if rid == "I4b":
        a, b = success_result(harness, "I4a").result, success_result(harness, "RB1").result
        (y, e), (ya, ea), (yb, eb) = read_key(out.result, key), read_key(a, key), read_key(b, key)
        y, ya, yb = Fr(float(y)), Fr(float(ya)), Fr(float(yb))
        need(abs(ya + y - yb) <= tolerance_sum(e, ea, eb) + 8*Fr(EPS)*abs(yb), f"{float(ya + y)!r} vs {float(yb)!r}")
        return True
    if rid == "Z2":
        other = success_result(harness, "RB1a").result
        (y, e), (yo, eo) = read_key(out.result, key), read_key(other, key)
        need(abs(Fr(float(y)) - Fr(float(yo))) <= tolerance_sum(e, eo), f"{y!r} vs {yo!r}")
        return True
    if rid == "Z1":
        region, quantity = split_region(key)
        value = read_key(out.result, key)[0]
        if quantity == "zero_loss":
            expected = zero_loss_expected(row["inputs"]["primary"], row["inputs"]["output_edges_eV"], region)
            need(bits(value) == bits(expected), f"zero loss {value!r}, expected {expected!r}")
        else:
            need(bits(value) == bits(read_key(out.result, f"{region}.zero_loss")[0]), "spectrum is not zero_loss")
        return True
    raise AssertionError(f"no identity rule for {rid}")


def check_zero(row, out, key):
    need(out.result is not None, f"no result: {out.summary()}")
    if key.startswith("inverse_window_area_interval_A."):
        need(G(out.result, "inverse_window_area_interval_A") is None and G(out.result, "inverse_window_area_A") is None,
             "inverse is not None")
        need(G(out.result, "inverse_reason") == "nonpositive_or_unresolved_area", "inverse_reason")
        return True
    value, error = read_key(out.result, key)
    need(float(value) == 0.0, f"value {value!r}")
    need(float(error) == 0.0, f"error {error!r} of an exact zero")
    return True


def candidates(d, family, end=None):
    return [c for c in d["candidates"] if c["family"] == family and (end is None or c["end"] == end)]


def near(cands, z, radius=1e-3):
    return [c for c in cands if abs(complex(c["z_re"], c["z_im"]) - z) <= radius]


_CE_ROOTS = {"CE5": (8.3811 - 0.5j, 10.0130 - 0.5j), "CE6": (8.3516 - 0.5j, 9.9875 - 0.5j)}
_CE_DISTANCE = {"CE1": ("surface", 0.0222, 0.05), "CE2": ("bulk", 0.02536, 0.01), "CE3": ("bulk", 0.00503, 0.01)}


def check_diagnostic(row, out, key, oracle, harness):
    rid, e = row["id"], row["expect"]
    d = diagnostics_of(out)
    c = out.cli
    # ---- counters
    if key == "kernel_evaluations":
        return need(evaluations(row, out) == 0, f"{evaluations(row, out)} evaluations")
    if key == "_call_kernel.calls":
        return need(len(out.calls) == 1, f"{len(out.calls)} calls")
    if key == "_call_kernel.calls >= 1":
        return need(c["calls_in_cli"] >= 1, "the normal C02 did not pass through _call_kernel")
    if key.startswith("_call_kernel."):
        name = key.split(".")[1]
        cond = conditions(row)
        expected = {"weight_length_A": cond["L"], "hole_charge": cond["charge"], "epsilon_momentum": cond["momentum"]}[name]
        need(out.calls, "no call of _call_kernel")
        # hole_charge has the upstream default 1.0; the other two are mandatory upstream
        got = {call["kwargs"].get(name, 1.0 if name == "hole_charge" else "missing") for call in out.calls}
        return need(got == {expected}, f"{name} passed upstream: {got}")
    # ---- exceptions
    if key == "exception.required_panels" or (key == "required_panels" and out.exc is not None):
        need(out.exc is not None and out.exc.required_panels == e["attributes"]["required_panels"],
             f"required_panels {getattr(out.exc, 'required_panels', None)!r}")
        return True
    if key == "exception.dominant_part":
        return need(out.exc.dominant_part == e["attributes"]["dominant_part"], f"{out.exc.dominant_part!r}")
    if key == "exception.diagnostics.required_panels":
        return need(out.exc.diagnostics["required_panels"] == 9, f"{out.exc.diagnostics['required_panels']!r}")
    if key == "exception.diagnostics.initial_cells":
        return need(out.exc.diagnostics["initial_cells"] == 89, f"{out.exc.diagnostics['initial_cells']!r}")
    if key == "exception.panel_width_ulp":
        return need(out.exc.panel_width_ulp == 1, f"{out.exc.panel_width_ulp!r}")
    if key == "exception.panel_eV":
        return check_panel(row, out)
    if key.startswith("exception.") and rid == "N28":
        name = key.split(".", 1)[1]
        return need(getattr(out.exc, name) == e["attributes"][name], f"{name} {getattr(out.exc, name)!r}")
    # ---- N29 and the failing CLI rows
    if key.startswith("stderr.") and rid == "N29":
        name = key.split(".", 1)[1]
        s = c["stderr_json"]
        if name in ("index", "loss_eV", "column", "inner_attributes"):
            return check_n29(out, name)
        expected = {"reason": e["reason"], "node_kind": e["attributes"]["node_kind"],
                    "inner_reason": e["attributes"]["inner_reason"]}[name]
        return need(s[name] == expected, f"stderr {name} {s.get(name)!r}")
    if key == "injection_fired":
        return need("injected" in out.extras, "the injection did not fire")
    if key == "exit_code":
        return need(c["exit_code"] == e["exit_code"], f"exit code {c['exit_code']}")
    if key == "output_dir_after":
        return need(c["output_dir_after"] == e["output_dir_after"], c["output_dir_after"])
    if key == "stderr.error_type":
        return need(c["stderr_json"]["error_type"] == e["stderr_json"]["error_type"], str(c["stderr_json"]))
    # ---- comparison display
    if key in ("comparison.areas", "comparison.spectrum"):
        group = key.split(".")[1]
        m = c["manifest"]
        need(m["comparison"][group] == e["comparison"][group], f"{m['comparison'][group]!r}")
        need(m["diagnostics"]["comparison"][group] == e["comparison"][group], "diagnostics.comparison differs")
        need(m["comparison"]["physical_validity"] == "not_established", "physical_validity")
        return True
    if key == "warnings[outside_matched_scope]":
        m = c["manifest"]
        expected = any(v != "matched" for v in (e["comparison"]["areas"], e["comparison"]["spectrum"]))
        return need(("outside_matched_scope" in warning_codes(m["diagnostics"])) == expected,
                    f"outside_matched_scope expected {expected}")
    if key.startswith("comparison.key_status"):
        name = "outside.low.single_loss_total"
        m = c["manifest"]
        need(m["comparison"]["key_status"]["spectrum"][name]["status"] == "no_reference", "status of the missing key")
        need(m["comparison"]["key_status"]["spectrum"][name]["reason"], "no reason for the missing key")
        second = out.extras["second_display"]
        need(second["exit_code"] == 0, "the run with test_table_C08_no_reference failed")
        for group in ("areas", "spectrum"):
            need(second["comparison"][group] == m["comparison"][group], f"display of {group} differs between the two tables")
        need(second["comparison"]["key_status"]["spectrum"][name]["status"] == "no_reference", "second table")
        return True
    # ---- run.json and CSV of C01, C02
    if key.startswith("run.json paths"):
        return check_manifest(row, out)
    if key == "run.json results":
        return check_results(row, out)
    if key.startswith("window_area.csv"):
        return check_window_csv(out, numeric="numeric" in key)
    if key.startswith("single_loss_spectrum.csv"):
        return check_spectrum_csv(out, numeric="numeric" in key)
    # ---- successful runs
    need(d is not None, f"no diagnostics: {out.summary()}")
    if key in ("diagnostics.required_panels", "required_panels"):
        return need(d["required_panels"] == e["attributes"]["required_panels"], f"{d['required_panels']!r}")
    if key == "diagnostics.stages":
        return need(d["stages"] >= 1 and d["final_panels"] > d["required_panels"], f"stages {d['stages']!r}")
    if key == "warnings":
        codes = warning_codes(d)
        expected = [oracle] if oracle.startswith("E0_") else e.get("warnings_include", [])
        for code in expected:
            need(code in codes, f"warning {code} missing: {codes}")
        if rid.startswith("CE") and "resolution_basis" in e["attributes"]:
            flag = e["attributes"]["resolution_basis"] == "computed_candidates_only"
            need(("resolution_basis_computed_candidates_only" in codes) == flag, "basis warning")
        return True
    if key == "warnings[conditional on the basis]":
        codes = warning_codes(d)
        return need((d["resolution_basis"] == "computed_candidates_only")
                    == ("resolution_basis_computed_candidates_only" in codes), "basis and warning disagree")
    if key == "diagnostics.resolution_basis":
        if rid == "Z1":
            return need(d["resolution_basis"] is None and d["resolution_basis_reason"] == "zero_response", "basis of Z1")
        return need(d["resolution_basis"] == e["attributes"]["resolution_basis"], f"{d['resolution_basis']!r}")
    if key == "diagnostics.candidates":
        return need(d["candidates"] == [], "candidates not empty")
    if key == "diagnostics.kernel_calls":
        if rid == "Z1":
            return need(d["kernel_calls"]["calls"] == 0 and not out.calls, "upstream was called")
        return check_pointwise(row, out)
    if key.startswith("diagnostics.candidates["):
        tag = key[len("diagnostics.candidates["):-1]
        if rid in _CE_DISTANCE:
            family, distance, tolerance = _CE_DISTANCE[rid]
            hits = [x for x in candidates(d, family, "lower")
                    if abs(x["distance_to_window_eV"] - distance) <= tolerance*distance]
            return need(hits, f"no {family}/lower candidate at {distance} eV from the window")
        if tag == "bulk,both":
            return need(candidates(d, "bulk", "both"), "no bulk pole with end 'both'")
        first, second = _CE_ROOTS[rid]
        h1, h2 = near(candidates(d, "surface"), first), near(candidates(d, "surface"), second)
        need(h1 and h2 and any(x is not y for x in h1 for y in h2), "two different surface candidates not found")
        if tag.startswith("surface,root2,unconfirmed"):
            need(any(x["confirmed"] is False for x in h2), "the second root is not marked unconfirmed")
        return True
    if key == "diagnostics.primary_support":
        s = d["primary_support"]
        if rid == "Z2":
            return need(all(s[k] is None for k in ("min_eV", "max_eV", "min_minus_E0_eV", "max_minus_E0_eV")), str(s))
        return need(s["min_minus_E0_eV"] == 8.0 and s["max_minus_E0_eV"] == 8.0, str(s))
    if key == "diagnostics.error_budget.outputs[*].ratio":
        outputs = d["error_budget"]["outputs"]
        need(outputs, "no outputs")
        flagged = 0
        for name, item in outputs.items():
            if item["estimate"] == 0 and item["budget"] == 0:
                need(item["ratio"] is None and item["flag"] == "zero_budget_zero_error", f"{name}: {item}")
                flagged += 1
            else:
                need(item["flag"] is None and item["ratio"] is not None, f"{name}: {item}")
        return need(flagged > 0, "no key with e = B = 0")
    raise AssertionError(f"no diagnostic rule for {rid} {key}")


def check_pointwise(row, out):
    """I8: one new interior point per call plus the sentinels a and b."""
    a, b = conditions(row)["window"]
    seen, sizes = set(), []
    for call in out.calls:
        loss = list(call["loss_eV"])
        need(a in loss and b in loss, "a call without the sentinels")
        need(len(loss) == len(set(loss)), "duplicated points in a call")
        new = [t for t in loss if t not in (a, b) and t not in seen]
        need(len(new) <= 1 and len(loss) <= 3, f"a call with {len(loss)} points, {len(new)} new")
        seen.update(new)
        sizes.append(len(loss))
    k = G(out.result, "diagnostics")["kernel_calls"]
    need(k["calls"] == len(out.calls) and list(k["batch_sizes"]) == sizes and k["points_sent_total"] == sum(sizes),
         "kernel_calls differs from the recorded calls")
    need(k["distinct_points"] == len(seen | {a, b}), "distinct_points")
    return True


def check_manifest(row, out):
    xs = impl()[1]
    c, m = out.cli, out.cli["manifest"]
    area = row["entry"] == "cli_area"
    for name in ("schema_version", "status", "command", "created_utc", "elapsed_seconds", "input_sha256", "input_file",
                 "package_version", "model_id", "runtime", "config", "units", "constants", "files", "warnings",
                 "dependent_model_ids", "spectrum_input", "spectrum_definition", "effective_settings",
                 "effective_settings_sources", "primary", "output_edges_eV", "results", "diagnostics", "comparison"):
        need(name in m, f"run.json lacks {name}")
    need(m["model_id"] == xs.SINGLE_LOSS_MODEL_ID == "simonsen1997_single_loss_spectrum_external_v0", "model_id")
    need(set(m["dependent_model_ids"]) == {"depth_averaged", "fixed_depth", "dielectric"}, "dependent_model_ids")
    need(m["input_sha256"] == c["config_sha256"], "input_sha256")
    d = m["spectrum_definition"]
    argv = row["inputs"]["argv"]
    need(d["loss_window_eV"] == [float(x) for x in argv["--loss-window-eV"]], "loss_window_eV")
    need(d["weight"]["length_A"] == float(argv["--weight-length-A"]) and d["weight"]["length_kind"] == "external_constant"
         and d["weight"]["length_source"] == argv.get("--weight-length-source"), f"weight {d['weight']}")
    need(d["epsilon_momentum"] == argv["--epsilon-momentum"] and d["hole_charge"] == 1.0 and d["theta_deg"] == 0.0
         and d["kernel_energy_eV"] == 1000.0 and d["config_loss_grid_used"] is False and "axis" in d,
         f"spectrum_definition {d}")
    names = AREA_SETTINGS if area else SPECTRUM_SETTINGS
    defaults = table()["defaults"]["settings"]
    need(set(m["effective_settings"]) == set(names) == set(m["effective_settings_sources"]), "effective_settings keys")
    need(m["diagnostics"]["effective_settings"] == m["effective_settings"], "diagnostics.effective_settings")
    for name in names:
        need(m["effective_settings"][name] == defaults[name], f"effective {name}")
        need(m["effective_settings_sources"][name] == ("config" if name.startswith("kernel_") else "default"),
             f"source of {name}: {m['effective_settings_sources'][name]}")
    if area:
        need(m["spectrum_input"] is None and m["primary"] is None and m["output_edges_eV"] is None, "area run fields")
        need("spectrum" not in m["comparison"], "an area run has a spectrum group")
    else:
        s = m["spectrum_input"]
        edges = np.array(m["output_edges_eV"], "<f8")
        need(s["sha256"] == c["spectrum_sha256"] and "file" in s, "spectrum_input")
        need(s["output_edges_sha256_float64_le"] == hashlib.sha256(edges.tobytes()).hexdigest(), "edge hash")
        need(m["primary"]["support"] == m["diagnostics"]["primary_support"], "primary.support")
    codes = warning_codes(m["diagnostics"])
    always = ["window_only", "weight_not_probability", "kernel_fixed_at_E0"] + ([] if area else ["primary_assumed_without_intrinsic_loss"])
    need(all(code in codes for code in always) and isinstance(m["warnings"], list) and m["warnings"], f"warnings {codes}")
    return True


def flatten(x):
    if isinstance(x, dict):
        return {k: flatten(v) for k, v in x.items()}
    if isinstance(x, (list, tuple, np.ndarray)):
        return [flatten(v) for v in x]
    if isinstance(x, (float, np.floating)):
        return float(x).hex()
    if isinstance(x, np.integer):
        return int(x)
    return x


def check_results(row, out):
    """run.json results: every field of prereg 4.1 equals the API replay bitwise."""
    results, api = out.cli["manifest"]["results"], out.extras["api_replay"]
    names = ["electron_energy_eV", "loss_window_eV", "weight_length_A", "epsilon_momentum", "hole_charge", "theta_deg",
             "window_area", "window_area_error", "window_area_budget", "single_loss_weight", "single_loss_weight_error",
             "inverse_window_area_A", "inverse_window_area_interval_A", "inverse_reason", "cutoff_sensitivity",
             "cutoff_sensitivity_error"]
    if row["entry"] == "cli_spectrum":
        five = ("zero_loss", "single_loss_total", "single_loss_no_hole", "single_loss_hole", "spectrum")
        names += ["output_edges_eV", "primary_integral", "primary_abs_integral", "outside", "cutoff_sensitivity_bins",
                  "cutoff_sensitivity_bins_error"]
        for q in five:
            names += [q, q + "_error", q + "_per_eV", q + "_per_eV_error"]
    for name in names:
        need(name in results, f"results lacks {name}")
        need(flatten(results[name]) == flatten(G(api, name)), f"results.{name} differs from the API replay")
    return True


def parse_csv(raw):
    need(raw.isascii(), "CSV is not ASCII")
    lines = raw.decode("ascii").splitlines()
    return lines[0], [line.split(",") for line in lines[1:]]


def check_window_csv(out, numeric):
    header, body = parse_csv(out.cli["files"]["window_area.csv"])
    api = out.extras["api_replay"]
    if not numeric:
        need(header == WINDOW_CSV, f"header {header!r}")
        return need([r[0] for r in body] == list(COLUMNS), "order of the rows")
    need(len(body) == 11, "eleven rows")
    for r in body:
        c = r[0]
        expected = [G(api, n)[c] for n in ("window_area", "window_area_error", "window_area_budget",
                                           "cutoff_sensitivity", "cutoff_sensitivity_error")]
        need([float(x).hex() for x in r[1:]] == [float(x).hex() for x in expected], f"row {c}")
    return True


def spectrum_csv_header():
    names = ["region", "edge_low_eV", "edge_high_eV"]
    for q in QUANTITIES:
        names += [q, q + "_error"]
    for s in ("plus", "minus"):
        for c in ("total", "no_hole", "hole"):
            names += [f"cutoff_{s}_{c}", f"cutoff_{s}_{c}_error"]
    return ",".join(names)


def check_spectrum_csv(out, numeric):
    header, body = parse_csv(out.cli["files"]["single_loss_spectrum.csv"])
    api = out.extras["api_replay"]
    edges = [float(x) for x in G(api, "output_edges_eV")]
    n = len(edges) - 1
    if not numeric:
        need(header == spectrum_csv_header(), f"header {header!r}")
        return need([r[0] for r in body] == ["low"] + [str(i) for i in range(n)] + ["high"], "order of the rows")
    for r, region in zip(body, region_names(n)):
        low = "" if region == "outside.low" else (edges[-1] if region == "outside.high" else edges[int(region[4:-1])])
        high = "" if region == "outside.high" else (edges[0] if region == "outside.low" else edges[int(region[4:-1]) + 1])
        for got, want in zip(r[1:3], (low, high)):
            need(got == "" if want == "" else (got != "" and float(got) == want), f"edges of {region}")
        expected = []
        for q in QUANTITIES:
            expected += list(read_key(api, f"{region}.{q}"))
        for s in ("plus", "minus"):
            for c in ("total", "no_hole", "hole"):
                expected += list(read_key(api, f"{region}.cutoff_sensitivity.{s}.{c}"))
        need([float(x).hex() for x in r[3:]] == [float(x).hex() for x in expected], f"row {region}")
    return True


def planned_counts(ids):
    out = {}
    for rid in ids:
        for _, kind, _ in plan(rows()[rid]):
            if kind in OWN_KINDS:
                out[kind] = out.get(kind, 0) + 1
    return out


def auxiliary_controls():
    """The two auxiliary controls of V07 (reported apart from the rows)."""
    v = v_records()
    if "unavailable" in v:
        return dict(status="not_run", reason=v["unavailable"][1])
    out = {}
    for name, kind, reason in (("aux_i", "FixedDepthGapDomainError", "surface_pole_not_excluded"),
                               ("aux_ii", "DielectricGapEdgeError", None)):
        new, old = v["current"].get(name), v["tag"].get(name)
        ok = bool(new) and new.get("type") == kind and (reason is None or new.get("reason") == reason) and new == old
        if name == "aux_ii" and ok:
            ok = new["attrs"].get("offending_indices") == [10]
        out[name] = "pass" if ok else "fail"
    return out


# ---------------------------------------------------------------- R-A comparison (prereg 7.2, 7.3)

def ra_references():
    if "ra" not in _CACHE:
        raw = RA_PATH.read_bytes()
        _CACHE["ra_sha256"] = hashlib.sha256(raw).hexdigest()
        _CACHE["ra"] = json.loads(raw)
    return _CACHE["ra"]


def ra_sha256():
    ra_references()
    return _CACHE["ra_sha256"]


def exact_float(x):
    x = float(x)
    if not math.isfinite(x):
        raise AssertionError(f"non-finite value {x!r}")
    return Fr(x)


def judge_direct(y, e_y, ref):
    """(status, detail) with the precision check (7.2) and the three-valued estimate check (7.3)."""
    if ref["status"] != "reference":
        return "no_reference", dict(rule=ref["no_reference_rule_applied"])
    r, u = Fr(ref["value"]), Fr(ref["u_r"])
    d_abs, d_rel = Fr(ref["d_abs"]), Fr(ref["d_rel"])
    width = d_abs + d_rel*abs(r) if r else d_abs
    y, e_y = exact_float(y), exact_float(e_y)
    d = abs(y - r)
    precision = d + u <= width
    e = e_y + 4*Fr(EPS)*max(abs(y), abs(r))
    estimate = "covered" if d + u <= e else ("underestimated" if d - u > e else "unresolved")
    detail = dict(value=float(y), reported_error=float(e_y), difference=float(d), width=float(width),
                  precision="pass" if precision else "fail", estimate=estimate)
    if not precision or estimate == "underestimated" or e_y < 0:
        return "fail", detail
    return ("pass" if estimate == "covered" else "unresolved"), detail


def round_down(x):
    f = float(x)
    return math.nextafter(f, -math.inf) if Fr(f) > x else f


def round_up(x):
    f = float(x)
    return math.nextafter(f, math.inf) if Fr(f) < x else f


def judge_inverse(result, ref_area, end):
    """Both checks of prereg 7.2 for one end of the interval of the inverse."""
    interval = G(result, "inverse_window_area_interval_A")
    if interval is None or G(result, "inverse_reason") is not None:
        return "fail", dict(reason=G(result, "inverse_reason"))
    lo, hi = (exact_float(v) for v in interval)
    r, u = Fr(ref_area["value"]), Fr(ref_area["u_r"])
    ref_lo, ref_hi = 1/(r + u), 1/(r - u)
    contained, disjoint = lo <= ref_lo and ref_hi <= hi, hi < ref_lo or ref_hi < lo
    a, e = exact_float(G(result, "window_area")["total"]), exact_float(G(result, "window_area_error")["total"])
    if a - e <= 0:
        return "fail", dict(reason="A - e_A <= 0 but an interval was returned")
    bitwise = (float(lo) == round_down(1/(a + e)) if end == "lower" else float(hi) == round_up(1/(a - e)))
    nearest = float(G(result, "inverse_window_area_A")) == float(1/a)
    detail = dict(interval=[float(lo), float(hi)], reference=[float(ref_lo), float(ref_hi)], contains_reference=contained,
                  directed_rounding=bitwise, inverse_is_nearest=nearest)
    if disjoint or not bitwise or not nearest:
        return "fail", detail
    return ("pass" if contained else "unresolved"), detail


def judge_reference_row(row, out):
    """{key: (status, detail)} for the R-A keys of a row, plus the row checks under '@...'."""
    rid = row["id"]
    refs = ra_references()["rows"][rid]["keys"]
    verdicts = {}

    def all_keys(status, detail):
        return {key: (status, detail) for key in refs}
    if out.status == "harness_error":
        return all_keys("not_run", "unregistered error (harness or implementation): " + str(out.harness_error).strip().splitlines()[-1][:400])
    if rid == "RA14":
        key, = refs
        ref = refs[key]
        try:
            y = exact_float(out.exc.diagnostics["error_budget"]["outputs"]["window_area.total"]["estimate"])
        except Exception as exc:
            return all_keys("not_run", f"no estimate in the exception snapshot ({type(exc).__name__}: {exc}); {out.summary()}")
        r, u = Fr(ref["value"]), Fr(ref["u_r"])
        ok = abs(y - r) + u <= Fr(ref["d_rel"])*abs(r)
        return {key: ("pass" if ok else "fail", dict(value=float(y), reference=float(r)))}
    if out.status != "success":
        return all_keys("not_run", f"the row did not succeed: {out.summary()}")
    if row["entry"] == "layer_b":
        outputs = G(out.folded, "outputs")
        for key, ref in refs.items():
            try:
                if ref["kind"] == "error_part":
                    name = key.rsplit(".", 1)[0]
                    y = exact_float(outputs[name]["kernel_point_error"])
                    r, u = Fr(ref["value"]), Fr(ref["u_r"])
                    verdicts[key] = ("pass" if abs(y - r) + u <= Fr(ref["d_rel"])*abs(r) else "fail",
                                     dict(value=float(y), reference=float(r)))
                else:
                    verdicts[key] = judge_direct(outputs[key]["value"], outputs[key]["estimate"], ref)
            except Exception as exc:
                verdicts[key] = ("fail", f"{type(exc).__name__}: {exc}")
        if row["expect"]["component_identity_applied"]:
            verdicts["@component_identity"] = _guard(lambda: component_identity(
                lambda name: (outputs[name]["value"], outputs[name]["estimate"]), ["bin[0]"], conditions(row)["charge"], False))
        return verdicts
    result = out.result
    for key, ref in refs.items():
        try:
            if ref["kind"] == "inverse_interval":
                verdicts[key] = judge_inverse(result, refs["window_area.total"], ref["end"])
            else:
                verdicts[key] = judge_direct(*read_key(result, key), ref)
        except Exception as exc:  # a missing field is a failure (never a skip)
            verdicts[key] = ("fail", f"{type(exc).__name__}: {str(exc)[:300]}")
    c, s = conditions(row), settings_of(row)
    spectrum = "primary" in row["inputs"]
    regions = region_names(result) if spectrum else []
    verdicts["@component_identity"] = _guard(lambda: component_identity(
        lambda name: read_key(result, name), regions, c["charge"], True))
    verdicts["@estimate_within_budget"] = _guard(lambda: within_budget(row, result, regions, s))
    if spectrum:
        verdicts["@density"] = _guard(lambda: density(result))
    return verdicts


def judge_n27(out):
    """The error_part key of N27 (oracle 'analytic'): the point-error part of the area, relative 1e-12."""
    key = "diagnostics.error_budget.outputs[window_area.total].kernel_point_error"
    ref = ra_references()["auxiliary"]["N27.window_area.total.kernel_point_error"]
    if out.status != "success":
        return {key: ("not_run", f"the row did not succeed: {out.summary()}")}
    try:
        y = exact_float(G(out.result, "diagnostics")["error_budget"]["outputs"]["window_area.total"]["kernel_point_error"])
    except Exception as exc:
        return {key: ("fail", f"{type(exc).__name__}: {exc}")}
    r = Fr(ref["value"])
    return {key: ("pass" if abs(y - r) <= Fr(ref["d_rel"])*abs(r) else "fail", dict(value=float(y), reference=float(r)))}


def _guard(check):
    try:
        check()
    except Exception as exc:
        return "fail", f"{type(exc).__name__}: {str(exc)[:400]}"
    return "pass", None


def component_identity(read, regions, charge, areas):
    """total = no_hole + Q_h*hole within the three errors and rounding (prereg 7.4)."""
    names = [("window_area.total", "window_area.no_hole", "window_area.hole")] if areas else []
    names += [tuple(f"{r}.single_loss_{c}" for c in ("total", "no_hole", "hole")) for r in regions]
    q = Fr(float(charge))
    for total, no_hole, hole in names:
        (t, et), (n, en), (h, eh) = (tuple(map(exact_float, read(k))) for k in (total, no_hole, hole))
        bound = et + en + abs(q)*eh + 8*Fr(EPS)*(abs(n) + abs(q*h) + abs(t))
        need(abs(t - (n + q*h)) <= bound, f"{total}: {float(t)!r} != {float(n)!r} + {float(q)}*{float(h)!r}")


def within_budget(row, result, regions, s):
    """e_y <= B_y for every judged quantity of a successful run (prereg 2, 7.3)."""
    rel, floor = Fr(float(s["epsrel"])), Fr(float(s["area_epsabs"]))
    for c in COLUMNS:
        y, e = map(exact_float, read_key(result, f"window_area.{c}"))
        budget = max(floor, rel*abs(y))
        need(e <= budget, f"window_area.{c}: error {float(e)!r} above the budget")
        need(abs(exact_float(G(result, "window_area_budget")[c]) - budget) <= 4*Fr(EPS)*budget, f"window_area_budget[{c}]")
    if regions:
        strength = Fr(float(G(result, "primary_abs_integral")))
        _, _, values = primary_arrays(row["inputs"]["primary"])
        need(abs(strength - sum(Fr(v) for v in values)) <= 4*Fr(EPS)*strength, "primary_abs_integral")
        eta = Fr(float(s["spectrum_epsabs_fraction"]))*strength
        for region in regions:
            for q in QUANTITIES[:4]:
                y, e = map(exact_float, read_key(result, f"{region}.{q}"))
                need(e <= max(eta, rel*abs(y)), f"{region}.{q}: error {float(e)!r} above the budget")
    outputs = G(result, "diagnostics")["error_budget"]["outputs"]
    for name, item in outputs.items():
        need(item["flag"] == "zero_budget_zero_error" or item["estimate"] <= item["budget"], f"error_budget[{name}]")
        parts = sum(Fr(float(item[p])) for p in ("interpolation", "verification", "kernel_point_error", "edge", "rounding"))
        need(abs(parts - Fr(float(item["estimate"]))) <= 8*Fr(EPS)*parts, f"error_budget[{name}]: estimate is not the sum of the parts")
    expected = {f"window_area.{c}" for c in COLUMNS} | {f"{r}.{q}" for r in regions for q in QUANTITIES[:4]}
    need(set(outputs) == expected, f"error_budget.outputs keys differ: {sorted(set(outputs) ^ expected)[:6]}")


def density(result):
    """The densities equal value/width within 2 eps (prereg 7.1)."""
    edges = [float(x) for x in G(result, "output_edges_eV")]
    for q in QUANTITIES:
        values, per = G(result, q), G(result, q + "_per_eV")
        need(len(values) == len(edges) - 1 == len(per), f"shape of {q}")
        for i, (y, p) in enumerate(zip(values, per)):
            want = Fr(float(y))/(Fr(edges[i+1]) - Fr(edges[i]))
            need(abs(Fr(float(p)) - want) <= 2*Fr(EPS)*abs(want), f"{q}_per_eV[{i}]")


if __name__ == "__main__" and sys.argv[1:2] == ["--v-snapshot"]:
    destination = Path(sys.argv[3])
    destination.write_text(json.dumps(v_snapshot(sys.argv[2], destination.parent)), encoding="ascii")
