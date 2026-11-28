"""Shared harness of the M4b-1 tests of R-A (prereg 4, 5, 7): rows, kernels, judgements.

The implementation is reached only through the surface of prereg section 4: the
public names of ``xpsloss``, the CLI, and the test ports
``xpsloss._loss_orders.integrate_orders``, ``convolve_panels``, ``_fold_order``,
``xpsloss._loss_reference_models.integrate_orders_reference``,
``xpsloss.xps_single_loss._call_kernel`` and ``_run`` (wrapped, never replaced).
A key that cannot be read is a failure, never a skip.

Conventions fixed here where the registration is silent (reported in ra_report.md):

* ``side`` of ``kernel(T, side)``: ``"left"`` is T-0, ``"right"`` is T+0; ``"interior"`` at a
  declared jump of a box kernel raises ``KernelContractError`` (as in the M4a harness).
* ``features`` carries the five keys of the M4a analytic layer.
* ``panels_by_order`` may be a mapping keyed by the order or a sequence (index n, or n - 1
  when its length equals the number of orders).
* the CLI configuration carries ``loss_grid_eV`` {1, 100, 1} (required by the existing parser,
  not used by the command).
* this module run as ``python tests/m4b1_states.py --baseline-child OUT IDS`` is the baseline
  child of the M4A rows; it imports only names that exist in the tag v0.3.0.
"""
from contextlib import redirect_stderr, redirect_stdout
import dataclasses
from fractions import Fraction as Fr
import hashlib
import importlib
import io
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import time
import traceback
import zipfile

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DIRECTORY = ROOT/"benchmarks/numerics/m4b1"
TABLE_PATH = DIRECTORY/"expected_states_m4b1.json"
TABLE_PY = DIRECTORY/"expected_states_m4b1.py"
RA_PATH = DIRECTORY/"references_ra_m4b1.json"
TABLE_SHA256 = "d8d0f892a49c3ed9933b08d6a7c92e5b3aaf05ee5d6c7575050fc164991035e2"
TABLE_PY_SHA256 = "7c06f0bab21b5f195696fcfe6fbd4bdb0232a7bcf7d9b5f994bd3e2bfbbc248b"
BASELINE_TAG, BASELINE_COMMIT = "v0.3.0", "c40a38940d52f5e637ee987b9c6aad911adca4dc"
EPS = 2.0**-52
FEPS = Fr(1, 2**52)
COLUMNS = ("term_13a", "term_13b", "term_14", "term_15_T1", "term_15_T2", "term_15_T3",
           "term_16_e", "term_16_h", "no_hole", "hole", "total")
NO_HOLE = ("term_13a", "term_13b", "term_15_T1", "term_15_T3", "term_16_e")
HOLE = ("term_14", "term_15_T2", "term_16_h")
COMPONENTS = ("total", "no_hole", "hole")
PARTS = ("inherited", "propagated", "local", "multiplier_interval")
EXTERNAL = ("direct_reference", "bound_reference", "error_part", "interval_contains")
OWN_KINDS = ("state", "identity", "exact_zero_or_null", "bitwise", "replay", "diagnostic",
             "recomputed_from_report", "inequality")
ALWAYS_WARNINGS = ("model_is_approximation", "excess_enters_once", "depth_correlation_dropped",
                   "no_elastic_scattering", "order_areas_not_probabilities")
FIRST_ORDER_NAMES = ("epsrel", "area_epsabs", "spectrum_epsabs_fraction", "panel_limit", "max_cells",
                     "kernel_epsabs", "kernel_epsrel", "kernel_limit")
ORDER_NAMES = ("order_epsrel", "order_epsabs_fraction", "order_limit", "max_convolution_points",
               "neglect_tolerance_fraction")
UNSUPPORTED = ["oblique_emission", "trajectory_dependent_weight_length", "self_consistent_length",
               "intrinsic_loss_separation", "losses_outside_window", "elastic_scattering", "instrument_broadening",
               "background_subtraction", "energy_dependent_kernel", "binding_energy_axis"]
# Removals from UNSUPPORTED registered by later frozen pre-registrations (post-freeze test change
# approved by the owner on 2026-11-28; docs/M4B2_PREREG.md section 9.1, frozen at 62c1153).
LATER_REMOVED_UNSUPPORTED = ("binding_energy_axis",)


def unsupported_ok(actual):
    """The frozen list, with any subset of the registered removals taken out; order is kept."""
    removed = [name for name in UNSUPPORTED if name not in actual]
    return (set(removed) <= set(LATER_REMOVED_UNSUPPORTED)
            and actual == [name for name in UNSUPPORTED if name not in removed])


_CACHE = {}


class KernelContractError(AssertionError):
    """The analytic kernel was asked for an undefined value (see the module docstring)."""


class Unreadable(AssertionError):
    """A field of the surface of prereg 4 is missing or has an unexpected form."""


# Post-freeze changes of the registration (prereg 11): none.  The frozen table keeps its bytes.
POST_FREEZE_EXPECT = {}


# ---------------------------------------------------------------- table and ledger
def table():
    if "table" not in _CACHE:
        raw = TABLE_PATH.read_bytes()
        _CACHE["table_sha256"] = hashlib.sha256(raw).hexdigest()
        _CACHE["table"] = json.loads(raw.decode("utf-8"))
    return _CACHE["table"]


def table_sha256():
    table()
    return _CACHE["table_sha256"]


def rows():
    if "rows" not in _CACHE:
        _CACHE["rows"] = {r["id"]: r for r in table()["rows"]}
    return _CACHE["rows"]


def ledger():
    """Concrete (key, kind, oracle) records per row, expanded by the frozen table script itself."""
    if "ledger" not in _CACHE:
        src = TABLE_PY.read_bytes()
        if hashlib.sha256(src).hexdigest() != TABLE_PY_SHA256:
            raise AssertionError("the table script is not the frozen one")
        scope = {"__file__": str(TABLE_PY), "__name__": "m4b1_expected_states_table"}
        saved = sys.argv
        sys.argv = [str(TABLE_PY), "--check"]
        try:
            with redirect_stdout(io.StringIO()):
                exec(compile(src, str(TABLE_PY), "exec"), scope)
        except SystemExit as stop:
            if stop.code != 0:
                raise AssertionError("the table script does not reproduce the frozen table")
        finally:
            sys.argv = saved
        _CACHE["ledger"] = {rid: [tuple(r) for r in records] for rid, records in scope["LEDGER"].items()}
    return _CACHE["ledger"]


def own_plan(rid):
    return [(k, kind, oracle) for k, kind, oracle in ledger()[rid] if kind in OWN_KINDS]


def planned_counts(ids):
    counts = {}
    for rid in ids:
        for _, kind, _ in own_plan(rid):
            counts[kind] = counts.get(kind, 0) + 1
    return counts


def full_run():
    return os.environ.get("XPSLOSS_M4B1_FULL") == "1"


def selected_ids():
    """Rows of this run in the registered execution order (CI: 90 rows; XPSLOSS_M4B1_FULL=1: all 109)."""
    only = os.environ.get("XPSLOSS_M4B1_ROWS")
    out = []
    for rid in table()["execution_order"]:
        if only:
            if rid in only.split(","):
                out.append(rid)
        elif rows()[rid]["ci"] or full_run():
            out.append(rid)
    return out


def num(obj):
    if isinstance(obj, dict):
        if "invalid" in obj:
            return float(obj["invalid"])
        value = float.fromhex(obj["hex"])
        if value != float(obj["value"]):
            raise AssertionError(obj)
        return value
    return obj


def G(obj, name):
    try:
        return obj[name] if isinstance(obj, dict) else getattr(obj, name)
    except (KeyError, AttributeError, TypeError, IndexError) as exc:
        raise Unreadable(f"missing field {name!r} ({type(exc).__name__})") from None


def settings_of(row):
    d = table()["defaults"]
    out = {**d["first_order_settings"], **d["order_settings"]}
    out.update({k: num(v) for k, v in row["settings"].items()})
    return out


def fr(x):
    return Fr(float(x))


def round_down(q):
    x = float(q)
    return math.nextafter(x, -math.inf) if Fr(x) > q else x


def round_up(q):
    x = float(q)
    return math.nextafter(x, math.inf) if Fr(x) < q else x


def bits(x):
    return float(x).hex()


# ---------------------------------------------------------------- inputs
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


def make_primary(name):
    import xpsloss
    kind, x, s = primary_arrays(name)
    return (xpsloss.DeltaPrimary if kind == "delta" else xpsloss.BinnedPrimary)(np.array(x, float), np.array(s, float))


def edge_values(name):
    return [float.fromhex(h) for h in table()["output_edges"][name]["hex"]]


def region_names(n_bins):
    return ["outside.low"] + [f"bin[{i}]" for i in range(n_bins)] + ["outside.high"]


def row_regions(row):
    return region_names(len(edge_values(row["inputs"]["output_edges_eV"])) - 1)


# ---------------------------------------------------------------- analytic kernels
_SIGN = {"left": -1, "right": 1, "interior": 0}


def legendre_exact(degree):
    p0, p1 = [Fr(1)], [Fr(0), Fr(1)]
    if degree == 0:
        return p0
    for n in range(1, degree):
        nxt = [Fr(0)]*(len(p1) + 1)
        for k, c in enumerate(p1):
            nxt[k + 1] += c*Fr(2*n + 1, n + 1)
        for k, c in enumerate(p0):
            nxt[k] -= c*Fr(n, n + 1)
        p0, p1 = p1, nxt
    return p1


class AnalyticKernel:
    """kernel(T, side) -> (values[11, n], errors[11, n]) of a kernel of the table; counts its calls."""

    def __init__(self, name, charge=1.0):
        self.spec = table()["kernels"][name]
        self.charge = float(charge)
        self.calls = self.points = 0
        self.jumps = {float(x) for x in self.spec["features"]["jumps_eV"]}

    def piece(self, p, t, sign):
        kind = p["type"]
        if kind == "box":
            lo, hi = (num(x) for x in p["interval_eV"])
            at_jump = ((t == lo) | (t == hi)) & (sign == 0)
            if np.any(at_jump):
                raise KernelContractError(f"side 'interior' at a jump of the kernel ({lo}, {hi})")
            inside = ((t > lo) & (t < hi)) | ((t == lo) & (sign > 0)) | ((t == hi) & (sign < 0))
            return np.where(inside, num(p["value"]), 0.0)
        if kind == "exponential":
            return num(p["amplitude"])*np.exp(-num(p["decay_per_eV"])*t)
        if kind == "lorentzian":
            return num(p["amplitude"])/(1.0 + ((t - num(p["center_eV"]))/num(p["half_width_eV"]))**2)
        if kind == "linear":
            return num(p["slope_per_eV"])*(t - num(p["zero_eV"]))
        if kind == "legendre":
            return np.array([legendre_value(p, x) for x in t])
        raise AssertionError(kind)

    def __call__(self, T, side):
        self.calls += 1
        t = np.array(T, float, ndmin=1)
        self.points += t.size
        sides = [side]*t.size if isinstance(side, str) else [str(s) for s in np.array(side, ndmin=1)]
        sign = np.array([_SIGN[s] for s in sides])
        values = np.zeros((11, t.size))
        for name, pieces in self.spec["columns"].items():
            for p in pieces:
                values[COLUMNS.index(name)] += self.piece(p, t, sign)
        values[8] = sum(values[COLUMNS.index(n)] for n in NO_HOLE)
        values[9] = sum(values[COLUMNS.index(n)] for n in HOLE)
        values[10] = values[8] + self.charge*values[9]
        return values, np.zeros_like(values)


def legendre_value(piece, x):
    """scale*(1 + P_d((2 T - 7)/3)) evaluated in exact rationals at the binary64 T, rounded once."""
    key = ("legendre", piece["degree"])
    if key not in _CACHE:
        _CACHE[key] = legendre_exact(piece["degree"])
    u = (2*Fr(float(x)) - 7)/3
    total = Fr(0)
    for c in reversed(_CACHE[key]):
        total = total*u + c
    return float(Fr(num(piece["scale"]))*(1 + total))


def features_of(spec):
    f = spec["features"]
    return dict(jumps_eV=tuple(float(x) for x in f.get("jumps_eV", ())),
                kinks_eV=tuple(float(x) for x in f.get("kinks_eV", ())),
                singularities=tuple((float(c), float(g)) for c, g in f.get("singularities", ())),
                min_half_width_eV=float(f["min_half_width_eV"]) if "min_half_width_eV" in f else None,
                min_variation_scale_eV=float(f["min_variation_scale_eV"]) if "min_variation_scale_eV" in f else None)


def make_inject(spec):
    """inject(role, T, values, errors) of prereg 4.6 item 2: point errors on the named columns and roles."""
    roles = set(spec["roles"])
    errors_by_column = {}
    for name, value in spec["errors"].items():
        for col in (COLUMNS if name == "all 11 columns" else (name,)):
            errors_by_column[COLUMNS.index(col)] = num(value)
    fired = {"interpolation": 0, "verification": 0}

    def inject(role, T, values, errors):
        fired[role] += 1
        values, errors = np.array(values, float), np.array(errors, float)
        if role in roles:
            for j, e in errors_by_column.items():
                errors[j] = e
        return values, errors
    inject.fired = fired
    return inject


def make_order_inject(row):
    """order_inject of prereg 4.6 item 3 for OI1 and OI2 (order 2, interpolation nodes, three columns)."""
    spec = row["inputs"]["order_inject"]
    cols = [COLUMNS.index(c) for c in spec["columns"]]
    seen = []

    def order_inject(order, role, panel_eV, nodes_eV, values, local_errors):
        values, local_errors = np.array(values, float), np.array(local_errors, float)
        if order == spec["order"] and role in spec["roles"]:
            width = float(panel_eV[1]) - float(panel_eV[0])
            e = 1e-12 if row["id"] == "OI2" or width <= 1.5 else 1e-3
            local_errors[cols] = e
            seen.append((float(panel_eV[0]), float(panel_eV[1]), e))
        return values, local_errors
    order_inject.seen = seen
    return order_inject


# ---------------------------------------------------------------- outcome of a row
class Outcome:
    def __init__(self, row_id):
        self.row_id, self.status, self.exc = row_id, None, None
        self.result = self.panels = self.kernel = self.cli = None
        self.calls, self.deliveries, self.run_calls = [], [], 0
        self.extras, self.seconds, self.harness_error = {}, 0.0, None

    def summary(self):
        out = dict(row=self.row_id, status=self.status, seconds=round(self.seconds, 4))
        if self.exc is not None:
            out.update(exception=type(self.exc).__name__, reason=getattr(self.exc, "reason", None),
                       message=str(self.exc)[:400])
        if self.harness_error:
            out["harness_error"] = self.harness_error[-800:]
        return out


def classify(exc):
    from xpsloss import NumericalConvergenceError
    return "not_converged" if isinstance(exc, NumericalConvergenceError) else "error"


class patched:
    """Record the calls of _call_kernel and the deliveries of on_final_panels (by wrapping _run)."""

    def __init__(self, out):
        self.out = out

    def __enter__(self):
        from xpsloss import xps_single_loss as xs
        self.xs, self.saved = xs, (xs._call_kernel, xs._run)
        original_kernel, original_run, out = xs._call_kernel, xs._run, self.out

        def call_kernel(material_, electron_energy_eV, loss_eV, **kw):
            out.calls.append(dict(loss_eV=np.array(loss_eV, float, ndmin=1).copy(), kwargs=dict(kw),
                                  electron_energy_eV=float(electron_energy_eV)))
            return original_kernel(material_, electron_energy_eV, loss_eV, **kw)

        def run(*args, **kw):
            out.run_calls += 1
            callback = kw.get("on_final_panels")
            if callback is not None:
                def observed(result, panels):
                    out.deliveries.append((result, panels))
                    return callback(result, panels)
                kw["on_final_panels"] = observed
            return original_run(*args, **kw)
        xs._call_kernel, xs._run = call_kernel, run
        return self

    def __exit__(self, *exc):
        self.xs._call_kernel, self.xs._run = self.saved
        return False


def by_order(panels_by_order, n, total=None):
    if isinstance(panels_by_order, dict):
        return panels_by_order[n]
    if total is not None and len(panels_by_order) == total:
        return panels_by_order[n - 1]
    return panels_by_order[n]


# ---------------------------------------------------------------- entries
def analytic_kwargs(row, inputs=None, settings=None):
    i = row["inputs"] if inputs is None else inputs
    s = settings_of(row) if settings is None else settings
    d = table()["defaults"]["analytic"]
    spec = table()["kernels"][i["kernel"]]
    window = tuple(num(x) for x in i.get("window_override_eV", spec["window_eV"]))
    charge = float(i.get("hole_charge", d["hole_charge"]))
    features = features_of(spec)
    if "window_override_eV" in i:
        features["jumps_eV"] = tuple(window)
    kw = dict(features=features, weight_length_A=float(d["weight_length_A"]), hole_charge=charge,
              model=i["model"], max_order=i["max_order"], epsrel=s["epsrel"], area_epsabs=s["area_epsabs"],
              spectrum_epsabs_fraction=s["spectrum_epsabs_fraction"], panel_limit=s["panel_limit"],
              max_cells=s["max_cells"], order_epsrel=s["order_epsrel"],
              order_epsabs_fraction=s["order_epsabs_fraction"], order_limit=s["order_limit"],
              max_convolution_points=s["max_convolution_points"])
    if "primary" in i:
        kw.update(primary=make_primary(i["primary"]), output_edges_eV=np.array(edge_values(i["output_edges_eV"])))
    if "neglect_tolerance_fraction" in i:
        kw["neglect_tolerance_fraction"] = num(i["neglect_tolerance_fraction"])
    for name in ("multiplier", "multiplier_A"):
        if name in i:
            kw[name] = i[name]
    if "initial_edges" in i:
        kw["initial_edges"] = [float(x) for x in i["initial_edges"]]
    for name in ("order_panel_limit", "order_min_panel_width_eV"):
        if name in s:
            kw[name] = s[name]
    return AnalyticKernel(i["kernel"], charge), window, kw


def run_analytic(row, out):
    lo = importlib.import_module("xpsloss._loss_orders")
    kernel, window, kw = analytic_kwargs(row)
    i = row["inputs"]
    if "inject" in i:
        kw["inject"] = out.extras["inject"] = make_inject(i["inject"])
    if "order_inject" in i:
        kw["order_inject"] = out.extras["order_inject"] = make_order_inject(row)
    out.kernel = kernel
    out.extras["call"] = (window, kw)
    out.result, out.panels = lo.integrate_orders(kernel, window, **kw)


def run_reference_model(row, out):
    ref = importlib.import_module("xpsloss._loss_reference_models")
    kernel, window, kw = analytic_kwargs(row)
    out.kernel = kernel
    got = ref.integrate_orders_reference(kernel, window, **kw)
    out.result, out.panels = (got[0], got[1]) if isinstance(got, tuple) else (got, None)


def run_layer_c(row, out):
    lo = importlib.import_module("xpsloss._loss_orders")
    from xpsloss._loss_panels import Panel, _nodes
    i = row["inputs"]
    piece, = table()["kernels"][i["panels"]["kernel"]]["columns"]["term_13a"]
    window = tuple(num(x) for x in i["window"])
    nodes = _nodes(*window)[0]
    values, errors = np.zeros((11, 17)), np.zeros((11, 17))
    values[[0, 8, 10]] = np.array([legendre_value(piece, x) for x in nodes])
    panel = Panel.from_lobatto(window[0], window[1], values, errors)
    features = dict(jumps_eV=window, kinks_eV=(), singularities=(), min_half_width_eV=None,
                    min_variation_scale_eV=None)
    out.panels = lo.convolve_panels([panel], window, model=i["model"], max_order=i["max_order"], features=features)
    out.result = out.panels


def api_kwargs(row):
    i, d = row["inputs"], table()["defaults"]["public"]
    kw = dict(loss_window_eV=tuple(num(x) for x in i["loss_window_eV"]),
              weight_length_A=i.get("weight_length_A", d["weight_length_A"]),
              epsilon_momentum=i.get("epsilon_momentum", d["epsilon_momentum"]))
    for name in ("hole_charge", "multiplier", "multiplier_A", "weight_length_kind", "extra_inverse_length_per_A"):
        if name in i:
            kw[name] = i[name]
    first = dict(kw)
    for name, value in row["settings"].items():
        value = num(value)
        kw[name] = value
        if name in FIRST_ORDER_NAMES:
            first[name] = value
    kw.update(model=i["model"], max_order=i["max_order"])
    if "neglect_tolerance_fraction" in i:
        kw["neglect_tolerance_fraction"] = num(i["neglect_tolerance_fraction"])
    return kw, first


def run_api(row, out):
    import xpsloss
    i = row["inputs"]
    kw, first = api_kwargs(row)
    mat = material(i["material"])
    energy = float(table()["defaults"]["public"]["electron_energy_eV"])
    out.extras["first_order_call"] = (mat, energy, first)
    with patched(out):
        if row["entry"] == "api_areas":
            out.result = xpsloss.multiple_loss_order_areas(mat, energy, **kw)
        else:
            out.result = xpsloss.multiple_loss_spectrum(mat, energy, make_primary(i["primary"]),
                                                        np.array(edge_values(i["output_edges_eV"])), **kw)


def old_entry(row, out):
    """The old entry called in the same process with the same first-order arguments (recorded separately)."""
    import xpsloss
    mat, energy, first = out.extras["first_order_call"]
    other = Outcome(row["id"] + ".old")
    with patched(other):
        if row["entry"] == "api_areas":
            first = {k: v for k, v in first.items() if k not in ("spectrum_epsabs_fraction", "max_cells")}
            other.result = xpsloss.loss_window_area(mat, energy, **first)
        else:
            i = row["inputs"]
            other.result = xpsloss.single_loss_spectrum(mat, energy, make_primary(i["primary"]),
                                                        np.array(edge_values(i["output_edges_eV"])), **first)
    return other


def cli_inputs(row, folder):
    i = row["inputs"]
    cfg = folder/"input.json"
    cfg.write_bytes(json.dumps(dict(material=material_dict(i["config_material"]),
                                    electron_energy_eV=i["config_electron_energy_eV"],
                                    loss_grid_eV={"min": 1, "max": 100, "step": 1},
                                    integration=i["config_integration"])).encode("ascii"))
    s = i["spectrum"]
    kind, x, strength = primary_arrays(s["primary"])
    primary = {"kind": kind, "energies_eV" if kind == "delta" else "edges_eV": x, "intensities": strength,
               "intensity_unit": "arb. unit", "description": None, "included_broadening": None, "truncation": None,
               "intrinsic_satellites": None}
    spectrum = folder/"spectrum.json"
    spectrum.write_bytes(json.dumps({"schema": s["schema"], "axis": s["axis"], "primary": primary,
                                     "output_edges_eV": edge_values(s["output_edges_eV"])}).encode("ascii"))
    out = folder/"output"
    argv = ["xps-multiple-loss", str(cfg), "--spectrum", str(spectrum), *i["argv"], "--out", str(out)]
    return argv, cfg, spectrum, out


def run_cli(row, out):
    from xpsloss import cli
    with tempfile.TemporaryDirectory() as tmp:
        if row["entry"] == "cli_capabilities":
            argv, cfg, spectrum, folder = list(row["inputs"]["argv"]), None, None, Path(tmp)/"none"
        else:
            argv, cfg, spectrum, folder = cli_inputs(row, Path(tmp))
        before = None
        if row["inputs"].get("output_dir") == "existing":
            folder.mkdir()
            (folder/"keep.txt").write_bytes(b"existing content\n")
            before = {p.name: p.read_bytes() for p in folder.iterdir()}
        stdout, stderr = io.StringIO(), io.StringIO()
        with patched(out), redirect_stdout(stdout), redirect_stderr(stderr):
            try:
                code = cli.main(argv)
            except SystemExit as exc:
                code = exc.code
        after = {p.name: p.read_bytes() for p in folder.iterdir()} if folder.exists() else None
        state = "absent" if after is None else ("unchanged" if before is not None and after == before else "created")
        record = dict(exit_code=code, stdout=stdout.getvalue(), stderr=stderr.getvalue(), output_dir_after=state,
                      argv=argv, files=after or {},
                      config_sha256=hashlib.sha256(cfg.read_bytes()).hexdigest() if cfg else None,
                      spectrum_sha256=hashlib.sha256(spectrum.read_bytes()).hexdigest() if spectrum else None)
        out.cli = record
        try:
            record["stderr_json"] = json.loads(record["stderr"]) if record["stderr"].strip() else None
        except ValueError:
            record["stderr_json"] = "not JSON"
        if code == 0 and "run.json" in record["files"]:
            record["manifest"] = json.loads(record["files"]["run.json"].decode("utf-8"))
            if row["id"] == "C01":
                try:
                    out.extras["api_replay"] = cli_replay(record["manifest"])
                except Exception:
                    out.extras["api_replay_error"] = traceback.format_exc()
    if code == 0:
        out.status = "complete"
    else:
        payload = record["stderr_json"] if isinstance(record["stderr_json"], dict) else {}
        out.status = payload.get("status", "error")


def cli_replay(manifest):
    """The API call rebuilt from run.json alone (prereg 4.3)."""
    import xpsloss
    d = manifest["spectrum_definition"]
    mat = manifest["config"]["material"]
    m = xpsloss.Material(mat["name"], mat["gap_eV"], tuple(xpsloss.Oscillator(**o) for o in mat["oscillators"]))
    kw = dict(loss_window_eV=tuple(d["loss_window_eV"]), weight_length_A=d["weight"]["length_A"],
              epsilon_momentum=d["epsilon_momentum"], hole_charge=d["hole_charge"], theta_deg=d["theta_deg"],
              model=manifest["model_id"], max_order=d["max_order"], **manifest["effective_settings"])
    p = manifest["primary"]
    primary = (xpsloss.DeltaPrimary(np.array(p["energies_eV"]), np.array(p["intensities"])) if "energies_eV" in p
               else xpsloss.BinnedPrimary(np.array(p["edges_eV"]), np.array(p["intensities"])))
    return xpsloss.multiple_loss_spectrum(m, d["kernel_energy_eV"], primary, np.array(manifest["output_edges_eV"]), **kw)


RUNNERS = {"analytic_spectrum": run_analytic, "analytic_areas": run_analytic, "layer_c": run_layer_c,
           "reference_model": run_reference_model, "api_spectrum": run_api, "api_areas": run_api,
           "cli": run_cli, "cli_capabilities": run_cli}


def run_row(rid):
    """Run a row once (cached: reuse of row.result by later rows reads the same outcome)."""
    runs = _CACHE.setdefault("runs", {})
    if rid in runs:
        return runs[rid]
    row = rows()[rid]
    out = Outcome(rid)
    start = time.perf_counter()
    try:
        RUNNERS[row["entry"]](row, out)
        if out.status is None:
            out.status = "success"
    except (ImportError, AttributeError, TypeError, Unreadable, KernelContractError) as exc:
        # a missing module or name of the surface, or a call that does not fit prereg 4: a failure of the row
        out.exc, out.status = exc, "harness_error"
        out.harness_error = traceback.format_exc()
    except Exception as exc:
        out.exc = exc
        try:
            out.status = classify(exc)
        except Exception:
            out.status, out.harness_error = "harness_error", traceback.format_exc()
    out.seconds = time.perf_counter() - start
    runs[rid] = out
    return out


# ---------------------------------------------------------------- reading the result
def order_of(result, n):
    for element in G(result, "orders"):
        if G(element, "order") == n:
            return element
    raise Unreadable(f"order {n} is not in result.orders")


def read_region(result, n, region, comp):
    o = order_of(result, n)
    if region.startswith("bin["):
        i = int(region[4:-1])
        return float(G(o, comp)[i]), float(G(o, comp + "_error")[i])
    side = G(G(o, "outside"), region.split(".")[1])
    return float(G(side, comp)), float(G(side, comp + "_error"))


def read_region_part(result, n, region, comp, part):
    o = order_of(result, n)
    if region.startswith("bin["):
        return float(G(G(o, comp + "_error_parts"), part)[int(region[4:-1])])
    return float(G(G(G(G(o, "outside"), region.split(".")[1]), comp + "_error_parts"), part))


def read_area(result, n, col, name="area"):
    o = order_of(result, n)
    return float(G(G(o, name), col)), float(G(G(o, name + "_error"), col))


def read_spectrum(result, region):
    if region.startswith("bin["):
        i = int(region[4:-1])
        return float(G(result, "spectrum")[i]), float(G(result, "spectrum_error")[i])
    side = G(G(result, "outside"), region.split(".")[1])
    return float(G(side, "spectrum")), float(G(side, "spectrum_error"))


def read_zero_loss(result, region):
    if region.startswith("bin["):
        i = int(region[4:-1])
        return float(G(result, "zero_loss")[i]), float(G(result, "zero_loss_error")[i])
    side = G(G(G(result, "single_loss"), "outside"), region.split(".")[1])
    return float(G(side, "zero_loss")), float(G(side, "zero_loss_error"))


def primary_sums(row):
    _, _, s = primary_arrays(row["inputs"]["primary"])
    return sum(Fr(v) for v in s), sum(abs(Fr(v)) for v in s)


def multiplier_of(result):
    single = G(result, "single_loss")
    m = float(G(single, "multiplier_A"))
    interval = G(single, "multiplier_interval_A")
    return m, (None if interval is None else (float(interval[0]), float(interval[1])))


def warning_codes(result):
    return [G(w, "code") for w in G(G(result, "diagnostics"), "warnings")]


def canonical(x):
    """Bit-exact, comparable form of a result (arrays by dtype, shape and bytes; floats by hex)."""
    if isinstance(x, np.ndarray):
        return ["ndarray", x.dtype.str, list(x.shape), hashlib.sha256(np.ascontiguousarray(x).tobytes()).hexdigest()]
    if isinstance(x, (bool, np.bool_)):
        return bool(x)
    if isinstance(x, (float, np.floating)):
        return ["float", float(x).hex()]
    if isinstance(x, (int, np.integer)):
        return int(x)
    if dataclasses.is_dataclass(x) and not isinstance(x, type):
        return {"@type": type(x).__name__, **{f.name: canonical(getattr(x, f.name)) for f in dataclasses.fields(x)}}
    if isinstance(x, dict):
        return {str(k): canonical(v) for k, v in x.items()}
    if isinstance(x, tuple):
        return ["tuple", *[canonical(v) for v in x]]
    if isinstance(x, list):
        return [canonical(v) for v in x]
    if isinstance(x, BaseException):
        return {"@exception": type(x).__name__, "message": str(x),
                "attributes": canonical({k: v for k, v in vars(x).items()})}
    return x


def first_difference(a, b, path=""):
    if type(a) is not type(b):
        return f"{path}: type {type(a).__name__} != {type(b).__name__}"
    if isinstance(a, dict):
        if list(a) != list(b):
            return f"{path}: keys {sorted(set(a) ^ set(b)) or 'order differs'}"
        for k in a:
            d = first_difference(a[k], b[k], f"{path}.{k}")
            if d:
                return d
        return None
    if isinstance(a, list):
        if len(a) != len(b):
            return f"{path}: length {len(a)} != {len(b)}"
        for i, (x, y) in enumerate(zip(a, b)):
            d = first_difference(x, y, f"{path}[{i}]")
            if d:
                return d
        return None
    return None if a == b else f"{path}: {a!r} != {b!r}"


def canonical_calls(calls):
    return [dict(electron_energy_eV=float(c["electron_energy_eV"]).hex(), loss_eV=canonical(c["loss_eV"]),
                 kwargs=canonical(dict(sorted(c["kwargs"].items())))) for c in calls]


# ---------------------------------------------------------------- baseline child (tag v0.3.0)
M4A_FIRST = {"M4A-o3": "M4A-c1", "M4A-o3s": "M4A-c4"}


def baseline_inputs(rid):
    """First-order arguments of an M4A row as plain data (no object of the package)."""
    row = rows()[rid]
    kw, first = api_kwargs(row)
    i = row["inputs"]
    return dict(entry=row["entry"], material=i["material"], first=first,
                primary=i.get("primary"), edges=i.get("output_edges_eV"))


def _baseline_child(out_path, payload_path):
    """Runs inside the unpacked tag: old entries only."""
    import xpsloss
    from xpsloss import cli, xps_single_loss as xs
    payload = json.loads(Path(payload_path).read_text(encoding="utf-8"))
    record = dict(file=str(Path(xpsloss.__file__).resolve()), module=str(Path(xs.__file__).resolve()),
                  version=xpsloss.__version__, executable=sys.executable, numpy=np.__version__,
                  has_new_api=hasattr(xpsloss, "multiple_loss_spectrum"), rows={})
    energy = payload["electron_energy_eV"]
    for rid, item in payload["rows"].items():
        m = item["material_data"]
        mat = xpsloss.Material(item["material"], float.fromhex(m["gap_hex"]),
                               tuple(xpsloss.Oscillator(*(float.fromhex(h) for h in o)) for o in m["oscillators_hex"]))
        calls, original = [], xs._call_kernel

        def recorder(material_, electron_energy_eV, loss_eV, **kw):
            calls.append(dict(loss_eV=np.array(loss_eV, float, ndmin=1).copy(), kwargs=dict(kw),
                              electron_energy_eV=float(electron_energy_eV)))
            return original(material_, electron_energy_eV, loss_eV, **kw)
        first = dict(item["first"])
        first["loss_window_eV"] = tuple(first["loss_window_eV"])
        xs._call_kernel = recorder
        try:
            if item["entry"] == "api_areas":
                first = {k: v for k, v in first.items() if k not in ("spectrum_epsabs_fraction", "max_cells")}
                result = xpsloss.loss_window_area(mat, energy, **first)
            else:
                kind, x, s = item["primary_data"]
                primary = (xpsloss.DeltaPrimary if kind == "delta" else xpsloss.BinnedPrimary)(
                    np.array(x, float), np.array(s, float))
                result = xpsloss.single_loss_spectrum(mat, energy, primary, np.array(item["edges_data"]), **first)
            snapshot = canonical(result)
        except Exception as exc:
            snapshot = canonical(exc)
        finally:
            xs._call_kernel = original
        record["rows"][rid] = dict(result=snapshot, calls=canonical_calls(calls))
    stdout = io.StringIO()
    with redirect_stdout(stdout):
        try:
            code = cli.main(["capabilities"])
        except SystemExit as exc:
            code = exc.code
    record["capabilities"] = dict(exit_code=code, stdout=stdout.getvalue())
    Path(out_path).write_text(json.dumps(record), encoding="utf-8")


def baseline():
    """Snapshots of the baseline (tag v0.3.0) for the selected M4A rows and C11: one child process."""
    if "baseline" in _CACHE:
        return _CACHE["baseline"]
    box = _CACHE["baseline"] = {"unavailable": None, "seconds": None}
    if not (ROOT/".git").exists():
        box["unavailable"] = ("skipped", "no git history (unpacked sdist): the baseline comparison is not counted")
        return box
    try:
        commit = subprocess.check_output(["git", "rev-parse", BASELINE_TAG + "^{commit}"], cwd=ROOT, text=True,
                                         stderr=subprocess.STDOUT).strip()
    except Exception as exc:
        box["unavailable"] = ("fail", f"tag {BASELINE_TAG} is not available in the source tree: {exc}")
        return box
    if commit != BASELINE_COMMIT:
        box["unavailable"] = ("fail", f"tag {BASELINE_TAG} is {commit}, not the registered commit")
        return box
    ids = [rid for rid in selected_ids() if rows()[rid]["group"] == "M4A"]
    wanted = sorted({M4A_FIRST.get(rid, rid) for rid in ids} | ({"M4A-c1"} if "C11" in selected_ids() else set()))
    payload = dict(electron_energy_eV=float(table()["defaults"]["public"]["electron_energy_eV"]), rows={})
    for rid in wanted:
        item = baseline_inputs(rid)
        item["material_data"] = {k: table()["materials"][item["material"]][k] for k in ("gap_hex", "oscillators_hex")}
        if item["primary"] is not None:
            item["primary_data"] = list(primary_arrays(item["primary"]))
            item["edges_data"] = edge_values(item["edges"])
        payload["rows"][rid] = item
    start = time.perf_counter()
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        raw = subprocess.check_output(["git", "archive", "--format=zip", BASELINE_COMMIT, "src/xpsloss"], cwd=ROOT)
        with zipfile.ZipFile(io.BytesIO(raw)) as z:
            z.extractall(folder)
        src = folder/"src"
        (folder/"payload.json").write_text(json.dumps(payload), encoding="utf-8")
        env = dict(os.environ, PYTHONPATH=str(src))
        done = subprocess.run([sys.executable, str(Path(__file__).resolve()), "--baseline-child",
                               str(folder/"out.json"), str(folder/"payload.json")], cwd=folder, env=env,
                              capture_output=True, text=True)
        if done.returncode != 0:
            box["unavailable"] = ("fail", f"baseline child failed: {done.stderr[-600:]}")
            return box
        record = json.loads((folder/"out.json").read_text(encoding="utf-8"))
        inside = Path(record["file"]).resolve().is_relative_to(src.resolve()) and \
            Path(record["module"]).resolve().is_relative_to(src.resolve())
    box.update(record=record, seconds=time.perf_counter() - start, from_archive=inside)
    if not inside or record["has_new_api"] or record["version"] != "0.3.0" or record["executable"] != sys.executable:
        box["unavailable"] = ("fail", f"baseline child environment: {record['file']}, version {record['version']}")
    return box


# ---------------------------------------------------------------- judgements
class Judgement:
    def __init__(self):
        self.items = {}

    def add(self, rid, key, kind, status, detail=""):
        if (rid, key, kind) in self.items:
            raise AssertionError(f"duplicate judgement {(rid, key, kind)}")
        self.items[(rid, key, kind)] = (status, detail)

    def counts(self):
        out = {}
        for (_, _, kind), (status, _) in self.items.items():
            out.setdefault(kind, {}).setdefault(status, 0)
            out[kind][status] += 1
        return out

    def failures(self, rid):
        return {k: v for k, v in self.items.items() if k[0] == rid and v[0] != "pass"}


def need(condition, message):
    if not condition:
        raise AssertionError(message)


def guard(check):
    """Run a check; (status, detail).  An unreadable field or any exception is a failure."""
    try:
        got = check()
        if got is None or got is True:
            return "pass", ""
        if isinstance(got, tuple):
            return got
        return ("pass", "") if got else ("fail", "condition is false")
    except AssertionError as exc:
        return "fail", str(exc)[:500]
    except Exception as exc:
        return "fail", f"{type(exc).__name__}: {str(exc)[:400]}"


def success_result(rid):
    out = run_row(rid)
    need(out.status == "success" and out.result is not None, f"{rid} did not succeed ({out.status}): reuse impossible")
    return out.result


def truncation(result):
    return G(result, "truncation")


def bound_triple(result, i):
    t = truncation(result)
    return (G(t, "neglected_at_or_above_region_bound")[i], G(t, "bound_reason")[i], G(t, "tau_per_eV")[i])


def is_complete_bound(result, i):
    b, reason, tau = bound_triple(result, i)
    need(isinstance(b, float) and b == 0.0 and reason == "complete" and tau is None,
         f"bin[{i}]: {(b, reason, tau)} is not (0.0, 'complete', None)")


def kernel_evaluations(row, out):
    if row["entry"].startswith("analytic") or row["entry"] == "reference_model":
        return None if out.kernel is None else out.kernel.calls
    return len(out.calls)


def check_state(row, out):
    e = row["expect"]
    need(out.status != "harness_error", f"harness error: {(out.harness_error or '')[-500:]}")
    if row["entry"] in ("cli", "cli_capabilities"):
        return check_cli_state(row, out)
    need(out.status == e["status"], f"status {out.status} ({out.summary().get('reason')}: "
                                    f"{out.summary().get('message')}), expected {e['status']}")
    count = kernel_evaluations(row, out)
    if "kernel_evaluations" in e:
        if e["kernel_evaluations"] == 0:
            need(count == 0, f"{count} kernel evaluations, expected none")
        else:
            need(count is not None and count >= 1, f"{count} kernel evaluations, expected at least 1")
    if e["status"] != "success":
        exc = out.exc
        need(type(exc).__name__ == e["exception"], f"exception {type(exc).__name__}, expected {e['exception']}")
        need(G(exc, "reason") == e["reason"], f"reason {G(exc, 'reason')}, expected {e['reason']}")
        for name in ("requested_order", "order_limit", "limit", "lower_bound", "resource", "order", "stage",
                     "dominant_part", "argument"):
            if name in e:
                need(G(exc, name) == e[name], f"{name} = {G(exc, name)!r}, expected {e[name]!r}")
        if e["exception"].startswith("MultipleLoss"):
            import xpsloss
            base = ValueError if e["status"] == "error" else xpsloss.NumericalConvergenceError
            need(isinstance(exc, getattr(xpsloss, e["exception"])) and isinstance(exc, base), "exception base class")
        if "required_order_estimate" in e:
            need(G(exc, "required_order_estimate") >= 5, f"required_order_estimate {G(exc, 'required_order_estimate')}")
        if e["status"] == "not_converged" and e["exception"] == "MultipleLossConvergenceError":
            need(G(exc, "stage") in ("after_first_order", "order_refinement", "order_fold", "truncation"), "stage")
            inner = G(exc, "inner")
            need(inner is None or {"type", "reason", "stage"} <= set(inner), "inner")
            json.dumps(inner)
            if e["reason"].startswith("truncation_bound"):
                need(G(exc, "stage") == "truncation", "truncation stage")
        if "exceeded_includes" in e:
            check_exceeded(row, exc, e["exceeded_includes"])
        if row["id"] == "N15":
            p1 = G(G(G(exc, "diagnostics"), "first_order"), "final_panels")
            need(p1 >= 9 and G(exc, "estimate") == 280.5*sum(n*(p1 + n - 1) for n in range(2, 7))
                 and G(exc, "limit") == 50000 and G(exc, "estimate") >= 70125, "N15: estimate and P1")
        if row["id"] == "N16":
            need(G(exc, "measured") > 4500, "N16: measured count")
            diag = G(G(exc, "diagnostics"), "first_order")
            if G(diag, "final_panels") != 6:
                return "not_run", "precondition: 6 first-order panels"
        if row["id"] in ("N18", "N21"):
            panel = G(exc, "panel_eV")
            need(panel is not None and len(panel) == 2 and panel[0] < panel[1], "panel_eV")
        return None
    result = out.result
    if row["entry"] == "layer_c":
        return None
    if e.get("reference_only"):
        need(G(result, "reference_only") is True, "reference_only")
        return None
    if "truncation_state" in e:
        need(G(truncation(result), "state") == e["truncation_state"],
             f"state {G(truncation(result), 'state')}, expected {e['truncation_state']}")
    for name in ("support_order", "computed"):
        if name in e:
            need(G(truncation(result), name) == e[name], f"{name} {G(truncation(result), name)}, expected {e[name]}")
    if e.get("orders_empty"):
        need(G(result, "orders") == (), "orders is not the empty tuple")
    codes = warning_codes(result)
    need(set(ALWAYS_WARNINGS) <= set(codes), f"missing constant warnings: {sorted(set(ALWAYS_WARNINGS) - set(codes))}")
    need(len(G(result, "warnings")) >= len(ALWAYS_WARNINGS) and all(isinstance(w, str) for w in G(result, "warnings")),
         "warnings is not a tuple of messages")
    if row["entry"] != "analytic_areas" and row["entry"] != "api_areas":
        state = G(truncation(result), "state")
        need(("truncated" in codes) == (state == "truncated"), "the warning 'truncated' does not follow the state")
    for code in e.get("warnings_include", ()):
        need(code in codes, f"warning {code} is missing")
    for code in e.get("warnings_exclude", ()):
        need(code not in codes, f"warning {code} must not appear")
    if "neglected_total_area" in e:
        got = G(truncation(result), "neglected_total_area")
        want = e["neglected_total_area"]
        need(G(got, "reason") == want["reason"], f"neglected_total_area reason {G(got, 'reason')}")
        value = G(got, "value")
        if want["value"] is None:
            need(value is None, "neglected_total_area value must be None")
        elif want["value"] == "finite":
            need(isinstance(value, float) and math.isfinite(value) and value > 0, "neglected_total_area value")
        else:
            need(isinstance(value, float) and value == 0.0, "neglected_total_area value must be 0.0")
        G(got, "interval")
    if "convention_ratio" in e:
        need(G(result, "convention_ratio") == e["convention_ratio"], f"convention_ratio {G(result, 'convention_ratio')}")
    if "beta_deliveries" in e:
        need(G(G(G(result, "diagnostics"), "first_order"), "beta_deliveries") == 1 and len(out.deliveries) == 1
             and out.run_calls == 1, f"deliveries {len(out.deliveries)}, _run calls {out.run_calls}")
    if "neglected_bound" in e:
        n_bins = len(G(truncation(result), "neglected_at_or_above_region_bound"))
        if e["neglected_bound"].startswith("exactly 0"):
            for i in range(n_bins):
                b = bound_triple(result, i)[0]
                need(isinstance(b, float) and b == 0.0, f"bound of bin[{i}] is {b!r}")
        else:
            for i in range(4):
                need(bound_triple(result, i) == (None, "no_usable_tau", None), f"T5n bin[{i}]: {bound_triple(result, i)}")
            for i in (4, 5):
                is_complete_bound(result, i)
    if row["entry"] in ("analytic_spectrum", "analytic_areas"):
        need(G(G(result, "single_loss"), "electron_energy_eV") is None, "analytic entry: electron_energy_eV")
    for name, value in (("model_id", table()["model_id"]), ("kernel_energy_policy", "fixed_E0"),
                        ("repeated_kernel", "term_13a"), ("hole_treatment", "first_event_only")):
        need(G(result, name) == value, f"{name} = {G(result, name)!r}")
    if row["id"] == "RB1a":
        for n in (1, 2, 3):
            mesh = G(order_of(result, n), "mesh")
            need(G(mesh, "cells") == 0 and G(mesh, "formal_folds") == 0, "area entry: mesh.cells and formal_folds")
    if row["id"] == "A3c":
        need("negative_total_kernel" in codes, "negative_total_kernel")
    if row["id"] == "A3d":
        need("multiplier_differs_from_weight_length" in codes, "multiplier_differs_from_weight_length")
    return None


def box_convolution_bin(height, lo, hi, s0, s1):
    """int_{s0}^{s1} (height*1[lo,hi] * 1[lo,hi])(s) ds, exact."""
    def cumulative(s):       # int of the triangle height*(min(s, lo+hi) - 2 lo ... ) by truncated powers
        f = lambda x: (max(x, Fr(0))**2)/2
        return height*(f(s - 2*lo) - 2*f(s - lo - hi) + f(s - 2*hi))
    return cumulative(s1) - cumulative(s0)


def check_exceeded(row, exc, want):
    listed = G(exc, "exceeded")
    hits = [x for x in listed if str(G(x, "key")).endswith(want["key"])]
    need(len(hits) == 1, f"exceeded has {len(hits)} entries for {want['key']}")
    hit = hits[0]
    need(G(hit, "side") == want["side"] and G(hit, "suggestion") == want["suggestion"],
         f"side {G(hit, 'side')}, suggestion {G(hit, 'suggestion')}")
    need(not any("area" in str(G(x, "key")) for x in listed), "an area key is listed in exceeded")
    # the numbers of the listed element (prereg 5.2, N19 and N20): K_SGN, |hole| = 0.04 on [2, 5]
    s = settings_of(row)
    e = Fr(num(row["inputs"]["inject"]["errors"]["term_13a"]))
    m = Fr(2)
    edges = [Fr(x) for x in edge_values(row["inputs"]["output_edges_eV"])]
    i = int(re.search(r"bin\[(\d+)\]", want["key"]).group(1))
    s0, s1 = 1000 - edges[i + 1], 1000 - edges[i]
    envelope = m*m*e*box_convolution_bin(Fr(0.04), Fr(2), Fr(5), s0, s1)
    need(Fr(G(hit, "inherited")) >= envelope*(1 - Fr(1, 10**9)), "inherited is below the envelope integral")
    phi2 = m*Fr(0.1)*3 + m*Fr(0.04)*3
    absolute = phi2*Fr(s["spectrum_epsabs_fraction"])
    need(abs(Fr(G(hit, "budget_absolute")) - absolute) <= Fr(1, 10**6)*absolute + Fr(1, 10**300), "budget_absolute")
    if row["id"] == "N19":
        y = m*m*Fr(0.1)*Fr(0.04)*_sgn_bin(s0, s1)
        relative = 2*Fr(s["epsrel"])*abs(y)
        need(abs(Fr(G(hit, "budget_relative")) - relative) <= Fr(1, 10**6)*relative, "budget_relative")
        need(Fr(G(hit, "inherited")) > relative, "inherited does not exceed the relative budget")
    else:
        need(Fr(G(hit, "inherited")) > absolute, "inherited does not exceed the absolute budget")
    need(G(hit, "side") == ("relative" if Fr(G(hit, "budget_relative")) >= Fr(G(hit, "budget_absolute")) else "absolute"),
         "the side is not the larger budget")


def _sgn_bin(s0, s1):
    """int_{s0}^{s1} (sgn-box * 1[2,5])(s) ds for the hole of K_SGN in units of its height: +1 on [2,3.5], -1 on [3.5,5]."""
    f = lambda x: (max(x, Fr(0))**2)/2
    cumulative = lambda s: (f(s - 4) - f(s - 7)) - 2*(f(s - Fr(11, 2)) - f(s - Fr(17, 2))) + (f(s - 7) - f(s - 10))
    return cumulative(s1) - cumulative(s0)


def check_cli_state(row, out):
    e, record = row["expect"], out.cli
    need(record is not None, "the CLI did not run")
    need(record["exit_code"] == e["exit_code"], f"exit code {record['exit_code']}: {record['stderr'][-300:]}")
    need(out.status == e["status"], f"status {out.status}, expected {e['status']}")
    if row["entry"] == "cli_capabilities":
        return None
    need(record["output_dir_after"] == e["output_dir_after"], f"output directory {record['output_dir_after']}")
    need(record["stderr"].isascii() and record["stdout"].isascii(), "non-ASCII CLI output")
    if e["exit_code"] == 0:
        manifest = record["manifest"]
        need(sorted(record["files"]) == ["multiple_loss_spectrum.csv", "order_areas.csv", "run.json"], "output files")
        need(G(G(G(manifest, "results"), "truncation"), "state") == e["truncation_state"], "truncation state")
        return None
    payload = record["stderr_json"]
    need(isinstance(payload, dict), "stderr is not a JSON object")
    need(payload.get("error_type") == e["error_type"], f"error_type {payload.get('error_type')}")
    need(payload.get("status") == e["status"] and "message" in payload, "status and message")
    need(payload.get("reason") == e.get("reason"), f"reason {payload.get('reason')}")
    if row["id"] in ("C03", "C04"):
        option = "--model" if row["id"] == "C03" else "--max-order"
        need(option in payload["message"] and "invalid choice" not in payload["message"],
             f"the usage error does not name the missing option {option}: {payload['message'][:200]}")
    if row["id"] == "C05":
        for name in ("order", "stage", "inner", "required_order_estimate"):
            need(name in payload, f"stderr JSON lacks {name}")
        need(payload["stage"] == "truncation", "stage")
    if row["id"] == "C10":
        need(payload.get("stage") == "after_first_order" and payload.get("order") is None
             and payload.get("inner") is None and payload.get("limit") == 50000 and payload.get("estimate") is not None,
             "C10 attributes")
    return None


# ---------------------------------------------------------------- identities (prereg 7.4)
def check_components(row, result, n, slot):
    charge = fr(row["inputs"].get("hole_charge", 1.0))
    if slot == "area":
        (t, et), (nh, enh), (h, eh) = (read_area(result, n, c) for c in COMPONENTS)
    else:
        (t, et), (nh, enh), (h, eh) = (read_region(result, n, slot, c) for c in COMPONENTS)
    diff = abs(Fr(t) - Fr(nh) - charge*Fr(h))
    bound = Fr(et) + Fr(enh) + Fr(eh) + 8*FEPS*(abs(Fr(t)) + abs(Fr(nh)) + abs(charge*Fr(h)))
    need(diff <= bound, f"order {n} {slot}: |total - no_hole - Q*hole| = {float(diff):.3e} > {float(bound):.3e}")


def check_region_sum(row, result, n, comp):
    S, _ = primary_sums(row)
    regions = row_regions(row)
    values = [read_region(result, n, r, comp) for r in regions]
    w, ew = read_area(result, n, comp)
    total = sum(Fr(v) for v, _ in values)
    bound = sum(Fr(e) for _, e in values) + Fr(ew)*abs(S) + 8*FEPS*(sum(abs(Fr(v)) for v, _ in values) + abs(S*Fr(w)))
    need(abs(total - S*Fr(w)) <= bound, f"order {n} {comp}: sum over the regions differs from S*w_n by "
                                        f"{float(abs(total - S*Fr(w))):.3e} > {float(bound):.3e}")


def check_spectrum_sum(result, region):
    orders = [G(o, "order") for o in G(result, "orders")]
    terms = [read_zero_loss(result, region)] + [read_region(result, n, region, "total") for n in orders]
    y, e = read_spectrum(result, region)
    exact = sum(Fr(v) for v, _ in terms)
    scale = sum(abs(Fr(v)) for v, _ in terms)
    need(abs(Fr(y) - exact) <= 8*FEPS*scale, f"{region}: spectrum differs from J0 + sum of the orders")
    errors = sum(Fr(v) for _, v in terms)
    need(Fr(e) >= errors*(1 - 4*FEPS) and Fr(e) <= errors*(1 + 8*FEPS) + 8*FEPS*scale,
         f"{region}: spectrum_error {e} against the sum of the returned errors {float(errors)}")


def close_pair(a, b, label):
    (x, ex), (y, ey) = a, b
    need(abs(Fr(x) - Fr(y)) <= Fr(ex) + Fr(ey), f"{label}: |{x} - {y}| > {ex} + {ey}")


def enumerated(key):
    m = re.search(r" \[(\d+)/(\d+)\]$", key)
    return (int(m.group(1)) - 1, int(m.group(2))) if m else (0, 1)


def check_identity(row, out, key, oracle):
    result, rid = out.result, row["id"]
    need(out.status == "success", f"the row did not succeed ({out.status})")
    index, count = enumerated(key)
    m = re.match(r"order\[(\d+)\.\.(\d+)\]: total = no_hole", key)
    if m:
        lo, hi = int(m.group(1)), int(m.group(2))
        slots = row_regions(row) + ["area"]
        need(count == (hi - lo + 1)*len(slots), "ledger count of the component identity")
        return check_components(row, result, lo + index // len(slots), slots[index % len(slots)])
    m = re.match(r"order\[(\d+)\.\.(\d+)\]: sum over all regions", key)
    if m:
        lo = int(m.group(1))
        return check_region_sum(row, result, lo + index // 3, COMPONENTS[index % 3])
    if key.startswith("spectrum."):
        return check_spectrum_sum(result, key[len("spectrum."):])
    if rid == "RB1a":
        other = success_result("RB1")
        if index == 66:
            a, b = G(result, "order_area_ratio"), G(other, "order_area_ratio")
            return close_pair((G(a, "value"), G(a, "error")), (G(b, "value"), G(b, "error")), "order_area_ratio")
        name = "area" if index < 33 else "cutoff_sensitivity"
        n, col = 1 + (index % 33)//11, COLUMNS[index % 11]
        return close_pair(read_area(result, n, col, name), read_area(other, n, col, name), f"order {n} {name} {col}")
    if rid == "I-areas":
        other = success_result("A3a")
        n, col = 1 + index//11, COLUMNS[index % 11]
        return close_pair(read_area(result, n, col), read_area(other, n, col), f"order {n} area {col}")
    if rid == "REF2":
        other = success_result("REF1")
        regions = row_regions(row)
        n = 2 + index//(len(regions)*3)
        region, comp = regions[(index // 3) % len(regions)], COMPONENTS[index % 3]
        (y, e), (y1, e1) = read_region(result, n, region, comp), read_region(other, n, region, comp)
        f = math.factorial(n)
        need(abs(Fr(y) - Fr(y1)/f) <= Fr(e) + Fr(e1)/f + 8*FEPS*abs(Fr(y)), f"order {n} {region} {comp}: REF2 != REF1/n!")
        return None
    if rid == "CT3":
        other = success_result("M4A-o3")
        regions = row_regions(row)
        per = len(regions)*3 + 11
        n, j = 1 + index//per, index % per
        if j < len(regions)*3:
            region, comp = regions[j//3], COMPONENTS[j % 3]
            return close_pair(read_region(result, n, region, comp), read_region(other, n, region, comp), f"{n} {region} {comp}")
        col = COLUMNS[j - len(regions)*3]
        return close_pair(read_area(result, n, col), read_area(other, n, col), f"order {n} area {col}")
    if rid == "CT7":
        n = int(key[5])
        base = success_result("A2")
        a = num(table()["kernels"]["K_A2"]["window_eV"][0])
        a2 = num(row["inputs"]["window_override_eV"][0])
        checked = 0
        for col in COLUMNS:
            h = read_area(base, n, col, "cutoff_sensitivity")[0]
            if h == 0.0:
                continue
            fd = (read_area(base, n, col)[0] - read_area(result, n, col)[0])/math.log(a2/a)
            need(abs(h - fd) <= 0.02*abs(h), f"h_{n} {col}: {h} against the difference quotient {fd}")
            checked += 1
        need(checked >= 1, "no non-zero sensitivity column")
        return None
    if rid == "A3f" and key.startswith("control_external."):
        control = a3f_controls(out)["external"]
        m = re.match(r"control_external\.order\[(\d)\]\.area_error_parts\[(\w+)\]\[total\]", key)
        n, part = int(m.group(1)), m.group(2)
        a = float(G(G(G(order_of(result, n), "area_error_parts"), part), "total"))
        b = float(G(G(G(order_of(control, n), "area_error_parts"), part), "total"))
        need(abs(a - b) <= 4*EPS*max(abs(a), abs(b)), f"order {n} {part}: {a} (S2) against {b} (external)")
        need(float(G(G(G(order_of(control, n), "area_error_parts"), "multiplier_interval"), "total")) == 0.0,
             "the control has an interval part")
        return None
    raise AssertionError(f"no identity check for {key!r}")


# ---------------------------------------------------------------- exact zeros
def support_zero_regions(row):
    """(order, region) without positive-width overlap with [n a, n b] (A1s; exact)."""
    spec = table()["kernels"][row["inputs"]["kernel"]]
    a, b = (Fr(num(x)) for x in spec["window_eV"])
    edges = [Fr(x) for x in edge_values(row["inputs"]["output_edges_eV"])]
    p = Fr(primary_arrays(row["inputs"]["primary"])[1][0])
    out = []
    for n in range(1, row["inputs"]["max_order"] + 1):
        for name, lo, hi in [("outside.low", None, edges[0])] + \
                [(f"bin[{i}]", x, y) for i, (x, y) in enumerate(zip(edges, edges[1:]))] + [("outside.high", edges[-1], None)]:
            l1 = None if hi is None else p - hi
            l2 = None if lo is None else p - lo
            overlap = (l2 is None or l2 > n*a) and (l1 is None or l1 < n*b)
            if not overlap:
                out.append((n, name))
    return out


def zero_pair(pair, label):
    y, e = pair
    need(y == 0.0 and e == 0.0, f"{label}: value {y!r}, error {e!r} (both must be exactly 0)")


def check_zero(row, out, key, oracle):
    result, rid = out.result, row["id"]
    need(out.status == "success", f"the row did not succeed ({out.status})")
    index, count = enumerated(key)
    if rid == "A1s":
        cells = [(n, r, c) for n, r in support_zero_regions(row) for c in COMPONENTS]
        need(len(cells) == count == 57, f"{len(cells)} zero cells")
        n, region, comp = cells[index]
        return zero_pair(read_region(result, n, region, comp), f"order {n} {region} {comp}")
    if rid == "REF4a":
        if index < 7:
            b = bound_triple(result, index)[0]
            need(isinstance(b, float) and b == 0.0, f"bound of bin[{index}] is {b!r}")
            return None
        got = G(truncation(result), "neglected_total_area")
        need(G(got, "value") == 0.0 and G(got, "reason") is None, "neglected_total_area")
        return None
    if rid == "M4A-Z":
        regions = row_regions(row)
        per = len(regions)*3 + 11
        n, j = 1 + index//per, index % per
        if j < len(regions)*3:
            return zero_pair(read_region(result, n, regions[j//3], COMPONENTS[j % 3]), f"order {n}")
        col = COLUMNS[j - len(regions)*3]
        return zero_pair(read_area(result, n, col), f"order {n} area {col}")
    if rid in ("T2a", "T2c", "T2e"):
        i = 4 + index//3
        b, reason, tau = bound_triple(result, i)
        want = (index % 3 == 0 and isinstance(b, float) and b == 0.0) or (index % 3 == 1 and reason == "complete") \
            or (index % 3 == 2 and tau is None)
        need(want, f"bin[{i}]: {(b, reason, tau)}")
        return None
    if rid == "T4u" or key.startswith("truncation.bound."):
        return is_complete_bound(result, int(re.search(r"bin\[(\d+)\]", key).group(1)))
    if rid == "DS1":
        return is_complete_bound(result, index)
    if rid == "DS0":
        region = key[len("spectrum."):-len(" = zero_loss")]
        (y, e), (z, ez) = read_spectrum(result, region), read_zero_loss(result, region)
        need(bits(y) == bits(z) and bits(e) == bits(ez), f"{region}: spectrum {y} is not zero_loss {z}")
        return None
    if rid == "CT1":
        regions = row_regions(row)
        if index < 3*len(regions)*3:
            n = 1 + index//(len(regions)*3)
            j = index % (len(regions)*3)
            return zero_pair(read_region(result, n, regions[j//3], COMPONENTS[j % 3]), f"order {n}")
        return zero_pair(read_spectrum(result, regions[index - 3*len(regions)*3]), "spectrum")
    raise AssertionError(f"no zero check for {key!r}")


# ---------------------------------------------------------------- bitwise
def spectrum_snapshot(result):
    """Everything of orders and sums that I-shift compares (a shift of the axis must not change a bit)."""
    out = {"spectrum": canonical(G(result, "spectrum")), "spectrum_error": canonical(G(result, "spectrum_error")),
           "zero_loss": canonical(G(result, "zero_loss")), "zero_loss_error": canonical(G(result, "zero_loss_error")),
           "outside": canonical({s: {k: G(G(G(result, "outside"), s), k) for k in ("spectrum", "spectrum_error")}
                                 for s in ("low", "high")})}
    for o in G(result, "orders"):
        n = G(o, "order")
        for c in COMPONENTS:
            out[f"{n}.{c}"] = canonical(G(o, c))
            out[f"{n}.{c}_error"] = canonical(G(o, c + "_error"))
            for side in ("low", "high"):
                box = G(G(o, "outside"), side)
                out[f"{n}.{side}.{c}"] = canonical((G(box, c), G(box, c + "_error")))
        out[f"{n}.area"] = canonical((dict(G(o, "area")), dict(G(o, "area_error"))))
    return out


def order1_mapping(result, areas_only=False):
    """Order 1 is the transcription of the legacy fields (prereg 4.2)."""
    single = G(result, "single_loss")
    o = order_of(result, 1)
    need(bits(G(G(o, "area"), "total")) == bits(G(single, "single_loss_weight")), "area.total is not single_loss_weight")
    need(bits(G(G(o, "area_error"), "total")) == bits(G(single, "single_loss_weight_error")),
         "area_error.total is not single_loss_weight_error")
    for col in COLUMNS:
        need(bits(G(G(o, "cutoff_sensitivity"), col)) == bits(G(G(single, "cutoff_sensitivity"), col)), "h_1")
        need(bits(G(G(o, "cutoff_sensitivity_error"), col)) == bits(G(G(single, "cutoff_sensitivity_error"), col)), "h_1 error")
    m = Fr(float(G(single, "multiplier_A")))
    for col in COLUMNS[:-1]:
        want = float(m*Fr(float(G(G(single, "window_area"), col))))
        need(bits(G(G(o, "area"), col)) == bits(want), f"area[{col}] is not m*A rounded once")
    if areas_only:
        return
    for c in COMPONENTS:
        need(canonical(G(o, c)) == canonical(G(single, "single_loss_" + c)), f"order 1 {c}")
        need(canonical(G(o, c + "_error")) == canonical(G(single, "single_loss_" + c + "_error")), f"order 1 {c} error")
        for side in ("low", "high"):
            mine, old = G(G(o, "outside"), side), G(G(single, "outside"), side)
            need(bits(G(mine, c)) == bits(G(old, "single_loss_" + c))
                 and bits(G(mine, c + "_error")) == bits(G(old, "single_loss_" + c + "_error")), f"order 1 outside {side} {c}")
    need(canonical(G(result, "zero_loss")) == canonical(G(single, "zero_loss")), "zero_loss copy")
    if len(G(result, "orders")) == 1:
        need(canonical(G(result, "spectrum")) == canonical(G(single, "spectrum"))
             and canonical(G(result, "spectrum_error")) == canonical(G(single, "spectrum_error")), "N = 1 spectrum copy")


def check_bitwise(row, out, key, oracle):
    result, rid = out.result, row["id"]
    need(out.status == "success", f"the row did not succeed ({out.status})")
    if rid == "OI2":
        other = success_result("A1s")
        a = np.asarray(G(G(order_of(result, 2), "mesh"), "panel_edges_eV"), float)
        b = np.asarray(G(G(order_of(other, 2), "mesh"), "panel_edges_eV"), float)
        need(a.tobytes() == b.tobytes(), "order-2 mesh differs from that of A1s")
        return None
    if rid == "I-shift":
        d = first_difference(spectrum_snapshot(result), spectrum_snapshot(success_result("A1s")))
        need(d is None, f"shifted run differs from A1s: {d}")
        return None
    if oracle.endswith("same process"):
        old = old_entry(row, out)
        d = first_difference(canonical(G(result, "single_loss")), canonical(old.result))
        need(d is None, f"single_loss differs from the old entry: {d}")
        d = first_difference(canonical_calls(out.calls), canonical_calls(old.calls))
        need(d is None, f"_call_kernel calls differ from the old entry: {d}")
        order1_mapping(result, areas_only=row["entry"] == "api_areas")
        return None
    if oracle.startswith("baseline"):
        base = baseline()
        if base["unavailable"]:
            return base["unavailable"]
        snap = base["record"]["rows"][M4A_FIRST.get(rid, rid)]
        if key == "single_loss.all_legacy_fields":
            d = first_difference(canonical(G(result, "single_loss")), snap["result"])
            need(d is None, f"single_loss differs from the baseline: {d}")
            if rid in M4A_FIRST:
                d = first_difference(canonical(G(result, "single_loss")),
                                     canonical(G(success_result(M4A_FIRST[rid]), "single_loss")))
                need(d is None, f"single_loss differs from the max_order = 1 row: {d}")
            return None
        return order1_mapping(result, areas_only=key == "order[1].area")
    raise AssertionError(f"no bitwise check for {key!r}")


# ---------------------------------------------------------------- diagnostics
def mandatory_knots(row, n):
    spec = table()["kernels"][row["inputs"]["kernel"]]
    a, b = (Fr(num(x)) for x in spec["window_eV"])
    return sorted({float(k*a + (n - k)*b) for k in range(n + 1)})


def first_order_edges(out):
    panels = by_order(out.panels, 1, len(G(out.result, "orders")))
    return [float(G(panels[0], "left_eV"))] + [float(G(p, "right_eV")) for p in panels]


def a3f_controls(out):
    """The two in-row controls of A3f (run once, cached on the outcome)."""
    if "controls" in out.extras:
        return out.extras["controls"]
    lo = importlib.import_module("xpsloss._loss_orders")
    row = rows()["A3f"]
    m, _ = multiplier_of(out.result)
    kernel, window, kw = analytic_kwargs(row, {**row["inputs"], "multiplier": "external", "multiplier_A": m})
    kw["inject"] = make_inject(row["inputs"]["inject"])
    external, _ = lo.integrate_orders(kernel, window, **kw)
    c = row["inputs"]["controls"]["control_point_path"]
    inputs = {"kernel": c["kernel"], "model": row["inputs"]["model"], "max_order": c["max_order"],
              "primary": c["primary"], "output_edges_eV": c["output_edges_eV"], "initial_edges": c["initial_edges"]}
    kernel, window, kw = analytic_kwargs(row, inputs)
    kw["inject"] = make_inject(c["inject"])
    captured, original = [], lo._fold_order

    def wrapped(panels, geometry, *, multiplier_power, error_parts):
        captured.append((panels, error_parts))
        return original(panels, geometry, multiplier_power=multiplier_power, error_parts=error_parts)
    lo._fold_order = wrapped
    try:
        lo.integrate_orders(kernel, window, **kw)
    finally:
        lo._fold_order = original
    out.extras["controls"] = {"external": external, "point_path": captured}
    return out.extras["controls"]


def node_part_at(captured, energy, part="inherited", column="total"):
    """Unmultiplied node parts handed to _fold_order at a panel end node at ``energy`` (all occurrences)."""
    j, found = COLUMNS.index(column), []
    for panels, parts in captured:
        for panel, box in zip(panels, parts):
            if float(G(panel, "left_eV")) == energy:
                found.append(float(G(box, part)[j, 0]))
            if float(G(panel, "right_eV")) == energy:
                found.append(float(G(box, part)[j, 16]))
    return found


def check_diagnostic(row, out, key, oracle):
    rid, result = row["id"], out.result
    index, count = enumerated(key)
    if oracle in ("R-INV", "R-A"):
        notes = ra_references()["rows"][rid]["notes"]
        name = key
        need(name in notes and notes[name]["passes"] is True, "the reference does not record the control")
        if rid == "T2e":
            need(all(b["outside_width"] for b in notes[name]["bins"]), "T2e reference bins")
        else:
            need(notes[name]["outside_width_without_tail"], "INV reference bins")
        return None
    if row["entry"] in ("cli", "cli_capabilities"):
        return check_cli_diagnostic(row, out, key, index)
    if row["group"] == "M4A":
        return check_m4a_diagnostic(row, out, key, oracle)
    need(out.status == "success", f"the row did not succeed ({out.status})")
    if rid == "A1" and key.startswith("order[2..6].mesh.knots_eV"):
        n = 2 + index
        got = [float(x) for x in G(G(order_of(result, n), "mesh"), "knots_eV")]
        need([bits(x) for x in got] == [bits(x) for x in mandatory_knots(row, n)], f"order {n}: knots_eV {got}")
        return None
    if rid == "A1":
        target = float(2*Fr(1.3) + Fr(7.7))          # 2a + b = 10.3 rounded once (k = 2)
        need(bits(target) in [bits(x) for x in mandatory_knots(row, 3)] and target == 10.3, "RN(2a+b)")
        name = "initial_panel_edges_eV" if "initial" in key else "panel_edges_eV"
        edges = [bits(x) for x in G(G(order_of(result, 3), "mesh"), name)]
        need(bits(target) in edges, f"RN(2a+b) is not an edge of mesh.{name}")
        if name == "panel_edges_eV":
            panels = by_order(out.panels, 3, 6)
            real = {bits(G(p, "left_eV")) for p in panels} | {bits(G(p, "right_eV")) for p in panels}
            need(bits(target) in real, "RN(2a+b) is not an edge of the returned order-3 panels")
        return None
    if key == "single_loss panels":
        want = [2.0, 4.0, 5.0] if rid == "AE1" else [1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 50.0]
        got = first_order_edges(out)
        return None if got == want else ("not_run", f"precondition: first-order panels {got}, expected {want}")
    if rid == "OI1":
        mesh = G(order_of(result, 2), "mesh")
        if key.endswith("mesh.stages"):
            need(G(mesh, "stages") >= 1, f"stages {G(mesh, 'stages')}")
        elif key.endswith("mesh.history"):
            need([float(x) for x in G(mesh, "initial_panel_edges_eV")] == [4.0, 7.0, 10.0], "initial order-2 mesh")
            history = G(mesh, "history")
            need(len(history) >= 2, "mesh.history does not record the failed initial mesh and a later stage")
            seen = out.extras["order_inject"].seen
            need(any(e == 1e-3 for _, _, e in seen) and any(e == 1e-12 for _, _, e in seen), "order_inject calls")
        else:
            widths = np.diff(np.asarray(G(mesh, "panel_edges_eV"), float))
            need(len(widths) >= 4 and widths.max() <= 1.5, f"final panel widths {widths}")
        return None
    if rid == "OI2":
        need(G(G(order_of(result, 2), "mesh"), "builds") == 1, "mesh.builds")
        return None
    if rid in ("A4b", "A5"):
        n = 2 + index
        need(G(G(order_of(result, n), "mesh"), "stages") >= 1, f"order {n}: no bisection stage")
        return None
    if rid in ("T1a", "T1b", "T1c"):
        name = key.split(".")[1]
        want = row["expect"]["support_order"] if name == "support_order" else row["expect"]["truncation_state"]
        edges = [Fr(x) for x in edge_values(row["inputs"]["output_edges_eV"])]
        exact = max(0, math.floor((Fr(1000) - edges[0])/Fr(1)))
        need(exact == row["expect"]["support_order"], "exact support order of the table")
        need(G(truncation(result), name) == want, f"{name} = {G(truncation(result), name)!r}, expected {want!r}")
        return None
    if rid == "T4u" and key.startswith("order_area_ratio"):
        lo, hi = G(G(result, "order_area_ratio"), "interval")
        need(lo < 1.0 < hi, f"interval {(lo, hi)} does not straddle 1")
        return None
    if rid == "T4u":
        i = int(re.search(r"bin\[(\d+)\]", key).group(1))
        b, reason, tau = bound_triple(result, i)
        need(isinstance(b, float) and math.isfinite(b) and b > 0 and reason is None and tau is not None and tau > 0,
             f"bin[{i}]: {(b, reason, tau)}")
        return None
    if rid == "T4d":
        lo, hi = G(G(result, "order_area_ratio"), "interval")
        need(lo >= 1.0, f"lower end {lo}")
        return None
    raise AssertionError(f"no diagnostic check for {key!r}")


def check_m4a_diagnostic(row, out, key, oracle):
    rid, result = row["id"], out.result
    need(out.status == "success", f"the row did not succeed ({out.status})")
    if key == "on_final_panels.deliveries":
        need(len(out.deliveries) == 1 and out.run_calls == 1, f"{len(out.deliveries)} deliveries, {out.run_calls} _run calls")
        need(G(G(G(result, "diagnostics"), "first_order"), "beta_deliveries") == 1, "beta_deliveries")
        return None
    if key.startswith("order[") and "formal_folds" in key:
        n = int(key[6])
        need(G(G(order_of(result, n), "mesh"), "formal_folds") >= 1, f"order {n}: formal_folds")
        return None
    if key == "on_final_panels.received_panels_are_order2_input":
        need(len(out.deliveries) == 1, "no single delivery")
        received = list(out.deliveries[0][1])
        counts = G(G(result, "diagnostics"), "first_order")
        need(G(counts, "final_panels") == len(received), "final_panels")
        # the same first-order arguments through the analytic port are not available for a material; the
        # hand-over is observed on the higher orders: the order-2 mesh starts from the delivered panels
        left = float(G(received[0], "left_eV"))
        right = float(G(received[-1], "right_eV"))
        knots = [float(x) for x in G(G(order_of(result, 2), "mesh"), "knots_eV")]
        need(knots == sorted({2*left, left + right, 2*right}), "order-2 knots do not come from the delivered panels")
        need(all(np.asarray(G(p, "values")).shape == (11, 17) for p in received), "delivered panels")
        # the delivered objects through layer C (multiplier 1) reproduce the order-2 area of the result
        lo = importlib.import_module("xpsloss._loss_orders")
        features = dict(jumps_eV=(left, right), kinks_eV=(), singularities=(), min_half_width_eV=None,
                        min_variation_scale_eV=None)
        again = lo.convolve_panels(received, (left, right), model=row["inputs"]["model"], max_order=2, features=features)
        m, _ = multiplier_of(result)
        area = m*m*panel_functional(by_order(again, 2, 2), "total", lambda s: np.ones_like(s))
        w, e = read_area(result, 2, "total")
        need(abs(area - w) <= e + 1e-6*abs(w), f"order-2 area from the delivered panels {area!r} against {w!r}")
        return None
    base = baseline()
    if base["unavailable"]:
        return base["unavailable"]
    if key == "_call_kernel.call_sequence":
        snap = base["record"]["rows"][M4A_FIRST.get(rid, rid)]
        d = first_difference(canonical_calls(out.calls), snap["calls"])
        need(d is None, f"_call_kernel calls differ from the baseline: {d}")
        return None
    raise AssertionError(f"no M4A diagnostic check for {key!r}")


def check_cli_diagnostic(row, out, key, index):
    rid, record = row["id"], out.cli
    need(record is not None and record["exit_code"] == 0, "the CLI did not complete")
    if rid == "C11":
        import xpsloss
        caps = json.loads(record["stdout"])
        base = baseline()
        if base["unavailable"]:
            return base["unavailable"]
        old = json.loads(base["record"]["capabilities"]["stdout"])
        if index == 0:
            need(caps["supported"] == old["supported"] + ["xps_multiple_loss"], "supported")
        elif index == 1:
            for name, value in old.items():
                if name not in ("supported", "version"):
                    need(caps.get(name) == value, f"block {name} differs from the baseline")
        elif index == 2:
            need(caps["version"] == xpsloss.__version__ and old["version"] == base["record"]["version"] == "0.3.0",
                 "each side reports its own version")
        else:
            block = caps["multiple_scattering"]
            need(block["model_id"] == table()["model_id"] and block["emission"] == "normal_only"
                 and block["kernel_energy_policy"] == "fixed_E0" and block["status"] == "approximate"
                 and unsupported_ok(block["unsupported_physics"]), "multiple_scattering block")
            need(not any(name in record["stdout"] for name in table()["reference_only_models"]),
                 "a reference-only model id appears in capabilities")
        return None
    manifest = record["manifest"]
    results = G(manifest, "results")
    if rid == "C02":
        t = G(results, "truncation")
        want = (("state", "complete_on_display"), ("computed", 2), ("support_order", 2))[index]
        need(G(t, want[0]) == want[1], f"truncation.{want[0]} = {G(t, want[0])!r}")
        return None
    import xpsloss
    checks = [
        lambda: G(manifest, "model_id") == table()["model_id"],
        lambda: G(manifest, "depends_on_model_id") == xpsloss.SINGLE_LOSS_MODEL_ID,
        lambda: {"depth_averaged", "fixed_depth", "dielectric", "single_loss"} <= set(G(manifest, "dependent_model_ids")),
        lambda: (G(manifest, "repeated_kernel"), G(manifest, "hole_treatment"), G(manifest, "kernel_energy_policy"))
        == ("term_13a", "first_event_only", "fixed_E0"),
        lambda: [G(o, "order") for o in G(manifest, "orders_definition")] == [1, 2, 3],
        lambda: G(G(manifest, "truncation"), "state") == "truncated",
        lambda: len(G(G(manifest, "order_area_ratio"), "interval")) == 2 and "convention_ratio" in manifest,
        lambda: set(G(manifest, "effective_settings")) == set(FIRST_ORDER_NAMES) | set(ORDER_NAMES)
        and set(G(manifest, "effective_settings_sources").values()) <= {"cli", "config", "default"}
        and len(G(manifest, "effective_settings")) == 13,
        lambda: G(G(manifest, "spectrum_definition"), "max_order") == 3,
        lambda: all(G(G(manifest, "comparison"), g) == "not_matched" or
                    G(G(G(manifest, "comparison"), g), "status") == "not_matched" for g in ("orders", "areas", "spectrum"))
        and G(manifest, "physical_validity") == "not_established",
        lambda: G(manifest, "package_version") == xpsloss.__version__ and "runtime" in manifest
        and len(G(manifest, "warnings")) >= len(ALWAYS_WARNINGS),
        lambda: G(manifest, "input_sha256") == record["config_sha256"]
        and record["spectrum_sha256"] in json.dumps(G(manifest, "spectrum_input")),
    ]
    need(index < len(checks), "ledger count of the run.json fields")
    need(checks[index](), f"run.json record field check {index + 1} of 12")
    return None


# ---------------------------------------------------------------- replay (C01)
def parse_csv(raw):
    text = raw.decode("ascii")
    need("\r" not in text, "CSV has carriage returns")
    lines = text.strip("\n").split("\n")
    return lines[0].split(","), [line.split(",") for line in lines[1:]]


def csv_float(text, value, label):
    need(text == "%.17g" % value, f"{label}: CSV {text!r} is not %.17g of {value!r}")
    need(bits(float(text)) == bits(value), f"{label}: CSV value does not read back bit for bit")


def check_replay(row, out, key):
    index, _ = enumerated(key)
    record = out.cli
    need(record is not None and record["exit_code"] == 0, "the CLI did not complete")
    need("api_replay" in out.extras, f"API replay failed: {out.extras.get('api_replay_error', '')[-400:]}")
    replay = out.extras["api_replay"]
    orders = [G(o, "order") for o in G(replay, "orders")]
    if index == 0:
        header, lines = parse_csv(record["files"]["multiple_loss_spectrum.csv"])
        want = ["region", "edge_low_eV", "edge_high_eV", "zero_loss", "zero_loss_error"]
        for n in orders:
            want += [f"order{n}_{c}{s}" for c in COMPONENTS for s in ("", "_error")]
        need(header == want + ["spectrum", "spectrum_error"], f"header {header}")
        edges = edge_values(row["inputs"]["spectrum"]["output_edges_eV"])
        regions = region_names(len(edges) - 1)
        need([line[0] for line in lines] == ["low"] + [str(i) for i in range(len(edges) - 1)] + ["high"], "row names")
        for line, region in zip(lines, regions):
            cells = iter(line[3:])
            values = list(read_zero_loss(replay, region))
            for n in orders:
                for c in COMPONENTS:
                    values += list(read_region(replay, n, region, c))
            values += list(read_spectrum(replay, region))
            for value in values:
                csv_float(next(cells), value, region)
        need(lines[0][1] == "" and lines[-1][2] == "", "the infinite ends are not empty")
        return None
    if index == 1:
        header, lines = parse_csv(record["files"]["order_areas.csv"])
        need(header == ["order", "column", "area", "area_error", "cutoff_sensitivity", "cutoff_sensitivity_error"], "header")
        need([(int(line[0]), line[1]) for line in lines] == [(n, c) for n in orders for c in COLUMNS], "rows")
        for line in lines:
            n, col = int(line[0]), line[1]
            a, ea = read_area(replay, n, col)
            h, eh = read_area(replay, n, col, "cutoff_sensitivity")
            for text, value in zip(line[2:], (a, ea, h, eh)):
                csv_float(text, value, f"order {n} {col}")
        return None
    manifest = record["manifest"]
    results = G(manifest, "results")
    for name in ("spectrum", "spectrum_error", "zero_loss", "zero_loss_error"):
        need([bits(x) for x in G(results, name)] == [bits(x) for x in G(replay, name)], f"results.{name}")
    for o_json, o in zip(G(results, "orders"), G(replay, "orders")):
        for c in COMPONENTS:
            need([bits(x) for x in G(o_json, c)] == [bits(x) for x in G(o, c)], f"results.orders {c}")
            need([bits(x) for x in G(o_json, c + "_error")] == [bits(x) for x in G(o, c + "_error")], f"orders {c}_error")
        need({k: bits(v) for k, v in G(o_json, "area").items()} == {k: bits(v) for k, v in G(o, "area").items()}, "area")
    need(G(G(results, "truncation"), "state") == G(truncation(replay), "state"), "state")
    need(list(G(manifest, "warnings")) == list(G(replay, "warnings")), "warnings")
    need(G(manifest, "model_id") == G(replay, "model_id"), "model id")
    return None


# ---------------------------------------------------------------- recomputed from the report, inequalities
def rho_n(m, interval, n):
    if interval is None:
        return Fr(0)
    m, lo, hi = Fr(m), Fr(interval[0]), Fr(interval[1])
    return max(1 - (lo/m)**n, (hi/m)**n - 1)


def outward(y, exact, label, rel=4*FEPS):
    need(Fr(y) >= exact and Fr(y) <= exact*(1 + rel) + Fr(5, 10**324),
         f"{label}: {y!r} against the rebuilt {float(exact)!r} (not below; above by at most {float(rel):.1e})")


def area_parts(result, n, col):
    box = G(order_of(result, n), "area_error_parts")
    return {p: Fr(float(G(G(box, p), col))) for p in PARTS + ("convention",)}


def local_terms_sum(result, n, col):
    terms = np.asarray(G(G(order_of(result, n), "mesh"), "local_area_terms"), float)
    need(terms.ndim == 2 and terms.shape[1] == 11, f"local_area_terms shape {terms.shape}")
    return sum((Fr(float(x)) for x in terms[:, COLUMNS.index(col)]), Fr(0))


def check_recomputed(row, out, key, oracle):
    rid, result = row["id"], out.result
    need(out.status == "success", f"the row did not succeed ({out.status})")
    index, count = enumerated(key)
    if rid == "DS0":
        single = G(result, "single_loss")
        m, interval = multiplier_of(result)
        m_lo, m_hi = (Fr(m), Fr(m)) if interval is None else (Fr(interval[0]), Fr(interval[1]))
        A = Fr(float(G(G(single, "window_area"), "total")))
        e = Fr(float(G(G(single, "window_area_error"), "total")))
        corners = [x*y for x in (m_lo, m_hi) for y in (A - e, A + e)]
        w_lo, w_hi = min(corners), max(corners)
        r_lo, r_hi = (Fr(float(x)) for x in G(G(result, "order_area_ratio"), "interval"))
        need(r_hi < 1, "the ratio interval reaches 1")
        S, _ = primary_sums(row)
        got = G(truncation(result), "neglected_total_area")
        lo, hi = round_down(S*w_lo/(1 - r_lo)), round_up(S*w_hi/(1 - r_hi))
        need((bits(G(got, "interval")[0]), bits(G(got, "interval")[1])) == (bits(lo), bits(hi)),
             f"interval {G(got, 'interval')} against the rebuilt {(lo, hi)}")
        need(lo <= G(got, "value") <= hi, "value outside its interval")
        return None
    if rid == "AE1":
        m, _ = multiplier_of(result)
        exact = Fr(m)**2*local_terms_sum(result, 2, "hole")
        outward(float(area_parts(result, 2, "hole")["local"]), exact, "order 2 local part of hole", Fr(1, 10**9))
        return None
    need(rid == "A3f", f"no recomputation for {rid}")
    single = G(result, "single_loss")
    m, interval = multiplier_of(result)
    seeds = G(G(G(result, "diagnostics"), "first_order"), "fixed_multiplier_area_error")
    if key.startswith("order[") and ".area_error." in key:
        mm = re.match(r"order\[(\d)\]\.area_error\.(\w+)$", key)
        n, col = int(mm.group(1)), mm.group(2)
        parts = area_parts(result, n, col)
        w, e = read_area(result, n, col)
        e_f = parts["inherited"] + parts["propagated"] + parts["local"]
        rho = rho_n(m, interval, n)
        need(rho > 0, "no multiplier interval")
        outward(e, e_f + rho*(abs(Fr(w)) + e_f), f"order {n} area_error {col}", 8*FEPS)
        outward(float(parts["multiplier_interval"]), rho*(abs(Fr(w)) + e_f), f"order {n} interval part {col}", 8*FEPS)
        need(parts["convention"] == 0, "convention part of a higher order")
        return None
    if key.startswith("order[1] area and error mapping"):
        A = {c: Fr(float(G(G(single, "window_area"), c))) for c in COLUMNS}
        eA = {c: Fr(float(G(G(single, "window_area_error"), c))) for c in COLUMNS}
        if index >= 11:
            got = G(single, "multiplier_interval_A")
            want = round_down(1/(A["total"] + eA["total"])) if index == 11 else round_up(1/(A["total"] - eA["total"]))
            need(bits(got[index - 11]) == bits(want), f"interval end {got[index - 11]!r}, rebuilt {want!r}")
            return None
        col = COLUMNS[index]
        parts = area_parts(result, 1, col)
        w, e = read_area(result, 1, col)
        seed = Fr(float(G(seeds, col)))
        outward(float(seed), Fr(m)*eA[col] + FEPS*abs(Fr(m)*A[col]), f"fixed_multiplier_area_error[{col}]")
        if col == "total":
            need(bits(w) == bits(G(single, "single_loss_weight")) and bits(e) == bits(G(single, "single_loss_weight_error")),
                 "order 1 total is not the legacy weight")
            need(G(G(G(result, "diagnostics"), "first_order"), "area_constraint") == "inverse_window_area", "area_constraint")
            need(parts["convention"] == Fr(e) and all(parts[p] == 0 for p in PARTS), f"S2 parts of total: {parts}")
            return None
        need(bits(w) == bits(float(Fr(m)*A[col])), f"area[{col}] is not m*A rounded once")
        e_f = Fr(m)*eA[col] + FEPS*abs(Fr(w))
        outward(float(parts["inherited"]), e_f, f"order 1 inherited {col}")
        rho = Fr(float(G(G(G(single, "diagnostics"), "multiplier"), "rho")))
        outward(e, parts["inherited"] + rho*(abs(Fr(w)) + parts["inherited"]), f"order 1 area_error {col}")
        need(parts["propagated"] == 0 and parts["local"] == 0 and parts["convention"] == 0, "order 1 parts")
        return None
    mm = re.match(r"order\[3\]\.area_error_parts\[(\w+)\]\[total\]$", key)
    need(mm is not None, f"no recomputation for {key!r}")
    part = mm.group(1)
    norms = ra_references()["rows"]["A3f"]["notes"]["norms_for_order_3"]
    n_r, n_c2 = Fr(norms["norm_R_per_A"]), Fr(norms["norm_C2_total_per_A2"])
    p2, p3 = area_parts(result, 2, "total"), area_parts(result, 3, "total")
    e_r = Fr(float(G(seeds, "term_13a")))
    M = Fr(m)
    if part == "propagated":
        outward(float(p3[part]), M*n_r*(p2["local"] + p2["propagated"]), "order 3 propagated", Fr(1, 10**9))
    elif part == "inherited":
        exact = M*n_r*p2["inherited"] + M*M*n_c2*e_r + p2["inherited"]*e_r + (p2["local"] + p2["propagated"])*e_r
        need(Fr(float(p3[part])) >= exact*(1 - Fr(1, 10**12)) and Fr(float(p3[part])) <= exact*(1 + Fr(1, 10**9)),
             f"order 3 inherited {float(p3[part])!r} against the rebuilt {float(exact)!r}")
    else:
        outward(float(p3[part]), M**3*local_terms_sum(result, 3, "total"), "order 3 local", Fr(1, 10**9))
    return None


def check_inequality(row, out, key):
    rid, result = row["id"], out.result
    need(out.status == "success", f"the row did not succeed ({out.status})")
    if rid == "AE1":
        s = settings_of(row)
        w, _ = read_area(result, 2, "hole")
        budget = max(Fr(s["order_epsabs_fraction"]), Fr(s["order_epsrel"])*abs(Fr(w)))
        need(area_parts(result, 2, "hole")["local"] <= budget, "the local part exceeds the local budget")
        return None
    m = re.match(r"order\[(\d)\]\.bin\[1\]\.total_error_parts\[(\w+)\]", key)
    n, part = int(m.group(1)), m.group(2)
    value = read_region_part(result, n, "bin[1]", "total", part)
    if rid == "OI1":
        y, _ = read_region(result, 2, "bin[1]", "total")
        if part == "local":
            need(value > 0, "the local part is not positive")
        elif part == "propagated":
            need(value == 0.0, f"the propagated part is {value!r}")
        else:
            need(value <= 1e-10*abs(y) + 1e-13, f"the inherited part is {value!r}")
        return None
    if key.endswith("> 0"):
        need(value > 0, f"{part} of order {n} is {value!r}")
    elif key.endswith("= 0"):
        need(value == 0.0, f"{part} of order {n} is {value!r}")
    else:
        other = read_region_part(success_result("A1s"), n, "bin[1]", "total", part)
        need(abs(value - other) <= 1e-6*abs(other) + 1e-18, f"inherited {value!r} against {other!r} of A1s")
    return None


# ---------------------------------------------------------------- the state test: one row
def dispatch(row, out, key, kind, oracle):
    if kind == "state":
        return check_state(row, out)
    if out.status == "harness_error":
        raise AssertionError(f"harness error: {(out.harness_error or '')[-300:]}")
    if kind == "identity":
        return check_identity(row, out, key, oracle)
    if kind == "exact_zero_or_null":
        return check_zero(row, out, key, oracle)
    if kind == "bitwise":
        return check_bitwise(row, out, key, oracle)
    if kind == "diagnostic":
        return check_diagnostic(row, out, key, oracle)
    if kind == "replay":
        return check_replay(row, out, key)
    if kind == "recomputed_from_report":
        return check_recomputed(row, out, key, oracle)
    if kind == "inequality":
        return check_inequality(row, out, key)
    raise AssertionError(kind)


def judge_row(row, judgement):
    """Run the row and judge every planned key of the kinds of the state test."""
    out = run_row(row["id"])
    for key, kind, oracle in own_plan(row["id"]):
        status, detail = guard(lambda: dispatch(row, out, key, kind, oracle))
        judgement.add(row["id"], key, kind, status, detail)
    return out


# ---------------------------------------------------------------- R-A references
def ra_references():
    if "ra" not in _CACHE:
        raw = RA_PATH.read_bytes()
        _CACHE["ra_sha256"] = hashlib.sha256(raw).hexdigest()
        _CACHE["ra"] = json.loads(raw.decode("ascii"))
    return _CACHE["ra"]


def ra_sha256():
    ra_references()
    return _CACHE["ra_sha256"]


def judge_accuracy(y, ref):
    """Prereg 7.2: |y - r| + u_r <= d_abs + d_rel*|r| (the floor alone for r = 0)."""
    if ref["status"] != "reference":
        return "no_reference", ref["no_value_reason"]
    if not math.isfinite(y):
        return "fail", f"non-finite value {y!r}"
    r, u = Fr(ref["value"]), Fr(ref["u_r"])
    width = Fr(ref["width"]["d_abs"]) + Fr(ref["width"]["d_rel"])*abs(r)
    d = abs(Fr(y) - r)
    if d + u <= width:
        return "pass", ""
    return "fail", f"|y - r| = {float(d):.3e} > width {float(width):.3e} (y = {y!r}, r = {float(r)!r})"


def judge_coverage(y, e, ref):
    """Prereg 7.3 (three values): covered, underestimate or unresolved."""
    if ref["status"] != "reference":
        return "no_reference", ""
    r, u = Fr(ref["value"]), Fr(ref["u_r"])
    d = abs(Fr(y) - r)
    allowed = Fr(e) + 4*FEPS*max(abs(Fr(y)), abs(r))
    if d + u <= allowed:
        return "covered", ""
    if d - u > allowed:
        return "underestimate", f"|y - r| = {float(d):.3e} > returned error {e!r}"
    return "unresolved", ""


def gauss(n):
    return np.polynomial.legendre.leggauss(n)


def panel_functional(panels, column, weight, points=24):
    """Sum over the panels of the integral of weight(s)*interpolant(s) with a Gauss rule of the test."""
    j = COLUMNS.index(column)
    x, w = gauss(points)
    k = np.arange(17)
    nodes = -np.cos(np.pi*k/16)
    bary = (-1.0)**k
    bary[[0, 16]] *= 0.5
    total = 0.0
    for p in panels:
        left, right = float(G(p, "left_eV")), float(G(p, "right_eV"))
        values = np.asarray(G(p, "values"), float)[j]
        diff = x[:, None] - nodes[None, :]
        terms = bary[None, :]/diff
        interpolant = (terms @ values)/terms.sum(axis=1)
        s = (left + right)/2 + (right - left)/2*x
        total += (right - left)/2*float(np.sum(w*weight(s)*interpolant))
    return total


def read_reference_value(row, out, key, ref):
    """(y, returned error or None, extra uncertainty of the test) of a direct_reference key."""
    klass, result = ref["class"], out.result
    if klass == "region":
        return (*read_region(result, ref["order"], ref["region"], ref["component"]), 0.0)
    if klass == "area":
        return (*read_area(result, ref["order"], ref["column"]), 0.0)
    if klass == "sensitivity":
        return (*read_area(result, ref["order"], ref["column"], "cutoff_sensitivity"), 0.0)
    if klass == "spectrum":
        return (*read_spectrum(result, ref["region"]), 0.0)
    if klass == "region_difference":
        y, _ = read_region(result, ref["order"], ref["region"], ref["component"])
        z, _ = read_region(success_result("A3a"), ref["order"], ref["region"], ref["component"])
        return float(Fr(y) - Fr(z)), None, 0.0
    if klass in ("panel_laplace", "panel_moment"):
        panels = by_order(out.panels, ref["order"], len(G(result, "orders")))
        if klass == "panel_laplace":
            tau = ref["tau_per_eV"]
            weight = lambda s: np.exp(-tau*s)
        else:
            weight = (lambda s: np.ones_like(s)) if ref["moment"] == 0 else (lambda s: s)
        y = panel_functional(panels, ref["column"], weight, 24)
        return y, None, abs(y - panel_functional(panels, ref["column"], weight, 32))
    if klass == "window_area":
        return float(G(G(G(result, "single_loss"), "window_area"), "total")), None, 0.0
    raise AssertionError(klass)


COVERAGE_CLASSES = ("region", "area", "sensitivity", "spectrum")


def judge_reference_key(row, out, key, ref):
    """(accuracy status, detail, coverage status or None) of one R-A key."""
    klass, kind, rid = ref["class"], ref["kind"], row["id"]
    if row["entry"] == "layer_c":
        need(out.status == "success", f"the row did not succeed ({out.status})")
        panels = by_order(out.panels, 2, 2)
        j = COLUMNS.index(ref["column"])
        ends = [float(np.asarray(G(p, "values"), float)[j, 16]) for p in panels if float(G(p, "right_eV")) == 7.0] + \
               [float(np.asarray(G(p, "values"), float)[j, 0]) for p in panels if float(G(p, "left_eV")) == 7.0]
        need(len(ends) == 2, f"{len(ends)} panel end nodes at 7 eV (a boundary at 7 eV is required)")
        verdicts = [judge_accuracy(y, ref) for y in ends]
        bad = [v for v in verdicts if v[0] != "pass"]
        return (*(bad[0] if bad else ("pass", "")), None)
    need(out.status == "success", f"the row did not succeed ({out.status}: {out.summary().get('reason')})")
    result = out.result
    if kind == "direct_reference" and klass != "window_area_error":
        y, e, extra = read_reference_value(row, out, key, ref)
        shown = dict(ref)
        if extra:
            shown["u_r"] = str(float(Fr(ref["u_r"]) + Fr(extra)))
        status, detail = judge_accuracy(y, shown)
        coverage = None
        if klass in COVERAGE_CLASSES and row["entry"] != "reference_model" and e is not None:
            coverage = judge_coverage(y, e, ref)[0]
        return status, detail, coverage
    r, u = Fr(ref["value"]), Fr(ref["u_r"])
    if klass == "window_area_error":
        e = Fr(float(G(G(G(result, "single_loss"), "window_area_error"), "total")))
        need(r <= e <= r*(1 + Fr(1, 10**9)), f"window_area_error {float(e)!r} outside [{float(r)!r}, (1 + 1e-9) times it]")
        return "pass", "", None
    if klass == "order_area_ratio":
        box = G(result, "order_area_ratio")
        lo, hi = (Fr(float(x)) for x in G(box, "interval"))
        single = G(result, "single_loss")
        m, interval = multiplier_of(result)
        m_lo, m_hi = (Fr(m), Fr(m)) if interval is None else (Fr(interval[0]), Fr(interval[1]))
        A = Fr(float(G(G(single, "window_area"), "term_13a")))
        e = Fr(float(G(G(single, "window_area_error"), "term_13a")))
        corners = [x*y for x in (m_lo, m_hi) for y in (A - e, A + e)]
        need((bits(float(lo)), bits(float(hi))) == (bits(round_down(min(corners))), bits(round_up(max(corners)))),
             "the interval is not the outward-rounded exact corner product of the reported m and A_13a +- e")
        need(bits(G(box, "value")) == bits(float(Fr(m)*A)), "value is not m*A_13a rounded once")
        if lo <= r - u and r + u <= hi:
            return "pass", "", None
        if hi < r - u or lo > r + u:
            return "fail", f"interval [{float(lo)!r}, {float(hi)!r}] does not meet r = {float(r)!r}", None
        return "unresolved", "", None
    if klass == "neglected_total_area":
        lo, hi = (Fr(float(x)) for x in G(G(truncation(result), "neglected_total_area"), "interval"))
        if lo <= r - u and r + u <= hi:
            return "pass", "", None
        if hi < r - u or lo > r + u:
            return "fail", f"interval [{float(lo)!r}, {float(hi)!r}] does not meet r = {float(r)!r}", None
        return "unresolved", "", None
    if klass == "bound_grid_minimum":
        if rid == "T2e" and first_order_edges(out) != ref["panel_edges_eV"]:
            return "not_run", "precondition: first-order panels", None
        b, reason, tau = bound_triple(result, ref["bin"])
        need(isinstance(b, float) and reason is None and tau is not None and tau > 0, f"bound triple {(b, reason, tau)}")
        lower = Fr(b) >= r*(1 - Fr(1, 10**9)) - u
        upper = Fr(b) <= r*(1 + Fr(1, 10**6)) + u
        need(lower and upper, f"bound {b!r} against the grid minimum {float(r)!r}: lower side "
                              f"{'pass' if lower else 'fail'}, upper side {'pass' if upper else 'fail'}")
        return "pass", "", None
    if klass in ("bound_covers_exact", "bound_ratio"):
        b = bound_triple(result, ref["bin"])[0]
        need(isinstance(b, float), f"bound {b!r}")
        if klass == "bound_covers_exact":
            need(Fr(b) >= r, f"bound {b!r} is below the exact neglected quantity {float(r)!r}")
        else:
            need(Fr(b) <= 10*r, f"bound {b!r} is more than 10 times the exact neglected quantity {float(r)!r}")
        return "pass", "", None
    if klass == "bound_model":
        if not ref["reached"]:
            is_complete_bound(result, ref["bin"])
        else:
            b = bound_triple(result, ref["bin"])[0]
            need(isinstance(b, float) and b > 0 and Fr(b) >= r, f"bound {b!r} against the exact {float(r)!r}")
        return "pass", "", None
    if klass == "young_area":
        if rid == "AE1" and first_order_edges(out) != [2.0, 4.0, 5.0]:
            return "not_run", "precondition: first-order panels [2,4], [4,5]", None
        col = ref["column"]
        m, _ = multiplier_of(result)
        seeds = G(G(G(result, "diagnostics"), "first_order"), "fixed_multiplier_area_error")
        e_j, e_r = Fr(float(G(seeds, col))), Fr(float(G(seeds, "term_13a")))
        rebuilt = Fr(m)*Fr(ref["norm_R_per_A"])*e_j + Fr(m)*Fr(ref["norm_column_per_A"])*e_r + e_j*e_r
        y = area_parts(result, 2, col)["inherited"]
        floor = max(rebuilt, r - u) if rid == "AE1" else rebuilt
        need(y >= floor, f"inherited {float(y)!r} is below the Young right-hand side {float(floor)!r}")
        need(y <= rebuilt*(1 + Fr(1, 10**9)) + 8*FEPS*rebuilt, f"inherited {float(y)!r} exceeds the Young right-hand side "
                                                               f"{float(rebuilt)!r} by more than 1e-9")
        need(abs(rebuilt - r) <= Fr(1, 10**6)*r, "the reported seeds are not the registered injection")
        return "pass", "", None
    if klass == "node_part":
        found = node_part_at(a3f_controls(out)["point_path"], ref["at_eV"])
        need(found, "no node at 7 eV was handed to _fold_order")
        for y in found:
            need(abs(Fr(y) - r) <= Fr(ref["width"]["d_abs"]), f"node part {y!r} against {float(r)!r}")
        return "pass", "", None
    raise AssertionError(f"no reference judgement for class {klass}")


def judge_reference_row(row, out):
    """{key: (status, detail)} of the accuracy and {key: status} of the coverage (prereg 7.2 and 7.3)."""
    accuracy, coverage = {}, {}
    for key, ref in ra_references()["rows"][row["id"]]["keys"].items():
        def check():
            return judge_reference_key(row, out, key, ref)
        try:
            status, detail, cover = check()
        except AssertionError as exc:
            status, detail, cover = "fail", str(exc)[:400], None
        except Exception as exc:
            status, detail, cover = "fail", f"{type(exc).__name__}: {str(exc)[:300]}", None
        accuracy[key] = (status, detail)
        if ref["class"] in COVERAGE_CLASSES and row["entry"] not in ("reference_model", "layer_c"):
            coverage[key] = cover if cover is not None else "not_run"
    return accuracy, coverage


def judge_budgets(row, out):
    """Prereg 7.3: on a successful row the inherited and local parts are within their budgets (rebuilt)."""
    if row["entry"] != "analytic_spectrum" or "inject" in row["inputs"] or out.status != "success":
        return None
    result, s = out.result, settings_of(row)
    _, s_abs = primary_sums(row)
    for o in G(result, "orders"):
        n = G(o, "order")
        if n < 2:
            continue
        for region in row_regions(row):
            for c in COMPONENTS:
                y, _ = read_region(result, n, region, c)
                local = read_region_part(result, n, region, c, "local")
                budget = max(Fr(s["order_epsabs_fraction"])*s_abs, Fr(s["order_epsrel"])*abs(Fr(y)))
                if Fr(local) > budget:
                    return "fail", f"order {n} {region} {c}: local part {local!r} above its budget {float(budget)!r}"
    return "pass", ""


if __name__ == "__main__":
    if len(sys.argv) == 4 and sys.argv[1] == "--baseline-child":
        _baseline_child(sys.argv[2], sys.argv[3])
    else:
        raise SystemExit("usage: m4b1_states.py --baseline-child OUT PAYLOAD")
