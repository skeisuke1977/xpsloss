"""Shared harness of the M4b-2 tests of R-C (prereg docs/M4B2_PREREG.md sections 4, 5, 7).

The implementation is reached only through the surface of prereg section 4: the public names
of ``xpsloss`` (``BindingEnergyAxis``, ``EnergyAxisReference``, ``BindingDeltaPrimary``,
``BindingBinnedPrimary``, ``multiple_loss_spectrum(..., energy_axis=)``), the CLI, the test port
``xpsloss._loss_orders.integrate_orders(..., energy_axis=)``, and the read-only wrap of
``_order_bounds.truncation`` of prereg 4.5 (``xpsloss.xps_single_loss._call_kernel`` is wrapped
to count the kernel calls, never replaced).  A key that cannot be read is a failure, never a skip.

Every ledger record ``(row, key, kind, oracle)`` of the frozen table gets exactly one judgement:
``pass``, ``fail``, ``not_run`` (with a reason, e.g. ``sdist_without_git_history``) or
``no_reference``.  A key of the reference that is not used, or used twice, is a fault of its row.

Conventions fixed here where the registration is silent (reported in rc_report.md):

* the measured-spectrum file of a CLI row carries only the required keys; ``primary`` has
  ``kind``, the positions (``binding_energies_eV`` or ``binding_edges_eV``) and ``intensities``;
* the CLI configuration carries ``loss_grid_eV`` {1, 100, 1} as in the M4b-1 harness;
* ``diagnostics`` is compared with the KE control key by key over the keys of the control,
  without ``warnings`` and without the tolerance entry of ``effective_settings``;
* this module run as ``python tests/m4b2_states.py --baseline-child OUT PAYLOAD`` is the
  baseline child of the BASE rows; it uses only names that exist in the baseline commit.
"""
from contextlib import redirect_stderr, redirect_stdout
import copy
import dataclasses
import datetime
from fractions import Fraction as Fr
import hashlib
import importlib
import io
import json
import math
import os
from pathlib import Path
import re
import struct
import subprocess
import sys
import tempfile
import time
import traceback
import types
import zipfile

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DIRECTORY = ROOT/"benchmarks/numerics/m4b2"
TABLE_PATH = DIRECTORY/"expected_states_m4b2.json"
TABLE_PY = DIRECTORY/"expected_states_m4b2.py"
RC_PATH = DIRECTORY/"references_m4b2.json"
RC_PY = DIRECTORY/"generate_references_m4b2.py"
TABLE_SHA256 = "9c97907dab69a00e110414d6862d0d061bc80212000b8f3180f905bcb28108ea"
TABLE_PY_SHA256 = "58fd1c3d35e2f6e205ae52039cecaf4a0c22ca23a3fe84ea2bde8d533a024d2f"
BASELINE_COMMIT = "8c33d2939dab866decc37fdf9663f5ac2e29b78d"
NO_HISTORY = "sdist_without_git_history"
FEPS = Fr(1, 2**52)
COLUMNS = ("term_13a", "term_13b", "term_14", "term_15_T1", "term_15_T2", "term_15_T3",
           "term_16_e", "term_16_h", "no_hole", "hole", "total")
NO_HOLE = ("term_13a", "term_13b", "term_15_T1", "term_15_T3", "term_16_e")
HOLE = ("term_14", "term_15_T2", "term_16_h")
COMPONENTS = ("total", "no_hole", "hole")
PARTS = ("inherited", "propagated", "local", "multiplier_interval")
BE_CODES = ("binding_axis_assumed_calibrated", "energy_axis_field_not_recorded",
            "recorded_emission_angle_differs_from_theta", "axis_conversion_support_disagreement")
SETTING_NAMES = ("epsrel", "area_epsabs", "spectrum_epsabs_fraction", "panel_limit", "max_cells", "kernel_epsabs",
                 "kernel_epsrel", "kernel_limit", "order_epsrel", "order_epsabs_fraction", "order_limit",
                 "max_convolution_points", "neglect_tolerance_fraction")
NINE_UNSUPPORTED = ["oblique_emission", "trajectory_dependent_weight_length", "self_consistent_length",
                    "intrinsic_loss_separation", "losses_outside_window", "elastic_scattering",
                    "instrument_broadening", "background_subtraction", "energy_dependent_kernel"]
MEASURED = "xpsloss-measured-spectrum-v0"
FILES = {"multiple_loss_spectrum.csv", "order_areas.csv", "run.json"}
_CACHE = {}

# Post-freeze changes of the registration (prereg 11): none.  The frozen table keeps its bytes.
POST_FREEZE_EXPECT = {}


class Unreadable(AssertionError):
    """A field of the surface of prereg 4 is missing or has an unexpected form."""


class KernelContractError(AssertionError):
    """The analytic kernel was asked for an undefined value."""


class NotRun(Exception):
    """A check whose premise does not hold (prereg 7.3): counted apart, never as a pass."""

    def __init__(self, reason, status="not_run"):
        super().__init__(reason)
        self.reason, self.status = reason, status


def need(condition, message):
    if not condition:
        raise AssertionError(message)


def G(obj, name):
    try:
        return obj[name] if isinstance(obj, (dict, list, tuple)) else getattr(obj, name)
    except (KeyError, AttributeError, TypeError, IndexError) as exc:
        raise Unreadable(f"missing field {name!r} ({type(exc).__name__})") from None


def bits(x):
    """Comparable form of a scalar: floats by their bits, None and strings as they are."""
    if x is None or isinstance(x, str):
        return x
    if isinstance(x, (bool, np.bool_)):
        return bool(x)
    if isinstance(x, (int, np.integer)):
        return int(x)
    if isinstance(x, (float, np.floating)):
        return ["float", float(x).hex()]
    raise Unreadable(f"not a scalar: {type(x).__name__}")


def canon(x):
    """Bit-exact, JSON-serialisable form of a result (arrays by dtype, shape and bytes; floats by hex)."""
    if isinstance(x, np.ndarray):
        return ["ndarray", x.dtype.str, list(x.shape), hashlib.sha256(np.ascontiguousarray(x).tobytes()).hexdigest()]
    if x is None or isinstance(x, str):
        return x
    if isinstance(x, (bool, np.bool_)):
        return bool(x)
    if isinstance(x, (float, np.floating)):
        return ["float", float(x).hex()]
    if isinstance(x, (int, np.integer)):
        return int(x)
    if dataclasses.is_dataclass(x) and not isinstance(x, type):
        return {"@type": type(x).__name__, **{f.name: canon(getattr(x, f.name)) for f in dataclasses.fields(x)}}
    if isinstance(x, dict):
        return {str(k): canon(v) for k, v in x.items()}
    if isinstance(x, tuple):
        return ["tuple", *[canon(v) for v in x]]
    if isinstance(x, list):
        return [canon(v) for v in x]
    if isinstance(x, (set, frozenset)):
        return ["set", *sorted(repr(v) for v in x)]
    if isinstance(x, BaseException):
        return {"@exception": type(x).__name__, "message": str(x), "attributes": canon(dict(vars(x)))}
    if isinstance(x, Fr):
        return ["fraction", str(x)]
    if callable(x):
        return ["callable", getattr(x, "__qualname__", type(x).__name__)]
    if hasattr(x, "__dict__"):
        return {"@object": type(x).__name__, **{str(k): canon(v) for k, v in vars(x).items()}}
    return ["repr", type(x).__name__, repr(x)]


def plain(x):
    """canon() passed through JSON so that a snapshot of the parent equals one read from the child."""
    return json.loads(json.dumps(canon(x)))


def first_difference(a, b, path=""):
    if type(a) is not type(b):
        return f"{path}: type {type(a).__name__} != {type(b).__name__}"
    if isinstance(a, dict):
        if set(a) != set(b):
            return f"{path}: keys {sorted(set(a) ^ set(b))}"
        if list(a) != list(b):
            return f"{path}: key order differs"
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


def same(a, b, label):
    d = first_difference(plain(a), plain(b), label)
    need(d is None, d)


# ---------------------------------------------------------------- table, ledger, reference
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
        scope = {"__file__": str(TABLE_PY), "__name__": "m4b2_expected_states_table"}
        saved = sys.argv
        sys.argv = [str(TABLE_PY), "--check"]      # --check never writes
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


def references():
    if "rc" not in _CACHE:
        raw = RC_PATH.read_bytes()
        _CACHE["rc_sha256"] = hashlib.sha256(raw).hexdigest()
        _CACHE["rc"] = json.loads(raw.decode("ascii"))
    return _CACHE["rc"]


def rc_sha256():
    references()
    return _CACHE["rc_sha256"]


def rc_row(rid):
    try:
        return references()["rows"][rid]
    except KeyError:
        raise NotRun(f"row {rid} is not in the reference", "no_reference") from None


def rc_ledger_keys(rid):
    return [key for key, _kind, oracle in ledger()[rid] if oracle.startswith("R-C")]


class RefUse:
    """Counts how often every reference key of a row is read (prereg 7.1: unused and duplicate use fail)."""

    def __init__(self, rid):
        self.rid, self.used = rid, {}

    def get(self, key):
        keys = rc_row(self.rid)["keys"]
        if key not in keys:
            raise NotRun(f"reference key {key!r} is missing", "no_reference")
        self.used[key] = self.used.get(key, 0) + 1
        return keys[key]

    def faults(self):
        if self.rid not in references()["rows"]:
            return {}
        keys = references()["rows"][self.rid]["keys"]
        out = {}
        unused = sorted(k for k in keys if k not in self.used)
        twice = sorted(k for k, n in self.used.items() if n > 1)
        if unused:
            out["unused"] = unused
        if twice:
            out["duplicate"] = twice
        return out


def full_run():
    return os.environ.get("XPSLOSS_M4B2_FULL") == "1"


def selected_ids():
    """Rows of this run in the registered execution order (CI: 70 rows; XPSLOSS_M4B2_FULL=1: all 82)."""
    only = os.environ.get("XPSLOSS_M4B2_ROWS")
    out = []
    for rid in table()["execution_order"]:
        if only:
            if rid in only.split(","):
                out.append(rid)
        elif rows()[rid]["ci"] or full_run():
            out.append(rid)
    return out


def planned_counts(ids):
    counts = {}
    for rid in ids:
        for _key, kind, _oracle in ledger()[rid]:
            counts[kind] = counts.get(kind, 0) + 1
    return counts


# ---------------------------------------------------------------- numbers of the table and of the reference
def num(obj):
    if isinstance(obj, dict):
        if "float" in obj:
            return float(obj["float"])
        value = float.fromhex(obj["hex"])
        if "value" in obj and value != float(obj["value"]):
            raise AssertionError(obj)
        return value
    return obj


def hexes(obj):
    return [float.fromhex(h) for h in obj["hex"]]


def rv(obj):
    """A value of the reference file: {'hex': ...} is a binary64, everything else is itself."""
    if isinstance(obj, dict) and set(obj) == {"hex"}:
        return float.fromhex(obj["hex"])
    return obj


def within_membership(y, r):
    """Prereg 7.2, bin_membership: reference 0 -> exactly 0; else |y - r| <= 4 eps |r| (exact arithmetic)."""
    y = float(y)
    if not math.isfinite(y):
        return False
    return y == 0.0 if r == 0 else abs(Fr(y) - r) <= 4*FEPS*abs(r)


# ---------------------------------------------------------------- inputs
def material_spec(name):
    m = table()["materials"][name]
    return dict(name=name, gap_hex=m["gap_hex"], oscillators_hex=m["oscillators_hex"], gap_eV=m["gap_eV"],
                oscillators=m["oscillators"])


def make_material(xp, spec):
    return xp.Material(spec["name"], float.fromhex(spec["gap_hex"]),
                       tuple(xp.Oscillator(*(float.fromhex(h) for h in o)) for o in spec["oscillators_hex"]))


def material_dict(spec):
    return dict(name=spec["name"], gap_eV=spec["gap_eV"], oscillators=[
        dict(zip(("energy_eV", "strength_eV2", "width_eV", "dispersion"), o)) for o in spec["oscillators"]])


def ke_primary_arrays(name):
    p = table()["kinetic_primaries"][name]
    key = "energies_eV" if p["kind"] == "delta" else "edges_eV"
    return p["kind"], hexes(p[key]), hexes(p["intensities"])


def be_primary_arrays(name):
    p = table()["binding_primaries"][name]
    key = "binding_energies_eV" if p["kind"] == "delta" else "binding_edges_eV"
    return p["kind"], hexes(p[key]), hexes(p["intensities"])


def make_ke_primary(xp, kind, x, s):
    return (xp.DeltaPrimary if kind == "delta" else xp.BinnedPrimary)(np.array(x, float), np.array(s, float))


def make_be_primary(xp, kind, x, s):
    cls = G(xp, "BindingDeltaPrimary" if kind == "delta" else "BindingBinnedPrimary")
    return cls(np.array(x, float), np.array(s, float))


def be_edges(spec):
    """BE output edges of a row: (values in input order, the JSON form for a CLI file)."""
    if isinstance(spec, dict):
        registered = table()["binding_output_edges"][spec["compact_form_of"]]
        return hexes(registered), dict(registered["compact"])
    values = hexes(table()["binding_output_edges"][spec])
    return values, values


def axis_plain(spec):
    """The axis of a row as nested plain data (prereg 3), before the defects of the row."""
    axis = table()["axes"][spec["axis"]]
    ref = {k: (num(v) if isinstance(v, dict) else v) for k, v in axis["reference"].items()}
    out = dict(photon_energy_eV=num(axis["photon_energy_eV"]), reference=ref,
               primary_reference_binding_eV=num(spec["primary_reference_binding_eV"]),
               primary_reference_choice=spec["primary_reference_choice"])
    out.update(copy.deepcopy(table()["records"][spec["record"]]))
    return out


def apply_defects(obj, defects, allow_delete=False):
    """The defect operations of the table on plain data; returns True for pass_plain_dict."""
    as_dict = False
    for d in defects:
        if d["op"] == "pass_plain_dict":
            as_dict = True
            continue
        parts = d["path"].split(".")
        target = obj
        for p in parts[:-1]:
            target = target[p]
        if d["op"] == "delete_key":
            need(allow_delete, "delete_key is a CLI operation")
            del target[parts[-1]]
        else:
            need(d["op"] == "set", d["op"])
            target[parts[-1]] = num(d["value"]) if isinstance(d["value"], dict) and "float" in d["value"] else d["value"]
    return as_dict


def make_axis(xp, plain_axis, as_dict=False):
    if as_dict:
        return copy.deepcopy(plain_axis)
    data = copy.deepcopy(plain_axis)
    reference = G(xp, "EnergyAxisReference")(**data.pop("reference"))
    return G(xp, "BindingEnergyAxis")(reference=reference, **data)


def not_recorded_names(plain_axis):
    """The list of prereg 3: registered order, a parent that is 'not_recorded' stands for its children."""
    out = []
    for name in table()["not_recorded_order"]:
        parts = name.split(".")
        node, found = plain_axis, True
        for depth, p in enumerate(parts):
            if not isinstance(node, dict) or p not in node:
                found = False
                break
            node = node[p]
        if found and isinstance(node, str) and node == "not_recorded":
            out.append(name)
    return tuple(out)


def record_fields(plain_axis):
    return {name: plain_axis[name] for name in table()["record_fields"]}


# ---------------------------------------------------------------- analytic kernels
_SIGN = {"left": -1, "right": 1, "interior": 0}


class AnalyticKernel:
    """kernel(T, side) -> (values[11, n], errors[11, n]) of a box kernel of the table; counts its calls."""

    def __init__(self, name, charge=1.0):
        self.spec = table()["kernels"][name]
        self.charge = float(charge)
        self.calls = self.points = 0

    def __call__(self, T, side):
        self.calls += 1
        t = np.array(T, float, ndmin=1)
        self.points += t.size
        sides = [side]*t.size if isinstance(side, str) else [str(s) for s in np.array(side, ndmin=1)]
        sign = np.array([_SIGN[s] for s in sides])
        values = np.zeros((11, t.size))
        for name, pieces in self.spec["columns"].items():
            for p in pieces:
                need(p["type"] == "box", p["type"])
                lo, hi = (num(x) for x in p["interval_eV"])
                if np.any(((t == lo) | (t == hi)) & (sign == 0)):
                    raise KernelContractError(f"side 'interior' at a jump of the kernel ({lo}, {hi})")
                inside = ((t > lo) & (t < hi)) | ((t == lo) & (sign > 0)) | ((t == hi) & (sign < 0))
                values[COLUMNS.index(name)] += np.where(inside, num(p["value"]), 0.0)
        values[8] = sum(values[COLUMNS.index(n)] for n in NO_HOLE)
        values[9] = sum(values[COLUMNS.index(n)] for n in HOLE)
        values[10] = values[8] + self.charge*values[9]
        return values, np.zeros_like(values)


def features_of(name):
    f = table()["kernels"][name]["features"]
    return dict(jumps_eV=tuple(float(x) for x in f.get("jumps_eV", ())),
                kinks_eV=tuple(float(x) for x in f.get("kinks_eV", ())),
                singularities=tuple((float(c), float(g)) for c, g in f.get("singularities", ())),
                min_half_width_eV=float(f["min_half_width_eV"]) if "min_half_width_eV" in f else None,
                min_variation_scale_eV=float(f["min_variation_scale_eV"]) if "min_variation_scale_eV" in f else None)


def kernel_window(name):
    return tuple(num(x) for x in table()["kernels"][name]["window_eV"])


# ---------------------------------------------------------------- outcome of a row
class Outcome:
    def __init__(self, row_id):
        self.row_id, self.status, self.exc = row_id, None, None
        self.result = self.panels = self.kernel = self.cli = None
        self.calls, self.truncation_calls = [], []
        self.extras, self.seconds, self.harness_error = {}, 0.0, None


def classify(exc):
    from xpsloss import NumericalConvergenceError
    return "not_converged" if isinstance(exc, NumericalConvergenceError) else "error"


class patched:
    """Record the calls of _call_kernel (wrapped, never replaced)."""

    def __init__(self, out):
        self.out = out

    def __enter__(self):
        from xpsloss import xps_single_loss as xs
        self.xs, self.saved = xs, xs._call_kernel
        original, out = xs._call_kernel, self.out

        def call_kernel(material_, electron_energy_eV, loss_eV, **kw):
            out.calls.append(float(electron_energy_eV))
            return original(material_, electron_energy_eV, loss_eV, **kw)
        xs._call_kernel = call_kernel
        return self

    def __exit__(self, *exc):
        self.xs._call_kernel = self.saved
        return False


class watched_truncation:
    """Read-only wrap of the function _loss_orders calls for the truncation record (prereg 4.5)."""

    def __init__(self, out):
        self.out = out

    def __enter__(self):
        self.ob = importlib.import_module("xpsloss._order_bounds")
        self.lo = importlib.import_module("xpsloss._loss_orders")
        original, out = self.ob.truncation, self.out
        self.saved = (original, getattr(self.lo, "truncation", None))

        def truncation(*args, **kw):
            out.truncation_calls.append((args, dict(kw)))
            return original(*args, **kw)
        self.ob.truncation = truncation
        if self.saved[1] is not None:
            need(self.saved[1] is original, "_loss_orders.truncation is not _order_bounds.truncation")
            self.lo.truncation = truncation
        self.original = original
        return self

    def __exit__(self, *exc):
        self.ob.truncation = self.saved[0]
        if self.saved[1] is not None:
            self.lo.truncation = self.saved[1]
        return False


# ---------------------------------------------------------------- entries
def analytic_call(row, *, control=None, tolerance="row"):
    """(kernel, window, kwargs) of an analytic row; ``control`` replaces the BE inputs by the KE ones."""
    import xpsloss
    i = row["inputs"]
    kw = dict(features=features_of(i["kernel"]), weight_length_A=float(i["weight_length_A"]),
              hole_charge=float(i["hole_charge"]), model=i["model"], max_order=i["max_order"])
    kw.update({k: num(v) for k, v in row["settings"].items()})
    extras = {}
    if control is not None:
        c = control
        kw["primary"] = make_ke_primary(xpsloss, c["kinetic_primary"]["kind"],
                                        [rv(x) for x in c["kinetic_primary"]["coordinates"]],
                                        [rv(x) for x in c["kinetic_primary"]["intensities"]])
        kw["output_edges_eV"] = np.array([rv(x) for x in c["kinetic_output_edges_eV"]], float)
        kw["max_order"] = c["max_order"]
        if tolerance != "row":
            kw.pop("neglect_tolerance_fraction", None)
    else:
        spec = i.get("energy_axis")
        if isinstance(i["primary"], dict):
            primary = make_ke_primary(xpsloss, *ke_primary_arrays(i["primary"]["kinetic_primary"]))
        elif i["primary"] in table()["kinetic_primaries"]:
            primary = make_ke_primary(xpsloss, *ke_primary_arrays(i["primary"]))
        else:
            primary = make_be_primary(xpsloss, *be_primary_arrays(i["primary"]))
        if i["output_edges_eV"] in table()["kinetic_output_edges"]:
            edges = hexes(table()["kinetic_output_edges"][i["output_edges_eV"]])
        else:
            edges = be_edges(i["output_edges_eV"])[0]
        kw["primary"], kw["output_edges_eV"] = primary, np.array(edges, float)
        if spec is not None:
            plain_axis = axis_plain(spec)
            as_dict = apply_defects(plain_axis, spec.get("defects", ()))
            kw["energy_axis"] = make_axis(xpsloss, plain_axis, as_dict)
            extras.update(axis_plain=plain_axis, axis_object=kw["energy_axis"])
        extras.update(primary=primary, edges=np.array(edges, float))
    return AnalyticKernel(i["kernel"], float(i["hole_charge"])), kernel_window(i["kernel"]), kw, extras


def run_analytic(row, out):
    lo = importlib.import_module("xpsloss._loss_orders")
    kernel, window, kw, extras = analytic_call(row)
    out.kernel = kernel
    out.extras.update(extras)
    out.result, out.panels = lo.integrate_orders(kernel, window, **kw)


def api_call(row, *, control=None, tolerance="row"):
    import xpsloss
    i = row["inputs"]
    kw = dict(loss_window_eV=tuple(num(x) for x in i["loss_window_eV"]), weight_length_A=float(i["weight_length_A"]),
              epsilon_momentum=i["epsilon_momentum"])
    if "multiplier" in i:
        kw["multiplier"] = i["multiplier"]
    kw.update({k: num(v) for k, v in row["settings"].items()})
    if row["entry"] in ("api_spectrum", "api_areas"):
        kw.update(model=i["model"], max_order=i["max_order"])
    mat = make_material(xpsloss, material_spec(i["material"]))
    energy = float(i["electron_energy_eV"])
    extras = {}
    if row["entry"] in ("api_areas", "api_single_area"):
        return mat, energy, None, None, kw, extras
    if control is not None:
        c = control
        primary = make_ke_primary(xpsloss, c["kinetic_primary"]["kind"],
                                  [rv(x) for x in c["kinetic_primary"]["coordinates"]],
                                  [rv(x) for x in c["kinetic_primary"]["intensities"]])
        edges = [rv(x) for x in c["kinetic_output_edges_eV"]]
        kw["max_order"] = c["max_order"]
        if tolerance != "row":
            kw.pop("neglect_tolerance_fraction", None)
    else:
        if i["primary"] in table()["kinetic_primaries"]:
            primary = make_ke_primary(xpsloss, *ke_primary_arrays(i["primary"]))
        else:
            primary = make_be_primary(xpsloss, *be_primary_arrays(i["primary"]))
        if i["output_edges_eV"] in table()["kinetic_output_edges"]:
            edges = hexes(table()["kinetic_output_edges"][i["output_edges_eV"]])
        else:
            edges = be_edges(i["output_edges_eV"])[0]
        spec = i.get("energy_axis")
        if spec is not None:
            plain_axis = axis_plain(spec)
            as_dict = apply_defects(plain_axis, spec.get("defects", ()))
            kw["energy_axis"] = make_axis(xpsloss, plain_axis, as_dict)
            extras.update(axis_plain=plain_axis, axis_object=kw["energy_axis"])
        extras.update(primary=primary, edges=np.array(edges, float))
    return mat, energy, primary, np.array(edges, float), kw, extras


def call_api(entry, mat, energy, primary, edges, kw):
    import xpsloss
    if entry == "api_spectrum":
        return xpsloss.multiple_loss_spectrum(mat, energy, primary, edges, **kw)
    if entry == "api_areas":
        return xpsloss.multiple_loss_order_areas(mat, energy, **kw)
    if entry == "api_single_spectrum":
        return xpsloss.single_loss_spectrum(mat, energy, primary, edges, **kw)
    need(entry == "api_single_area", entry)
    return xpsloss.loss_window_area(mat, energy, **kw)


def run_api(row, out):
    mat, energy, primary, edges, kw, extras = api_call(row)
    out.extras.update(extras)
    out.extras["api"] = (mat, energy, primary, edges, kw)
    with patched(out):
        out.result = call_api(row["entry"], mat, energy, primary, edges, kw)


def spectrum_json(row):
    """The --spectrum file of a CLI row as plain data, after the defects of the row."""
    s = row["inputs"]["spectrum"]
    if s["schema"] != MEASURED:
        kind, x, strength = ke_primary_arrays(s["primary"])
        primary = {"kind": kind, "energies_eV" if kind == "delta" else "edges_eV": x, "intensities": strength,
                   "intensity_unit": "arb. unit", "description": None, "included_broadening": None,
                   "truncation": None, "intrinsic_satellites": None}
        return {"schema": s["schema"], "axis": s["axis"], "primary": primary,
                "output_edges_eV": hexes(table()["kinetic_output_edges"][s["output_edges_eV"]])}
    kind, x, strength = be_primary_arrays(s["primary"])
    primary = {"kind": kind, "binding_energies_eV" if kind == "delta" else "binding_edges_eV": x,
               "intensities": strength}
    data = {"schema": s["schema"], "axis": s["axis"], "measurement": axis_plain(s["measurement"]),
            "primary": primary, "output_binding_edges_eV": be_edges(s["output_binding_edges_eV"])[1]}
    apply_defects(data, s.get("defects", ()), allow_delete=True)
    return data


def cli_files(row, folder, spec=None):
    i = row["inputs"]
    cfg = folder/"input.json"
    cfg.write_bytes(json.dumps(dict(material=material_dict(spec or material_spec(i["config_material"])),
                                    electron_energy_eV=i["config_electron_energy_eV"],
                                    loss_grid_eV={"min": 1, "max": 100, "step": 1},
                                    integration=i["config_integration"])).encode("ascii"))
    spectrum = folder/"spectrum.json"
    data = spectrum_json(row)
    spectrum.write_bytes(json.dumps(data).encode("ascii"))
    out = folder/"output"
    argv = [i["command"], str(cfg), "--spectrum", str(spectrum), *i["argv"], "--out", str(out)]
    return argv, cfg, spectrum, out, data


def cli_run(cli, argv, folder):
    stdout, stderr = io.StringIO(), io.StringIO()
    with redirect_stdout(stdout), redirect_stderr(stderr):
        try:
            code = cli.main(argv)
        except SystemExit as exc:
            code = exc.code
    files = {p.name: p.read_bytes() for p in folder.iterdir()} if folder.exists() else None
    return dict(exit_code=0 if code is None else code, stdout=stdout.getvalue(), stderr=stderr.getvalue(),
                output_dir_after="absent" if files is None else "created", files=files or {})


def run_cli(row, out):
    from xpsloss import cli
    with tempfile.TemporaryDirectory() as tmp:
        if row["entry"] == "cli_capabilities":
            with patched(out):
                record = cli_run(cli, [row["inputs"]["command"], *row["inputs"]["argv"]], Path(tmp)/"none")
        else:
            argv, cfg, spectrum, folder, data = cli_files(row, Path(tmp))
            with patched(out):
                record = cli_run(cli, argv, folder)
            record.update(config_path=str(cfg), spectrum_path=str(spectrum), spectrum_data=data,
                          config_sha256=hashlib.sha256(cfg.read_bytes()).hexdigest(),
                          spectrum_sha256=hashlib.sha256(spectrum.read_bytes()).hexdigest())
    out.cli = record
    try:
        record["stderr_json"] = json.loads(record["stderr"]) if record["stderr"].strip() else None
    except ValueError:
        record["stderr_json"] = "not JSON"
    if record["exit_code"] == 0:
        out.status = "complete"
        if "run.json" in record["files"]:
            record["manifest"] = json.loads(record["files"]["run.json"].decode("utf-8"))
            if row["inputs"]["spectrum"] and row["inputs"]["spectrum"]["schema"] == MEASURED:
                try:
                    replay_out = Outcome(row["id"] + ".replay")
                    with patched(replay_out):
                        out.extras["api_replay"] = cli_replay(record["manifest"])
                except Exception:
                    out.extras["api_replay_error"] = traceback.format_exc()
    else:
        payload = record["stderr_json"] if isinstance(record["stderr_json"], dict) else {}
        out.status = payload.get("status", "error")


def cli_replay(manifest):
    """The API call rebuilt from run.json alone (prereg 4.6)."""
    import xpsloss
    d = G(manifest, "spectrum_definition")
    mat = G(G(manifest, "config"), "material")
    m = xpsloss.Material(mat["name"], mat["gap_eV"], tuple(xpsloss.Oscillator(**o) for o in mat["oscillators"]))
    ea, meas = G(manifest, "energy_axis"), G(manifest, "measurement")
    # the fields that energy_axis and measurement share are compared before use
    need(bits(G(ea, "photon_energy_eV")) == bits(G(meas, "photon_energy_eV")), "photon energy of the two records")
    ref = G(meas, "reference")
    need(G(ea, "reference_kind") == G(ref, "kind"), "reference kind of the two records")
    for name, target in (("kinetic_energy_offset_eV", "kinetic_energy_offset_eV"), ("binding_anchor_eV", "binding_anchor_eV"),
                         ("kinetic_anchor_eV", "kinetic_anchor_eV"), ("line", "anchor_line"), ("citation", "anchor_citation")):
        need(bits(ref.get(name)) == bits(G(ea, target)), f"reference.{name} of the two records")
    for name in ("primary_reference_binding_eV", "primary_reference_choice"):
        need(bits(G(ea, name)) == bits(G(meas, name)), f"{name} of the two records")
    need(G(ea, "record") == {k: G(meas, k) for k in table()["record_fields"]}, "the six record fields of the two records")
    axis = make_axis(xpsloss, meas)
    p = G(ea, "input_primary")
    need(set(p) <= {"kind", "binding_energies_eV", "binding_edges_eV", "intensities", "intensity_unit", "description",
                    "included_broadening", "truncation", "intrinsic_satellites"}, f"input_primary keys {sorted(p)}")
    optional = {k: p[k] for k in ("intensity_unit", "description", "included_broadening", "truncation",
                                  "intrinsic_satellites") if k in p}
    if G(p, "kind") == "delta":
        primary = xpsloss.BindingDeltaPrimary(np.array(p["binding_energies_eV"]), np.array(p["intensities"]), **optional)
    else:
        primary = xpsloss.BindingBinnedPrimary(np.array(p["binding_edges_eV"]), np.array(p["intensities"]), **optional)
    weight, multiplier = G(d, "weight"), G(d, "multiplier")
    need(G(weight, "length_kind") == "external_constant" and G(multiplier, "convention") == "weight_length",
         "the replay of this harness covers the external constant length with the weight-length multiplier")
    kw = dict(loss_window_eV=tuple(d["loss_window_eV"]), weight_length_A=weight["length_A"],
              epsilon_momentum=d["epsilon_momentum"], hole_charge=d["hole_charge"], theta_deg=d["theta_deg"],
              model=G(manifest, "model_id"), max_order=d["max_order"], **G(manifest, "effective_settings"))
    return xpsloss.multiple_loss_spectrum(m, d["kernel_energy_eV"], primary,
                                          np.array(G(ea, "input_binding_edges_eV"), float), energy_axis=axis, **kw)


RUNNERS = {"analytic_spectrum": run_analytic, "api_spectrum": run_api, "api_areas": run_api,
           "api_single_spectrum": run_api, "api_single_area": run_api, "cli": run_cli, "cli_capabilities": run_cli}


def execute(out, call):
    start = time.perf_counter()
    try:
        call()
        if out.status is None:
            out.status = "success"
    except (ImportError, AttributeError, TypeError, AssertionError) as exc:
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
    return out


def run_row(rid):
    """Run a row once (cached: a later row that reuses row.result reads the same outcome)."""
    runs = _CACHE.setdefault("runs", {})
    if rid not in runs:
        row = rows()[rid]
        out = Outcome(rid)
        runs[rid] = execute(out, lambda: RUNNERS[row["entry"]](row, out))
    return runs[rid]


def run_control(rid, tolerance=None):
    """The physical control of prereg 7.1: the same entry on the KE coordinates of R-C.

    ``tolerance=None`` removes the tolerance; ``tolerance='row'`` keeps the tolerance of the row."""
    runs = _CACHE.setdefault("controls", {})
    if (rid, tolerance) not in runs:
        row = rows()[rid]
        out = Outcome(f"{rid}.control")
        control = rc_row(rid)["control"]

        def call():
            if row["entry"] == "analytic_spectrum":
                lo = importlib.import_module("xpsloss._loss_orders")
                kernel, window, kw, _ = analytic_call(row, control=control, tolerance=tolerance)
                out.kernel = kernel
                with watched_truncation(out):
                    out.result, out.panels = lo.integrate_orders(kernel, window, **kw)
            else:
                mat, energy, primary, edges, kw, _ = api_call(row, control=control, tolerance=tolerance)
                with patched(out), watched_truncation(out):
                    out.result = call_api(row["entry"], mat, energy, primary, edges, kw)
        runs[(rid, tolerance)] = execute(out, call)
    return runs[(rid, tolerance)]


def success(out):
    need(out.status in ("success", "complete") and (out.result is not None or out.cli is not None),
         f"upstream: {out.row_id} did not succeed ({out.status}: {type(out.exc).__name__}: {str(out.exc)[:200]})")
    return out.result


def control_result(rid, tolerance=None):
    out = run_control(rid, tolerance)
    if out.status == "not_converged":
        raise NotRun(f"control_not_converged: {getattr(out.exc, 'reason', None)}")
    need(out.status == "success", f"control: the KE control of {rid} failed ({out.status}: "
                                  f"{type(out.exc).__name__}: {str(out.exc)[:200]})")
    return out.result


def failure(out):
    need(out.exc is not None and out.status in ("error", "not_converged"),
         f"upstream: {out.row_id} raised no registered exception (status {out.status}: "
         f"{type(out.exc).__name__ if out.exc is not None else None}: {str(out.exc)[:200]})")
    return out.exc


# ---------------------------------------------------------------- reading results
def region_of(key_region):
    """'bin[3]' -> ('bin', 3); 'outside.low_binding' -> ('outside', 'low_binding')."""
    if key_region.startswith("bin["):
        return "bin", int(key_region[4:-1])
    return "outside", key_region.split(".", 1)[1]


def to_kinetic(rid, region):
    kind, where = region
    m = rc_row(rid)["map"]
    if kind == "bin":
        return "bin", m["kinetic_bin_of_binding_bin"][where]
    return "outside", m["outside"][where]


def order_of(result, n):
    for element in G(result, "orders"):
        if G(element, "order") == n:
            return element
    raise Unreadable(f"order {n} is not in result.orders")


_KEY_ORDER_REGION = re.compile(r"^order\[(\d+)\]\.(bin\[\d+\]|outside\.\w+)\.(total|no_hole|hole)(_error)?(?:_parts\[(\w+)\])?$")
_KEY_ORDER_AREA = re.compile(r"^order\[(\d+)\]\.(area|area_error|area_zero_kind|cutoff_sensitivity|cutoff_sensitivity_error)\.(\w+)$")
_KEY_ORDER_PARTS = re.compile(r"^order\[(\d+)\]\.area_error_parts\[(\w+)\]\[(\w+)\]$")
_KEY_TOP = re.compile(r"^(zero_loss|zero_loss_error|spectrum|spectrum_error)\.(bin\[\d+\]|outside\.\w+)$")
_KEY_TRUNC = re.compile(r"^truncation\.(bound_reason|neglected_at_or_above_region_bound|tau_per_eV)\[bin\[(\d+)\]\]$")


def read_value(result, key, region_map=None, kinetic=False):
    """The scalar a physical ledger key names.  ``region_map`` maps the BE region of the key to the one read;
    ``kinetic`` reads a KE result (outside zero loss inside single_loss, prereg 4.3)."""
    mapped = (lambda r: r) if region_map is None else region_map
    m = _KEY_ORDER_REGION.match(key)
    if m:
        n, region, comp, err, part = int(m[1]), mapped(region_of(m[2])), m[3], m[4], m[5]
        o = order_of(result, n)
        name = comp + ("_error_parts" if part else ("_error" if err else ""))
        holder = o if region[0] == "bin" else G(G(o, "outside"), region[1])
        value = G(holder, name)
        if part:
            value = G(value, part)
        return float(value[region[1]]) if region[0] == "bin" else float(value)
    m = _KEY_ORDER_AREA.match(key)
    if m:
        return G(G(order_of(result, int(m[1])), m[2]), m[3])
    m = _KEY_ORDER_PARTS.match(key)
    if m:
        return G(G(G(order_of(result, int(m[1])), "area_error_parts"), m[2]), m[3])
    m = _KEY_TOP.match(key)
    if m:
        name, region = m[1], mapped(region_of(m[2]))
        if region[0] == "bin":
            return float(G(result, name)[region[1]])
        if kinetic and name.startswith("zero_loss"):
            return float(G(G(G(G(result, "single_loss"), "outside"), region[1]), name))
        return float(G(G(G(result, "outside"), region[1]), name))
    m = _KEY_TRUNC.match(key)
    if m:
        index = mapped(("bin", int(m[2])))[1]
        return G(G(G(result, "truncation"), m[1]), index)
    if key == "convention_ratio":
        return G(result, "convention_ratio")
    if key in ("order_area_ratio.value", "order_area_ratio.error"):
        return G(G(result, "order_area_ratio"), key.split(".")[1])
    m = re.match(r"^order_area_ratio\.interval\[(\d)\]$", key)
    if m:
        return G(G(result, "order_area_ratio"), "interval")[int(m[1])]
    raise Unreadable(f"no reader for the key {key!r}")


def triple(result, i):
    t = G(result, "truncation")
    return (G(t, "neglected_at_or_above_region_bound")[i], G(t, "bound_reason")[i], G(t, "tau_per_eV")[i])


def is_complete(tr):
    b, reason, tau = tr
    return isinstance(b, float) and b == 0.0 and reason == "complete" and tau is None


def finite_positive(tr, label):
    b, reason, tau = tr
    need(isinstance(b, float) and math.isfinite(b) and b > 0.0 and reason is None
         and isinstance(tau, float) and math.isfinite(tau), f"{label}: {tr} is not (finite > 0, None, finite tau)")


def warning_codes(result):
    return [G(w, "code") for w in G(G(result, "diagnostics"), "warnings")]


def energy_axis_of(result):
    ea = G(result, "energy_axis")
    need(isinstance(ea, dict), f"result.energy_axis is {type(ea).__name__}, not a dict")
    return ea


def kinetic_coordinates(ea):
    p = G(ea, "kinetic_primary")
    return np.asarray(G(p, "energies_eV") if hasattr(p, "energies_eV") else G(p, "edges_eV"), float)


def n_bins_of(row):
    i = row["inputs"]
    spec = i["spectrum"]["output_binding_edges_eV"] if row["entry"] == "cli" else i["output_edges_eV"]
    return len(be_edges(spec)[0]) - 1


def physical_snapshot(result, permutation=None, lists=False):
    """The physical keys of prereg 7.1 of a BE result as one comparable record.

    ``permutation[i]`` is the display bin of this result placed at position i."""
    n = len(G(result, "zero_loss"))
    idx = list(range(n)) if permutation is None else list(permutation)

    def arr(a):
        a = np.asarray(a, float)
        return [float(a[i]).hex() for i in idx]
    orders = []
    for o in G(result, "orders"):
        rec = dict(order=G(o, "order"), outside=canon(G(o, "outside")), mesh=canon(G(o, "mesh")))
        for c in COMPONENTS:
            rec[c], rec[c + "_error"] = arr(G(o, c)), arr(G(o, c + "_error"))
            rec[c + "_error_parts"] = {p: arr(G(G(o, c + "_error_parts"), p)) for p in PARTS}
        for name in ("area", "area_error", "area_error_parts", "area_zero_kind", "cutoff_sensitivity",
                     "cutoff_sensitivity_error"):
            rec[name] = canon(G(o, name))
        orders.append(rec)
    ea = energy_axis_of(result)
    t = G(result, "truncation")
    snap = dict(orders=orders, outside=canon(G(result, "outside")),
                **{name: arr(G(result, name)) for name in ("zero_loss", "zero_loss_error", "spectrum", "spectrum_error")},
                order_area_ratio=canon(G(result, "order_area_ratio")), convention_ratio=G(result, "convention_ratio"),
                single_loss=canon(G(result, "single_loss")), neglected_total_area=canon(G(t, "neglected_total_area")),
                kinetic_output_edges_eV=canon(np.asarray(G(ea, "kinetic_output_edges_eV"), float)),
                kinetic_primary=canon(G(ea, "kinetic_primary")))
    if lists:
        for name in ("neglected_at_or_above_region_bound", "bound_reason", "tau_per_eV"):
            snap[name] = [canon(G(t, name)[i]) for i in idx]
    return json.loads(json.dumps(snap))


def kernel_evaluations(row, out):
    if row["entry"].startswith("analytic"):
        need(out.kernel is not None, "the kernel of the row was not built (the call did not reach the entry)")
        return out.kernel.calls
    return len(out.calls)


# ---------------------------------------------------------------- baseline child (commit 8c33d29)
def base_item(rid):
    """A BASE row as plain data (no object of the package)."""
    row = rows()[rid]
    i = row["inputs"]
    item = dict(entry=row["entry"], inputs=i, settings=row["settings"])
    name = i.get("material") or i.get("config_material")
    if name:
        item["material"] = material_spec(name)
    if row["entry"] == "cli":
        item["spectrum"] = spectrum_json(row)
    else:
        if "primary" in i:
            item["primary"] = list(ke_primary_arrays(i["primary"]))
            item["edges"] = hexes(table()["kinetic_output_edges"][i["output_edges_eV"]])
        if "kernel" in i:
            item["window"] = list(kernel_window(i["kernel"]))
            item["features"] = features_of(i["kernel"])
    return item


def base_call(item, explicit_none=False):
    """Run a BASE row from its plain data; the same code runs in the parent and in the baseline child."""
    import xpsloss
    from xpsloss import cli, xps_single_loss as xs
    i, entry = item["inputs"], item["entry"]
    extra = {"energy_axis": None} if explicit_none else {}
    record = {}
    if entry == "cli":
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            cfg, spectrum, target = folder/"input.json", folder/"spectrum.json", folder/"output"
            cfg.write_bytes(json.dumps(dict(material=material_dict(item["material"]),
                                            electron_energy_eV=i["config_electron_energy_eV"],
                                            loss_grid_eV={"min": 1, "max": 100, "step": 1},
                                            integration=i["config_integration"])).encode("ascii"))
            spectrum.write_bytes(json.dumps(item["spectrum"]).encode("ascii"))
            got = cli_run(cli, [i["command"], str(cfg), "--spectrum", str(spectrum), *i["argv"], "--out", str(target)],
                          target)
            record.update(exit_code=got["exit_code"], stdout=got["stdout"], stderr=got["stderr"],
                          files={k: v.decode("utf-8") for k, v in got["files"].items()},
                          config_path=str(cfg), spectrum_path=str(spectrum), output_path=str(target),
                          config_sha256=hashlib.sha256(cfg.read_bytes()).hexdigest(),
                          spectrum_sha256=hashlib.sha256(spectrum.read_bytes()).hexdigest())
        return record
    calls, original = [], xs._call_kernel

    def recorder(material_, electron_energy_eV, loss_eV, **kw):
        calls.append(float(electron_energy_eV))
        return original(material_, electron_energy_eV, loss_eV, **kw)
    kernel = None
    xs._call_kernel = recorder
    try:
        settings = {k: num(v) for k, v in item["settings"].items()}
        if entry == "analytic_spectrum":
            lo = importlib.import_module("xpsloss._loss_orders")
            kernel = AnalyticKernel(i["kernel"], float(i["hole_charge"]))
            f = item["features"]
            features = dict(jumps_eV=tuple(f["jumps_eV"]), kinks_eV=tuple(f["kinks_eV"]),
                            singularities=tuple(tuple(s) for s in f["singularities"]),
                            min_half_width_eV=f["min_half_width_eV"], min_variation_scale_eV=f["min_variation_scale_eV"])
            result, panels = lo.integrate_orders(
                kernel, tuple(item["window"]), features=features, weight_length_A=float(i["weight_length_A"]),
                hole_charge=float(i["hole_charge"]), primary=make_ke_primary(xpsloss, *item["primary"]),
                output_edges_eV=np.array(item["edges"], float), model=i["model"], max_order=i["max_order"],
                **settings, **extra)
            record["panels"] = canon(panels)
        else:
            kw = dict(loss_window_eV=tuple(num(x) for x in i["loss_window_eV"]),
                      weight_length_A=float(i["weight_length_A"]), epsilon_momentum=i["epsilon_momentum"], **settings)
            if "multiplier" in i:
                kw["multiplier"] = i["multiplier"]
            mat = make_material(xpsloss, item["material"])
            energy = float(i["electron_energy_eV"])
            primary = make_ke_primary(xpsloss, *item["primary"]) if "primary" in item else None
            edges = np.array(item["edges"], float) if "edges" in item else None
            if entry == "api_spectrum":
                result = xpsloss.multiple_loss_spectrum(mat, energy, primary, edges, model=i["model"],
                                                        max_order=i["max_order"], **kw, **extra)
            elif entry == "api_areas":
                result = xpsloss.multiple_loss_order_areas(mat, energy, model=i["model"], max_order=i["max_order"],
                                                           **kw, **extra)
            elif entry == "api_single_spectrum":
                result = xpsloss.single_loss_spectrum(mat, energy, primary, edges, **kw, **extra)
            else:
                result = xpsloss.loss_window_area(mat, energy, **kw, **extra)
        record["result"] = canon(result)
    except Exception as exc:
        record["result"] = canon(exc)
        record["raised"] = type(exc).__name__
    finally:
        xs._call_kernel = original
    record["kernel_calls"] = kernel.calls if kernel is not None else len(calls)
    return json.loads(json.dumps(record))


def _baseline_child(out_path, payload_path):
    """Runs inside the unpacked baseline commit: old names only."""
    import scipy
    import xpsloss
    from xpsloss import cli, xps_single_loss as xs
    payload = json.loads(Path(payload_path).read_text(encoding="utf-8"))
    record = dict(file=str(Path(xpsloss.__file__).resolve()), module=str(Path(xs.__file__).resolve()),
                  version=xpsloss.__version__, executable=sys.executable, numpy=np.__version__, scipy=scipy.__version__,
                  has_new_api=hasattr(xpsloss, "BindingEnergyAxis"), rows={})
    for rid, item in payload["rows"].items():
        start = time.perf_counter()
        record["rows"][rid] = base_call(item)
        record["rows"][rid]["seconds"] = time.perf_counter() - start
    record["capabilities"] = cli_run(cli, ["capabilities"], Path(payload_path).parent/"none")
    record["capabilities"].pop("files")
    Path(out_path).write_text(json.dumps(record), encoding="utf-8")


def baseline():
    """Snapshots of the baseline commit for the selected BASE rows and BC11: one child process."""
    if "baseline" in _CACHE:
        return _CACHE["baseline"]
    box = _CACHE["baseline"] = {"unavailable": None, "seconds": None}
    if not (ROOT/".git").exists():
        box["unavailable"] = ("not_run", NO_HISTORY)
        return box
    try:
        commit = subprocess.check_output(["git", "rev-parse", "--verify", BASELINE_COMMIT + "^{commit}"], cwd=ROOT,
                                         text=True, stderr=subprocess.STDOUT).strip()
    except Exception as exc:
        box["unavailable"] = ("fail", f"baseline: commit {BASELINE_COMMIT} is not available in the checkout: {exc}")
        return box
    if commit != BASELINE_COMMIT or commit != table()["baseline"]["commit"]:
        box["unavailable"] = ("fail", f"baseline: resolved {commit}, not the registered commit")
        return box
    ids = [rid for rid in selected_ids() if rows()[rid]["group"] == "BASE"]
    payload = dict(rows={rid: base_item(rid) for rid in ids})
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
            box["unavailable"] = ("fail", f"baseline: the child failed: {done.stderr[-600:]}")
            return box
        record = json.loads((folder/"out.json").read_text(encoding="utf-8"))
        inside = Path(record["file"]).resolve().is_relative_to(src.resolve()) and \
            Path(record["module"]).resolve().is_relative_to(src.resolve())
    import scipy
    box.update(record=record, seconds=time.perf_counter() - start, from_archive=inside)
    if (not inside or record["has_new_api"] or record["version"] != table()["baseline"]["package_version"]
            or record["executable"] != sys.executable or record["numpy"] != np.__version__
            or record["scipy"] != scipy.__version__):
        box["unavailable"] = ("fail", f"baseline: child environment: {record['file']}, version {record['version']}, "
                                      f"numpy {record['numpy']}, scipy {record['scipy']}")
    return box


def baseline_row(rid):
    box = baseline()
    if box["unavailable"]:
        status, detail = box["unavailable"]
        if status == "not_run":
            raise NotRun(detail)
        raise AssertionError(detail)
    return box["record"]["rows"][rid]


def base_current(rid, explicit_none=False):
    """The BASE row on the code under test (cached)."""
    runs = _CACHE.setdefault("base_current", {})
    if (rid, explicit_none) not in runs:
        try:
            runs[(rid, explicit_none)] = base_call(base_item(rid), explicit_none)
        except Exception:
            runs[(rid, explicit_none)] = {"harness_error": traceback.format_exc()}
    got = runs[(rid, explicit_none)]
    need("harness_error" not in got, f"upstream: the row could not be run: {got.get('harness_error', '')[-400:]}")
    return got


def without_axis_field(snapshot):
    """Drop the one allowed addition: 'energy_axis' directly under a spectrum result, and only when it is None."""
    if isinstance(snapshot, dict) and snapshot.get("@type") == "MultipleLossSpectrumResult" and "energy_axis" in snapshot:
        need(snapshot["energy_axis"] is None, "result.energy_axis is not None on the KE path")
        return {k: v for k, v in snapshot.items() if k != "energy_axis"}
    return snapshot


def has_key(obj, name):
    if isinstance(obj, dict):
        return name in obj or any(has_key(v, name) for v in obj.values())
    if isinstance(obj, list):
        return any(has_key(v, name) for v in obj)
    return False


def run_base(row, out):
    """BASE rows run through base_call so that both sides execute the same code."""
    got = base_current(row["id"])
    out.extras["base"] = got
    if row["entry"] == "cli":
        out.status = "complete" if got["exit_code"] == 0 else "error"
    elif "raised" in got:
        import xpsloss
        cls = getattr(xpsloss, got["raised"], None)
        convergence = isinstance(cls, type) and issubclass(cls, xpsloss.NumericalConvergenceError)
        out.status = "not_converged" if convergence else "error"
    else:
        out.status = "success"


# ---------------------------------------------------------------- judgements
class Judgement:
    def __init__(self):
        self.items, self.faults, self.seconds = {}, {}, {}

    def add(self, rid, key, kind, oracle, status, detail=""):
        if (rid, key, kind, oracle) in self.items:
            raise AssertionError(f"duplicate judgement {(rid, key, kind, oracle)}")
        self.items[(rid, key, kind, oracle)] = (status, detail)

    def counts(self):
        out = {}
        for (_, _, kind, _), (status, _) in self.items.items():
            out.setdefault(kind, {}).setdefault(status, 0)
            out[kind][status] += 1
        return out

    def by_status(self):
        out = {}
        for status, _ in self.items.values():
            out[status] = out.get(status, 0) + 1
        return out

    def categories(self):
        """Failures by cause (prereg 7.3): upstream, control, baseline, reference, other."""
        out = {}
        for (status, detail) in self.items.values():
            if status == "pass":
                continue
            if status != "fail":
                name = f"{status}: {detail.split(':')[0]}"
            else:
                head = detail.split(":")[0]
                name = head if head in ("upstream", "control", "baseline", "reference") else "comparison"
            out[name] = out.get(name, 0) + 1
        return out

    def not_passed(self, rid):
        return {k: v for k, v in self.items.items() if k[0] == rid and v[0] != "pass"}


def guard(check):
    """Run a check; (status, detail).  An unreadable field or any exception is a failure."""
    try:
        got = check()
        if got is None or got is True:
            return "pass", ""
        return ("pass", "") if got else ("fail", "condition is false")
    except NotRun as exc:
        return exc.status, exc.reason
    except AssertionError as exc:
        return "fail", str(exc)[:500]
    except Exception as exc:
        return "fail", f"{type(exc).__name__}: {str(exc)[:400]}"


def expect_of(row):
    return {**row["expect"], **POST_FREEZE_EXPECT.get(row["id"], {})}


# ---- state (oracle expect)
def check_state(row, out, key, use):
    e, rid = expect_of(row), row["id"]
    if key == "terminal_state":
        need(out.status == e["status"], f"status {out.status!r}, expected {e['status']!r}"
             + (f" ({type(out.exc).__name__}: {str(out.exc)[:300]})" if out.exc is not None else "")
             + (f" [{(out.harness_error or '')[-300:]}]" if out.status == "harness_error" else ""))
        return
    if row["group"] == "BASE":
        return check_base_state(row, out, key)
    if row["entry"] in ("cli", "cli_capabilities"):
        return check_cli_state(row, out, key)
    exc = failure(out)
    if key == "exception_type":
        import xpsloss
        need(type(exc).__name__ == e["exception"] and isinstance(exc, getattr(xpsloss, e["exception"])),
             f"{type(exc).__name__}: {str(exc)[:200]}")
    elif key == "reason":
        need(G(exc, "reason") == e["reason"], f"reason {G(exc, 'reason')!r}")
    elif key in ("field", "detail", "resource", "stage"):
        need(G(exc, key) == e[key], f"{key} {G(exc, key)!r}, expected {e[key]!r}")
    elif key in ("limit", "lower_bound", "order_limit", "requested_order"):
        need(type(G(exc, key)) is int and G(exc, key) == e[key], f"{key} {G(exc, key)!r}, expected {e[key]!r}")
    elif key == "kernel_evaluations":
        need(kernel_evaluations(row, out) == e["kernel_evaluations"], f"{kernel_evaluations(row, out)} kernel calls")
    elif key in ("binding_eV", "kinetic_eV"):
        got = G(exc, key)
        need(isinstance(got, tuple) and len(got) == 2, f"{key} is not a tuple of two: {got!r}")
        need([float(x).hex() for x in got] == [x.hex() for x in hexes(e[key])], f"{key} {got!r}")
        need([float(x).hex() for x in got] == [rv(x).hex() for x in rc_row(rid)["collision"][key]], f"{key} against R-C")
    elif key == "roles":
        got = G(exc, "roles")
        need(isinstance(got, tuple) and list(got) == e["roles"] == rc_row(rid)["collision"]["roles"], f"roles {got!r}")
    elif key.startswith("required_order_estimate =="):
        need(G(exc, "required_order_estimate") == e["required_order_estimate"] == int(key.split("==")[1]),
             f"required_order_estimate {G(exc, 'required_order_estimate')!r}")
    elif key.startswith("required_order_estimate >="):
        got = G(exc, "required_order_estimate")
        need(type(got) is int and got >= e["required_order_estimate_at_least"] == int(key.split(">=")[1]),
             f"required_order_estimate {got!r}")
    else:
        raise Unreadable(f"no state check for the key {key!r}")


def check_cli_state(row, out, key):
    e, record = expect_of(row), out.cli
    need(record is not None, f"upstream: the CLI row was not run ({out.status}: {(out.harness_error or '')[-300:]})")
    payload = record["stderr_json"]
    if key == "exit_code":
        need(record["exit_code"] == e["exit_code"], f"exit code {record['exit_code']!r}; stderr {record['stderr'][-300:]}")
    elif key == "output_dir_after == absent":
        need(record["output_dir_after"] == "absent" == e["output_dir_after"], record["output_dir_after"])
    elif key == "stderr is one JSON object":
        need(isinstance(payload, dict) and record["stderr"].isascii(), f"stderr {record['stderr'][:200]!r}")
        need(set(payload) >= {"status", "error_type", "reason", "message"}, f"keys {sorted(payload)}")
    else:
        need(isinstance(payload, dict), f"stderr is not one JSON object: {record['stderr'][:200]!r}")
        if key == "error_type":
            need(G(payload, "error_type") == e["error_type"], f"error_type {payload.get('error_type')!r}")
        elif key in ("reason", "field", "detail"):
            need(G(payload, key) == e[key], f"{key} {payload.get(key)!r}, expected {e[key]!r}")
        elif key == "detail attribute is a non-empty ASCII string":
            d = G(payload, "detail")
            need(isinstance(d, str) and d.strip() != "" and d.isascii(), f"detail {d!r}")
        elif key in ("binding_eV", "kinetic_eV"):
            got = G(payload, key)
            need(isinstance(got, list) and len(got) == 2, f"{key} is not an array of two: {got!r}")
            need([float(x).hex() for x in got] == [x.hex() for x in hexes(e[key])], f"{key} {got!r}")
            need([float(x).hex() for x in got] == [rv(x).hex() for x in rc_row(row["id"])["collision"][key]],
                 f"{key} against R-C")
        elif key == "roles":
            need(G(payload, "roles") == e["roles"] == rc_row(row["id"])["collision"]["roles"], f"roles {payload.get('roles')!r}")
        elif key == "kernel_evaluations == 0":
            need(len(out.calls) == 0, f"{len(out.calls)} kernel calls")
        else:
            raise Unreadable(f"no CLI state check for the key {key!r}")


def check_base_state(row, out, key):
    e = expect_of(row)
    got = base_current(row["id"])
    if key == "exit_code == 0 on both sides":
        need(got["exit_code"] == 0, f"exit code {got['exit_code']}: {got['stderr'][-300:]}")
        need(baseline_row(row["id"])["exit_code"] == 0, "baseline exit code")
        return
    need("raised" in got, "upstream: the row raised nothing")
    attributes = got["result"]["attributes"]
    if key == "exception_type":
        need(got["raised"] == e["exception"], got["raised"])
    elif key == "kernel_evaluations":
        need(got["kernel_calls"] == e["kernel_evaluations"], f"{got['kernel_calls']} kernel calls")
    elif key in ("reason", "resource", "limit", "lower_bound"):
        need(attributes.get(key) == e[key], f"{key} {attributes.get(key)!r}")
    else:
        raise Unreadable(f"no BASE state check for the key {key!r}")


# ---- exact conversion (oracle R-C, or None for the analytic entry)
def check_conversion(row, out, key, oracle, use):
    ref = use.get(key) if oracle == "R-C" else None      # read first: an upstream failure is not an unused key
    if key.startswith("run.json: "):
        ea = G(manifest_of(out), "energy_axis")
        name = key[len("run.json: energy_axis."):]
    else:
        ea = energy_axis_of(success(out))
        name = key[len("energy_axis."):]
    if oracle == "None (analytic)":
        need(name == "electron_energy_minus_reference_eV" and G(ea, name) is None, f"{name} is {ea.get(name)!r}, not None")
        return
    need(oracle == "R-C", oracle)
    m = re.match(r"^(kinetic_output_edges_eV|kinetic_primary\.coordinates)\[(\d+)\]$", name)
    if m:
        values = kinetic_coordinates(ea) if m[1].startswith("kinetic_primary") else G(ea, "kinetic_output_edges_eV")
        got = float(values[int(m[2])])
        need(got.hex() == rv(ref).hex(), f"{got.hex()} != {rv(ref).hex()}")
    elif name == "conversion_constant_exact":
        got = G(ea, name)
        need(isinstance(got, str) and re.fullmatch(r"-?\d+/[1-9]\d*", got) is not None, f"form of {got!r}")
        need(Fr(got) == Fr(ref) and Fr(got).denominator == int(got.split("/")[1]), f"{got} != {ref} (irreducible)")
    elif name in ("support_order_binding", "support_order_kinetic"):
        need(type(G(ea, name)) is int and G(ea, name) == ref, f"{name} {ea.get(name)!r} != {ref}")
    else:
        got, want = G(ea, name), rv(ref)
        if want is None:
            need(got is None, f"{name} is {got!r}, not None")
        else:
            need(isinstance(got, float) and got.hex() == want.hex(), f"{name} {got!r} != {want.hex()}")


# ---- zero loss membership (oracle R-C, KE widths)
def check_membership(row, out, key, use):
    r = Fr(use.get(key))
    if key.startswith("CSV "):
        table_ = spectrum_csv(out)
        region = region_of(key[len("CSV zero_loss."):])
        name = str(region[1])
        line = table_["by_region"].get(name)
        need(line is not None, f"no CSV row {name!r}")
        y, err = float(line["zero_loss"]), float(line["zero_loss_error"])
    else:
        result = success(out)
        y = read_value(result, key)
        err = read_value(result, key.replace("zero_loss.", "zero_loss_error.", 1))
    need(within_membership(y, r), f"zero loss {y!r} against {r}")
    if r == 0:
        need(err == 0.0, f"zero-loss error {err!r} where the reference is 0")
    else:
        need(math.isfinite(err) and err >= 0.0, f"zero-loss error {err!r}")


# ---- bit identity
def check_bit_identity(row, out, key, oracle):
    rid = row["id"]
    result = success(out)
    if oracle == "KE control":
        control = control_result(rid)
        got = read_value(result, key)
        want = read_value(control, key, region_map=lambda region: to_kinetic(rid, region), kinetic=True)
        need(bits(got) == bits(want), f"{bits(got)} != control {bits(want)}")
        return
    need(oracle.endswith(".result"), oracle)
    other = success(run_row(oracle[:-len(".result")]))
    name = key.split(" == ")[0]
    got, want = read_value(result, name), read_value(other, name)
    need(bits(got) == bits(want), f"{bits(got)} != {oracle} {bits(want)}")


# ---- records
def strip_diagnostics(diag):
    d = {k: v for k, v in dict(diag).items() if k != "warnings"}
    if "effective_settings" in d:
        d["effective_settings"] = {k: v for k, v in dict(d["effective_settings"]).items()
                                   if k != "neglect_tolerance_fraction"}
    return d


def check_record(row, out, key, oracle, use):
    rid = row["id"]
    ref = use.get(key) if oracle.startswith("R-C") else None
    if row["group"] == "BASE":
        return check_base_record(row, out, key, oracle)
    if oracle == "replay":
        return check_replay(row, out, key)
    if oracle in ("BS4.result", "BT2c.result") and key.startswith("the three per-bin lists of the exception"):
        exc = failure(out)
        t = G(G(exc, "diagnostics"), "truncation")
        other = G(success(run_row(oracle[:-len(".result")])), "truncation")
        for name in ("neglected_at_or_above_region_bound", "bound_reason", "tau_per_eV"):
            same(list(G(t, name)), list(G(other, name)), name)
        return
    result = success(out)
    if oracle == "KE control":
        control = control_result(rid)
        if key.startswith("order["):
            n = int(key[6:key.index("]")])
            same(G(order_of(result, n), "mesh"), G(order_of(control, n), "mesh"), "mesh")
        elif key.startswith("single_loss"):
            same(G(result, "single_loss"), G(control, "single_loss"), "single_loss")
        elif key == "truncation.neglected_total_area":
            same(G(G(result, "truncation"), "neglected_total_area"), G(G(control, "truncation"), "neglected_total_area"), key)
        elif key.startswith("diagnostics"):
            mine, theirs = strip_diagnostics(G(result, "diagnostics")), strip_diagnostics(G(control, "diagnostics"))
            need(set(theirs) <= set(mine), f"diagnostics lacks {sorted(set(theirs) - set(mine))}")
            same({k: mine[k] for k in theirs}, theirs, "diagnostics")
        elif key.startswith("warnings"):
            drop = set(BE_CODES) | {"truncated"}
            mine = [(G(w, "source"), G(w, "code")) for w in G(G(result, "diagnostics"), "warnings") if G(w, "code") not in drop]
            theirs = [(G(w, "source"), G(w, "code")) for w in G(G(control, "diagnostics"), "warnings") if G(w, "code") not in drop]
            need(mine == theirs, f"warning identifiers {mine} != control {theirs}")
        else:
            raise Unreadable(f"no record check for {key!r}")
        return
    if oracle == "R-C (both-sided rule)":
        got = G(G(energy_axis_of(result), "support_comparison"), "first_mismatches")
        got = [dict(x) for x in got]
        want = [{k: rv(v) for k, v in x.items()} for x in ref]
        same(got, want, "first_mismatches")
        return
    if oracle == "inputs":
        return check_input_record(row, out, key, result)
    if oracle.endswith(".result"):
        other = success(run_row(oracle[:-len(".result")]))
        n = n_bins_of(row)
        if "bin[i] -> bin[" in key:
            mine = physical_snapshot(result, permutation=[n - 1 - i for i in range(n)])
            theirs = physical_snapshot(other)
        else:
            mine, theirs = physical_snapshot(result, lists=True), physical_snapshot(other, lists=True)
        d = first_difference(mine, theirs, "physical")
        need(d is None, d)
        return
    raise Unreadable(f"no record check for {key!r} with the oracle {oracle!r}")


def check_input_record(row, out, key, result):
    ea = energy_axis_of(result)
    axis = out.extras["axis_plain"]
    ref = axis["reference"]
    if key.startswith("energy_axis.axis"):
        need(G(ea, "axis") == "binding_energy_eV" and G(ea, "binding_reference") == "fermi_level", "axis names")
        need(G(ea, "reference_kind") == ref["kind"], "reference_kind")
        pairs = (("photon_energy_eV", axis["photon_energy_eV"]),
                 ("kinetic_energy_offset_eV", ref.get("kinetic_energy_offset_eV")),
                 ("binding_anchor_eV", ref.get("binding_anchor_eV")), ("kinetic_anchor_eV", ref.get("kinetic_anchor_eV")),
                 ("anchor_line", ref.get("line")), ("anchor_citation", ref.get("citation")),
                 ("primary_reference_binding_eV", axis["primary_reference_binding_eV"]),
                 ("primary_reference_choice", axis["primary_reference_choice"]))
        for name, want in pairs:
            need(bits(G(ea, name)) == bits(want), f"{name}: {ea.get(name)!r} != {want!r}")
    elif key.startswith("energy_axis.direction"):
        edges = out.extras["edges"]
        kind, x, _ = be_primary_arrays(row["inputs"]["primary"])
        want = {"output_edges": "increasing" if edges[0] < edges[-1] else "decreasing",
                "primary": None if kind == "delta" else ("increasing" if x[0] < x[-1] else "decreasing")}
        need(dict(G(ea, "direction")) == want, f"direction {ea.get('direction')!r} != {want}")
        need(want == rc_row(row["id"])["map"]["direction"], "direction against R-C")
    elif key.startswith("energy_axis.input_binding_edges_eV"):
        got = G(ea, "input_binding_edges_eV")
        need(isinstance(got, np.ndarray) and got.dtype == np.float64 and got.tobytes() == out.extras["edges"].tobytes(),
             "input_binding_edges_eV differs from the input (bits or order)")
    elif key.startswith("energy_axis.input_primary"):
        got, given = G(ea, "input_primary"), out.extras["primary"]
        need(type(got) is type(given), f"input_primary is {type(got).__name__}")
        same(got, given, "input_primary")
    elif key.startswith("energy_axis.not_recorded"):
        got = G(ea, "not_recorded")
        need(isinstance(got, tuple) and got == not_recorded_names(axis), f"not_recorded {got!r} != {not_recorded_names(axis)}")
    elif key.startswith("energy_axis.record"):
        got, want = G(ea, "record"), record_fields(axis)
        need(isinstance(got, dict) and list(got) == list(want), f"record keys {list(got)}")
        same(got, want, "record")
        given = out.extras["axis_object"]
        for name, value in want.items():
            if isinstance(value, dict):
                need(got[name] is not getattr(given, name), f"record.{name} is the caller's object, not a copy")
    else:
        raise Unreadable(f"no input record check for {key!r}")


# ---- truncation
def check_truncation(row, out, key, oracle, use):
    rid, e = row["id"], expect_of(row)
    if oracle == "R-C (both-sided rule)":
        ref = use.get(key)
        if key.startswith("R-C support orders"):
            need(ref[0] == ref[1] == 50, f"R-C support orders {ref}")
            return
        if key.startswith("requested_order == R-C"):
            need(G(failure(out), "requested_order") == ref, f"requested_order {G(failure(out), 'requested_order')!r} != {ref}")
            return
        if key.startswith("required_order_estimate == R-C"):
            got = G(failure(out), "required_order_estimate")
            need(got == ref, f"required_order_estimate {got!r} != {ref}")
            return
        result = success(out)
        if key.startswith("warning axis_conversion_support_disagreement"):
            codes = warning_codes(result)
            need(codes.count("axis_conversion_support_disagreement") == (1 if ref else 0),
                 f"warning present {codes.count('axis_conversion_support_disagreement')} times, R-C says {ref}")
            return
        m = re.match(r"^truncation\.bound_class\[bin\[(\d+)\]\]$", key)
        if m:
            got = "complete" if is_complete(triple(result, int(m[1]))) else "reached"
            need(got == ref, f"{triple(result, int(m[1]))} is {got}, R-C says {ref}")
            return
        m = re.match(r"^energy_axis\.region_loss_upper_eV\[bin\[(\d+)\]\]$", key)
        if m:
            got = G(energy_axis_of(result), "region_loss_upper_eV")[int(m[1])]
            need(bits(got) == bits(rv(ref)), f"{got!r} != {rv(ref)!r}")
            return
        if key.startswith("truncation."):
            got = G(G(result, "truncation"), key.split(".", 1)[1])
            need(type(got) is type(ref) and got == ref, f"{key} {got!r} != {ref!r}")
            return
        ea = energy_axis_of(result)
        if key.startswith("energy_axis.support_comparison."):
            got = G(G(ea, "support_comparison"), key.rsplit(".", 1)[1])
            got = list(got) if isinstance(ref, list) else got
            need(got == ref and (isinstance(ref, list) or type(got) is int), f"{key} {got!r} != {ref!r}")
            return
        got = G(ea, key.split(".", 1)[1])
        need(type(got) is int and got == ref, f"{key} {got!r} != {ref!r}")
        return
    if oracle == "same result":
        result = success(out)
        present = warning_codes(result).count("truncated")
        state = G(G(result, "truncation"), "state")
        need(present == (1 if state == "truncated" else 0), f"'truncated' appears {present} times with state {state!r}")
        return
    if oracle == "KE control with the same tolerance":
        control = run_control(rid, "row")
        if key == "required_order_estimate":
            exc = failure(out)
            need(control.exc is not None and type(control.exc) is type(exc)
                 and G(control.exc, "reason") == G(exc, "reason"),
                 f"control: the KE control with the tolerance ended as {control.status}: {control.exc!r}")
            need(G(exc, "required_order_estimate") == G(control.exc, "required_order_estimate"),
                 f"{G(exc, 'required_order_estimate')!r} != control {G(control.exc, 'required_order_estimate')!r}")
        else:
            need(control.status == "success", f"control: the KE control with the tolerance ended as {control.status}")
            got = G(G(success(out), "truncation"), "required_order_estimate")
            want = G(G(control.result, "truncation"), "required_order_estimate")
            need(type(got) is type(want) and got == want, f"{got!r} != control {want!r}")
        return
    if oracle == "inputs":
        result = success(out)
        tol = num(row["settings"]["neglect_tolerance_fraction"])
        need(bits(G(G(result, "truncation"), "neglect_tolerance_fraction")) == bits(tol), "truncation tolerance")
        need(bits(G(G(G(result, "diagnostics"), "effective_settings"), "neglect_tolerance_fraction")) == bits(tol),
             "effective_settings tolerance")
        return
    need(oracle == "expect", oracle)
    if key.startswith("run.json: truncation three fields") or key.startswith("run.json: warning"):
        manifest = manifest_of(out)
        if "warning" in key:
            codes = [G(w, "code") for w in G(G(manifest, "diagnostics"), "warnings")]
            need(codes.count("axis_conversion_support_disagreement") == 1, f"codes {codes}")
        else:
            t = G(manifest, "truncation")
            finite_positive((G(t, "neglected_at_or_above_region_bound")[0], G(t, "bound_reason")[0], G(t, "tau_per_eV")[0]),
                            "run.json bin[0]")
        return
    if key.startswith("API replay:"):
        replay = replay_of(out)
        if "warning" in key:
            need(warning_codes(replay).count("axis_conversion_support_disagreement") == 1, "replay warning")
        else:
            finite_positive(triple(replay, 0), "replay bin[0]")
        return
    result = success(out)
    if key.startswith("truncation three fields of bin["):
        i = int(key[len("truncation three fields of bin["):].split("]")[0])
        finite_positive(triple(result, i), f"bin[{i}]")
    elif key.startswith("truncation.required_order_estimate =="):
        got = G(G(result, "truncation"), "required_order_estimate")
        need(type(got) is int and got == e["required_order_estimate"] == int(key.split("==")[1]), f"{got!r}")
    elif key == "truncation.state":
        need(G(G(result, "truncation"), "state") == e["truncation_state"], f"state {G(G(result, 'truncation'), 'state')!r}")
    else:
        raise Unreadable(f"no truncation check for {key!r}")


def check_bound_bracket(row, out, key):
    """Prereg 7.2: b_K <= b <= the candidate at tau_K with the R-C width, evaluated by the product function."""
    rid = row["id"]
    i = int(_KEY_TRUNC.match(key)[2])
    result = success(out)
    control_out = run_control(rid)
    control = control_result(rid)
    j = to_kinetic(rid, ("bin", i))[1]
    b, reason, tau = triple(result, i)
    b_k, reason_k, tau_k = triple(control, j)
    width = rv(rc_row(rid)["support"]["region_loss_upper_eV"][i])
    need(not rc_row(rid)["support"]["width_equals_kinetic"][i], "the R-C width equals the kinetic width")
    need(len(control_out.truncation_calls) == 1, f"control: _order_bounds.truncation was called "
                                                 f"{len(control_out.truncation_calls)} times")
    args, kw = control_out.truncation_calls[0]
    need(len(args) + len(kw) >= 10, "control: unexpected arguments of _order_bounds.truncation")
    names = ("single", "panels", "geometry", "requested", "n", "support", "ratio", "tol", "negative", "diagnostics")
    a = dict(zip(names, args))
    a.update({k: v for k, v in kw.items() if k in names})
    ob = importlib.import_module("xpsloss._order_bounds")
    if tau_k is None or b_k is None:
        need(reason == reason_k and b == b_k, f"no tau on the KE side ({reason_k}); BE side {reason}, {b!r}")
        return
    geometry = types.SimpleNamespace(points=(width,), edges=(0.0, 0.0))

    def tail_model(taus, count, logk, logr):
        available = logr < 0
        log_tail = np.full(len(taus), np.inf)
        power = ob.log_power(logr[available], count)
        log_tail[available] = ob.log_up(ob.log_up(logk[available] + power) + ob.log_geometric_factor(logr[available]))
        keep = np.asarray(taus) == tau_k
        need(int(keep.sum()) == 1, "tau_K is not a point of the tau grid")
        return np.where(keep, log_tail, np.inf)

    def call(model):
        diag = dict(a["diagnostics"]) if isinstance(a["diagnostics"], dict) else a["diagnostics"]
        t = ob.truncation(a["single"], a["panels"], geometry, a["requested"], a["n"], a["support"], a["ratio"], None,
                          a["negative"], diag, tail_model=model)
        return t["neglected_at_or_above_region_bound"][0], t["bound_reason"][0]
    upper, upper_reason = call(tail_model)
    if upper is None:
        upper, upper_reason = call(None)       # the smallest of the other finite candidates
    if upper is None:
        need(b is None and reason == upper_reason, f"no finite candidate ({upper_reason}); BE side {b!r}, {reason}")
        return
    need(isinstance(b, float) and math.isfinite(b), f"bound {b!r} ({reason})")
    need(b_k <= b <= upper, f"b_K {b_k!r} <= b {b!r} <= upper {upper!r} does not hold")


# ---- diagnostics
def manifest_of(out):
    need(out.cli is not None and out.cli["exit_code"] == 0 and "manifest" in out.cli,
         f"upstream: the CLI did not complete (exit {None if out.cli is None else out.cli['exit_code']}: "
         f"{'' if out.cli is None else out.cli['stderr'][-300:]})")
    return out.cli["manifest"]


def replay_of(out):
    manifest_of(out)
    need("api_replay" in out.extras, f"API replay failed: {out.extras.get('api_replay_error', '')[-400:]}")
    return out.extras["api_replay"]


def parse_csv(raw):
    text = raw.decode("ascii")
    need("\r" not in text, "CSV has carriage returns")
    lines = text.strip("\n").split("\n")
    return lines[0].split(","), [line.split(",") for line in lines[1:]]


def spectrum_csv(out):
    manifest_of(out)
    need("multiple_loss_spectrum.csv" in out.cli["files"], "multiple_loss_spectrum.csv is missing")
    header, lines = parse_csv(out.cli["files"]["multiple_loss_spectrum.csv"])
    need(all(len(line) == len(header) for line in lines), "ragged CSV")
    records = [dict(zip(header, line)) for line in lines]
    return dict(header=header, lines=lines, records=records, by_region={r.get("region"): r for r in records})


def csv_float(text, value, label):
    need(text == "%.17g" % value, f"{label}: CSV {text!r} is not %.17g of {value!r}")
    need(float(text).hex() == float(value).hex(), f"{label}: CSV value does not read back bit for bit")


def edge_hash(values):
    return hashlib.sha256(b"".join(struct.pack("<d", float(x)) for x in values)).hexdigest()


def cli_region_names(row):
    edges = be_edges(row["inputs"]["spectrum"]["output_binding_edges_eV"])[0]
    n = len(edges) - 1
    inner = [str(i) for i in range(n)]
    return (["low_binding"] + inner + ["high_binding"]) if edges[0] < edges[-1] else (["high_binding"] + inner + ["low_binding"])


def check_diagnostic(row, out, key, oracle, use):
    rid, e = row["id"], expect_of(row)
    if row["group"] == "BASE":
        return check_base_diagnostic(row, out, key, oracle)
    if rid == "BC11":
        return check_capabilities(row, out, key)
    if oracle == "prereg 4.4":
        result = success(out)
        records = [(G(w, "source"), G(w, "code"), G(w, "message")) for w in G(G(result, "diagnostics"), "warnings")]
        codes = [c for _s, c, _m in records]
        axis = out.extras["axis_plain"]
        want = ["binding_axis_assumed_calibrated"]
        if not_recorded_names(axis):
            want.append("energy_axis_field_not_recorded")
        angle = axis["emission_angle_deg"]
        if isinstance(angle, (int, float)) and not isinstance(angle, bool) and float(angle) != 0.0:
            want.append("recorded_emission_angle_differs_from_theta")
        if "axis_conversion_support_disagreement" in e["warnings_include"]:
            want.append("axis_conversion_support_disagreement")
        got = [c for c in codes if c in BE_CODES]
        need(got == want, f"BE identifiers {got}, expected {want}")
        need(codes[len(codes) - len(got):] == got, f"the BE identifiers do not follow the M4b-1 identifiers: {codes}")
        need(all(s == "xps_multiple_loss" for s, c, _m in records if c in BE_CODES), "source of a BE warning")
        need(all(c in codes for c in e["warnings_include"]) and not any(c in codes for c in e.get("warnings_exclude", ())),
             f"registered warnings; codes {codes}")
        outer = list(G(result, "warnings"))
        need(outer == [m for _s, _c, m in records], "result.warnings and diagnostics.warnings do not pair in order")
        need(all(isinstance(m, str) and m.isascii() and m.strip() for m in outer), "a warning message is not ASCII text")
        if "energy_axis_field_not_recorded" in want:
            message = next(m for _s, c, m in records if c == "energy_axis_field_not_recorded")
            need(all(name in message for name in not_recorded_names(axis)), f"the message lacks a field name: {message!r}")
        return
    if key == "energy_axis.not_recorded":
        got = G(energy_axis_of(success(out)), "not_recorded")
        need(isinstance(got, tuple) and list(got) == e["not_recorded"], f"not_recorded {got!r}")
        return
    if key == "run.json: energy_axis.not_recorded":
        got = G(G(manifest_of(out), "energy_axis"), "not_recorded")
        want = list(not_recorded_names(out.cli["spectrum_data"]["measurement"]))
        need(got == want and len(want) > 0, f"not_recorded {got!r} != {want}")
        return
    if key.startswith("kernel energy argument"):
        success(out)
        want = float(row["inputs"]["electron_energy_eV"])
        need(len(out.calls) > 0 and all(x.hex() == want.hex() for x in out.calls),
             f"kernel energies {sorted(set(out.calls))} in {len(out.calls)} calls, expected {want}")
        return
    return check_cli_diagnostic(row, out, key)


def check_cli_diagnostic(row, out, key):
    rid, e = row["id"], expect_of(row)
    manifest = manifest_of(out)
    record, ref = out.cli, rc_row(rid)
    data = record["spectrum_data"]
    edges = be_edges(row["inputs"]["spectrum"]["output_binding_edges_eV"])[0]
    ke_edges = [rv(x) for x in ref["control"]["kinetic_output_edges_eV"]]
    if key == "files == the three names":
        need(set(record["files"]) == FILES, f"files {sorted(record['files'])}")
    elif key == "header of multiple_loss_spectrum.csv":
        orders = range(1, int(G(G(manifest, "truncation"), "computed")) + 1)
        want = ["region", "binding_edge_low_eV", "binding_edge_high_eV", "kinetic_edge_low_eV", "kinetic_edge_high_eV",
                "zero_loss", "zero_loss_error"]
        for n in orders:
            want += [f"order{n}_{c}{s}" for c in COMPONENTS for s in ("", "_error")]
        need(spectrum_csv(out)["header"] == want + ["spectrum", "spectrum_error"], f"header {spectrum_csv(out)['header']}")
    elif key == "row order of multiple_loss_spectrum.csv":
        got = [line[0] for line in spectrum_csv(out)["lines"]]
        need(got == cli_region_names(row), f"row order {got}")
    elif key.startswith("the four edge columns"):
        by = spectrum_csv(out)["by_region"]
        kinetic = {b: rv(c["kinetic_eV"]) for c in ref["conversion"]["coordinates"] for b in [rv(c["binding_eV"])]}
        names = ("binding_edge_low_eV", "binding_edge_high_eV", "kinetic_edge_low_eV", "kinetic_edge_high_eV")
        for i in range(len(edges) - 1):
            lo, hi = sorted((edges[i], edges[i + 1]))
            for name, value in zip(names, (lo, hi, kinetic[hi], kinetic[lo])):
                csv_float(G(by[str(i)], name), value, f"bin {i} {name}")
        low, high = min(edges), max(edges)
        cells = [G(by["low_binding"], n) for n in names]
        need(cells[0] == "" and cells[3] == "", f"low_binding: the infinite sides are not blank: {cells}")
        csv_float(cells[1], low, "low_binding upper BE")
        csv_float(cells[2], kinetic[low], "low_binding lower KE")
        cells = [G(by["high_binding"], n) for n in names]
        need(cells[1] == "" and cells[2] == "", f"high_binding: the infinite sides are not blank: {cells}")
        csv_float(cells[0], high, "high_binding lower BE")
        csv_float(cells[3], kinetic[high], "high_binding upper KE")
    elif key.startswith("run.json: effective_settings has the"):
        need(set(G(manifest, "effective_settings")) == set(SETTING_NAMES) and len(SETTING_NAMES) == 13,
             f"effective_settings {sorted(G(manifest, 'effective_settings'))}")
    elif key == "run.json: energy_axis == results.energy_axis":
        need(isinstance(G(manifest, "energy_axis"), dict), "energy_axis is not an object")
        same(G(manifest, "energy_axis"), G(G(manifest, "results"), "energy_axis"), "energy_axis")
    elif key == "run.json: measurement == the file":
        same(G(manifest, "measurement"), data["measurement"], "measurement")
    elif key == "run.json: model_id":
        need(G(manifest, "model_id") == table()["model_id"], f"model_id {manifest.get('model_id')!r}")
    elif key.startswith("run.json: primary and output_edges_eV are the BE input record"):
        p, oe = G(manifest, "primary"), G(manifest, "output_edges_eV")
        position = "binding_energies_eV" if data["primary"]["kind"] == "delta" else "binding_edges_eV"
        need(isinstance(p, dict) and G(p, "kind") == data["primary"]["kind"], "primary.kind")
        need([float(x).hex() for x in G(p, position)] == [float(x).hex() for x in data["primary"][position]],
             "the BE positions of the primary differ from the input (bits or order)")
        need("energies_eV" not in p and "edges_eV" not in p, "the top-level primary carries KE positions")
        need([float(x).hex() for x in G(p, "intensities")] == [float(x).hex() for x in data["primary"]["intensities"]],
             "primary intensities")
        if isinstance(oe, dict):
            lists = [v for v in oe.values() if isinstance(v, list)]
            need(len(lists) == 1, f"output_edges_eV object without exactly one list: {sorted(oe)}")
            values, axis_names = lists[0], [oe.get("axis")]
        else:
            values, axis_names = oe, []
        need([float(x).hex() for x in values] == [x.hex() for x in edges], "output_edges_eV differs from the BE input edges")
        axis_names.append(p.get("axis"))
        need("binding_energy_eV" in axis_names and not any(a == "kinetic_energy_eV" for a in axis_names),
             f"the axis of the top-level input record is not explicit: {axis_names}")
    elif key.startswith("run.json: results.single_loss.primary and edges are KE"):
        single = G(G(manifest, "results"), "single_loss")
        p = G(single, "primary")
        c = ref["control"]["kinetic_primary"]
        position = "energies_eV" if c["kind"] == "delta" else "edges_eV"
        need([float(x).hex() for x in G(p, position)] == [rv(x).hex() for x in c["coordinates"]], "KE primary positions")
        need([float(x).hex() for x in G(p, "intensities")] == [rv(x).hex() for x in c["intensities"]], "KE primary intensities")
        need([float(x).hex() for x in G(single, "output_edges_eV")] == [x.hex() for x in ke_edges], "KE edges of single_loss")
    elif key == "run.json: spectrum_definition.axis":
        need(G(G(manifest, "spectrum_definition"), "axis") == "binding_energy_eV", "spectrum_definition.axis")
    elif key.startswith("run.json: spectrum_input.binding_output_edges"):
        got = G(G(manifest, "spectrum_input"), "binding_output_edges_sha256_float64_le")
        need(got == edge_hash(edges) == ref["hashes"]["binding_output_edges_sha256_float64_le"], f"hash {got!r}")
        need(got == edge_hash(G(G(manifest, "energy_axis"), "input_binding_edges_eV")), "hash of the recorded BE edges")
    elif key.startswith("run.json: spectrum_input.kinetic_output_edges"):
        got = G(G(manifest, "spectrum_input"), "kinetic_output_edges_sha256_float64_le")
        need(got == edge_hash(ke_edges) == ref["hashes"]["kinetic_output_edges_sha256_float64_le"], f"hash {got!r}")
        need(got == edge_hash(G(G(manifest, "energy_axis"), "kinetic_output_edges_eV")), "hash of the recorded KE edges")
    elif key == "run.json: spectrum_input.schema":
        need(G(G(manifest, "spectrum_input"), "schema") == MEASURED, "spectrum_input.schema")
        need(set(G(manifest, "spectrum_input")) == {"file", "sha256", "schema", "binding_output_edges_sha256_float64_le",
                                                    "kinetic_output_edges_sha256_float64_le"}, "spectrum_input keys")
    elif key.startswith("run.json: spectrum_input.sha"):
        need(G(G(manifest, "spectrum_input"), "sha256") == record["spectrum_sha256"], "spectrum_input.sha256")
        need(G(G(manifest, "spectrum_input"), "file") == record["spectrum_path"], "spectrum_input.file")
    elif key == "run.json: truncation":
        t = G(manifest, "truncation")
        need(G(t, "state") == e["truncation_state"] and G(t, "support_order") == e["support_order"]
             == ref["support"]["support_order"] and G(t, "computed") == ref["support"]["computed"],
             f"truncation {t.get('state')!r}, {t.get('support_order')!r}, {t.get('computed')!r}")
        got = ["complete" if is_complete((G(t, "neglected_at_or_above_region_bound")[i], G(t, "bound_reason")[i],
                                          G(t, "tau_per_eV")[i])) else "reached" for i in range(len(edges) - 1)]
        need(got == ref["support"]["bound_class"], f"bound classes {got} != R-C {ref['support']['bound_class']}")
        same(t, G(G(manifest, "results"), "truncation"), "truncation")
        need(G(manifest, "physical_validity") == "not_established"
             and all(G(G(manifest, "comparison"), g) == "not_matched" for g in G(manifest, "comparison")),
             "comparison / physical_validity")
    elif key == "run.json: units.binding_energy":
        need(G(G(manifest, "units"), "binding_energy") == "eV", "units.binding_energy")
    elif key == "run.json: warnings":
        records = [(G(w, "code"), G(w, "message")) for w in G(G(manifest, "diagnostics"), "warnings")]
        codes = [c for c, _m in records]
        need(all(c in codes for c in e["warnings_include"]) and not any(c in codes for c in e.get("warnings_exclude", ())),
             f"registered warnings; codes {codes}")
        need(list(G(manifest, "warnings")) == [m for _c, m in records], "warnings and diagnostics.warnings do not pair")
        need(("axis_conversion_support_disagreement" in codes) == ref["support"]["warning"], "disagreement warning against R-C")
        need(all(m.isascii() for _c, m in records), "a warning message is not ASCII")
    else:
        raise Unreadable(f"no CLI diagnostic check for {key!r}")


def check_capabilities(row, out, key):
    import xpsloss
    record = out.cli
    need(record is not None and record["exit_code"] == 0, "upstream: capabilities did not complete")
    caps = json.loads(record["stdout"])
    block = G(caps, "multiple_scattering")
    if key == "multiple_scattering.energy_axes":
        need(G(block, "energy_axes") == ["kinetic_energy_eV", "binding_energy_eV"], f"{block.get('energy_axes')!r}")
    elif key == "multiple_scattering.measured_spectrum_schema":
        need(G(block, "measured_spectrum_schema") == MEASURED, f"{block.get('measured_spectrum_schema')!r}")
    elif key.startswith("multiple_scattering.unsupported_physics =="):
        need(G(block, "unsupported_physics") == NINE_UNSUPPORTED, f"{block.get('unsupported_physics')!r}")
    else:
        box = baseline()
        if box["unavailable"]:
            status, detail = box["unavailable"]
            if status == "not_run":
                raise NotRun(detail)
            raise AssertionError(detail)
        old = json.loads(box["record"]["capabilities"]["stdout"])
        if key == "supported == baseline":
            need(G(caps, "supported") == old["supported"], "supported")
        elif key == "each side reports its own version":
            need(G(caps, "version") == xpsloss.__version__ and old["version"] == box["record"]["version"], "versions")
        elif key.startswith("baseline capabilities with one removal and two additions"):
            want = copy.deepcopy(old)
            want["version"] = caps.get("version")
            ms = want["multiple_scattering"]
            need("binding_energy_axis" in ms["unsupported_physics"], "baseline: binding_energy_axis is not in the baseline")
            ms["unsupported_physics"].remove("binding_energy_axis")
            need("energy_axes" not in ms and "measured_spectrum_schema" not in ms, "baseline: the additions exist already")
            mine = copy.deepcopy(caps)
            added = {k: mine["multiple_scattering"].pop(k, None) for k in ("energy_axes", "measured_spectrum_schema")}
            need(added == {"energy_axes": ["kinetic_energy_eV", "binding_energy_eV"], "measured_spectrum_schema": MEASURED},
                 f"additions {added}")
            d = first_difference(mine, want, "capabilities")
            need(d is None, d)
        else:
            raise Unreadable(f"no capabilities check for {key!r}")


# ---- replay (prereg 4.6)
def jsonable(x):
    """A result as the JSON a manifest would hold: floats by bits, tuples as arrays."""
    if isinstance(x, np.ndarray):
        return [jsonable(v) for v in x.tolist()]
    if isinstance(x, (float, np.floating)):
        return ["float", float(x).hex()]
    if isinstance(x, (bool, np.bool_)) or x is None or isinstance(x, str):
        return x if not isinstance(x, np.bool_) else bool(x)
    if isinstance(x, (int, np.integer)):
        return int(x)
    if isinstance(x, dict):
        return {str(k): jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [jsonable(v) for v in x]
    if dataclasses.is_dataclass(x):
        return {f.name: jsonable(getattr(x, f.name)) for f in dataclasses.fields(x)}
    raise Unreadable(f"not JSON-like: {type(x).__name__}")


def same_json(a, b, label):
    d = first_difference(jsonable(a), jsonable(b), label)
    need(d is None, d)


def check_replay(row, out, key):
    manifest = manifest_of(out)
    replay = replay_of(out)
    record = out.cli
    orders = [G(o, "order") for o in G(replay, "orders")]
    if "multiple_loss_spectrum.csv" in key:
        table_ = spectrum_csv(out)
        names = cli_region_names(row)
        need([r["region"] for r in table_["records"]] == names, "row names")
        expected_columns = 7 + 6*len(orders) + 2
        need(len(table_["header"]) == expected_columns, f"{len(table_['header'])} columns, expected {expected_columns}")
        for name, cells in zip(names, table_["records"]):
            region = ("bin", int(name)) if name.isdigit() else ("outside", name)

            def top(field):
                return (float(G(replay, field)[region[1]]) if region[0] == "bin"
                        else float(G(G(G(replay, "outside"), region[1]), field)))
            want = {"zero_loss": top("zero_loss"), "zero_loss_error": top("zero_loss_error"),
                    "spectrum": top("spectrum"), "spectrum_error": top("spectrum_error")}
            for n in orders:
                o = order_of(replay, n)
                for c in COMPONENTS:
                    for suffix in ("", "_error"):
                        holder = o if region[0] == "bin" else G(G(o, "outside"), region[1])
                        value = G(holder, c + suffix)
                        want[f"order{n}_{c}{suffix}"] = float(value[region[1]]) if region[0] == "bin" else float(value)
            need(set(want) | {"region", "binding_edge_low_eV", "binding_edge_high_eV", "kinetic_edge_low_eV",
                              "kinetic_edge_high_eV"} == set(cells), "the column map does not cover the header")
            for column, value in want.items():
                csv_float(cells[column], value, f"{name} {column}")
        return
    if "order_areas.csv" in key:
        header, lines = parse_csv(G(record["files"], "order_areas.csv"))
        need(header == ["order", "column", "area", "area_error", "cutoff_sensitivity", "cutoff_sensitivity_error"], "header")
        need([(int(line[0]), line[1]) for line in lines] == [(n, c) for n in orders for c in COLUMNS], "rows")
        for line in lines:
            o = order_of(replay, int(line[0]))
            values = [G(G(o, name), line[1]) for name in ("area", "area_error", "cutoff_sensitivity", "cutoff_sensitivity_error")]
            for text, value in zip(line[2:], values):
                csv_float(text, float(value), f"order {line[0]} {line[1]}")
        return
    results = G(manifest, "results")
    for name in ("zero_loss", "zero_loss_error", "spectrum", "spectrum_error", "outside", "truncation",
                 "order_area_ratio", "convention_ratio"):
        same_json(G(results, name), G(replay, name), f"results.{name}")
    need(set(G(results, "outside")) == {"low_binding", "high_binding"}, "results.outside names")
    need(len(G(results, "orders")) == len(orders), "number of orders")
    for o_json, o in zip(G(results, "orders"), G(replay, "orders")):
        need(set(G(o_json, "outside")) == {"low_binding", "high_binding"}, "orders outside names")
        for name in [c + s for c in COMPONENTS for s in ("", "_error", "_error_parts")] + [
                "outside", "area", "area_error", "area_error_parts", "area_zero_kind", "cutoff_sensitivity",
                "cutoff_sensitivity_error", "mesh", "order"]:
            same_json(G(o_json, name), G(o, name), f"results.orders[{G(o, 'order')}].{name}")
    need(list(G(manifest, "warnings")) == list(G(replay, "warnings")), "warnings")
    need([G(w, "code") for w in G(G(manifest, "diagnostics"), "warnings")] == warning_codes(replay), "warning identifiers")
    need(G(G(manifest, "truncation"), "state") == G(G(replay, "truncation"), "state"), "state")
    ea_json, ea = G(manifest, "energy_axis"), energy_axis_of(replay)
    for name in ("conversion_constant_eV", "conversion_constant_exact", "implied_offset_eV", "rounding_max_eV",
                 "primary_reference_kinetic_eV", "electron_energy_minus_reference_eV", "kinetic_output_edges_eV",
                 "input_binding_edges_eV", "support_order_binding", "support_order_kinetic", "region_loss_upper_eV",
                 "support_comparison", "record", "not_recorded", "direction"):
        same_json(G(ea_json, name), G(ea, name), f"energy_axis.{name}")
    need(G(G(manifest, "spectrum_input"), "binding_output_edges_sha256_float64_le")
         == edge_hash(G(ea, "input_binding_edges_eV")), "binding edge hash against the replay")
    need(G(G(manifest, "spectrum_input"), "kinetic_output_edges_sha256_float64_le")
         == edge_hash(G(ea, "kinetic_output_edges_eV")), "kinetic edge hash against the replay")


# ---- BASE rows
def check_base_record(row, out, key, oracle):
    rid = row["id"]
    current = base_current(rid)
    if oracle == "same row":
        explicit = base_current(rid, explicit_none=True)
        d = first_difference(explicit, current, "explicit energy_axis=None")
        need(d is None, d)
        return
    old = baseline_row(rid)
    if key == "every field of the result or exception":
        mine = dict(current)
        mine["result"] = without_axis_field(mine["result"])
        theirs = {k: v for k, v in old.items() if k != "seconds"}
        need(not has_key(mine, "energy_axis"), "an energy_axis key remains outside the registered place")
        d = first_difference(mine, theirs, "result")
        need(d is None, d)
    elif key.startswith("multiple_loss_spectrum.csv or single_loss_spectrum.csv bytes"):
        name = "multiple_loss_spectrum.csv" if "multiple_loss_spectrum.csv" in old["files"] else "single_loss_spectrum.csv"
        need(G(current["files"], name) == old["files"][name], f"{name} differs from the baseline")
    elif key == "second CSV bytes":
        names = sorted(n for n in old["files"] if n.endswith(".csv"))
        need(len(names) == 2, f"baseline CSV files {names}")
        second = [n for n in names if n not in ("multiple_loss_spectrum.csv", "single_loss_spectrum.csv")]
        need(len(second) == 1 and G(current["files"], second[0]) == old["files"][second[0]], f"{second} differs from the baseline")
        need(sorted(current["files"]) == sorted(old["files"]), "file names")
    elif key.startswith("run.json without time"):
        mine, theirs = (json.loads(x["files"]["run.json"]) for x in (current, old))
        need(not has_key(mine, "energy_axis"), "the KE run.json carries an energy_axis key")
        d = first_difference(plain(strip_manifest(mine)), plain(strip_manifest(theirs)), "run.json")
        need(d is None, d)
        need(current["stderr"] == old["stderr"], "stderr differs from the baseline")
    else:
        raise Unreadable(f"no BASE record check for {key!r}")


def strip_manifest(m):
    m = copy.deepcopy(m)
    for name in ("created_utc", "elapsed_seconds", "input_file", "package_version", "output_dir"):
        m.pop(name, None)
    if isinstance(m.get("spectrum_input"), dict):
        m["spectrum_input"].pop("file", None)
    return m


def check_base_diagnostic(row, out, key, oracle):
    rid = row["id"]
    if key == "result.energy_axis is None":
        current = base_current(rid)["result"]
        need(current.get("@type") == "MultipleLossSpectrumResult" and "energy_axis" in current
             and current["energy_axis"] is None, "the spectrum result has no field energy_axis equal to None")
        return
    if key == "energy_axis keyword raises TypeError":
        got = base_current(rid, explicit_none=True)
        need(got.get("raised") == "TypeError", f"energy_axis=None raised {got.get('raised')!r}")
        return
    import xpsloss
    current, old = base_current(rid), baseline_row(rid)
    mine, theirs = (json.loads(x["files"]["run.json"]) for x in (current, old))
    if key == "created_utc is valid":
        for m in (mine, theirs):
            datetime.datetime.fromisoformat(G(m, "created_utc"))
    elif key == "elapsed_seconds is valid":
        for m in (mine, theirs):
            need(isinstance(G(m, "elapsed_seconds"), float) and 0.0 <= m["elapsed_seconds"] < 3600.0, "elapsed_seconds")
    elif key == "paths are valid":
        for m, side in ((mine, current), (theirs, old)):
            need(G(m, "input_file") == side["config_path"] and G(G(m, "spectrum_input"), "file") == side["spectrum_path"],
                 "paths of the manifest")
    elif key == "each side reports its own version":
        need(G(mine, "package_version") == xpsloss.__version__
             and G(theirs, "package_version") == baseline()["record"]["version"], "package_version")
    else:
        raise Unreadable(f"no BASE diagnostic check for {key!r}")


def check_inequality(row, out, key):
    need(row["id"] == "BS6", row["id"])
    mine = triple(success(out), 0)[0]
    other = triple(control_result("BS5"), to_kinetic("BS5", ("bin", 0))[1])[0]
    need(isinstance(mine, float) and isinstance(other, float) and mine > other, f"{mine!r} > {other!r} does not hold")


def dispatch(row, out, key, kind, oracle, use):
    if kind == "state":
        return check_state(row, out, key, use)
    if kind == "exact_conversion":
        return check_conversion(row, out, key, oracle, use)
    if kind == "bin_membership":
        return check_membership(row, out, key, use)
    if kind == "bit_identity":
        return check_bit_identity(row, out, key, oracle)
    if kind == "record_identity":
        return check_record(row, out, key, oracle, use)
    if kind == "truncation":
        return check_truncation(row, out, key, oracle, use)
    if kind == "bound_bracket":
        return check_bound_bracket(row, out, key)
    if kind == "diagnostic":
        return check_diagnostic(row, out, key, oracle, use)
    if kind == "inequality":
        return check_inequality(row, out, key)
    raise Unreadable(f"unknown kind {kind}")


def judge_row(row, judgement):
    """Run a row and judge every ledger record of it; returns the outcome."""
    rid = row["id"]
    start = time.perf_counter()
    if row["group"] == "BASE":
        runs = _CACHE.setdefault("runs", {})
        if rid not in runs:
            out = Outcome(rid)
            runs[rid] = execute(out, lambda: run_base(row, out))
        out = runs[rid]
    else:
        out = run_row(rid)
    use = RefUse(rid)
    for key, kind, oracle in ledger()[rid]:
        status, detail = guard(lambda: dispatch(row, out, key, kind, oracle, use))
        judgement.add(rid, key, kind, oracle, status, detail)
    faults = use.faults()
    # a reference key is counted as unused only when its consumer could run at all
    if faults:
        judgement.faults[rid] = faults
    judgement.seconds[rid] = time.perf_counter() - start
    return out


if __name__ == "__main__":
    if len(sys.argv) == 4 and sys.argv[1] == "--baseline-child":
        _baseline_child(sys.argv[2], sys.argv[3])
    else:
        raise SystemExit("usage: m4b2_states.py --baseline-child OUT PAYLOAD")
