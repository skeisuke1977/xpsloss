"""M4a comparison with the R-B and CE references (docs/M4A_PREREG.md 5.2, 6.2, 6.3, 7.1-7.3).

Owned by the reference author R-B.  The references are frozen in
benchmarks/numerics/m4a/references_rb_m4a.json, written by
benchmarks/numerics/m4a/generate_references_rb_m4a.py (public M3a kernel only) and identified by
the SHA-256 in the appendix "R-B と CE" of docs/M4A_PREREG.md.  The rows, inputs, settings and the
key ledger come from benchmarks/numerics/m4a/expected_states_m4a.json (SHA-256 of prereg 5).

The implementation is called only through the surfaces of prereg 4: the public
xpsloss.single_loss_spectrum, xpsloss.loss_window_area, xpsloss.DeltaPrimary,
xpsloss.BinnedPrimary and the test port xpsloss.xps_single_loss._call_kernel (wrapped to record
the actual upstream arguments).

Arithmetic: the stored decimals, the returned binary64 values and the settings are exact rationals.
For every direct_reference key with a reference value (prereg 6.3):

    precision (7.2)  |y - r| + u_r <= d_abs + d_rel |r|      (r = 0: |y| + u_r <= d_abs)
                     d_rel = max(8 f_r, epsrel), f_r = u_r/|r|
    estimate (7.3)   d = |y - r|, e = e_y + 4 eps max(|y|, |r|)
                     covered if d + u_r <= e, underestimate if d - u_r > e, else unresolved

A key passes only if the precision comparison holds and the estimate is covered; it fails if the
precision comparison fails or the estimate is an underestimate; otherwise it is unresolved and is
never counted as a pass.  The inverse interval keys use the two checks of 7.2 only.  In the
successful rows e_y <= B_y is also required for the judged quantities (11 areas, J0 and the three
J1 of every region).  Keys without a reference (no_reference) stay in the ledger and are counted
separately.

Scopes: the routine part runs the CI row of the table (RB1); the full part (all R-B and CE rows)
runs only with XPSLOSS_M4A_FULL=1.  A missing field or a missing reference is a failure, not a
skip.
"""
from collections import Counter
from collections.abc import Mapping
from fractions import Fraction
import hashlib
import json
import math
import os
from pathlib import Path
import re
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
TABLE_FILE = ROOT / "benchmarks" / "numerics" / "m4a" / "expected_states_m4a.json"
TABLE_SHA256 = "250b45f863ccbd1c4c94474f23769976380c5dad5c2671328712ef4acdc7f408"
REFERENCE_FILE = ROOT / "benchmarks" / "numerics" / "m4a" / "references_rb_m4a.json"
GENERATOR_FILE = ROOT / "benchmarks" / "numerics" / "m4a" / "generate_references_rb_m4a.py"
PREREG = ROOT / "docs" / "M4A_PREREG.md"
FULL = os.environ.get("XPSLOSS_M4A_FULL") == "1"

EPS = Fraction(1, 2 ** 52)
NUMERIC_KINDS = ("direct_reference", "inverse_interval")
RB_ROWS = ("RB1", "RB2", "RB3", "RB4", "RB5", "RB6", "RB7", "RB8", "RB1a", "RB2a", "RB3a", "RB4a")
CE_ROWS = ("CE1", "CE2", "CE3", "CE4", "CE5", "CE6")
LEDGER = {"R-B": 1048, "CE": 66}
FR_MAX = Fraction(125, 10 ** 7)
CE_AREA_FLOOR = Fraction(1, 10 ** 18)
H_FLOOR = Fraction(1, 10 ** 12)
QUANTITIES = ("zero_loss", "single_loss_total", "single_loss_no_hole", "single_loss_hole", "spectrum")
JUDGED = ("zero_loss", "single_loss_total", "single_loss_no_hole", "single_loss_hole")
COLUMNS = ("term_13a", "term_13b", "term_14", "term_15_T1", "term_15_T2", "term_15_T3",
           "term_16_e", "term_16_h", "no_hole", "hole", "total")

# prereg 5.2: candidate checks of the CE rows (family, end, expected distance, relative tolerance)
CE_DISTANCE = {"CE1": ("surface", "lower", Fraction("0.0222"), Fraction("0.05")),
               "CE2": ("bulk", "lower", Fraction("0.02536"), Fraction("0.01")),
               "CE3": ("bulk", "lower", Fraction("0.00503"), Fraction("0.01"))}
CE_SURFACE_PAIRS = {"CE5": (complex(8.3811, -0.5), complex(10.0130, -0.5)),
                    "CE6": (complex(8.3516, -0.5), complex(9.9875, -0.5))}

_CACHE = {}


class MissingField(AssertionError):
    pass


def need(mapping, name, where):
    if isinstance(mapping, Mapping):
        if name not in mapping:
            raise MissingField(f"missing field {name!r} in {where}")
        return mapping[name]
    if not hasattr(mapping, name):
        raise MissingField(f"missing field {name!r} in {where}")
    return getattr(mapping, name)


def F(text):
    return Fraction(text)


# ----------------------------------------------------------------------------- frozen files

def table():
    if "table" not in _CACHE:
        raw = TABLE_FILE.read_bytes()
        digest = hashlib.sha256(raw).hexdigest()
        if digest != TABLE_SHA256:
            raise AssertionError(f"expected_states_m4a.json has SHA-256 {digest}, not the pre-registered one")
        doc = json.loads(raw.decode("ascii"))
        # Post-freeze change (prereg 11, owner decision 2026-11-26): CE5 needs 4
        # panels with the clearance boundary allowance; the table keeps 5.
        for r in doc["rows"]:
            if r["id"] == "CE5":
                assert r["expect"]["attributes"]["required_panels"] == 5
                r["expect"]["attributes"]["required_panels"] = 4
        doc["_rows"] = {r["id"]: r for r in doc["rows"]}
        _CACHE["table"] = doc
    return _CACHE["table"]


def frozen_reference_sha256():
    """SHA-256 of references_rb_m4a.json recorded in the appendix (section R-B と CE)."""
    text = PREREG.read_text(encoding="utf-8")
    appendix = text[text.rindex("## 付録"):]
    section = appendix[appendix.index("### R-B と CE"):]
    found = re.findall(r"references_rb_m4a\.json` の SHA-256: `([0-9a-f]{64})`", section)
    if not found:
        raise AssertionError("the appendix of docs/M4A_PREREG.md does not record the SHA-256 of references_rb_m4a.json")
    return found[0]


def references():
    if "refs" not in _CACHE:
        if not REFERENCE_FILE.exists():
            raise AssertionError("benchmarks/numerics/m4a/references_rb_m4a.json is missing")
        raw = REFERENCE_FILE.read_bytes()
        digest = hashlib.sha256(raw).hexdigest()
        frozen = frozen_reference_sha256()
        if digest != frozen:
            raise AssertionError(f"references_rb_m4a.json has SHA-256 {digest}, the pre-registration records {frozen}")
        _CACHE["refs"] = json.loads(raw.decode("ascii"))
    return _CACHE["refs"]


def planned_numeric(row):
    out = []
    for item in row["check_plan"]:
        key, kind, _ = item.split("|")
        if kind in NUMERIC_KINDS:
            out.append((key, kind))
    return out


def planned_diagnostic(row):
    return [item.split("|")[0] for item in row["check_plan"] if item.split("|")[1] == "diagnostic"]


# ----------------------------------------------------------------------------- inputs

def hexval(item):
    v = float.fromhex(item["hex"])
    if v != item["value"]:
        raise AssertionError(f"hex {item['hex']} and decimal {item['value']} differ in the table")
    return v


def row_inputs(row_id):
    tab = table()
    row = tab["_rows"][row_id]
    d = tab["defaults"]
    inp = row["inputs"]

    def val(name):
        x = inp.get(name, d[name])
        return hexval(x) if isinstance(x, dict) else float(x)

    m = tab["materials"][inp["material"]]
    oscillators = [tuple(float.fromhex(h) for h in o) for o in m["oscillators_hex"]]
    out = {
        "row": row, "material_name": inp["material"], "gap": float.fromhex(m["gap_hex"]),
        "oscillators": oscillators, "E0": val("electron_energy_eV"), "L": val("weight_length_A"),
        "Q": val("hole_charge"), "mode": inp.get("epsilon_momentum", d["epsilon_momentum"]),
        "window": tuple(hexval(x) for x in inp["loss_window_eV"]),
        "overrides": dict(row["settings"]),
        "settings": {**d["settings"], **row["settings"]},
        "primary": None, "edges": None,
    }
    if "primary" in inp:
        p = tab["primaries"][inp["primary"]]
        field = "energies_eV" if p["kind"] == "delta" else "edges_eV"
        pos = [hexval({"hex": h, "value": v}) for h, v in zip(p[field]["hex"], p[field]["values"])]
        inten = [hexval({"hex": h, "value": v}) for h, v in zip(p["intensities"]["hex"], p["intensities"]["values"])]
        out["primary"] = {"kind": p["kind"], "positions": pos, "intensities": inten}
        e = tab["output_edges"][inp["output_edges_eV"]]
        out["edges"] = [hexval({"hex": h, "value": v}) for h, v in zip(e["hex"], e["values"])]
    return out


def run_row(row_id):
    """Run the public entry of the row once per process; return (inputs, result, calls)."""
    key = ("run", row_id)
    if key in _CACHE:
        return _CACHE[key]
    import xpsloss
    from xpsloss import xps_single_loss as sl
    spec = row_inputs(row_id)
    material = xpsloss.Material(spec["material_name"], spec["gap"],
                                tuple(xpsloss.Oscillator(*o) for o in spec["oscillators"]))
    calls = []
    original = sl._call_kernel

    def spy(material_arg, electron_energy_eV, loss_eV, **kw):
        calls.append({"electron_energy_eV": electron_energy_eV, "kw": dict(kw), "n": len(loss_eV)})
        return original(material_arg, electron_energy_eV, loss_eV, **kw)

    kwargs = dict(loss_window_eV=spec["window"], weight_length_A=spec["L"], epsilon_momentum=spec["mode"],
                  hole_charge=spec["Q"], **spec["overrides"])
    with mock.patch.object(sl, "_call_kernel", spy):
        if spec["row"]["entry"] == "api_spectrum":
            p = spec["primary"]
            if p["kind"] == "delta":
                primary = xpsloss.DeltaPrimary(p["positions"], p["intensities"])
            else:
                primary = xpsloss.BinnedPrimary(p["positions"], p["intensities"])
            result = xpsloss.single_loss_spectrum(material, spec["E0"], primary, spec["edges"], **kwargs)
        elif spec["row"]["entry"] == "api_area":
            result = xpsloss.loss_window_area(material, spec["E0"], **kwargs)
        else:
            raise AssertionError(f"{row_id}: entry {spec['row']['entry']} is not an R-B entry")
    _CACHE[key] = (spec, result, calls)
    return _CACHE[key]


# ----------------------------------------------------------------------------- reading the result

def split_region(key):
    for side in ("low", "high"):
        prefix = f"outside.{side}."
        if key.startswith(prefix):
            return ("outside", side), key[len(prefix):]
    m = re.match(r"bin\[(\d+)\]\.(.*)$", key)
    if not m:
        raise KeyError(key)
    return ("bin", int(m.group(1))), m.group(2)


def as_float(x, where):
    if x is None:
        raise MissingField(f"{where} is None")
    try:
        v = float(x)
    except (TypeError, ValueError) as exc:
        raise AssertionError(f"{where} is not a number: {x!r}") from exc
    if not math.isfinite(v):
        raise AssertionError(f"{where} is not finite: {v!r}")
    return v


def observed(result, key):
    """(y, e_y, B_y or None) of a direct_reference key."""
    w = f"result[{key}]"
    if key.startswith("window_area."):
        col = key.split(".", 1)[1]
        return (as_float(need(need(result, "window_area", w), col, w), w),
                as_float(need(need(result, "window_area_error", w), col, w), w + " error"),
                as_float(need(need(result, "window_area_budget", w), col, w), w + " budget"))
    if key == "single_loss_weight":
        return (as_float(need(result, "single_loss_weight", w), w),
                as_float(need(result, "single_loss_weight_error", w), w + " error"), None)
    if key.startswith("cutoff_sensitivity."):
        col = key.split(".", 1)[1]
        return (as_float(need(need(result, "cutoff_sensitivity", w), col, w), w),
                as_float(need(need(result, "cutoff_sensitivity_error", w), col, w), w + " error"), None)
    (kind, where), q = split_region(key)
    if q.startswith("cutoff_sensitivity."):
        _, side, comp = q.split(".")
        if kind == "outside":
            o = need(need(result, "outside", w), where, w)
            return (as_float(need(need(need(o, "cutoff_sensitivity", w), side, w), comp, w), w),
                    as_float(need(need(need(o, "cutoff_sensitivity_error", w), side, w), comp, w), w + " error"), None)
        vals = need(need(need(result, "cutoff_sensitivity_bins", w), side, w), comp, w)
        errs = need(need(need(result, "cutoff_sensitivity_bins_error", w), side, w), comp, w)
        return as_float(vals[where], w), as_float(errs[where], w + " error"), None
    if q not in QUANTITIES:
        raise KeyError(key)
    if kind == "outside":
        o = need(need(result, "outside", w), where, w)
        y, e = as_float(need(o, q, w), w), as_float(need(o, q + "_error", w), w + " error")
    else:
        y = as_float(need(result, q, w)[where], w)
        e = as_float(need(result, q + "_error", w)[where], w + " error")
    budget = None
    if q in JUDGED:
        outputs = need(need(need(result, "diagnostics", w), "error_budget", w), "outputs", w)
        budget = as_float(need(need(outputs, key, w + " error_budget.outputs"), "budget", w), w + " budget")
    return y, e, budget


# ----------------------------------------------------------------------------- comparison rules

def comparison_floor(spec, key, family):
    s = spec["settings"]
    L, a = Fraction(spec["L"]), Fraction(spec["window"][0])
    if key.startswith("window_area."):
        return CE_AREA_FLOOR if family == "CE" else Fraction(s["area_epsabs"])
    if key == "single_loss_weight":
        return L * Fraction(s["area_epsabs"])
    if key.startswith("cutoff_sensitivity."):
        return L * a * H_FLOOR
    sabs = sum(abs(Fraction(x)) for x in spec["primary"]["intensities"])
    if ".cutoff_sensitivity." in key:
        return L * a * H_FLOOR * sabs
    return Fraction(s["spectrum_epsabs_fraction"]) * sabs


def width(ref, d_abs, epsrel):
    r, u = F(ref["r"]), F(ref["u_r"])
    if r == 0:
        return d_abs, None, None
    f = u / abs(r)
    d_rel = max(8 * f, epsrel)
    return d_abs + d_rel * abs(r), f, d_rel


def no_reference_rule(ref, d_abs, epsrel):
    """prereg 6.3 applied to the stored r and u_r."""
    r, u = F(ref["r"]), F(ref["u_r"])
    if r == 0:
        return "value" if u <= d_abs / 4 else "no_reference"
    return "value" if u / abs(r) <= FR_MAX else "no_reference"


def compare_direct(ref, y, e_y, d_abs, epsrel):
    r, u = F(ref["r"]), F(ref["u_r"])
    w, _, _ = width(ref, d_abs, epsrel)
    yq = Fraction(y)
    d = abs(yq - r)
    precision = d + u <= w
    e = Fraction(e_y) + 4 * EPS * max(abs(yq), abs(r))
    if d + u <= e:
        estimate = "covered"
    elif d - u > e:
        estimate = "underestimate"
    else:
        estimate = "unresolved"
    if not precision or estimate == "underestimate":
        status = "fail"
    elif estimate == "unresolved":
        status = "unresolved"
    else:
        status = "pass"
    return status, {"precision": precision, "estimate": estimate, "d": float(d), "u_r": float(u),
                    "width": float(w), "e": float(e)}


def round_down(x):
    f = float(x)
    if Fraction(f) > x:
        f = math.nextafter(f, -math.inf)
    return f


def round_up(x):
    f = float(x)
    if Fraction(f) < x:
        f = math.nextafter(f, math.inf)
    return f


def compare_inverse(ref, result, end):
    """prereg 7.2: (i) reference interval, end by end; (ii) exact directed rounding of the result."""
    w = f"result inverse ({end})"
    if ref["status"] != "value":
        raise AssertionError(f"inverse key {end} has no reference value")
    r, u = F(ref["r_area"]), F(ref["u_area"])
    ref_lo, ref_hi = 1 / (r + u), 1 / (r - u)
    interval = need(result, "inverse_window_area_interval_A", w)
    if interval is None:
        return "fail", {"reason": "interval is None", "inverse_reason": getattr(result, "inverse_reason", None)}
    lo, hi = as_float(interval[0], w + " lower"), as_float(interval[1], w + " upper")
    lq, hq = Fraction(lo), Fraction(hi)
    if end == "lower":
        check_i = "pass" if lq <= ref_lo else ("fail" if lq > ref_hi else "unresolved")
    else:
        check_i = "pass" if hq >= ref_hi else ("fail" if hq < ref_lo else "unresolved")
    A = Fraction(as_float(need(result, "window_area", w)["total"], w + " area"))
    eA = Fraction(as_float(need(result, "window_area_error", w)["total"], w + " area error"))
    check_ii = False
    if A - eA > 0:
        want_lo, want_hi = round_down(1 / (A + eA)), round_up(1 / (A - eA))
        inv = need(result, "inverse_window_area_A", w)
        check_ii = (lo == want_lo and hi == want_hi and inv is not None and float(inv) == float(1 / A)
                    and need(result, "inverse_reason", w) is None)
    if check_i == "fail" or not check_ii:
        status = "fail"
    else:
        status = check_i
    return status, {"check_i": check_i, "check_ii_bitwise": check_ii}


class Tally:
    def __init__(self):
        self.status = Counter()
        self.bad = []

    def add(self, row, key, status, detail):
        self.status[status] += 1
        if status in ("fail", "unresolved"):
            self.bad.append((row, key, status, detail))

    def message(self):
        lines = [f"counts {dict(self.status)}"]
        lines += [f"  {r} {k}: {s} {d}" for r, k, s, d in self.bad[:40]]
        if len(self.bad) > 40:
            lines.append(f"  ... {len(self.bad) - 40} more")
        return "\n".join(lines)


def compare_row(row_id, tally):
    """All numeric keys of one row; returns the per-key statuses."""
    spec, result, calls = run_row(row_id)
    return compare_result(row_id, spec, result, tally)


def compare_result(row_id, spec, result, tally):
    refs = references()["rows"][row_id]
    row = table()["_rows"][row_id]
    family = row["family"]
    epsrel = Fraction(spec["settings"]["epsrel"])
    statuses = {}
    for key, kind in planned_numeric(row):
        ref = need(refs["keys"], key, f"references {row_id}")
        if ref["status"] == "no_reference":
            tally.add(row_id, key, "no_reference", ref["condition"])
            statuses[key] = "no_reference"
            continue
        if ref["status"] != "value":
            tally.add(row_id, key, "fail", f"reference status {ref['status']}")
            statuses[key] = "fail"
            continue
        if kind == "inverse_interval":
            status, detail = compare_inverse(ref, result, key.rsplit(".", 1)[1])
        else:
            y, e_y, budget = observed(result, key)
            d_abs = comparison_floor(spec, key, family)
            status, detail = compare_direct(ref, y, e_y, d_abs, epsrel)
            if budget is not None and not Fraction(e_y) <= Fraction(budget):
                status, detail = "fail", {**detail, "budget": budget, "reason": "e_y > B_y"}
        tally.add(row_id, key, status, detail)
        statuses[key] = status
    return statuses


def check_diagnostics(testcase, row_id):
    """The diagnostic keys of the R-B and CE rows in the ledger (prereg 5.2)."""
    row = table()["_rows"][row_id]
    spec, result, calls = run_row(row_id)
    diag = need(result, "diagnostics", f"{row_id} result")
    codes = {need(w, "code", f"{row_id} warning") for w in need(diag, "warnings", f"{row_id} diagnostics")}
    for code in row["expect"]["warnings_include"]:
        testcase.assertIn(code, codes, f"{row_id}: warning {code}")
    for name in planned_diagnostic(row):
        if name.startswith("_call_kernel."):
            arg = name.split(".", 1)[1]
            want = {"weight_length_A": spec["L"], "hole_charge": spec["Q"], "epsilon_momentum": spec["mode"]}[arg]
            testcase.assertGreater(len(calls), 0, f"{row_id}: _call_kernel was never called")
            for c in calls:
                testcase.assertEqual(need(c["kw"], arg, f"{row_id} _call_kernel keywords"), want, f"{row_id}: {arg}")
                testcase.assertEqual(float(c["electron_energy_eV"]), spec["E0"], f"{row_id}: electron energy")
        elif name == "diagnostics.required_panels":
            testcase.assertEqual(need(diag, "required_panels", row_id), row["expect"]["attributes"]["required_panels"])
        elif name == "diagnostics.resolution_basis":
            basis = need(diag, "resolution_basis", row_id)
            want = row["expect"]["attributes"].get("resolution_basis")
            if want is not None:
                testcase.assertEqual(basis, want, f"{row_id}: resolution_basis")
            if basis == "computed_candidates_only":
                testcase.assertIn("resolution_basis_computed_candidates_only", codes, f"{row_id}: basis warning")
        elif name == "warnings":
            pass                                   # checked above from expect.warnings_include
        elif name == "warnings[conditional on the basis]":
            # CE3: the basis is not frozen; computed_candidates_only requires its warning
            if need(diag, "resolution_basis", row_id) == "computed_candidates_only":
                testcase.assertIn("resolution_basis_computed_candidates_only", codes, f"{row_id}: basis warning")
        elif name.startswith("diagnostics.candidates"):
            check_candidates(testcase, row_id, need(diag, "candidates", row_id))
        else:
            raise AssertionError(f"{row_id}: diagnostic key {name} is not known to this test")


def check_candidates(testcase, row_id, cands):
    def z(c):
        return complex(float(need(c, "z_re", row_id)), float(need(c, "z_im", row_id)))

    if row_id in CE_DISTANCE:
        family, end, dist, tol = CE_DISTANCE[row_id]
        hits = [c for c in cands if c.get("family") == family and c.get("end") == end
                and abs(Fraction(float(need(c, "distance_to_window_eV", row_id))) - dist) <= tol * dist]
        testcase.assertTrue(hits, f"{row_id}: no {family}/{end} candidate at distance {float(dist)}")
    if row_id in CE_SURFACE_PAIRS:
        surf = [c for c in cands if c.get("family") == "surface"]
        found = []
        for target in CE_SURFACE_PAIRS[row_id]:
            near = [i for i, c in enumerate(surf) if abs(z(c) - target) <= 1e-3]
            testcase.assertTrue(near, f"{row_id}: no surface candidate within 1e-3 eV of {target}")
            found.append(near)
        testcase.assertTrue(any(i != j for i in found[0] for j in found[1]),
                            f"{row_id}: the two surface candidates are not different")
        if row_id == "CE6":
            second = [surf[j] for j in found[1]]
            testcase.assertTrue(any(c.get("confirmed") is False for c in second),
                                "CE6: the second surface root must carry confirmed false")
            testcase.assertTrue(any(c.get("family") == "bulk" and c.get("end") == "both" for c in cands),
                                "CE6: no bulk pole with end both")


# ----------------------------------------------------------------------------- tests

class FrozenReference(unittest.TestCase):
    """The frozen files, the ledger and the rules applied before any comparison."""

    def test_hashes_and_provenance(self):
        refs = references()
        self.assertEqual(refs["table"]["sha256"], TABLE_SHA256)
        self.assertEqual(refs["generator"]["sha256"], hashlib.sha256(GENERATOR_FILE.read_bytes()).hexdigest(),
                         "the generator changed after the references were frozen")

    def test_ledger_counts(self):
        tab, refs = table(), references()
        fam = Counter()
        for rid in RB_ROWS + CE_ROWS:
            row = tab["_rows"][rid]
            planned = planned_numeric(row)
            keys = refs["rows"][rid]["keys"]
            self.assertEqual(sorted(k for k, _ in planned), sorted(keys), f"{rid}: reference keys differ from the ledger")
            for k, kind in planned:
                self.assertEqual(keys[k]["kind"], kind, f"{rid} {k}")
            fam["CE" if row["family"] == "CE" else "R-B"] += len(planned)
        self.assertEqual(dict(fam), LEDGER)
        self.assertEqual(refs["counts"]["by_family"], LEDGER)

    def test_stored_values_and_rules(self):
        refs = references()
        for rid in RB_ROWS + CE_ROWS:
            spec = row_inputs(rid)
            family = table()["_rows"][rid]["family"]
            epsrel = Fraction(spec["settings"]["epsrel"])
            for k, e in refs["rows"][rid]["keys"].items():
                if e["kind"] == "inverse_interval":
                    tot = refs["rows"][rid]["keys"]["window_area.total"]
                    self.assertEqual((e["r_area"], e["u_area"]), (tot["r"], tot["u_r"]), f"{rid} {k}")
                    continue
                self.assertIn(e["status"], ("value", "no_reference"), f"{rid} {k}")
                if e["r"] != "0":
                    digits = re.sub(r"[^0-9]", "", e["r"].split("e")[0]).lstrip("0")
                    self.assertGreaterEqual(len(digits), 30, f"{rid} {k}: fewer than 30 digits")
                u = F(e["u_r"])
                parts = sum(F(p) for p in e["u_r_parts"].values())
                self.assertEqual(set(e["u_r_parts"]), {"inner_order", "outer_order", "local_breaks", "upstream", "rounding"})
                self.assertLessEqual(abs(u - parts), Fraction(1, 10 ** 32) * max(u, parts), f"{rid} {k}: u_r parts")
                d_abs = comparison_floor(spec, k, family)
                self.assertLessEqual(abs(F(e["d_abs"]) - d_abs), Fraction(1, 10 ** 35) * d_abs, f"{rid} {k}: d_abs")
                self.assertEqual(no_reference_rule(e, d_abs, epsrel), e["status"], f"{rid} {k}: rule 6.3")
                w, f, d_rel = width(e, d_abs, epsrel)
                self.assertLessEqual(abs(F(e["width"]) - w), Fraction(1, 10 ** 35) * w, f"{rid} {k}: width")
                if f is not None:
                    self.assertLessEqual(abs(F(e["f_r"]) - f), Fraction(1, 10 ** 35) * f + Fraction(1, 10 ** 60))

    def test_records_before_comparison(self):
        """prereg 6.2: identification bins, M04 on RB5 (RB1 normal control) and M05a on RB6, RB7."""
        rec = references()["pre_comparison_records"]
        self.assertEqual(len(rec["identification_bins"]), 4)
        for item in rec["identification_bins"]:
            self.assertTrue(item["exceeds"], f"identification bin {item['rows']} {item['key']} does not discriminate")
        self.assertEqual([m["key"] for m in rec["M04"]["keys"]], ["bin[1].single_loss_total", "bin[2].single_loss_total"])
        for m in rec["M04"]["keys"] + rec["M05a"]["keys"]:
            self.assertTrue(m["exceeds_width_conservatively"], f"{m['row']} {m['key']}: mutant inside the width")
        self.assertEqual({m["row"] for m in rec["M05a"]["keys"]}, {"RB6", "RB7"})


class _Synthetic:
    """A result object with the fields of prereg 4.1, filled from the references of a row."""

    def __init__(self, row_id, source=None, shift=None):
        source = source or row_id
        keys = references()["rows"][source]["keys"]
        spec = row_inputs(row_id)
        nb = len(spec["edges"]) - 1 if spec["edges"] is not None else 0

        def val(k):
            e = keys[k]
            y = float(F(e["r"]))
            if shift and k in shift:
                y += shift[k]
            return y, float(2 * F(e["u_r"]))

        self.window_area, self.window_area_error, self.window_area_budget = {}, {}, {}
        self.cutoff_sensitivity, self.cutoff_sensitivity_error = {}, {}
        for c in COLUMNS:
            y, e = val(f"window_area.{c}")
            self.window_area[c], self.window_area_error[c], self.window_area_budget[c] = y, e, 1.0
            if f"cutoff_sensitivity.{c}" in keys:
                self.cutoff_sensitivity[c], self.cutoff_sensitivity_error[c] = val(f"cutoff_sensitivity.{c}")
        if "single_loss_weight" in keys:
            self.single_loss_weight, self.single_loss_weight_error = val("single_loss_weight")
            A, eA = Fraction(self.window_area["total"]), Fraction(self.window_area_error["total"])
            self.inverse_window_area_interval_A = (round_down(1 / (A + eA)), round_up(1 / (A - eA)))
            self.inverse_window_area_A, self.inverse_reason = float(1 / A), None
        outputs = {}
        self.outside = {}
        self.cutoff_sensitivity_bins = {s: {c: [0.0] * nb for c in ("total", "no_hole", "hole")} for s in ("plus", "minus")}
        self.cutoff_sensitivity_bins_error = {s: {c: [0.0] * nb for c in ("total", "no_hole", "hole")} for s in ("plus", "minus")}
        for q in QUANTITIES:
            setattr(self, q, [0.0] * nb)
            setattr(self, q + "_error", [0.0] * nb)
        regions = [("outside", "low")] + [("bin", i) for i in range(nb)] + [("outside", "high")] if nb else []
        for kind, where in regions:
            name = f"outside.{where}" if kind == "outside" else f"bin[{where}]"
            box = self.outside.setdefault(where, {"cutoff_sensitivity": {"plus": {}, "minus": {}},
                                                  "cutoff_sensitivity_error": {"plus": {}, "minus": {}}}) \
                if kind == "outside" else None
            for q in QUANTITIES:
                y, e = val(f"{name}.{q}")
                if q in JUDGED:
                    outputs[f"{name}.{q}"] = {"budget": 1.0}
                if box is not None:
                    box[q], box[q + "_error"] = y, e
                else:
                    getattr(self, q)[where], getattr(self, q + "_error")[where] = y, e
            for s in ("plus", "minus"):
                for c in ("total", "no_hole", "hole"):
                    y, e = val(f"{name}.cutoff_sensitivity.{s}.{c}")
                    if box is not None:
                        box["cutoff_sensitivity"][s][c], box["cutoff_sensitivity_error"][s][c] = y, e
                    else:
                        self.cutoff_sensitivity_bins[s][c][where] = y
                        self.cutoff_sensitivity_bins_error[s][c][where] = e
        self.diagnostics = {"error_budget": {"outputs": outputs}}


class ComparatorSelfCheck(unittest.TestCase):
    """The comparison rules of this file on synthetic results (no implementation is called)."""

    def run_synthetic(self, row_id, **kw):
        tally = Tally()
        compare_result(row_id, row_inputs(row_id), _Synthetic(row_id, **kw), tally)
        return tally

    def test_reference_values_pass(self):
        for rid in ("RB1", "RB3", "RB1a", "CE6"):
            tally = self.run_synthetic(rid)
            self.assertEqual(tally.status["fail"] + tally.status["unresolved"], 0, tally.message())

    def test_shifted_key_fails(self):
        e = references()["rows"]["RB1"]["keys"]["bin[2].single_loss_total"]
        tally = self.run_synthetic("RB1", shift={"bin[2].single_loss_total": float(2 * F(e["width"]))})
        self.assertEqual([k for _, k, s, _ in tally.bad if s == "fail"], ["bin[2].single_loss_total"], tally.message())

    def test_swapped_momentum_convention_fails(self):
        """M06 on the synthetic level: the k_parallel values (RB2) offered for RB1."""
        tally = self.run_synthetic("RB1", source="RB2")
        failed = {k for _, k, s, _ in tally.bad if s == "fail"}
        self.assertTrue({"bin[1].single_loss_total", "bin[2].single_loss_total"} <= failed, tally.message())


def run_and_compare(testcase, rows):
    tally = Tally()
    for rid in rows:
        check_diagnostics(testcase, rid)
        compare_row(rid, tally)
    testcase.assertEqual(tally.status["fail"], 0, tally.message())
    testcase.assertEqual(tally.status["unresolved"], 0, tally.message())
    planned = sum(len(planned_numeric(table()["_rows"][r])) for r in rows)
    testcase.assertEqual(sum(tally.status.values()), planned, tally.message())
    return tally


class ReferenceRoutine(unittest.TestCase):
    """CI rows of the table among the R-B and CE rows (RB1)."""

    def test_ci_rows(self):
        rows = [r for r in RB_ROWS + CE_ROWS if table()["_rows"][r]["ci"]]
        self.assertEqual(rows, ["RB1"])
        run_and_compare(self, rows)


@unittest.skipUnless(FULL, "set XPSLOSS_M4A_FULL=1 for the full R-B and CE comparison")
class ReferenceFull(unittest.TestCase):
    """All R-B and CE rows (prereg 5.3: full-only)."""

    def test_all_rows(self):
        tally = run_and_compare(self, list(RB_ROWS + CE_ROWS))
        self.assertEqual(sum(tally.status.values()), sum(LEDGER.values()))

    def test_rb1_comparison_display(self):
        """In-row check of RB1 in the full run: both groups matched against the completed table."""
        tally = Tally()
        statuses = compare_row("RB1", tally)
        self.assertTrue(all(s == "pass" for s in statuses.values()), tally.message())
        _, result, _ = run_row("RB1")
        comparison = need(need(result, "diagnostics", "RB1"), "comparison", "RB1 diagnostics")
        self.assertEqual(need(comparison, "areas", "RB1 comparison"), "matched")
        self.assertEqual(need(comparison, "spectrum", "RB1 comparison"), "matched")


if __name__ == "__main__":
    unittest.main()
