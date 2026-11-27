"""M4b-1 comparison with the R-B references (docs/M4B1_PREREG.md 5.2, 5.3, 6.2, 6.3, 7.1-7.3).

Owned by the reference author R-B.  The references are frozen in
benchmarks/numerics/m4b1/references_rb_m4b1.json, written by
benchmarks/numerics/m4b1/generate_references_rb_m4b1.py (public M3a point values only) and
identified by the SHA-256 in the appendix of docs/M4B1_PREREG.md (R-B entry).  The rows, inputs,
settings and the key ledger come from benchmarks/numerics/m4b1/expected_states_m4b1.json.

The implementation is called only through the public surface of prereg 4.1 and 4.2:
xpsloss.multiple_loss_spectrum, xpsloss.multiple_loss_order_areas, xpsloss.DeltaPrimary,
xpsloss.BinnedPrimary, xpsloss.Material, xpsloss.Oscillator and the fields of the results.

Arithmetic: the stored decimals, the returned binary64 values and the settings are exact rationals.
For every direct_reference key with a reference value (prereg 6.3):

    precision (7.2)  |y - r| + u_r <= d_abs + d_rel |r|      (r = 0: |y| + u_r <= d_abs)
                     d_rel = min(max(8 f_r, n epsrel + (n-1) order_epsrel), 1e-4), stored per key
    estimate (7.3)   d = |y - r|, e = e_y + 4 eps max(|y|, |r|)
                     covered if d + u_r <= e, underestimate if d - u_r > e, else unresolved

A key passes only if the precision comparison holds and the estimate is covered; it fails if the
precision comparison fails or the estimate is an underestimate; otherwise it is unresolved and is
never counted as a pass.  For orders n >= 2 of bins, outside regions and areas the returned
`local` part must be within its budget (iii) and the `inherited` part within its budget (i)
(prereg 7.3, design 6.2).  Budget (i) is rebuilt with the reference norms; an excess on the
absolute side of less than 1e-3 of the budget is unresolved, not a failure, because the
implementation may hold slightly different norm estimates (decision of the reference author).
The key order_area_ratio (interval_contains, 7.2): the returned interval must contain the reference
interval (pass), is a failure if disjoint, unresolved otherwise; both ends and the value must also
equal the directed rounding of the corners rebuilt from the reported first-order fields.
Keys without a reference (no_reference) stay in the ledger and are counted separately.  A row is
"reference insufficient" (never a pass) if more than one tenth of its keys have no reference or an
identification key (largest order-2 total bin and its neighbours) has none (6.3 item 5).

Scopes: the routine part runs the CI row of the table (RB1); the full part (all nine RM-B rows)
runs only with XPSLOSS_M4B1_FULL=1 (prereg 5.3).  A missing field, a missing public name or a
missing reference is a failure, not a skip.
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

ROOT = Path(__file__).resolve().parents[1]
TABLE_FILE = ROOT / "benchmarks" / "numerics" / "m4b1" / "expected_states_m4b1.json"
TABLE_SHA256 = "d8d0f892a49c3ed9933b08d6a7c92e5b3aaf05ee5d6c7575050fc164991035e2"
REFERENCE_FILE = ROOT / "benchmarks" / "numerics" / "m4b1" / "references_rb_m4b1.json"
GENERATOR_FILE = ROOT / "benchmarks" / "numerics" / "m4b1" / "generate_references_rb_m4b1.py"
PREREG = ROOT / "docs" / "M4B1_PREREG.md"
FULL = os.environ.get("XPSLOSS_M4B1_FULL") == "1"

EPS = Fraction(1, 2 ** 52)
RB_ROWS = ("RB1", "RB2", "RB3", "RB4", "RB5", "RB6", "RB7", "RB8", "RB1a")
LEDGER = 6096
FR_MAX = Fraction(125, 10 ** 7)
DREL_CAP = Fraction(1, 10 ** 4)
H_FLOOR = Fraction(1, 10 ** 12)
NORM_SLACK = Fraction(1, 1000)
STORE_TOL = Fraction(1, 10 ** 30)        # stored decimals carry 40 digits
COMPONENTS = ("total", "no_hole", "hole")
COLUMNS = ("term_13a", "term_13b", "term_14", "term_15_T1", "term_15_T2", "term_15_T3",
           "term_16_e", "term_16_h", "no_hole", "hole", "total")
NUMERIC_KINDS = ("direct_reference", "interval_contains")

_CACHE = {}


class MissingField(AssertionError):
    pass


def need(obj, name, where):
    if isinstance(obj, Mapping):
        if name not in obj:
            raise MissingField(f"missing field {name!r} in {where}")
        return obj[name]
    if not hasattr(obj, name):
        raise MissingField(f"missing field {name!r} in {where}")
    return getattr(obj, name)


def F(text):
    return Fraction(text)


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


# ----------------------------------------------------------------------------- frozen files

def table():
    if "table" not in _CACHE:
        raw = TABLE_FILE.read_bytes()
        digest = hashlib.sha256(raw).hexdigest()
        if digest != TABLE_SHA256:
            raise AssertionError(f"expected_states_m4b1.json has SHA-256 {digest}, not the frozen one")
        doc = json.loads(raw.decode("utf-8"))
        doc["_rows"] = {r["id"]: r for r in doc["rows"]}
        _CACHE["table"] = doc
    return _CACHE["table"]


def frozen_reference_sha256():
    """SHA-256 of references_rb_m4b1.json recorded in the appendix (R-B entry)."""
    text = PREREG.read_text(encoding="utf-8")
    appendix = text[text.rindex("## 付録"):]
    found = re.findall(r"references_rb_m4b1\.json` の SHA-256: `([0-9a-f]{64})`", appendix)
    if not found:
        raise AssertionError("the appendix of docs/M4B1_PREREG.md does not record the SHA-256 of references_rb_m4b1.json")
    return found[0]


def references():
    if "refs" not in _CACHE:
        if not REFERENCE_FILE.exists():
            raise AssertionError("benchmarks/numerics/m4b1/references_rb_m4b1.json is missing")
        raw = REFERENCE_FILE.read_bytes()
        digest = hashlib.sha256(raw).hexdigest()
        frozen = frozen_reference_sha256()
        if digest != frozen:
            raise AssertionError(f"references_rb_m4b1.json has SHA-256 {digest}, the pre-registration records {frozen}")
        doc = json.loads(raw.decode("ascii"))
        for row in doc["rows"].values():
            row["_by_key"] = {}
            for rec in row["records"]:
                if (rec["key"], rec["kind"]) in row["_by_key"]:
                    raise AssertionError(f"duplicate reference key {rec['key']}")
                row["_by_key"][(rec["key"], rec["kind"])] = rec
        _CACHE["refs"] = doc
    return _CACHE["refs"]


# ----------------------------------------------------------------------------- inputs and the ledger

def row_inputs(row_id):
    tab = table()
    row = tab["_rows"][row_id]
    inp = row["inputs"]
    pub = tab["defaults"]["public"]
    first = {**tab["defaults"]["first_order_settings"]}
    order = {**tab["defaults"]["order_settings"]}
    for k, v in row["settings"].items():
        (first if k in first else order)[k] = v
    m = tab["materials"][inp["material"]]
    spec = {
        "row": row, "material_name": inp["material"], "gap": float.fromhex(m["gap_hex"]),
        "oscillators": [tuple(float.fromhex(h) for h in o) for o in m["oscillators_hex"]],
        "E0": float(pub["electron_energy_eV"]), "L": float(inp.get("weight_length_A", pub["weight_length_A"])),
        "Q": float(inp.get("hole_charge", pub["hole_charge"])),
        "mode": inp.get("epsilon_momentum", pub["epsilon_momentum"]),
        "window": tuple(float.fromhex(x["hex"]) for x in inp["loss_window_eV"]),
        "model": inp["model"], "max_order": inp["max_order"],
        "multiplier": inp.get("multiplier", pub["multiplier"]), "multiplier_A": inp.get("multiplier_A"),
        "overrides": dict(row["settings"]), "first": first, "order": order,
        "primary": None, "edges": None,
    }
    spec["m"] = spec["L"] if spec["multiplier"] == "weight_length" else float(spec["multiplier_A"])
    if "primary" in inp:
        p = tab["primaries"][inp["primary"]]
        field = "energies_eV" if p["kind"] == "delta" else "edges_eV"
        spec["primary"] = {"kind": p["kind"], "positions": [float.fromhex(h) for h in p[field]["hex"]],
                           "intensities": [float.fromhex(h) for h in p["intensities"]["hex"]]}
        spec["edges"] = [float.fromhex(h) for h in tab["output_edges"][inp["output_edges_eV"]]["hex"]]
    return spec


def planned_numeric(row_id):
    """The R-B key groups of the frozen ledger expanded to concrete (key, kind)."""
    key = ("planned", row_id)
    if key in _CACHE:
        return _CACHE[key]
    spec = row_inputs(row_id)
    out = []
    for item in spec["row"]["check_plan"]:
        if item["oracle"] != "R-B":
            continue
        if item["kind"] not in NUMERIC_KINDS:
            raise AssertionError(f"{row_id}: unexpected R-B kind {item['kind']}")
        name = item["keys"]
        if name == "order[2..3].<outside.low|bin[i]|outside.high>.<total|no_hole|hole>":
            nb = len(spec["edges"]) - 1
            regions = ["outside.low"] + [f"bin[{i}]" for i in range(nb)] + ["outside.high"]
            keys = [f"order[{n}].{r}.{c}" for n in (2, 3) for c in COMPONENTS for r in regions]
        elif name == "order[2..3].area.<11 columns>":
            keys = [f"order[{n}].area.{c}" for n in (2, 3) for c in COLUMNS]
        elif name == "order[1..3].cutoff_sensitivity.<11 columns>":
            keys = [f"order[{n}].cutoff_sensitivity.{c}" for n in (1, 2, 3) for c in COLUMNS]
        elif name == "order_area_ratio":
            keys = ["order_area_ratio"]
        else:
            raise AssertionError(f"{row_id}: unknown R-B key group {name!r}")
        if len(keys) != item["count"]:
            raise AssertionError(f"{row_id}: group {name!r} expands to {len(keys)} keys, the table has {item['count']}")
        out += [(k, item["kind"]) for k in keys]
    if len(set(out)) != len(out) or not out:
        raise AssertionError(f"{row_id}: empty or duplicate expansion")
    _CACHE[key] = out
    return out


def parse_key(key):
    """(order, family, where, column)."""
    m = re.fullmatch(r"order\[(\d+)\]\.(.*)", key)
    if not m:
        raise KeyError(key)
    n, rest = int(m.group(1)), m.group(2)
    if rest.startswith("area."):
        return n, "area", None, rest[5:]
    if rest.startswith("cutoff_sensitivity."):
        return n, "sensitivity", None, rest[19:]
    for side in ("low", "high"):
        if rest.startswith(f"outside.{side}."):
            return n, "outside", side, rest[len(side) + 9:]
    b = re.fullmatch(r"bin\[(\d+)\]\.(\w+)", rest)
    if not b:
        raise KeyError(key)
    return n, "bin", int(b.group(1)), b.group(2)


# ----------------------------------------------------------------------------- running the rows

def run_row(row_id):
    """Run the public entry of the row once per process."""
    key = ("run", row_id)
    if key in _CACHE:
        return _CACHE[key]
    import xpsloss
    spec = row_inputs(row_id)
    material = xpsloss.Material(spec["material_name"], spec["gap"],
                                tuple(xpsloss.Oscillator(*o) for o in spec["oscillators"]))
    kwargs = dict(loss_window_eV=spec["window"], weight_length_A=spec["L"], epsilon_momentum=spec["mode"],
                  model=spec["model"], max_order=spec["max_order"], hole_charge=spec["Q"], **spec["overrides"])
    if spec["multiplier"] != "weight_length":
        kwargs.update(multiplier=spec["multiplier"], multiplier_A=float(spec["multiplier_A"]))
    entry = spec["row"]["entry"]
    if entry == "api_spectrum":
        p = spec["primary"]
        cls = xpsloss.DeltaPrimary if p["kind"] == "delta" else xpsloss.BinnedPrimary
        primary = cls(p["positions"], p["intensities"])
        function = need(xpsloss, "multiple_loss_spectrum", "xpsloss")
        result = function(material, spec["E0"], primary, spec["edges"], **kwargs)
    elif entry == "api_areas":
        function = need(xpsloss, "multiple_loss_order_areas", "xpsloss")
        result = function(material, spec["E0"], **kwargs)
    else:
        raise AssertionError(f"{row_id}: entry {entry} is not an R-B entry")
    _CACHE[key] = (spec, result)
    return _CACHE[key]


def order_element(result, n, where):
    found = [o for o in need(result, "orders", where) if need(o, "order", where) == n]
    if len(found) != 1:
        raise MissingField(f"{where}: {len(found)} elements of orders have order == {n}")
    return found[0]


def observed(result, key):
    """(y, e_y, parts) of a direct_reference key; parts is {'inherited', 'local'} or None."""
    w = f"result[{key}]"
    n, family, where, column = parse_key(key)
    o = order_element(result, n, w)

    def parts_of(read):
        if n == 1:
            return None
        return {p: as_float(read(p), f"{w} {p}") for p in ("inherited", "local")}

    if family == "area":
        box = need(o, "area_error_parts", w)
        return (as_float(need(need(o, "area", w), column, w), w),
                as_float(need(need(o, "area_error", w), column, w), w + " error"),
                parts_of(lambda p: need(need(box, p, w), column, w)))
    if family == "sensitivity":
        return (as_float(need(need(o, "cutoff_sensitivity", w), column, w), w),
                as_float(need(need(o, "cutoff_sensitivity_error", w), column, w), w + " error"), None)
    if family == "outside":
        box = need(need(o, "outside", w), where, w)
        ep = need(box, column + "_error_parts", w)
        return (as_float(need(box, column, w), w), as_float(need(box, column + "_error", w), w + " error"),
                parts_of(lambda p: need(ep, p, w)))
    ep = need(o, column + "_error_parts", w)
    return (as_float(need(o, column, w)[where], w), as_float(need(o, column + "_error", w)[where], w + " error"),
            parts_of(lambda p: need(ep, p, w)[where]))


# ----------------------------------------------------------------------------- comparison rules

def base_d_rel(spec, n):
    return n * Fraction(spec["first"]["epsrel"]) + (n - 1) * Fraction(spec["order"]["order_epsrel"])


def phi_psi(spec, norms, n, column):
    r_abs, w_abs = F(norms["r_abs"]), F(norms["w_abs"][column])
    phi = r_abs ** (n - 1) + (n - 1) * w_abs * r_abs ** (n - 2)
    psi = sum(r_abs ** (n - k) for k in range(2, n + 1))
    return phi, psi, r_abs


def s_abs_of(spec):
    return sum(abs(Fraction(x)) for x in spec["primary"]["intensities"])


def recomputed_width(spec, norms, key):
    """(d_abs, base d_rel) of prereg 7.2 from the stored norms and the settings of the row."""
    n, family, _, column = parse_key(key)
    m, a = Fraction(spec["m"]), Fraction(spec["window"][0])
    phi, psi, r_abs = phi_psi(spec, norms, n, column)
    eta_o = Fraction(spec["order"]["order_epsabs_fraction"])
    if family == "area":
        d_abs = phi * m * Fraction(spec["first"]["area_epsabs"]) + psi * eta_o
    elif family == "sensitivity":
        d_abs = n * m * a * H_FLOOR * max(Fraction(1), r_abs) ** (n - 1)
    else:
        d_abs = (phi * Fraction(spec["first"]["spectrum_epsabs_fraction"]) + psi * eta_o) * s_abs_of(spec)
    return d_abs, base_d_rel(spec, n)


def has_value(rec):
    return rec["value"] is not None and rec["no_value_reason"] is None


def compare_direct(rec, y, e_y):
    r, u = F(rec["value"]), F(rec["u_r"])
    d_abs, d_rel = F(rec["width"]["d_abs"]), F(rec["width"]["d_rel"])
    w = d_abs if r == 0 else d_abs + d_rel * abs(r)
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


def budgets(spec, norms, key, y):
    """(budget (i) absolute, (i) relative, budget (iii)) of design 6.2 for a key of order n >= 2."""
    n, family, _, column = parse_key(key)
    phi, _, _ = phi_psi(spec, norms, n, column)
    eps, eps_o = Fraction(spec["first"]["epsrel"]), Fraction(spec["order"]["order_epsrel"])
    eta_o = Fraction(spec["order"]["order_epsabs_fraction"])
    ay = abs(Fraction(y))
    if family == "area":
        return (phi * Fraction(spec["m"]) * Fraction(spec["first"]["area_epsabs"]), n * eps * ay,
                max(eta_o, eps_o * ay))
    s_abs = s_abs_of(spec)
    return (phi * Fraction(spec["first"]["spectrum_epsabs_fraction"]) * s_abs, n * eps * ay,
            max(eta_o * s_abs, eps_o * ay))


def check_budgets(spec, norms, key, y, parts):
    b_abs, b_rel, b_local = budgets(spec, norms, key, y)
    inherited, local = Fraction(parts["inherited"]), Fraction(parts["local"])
    detail = {"inherited": float(inherited), "budget_i": float(max(b_abs, b_rel)), "local": float(local),
              "budget_iii": float(b_local)}
    if inherited < 0 or local < 0 or local > b_local:
        return "fail", detail
    if inherited <= max(b_abs, b_rel):
        return "pass", detail
    if inherited <= max(b_abs * (1 + NORM_SLACK), b_rel):
        return "unresolved", detail
    return "fail", detail


def compare_ratio(rec, result):
    w = "result order_area_ratio"
    ratio = need(result, "order_area_ratio", w)
    interval = need(ratio, "interval", w)
    lo, hi = as_float(interval[0], w + " lower"), as_float(interval[1], w + " upper")
    value = as_float(need(ratio, "value", w), w + " value")
    r, u = F(rec["value"]), F(rec["u_r"])
    if Fraction(lo) <= r - u and r + u <= Fraction(hi):
        contains = "pass"
    elif Fraction(hi) < r - u or Fraction(lo) > r + u:
        contains = "fail"
    else:
        contains = "unresolved"
    single = need(result, "single_loss", w)
    area = Fraction(as_float(need(need(single, "window_area", w), "term_13a", w), w + " area"))
    err = Fraction(as_float(need(need(single, "window_area_error", w), "term_13a", w), w + " area error"))
    m_hat = Fraction(as_float(need(single, "multiplier_A", w), w + " multiplier"))
    m_int = need(single, "multiplier_interval_A", w)
    ms = (m_hat, m_hat) if m_int is None else (Fraction(as_float(m_int[0], w)), Fraction(as_float(m_int[1], w)))
    corners = [m * x for m in ms for x in (area - err, area + err)]
    rebuilt = (lo == round_down(min(corners)) and hi == round_up(max(corners)) and value == float(m_hat * area))
    status = "fail" if contains == "fail" or not rebuilt else contains
    return status, {"contains": contains, "rebuilt_bitwise": rebuilt, "interval": (lo, hi), "value": value}


class Tally:
    def __init__(self):
        self.status = Counter()
        self.bad = []
        self.used = Counter()

    def add(self, row, key, status, detail):
        self.status[status] += 1
        self.used[(row, key)] += 1
        if status in ("fail", "unresolved"):
            self.bad.append((row, key, status, detail))

    def message(self):
        lines = [f"counts {dict(self.status)}"]
        lines += [f"  {r} {k}: {s} {d}" for r, k, s, d in self.bad[:40]]
        if len(self.bad) > 40:
            lines.append(f"  ... {len(self.bad) - 40} more")
        return "\n".join(lines)


def insufficiency(row_id):
    """Reasons that make the row 'reference insufficient' (prereg 6.3 item 5); empty if none."""
    refs = references()["rows"][row_id]
    planned = planned_numeric(row_id)
    recs = [refs["_by_key"].get(k) for k in planned]
    missing = [k for k, rec in zip(planned, recs) if rec is None or not has_value(rec)]
    reasons = []
    if 10 * len(missing) > len(planned):
        reasons.append(f"{len(missing)} of {len(planned)} keys have no reference")
    if table()["_rows"][row_id]["entry"] == "api_spectrum":
        ident = refs["identification"]
        nb = len(row_inputs(row_id)["edges"]) - 1
        values = [refs["_by_key"][(f"order[2].bin[{i}].total", "direct_reference")] for i in range(nb)]
        if any(not has_value(v) for v in values):
            best = None
        else:
            best = max(range(nb), key=lambda i: (F(values[i]["value"]), -i))
        want = None if best is None else [i for i in (best - 1, best, best + 1) if 0 <= i < nb]
        if want is None or ident["bins"] != want or any(not has_value(values[i]) for i in want):
            reasons.append("an identification key (largest order-2 total bin and its neighbours) has no reference")
    else:
        must = [f"order[{n}].area.total" for n in (2, 3)] + [f"order[{n}].cutoff_sensitivity.total" for n in (1, 2, 3)]
        if any(not has_value(refs["_by_key"][(k, "direct_reference")]) for k in must):
            reasons.append("a required area or sensitivity key has no reference")
    return reasons


def compare_result(row_id, spec, result, tally):
    refs = references()["rows"][row_id]
    statuses = {}
    for key, kind in planned_numeric(row_id):
        rec = need(refs["_by_key"], (key, kind), f"references {row_id}")
        if not has_value(rec):
            tally.add(row_id, key, "no_reference", rec["no_value_reason"])
            statuses[key] = "no_reference"
            continue
        if kind == "interval_contains":
            status, detail = compare_ratio(rec, result)
        else:
            y, e_y, parts = observed(result, key)
            status, detail = compare_direct(rec, y, e_y)
            if parts is not None:
                b_status, b_detail = check_budgets(spec, refs["norms"], key, y, parts)
                detail = {**detail, **b_detail}
                if b_status == "fail" or (b_status == "unresolved" and status == "pass"):
                    status = b_status
        tally.add(row_id, key, status, detail)
        statuses[key] = status
    return statuses


def run_and_compare(testcase, rows):
    tally = Tally()
    for rid in rows:
        testcase.assertEqual(insufficiency(rid), [], f"{rid}: reference insufficient")
        spec, result = run_row(rid)
        compare_result(rid, spec, result, tally)
    testcase.assertEqual(tally.status["fail"], 0, tally.message())
    testcase.assertEqual(tally.status["unresolved"], 0, tally.message())
    planned = sum(len(planned_numeric(r)) for r in rows)
    testcase.assertEqual(sum(tally.status.values()), planned, tally.message())
    testcase.assertTrue(all(v == 1 for v in tally.used.values()), "a reference key was used more than once")
    return tally


# ----------------------------------------------------------------------------- tests without the implementation

class ReferenceFile(unittest.TestCase):
    """The frozen reference file against the ledger and the rules of prereg 6.3 and 7.2."""

    def test_identity_and_counts(self):
        refs = references()
        self.assertEqual(refs["table"]["sha256"], TABLE_SHA256)
        self.assertEqual(refs["generator"]["sha256"], hashlib.sha256(GENERATOR_FILE.read_bytes()).hexdigest())
        self.assertEqual(sorted(refs["rows"]), sorted(RB_ROWS))
        total = 0
        for rid in RB_ROWS:
            planned = planned_numeric(rid)
            stored = refs["rows"][rid]["_by_key"]
            self.assertEqual(set(stored), set(planned), f"{rid}: unused or missing reference keys")
            total += len(planned)
        self.assertEqual(total, LEDGER)
        self.assertEqual(total, table()["counts"]["planned_numeric_keys_by_group"]["RM-B"])
        self.assertEqual(refs["counts"]["total"], LEDGER)

    def test_no_reference_rule_and_widths(self):
        """6.3 and 7.2 re-applied to the stored decimals (30 digits or more)."""
        refs = references()
        listed = {(x["row"], x["key"]) for x in refs["no_reference"]}
        seen = set()
        for rid in RB_ROWS:
            spec = row_inputs(rid)
            row = refs["rows"][rid]
            for (key, kind), rec in row["_by_key"].items():
                self.assertEqual(rec["oracle"], "R-B")
                if not has_value(rec):
                    seen.add((rid, key))
                    self.assertIn(rec["no_value_reason"], ("f_r_exceeds_1.25e-5", "zero_value_u_r_exceeds_quarter_floor"))
                    r, u = F(rec["rejected"]["value"]), F(rec["rejected"]["u_r"])
                    if r != 0:
                        self.assertGreater(u / abs(r), FR_MAX, f"{rid} {key}")
                    continue
                r, u = F(rec["value"]), F(rec["u_r"])
                self.assertGreaterEqual(u, 0)
                if r != 0:
                    digits = len(rec["value"].lstrip("-").split("e")[0].replace(".", ""))
                    self.assertGreaterEqual(digits, 30, f"{rid} {key}")
                    self.assertGreaterEqual(u, sum(F(x) for x in rec["u_r_parts"].values()) * (1 - STORE_TOL))
                    self.assertLessEqual(u / abs(r), FR_MAX, f"{rid} {key}")
                if kind == "interval_contains":
                    lo, hi = (F(x) for x in rec["interval"])
                    self.assertLessEqual(lo, r - u)
                    self.assertGreaterEqual(hi, r + u)
                    continue
                d_abs, base = recomputed_width(spec, row["norms"], key)
                self.assertLessEqual(abs(F(rec["width"]["d_abs"]) - d_abs), STORE_TOL * d_abs, f"{rid} {key}")
                self.assertLessEqual(F(rec["width"]["d_abs"]), d_abs * (1 + STORE_TOL))
                if r == 0:
                    self.assertLessEqual(u, d_abs / 4, f"{rid} {key}")
                    want = base
                else:
                    want = min(max(8 * u / abs(r), base), DREL_CAP)
                self.assertLessEqual(abs(F(rec["width"]["d_rel"]) - want), STORE_TOL * want, f"{rid} {key}")
                self.assertLessEqual(F(rec["width"]["d_rel"]), DREL_CAP)
        self.assertEqual(seen, listed)
        self.assertEqual(refs["counts"]["no_reference"], len(listed))
        self.assertEqual(refs["reference_path_failures"], [])

    def test_rows_are_sufficient(self):
        refs = references()
        for rid in RB_ROWS:
            self.assertEqual(insufficiency(rid), [], f"{rid}: reference insufficient")
        self.assertEqual(refs["insufficient_rows"], [])

    def test_records_made_before_the_comparison(self):
        """6.2: the two momentum prescriptions differ by more than the widths; X01 leaves the RB1 width."""
        pre = references()["pre_comparison_records"]
        for item in pre["momentum_prescriptions"]:
            self.assertEqual(len(item["bins"]), 3, item["rows"])
            self.assertTrue(all(b["exceeds"] for b in item["bins"]), item["rows"])
        x01 = pre["X01_wrong_answer_vs_RB1_width"]
        self.assertEqual(x01["row"], "RB1")
        order2 = x01["orders"]["2"]
        self.assertEqual(len(order2["identification_bins"]), 3)
        self.assertTrue(all(b["exceeds_width"] for b in order2["identification_bins"]))
        self.assertGreater(order2["bins_exceeding_count"], 0)


class _Synthetic(dict):
    """A result built from the reference values themselves (mappings only)."""

    def __init__(self, row_id, source=None, shift=None):
        super().__init__()
        refs = references()["rows"][source or row_id]["_by_key"]
        shift = shift or {}

        def val(key):
            rec = refs[(key, "direct_reference")]
            text = rec["value"] if has_value(rec) else rec["rejected"]["value"]
            u = F(rec["u_r"]) if has_value(rec) else Fraction(0)
            y = float(F(text)) + shift.get(key, 0.0)
            return y, round_up(2 * (u + abs(Fraction(y) - F(text))))

        spec = row_inputs(row_id)
        orders = []
        for n in (1, 2, 3):
            o = {"order": n, "area": {}, "area_error": {}, "cutoff_sensitivity": {}, "cutoff_sensitivity_error": {},
                 "area_error_parts": {"inherited": {}, "local": {}}}
            for c in COLUMNS:
                if (f"order[{n}].area.{c}", "direct_reference") in refs:
                    o["area"][c], o["area_error"][c] = val(f"order[{n}].area.{c}")
                    o["area_error_parts"]["inherited"][c] = o["area_error_parts"]["local"][c] = 0.0
                if (f"order[{n}].cutoff_sensitivity.{c}", "direct_reference") in refs:
                    o["cutoff_sensitivity"][c], o["cutoff_sensitivity_error"][c] = val(f"order[{n}].cutoff_sensitivity.{c}")
            if spec["edges"] and n >= 2:
                nb = len(spec["edges"]) - 1
                o["outside"] = {"low": {}, "high": {}}
                for c in COMPONENTS:
                    pairs = [val(f"order[{n}].bin[{i}].{c}") for i in range(nb)]
                    o[c], o[c + "_error"] = [p[0] for p in pairs], [p[1] for p in pairs]
                    o[c + "_error_parts"] = {"inherited": [0.0] * nb, "local": [0.0] * nb}
                    for side in ("low", "high"):
                        box = o["outside"][side]
                        box[c], box[c + "_error"] = val(f"order[{n}].outside.{side}.{c}")
                        box[c + "_error_parts"] = {"inherited": 0.0, "local": 0.0}
            orders.append(o)
        self["orders"] = orders
        rec = refs[("order_area_ratio", "interval_contains")]
        m = Fraction(spec["m"])
        area = float(F(rec["value"]) / m)
        err = round_up(2 * (F(rec["u_r"]) / m + abs(Fraction(area) - F(rec["value"]) / m)))
        corners = [m * (Fraction(area) - Fraction(err)), m * (Fraction(area) + Fraction(err))]
        self["order_area_ratio"] = {"value": float(m * Fraction(area)),
                                    "interval": (round_down(min(corners)), round_up(max(corners)))}
        self["single_loss"] = {"window_area": {"term_13a": area}, "window_area_error": {"term_13a": err},
                               "multiplier_A": spec["m"], "multiplier_interval_A": None}


class ComparatorSelfCheck(unittest.TestCase):
    """The comparison rules of this file on synthetic results (no implementation is called)."""

    def run_synthetic(self, row_id, **kw):
        tally = Tally()
        compare_result(row_id, row_inputs(row_id), _Synthetic(row_id, **kw), tally)
        return tally

    def test_reference_values_pass(self):
        for rid in ("RB1", "RB3", "RB7", "RB1a"):
            tally = self.run_synthetic(rid)
            self.assertEqual(tally.status["fail"] + tally.status["unresolved"], 0, tally.message())
            self.assertEqual(sum(tally.status.values()), len(planned_numeric(rid)))

    def test_shifted_key_fails(self):
        refs = references()["rows"]["RB1"]
        i = refs["identification"]["largest_order2_total_bin"]
        key = f"order[2].bin[{i}].total"
        rec = refs["_by_key"][(key, "direct_reference")]
        width = F(rec["width"]["d_abs"]) + F(rec["width"]["d_rel"]) * abs(F(rec["value"]))
        tally = self.run_synthetic("RB1", shift={key: float(2 * width)})
        self.assertEqual([k for _, k, s, _ in tally.bad if s == "fail"], [key], tally.message())

    def test_swapped_momentum_prescription_fails(self):
        """The k_parallel values (RB2) offered for RB1 leave the widths in the identification bins."""
        tally = self.run_synthetic("RB1", source="RB2")
        failed = {k for _, k, s, _ in tally.bad if s == "fail"}
        self.assertTrue(set(references()["rows"]["RB1"]["identification"]["keys"]) <= failed, tally.message())

    def test_budget_rules(self):
        spec = row_inputs("RB1")
        norms = references()["rows"]["RB1"]["norms"]
        key = "order[2].bin[90].total"
        b_abs, b_rel, b_local = budgets(spec, norms, key, 0.0)
        self.assertEqual(check_budgets(spec, norms, key, 0.0, {"inherited": float(b_abs) * 0.5, "local": 0.0})[0], "pass")
        self.assertEqual(check_budgets(spec, norms, key, 0.0, {"inherited": float(b_abs) * 1.0005, "local": 0.0})[0],
                         "unresolved")
        self.assertEqual(check_budgets(spec, norms, key, 0.0, {"inherited": float(b_abs) * 1.1, "local": 0.0})[0], "fail")
        self.assertEqual(check_budgets(spec, norms, key, 0.0, {"inherited": 0.0, "local": float(b_local) * 1.1})[0], "fail")


# ----------------------------------------------------------------------------- tests of the implementation

class ReferenceRoutine(unittest.TestCase):
    """CI rows of the table among the RM-B rows (RB1)."""

    def test_ci_rows(self):
        rows = [r for r in RB_ROWS if table()["_rows"][r]["ci"]]
        self.assertEqual(rows, ["RB1"])
        run_and_compare(self, rows)


@unittest.skipUnless(FULL, "set XPSLOSS_M4B1_FULL=1 for the full RM-B comparison (prereg 5.3)")
class ReferenceFull(unittest.TestCase):
    """All RM-B rows (prereg 5.3: RB2-RB8 and RB1a are full-only)."""

    def test_all_rows(self):
        tally = run_and_compare(self, list(RB_ROWS))
        self.assertEqual(sum(tally.status.values()), LEDGER)


if __name__ == "__main__":
    unittest.main()
