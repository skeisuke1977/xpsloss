#!/usr/bin/env python3
"""Mapped reference table for M4a-2 (docs/M4A2_PREREG.md, sections 4.1-4.5).

Reads only frozen JSON files:

* benchmarks/numerics/m4a/references_rb_m4a.json   (the independent R-B / CE references of M4a)
* benchmarks/numerics/m4a/expected_states_m4a.json (hash check only)
* benchmarks/numerics/m4a2/expected_states_m4a2.json (planned rows and keys)

and writes benchmarks/numerics/m4a2/references_m4a2.json.  Only the standard
library is used and the package under test is never imported.  All arithmetic
is exact (fractions.Fraction); decimals are written with 40 significant digits.

Usage:
    generate_references_m4a2.py            write the table
    generate_references_m4a2.py --check    regenerate in memory and compare bytes
    generate_references_m4a2.py --evaluate-fixed-m ROW_ID M_HAT_HEX
                                           evaluate the frozen V4 fixed_m formulas
                                           of one row for a reported binary64 m_hat
"""
import argparse
import hashlib
import json
import re
import sys
from fractions import Fraction as F
from pathlib import Path

HERE = Path(__file__).resolve().parent
RB_PATH = HERE.parent / "m4a" / "references_rb_m4a.json"
M4A_STATES_PATH = HERE.parent / "m4a" / "expected_states_m4a.json"
STATES_PATH = HERE / "expected_states_m4a2.json"
OUT_PATH = HERE / "references_m4a2.json"

RB_SHA = "942353c748b4959220ba4da15f078a8ebdbdd9833d3abdc2e80cb8acfd6e8d94"
M4A_STATES_SHA = "250b45f863ccbd1c4c94474f23769976380c5dad5c2671328712ef4acdc7f408"
STATES_SHA = "4a8264a58b9902efee3d843c10be812569ff6f46d720c9f23ebf37a4f263073f"

SCHEMA = "xpsloss-m4a2-references-v1"
DIGITS = 40
D_REL_LIMIT = F(1, 10 ** 4)       # prereg 4.4(1): above this the key is no_reference
H_FLOOR = F(1, 10 ** 12)          # prereg 4.2: the 1e-12 of the sensitivity floors (exact decimal)
CE_AREA_FLOOR = F(1, 10 ** 18)    # prereg 4.2: B_A^0 of the CE rows
C_VALUES = {"z": F(0), "p": F(1, 32)}
RB_IDS = ["RB1", "RB2", "RB3", "RB4", "RB5", "RB6", "RB7", "RB8", "RB1a", "RB2a", "RB3a", "RB4a"]
EXTERNAL_Q = 2                    # V3: m = 2 * L_ref

# Planned keys that are checked without a mapped reference (bit / rational
# checks of the run itself); they are listed in the metadata, not as records.
UNMAPPED = {
    "s2_map": {
        "denominator": "4.4(5), 4.5(2): beta-stage denominator of the run; no mapped reference",
        "multiplier_A": "4.4(5): bit check RN(1/A_total^beta) against the run's own area",
        "multiplier_interval_A": "4.5(2): bit check against the run's own area and error (c = 0, total)",
        "returned_errors": "4.4(4): bit checks of the returned errors",
        "rho": "4.4(4): bit check of diagnostics.multiplier.rho",
    },
    "seed_reference": {
        "seed.error13": "4.5(1): e_y of the area comparison (an estimate, not a referenced value)",
        "weight_length_A": "4.5(2): bit check RN(1/D) against the run's own area and c",
    },
}


# ----------------------------------------------------------------------------
# exact numbers and decimals
# ----------------------------------------------------------------------------

def fr_hex(text):
    """binary64 hex string -> exact Fraction."""
    return F(float.fromhex(text))


def dec(x, mode="n"):
    """Decimal string with DIGITS significant digits.

    mode "n": nearest, ties to even; "d": toward -inf; "u": toward +inf.
    An exact zero is "0".
    """
    x = F(x)
    if x == 0:
        return "0"
    neg = x < 0
    a = -x if neg else x
    if neg and mode != "n":
        mode = "u" if mode == "d" else "d"
    e = len(str(a.numerator)) - len(str(a.denominator))
    while a >= F(10) ** (e + 1):
        e += 1
    while a < F(10) ** e:
        e -= 1
    m = a / F(10) ** (e - (DIGITS - 1))
    n, rem = divmod(m.numerator, m.denominator)
    if rem:
        if mode == "u":
            n += 1
        elif mode == "n":
            twice = 2 * rem
            if twice > m.denominator or (twice == m.denominator and n % 2):
                n += 1
    if n == 10 ** DIGITS:
        n //= 10
        e += 1
    s = str(n)
    return ("-" if neg else "") + s[0] + "." + s[1:] + ("" if e == 0 else "e%+d" % e)


def plain(x):
    """Exact decimal string of a Fraction with a power-of-two/five denominator."""
    x = F(x)
    if x.denominator == 1:
        return str(x.numerator)
    d, k = x.denominator, 0
    while d % 10 == 0:
        d //= 10
        k += 1
    while d % 2 == 0:
        d //= 2
        k += 1
    while d % 5 == 0:
        d //= 5
        k += 1
    if d != 1:
        raise ValueError("not a finite decimal: %s" % x)
    # k is an upper bound of the number of decimal places
    n = x * 10 ** k
    assert n.denominator == 1
    digits = str(abs(n.numerator)).rjust(k + 1, "0")
    text = digits[:-k] + "." + digits[-k:]
    text = text.rstrip("0").rstrip(".")
    return ("-" if x < 0 else "") + text


def frac(x):
    """Exact rational string "num/den" (or "num")."""
    x = F(x)
    return str(x.numerator) if x.denominator == 1 else "%d/%d" % (x.numerator, x.denominator)


def sha_file(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def sha_generator():
    """SHA-256 of this file with CRLF normalised to LF."""
    return hashlib.sha256(Path(__file__).read_bytes().replace(b"\r\n", b"\n")).hexdigest()


# ----------------------------------------------------------------------------
# judgement of one mapped key (prereg 4.4(1))
# ----------------------------------------------------------------------------

def finalize(center, lo, hi, d_abs, epsrel, want_interval=False):
    """Stored fields of a mapped key from exact center, exact interval and floor."""
    center, lo, hi, d_abs = F(center), F(lo), F(hi), F(d_abs)
    assert lo <= center <= hi
    value = dec(center, "n")
    vs = F(value)
    u_exact = max(center - lo, hi - center)
    u_r = dec(max(vs - lo, hi - vs), "u")
    us = F(u_r)
    out = {
        "value": value,
        "u_r": u_r,
        "d_abs": dec(d_abs, "d"),
        "d_abs_exact": frac(d_abs),
        "d_rel": None,
        "interval": [dec(lo, "d"), dec(hi, "u")] if want_interval else None,
        "status": "value",
        "reason": None,
    }
    if center == 0:
        if d_abs > 0 and u_exact <= d_abs / 4 and us <= d_abs / 4:
            out["reason"] = "r = 0, u_r <= d_abs/4 (4.4(1)); width d_abs"
        else:
            out["status"] = "no_reference"
            out["reason"] = "r = 0 and u_r > d_abs/4 (4.4(1))"
        return out
    d_rel = max(8 * u_exact / abs(center), epsrel)
    d_rel_stored = max(8 * us / abs(vs), epsrel)
    out["d_rel"] = dec(d_rel, "d")
    if d_rel > D_REL_LIMIT or d_rel_stored > D_REL_LIMIT:
        out["status"] = "no_reference"
        out["reason"] = "d_rel = max(8 u_r/|r|, epsrel) exceeds 1e-4 after the mapping (4.4(1))"
    else:
        out["reason"] = "r != 0, d_rel = max(8 u_r/|r|, epsrel) <= 1e-4 before and after decimal rounding (4.4(1))"
    return out


# ----------------------------------------------------------------------------
# source rows
# ----------------------------------------------------------------------------

class Source:
    """One row of the frozen R-B / CE table."""

    def __init__(self, row_id, row):
        self.id = row_id
        self.row = row
        self.failure = row.get("reference_path_failure")
        cond = row["conditions"]
        self.cond = cond
        self.family = row["family"]
        self.L = fr_hex(cond["weight_length_A"])
        self.a = fr_hex(cond["loss_window_eV"][0])
        cs = cond["comparison_settings"]
        self.epsrel = F(cs["epsrel"])
        self.area_floor = CE_AREA_FLOOR if self.family == "CE" else F(cs["area_epsabs"])
        self.area_floor_name = "1e-18 (CE)" if self.family == "CE" else "Fr(area_epsabs)"
        self.eta = F(cs["spectrum_epsabs_fraction"])
        self.sabs = None
        self.widths = []
        if "primary" in cond:
            self.sabs = sum(abs(fr_hex(x)) for x in cond["primary"]["intensities"])
            edges = [fr_hex(x) for x in cond["output_edges_eV"]]
            self.widths = [edges[i + 1] - edges[i] for i in range(len(edges) - 1)]
            assert all(w > 0 for w in self.widths)

    def get(self, key):
        """(r, u) exact, or a reason string when the key cannot be used."""
        k = self.row["keys"].get(key)
        if k is None:
            return "source key %s is missing in %s" % (key, self.id)
        if k.get("status") != "value":
            return "source key %s of %s has status %s" % (key, self.id, k.get("status"))
        if "r" not in k:
            return "source key %s of %s has no center" % (key, self.id)
        r, u = F(k["r"]), F(k["u_r"])
        if u < 0:
            return "source key %s of %s has a negative u_r" % (key, self.id)
        return r, u


REGION = r"(bin\[(\d+)\]|outside\.low|outside\.high)"
ZERO = (F(0), F(0))


def spec_for(src, key):
    """Linear form of a planned key in the multiplier m.

    value(m) = (const + (m/L) * prop) / width
    d_abs(m) = (floor_const + m * floor_per_m) / width
    Returns a dict, or a string (reason for no_reference).
    """
    width, width_note, base = F(1), None, key
    m = re.fullmatch(r"bin\[(\d+)\]\.(\w+)_per_eV", key)
    if m:
        i = int(m.group(1))
        if i >= len(src.widths):
            return "bin %d is outside the output edges of %s" % (i, src.id)
        width = src.widths[i]
        width_note = "conditions.output_edges_eV[%d..%d]" % (i, i + 1)
        base = "bin[%d].%s" % (i, m.group(2))
    sp = {"const": None, "prop": None, "floor_const": F(0), "floor_per_m": F(0), "width": width,
          "num_keys": [], "add_keys": [], "cond": ["conditions.weight_length_A"], "cls": None}
    if width_note:
        sp["cond"].append(width_note)

    def need(k):
        got = src.get(k)
        return got

    if base.startswith("window_area."):
        got = need(base)
        if isinstance(got, str):
            return got
        sp.update(const=got, floor_const=src.area_floor, add_keys=[base], cls="area",
                  transform="r", unit="A^-1", floor="B_A0", cond=["conditions.comparison_settings.area_epsabs"])
        return sp
    if base == "single_loss_weight":
        got = need(base)
        if isinstance(got, str):
            return got
        sp.update(prop=got, floor_per_m=src.area_floor, num_keys=[base], cls="weight",
                  transform="q*r", unit="1", floor="m*B_A0")
        sp["cond"].append("conditions.comparison_settings.area_epsabs")
        return sp
    if base.startswith("cutoff_sensitivity."):
        got = need(base)
        if isinstance(got, str):
            return got
        sp.update(prop=got, floor_per_m=src.a * H_FLOOR, num_keys=[base], cls="H",
                  transform="q*r", unit="1", floor="m*a*1e-12")
        sp["cond"].append("conditions.loss_window_eV[0]")
        return sp
    m = re.fullmatch(REGION + r"\.(.+)", base)
    if not m or src.sabs is None:
        return "key %s has no mapping for source row %s" % (key, src.id)
    region, what = m.group(1), m.group(3)
    per = "" if width == 1 and width_note is None else "/dE"
    unit = "I" if not per else "I/eV"
    eta_floor = src.eta * src.sabs
    if what == "zero_loss":
        got = need(base)
        if isinstance(got, str):
            return got
        sp.update(const=got, floor_const=eta_floor, add_keys=[base], cls="zero",
                  transform="r0" + per, unit=unit, floor="eta*S_abs" + per)
        sp["cond"] = [c for c in sp["cond"] if c != "conditions.weight_length_A"]
    elif what in ("single_loss_total", "single_loss_no_hole", "single_loss_hole"):
        got = need(base)
        if isinstance(got, str):
            return got
        sp.update(prop=got, floor_const=eta_floor, num_keys=[base], cls="satellite",
                  transform="q*r1" + per, unit=unit, floor="eta*S_abs" + per)
    elif what == "spectrum":
        k0, k1 = region + ".zero_loss", region + ".single_loss_total"
        g0, g1 = need(k0), need(k1)
        if isinstance(g0, str):
            return g0
        if isinstance(g1, str):
            return g1
        sp.update(const=g0, prop=g1, floor_const=eta_floor, add_keys=[k0], num_keys=[k1], cls="spectrum",
                  transform="(r0+q*r1)" + per if per else "r0+q*r1", unit=unit, floor="eta*S_abs" + per)
    elif re.fullmatch(r"cutoff_sensitivity\.(plus|minus)\.(total|no_hole|hole)", what) and not per:
        got = need(base)
        if isinstance(got, str):
            return got
        sp.update(prop=got, floor_per_m=src.a * H_FLOOR * src.sabs, num_keys=[base], cls="region_sensitivity",
                  transform="q*r", unit="I", floor="m*a*1e-12*S_abs")
        sp["cond"].append("conditions.loss_window_eV[0]")
    else:
        return "key %s has no mapping for source row %s" % (key, src.id)
    sp["cond"] += ["conditions.comparison_settings.spectrum_epsabs_fraction", "conditions.primary.intensities"] \
        if "eta" in sp["floor"] else ["conditions.primary.intensities"]
    return sp


def base_record(row_id, kind, key, src):
    return {
        "row_id": row_id, "check_kind": kind, "key": key,
        "source_sha256": RB_SHA, "source_row": src.id,
        "source_keys": None, "transform": None, "unit": None,
        "value": None, "u_r": None, "d_abs": None, "d_abs_exact": None, "d_abs_formula": None,
        "d_rel": None, "interval": None, "interval_exact": None,
        "status": None, "reason": None, "parameter": None, "terms": None, "counted_as": None,
    }


def fail_record(rec, src, reason):
    if src.failure is not None:
        rec["status"] = "reference_path_failure"
        rec["reason"] = "reference path of source row %s failed: %s" % (src.id, json.dumps(src.failure, sort_keys=True))
    else:
        rec["status"] = "no_reference"
        rec["reason"] = reason
    return rec


def key_class(key):
    if key.startswith("inverse_window_area_interval_A."):
        return "inverse_interval"
    if key.endswith("_per_eV"):
        return "derived_density"
    return "direct_reference"


def inverse_record(rec, src, area_key, c, endpoint):
    """Reference interval [1/(r+u+c), 1/(r-u+c)] of an area key (prereg 4.2, 4.5(3))."""
    rec["source_keys"] = {"numerator": [], "added": [], "denominator": [area_key], "conditions": []}
    rec["transform"] = "1/([r-u_r, r+u_r]+c)" if c else "1/[r-u_r, r+u_r]"
    rec["unit"] = "A"
    rec["counted_as"] = "inverse_interval"
    got = src.get(area_key)
    if src.failure is not None or isinstance(got, str):
        return fail_record(rec, src, got if isinstance(got, str) else None)
    r, u = got
    lo_den, hi_den = r - u + c, r + u + c
    if lo_den <= 0:
        return fail_record(rec, src, "the reference interval of the denominator contains 0 or is negative (4.4(1))")
    lo, hi = 1 / hi_den, 1 / lo_den
    rec["interval"] = [dec(lo, "d"), dec(hi, "u")]
    rec["interval_exact"] = [frac(lo), frac(hi)]
    if endpoint == "lower":
        rec["value"] = rec["interval"][0]
    elif endpoint == "upper":
        rec["value"] = rec["interval"][1]
    rec["status"] = "value"
    rec["reason"] = ("4.5(3): reference interval I compared exactly with the reported interval P; "
                     "I subset of P passes, disjoint fails, otherwise unresolved; no width is added to the ends")
    return rec


def mapped_record(row_id, kind, key, src):
    """One V3 / V4 record. kind: external, fixed_m, ratio."""
    rec = base_record(row_id, kind, key, src)
    cls = key_class(key)
    rec["counted_as"] = cls
    L = src.L
    if cls == "inverse_interval":
        rec = inverse_record(rec, src, "window_area.total", F(0), key.rsplit(".", 1)[1])
        rec["parameter"] = {"note": "independent of the multiplier"}
        return rec
    if kind == "ratio" and key == "single_loss_weight":
        # prereg 4.3 / 4.4(5): the same area cancels; counted as the rounding identity with target 1
        rec["source_keys"] = {"numerator": ["window_area.total"], "added": [], "denominator": ["window_area.total"],
                              "conditions": ["conditions.comparison_settings.area_epsabs"]}
        rec["transform"] = "identity: m_ref*A = 1"
        rec["unit"] = "1"
        rec["counted_as"] = "identity_4.4(5)"
        got = src.get("window_area.total")
        if src.failure is not None or isinstance(got, str):
            return fail_record(rec, src, got if isinstance(got, str) else None)
        A, uA = got
        if A - uA <= 0:
            return fail_record(rec, src, "the reference interval of the denominator contains 0 or is negative (4.4(1))")
        rec.update(finalize(1, 1, 1, src.area_floor / A, src.epsrel, want_interval=True))
        rec["d_abs_formula"] = "m_ref*B_A0"
        rec["parameter"] = ratio_parameter(src, A, uA)
        rec["reason"] = ("4.3, 4.4(5): the same area cancels exactly; target 1 is counted as the rounding identity "
                         "w = RN(m_hat*A_total), error = up(|w-1|), not as an independent numerical comparison")
        return rec
    sp = spec_for(src, key)
    if src.failure is not None or isinstance(sp, str):
        return fail_record(rec, src, sp if isinstance(sp, str) else None)
    width = sp["width"]
    const = sp["const"] or ZERO
    prop = sp["prop"] or ZERO
    m_dependent = sp["prop"] is not None or sp["floor_per_m"] != 0
    rec["unit"] = sp["unit"]
    rec["d_abs_formula"] = sp["floor"]
    rec["source_keys"] = {"numerator": sp["num_keys"], "added": sp["add_keys"], "denominator": [],
                          "conditions": sp["cond"]}
    rec["transform"] = sp["transform"]
    want_interval = False
    if kind == "external" or not m_dependent:
        mval = EXTERNAL_Q * L
        q = F(EXTERNAL_Q)
        center = (const[0] + q * prop[0]) / width
        half = (const[1] + q * prop[1]) / width
        lo, hi = center - half, center + half
        d_abs = (sp["floor_const"] + mval * sp["floor_per_m"]) / width
        if m_dependent:
            rec["parameter"] = {"m_A": plain(mval), "m_hex": float(mval).hex(), "reference_length_A": plain(L),
                                "q": plain(q)}
        else:
            rec["parameter"] = {"note": "independent of the multiplier"}
    elif kind == "fixed_m":
        rec["status"] = "parametric"
        rec["reason"] = ("m_hat is not known at freezing (4.1); evaluate once with the reported binary64 m_hat "
                         "by metadata.fixed_m_evaluation")
        rec["parameter"] = {"symbol": "m_hat", "reference_length_A": plain(L), "q": "Fr(m_hat)/reference_length_A"}
        rec["terms"] = {
            "constant": {"value": dec_exact(const[0]), "u_r": dec_exact(const[1])},
            "proportional": {"value": dec_exact(prop[0]), "u_r": dec_exact(prop[1])},
            "divide_by": plain(width),
            "d_abs": {"constant_exact": frac(sp["floor_const"]), "per_m_hat_exact": frac(sp["floor_per_m"])},
            "epsrel_exact": frac(src.epsrel),
        }
        return rec
    else:  # ratio, m-dependent
        rec["source_keys"]["denominator"] = ["window_area.total"]
        rec["transform"] = sp["transform"].replace("q*", "").replace("r1", "r1/(L*A)").replace(
            "r)", "r/(L*A))") if sp["cls"] != "weight" else sp["transform"]
        if sp["transform"] == "q*r":
            rec["transform"] = "r/(L*A)"
        got = src.get("window_area.total")
        if isinstance(got, str):
            return fail_record(rec, src, got)
        A, uA = got
        if A - uA <= 0:
            return fail_record(rec, src, "the reference interval of the denominator contains 0 or is negative (4.4(1))")
        a_lo, a_hi = A - uA, A + uA
        n, n_lo, n_hi = prop[0] / L, (prop[0] - prop[1]) / L, (prop[0] + prop[1]) / L
        corners = [n_lo / a_lo, n_lo / a_hi, n_hi / a_lo, n_hi / a_hi]
        center = (const[0] + n / A) / width
        lo = (const[0] - const[1] + min(corners)) / width
        hi = (const[0] + const[1] + max(corners)) / width
        d_abs = (sp["floor_const"] + sp["floor_per_m"] / A) / width
        rec["d_abs_formula"] = sp["floor"].replace("m*", "m_ref*")
        rec["parameter"] = ratio_parameter(src, A, uA)
        want_interval = True
    rec.update(finalize(center, lo, hi, d_abs, src.epsrel, want_interval=want_interval))
    return rec


def dec_exact(x):
    """Exact decimal of a source number (a 40-digit decimal read exactly)."""
    return plain(x)


def ratio_parameter(src, A, uA):
    return {
        "reference_length_A": plain(src.L),
        "A_total_per_A": dec_exact(A),
        "A_total_interval_per_A": [dec_exact(A - uA), dec_exact(A + uA)],
        "m_ref_A": dec(1 / A, "n"),
        "m_ref_exact": frac(1 / A),
        "q": dec(1 / (src.L * A), "n"),
        "q_exact": frac(1 / (src.L * A)),
    }


def seed_records(row, src):
    """V5 seed row: seed.area13 (seed_area) and weight_length_interval_A (inverse_interval)."""
    inp = row["inputs"]
    c = F(inp["c"])
    param = {"c_per_A": plain(c), "c_hex": float(c).hex(), "seed_A": inp.get("seed_A", inp.get("seeds_A"))}
    out = []
    for key in row["keys"]:
        if key == "seed.area13":
            rec = base_record(row["id"], "seed_area", key, src)
            rec["source_keys"] = {"numerator": [], "added": ["window_area.term_13a"], "denominator": [],
                                  "conditions": ["conditions.comparison_settings.area_epsabs"]}
            rec["transform"] = "r"
            rec["unit"] = "A^-1"
            rec["counted_as"] = "seed_area"
            rec["d_abs_formula"] = "B_A0"
            got = src.get("window_area.term_13a")
            if src.failure is not None or isinstance(got, str):
                fail_record(rec, src, got if isinstance(got, str) else None)
            else:
                r, u = got
                rec.update(finalize(r, r - u, r + u, src.area_floor, src.epsrel))
            rec["parameter"] = dict(param, note="c and the seed do not enter the area comparison (4.5(1))")
            out.append(rec)
        elif key in ("weight_length_interval_A", "common_interval"):
            rec = base_record(row["id"], "inverse_interval", key, src)
            rec = inverse_record(rec, src, "window_area.term_13a", c, None)
            rec["parameter"] = dict(param)
            if key == "common_interval" and rec["status"] == "value":
                rec["reason"] = ("4.5(3): the three reported intervals must share a point (max L_- <= min L_+); "
                                 "each seed is also judged on its own row against this same reference interval I, "
                                 "and I is inside the common interval when all three pass")
            out.append(rec)
        elif key not in UNMAPPED["seed_reference"]:
            raise SystemExit("unexpected V5 key %s in %s" % (key, row["id"]))
    return out


# ----------------------------------------------------------------------------
# evaluation of the frozen fixed_m formulas
# ----------------------------------------------------------------------------

FIXED_M_EVALUATION = [
    "m = Fr(m_hat): the reported binary64 multiplier read as an exact rational; m must be finite and > 0",
    "q = m / reference_length_A (exact)",
    "all strings in terms are exact rationals (decimals or num/den)",
    "r = (constant.value + q*proportional.value) / divide_by",
    "u = (constant.u_r + q*proportional.u_r) / divide_by",
    "d = (d_abs.constant_exact + m*d_abs.per_m_hat_exact) / divide_by",
    "value = r rounded to 40 significant digits, nearest, ties to even (0 stays \"0\")",
    "u_r = (u + |value - r|) rounded up to 40 significant digits",
    "d_abs = d rounded down to 40 significant digits; the exact d is used in the comparison",
    "if r == 0: d_rel = null; the key is no_reference unless u <= d/4 and u_r <= d/4",
    "if r != 0: d_rel = max(8*u/|r|, epsrel_exact) rounded down; the key is no_reference when this value or "
    "max(8*u_r/|value|, epsrel_exact) exceeds 1e-4",
    "substitute once per run; record the hex of m_hat and value, u_r, d_abs, d_rel in the run ledger",
]


def evaluate_fixed_m(rec, m_hat_hex):
    """value, u_r, d_abs, d_rel, status of one parametric record for a reported m_hat."""
    m = F(float.fromhex(m_hat_hex))
    if m <= 0:
        raise ValueError("m_hat must be positive")
    t = rec["terms"]
    q = m / F(rec["parameter"]["reference_length_A"])
    w = F(t["divide_by"])
    r = (F(t["constant"]["value"]) + q * F(t["proportional"]["value"])) / w
    u = (F(t["constant"]["u_r"]) + q * F(t["proportional"]["u_r"])) / w
    d = (F(t["d_abs"]["constant_exact"]) + m * F(t["d_abs"]["per_m_hat_exact"])) / w
    out = finalize(r, r - u, r + u, d, F(t["epsrel_exact"]))
    out.pop("interval")
    return out


# ----------------------------------------------------------------------------
# build
# ----------------------------------------------------------------------------

def load():
    for path, sha in ((RB_PATH, RB_SHA), (M4A_STATES_PATH, M4A_STATES_SHA), (STATES_PATH, STATES_SHA)):
        got = sha_file(path)
        if got != sha:
            raise SystemExit("frozen input changed: %s has SHA-256 %s, expected %s" % (path.name, got, sha))
    rb = json.loads(RB_PATH.read_text(encoding="utf-8"))
    states = json.loads(STATES_PATH.read_text(encoding="utf-8"))
    return rb, states


def self_checks(sources, states, records):
    """Checks of the mapping against the frozen source (no observed values are involved)."""
    out = {}
    # 1. with m = L (q = 1) the floors and d_rel of 4.2/4.4 reproduce the source table
    n = bad = 0
    for sid, src in sources.items():
        for key, k in src.row["keys"].items():
            if k.get("kind") != "direct_reference":
                continue
            sp = spec_for(src, key)
            assert not isinstance(sp, str), sp
            if key.endswith(".spectrum"):
                continue  # mapped from zero_loss and single_loss_total, checked below
            r, u = src.get(key)
            d_abs = sp["floor_const"] + src.L * sp["floor_per_m"]
            fin = finalize(r, r - u, r + u, d_abs, src.epsrel)
            n += 1
            ok = dec(d_abs, "n") == k["d_abs"] and fin["status"] == k["status"]
            if r != 0:
                ok = ok and dec(max(8 * u / abs(r), src.epsrel), "n") == k["d_rel"]
            else:
                ok = ok and k["d_rel"] is None
            bad += 0 if ok else 1
    out["source_floor_and_d_rel_reproduced_at_q_1"] = {"checked": n, "mismatch": bad}
    # 2. the source spectrum center equals r0 + r1 and its u_r is covered by u0 + u1 (+ decimal rounding)
    n = bad = 0
    worst = F(0)
    for sid, src in sources.items():
        for key, k in src.row["keys"].items():
            if not key.endswith(".spectrum"):
                continue
            region = key[:-len(".spectrum")]
            r, u = src.get(key)
            r0, u0 = src.get(region + ".zero_loss")
            r1, u1 = src.get(region + ".single_loss_total")
            n += 1
            diff = abs(r - (r0 + r1))
            worst = max(worst, diff)
            if diff > u + u0 + u1:
                bad += 1
    out["source_spectrum_equals_zero_loss_plus_single_loss_total"] = {
        "checked": n, "outside_combined_u_r": bad, "max_center_difference": dec(worst, "u")}
    # 3. the inverse interval of window_area.total reproduces the source interval
    #    (the source widens its ends by a relative 1e-37 before rounding, so it contains the tight interval)
    n = bad = 0
    worst = F(0)
    for sid in RB_IDS:
        src = sources[sid]
        r, u = src.get("window_area.total")
        lo, hi = 1 / (r + u), 1 / (r - u)
        for end in ("lower", "upper"):
            n += 1
            s_lo, s_hi = (F(x) for x in src.row["keys"]["inverse_window_area_interval_A." + end]["reference_interval_A"])
            rel = max((lo - s_lo) / lo, (s_hi - hi) / hi)
            worst = max(worst, rel)
            if not (s_lo <= lo and hi <= s_hi and rel <= F(1, 10 ** 36)):
                bad += 1
    out["source_inverse_interval_contains_exact_interval"] = {
        "checked": n, "mismatch": bad, "max_relative_excess_of_source": dec(worst, "u")}
    # 4. 13a areas inside one V5 condition group have overlapping reference intervals
    groups = []
    for group in states["reference_inventory"]["v5_area_groups"]:
        ivs = [sources[g].get("window_area.term_13a") for g in group]
        lo = max(r - u for r, u in ivs)
        hi = min(r + u for r, u in ivs)
        groups.append({"rows": group, "intervals_overlap": lo <= hi,
                       "identical_strings": len({(sources[g].row["keys"]["window_area.term_13a"]["r"],
                                                  sources[g].row["keys"]["window_area.term_13a"]["u_r"])
                                                 for g in group}) == 1})
    out["v5_area_groups"] = groups
    # 5. the two values named in prereg 4.2
    by = {(r["row_id"], r["check_kind"], r["key"]): r for r in records}
    a = by[("V3-RB1", "external", "bin[6].spectrum")]
    b = by[("V3-RB3", "external", "outside.high.spectrum")]
    out["prereg_4.2_named_values"] = {
        "V3-RB1 bin[6].spectrum": a["value"], "V3-RB1 bin[6].single_loss_total":
            by[("V3-RB1", "external", "bin[6].single_loss_total")]["value"],
        "V3-RB3 outside.high.spectrum": b["value"], "V3-RB3 outside.high.single_loss_total":
            by[("V3-RB3", "external", "outside.high.single_loss_total")]["value"],
        "as_stated": F(a["value"]) == 1 and abs(F(b["value"]) - F(35, 100)) < F(1, 10 ** 15),
    }
    # 6. unique primary keys
    out["primary_key_unique"] = len(by) == len(records)
    return out


def build():
    rb, states = load()
    sources = {sid: Source(sid, row) for sid, row in rb["rows"].items()}
    rows = {r["id"]: r for r in states["rows"]}
    inventory = states["reference_inventory"]["rb_keys_per_map"]
    records, planned, unmapped = [], {}, {}

    def plan(row_id, kind, cls):
        planned.setdefault(row_id, {}).setdefault(kind, {}).setdefault(cls, 0)
        planned[row_id][kind][cls] += 1

    for row_id in states["execution_order"]:
        row = rows[row_id]
        rule = row.get("rule")
        if rule not in ("external_map", "s2_map", "seed_reference", "seed_intersection"):
            continue
        fixture = row["inputs"]["fixture"]
        if fixture not in sources:
            raise SystemExit("row %s refers to unknown source row %s" % (row_id, fixture))
        src = sources[fixture]
        if rule in ("external_map", "s2_map"):
            if rule == "external_map":
                assert row["inputs"]["m_over_Lref"] == EXTERNAL_Q
            kinds = ["external"] if rule == "external_map" else ["fixed_m", "ratio"]
            skip = UNMAPPED.get(rule, {})
            for key in row["keys"]:
                if key in skip:
                    unmapped.setdefault(row_id, []).append(key)
            for kind in kinds:
                for key in row["keys"]:
                    if key in skip:
                        continue
                    plan(row_id, kind, key_class(key))
                    records.append(mapped_record(row_id, kind, key, src))
        else:
            for key in row["keys"]:
                if key in UNMAPPED["seed_reference"]:
                    unmapped.setdefault(row_id, []).append(key)
            recs = seed_records(row, src)
            for rec in recs:
                plan(row_id, rec["check_kind"], rec["check_kind"])
            records.extend(recs)

    # planned totals against the frozen inventory (prereg 4.4(6))
    totals = {}
    for row_id, kinds in planned.items():
        for kind, classes in kinds.items():
            t = totals.setdefault(kind, {})
            for cls, n in classes.items():
                t[cls] = t.get(cls, 0) + n
    for kind in ("external", "fixed_m", "ratio"):
        if totals.get(kind) != inventory:
            raise SystemExit("planned keys of %s %s differ from the frozen inventory %s" % (kind, totals.get(kind), inventory))

    produced, produced_totals, no_ref, path_fail = {}, {}, [], []
    for rec in records:
        p = produced.setdefault(rec["row_id"], {}).setdefault(rec["check_kind"], {})
        p[rec["status"]] = p.get(rec["status"], 0) + 1
        t = produced_totals.setdefault(rec["check_kind"], {})
        t[rec["status"]] = t.get(rec["status"], 0) + 1
        ident = {"row_id": rec["row_id"], "check_kind": rec["check_kind"], "key": rec["key"], "reason": rec["reason"]}
        if rec["status"] == "no_reference":
            no_ref.append(ident)
        elif rec["status"] == "reference_path_failure":
            path_fail.append(ident)
    for row_id, kinds in planned.items():
        for kind, classes in kinds.items():
            if sum(classes.values()) != sum(produced[row_id][kind].values()):
                raise SystemExit("planned and produced counts differ in %s %s" % (row_id, kind))

    checks = self_checks(sources, states, records)
    metadata = {
        "schema": SCHEMA,
        "rules_document": "docs/M4A2_PREREG.md (4.1-4.5), frozen at commit 4f53580",
        "inputs": {
            "benchmarks/numerics/m4a/references_rb_m4a.json": RB_SHA,
            "benchmarks/numerics/m4a/expected_states_m4a.json": M4A_STATES_SHA,
            "benchmarks/numerics/m4a2/expected_states_m4a2.json": STATES_SHA,
        },
        "generator": {"file": "benchmarks/numerics/m4a2/generate_references_m4a2.py",
                      "sha256_lf_normalised": sha_generator(),
                      "imports": "standard library only; the package under test is not imported"},
        "primary_key": ["row_id", "check_kind", "key"],
        "precision": {
            "arithmetic": "exact rationals; source decimals and binary64 inputs are read exactly",
            "stored_significant_digits": DIGITS,
            "value": "nearest, ties to even; exact zero is \"0\"",
            "interval": "outward (lower end down, upper end up); interval_exact holds the exact ends",
            "u_r": "max distance from the stored value to the exact interval ends, rounded up "
                   "(covers the source interval and the decimal rounding of the center; an exact center adds nothing)",
            "d_abs": "rounded down; d_abs_exact holds the exact rational that the comparison uses",
            "d_rel": "max(8*u/|r|, epsrel) with the exact mapped u and r, rounded down; epsrel = Fr(1e-6); "
                     "null when r = 0; the 1e-4 limit is tested on the exact value and on the stored value and u_r",
            "floors": {"B_A0": "Fr(area_epsabs) of the source row; 1e-18 (exact) for CE rows",
                       "eta": "Fr(spectrum_epsabs_fraction)", "S_abs": "sum of |Fr(intensity)| of the primary",
                       "1e-12": "exact decimal", "a": "Fr(lower window end)", "dE": "exact difference of the bin edges",
                       "m": "external: 2*L_ref; fixed_m: Fr(m_hat); ratio: m_ref = 1/A_total (center)"},
        },
        "check_kinds": {
            "external": "V3 (4.2): m = 2*L_ref, q = 2",
            "fixed_m": "V4 (1) (4.2): q = m_hat/L_ref; m-dependent keys are parametric (terms + fixed_m_evaluation)",
            "ratio": "V4 (2) (4.3): m_ref = 1/A_total of the same source row; four-corner interval; "
                     "single_loss_weight is the identity of 4.4(5)",
            "seed_area": "V5 (4.5(1)): window_area.term_13a of the source row in area units",
            "inverse_interval": "V5 (4.5(3)): exact reference interval [1/(r+u_r+c), 1/(r-u_r+c)]",
        },
        "status_values": ["value", "parametric", "no_reference", "reference_path_failure"],
        "fixed_m_evaluation": FIXED_M_EVALUATION,
        "planned_keys": {"by_row_and_kind": planned, "by_kind": totals,
                         "frozen_inventory_rb_keys_per_map": inventory},
        "produced_records": {"total": len(records), "by_row_and_kind": produced, "by_kind_and_status": produced_totals},
        "planned_keys_without_mapped_reference": {"reasons": UNMAPPED, "by_row": unmapped},
        "no_reference": no_ref,
        "reference_path_failures": path_fail,
        "self_checks": checks,
        "file_hash_note": "the SHA-256 of this file is recorded separately before any comparison (prereg 4.1)",
    }
    return metadata, records


def render(metadata, records):
    lines = ["{", '"metadata": ' + json.dumps(metadata, indent=1, sort_keys=True, ensure_ascii=True) + ",",
             '"records": [']
    for i, rec in enumerate(records):
        lines.append(json.dumps(rec, sort_keys=True, ensure_ascii=True) + ("," if i + 1 < len(records) else ""))
    lines += ["]", "}"]
    return "\n".join(lines) + "\n"


def summary(metadata):
    p = metadata["produced_records"]
    print("records: %d" % p["total"])
    for kind in sorted(p["by_kind_and_status"]):
        print("  %-17s %s" % (kind, json.dumps(p["by_kind_and_status"][kind], sort_keys=True)))
    for kind in sorted(metadata["planned_keys"]["by_kind"]):
        print("  planned %-17s %s" % (kind, json.dumps(metadata["planned_keys"]["by_kind"][kind], sort_keys=True)))
    print("no_reference: %d, reference_path_failure: %d" % (
        len(metadata["no_reference"]), len(metadata["reference_path_failures"])))
    print("self_checks: " + json.dumps({k: v for k, v in metadata["self_checks"].items() if k != "v5_area_groups"},
                                       sort_keys=True))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--check", action="store_true", help="regenerate and compare with the stored file")
    ap.add_argument("--evaluate-fixed-m", nargs=2, metavar=("ROW_ID", "M_HAT_HEX"),
                    help="evaluate the fixed_m records of one V4 row for a reported m_hat (hex)")
    args = ap.parse_args(argv)
    metadata, records = build()
    text = render(metadata, records)
    json.loads(text)
    if args.evaluate_fixed_m:
        row_id, m_hex = args.evaluate_fixed_m
        hits = [r for r in records if r["row_id"] == row_id and r["check_kind"] == "fixed_m"]
        if not hits:
            print("no fixed_m records for row %s" % row_id)
            return 2
        for rec in hits:
            out = evaluate_fixed_m(rec, m_hex) if rec["status"] == "parametric" else {
                k: rec[k] for k in ("value", "u_r", "d_abs", "d_abs_exact", "d_rel", "status", "reason")}
            print(json.dumps(dict(out, row_id=row_id, key=rec["key"], m_hat_hex=m_hex), sort_keys=True))
        return 0
    data = text.encode("utf-8")
    if args.check:
        if not OUT_PATH.exists():
            print("check: FAILED (%s does not exist)" % OUT_PATH.name)
            return 1
        stored = OUT_PATH.read_bytes()
        same = stored == data
        print("check: %s (%s, sha256 %s)" % ("ok" if same else "FAILED: stored file differs from a fresh build",
                                               OUT_PATH.name, hashlib.sha256(stored).hexdigest()))
        summary(metadata)
        return 0 if same else 1
    with open(OUT_PATH, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(text)
    print("wrote %s (%d bytes, sha256 %s)" % (OUT_PATH.name, len(data), hashlib.sha256(data).hexdigest()))
    summary(metadata)
    return 0


if __name__ == "__main__":
    sys.exit(main())
