"""Frozen M4a-2 rows, actual API calls and the per-key execution ledger."""

from collections import Counter
from contextlib import contextmanager
from copy import deepcopy
from dataclasses import dataclass, is_dataclass
from decimal import Context, Decimal, ROUND_CEILING, ROUND_FLOOR, ROUND_HALF_EVEN
from fractions import Fraction as F
from fnmatch import fnmatchcase
from functools import lru_cache
import hashlib
import json
import math
import os
from pathlib import Path
import traceback
from unittest.mock import patch

import numpy as np
import m4a_states as old
from xpsloss import xps_single_loss as api
from xpsloss import _fold, _loss_panels as panels
from xpsloss._single_loss_types import COMPONENTS, QUANTITIES, EPS

ROOT = Path(__file__).resolve().parents[1]
DIRECTORY = ROOT / "benchmarks/numerics/m4a2"
RAW = (DIRECTORY / "expected_states_m4a2.json").read_bytes()
STATES_SHA256 = "4a8264a58b9902efee3d843c10be812569ff6f46d720c9f23ebf37a4f263073f"
assert hashlib.sha256(RAW).hexdigest() == STATES_SHA256
TABLE = json.loads(RAW)
ROWS = {r["id"]: r for r in TABLE["rows"]}
DERIVED = dict(weight_length_kind="inverse_13a_window_area", weight_length_A=None)
S2 = dict(multiplier="inverse_window_area")
NEW = {
    "model_id",
    "weight_length_kind",
    "weight_length_interval_A",
    "extra_inverse_length_per_A",
    "multiplier_convention",
    "multiplier_A",
    "multiplier_interval_A",
}


def bits(v):
    if is_dataclass(v):
        return (type(v).__name__, bits(vars(v)))
    if isinstance(v, dict):
        return {k: bits(x) for k, x in v.items()}
    if isinstance(v, np.ndarray):
        return (v.dtype.str, v.shape, v.tobytes())
    if isinstance(v, (list, tuple)):
        return (type(v).__name__, [bits(x) for x in v])
    if isinstance(v, (float, np.floating)):
        return type(v).__name__, float(v).hex()
    return type(v).__name__, v


class Ledger:
    def __init__(self, row):
        self.row = row
        self.items = {}

    def record(self, kind, key, status, **details):
        self.items[kind, key] = dict(status=status, **details)

    def check(self, key, ok, kind="identity", **details):
        self.record(kind, key, "pass" if ok else "fail", **details)

    def equal(self, key, value, target, kind="identity"):
        self.check(key, bits(value) == bits(target), kind)

    def status(self, key=None):
        statuses = {v["status"] for (_, name), v in self.items.items() if key is None or name == key}
        return next(
            (
                s
                for s in ("fail", "unresolved", "reference_path_failure", "no_reference", "not_run", "pass")
                if s in statuses
            ),
            "not_run",
        )


# A registered key with the same spelling as a check maps directly. These
# groups name the finer ledger checks used for the other registered keys.
# Patterns address "kind/key"; no execution-completed marker is evidence.
REGION_VALUES = tuple("identity/*." + q for q in QUANTITIES)
SENSITIVITIES = {
    "w": ("identity/single_loss_weight",),
    "H": ("identity/cutoff_sensitivity.*",),
    "Hpm": ("identity/*.cutoff_sensitivity.*",),
    "LbK_b": ("identity/diagnostics.LbK_b.*",),
}
REGISTERED_CHECKS = {
    "paired_baseline": {"legacy": ("identity/fields", "identity/diagnostics", "identity/legacy.keys",
                                   "identity/manifest.*", "identity/files.*")},
    "m4a_ledger": {"original_check_plan": ("old.*/*",)},
    "beta_replay": {"beta_numeric": ("identity/fixed.*",), "beta_diag": ("identity/beta.*",),
                    "returned_errors": ("rounding/*_error",)},
    "area_invariance": {"area11": ("identity/window_area",), "area_error11": ("identity/window_area_error",),
                        "area_budget11": ("identity/window_area_budget",)},
    "s2_map": {"returned_errors": ("rounding/*_error",)},
    "analytic_map": {**SENSITIVITIES, "regions5": REGION_VALUES, "densities": ("identity/*_per_eV",)},
    "direct_tails": {"outside.low": ("identity/outside.low.*",), "outside.high": ("identity/outside.high.*",)},
    "refold": {"budget": ("identity/*.budget",)},
    "error_formula": {
        "spectrum_error": ("identity/*.spectrum_error",), "density_error": ("identity/*_per_eV_error",),
        "satellite_errors": ("rounding/bin*.single_loss_*_error",), "outside_errors": ("rounding/outside.*_error",),
        "H_errors": ("rounding/cutoff_sensitivity.*_error",), "Hpm_errors": ("rounding/*.cutoff_sensitivity.*_error",),
    },
    "wiring": {"kernel_calls": ("identity/*.call_kwargs.weight_length_A",)},
    "partial_sensitivity": {**SENSITIVITIES, "balance": ("identity/diagnostics.balance.*",)},
    "structural_zero": {"kernel_calls": ("identity/calls",), "interval": ("identity/weight_length_interval_A",),
                        "w": ("identity/J",), "J0": ("identity/J",)},
    "zero_primary": {"w_error": ("identity/w",), "regions5": REGION_VALUES},
    "signed_s2": {"single_loss_hole": ("identity/*.single_loss_hole",)},
    "denominator_progress": {"stage": ("identity/denominator", "identity/final_panels")},
    "old_failure": {"worst_or_geometry": ("identity/worst_key", "identity/cells", "identity/panel_eV")},
    "convention_record": {
        "stdout": ("identity/stdout.model_id",), "manifest": ("identity/manifest.model_id", "identity/replay.spectrum_definition"),
        "csv": ("identity/replay.csv",), "replay": ("identity/replay.*",),
        "comparison": ("identity/comparison.*",), "units": ("identity/units.*",),
    },
    "old_manifest": {"legacy": ("identity/legacy.keys", "identity/window_area.*"), "replay": ("identity/csv",)},
    "old_constructor": {"legacy": ("identity/field_order", "identity/old_args")},
    "capabilities": {"supported": ("identity/other_blocks",)},
    "kernel_failure": {"stdout": ("identity/exit",), "outdir": ("identity/exit",)},
}


def registered_checks(log):
    """Resolve each frozen key to actual checks. Keep absent evidence visible."""
    aliases = REGISTERED_CHECKS.get(log.row["rule"], {})
    by_key = {}
    for kind, key in log.items:
        if kind not in ("execution", "registration", "unavailable"):
            by_key.setdefault(key, []).append(f"{kind}/{key}")
    checks = [check for matches in by_key.values() for check in matches]
    mapping = {}
    for key in log.row["keys"]:
        patterns = ("identity/exception.*",) if key == "exception" else aliases.get(key, ())
        matches = set(by_key.get(key, ()))
        matches.update(check for pattern in patterns for check in checks if fnmatchcase(check, pattern))
        timed_out = log.items.get(("unavailable", key), {}).get("reason") == "time_limit"
        if not matches and not timed_out:
            log.record("registration", key, "not_run", reason="unmapped_registered_key")
        mapping[key] = sorted(matches)
    return mapping


def fixture(name="short", spectrum=False):
    if name in old.rows():
        row = old.rows()[name]
        c, i = old.conditions(row), row["inputs"]
        kw = dict(
            loss_window_eV=c["window"],
            weight_length_A=c["L"],
            epsilon_momentum=c["momentum"],
            hole_charge=c["charge"],
            **row["settings"],
        )
        mat, energy = old.material(i["material"]), c["E0"]
    else:
        mat, energy = old.material("synthetic_wide"), 1000.0
        kw = dict(loss_window_eV=(5.0, 6.0), weight_length_A=20.0, epsilon_momentum="k_total")
        i = {}
    args = [mat, energy]
    if spectrum:
        args.extend(
            [old.make_primary(i.get("primary", "D1000")), old.make_edges(i.get("output_edges_eV", "E_RB1"))]
        )
    return args, kw


def call(args, kw):
    fn = api.single_loss_spectrum if len(args) == 4 else api.loss_window_area
    return fn(*args, **kw)


def rounded(q, direction):
    n = float(q)
    if (direction > 0 and n < q) or (direction < 0 and n > q):
        return math.nextafter(n, direction)
    return n


def reciprocal(a, e, c=0):
    d, err = F(a) + F(c), F(e)
    return float(1 / d), (rounded(1 / (d + err), -math.inf), rounded(1 / (d - err), math.inf))


@contextmanager
def fixed_errors():
    observations = []
    propagate = _fold.propagate_multiplier

    def observe(result):
        observations.append(deepcopy(result))
        propagate(result)

    with patch.object(_fold, "propagate_multiplier", observe):
        yield observations


def check_errors(log, result, fixed, *, check_balance=True):
    m = F(result.multiplier_A)
    lo, hi = map(F, result.multiplier_interval_A)
    rho = rounded(max(m - lo, hi - m) / m, math.inf)
    log.equal("rho", result.diagnostics["multiplier"]["rho"], rho)
    targets = list(("cutoff_sensitivity." + c) for c in old.COLUMNS)
    regions = old.region_names(result) if hasattr(result, "outside") else []
    for region in regions:
        targets.extend(region + ".single_loss_" + c for c in COMPONENTS)
        targets.extend(f"{region}.cutoff_sensitivity.{s}.{c}" for s in ("plus", "minus") for c in COMPONENTS)
    if result.multiplier_convention == "weight_length":
        targets.append("single_loss_weight")
    else:
        log.equal("fixed_w", result.single_loss_weight, float(m * F(result.window_area["total"])))
        log.equal(
            "w_error",
            result.single_loss_weight_error,
            rounded(abs(F(result.single_loss_weight) - 1), math.inf),
        )
        if regions and check_balance:
            mass = F(result.primary_integral)
            total = sum(F(old.read_key(result, r + ".single_loss_total")[0]) for r in regions)
            bound = F(result.diagnostics["balance"]["single_loss_total"]["rounding_bound"])
            log.check(
                "normalized_balance",
                abs(total - mass) <= bound + mass * abs(m * F(result.window_area["total"]) - 1),
            )
    for key in targets:
        y, ef = map(F, old.read_key(fixed, key))
        value, err = old.read_key(result, key)
        want = rounded(ef + F(rho) * (abs(y) + ef), math.inf)
        corners = [p * q / m for p in (y - ef, y + ef) for q in (lo, hi)]
        log.check(
            key + "_error",
            float(err).hex() == want.hex()
            and F(value) == y
            and y - F(err) <= min(corners) <= max(corners) <= y + F(err),
            "rounding",
        )
    for r in regions:
        y, e = old.read_key(result, r + ".spectrum")
        _, e0 = old.read_key(result, r + ".zero_loss")
        _, e1 = old.read_key(result, r + ".single_loss_total")
        log.equal(r + ".zero_loss_error", float(e0), float(old.read_key(fixed, r + ".zero_loss")[1]))
        log.equal(r + ".spectrum_error", float(e), rounded(F(e0) + F(e1) + F(EPS) * abs(F(y)), math.inf))
        if r.startswith("bin"):
            i = int(r[4:-1])
            width = F(result.output_edges_eV[i + 1]) - F(result.output_edges_eV[i])
            for k in QUANTITIES:
                y, e = old.read_key(result, r + "." + k)
                log.equal(
                    r + "." + k + "_per_eV_error",
                    float(getattr(result, k + "_per_eV_error")[i]),
                    float((F(e) + F(EPS) * abs(F(y))) / width),
                )
    log.equal("window_area_error", result.window_area_error, fixed.window_area_error)


def beta(row, log):
    i = row["inputs"]
    args, kw = fixture(i["fixture"], row["entry"] == "api_spectrum")
    if row["rule"] == "area_invariance":
        a, b = call(args, kw | dict(multiplier="external", multiplier_A=7.0)), call(args, kw)
        log.equal("window_area.total", a.window_area["total"], b.window_area["total"])
        for k in ("window_area", "window_area_error", "window_area_budget"):
            log.equal(k, getattr(a, k), getattr(b, k))
        log.equal("area_parts11", a.diagnostics["error_budget"], b.diagnostics["error_budget"])
        return
    kw["multiplier"] = i["multiplier"]
    if kw["multiplier"] == "external":
        kw["multiplier_A"] = i["external_m_A"]
    with fixed_errors() as observed:
        r = call(args, kw | DERIVED | dict(extra_inverse_length_per_A=i["c"]))
    with fixed_errors() as control:
        b = call(args, kw | dict(weight_length_A=r.weight_length_A))
    if observed:
        check_errors(log, r, observed[0])
    lhs, rhs = (observed[0], control[0] if control else b) if observed else (r, b)
    for k in vars(lhs).keys() - NEW - {"diagnostics", "warnings"}:
        log.equal("fixed." + k, getattr(lhs, k), getattr(rhs, k))
        if kw["multiplier"] != "weight_length" or "error" not in k and k != "outside":
            log.equal(k, getattr(r, k), getattr(b, k))
    log.equal("cutoff_sensitivity.total", r.cutoff_sensitivity["total"], b.cutoff_sensitivity["total"])
    for k in b.diagnostics.keys() - {"comparison", "warnings", "multiplier"}:
        log.equal("beta." + k, r.diagnostics[k], b.diagnostics[k])


@dataclass
class ReferenceTable:
    digest: str
    records: dict

    def audit_usage(self, rows, used):
        """Count assigned comparisons only. Mutation witnesses are separate executions."""
        expected = Counter(k for k in self.records if k[0] in rows)
        if used != expected:
            raise ValueError(f"reference usage: missing {expected - used}, unexpected/repeated {used - expected}")


def index_references(doc):
    meta = doc["metadata"]
    paths = {DIRECTORY.parent / "m4a" / name: sha for name, sha in TABLE["frozen_inputs"].items()}
    paths[DIRECTORY / "expected_states_m4a2.json"] = STATES_SHA256
    if meta["inputs"] != {p.relative_to(ROOT).as_posix(): sha for p, sha in paths.items()}:
        raise ValueError("mapped reference input hashes differ from the registration")
    for p, sha in paths.items():
        if hashlib.sha256(p.read_bytes()).hexdigest() != sha:
            raise ValueError(f"frozen input hash mismatch: {p.name}")

    index = {}
    by_row, by_kind = Counter(), Counter()
    for rec in doc["records"]:
        rid, kind, key = pk = (rec["row_id"], rec["check_kind"], rec["key"])
        if pk in index:
            raise ValueError(f"duplicate mapped reference {pk}")
        if rec["source_sha256"] != TABLE["frozen_inputs"]["references_rb_m4a.json"]:
            raise ValueError(f"mapped reference source hash mismatch: {pk}")
        if rid not in ROWS or key not in ROWS[rid]["keys"] or rec["source_row"] != ROWS[rid]["inputs"]["fixture"]:
            raise ValueError(f"unregistered reference key or source row: {pk}")
        if kind not in meta["check_kinds"]:
            raise ValueError(f"unregistered reference kind: {pk}")
        category = rec["counted_as"]
        if kind == "ratio" and key == "single_loss_weight":
            if category != "identity_4.4(5)":
                raise ValueError(f"missing rounding identity: {pk}")
            # planned_keys includes these identities among the 1,024 direct keys.
            category = "direct_reference"
        by_row[rid, kind, category] += 1
        by_kind[kind, category] += 1
        index[pk] = rec
    plan = meta["planned_keys"]
    planned_rows = {(rid, kind, cat): n for rid, kinds in plan["by_row_and_kind"].items()
                    for kind, classes in kinds.items() for cat, n in classes.items()}
    planned_kinds = {(kind, cat): n for kind, classes in plan["by_kind"].items() for cat, n in classes.items()}
    if by_row != planned_rows or by_kind != planned_kinds:
        raise ValueError("mapped reference counts differ from planned_keys by row, kind or counted_as")
    return index


@lru_cache(None)
def references():
    path = DIRECTORY / "references_m4a2.json"
    if not path.exists():
        return None
    raw = path.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    return ReferenceTable(digest, index_references(json.loads(raw)))


def decimal40(x, rounding):
    if not x:
        return "0"
    value = Context(prec=40, rounding=rounding).divide(Decimal(x.numerator), Decimal(x.denominator))
    digits, exponent = format(value, ".39e").split("e")
    return digits if int(exponent) == 0 else f"{digits}e{int(exponent):+d}"


def fixed_reference(record, multiplier):
    """Evaluate metadata.fixed_m_evaluation using only the frozen operands."""
    if not math.isfinite(multiplier) or multiplier <= 0:
        raise ValueError("fixed_m requires a finite positive binary64 multiplier")
    m = F(multiplier)
    q = m / F(record["parameter"]["reference_length_A"])
    t = record["terms"]
    width = F(t["divide_by"])
    r = (F(t["constant"]["value"]) + q * F(t["proportional"]["value"])) / width
    u = (F(t["constant"]["u_r"]) + q * F(t["proportional"]["u_r"])) / width
    floor = (F(t["d_abs"]["constant_exact"]) + m * F(t["d_abs"]["per_m_hat_exact"])) / width
    value = decimal40(r, ROUND_HALF_EVEN)
    stored_u = decimal40(u + abs(F(value) - r), ROUND_CEILING)
    rel = max(8 * u / abs(r), F(t["epsrel_exact"])) if r else None
    if r:
        eligible = max(rel, 8 * F(stored_u) / abs(F(value))) <= F("1e-4")
    else:
        eligible = max(u, F(stored_u)) <= floor / 4
    return dict(value=value, u_r=stored_u, d_abs=decimal40(floor, ROUND_FLOOR), d_abs_exact=str(floor),
                d_rel=decimal40(rel, ROUND_FLOOR) if r else None, status="value" if eligible else "no_reference",
                reason=None if eligible else "reference_width")


def judge_numeric(value, error, ref, budget=None):
    r, u, floor = (F(ref[k]) for k in ("value", "u_r", "d_abs_exact"))
    rel = F(ref["d_rel"]) if r else 0
    if min(u, floor, rel) < 0:
        return dict(status="fail", reason="invalid_reference")
    if rel > F("0.0001") or (r and 8 * u > F("0.0001") * abs(r)) or (not r and u > floor / 4):
        return dict(status="no_reference", reason="reference_width")
    if not math.isfinite(value) or not math.isfinite(error) or error < 0:
        return dict(status="fail", reason="nonfinite_or_negative_error")
    y, e = F(value), F(error)
    delta, cover = abs(y - r), e + 4 * F(EPS) * max(abs(y), abs(r))
    precision = delta + u <= floor + rel * abs(r)
    coverage = "covered" if delta + u <= cover else "underestimate" if delta - u > cover else "partial"
    within = e <= F(budget) if budget is not None else None
    status = (
        "fail"
        if not precision or coverage == "underestimate" or within is False
        else ("pass" if coverage == "covered" else "unresolved")
    )
    return dict(status=status, precision=precision, coverage=coverage, budget=within)


def mapped(log, kind, key, value, error=None, *, interval=None, budget=None, multiplier=None, identity_value=None):
    if (kind, key) in log.items:
        raise ValueError(f"repeated comparison: {(log.row['id'], kind, key)}")
    table = references()
    if table is None:
        log.record(kind, key, "not_run", reason="references_absent")
        return
    log.reference_sha256 = table.digest
    rec = table.records.get((log.row["id"], kind, key))
    if rec is None:
        log.record(kind, key, "fail", reason="missing_planned_reference_record")
        return
    details = dict(counted_as=rec["counted_as"])
    if interval is None:
        details.update(observed=float(value).hex(), error=float(error).hex())
        if budget is not None:
            details["budget_value"] = float(budget).hex()
    else:
        details["reported_interval"] = [float(x).hex() for x in interval]
    if rec["status"] == "parametric":
        rec = fixed_reference(rec, multiplier)
        details.update(multiplier_hex=float(multiplier).hex(),
                       evaluated={k: rec[k] for k in ("value", "u_r", "d_abs", "d_rel")})
        details["d_abs_exact"] = rec["d_abs_exact"]
    if rec["status"] in ("no_reference", "reference_path_failure"):
        log.record(kind, key, rec["status"], reason=rec["reason"], **details)
        return
    if interval is not None:
        lower, upper = map(F, rec["interval_exact"])
        details["interval_exact"] = rec["interval_exact"]
        stored_lo, stored_hi = map(F, rec["interval"])
        if not stored_lo <= lower <= upper <= stored_hi:
            log.record(kind, key, "fail", reason="mapped_interval_rounding", **details)
            return
        a, b = map(F, interval)
        if a > b or b < lower or a > upper:
            status = "fail"
        elif a <= lower and upper <= b:
            status = "pass"
        else:
            status = "unresolved"
        log.record(kind, key, status, coverage="interval", **details)
        return
    if identity_value is not None:
        target = F(rec["value"])
        ok = (target == 1 and F(rec["u_r"]) == 0 and bits(value) == bits(identity_value)
              and bits(error) == bits(rounded(abs(F(value) - target), math.inf)))
        log.check(key, ok, kind, **details)
        return
    verdict = judge_numeric(value, error, rec, budget)
    if verdict["status"] != "pass":
        details["reference"] = {k: rec[k] for k in ("value", "u_r", "d_abs", "d_abs_exact", "d_rel")}
    log.record(kind, key, **verdict, **details)


def seed_observation(args, kw):
    class SeedComplete(Exception):
        pass

    store = {}
    original = api._seed_length
    raw = panels._areas

    def observe(p):
        outputs, local = raw(p)
        store["actual13"] = outputs["window_area.term_13a"].copy()
        return outputs, local

    def stop(*a):
        with patch.object(panels, "_areas", observe):
            store.update(original(*a))
        store["diag"] = a[-1]
        raise SeedComplete

    with patch.object(api, "_seed_length", stop):
        try:
            call(args, kw | DERIVED)
        except SeedComplete:
            return store
    raise AssertionError("seed entry not reached")


@lru_cache(None)
def cached_seed(mat, energy, window, mode, charge, settings, seed, c):
    return seed_observation(
        [mat, energy],
        dict(
            settings,
            loss_window_eV=window,
            epsilon_momentum=mode,
            hole_charge=charge,
            weight_length_seed_A=seed,
            extra_inverse_length_per_A=c,
        ),
    )


def seed_test(row, log, member_statuses):
    i = row["inputs"]
    args, kw = fixture(i["fixture"])
    args = (args[0], args[1], kw.pop("loss_window_eV"), kw.pop("epsilon_momentum"), kw.pop("hole_charge"))
    del kw["weight_length_A"]

    def run(s):
        return cached_seed(*args, tuple(sorted(kw.items())), s, i["c"])

    if "seeds_A" in i:
        intervals = [run(s)["interval_A"] for s in i["seeds_A"]]
        members = [row["id"].removesuffix("intersection") + str(s) for s in i["seeds_A"]]
        common = max(x[0] for x in intervals), min(x[1] for x in intervals)
        log.check("common_interval", common[0] <= common[1], member_rows=members)
        mapped(log, "inverse_interval", "common_interval", None, interval=common)
        if references() is not None:
            for member in members:
                outcome = member_statuses.get(member)
                if outcome is None:  # A standalone intersection still checks each seed.
                    outcome = execute(ROWS[member]).status()
                log.record("member", member, outcome)
        return
    v = run(i["seed_A"])
    raw = v["actual13"]
    length, interval = reciprocal(raw["value"], raw["estimate"], i["c"])
    log.equal("weight_length_A", v["length_A"], length)
    log.equal("weight_length_interval_A", v["interval_A"], interval)
    log.equal("seed.error13", v["error13"], raw["estimate"])
    log.check(
        "denominator",
        F(v["error13"]) <= F(v["diag"]["effective_settings"]["epsrel"]) * abs(F(v["area13"]) + F(i["c"]))
        and F(v["area13"]) + F(i["c"]) > F(v["error13"]),
    )
    budget = v["diag"]["error_budget"]["outputs"]["window_area.term_13a"]["budget"]
    mapped(log, "seed_area", "seed.area13", v["area13"], v["error13"], budget=budget)
    mapped(log, "inverse_interval", "weight_length_interval_A", None, interval=v["interval_A"])


def real_map(row, log):
    i = row["inputs"]
    args, kw = fixture(i["fixture"], row["entry"] == "api_spectrum")
    s2 = row["rule"] == "s2_map"
    with fixed_errors() as before:
        r = call(
            args, kw | (S2 if s2 else dict(multiplier="external", multiplier_A=2 * kw["weight_length_A"]))
        )
    center, interval = reciprocal(r.window_area["total"], r.window_area_error["total"])
    log.equal("inverse_window_area_A", r.inverse_window_area_A, center)
    log.equal("inverse_window_area_interval_A", r.inverse_window_area_interval_A, interval)
    if s2:
        log.equal("multiplier_A", r.multiplier_A, center)
        log.equal("multiplier_interval_A", r.multiplier_interval_A, interval)
        log.equal("denominator", r.diagnostics["denominator"]["value"], r.window_area["total"])
        check_errors(log, r, before[0])
    if not s2 and i["fixture"] in ("RB1", "RB3"):
        region = "bin[6]" if i["fixture"] == "RB1" else "outside.high"
        zero, _ = old.read_key(r, region + ".zero_loss")
        log.equal(region + ".spectrum", float(old.read_key(r, region + ".spectrum")[0]), float(zero))
        log.equal(
            region + ".single_loss_total", float(old.read_key(r, region + ".single_loss_total")[0]), 0.0
        )
    for kind in ("fixed_m", "ratio") if s2 else ("external",):
        obj = before[0] if kind == "fixed_m" else r
        for key in row["keys"]:
            if key in {"multiplier_A", "multiplier_interval_A", "returned_errors", "rho", "denominator"}:
                continue
            if key.startswith("inverse_window_area_interval_A."):
                mapped(log, kind, key, None, interval=r.inverse_window_area_interval_A)
            elif key == "single_loss_weight" and kind == "ratio":
                mapped(log, kind, key, r.single_loss_weight, r.single_loss_weight_error,
                       identity_value=float(F(center) * F(r.window_area["total"])))
            else:
                y, e = old.read_key(obj, key)
                rec = r.diagnostics["error_budget"]["outputs"].get(key)
                mapped(
                    log,
                    kind,
                    key,
                    y,
                    e,
                    multiplier=r.multiplier_A,
                    budget=rec["budget"] if rec and kind != "ratio" else None,
                )


def historical(row, log):
    from xpsloss.config import load_config

    i = row["inputs"]
    mat = i["fixture"].split("-")[-1].lower()
    cfg = load_config(ROOT / "benchmarks/inputs" / f"{mat}_elfdb_v1.0.bulk.json")
    args, kw = [cfg.material, cfg.electron_energy_eV], dict(epsilon_momentum="k_total")
    r = call(
        args,
        kw
        | DERIVED
        | dict(
            loss_window_eV=i["window"],
            weight_length_seed_A=i["seed_A"],
            multiplier="external",
            multiplier_A=i["external_m_A"],
        ),
    )
    for k, ref in (("weight_length_A", "reference_length_A"), ("single_loss_weight", "reference_weight")):
        target = F(i[ref])
        log.check(k, abs(F(getattr(r, k)) - target) <= F("4e-6") * abs(target))


def inherited(row, log):
    rid = row["inputs"]["fixture"]
    j = old.Judgement()
    out = old.judge_row(old.rows()[rid], old.run_row, j)
    for (_, key, kind), (status, detail) in j.items.items():
        log.record(
            "old." + kind,
            key,
            {"skipped": "not_run", "recorded": "recorded"}.get(status, status),
            detail=detail,
        )
    if rid == "N27" or (rid.startswith("RA") and rid in old.ra_references()["rows"]):
        judged = old.judge_n27(out) if rid == "N27" else old.judge_reference_row(old.rows()[rid], out)
        for key, (status, detail) in judged.items():
            log.record("old.reference", key, status, detail=detail)
    if rid.startswith(("RB", "CE")):
        import test_m4a_reference_rb as rb

        class Tally:
            def add(self, rid, key, status, detail):
                log.record("old.reference", key, status, detail=detail)

        rb.compare_result(rid, rb.row_inputs(rid), out.result, Tally())
    if "matched_display" in out.extras:
        status, detail = out.extras["matched_display"]
        log.record("old.display", "comparison", status, detail=detail)
    if rid == "N25":
        args, kw = fixture(rid)
        if (ROOT / ".git").exists():
            from test_m4a2_records import pairs

            before, _ = pairs(("V1-old-N25",))
            failure = before["V1-old-N25"]
            log.equal("legacy.failure_message", str(out.exc), failure["message"])
            log.equal(
                "legacy.failure_attributes",
                {k: v for k, v in vars(out.exc).items() if k != "phase"},
                failure["attributes"],
            )
        try:
            call(args, kw | dict(multiplier="external", multiplier_A=kw["weight_length_A"]))
        except api.SingleLossConvergenceError as exc:
            log.equal("N25.equal_external", vars(exc), vars(out.exc))
            log.check("exception.phase", exc.phase is None and out.exc.phase is None)
            log.check("exception.stage", exc.stage is None and out.exc.stage is None)
        else:
            log.check("exception.reason", False)


def execute(row, member_statuses=None):
    log = Ledger(row)
    try:
        rule = row["rule"]
        if rule == "m4a_ledger":
            inherited(row, log)
        elif rule in ("beta_replay", "area_invariance"):
            beta(row, log)
        elif rule in ("seed_reference", "seed_intersection"):
            seed_test(row, log, member_statuses if member_statuses is not None else {})
        elif rule in ("external_map", "s2_map"):
            real_map(row, log)
        elif rule == "historical_replay":
            historical(row, log)
        elif rule == "mutation":
            from test_m4a2_mutations import mutation

            mutation(row, log)
        elif rule in (
            "paired_baseline",
            "convention_record",
            "old_manifest",
            "old_constructor",
            "usage_error",
            "capabilities",
        ) or row["id"].startswith("V10-stderr"):
            from test_m4a2_records import records

            records(row, log)
        else:
            from test_m4a2_driver import driver

            driver(row, log)
        if not log.items:
            raise AssertionError("row produced no checks")
        log.record("execution", "completed_checks", "pass")
    except Exception as exc:
        log.record(
            "execution",
            "unexpected_exception",
            "fail",
            error=type(exc).__name__,
            reason=getattr(exc, "reason", None),
            traceback=traceback.format_exc(),
        )
        for key in row["keys"]:
            if log.status(key) == "not_run":
                log.record("unavailable", key, "not_run", reason="execution_failed")
    return log


def profile():
    full = os.environ.get("XPSLOSS_M4A2_FULL") == "1"
    mutations = os.environ.get("XPSLOSS_M4A2_MUTATIONS") == "1"
    only = os.environ.get("XPSLOSS_M4A2_ROWS", "").split(",")
    return [
        r
        for rid in TABLE["execution_order"]
        if (r := ROWS[rid])
        and (rid in only if only != [""] else (mutations if r["rule"] == "mutation" else full or r["ci"]))
    ]
