"""Frozen M4a-2 rows, actual API calls and the per-key execution ledger."""

from contextlib import contextmanager
from copy import deepcopy
from dataclasses import is_dataclass
from fractions import Fraction as F
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
assert hashlib.sha256(RAW).hexdigest() == "4a8264a58b9902efee3d843c10be812569ff6f46d720c9f23ebf37a4f263073f"
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

    def status(self, key):
        statuses = [v["status"] for (k, name), v in self.items.items() if name == key]
        return next(
            (
                s
                for s in ("fail", "unresolved", "reference_path_failure", "no_reference", "not_run", "pass")
                if s in statuses
            ),
            "not_run",
        )


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


@lru_cache(None)
def references():
    path = DIRECTORY / "references_m4a2.json"
    if not path.exists():
        return None
    raw = path.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    doc = json.loads(raw)
    index = {}
    for r in doc["records"]:
        key = (r["row_id"], r["check_kind"], r["key"])
        if key in index:
            raise ValueError(f"duplicate mapped reference {key}")
        index[key] = r
    return digest, index


@lru_cache(None)
def rb_sources():
    path = DIRECTORY.parent / "m4a/references_rb_m4a.json"
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != TABLE["frozen_inputs"][path.name]:
        raise ValueError("R-B source hash differs from the registration")
    return json.loads(raw)["rows"]


def fixed_reference(record, key, multiplier):
    """Substitute one binary64 m in the affine expression of preregistration 4.1-4.2."""
    p = record["parameter"]
    if isinstance(p, dict) and "proportional" in p:
        q = F(multiplier) / F(p["L_ref"])
        r = F(p["constant"]["value"]) + q * F(p["proportional"]["value"])
        u = F(p["constant"]["u_r"]) + q * F(p["proportional"]["u_r"])
        floor = F(p["d_abs_constant"]) + q * F(p["d_abs_proportional"])
        epsrel = F(p["epsrel"])
    else:
        # A symbolic parameter can identify the frozen operands by source row
        # and keys instead of duplicating their decimal coefficients.
        source = rb_sources()[record["source_row"]]
        conditions = source["conditions"]
        q = F(multiplier) / F(float.fromhex(conditions["weight_length_A"]))
        epsrel = F(conditions["comparison_settings"]["epsrel"])
        raw_key = key.removesuffix("_per_eV")
        width = F(1)
        if key != raw_key:
            n = int(key[4 : key.index("]")])
            edges = conditions["output_edges_eV"]
            width = F(float.fromhex(edges[n + 1])) - F(float.fromhex(edges[n]))
        terms = [(raw_key, q)]
        if raw_key.startswith("window_area.") or raw_key.endswith(".zero_loss"):
            terms = [(raw_key, F(1))]
        elif raw_key.endswith(".spectrum"):
            region = raw_key.removesuffix(".spectrum")
            terms = [(region + ".zero_loss", F(1)), (region + ".single_loss_total", q)]
        refs = source["keys"]
        r = sum(factor * F(refs[k]["r"]) for k, factor in terms) / width
        u = sum(factor * F(refs[k]["u_r"]) for k, factor in terms) / width
        scale_floor = raw_key == "single_loss_weight" or "cutoff_sensitivity." in raw_key
        floor = F(refs[raw_key]["d_abs"]) * (q if scale_floor else 1) / width
    rel = max(epsrel, 8 * u / abs(r)) if r else None
    return dict(value=r, u_r=u, d_abs=floor, d_rel=rel)


def judge_numeric(value, error, ref, budget=None):
    r, u, floor = (F(ref[k]) for k in ("value", "u_r", "d_abs"))
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


def mapped(log, kind, key, value, error=None, *, interval=None, budget=None, multiplier=None):
    table = references()
    if table is None:
        log.record(kind, key, "not_run", reason="references_absent")
        return
    log.reference_sha256, refs = table
    rec = refs.get((log.row["id"], kind, key))
    if rec is None:
        log.record(kind, key, "fail", reason="missing_planned_reference_record")
        return
    if rec["status"] in ("no_reference", "reference_path_failure"):
        log.record(kind, key, rec["status"], reason=rec["reason"])
        return
    if interval is not None:
        # Test containment against the exact reciprocal of the source
        # interval; a 40-digit outward-rounded endpoint is a display value.
        source = rb_sources()[rec["source_row"]]["keys"]
        seed = kind == "inverse_interval" and log.row["rule"] == "seed_reference"
        area = source["window_area.term_13a" if seed else "window_area.total"]
        d = F(area["r"]) + (F(log.row["inputs"]["c"]) if seed else 0)
        u = F(area["u_r"])
        if d <= u:
            log.record(kind, key, "no_reference", reason="reference_denominator")
            return
        lower, upper = 1 / (d + u), 1 / (d - u)
        bounds = rec["interval"]
        if isinstance(bounds, dict):
            bounds = (bounds["lower"], bounds["upper"])
        stored_lo, stored_hi = map(F, bounds)
        if not stored_lo <= lower <= upper <= stored_hi:
            log.record(kind, key, "fail", reason="mapped_interval_rounding")
            return
        a, b = map(F, interval)
        status = "pass" if a <= lower <= upper <= b else "fail" if b < lower or a > upper else "unresolved"
        log.record(kind, key, status, coverage="interval")
        return
    extra = {}
    if rec["status"] == "parametric":
        rec = fixed_reference(rec, key, multiplier)
        extra = dict(multiplier_hex=float(multiplier).hex(), evaluated={k: str(v) for k, v in rec.items()})
    log.record(
        kind,
        key,
        **judge_numeric(value, error, rec, budget),
        observed=float(value).hex(),
        error=float(error).hex(),
        **extra,
    )


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


def seed_test(row, log):
    i = row["inputs"]
    args, kw = fixture(i["fixture"])
    args = (args[0], args[1], kw.pop("loss_window_eV"), kw.pop("epsilon_momentum"), kw.pop("hole_charge"))
    del kw["weight_length_A"]

    def run(s):
        return cached_seed(*args, tuple(sorted(kw.items())), s, i["c"])

    if "seeds_A" in i:
        intervals = [run(s)["interval_A"] for s in i["seeds_A"]]
        members = [row["id"].removesuffix("intersection") + str(s) for s in i["seeds_A"]]
        log.check("common_interval", max(x[0] for x in intervals) <= min(x[1] for x in intervals),
                  member_rows=members)
        if references() is not None:
            for member in members:
                single = Ledger(ROWS[member])
                seed_test(ROWS[member], single)
                statuses = {v["status"] for v in single.items.values()}
                outcome = next(s for s in ("fail", "unresolved", "reference_path_failure", "no_reference",
                                            "not_run", "pass") if s in statuses)
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
                mapped(log, kind, key, None, interval=interval)
            elif key == "single_loss_weight" and kind == "ratio":
                log.equal(key, r.single_loss_weight, float(F(center) * F(r.window_area["total"])), kind)
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


def execute(row):
    log = Ledger(row)
    try:
        rule = row["rule"]
        if rule == "m4a_ledger":
            inherited(row, log)
        elif rule in ("beta_replay", "area_invariance"):
            beta(row, log)
        elif rule in ("seed_reference", "seed_intersection"):
            seed_test(row, log)
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
