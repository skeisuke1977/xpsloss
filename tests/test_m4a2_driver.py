"""Controlled observations of actual seed, refinement, folding and finish stages."""

from contextlib import ExitStack, contextmanager
from fractions import Fraction as F
import math
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
from xpsloss import DepthAverageConvergenceError
from xpsloss._single_loss_types import COLUMNS, COMPONENTS, PARTS, QUANTITIES, EPS
import m4a2_states as m

api, panels, fold = m.api, m.panels, m._fold


def columns(a=0.5, h=-0.25, b=0.0):
    return np.array([a, b, h, 0.0, 0.0, 0.0, 0.0, 0.0, a + b, h, a + b + h])


def kernel(t, side):
    return np.tile(columns()[:, None], (1, len(t))), np.zeros((11, len(t)))


def near(log, key, y, r):
    log.check(key, abs(F(y) - F(r)) <= 8 * F(EPS) * abs(F(r)))


def layer(row, log):
    i = row["inputs"]
    a, b = i.get("a", 2), 4
    tails = row["rule"] == "direct_tails"
    p = api.DeltaPrimary([10, 100, 200], [1, 2**60, 2]) if tails else api.DeltaPrimary([10], [1])
    edges = [90, 110] if tails else [6, 8, 9, 11]
    s2 = i.get("multiplier") == "inverse_window_area"
    r, mesh = panels.integrate_window(
        kernel,
        (a, b),
        features={},
        weight_length_A=8.0,
        primary=p,
        output_edges_eV=edges,
        epsrel=1e-6,
        area_epsabs=1e-10,
        spectrum_epsabs_fraction=1e-9,
        panel_limit=400,
        max_cells=2000000,
        **(m.S2 if s2 else dict(multiplier="external", multiplier_A=4.0)),
    )
    if tails:
        for side, zero, sat in (("low", 1, 2), ("high", 2, 4)):
            for name, v in (("zero_loss", zero), ("single_loss_total", sat)):
                near(log, f"outside.{side}.{name}", r.outside[side][name], v)
        return
    mult = 4 / (b - a) if s2 else 4.0
    near(log, "single_loss_weight", r.single_loss_weight, mult * (b - a) / 4)
    for c, k in zip(COLUMNS, columns()):
        near(log, "window_area." + c, r.window_area[c], (b - a) * k)
        near(log, "cutoff_sensitivity." + c, r.cutoff_sensitivity[c], mult * a * k)
        near(log, "diagnostics.LbK_b." + c, r.diagnostics["LbK_b"][c], mult * b * k)
    names = m.old.region_names(r)
    extended = [-math.inf, *edges, math.inf]
    for n, l, u in zip(names, extended, extended[1:]):
        width = max(0, min(b, 10 - l) - max(a, 10 - u))
        zero = float(l <= 10 < u)
        sat = [width * mult * k for k in (0.25, 0.5, -0.25)]
        for q, target in zip(QUANTITIES, [zero, *sat, zero + sat[0]]):
            near(log, n + "." + q, m.old.read_key(r, n + "." + q)[0], target)
            if n.startswith("bin"):
                near(
                    log,
                    n + "." + q + "_per_eV",
                    m.old.read_key(r, n + "." + q + "_per_eV")[0],
                    target / (u - l),
                )
        for s, here in (("plus", l < 10 - a <= u), ("minus", l <= 10 - a < u)):
            for c, k in zip(COMPONENTS, (0.25, 0.5, -0.25)):
                key = f"{n}.cutoff_sensitivity.{s}.{c}"
                near(log, key, m.old.read_key(r, key)[0], mult * a * k * here)
    fixed = fold.fold(mesh, p, edges, multiplier=r.multiplier_A, hole_charge=1.0)
    for q in QUANTITIES:
        rows = [fixed.outputs[n + "." + q] for n in names]
        bound = math.fsum(x["rounding"] + x["edge"] for x in rows)
        if q == "zero_loss":
            target = 1.0
        else:
            c = "total" if q == "spectrum" else q.removeprefix("single_loss_")
            target = (
                1 + r.single_loss_weight
                if q == "spectrum"
                else float(F(r.multiplier_A) * F(r.window_area[c]))
            )
            bound += float(
                F(r.multiplier_A)
                * F(r.diagnostics["error_budget"]["outputs"]["window_area." + c]["rounding"])
            )
        terms = [x["value"] for x in rows]
        log.equal(
            "diagnostics.balance." + q,
            r.diagnostics["balance"][q],
            dict(
                difference=math.fsum([*terms, -target]),
                rounding_bound=bound + 8 * EPS * (math.fsum(map(abs, terms)) + abs(target)),
            ),
        )


def inputs(row, log):
    i = row["inputs"]
    spec = row["entry"] == "api_spectrum"
    args, kw = m.fixture(i["fixture"], spec)
    if (
        i["fixture"] == "derived"
        or i.get("derived")
        or i.get("weight") == "derived"
        or (i["fixture"] == "Z1" and i.get("weight") != "external")
    ):
        kw.update(m.DERIVED)
    if i["fixture"] == "external_m":
        kw.update(multiplier="external", multiplier_A=7.0)
    if i["fixture"] == "s2":
        kw.update(m.S2)
    if "window" in i:
        kw["loss_window_eV"] = i["window"]
    if "c" in i:
        kw["extra_inverse_length_per_A"] = math.ulp(0.0) if i["c"] == "min_subnormal" else i["c"]
    for k in ("weight_length_kind", "weight_length_A", "multiplier"):
        if k in i:
            kw[k] = i[k]
    if i.get("external_m_A") is not None:
        kw["multiplier_A"] = i["external_m_A"]
    if "argument" in i:
        tokens = dict(
            none=None,
            bool=True,
            np_bool=np.bool_(True),
            array2=np.array([0, 1]),
            zero=0,
            text_number="2",
            unknown="bad",
            nan=math.nan,
            inf=math.inf,
            negative=-1,
            overflow_int=2**2000,
        )
        kw[i["argument"]] = tokens[i["value_token"]] if "value_token" in i else i["value"]
    if "primary" in i:
        args[2] = m.old.make_primary(i["primary"])
    with (
        patch.object(api, "_call_kernel", wraps=api._call_kernel) as spy,
        patch.object(api, "_seed_length", wraps=api._seed_length) as seed,
    ):
        try:
            r = m.call(args, kw)
        except Exception as exc:
            exception(row, log, exc)
            for k in ("panel_eV", "panel_width_ulp", "required_panels"):
                if k in row["expect"]:
                    log.check(k, np.array_equal(getattr(exc, k), row["expect"][k]))
            if type(exc) is ValueError:
                from test_m4a2_records import old_length_message

                log.equal("message", str(exc), old_length_message(row))
            if row["expect"].get("exception") == "SingleLossInputError" and i.get("argument"):
                log.check("argument", i["argument"] in str(exc))
            r = None
        else:
            log.check("success", row["expect"]["status"] == "success")
        if row["expect"].get("kernel_evaluations") == 0 or row["rule"] == "structural_zero":
            log.check("calls", spy.call_count == sum(len(v.args[2]) for v in spy.call_args_list) == 0)
        if "seed_entries" in row["expect"]:
            log.equal("seed_entries", seed.call_count, 0)
    if r is not None and row["rule"] == "structural_zero":
        center, interval = m.reciprocal(0.0, 0.0, 3)
        d = r.diagnostics["weight_length"]
        log.equal("weight_length_A", r.weight_length_A, center)
        log.equal("weight_length_interval_A", r.weight_length_interval_A, interval)
        log.check(
            "seed_skipped",
            d["seed_skipped"] == "structural_zero"
            and d["consistency"] == "not_applicable"
            and d["area13"] == d["error13"] == d["kernel_calls"]["calls"] == seed.call_count == 0,
        )
        log.check(
            "J",
            r.single_loss_weight == 0
            and np.array_equal(r.spectrum, r.zero_loss)
            and all(side["spectrum"] == side["zero_loss"] for side in r.outside.values()),
        )
    if r is not None and row["rule"] == "zero_primary":
        log.check("area11", len(r.window_area) == 11 and r.window_area["total"] > 0)
        log.check(
            "w",
            abs(r.single_loss_weight - 1) <= EPS
            and r.single_loss_weight_error == abs(r.single_loss_weight - 1),
        )
        for n in m.old.region_names(r):
            for q in QUANTITIES:
                log.check(n + "." + q, m.old.read_key(r, n + "." + q) == (0.0, 0.0))
                if q != "spectrum":
                    record = r.diagnostics["error_budget"]["outputs"][n + "." + q]
                    log.check(n + "." + q + ".budget", record["budget"] == 0 and record["ratio"] is None)
        log.check(
            "balance",
            all(
                v["difference"] == 0
                for v in r.diagnostics["balance"].values()
                if v is not r.diagnostics["balance"]["window_area.components"]
            ),
        )


def exception(row, log, exc):
    e = row["expect"]
    for k in ("exception", "reason", "phase", "stage"):
        target = e.get(k)
        value = type(exc).__name__ if k == "exception" else getattr(exc, k, None)
        log.check("exception." + k, value == target, actual=value, expected=target)
    if type(exc) in (api.SingleLossConvergenceError, api.SingleLossInputError):
        log.check("exception.attributes", hasattr(exc, "phase") and hasattr(exc, "stage"))


class Simulation:
    """Observe phase boundaries and change only the registered kernel/aggregate."""

    def __init__(self, row, local_first=False):
        self.row, self.i = row, row["inputs"]
        self.name = self.i["fixture"]
        self.target = self.i.get("phase", "production")
        self.phase, self.stage = None, None
        self.local_first = local_first
        self.calls, self.folds, self.meshes = [], [], []
        self.inner = None
        self.original_areas = panels._areas
        self.original_fold = fold._fold

    @contextmanager
    def patches(self):
        integrate = panels._integrate
        sample = api._DepthSampler.__call__

        def drive(*a, **kw):
            self.phase = "weight_length_seed" if kw.get("seed") else "production"
            return integrate(*a, **kw)

        def sampling(sampler, points, stage):
            self.stage = stage
            return sample(sampler, points, stage)

        with ExitStack() as stack:
            for owner, name, f in (
                (api, "_call_kernel", self.kernel),
                (api._DepthSampler, "__call__", sampling),
                (panels, "_integrate", drive),
                (panels, "_areas", self.areas),
                (fold, "_fold", self.fold),
            ):
                stack.enter_context(patch.object(owner, name, f))
            original_mesh = panels._initial_mesh
            if self.name in ("ties", "two_stages"):
                edges = [2.0, 3.0, 5.0, 7.0] if self.name == "ties" else [2.0, 3.0, 4.0]

                def mesh(*a, **kw):
                    return original_mesh(*a, **(kw | dict(initial_edges=edges)))

                stack.enter_context(patch.object(panels, "_initial_mesh", mesh))
            yield

    def kernel(self, material, energy, grid, **kw):
        self.calls.append((self.phase, self.stage, kw.copy(), grid.copy(), energy, material))
        if self.name == "upstream" and self.phase == self.target and self.stage == self.i["fail_stage"]:
            self.inner = DepthAverageConvergenceError(
                "budget_not_met",
                grid[1],
                1,
                column="total",
                estimate=1.25e-3,
                budget=5e-4,
                subdivisions=3,
                panels=4,
            )
            raise self.inner
        v = columns(1, 0) if self.name == "two_stages" else columns()
        if self.name == "length_kernel":
            v = columns(0.125, -1 / 32, kw["weight_length_A"] / 64)
        values = np.tile(v[:, None], (1, len(grid)))
        errors = np.zeros_like(values)
        if self.name == "length_kernel" and self.phase == "weight_length_seed":
            errors[1] = 2.0
        return SimpleNamespace(
            **dict(zip(COLUMNS, values)), component_errors=dict(zip(COLUMNS, errors)), warnings=[]
        )

    def fold(self, ps, geom, mval, lower):
        folded = self.original_fold(ps, geom, mval, lower)
        self.folds.append((self.phase, self.stage, mval, ps, geom, lower, folded))
        return folded

    def areas(self, ps):
        out, local = self.original_areas(ps)
        self.meshes.append((self.phase, self.stage, [(p.left_eV, p.right_eV) for p in ps]))
        seed = self.phase == "weight_length_seed"

        def put(column, y, e=0.0, part="verification"):
            k = "window_area." + column
            parts = [e if x == part else 0.0 for x in PARTS]
            out[k] = fold._entry(float(y), parts)
            local[k] = np.array([e] + [0.0] * (len(ps) - 1))

        if self.name == "consistency":
            x = 0.25 + 2**-23
            b = {"inside": math.nextafter(x, 0.0), "equal": x, "outside": math.nextafter(x, math.inf)}[
                self.i["position"]
            ]
            put("term_13a", 0.25 if seed else b, 2**-24)
            if not seed and not self.stage and self.i["position"] == "outside":
                put("term_13b", 0.0, 2e-10)
        if self.name == "two_stages" and self.stage == 0:
            put("total", 3.0)
            put("term_13b", 0.0, 2e-10)
        if self.name == "s2_equals_weight":
            put("total", 1 / 20)
        if self.phase != self.target:
            return out, local
        held = (2**-12, 2**-13)
        values = {
            "negative": (-1, 2**-24),
            "positive-limit": held,
            "zero-limit": (0, 0),
            "inverse-overflow": (math.ulp(0.0), 0),
            "end-overflow": (2**-1023, 2**-1024),
            "release": held,
            "negative-hold": (-held[0], held[1]),
            "ties": held if self.local_first else (0, 0),
            "refined_width": held,
            "refined_cells": held,
            "small_area_plus_c": (2**-40, 2**-39),
            "exact_sum": (1.0, 0.0),
            "budget_point": (1, 0),
            "budget_loss": (1, 0),
        }
        if self.name == "upstream" and self.i["fail_stage"] == 1:
            values[self.name] = held
        if self.name in values:
            y, e = (1, 0) if self.stage and self.name != "zero-limit" else values[self.name]
            put("term_13a" if seed else "total", y, e)
        if self.name.startswith("budget_"):
            put(
                "term_13b",
                0.0,
                2e-10,
                "kernel_point_error" if self.name == "budget_point" else "verification",
            )
        return out, local

    def run(self):
        name = self.name
        args, kw = m.fixture()
        kw["loss_window_eV"] = (2.0, 7.0) if name == "ties" else (2.0, 4.0)
        if name == "refined_width":
            kw["loss_window_eV"] = (2.0, 2.0 + 384 * math.ulp(2.0))
        seed = self.target == "weight_length_seed"
        if seed or name in ("upstream", "consistency", "length_kernel", "constant"):
            kw.update(m.DERIVED)
        if not seed or name == "upstream":
            kw.update(m.S2)
        if name in ("consistency", "length_kernel"):
            kw["multiplier"] = "weight_length" if name == "consistency" else "external"
        if name == "length_kernel":
            kw["multiplier_A"] = 7.0
        if name in ("small_area_plus_c", "exact_sum", "constant"):
            kw["extra_inverse_length_per_A"] = {
                "small_area_plus_c": 1.0,
                "exact_sum": 3 * 2**-54,
                "constant": 0.125,
            }[name]
        if name in (
            "negative",
            "negative-hold",
            "release",
            "ties",
            "positive-limit",
            "zero-limit",
            "upstream",
            "refined_width",
            "refined_cells",
            "inverse-overflow",
            "end-overflow",
        ):
            kw["area_epsabs"] = 1.0
        if name == "end-overflow":
            kw["epsrel"] = 0.5
        if name in ("positive-limit", "budget_point", "budget_loss"):
            kw["panel_limit"] = 1
        elif name in (
            "zero-limit",
            "negative-hold",
            "release",
            "consistency",
            "upstream",
            "refined_width",
            "refined_cells",
        ):
            kw["panel_limit"] = 2
        elif name == "ties":
            kw["panel_limit"] = 4
        if not seed or name == "upstream":
            edges = (
                [0, 20]
                if name == "refined_cells"
                else [6, 7, 8, 9, 11]
                if name == "two_stages"
                else [6, 8, 9, 11]
            )
            args += [api.DeltaPrimary([10], [1]), edges]
        if name == "refined_cells":
            kw["max_cells"] = 3
        with self.patches():
            if seed and name != "upstream":
                return m.seed_observation(args, kw)
            return m.call(args, kw)


def upstream_checks(log, sim, exc):
    log.equal("inner_attributes", exc.inner_attributes, vars(sim.inner))
    log.check(
        "inner_exception",
        exc.inner_exception == "DepthAverageConvergenceError" and exc.inner_reason == "budget_not_met",
    )
    log.check(
        "index",
        exc.index == sim.inner.index
        and exc.loss_eV == sim.inner.loss_eV
        and exc.column == sim.inner.column
        and exc.node_kind == "interior",
    )
    if sim.target == "weight_length_seed":
        log.check("beta_calls", not any(x[0] == "production" for x in sim.calls))
    if exc.stage == 0:
        log.equal("first_budget", exc.diagnostics["error_budget"], dict(outputs={}, max_ratio=None))


def simulation(row, log, local_first=False):
    sim = Simulation(row, local_first)
    name = sim.name
    try:
        result = sim.run()
    except (api.SingleLossInputError, api.SingleLossConvergenceError) as exc:
        exception(row, log, exc)
        if name == "ties":
            log.check(
                "split_panel",
                len(sim.meshes) == 2
                and set(sim.meshes[0][2]) - set(sim.meshes[1][2])
                == {(2.0, 3.0) if local_first else (3.0, 5.0)},
            )
        d = exc.diagnostics
        if name == "consistency":
            d = d["weight_length"]["consistency"]
            log.check("consistency.delta", d["delta"] == math.nextafter(0.25 + 2**-23, math.inf) - 0.25)
            log.equal("consistency.allowance", d["allowance"], 2**-23)
        elif name == "upstream":
            upstream_checks(log, sim, exc)
        else:
            log.check("kernel_calls", len(sim.calls) > 0)
            if "denominator" in d:
                den = d["denominator"]
                log.check("denominator", den["stage"] == (exc.stage if exc.stage is not None else 0))
            if row["rule"] == "denominator_failure":
                log.check("fold_calls", not sim.folds)
                injected = {
                    "negative": (-1.0, 2**-24),
                    "positive-limit": (2**-12, 2**-13),
                    "zero-limit": (0.0, 0.0),
                    "inverse-overflow": (math.ulp(0.0), 0.0),
                    "end-overflow": (2**-1023, 2**-1024),
                }
                value, error = injected[name]
                limit = F(d["effective_settings"]["epsrel"]) * abs(F(value))
                log.equal("denominator.value", den["value"], value)
                log.equal("denominator.error", den["error"], error)
                log.equal("denominator.threshold", den["threshold"], float(limit))
                log.check(
                    "denominator.resolved", den["resolved"] == (F(value) > F(error) and F(error) <= limit)
                )
                if exc.reason.endswith("unresolved"):
                    worst = d["error_budget"]["max_ratio"]
                    column = "term_13a" if sim.target == "weight_length_seed" else "total"
                    entry = d["error_budget"]["outputs"]["window_area." + column]
                    log.check("raw_budget", entry["estimate"] == error <= entry["budget"])
                    log.check(
                        "worst_key",
                        worst is None
                        or worst["value"] == d["error_budget"]["outputs"][worst["key"]]["ratio"],
                    )
            if row["rule"] == "old_failure":
                log.equal("diagnostics.stages", d["stages"], 0)
                if name in ("budget_point", "budget_loss"):
                    log.equal("worst_key", d["error_budget"]["max_ratio"]["key"], "window_area.term_13b")
                    entry = d["error_budget"]["outputs"]["window_area.term_13b"]
                    dominant = "kernel_point_error" if name == "budget_point" else "verification"
                    log.equal("worst_part", entry[dominant], 2e-10)
                elif name == "refined_width":
                    log.equal("panel_width_ulp", exc.panel_width_ulp, 192.0)
                    log.equal("panel_eV", exc.panel_eV, (2.0, 2.0 + 192 * math.ulp(2.0)))
                else:
                    log.check("cells", "max_cells=3" in str(exc) and "cells_at_least=6" in str(exc))
        return
    log.check("exception.reason", row["expect"]["status"] == "success")
    if row["rule"] == "denominator_progress":
        seed = sim.target == "weight_length_seed"
        d = result["diag"] if seed else result.diagnostics
        denom = d["denominator"]
        if name in ("small_area_plus_c", "exact_sum"):
            a, e, c = (2**-40, 2**-39, 1.0) if name == "small_area_plus_c" else (1.0, 0.0, 3 * 2**-54)
            length, interval = m.reciprocal(a, e, c)
            log.equal("weight_length_A", result["length_A"], length)
            log.equal("weight_length_interval_A", result["interval_A"], interval)
            raw = d["error_budget"]["outputs"]["window_area.term_13a"]
            log.check("raw_budget", raw["estimate"] == e and raw["budget"] == max(1e-10, 1e-6 * a))
            log.check(
                "denominator",
                denom["resolved"] and denom["stage"] == 0 and denom["value"] == float(F(a) + F(c)),
            )
        else:
            meshes = [ps for phase, stage, ps in sim.meshes if phase == sim.target]
            left = (2.0, 3.0) if local_first else (3.0, 5.0) if name == "ties" else (2.0, 4.0)
            log.check("split_panel", len(meshes) == 2 and set(meshes[0]) - set(meshes[1]) == {left})
            log.check("final_panels", d["final_panels"] == (4 if name == "ties" else 2) and d["stages"] == 1)
            log.check("denominator", denom["resolved"] and denom["stage"] == 1 and denom["value"] == 1.0)
            log.check(
                "fold_calls", [(p, s) for p, s, *_ in sim.folds] == ([] if seed else [("production", 1)])
            )
            if name == "ties" and not local_first:
                other = m.Ledger(row)
                simulation(row, other, True)
                log.check("local_error_first", all(x["status"] == "pass" for x in other.items.values()))
    elif name == "length_kernel":
        d = result.diagnostics
        for phase, length, count in (
            ("weight_length_seed", 20.0, d["weight_length"]["kernel_calls"]["calls"]),
            ("production", result.weight_length_A, d["kernel_calls"]["calls"]),
        ):
            calls = [x for x in sim.calls if x[0] == phase]
            label = "alpha" if phase == "weight_length_seed" else "beta"
            log.check(
                label + ".call_kwargs.weight_length_A",
                len(calls) == count == 1 and all(x[2]["weight_length_A"] == length for x in calls),
            )
            log.check(
                label + ".call_kwargs",
                all(
                    x[2]["epsilon_momentum"] == "k_total"
                    and x[2]["epsrel"] == 1e-9
                    and x[2]["epsabs"] == 1e-12
                    and x[2]["limit"] == 200
                    and x[4] == 1000
                    and tuple(x[3][[0, -1]]) == (2, 4)
                    for x in calls
                ),
            )
        near(log, "seed.area13", d["weight_length"]["area13"], 0.25)
        near(log, "window_area.term_13b", result.window_area["term_13b"], 2 * result.weight_length_A / 64)
        near(log, "LbK_b", d["LbK_b"]["term_13b"], 7 * 4 * result.weight_length_A / 64)
    elif name == "consistency":
        d = result.diagnostics["weight_length"]["consistency"]
        x = 0.25 + 2**-23
        log.equal(
            "consistency.delta",
            d["delta"],
            (math.nextafter(x, 0) if sim.i["position"] == "inside" else x) - 0.25,
        )
        log.equal("consistency.allowance", d["allowance"], 2**-23)
    elif name == "two_stages":
        log.check("stage_m", [f[2] for f in sim.folds] == [float(F(1, 3)), 1 / result.window_area["total"]])
        _, _, _, ps, geom, lower, last = sim.folds[-1]
        reference = sim.original_fold(ps, geom, result.multiplier_A, lower)
        near(log, "retained_panel.folded_values", result.single_loss_total[0], 0.5)
        log.equal("retained_panel.parts5", last.outputs, reference.outputs)
        log.equal("local", last.panel_errors, reference.panel_errors)
        log.equal("support", last.support_widths, reference.support_widths)
        for key, rec in result.diagnostics["error_budget"]["outputs"].items():
            if key in reference.outputs:
                log.equal(
                    key + ".budget", rec["budget"], max(1e-9, 1e-6 * abs(reference.outputs[key]["value"]))
                )
        log.check(
            "balance",
            all(abs(x["difference"]) <= x["rounding_bound"] for x in result.diagnostics["balance"].values()),
        )
    elif name == "constant":
        near(log, "weight_length_A", result.weight_length_A, 8 / 9)
        near(log, "multiplier_A", result.multiplier_A, 2.0)
        near(log, "window_area.hole", result.window_area["hole"], -0.5)
        near(log, "w", result.single_loss_weight, 1.0)
        for n, v in zip(m.old.region_names(result), (0.0, -1.0, 0.0, 0.0, 0.0)):
            near(log, n + ".single_loss_hole", m.old.read_key(result, n + ".single_loss_hole")[0], v)


def finish_error(row, log):
    # Fixed-m outputs enter the same finish step used after the final area/region budgets.
    args, kw = m.fixture(spectrum=True)
    args[3] = [6.0, 8.0, 10.0]
    r = m.call(args, kw)
    r.multiplier_A, r.multiplier_convention = 0.5, "inverse_window_area"
    r.multiplier_interval_A = m.reciprocal(2.0, 1.0)[1]
    r.window_area["total"], r.window_area_error["total"] = 2.0, 1.0
    r.diagnostics["multiplier"] = {}
    r.single_loss_weight, r.single_loss_weight_error = 1.0, 0.5 + EPS
    for name, value in zip(QUANTITIES, (3.0, 2.0, 0.0, -2.0, 5.0)):
        err = 1 / 16 if name == "zero_loss" else 1 / 8
        getattr(r, name)[:] = value
        getattr(r, name + "_error")[:] = err
        for side in r.outside.values():
            side[name], side[name + "_error"] = value, err
    r.single_loss_no_hole[1], r.single_loss_no_hole_error[1] = 1.0, 2**-54
    r.zero_loss_per_eV_error[:] = float((F(1, 16) + F(EPS) * 3) / 2)
    for s in ("plus", "minus"):
        for c in COMPONENTS:
            r.cutoff_sensitivity_bins[s][c][:] = -2.0
            r.cutoff_sensitivity_bins_error[s][c][:] = 0.125
            for side in r.outside.values():
                side["cutoff_sensitivity"][s][c], side["cutoff_sensitivity_error"][s][c] = -2.0, 0.125
    from copy import deepcopy

    fixed = deepcopy(r)
    fold.propagate_multiplier(r)
    m.check_errors(log, r, fixed, check_balance=False)
    log.equal("signed_error", float(r.single_loss_hole_error[0]), 2.25)
    log.equal("zero_loss_error", float(r.zero_loss_error[0]), 1 / 16)
    log.equal("zero_value_error", float(r.single_loss_no_hole_error[0]), 0.25)
    log.equal("upward_rounding", float(r.single_loss_no_hole_error[1]), math.nextafter(1.0, math.inf))


def driver(row, log):
    if row["rule"] in ("input_before_kernel", "structural_zero", "zero_primary"):
        inputs(row, log)
    elif row["rule"] in ("analytic_map", "direct_tails", "partial_sensitivity"):
        layer(row, log)
    elif row["rule"] == "error_formula":
        finish_error(row, log)
    else:
        simulation(row, log)
