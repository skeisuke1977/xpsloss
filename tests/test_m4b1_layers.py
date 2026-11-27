"""Analytical checks of the numerical layers, independent of the comparison tables."""

import math
from dataclasses import replace
from decimal import Decimal, localcontext
from fractions import Fraction as Q
import unittest
from unittest.mock import patch

import numpy as np
from numpy.polynomial.legendre import Legendre

from xpsloss import (
    DeltaPrimary,
    MULTIPLE_LOSS_MODEL_ID as MODEL,
    MultipleLossConvergenceError,
    MultipleLossInputError,
)
from xpsloss import _loss_orders as orders, _loss_panels as lp, _loss_reference_models as models
from xpsloss._order_bounds import log_laplace, log_power
from xpsloss._order_quadrature import area_terms, mandatory_knots, tau_grid


def box(t, side):
    v = np.zeros((11, len(t)))
    v[[0, 8, 10]] = 0.1
    return v, np.zeros_like(v)


def calculate(**kwargs):
    args = dict(
        features={},
        weight_length_A=2.0,
        primary=DeltaPrimary([1000], [1]),
        output_edges_eV=[988, 990, 994, 995, 996, 998, 1000, 1001],
        max_order=3,
        model=MODEL,
    )
    return orders.integrate_orders(kwargs.pop("kernel", box), kwargs.pop("window", (2.0, 5.0)), **(args | kwargs))


class Layers(unittest.TestCase):
    def test_exact_knots_and_geometric_grid(self):
        a, b = 1.3, 7.7
        self.assertEqual(
            mandatory_knots({}, (a, b), 3, 10), sorted(float(k * Q(a) + (3 - k) * Q(b)) for k in range(4))
        )
        tau = tau_grid(1.0)
        self.assertEqual(len(tau), 93)
        np.testing.assert_allclose(tau[1:] / tau[:-1], 10**0.05)

    def test_degree_32_gauss_product(self):
        t = lp._nodes(2.0, 5.0)[0]
        x = (2 * t - 7) / 3
        for degree in (9, 16):
            v, _ = box(t, None)
            v[[0, 8, 10]] = 0.05 * (1 + Legendre.basis(degree)(x))
            p = lp.Panel.from_lobatto(2.0, 5.0, v, 0 * v)
            panels = orders.convolve_panels([p], (2.0, 5.0), features={}, model=MODEL, max_order=2)
            expected = 0.05**2 * 3 * (1 + (-1) ** degree / (2 * degree + 1))
            at_knot = next(p for p in panels[2] if p.right_eV == 7.0)
            self.assertAlmostEqual(at_knot.values[0, -1] / expected, 1.0, places=13)

    def test_tau_grid_upper_edge(self):
        for a, count in ((1.0, 93), (4.0, 81), (0.4, 100)):
            with self.subTest(a=a):
                grid = tau_grid(a)
                self.assertEqual(len(grid), count)
                self.assertLessEqual(Q(grid[-1]), 40 / Q(a))
                self.assertGreater(Q(10 ** (-3 + count / 20)), 40 / Q(a))
        self.assertEqual(tau_grid(4.0)[-1], 10.0)

    def test_touching_and_negative_interior_minima(self):
        for center, depth in ((2.0, 0.0), (2.04, 1e-5)):
            def quadratic(t, side):
                v, e = box(t, side)
                v[[0, 8, 10]] = 0.1 * (t - center) ** 2 - depth
                return v, e

            r, _ = calculate(kernel=quadratic, window=(1., 3.), initial_edges=[1., 3.], max_order=1)
            codes = {w["code"] for w in r.diagnostics["warnings"]}
            for key, code in (("term_13a", "negative_13a_kernel"), ("total", "negative_total_kernel")):
                minimum = r.diagnostics["negative_values"][key]
                self.assertGreaterEqual(minimum["sampled_min"], 0)
                self.assertAlmostEqual(minimum["continuous_min"], -depth, places=14)
                self.assertEqual(code in codes, depth > 0)
            tail = r.truncation["neglected_total_area"]
            self.assertEqual(tail["reason"], "repeated_kernel_not_nonnegative" if depth else None)
            self.assertEqual(tail["value"] is None, depth > 0)

    def test_order_budgets_below_binary64_floor(self):
        tiny = math.ulp(0.0)
        for m in (1e-200, tiny):
            r, _ = calculate(multiplier="external", multiplier_A=m, order_epsabs_fraction=tiny,
                             primary=DeltaPrimary([1000], [0.01]))
            for order in r.orders[1:]:
                self.assertEqual(order.area["total"], 0.)
                self.assertEqual(order.area_error_parts["inherited"]["total"], tiny)
                self.assertEqual(order.area_error_parts["local"]["total"], tiny)
                self.assertTrue(np.any(order.total_error_parts["local"] == tiny))
                self.assertEqual(order.mesh["stages"], 0)
        # Both sides round to zero, but a real excess must still nominate the panel.
        budget = Q(tiny) / 8
        split = np.array([False])
        orders._choose_split(np.array([2 * budget]), 2 * budget, Q(0), budget, 1., split)
        self.assertTrue(split[0])

        tolerances = dict(epsrel=np.float32(1e-6), area_epsabs=np.float32(1e-10),
                          spectrum_epsabs_fraction=np.float32(1e-9), order_epsrel=np.float32(1e-6),
                          order_epsabs_fraction=np.float32(1e-9))
        scalars, _ = calculate(max_order=2, **tolerances)
        binary64, _ = calculate(max_order=2, **{k: float(v) for k, v in tolerances.items()})
        np.testing.assert_array_equal(scalars.orders[1].total, binary64.orders[1].total)
        self.assertEqual(scalars.orders[1].area_error_parts, binary64.orders[1].area_error_parts)

    def test_outward_first_total_with_multiplier_interval(self):
        r, panels = calculate(max_order=1, primary=None, output_edges_eV=None)
        # The exact product is halfway between two floats; nearest rounds down.
        e = 0.25 + 3 * math.ulp(0.25)
        s = replace(r.single_loss, multiplier_A=1.5, multiplier_interval_A=(1.25, 1.75),
                    window_area=r.single_loss.window_area | {"total": 0.0},
                    window_area_error=r.single_loss.window_area_error | {"total": e},
                    single_loss_weight=-0.0, single_loss_weight_error=0.5)
        o = orders.first_order(s, panels[1], None, {}, orders.seed_errors(s), False)
        exact = Q(3, 2) * Q(e)
        base = o.area_error_parts["inherited"]["total"]
        self.assertLess(Q(math.nextafter(base, -math.inf)), exact)
        self.assertGreaterEqual(Q(base), exact)
        self.assertEqual(o.area_error_parts["multiplier_interval"]["total"], 0.5 - base)
        self.assertEqual(o.area["total"].hex(), s.single_loss_weight.hex())
        self.assertEqual(o.area_error["total"].hex(), s.single_loss_weight_error.hex())

    def test_nan_in_order_arithmetic_has_stage_and_reason(self):
        values, errors = box(lp._nodes(2., 5.)[0], None)
        panel = lp.Panel.from_lobatto(2., 5., values, errors)
        with patch.object(orders, "area_terms", return_value=np.full(11, np.nan)):
            for entry in (lambda: calculate(max_order=2),
                          lambda: orders.convolve_panels([panel], (2., 5.), features={}, model=MODEL, max_order=2)):
                with self.assertRaises(MultipleLossConvergenceError) as why:
                    entry()
                self.assertEqual((why.exception.reason, why.exception.order, why.exception.stage),
                                 ("nonfinite_output", 2, "order_refinement"))

        cell = orders._fold._cell_weights
        order_cell = orders._order_cell

        def weights_with_nan(*args):
            weights, *rest = cell(*args)
            return np.full_like(weights, np.nan), *rest

        def invalid(*args):
            with patch.object(orders._fold, "_cell_weights", weights_with_nan):
                return order_cell(*args)

        with patch.object(orders, "_order_cell", invalid):
            with self.assertRaises(MultipleLossConvergenceError) as why:
                calculate(max_order=2)
        self.assertEqual((why.exception.reason, why.exception.order, why.exception.stage),
                         ("nonfinite_output", 2, "order_fold"))
        self.assertEqual(why.exception.inner["type"], "ValueError")

        with self.assertRaises(MultipleLossInputError) as why:
            calculate(order_epsrel=math.nan)
        self.assertEqual((why.exception.reason, why.exception.argument), ("invalid_settings", "order_epsrel"))

    def test_box_areas_sensitivity_and_direct_bins(self):
        res, panels = calculate()
        for n, o in enumerate(res.orders, 1):
            self.assertAlmostEqual(o.area["total"], 0.6**n)
            self.assertAlmostEqual(o.cutoff_sensitivity["total"], n * 0.4 * 0.6 ** (n - 1))
            self.assertAlmostEqual(sum(o.total) + sum(side["total"] for side in o.outside.values()), 0.6**n)
            self.assertEqual(len(panels[n]), o.mesh["panels"])
        self.assertEqual(res.orders[1].total[-1], 0.0)
        self.assertEqual(res.orders[1].total_error[-1], 0.0)
        self.assertEqual(res.spectrum[-1], 1.0)

    def test_inherited_point_path_and_local_propagation(self):
        captured = []
        fold = orders._fold_order

        def inspect(panels, geom, **kw):
            for p, parts in zip(panels, kw["error_parts"]):
                if p.right_eV == 7:
                    captured.append(parts["inherited"][-1, -1])
            return fold(panels, geom, **kw)

        def inject(role, t, v, e):
            if role == "interpolation":
                e[[0, 8, 10]] += 0.001
            return v, e

        with patch.object(orders, "_fold_order", inspect):
            calculate(
                max_order=2,
                inject=inject,
                epsrel=0.5,
                area_epsabs=1.0,
                spectrum_epsabs_fraction=1.0,
                initial_edges=[2, 5],
            )
        self.assertAlmostEqual(captured[-1], 0.000603, places=12)

        def local(n, role, bounds, t, v, e):
            if role == "interpolation" and n == 2:
                e[[0, 8, 10]] += 1e-3 if bounds[1] - bounds[0] >= 3 else 1e-12
            return v, e

        res, _ = calculate(order_inject=local)
        self.assertGreater(res.orders[1].mesh["stages"], 0)
        self.assertEqual(res.orders[1].total_error_parts["propagated"][1], 0.0)
        self.assertGreater(res.orders[2].total_error_parts["propagated"][1], 0.0)

    def test_absolute_norm_and_s2_seed(self):
        def linear(t, s):
            v, e = box(t, s)
            v[[2, 9]] = 0.02 * (t - 3)
            v[10] += v[9]
            return v, e

        def inject(role, t, v, e):
            if role == "interpolation":
                e[[0, 8, 10]] += 0.001
            return v, e

        kw = dict(
            kernel=linear,
            primary=None,
            output_edges_eV=None,
            inject=inject,
            initial_edges=[2, 4, 5],
            epsrel=0.5,
            area_epsabs=1.0,
            spectrum_epsabs_fraction=1.0,
        )
        result, _ = calculate(**kw)
        value = result.orders[1].area_error_parts["inherited"]["hole"]
        self.assertGreaterEqual(value, 0.0006)
        self.assertLess(value, 0.0006 * (1 + 1e-9))
        self.assertEqual(result.orders[0].mesh["panel_edges_eV"], [2, 4, 5])
        for run in (result, calculate(**kw, multiplier="external", multiplier_A=1.3)[0]):
            for order in run.orders[1:]:
                power = Q(run.single_loss.multiplier_A) ** order.order
                for column, terms in zip(lp.COLUMNS, zip(*order.mesh["local_area_terms"])):
                    with self.subTest(order=order.order, column=column, power=power):
                        exact = power * sum(map(Q, terms))
                        reported = order.area_error_parts["local"][column]
                        self.assertGreaterEqual(Q(reported), exact)
                        self.assertLess(Q(math.nextafter(reported, -math.inf)), exact)
        result, _ = calculate(**kw, multiplier="inverse_window_area")
        one = result.orders[0]
        self.assertEqual(one.area_error["total"], result.single_loss.single_loss_weight_error)
        self.assertEqual(one.area_error_parts["convention"]["total"], one.area_error["total"])
        self.assertEqual(one.area_error_parts["inherited"]["total"], 0.0)
        self.assertGreater(result.diagnostics["first_order"]["fixed_multiplier_area_error"]["total"], 0.001)

    def test_region_error_parts_survive_subnormal_scaling(self):
        tiny = math.ulp(0.0)
        zeros = np.zeros((11, 17))
        errors = np.full_like(zeros, tiny)
        panel = lp.Panel.from_lobatto(2.0, 5.0, zeros, errors)
        geometry = orders._fold.Geometry(DeltaPrimary([1000], [1]), [996, 997])
        parts = dict.fromkeys(("inherited", "propagated", "local"), errors)
        rows = orders._fold_order([panel], geometry, multiplier_power=Q(1, 16), error_parts=[parts])
        # Each of the three regions receives one eV of constant nodal uncertainty.
        for key, row in rows.items():
            for part in parts:
                with self.subTest(region_component=key, part=part):
                    self.assertEqual(row[part], tiny)

    def test_higher_order_errors_cover_the_reported_parts(self):
        result, _ = calculate(multiplier="inverse_window_area", output_edges_eV=[991.5, 993, 995])

        def check(total, parts):
            exact = sum(map(Q, parts))
            self.assertGreaterEqual(Q(total), exact)
            self.assertLess(Q(math.nextafter(total, -math.inf)), exact)

        for order in result.orders[1:]:
            for c, error in order.area_error.items():
                check(error, (p[c] for p in order.area_error_parts.values()))
            for c in lp.COMPONENTS:
                errors = getattr(order, c + "_error")
                parts = getattr(order, c + "_error_parts")
                for i, error in enumerate(errors):
                    check(error, (p[i] for p in parts.values()))
                for side in order.outside.values():
                    check(side[c + "_error"], side[c + "_error_parts"].values())

    def test_constant_panel_laplace_rounds_outward(self):
        taus = [0.0, 0.05, 12.0]
        with localcontext() as ctx:
            ctx.prec = 80
            for height, m in ((1e-240, 1e240), (1e240, 1e-240)):
                panel = lp.Panel.from_lobatto(2.0, 5.0, np.full((11, 17), height), np.zeros((11, 17)))
                e = Decimal.from_float(area_terms(panel)[0])
                logs = log_laplace([panel], 0, np.array(taus), m)
                for tau, reported in zip(taus, logs):
                    t = Decimal.from_float(tau)
                    weight = (1 - (-3 * t).exp()) / t if t else Decimal(3)
                    exact = Decimal.from_float(m) * (-2 * t).exp() * (Decimal.from_float(height) * weight + e)
                    self.assertGreaterEqual(Decimal.from_float(reported), exact.ln())

    def test_positive_all_axis_bound_underflow(self):
        def small(t, side):
            values, errors = box(t, side)
            return values * 1e-150, errors

        for model in (models.GEOMETRIC, models.FACTORIAL, models.EXCESS):
            result = models.integrate_orders_reference(
                small, (2.0, 5.0), features={}, weight_length_A=2.0, max_order=2, model=model,
                primary=DeltaPrimary([1000], [1]), output_edges_eV=[990, 1001],
            )
            self.assertEqual(result.truncation["neglected_total_area_bound"], math.ulp(0.0))

    def test_large_order_log_bounds(self):
        n = 2**53 + 3
        logs = np.array([-1.3, 1.3])
        for v, bound in zip(logs, log_power(logs, n)):
            self.assertGreaterEqual(Q(bound), n * Q(v))
        n = 300
        factorial = models._log_tail(models.FACTORIAL, [], 1.0)
        bound = factorial(None, n, np.array([-1.0]), None)[0]
        with localcontext() as ctx:
            ctx.prec = 80
            exact = -(n + 1) - Decimal(math.factorial(n + 1)).ln() + Decimal(-1).exp()
            self.assertGreaterEqual(Decimal.from_float(bound), exact)
            self.assertLess(Decimal.from_float(bound) - exact, Decimal("1e-9"))

    def test_resources_and_failure_timing(self):
        def forbidden(*args):
            self.fail("preflight called the kernel")

        with self.assertRaises(MultipleLossInputError) as cm:
            calculate(kernel=forbidden, max_order=4, max_convolution_points=4500)
        self.assertEqual((cm.exception.resource, cm.exception.lower_bound), ("max_convolution_points", 4743))
        with self.assertRaises(MultipleLossInputError) as cm:
            calculate(kernel=forbidden, max_order=40, max_cells=300)
        self.assertEqual((cm.exception.order, cm.exception.lower_bound), (34, 306))
        with self.assertRaises(MultipleLossConvergenceError) as cm:
            calculate(max_order=2, max_convolution_points=1054)
        self.assertEqual(cm.exception.stage, "after_first_order")
        self.assertEqual(
            cm.exception.estimate, 280.5 * 2 * (cm.exception.diagnostics["first_order"]["final_panels"] + 1)
        )
        jumps = np.arange(2.0, 5.5, 0.5)
        with self.assertRaises(MultipleLossConvergenceError) as cm:
            calculate(features={"jumps_eV": jumps}, max_order=2, max_convolution_points=4500)
        self.assertEqual(cm.exception.reason, "convolution_work_limit_reached")
        self.assertGreater(cm.exception.measured, 4500)

    def test_equality_support_and_cumulative_bounds(self):
        for edge, state in ((994.0, "truncated"), (math.nextafter(994.0, math.inf), "complete_on_display")):
            r, _ = calculate(max_order=2, output_edges_eV=[edge, 998, 1001])
            self.assertEqual(r.truncation["state"], state)
            self.assertEqual(r.truncation["bound_reason"][-1], "complete")
        r, _ = calculate(max_order=2, output_edges_eV=[988, 999, 1001])
        tau = np.r_[0.0, tau_grid(2.0)]
        lap = np.array([0.6 if t == 0 else 0.2 * math.exp(-2 * t) * (-math.expm1(-3 * t)) / t for t in tau])
        expected = np.min(np.exp(12 * tau) * lap**3 / (1 - lap))
        self.assertAlmostEqual(r.truncation["neglected_at_or_above_region_bound"][0] / expected, 1.0, places=12)

    def test_signed_inherited_budget_and_error_side(self):
        def signed(t, side):
            v, e = box(t, side)
            v[[2, 9]] = np.where((t < 3.5) | ((t == 3.5) & (side == "left")), 0.04, -0.04)
            v[10] += v[9]
            return v, e

        for edges, eta, err, key, side in (
            ([989, 990, 992, 993, 993.125, 994, 996, 1000, 1001], 1e-300, 1.5e-8, "bin[3].hole", "relative"),
            ([989, 990, 992.5, 993.5, 996, 1000, 1001], 1e-9, 1e-8, "bin[2].hole", "absolute"),
        ):

            def inject(role, t, v, e):
                e[[0, 8, 10]] += err
                return v, e

            with self.assertRaises(MultipleLossConvergenceError) as why:
                calculate(
                    kernel=signed,
                    max_order=2,
                    output_edges_eV=edges,
                    features={"jumps_eV": [3.5]},
                    area_epsabs=1e-6,
                    spectrum_epsabs_fraction=eta,
                    inject=inject,
                )
            exc = why.exception
            self.assertEqual(exc.reason, "inherited_error_exceeds_order_budget")
            entry = next(x for x in exc.exceeded if x["key"] == key)
            self.assertEqual(entry["side"], side)
            self.assertEqual(entry["suggestion"], "epsrel" if side == "relative" else "spectrum_epsabs_fraction")
            self.assertIn("lowering epsrel", str(exc))
            self.assertIn("1e-7", str(exc))

    def test_reference_models_with_zero_repeat(self):
        def holes(t, side):
            v, e = box(t, side)
            v[0] = 0.0
            v[2] = v[9] = 0.02
            v[3] = v[8]
            v[10] += v[9]
            return v, e

        for model in (models.GEOMETRIC, models.FACTORIAL, models.EXCESS):
            r = models.integrate_orders_reference(
                holes,
                (2.0, 5.0),
                features={},
                weight_length_A=2.0,
                primary=DeltaPrimary([1000], [1]),
                output_edges_eV=[988, 994, 1000],
                model=model,
                max_order=3,
            )
            self.assertTrue(r.reference_only)
            for n, o in enumerate(r.orders, 1):
                scale = 1 if model == models.GEOMETRIC else math.factorial(n)
                self.assertAlmostEqual(o.area["hole"], (0.72**n - 0.6**n) / scale)
            self.assertGreater(r.truncation["neglected_at_or_above_region_bound"][0], 0.0)

    def test_nonfinite_output_and_large_exact_coefficient(self):
        with self.assertRaises(MultipleLossConvergenceError) as why:
            calculate(max_order=2, multiplier="external", multiplier_A=1e200)
        self.assertEqual((why.exception.reason, why.exception.order), ("nonfinite_output", 2))

        def small(t, side):
            v, e = box(t, side)
            return v * 1e-99, e

        r, _ = calculate(kernel=small, max_order=2, multiplier="external", multiplier_A=1e160)
        self.assertAlmostEqual(r.orders[1].area["total"] / 9e120, 1.0, places=12)

    def test_mesh_limits_and_unavailable_tail(self):
        def perturb(n, role, bounds, t, v, e):
            e[-1] += 0.1
            return v, e

        for override, reason in (
            (dict(order_panel_limit=2), "order_budget_not_met"),
            (dict(order_min_panel_width_eV=2.0), "order_mesh_unresolved"),
        ):
            with self.assertRaises(MultipleLossConvergenceError) as why:
                calculate(order_inject=perturb, max_order=2, **override)
            self.assertEqual((why.exception.reason, why.exception.order), (reason, 2))
            self.assertIsNotNone(why.exception.panel_eV)
        with self.assertRaises(MultipleLossConvergenceError) as why:
            calculate(max_order=1, multiplier="external", multiplier_A=1e25, neglect_tolerance_fraction=0.001)
        self.assertEqual(why.exception.reason, "truncation_bound_unavailable")
        with self.assertRaises(MultipleLossConvergenceError) as why:
            calculate(max_order=1, neglect_tolerance_fraction=1e-30)
        self.assertEqual(why.exception.reason, "truncation_bound_exceeds_tolerance")
        self.assertGreater(why.exception.required_order_estimate, 1)
