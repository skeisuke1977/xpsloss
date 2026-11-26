"""Analytic probes of the loss mesh, polynomial quadrature and detector cells."""

from fractions import Fraction
import math
import unittest
from unittest.mock import Mock

import numpy as np

from xpsloss import BinnedPrimary, DeltaPrimary, SingleLossConvergenceError, SingleLossInputError
from xpsloss._fold import Geometry, fold
from xpsloss._loss_panels import Panel, _diagnostics, _initial_mesh, _nodes, integrate_window
from xpsloss._single_loss_types import inverse_area


SETTINGS = dict(features={}, weight_length_A=2., epsrel=1e-6, area_epsabs=1e-10,
                spectrum_epsabs_fraction=1e-9, panel_limit=400, max_cells=2_000_000)


def solve(kernel, window=(2., 5.), **kwargs):
    return integrate_window(kernel, window, **{**SETTINGS, **kwargs})


def columns(values):
    result = np.zeros((11, len(values)))
    result[[0, 8, 10]] = values
    return result


def constant(t, side):
    v = columns(np.full(t.size, .1))
    return v, np.zeros_like(v)


class LossMeshTests(unittest.TestCase):
    def test_exact_panel_counts_before_any_evaluation(self):
        kernel = Mock(side_effect=AssertionError("unexpected evaluation"))
        for width, expected in [(1e-6, 1_500_002), (.0625, 24)]:
            with self.assertRaises(SingleLossInputError) as err:
                solve(kernel, features=dict(min_half_width_eV=width), panel_limit=23)
            self.assertEqual(err.exception.reason, "resolution_exceeds_panel_limit")
            self.assertEqual(err.exception.required_panels, expected)
        kernel.assert_not_called()
        r, _ = solve(constant, features=dict(singularities=[(3.137, 1e-6)]))
        self.assertEqual(r.diagnostics["required_panels"], 41)
        r, _ = solve(constant, features=dict(min_half_width_eV=.0625), panel_limit=24)
        self.assertEqual(r.diagnostics["required_panels"], 24)
        self.assertIsNone(r.diagnostics["resolution_basis"])
        self.assertEqual(r.diagnostics["resolution_basis_reason"], "analytic_layer")

    def test_node_separation_and_exact_endpoints(self):
        h = math.ulp(2.)
        for width, valid in [(100, False), (256, True), (1, False)]:
            nodes, _, separated = _nodes(2., 2.+width*h)
            self.assertEqual(separated, valid)
            self.assertEqual((nodes[0], nodes[-1]), (2., 2.+width*h))
        kernel = Mock(side_effect=AssertionError)
        with self.assertRaises(SingleLossInputError) as err:
            solve(kernel, (2., 2.+100*h))
        self.assertEqual(err.exception.reason, "invalid_window")
        with self.assertRaises(SingleLossInputError) as err:
            solve(kernel, (2., 2.+256*h), features={"kinks_eV": [2.+64*h, 2.+65*h],
                                                     "min_variation_scale_eV": h/2})
        self.assertEqual(err.exception.reason, "panel_width_unresolved")
        self.assertEqual(err.exception.panel_width_ulp, 1)
        kernel.assert_not_called()

    def test_count_zero_mass_cells_and_both_outsides(self):
        bounds = _initial_mesh((1., 50.), {"min_half_width_eV": 4.}, 400, _diagnostics({}), energy=1000.)
        edges = [bounds[0][0], *(r for _, r in bounds)]
        self.assertEqual(edges, [1, 2, 4, 8, 16, 24, 32, 38, 44, 50])
        geom = Geometry(DeltaPrimary([1000], [1]), [960., 980., 984., 988., 995., 999., 1000., 1001.])
        self.assertEqual(geom.cell_count(edges, 89), 89)
        counts = [9+len(set(geom.breaks(i, Fraction(1), Fraction(50)))-set(edges))
                  for i in range(len(geom.regions))]
        self.assertEqual(counts, [10, 11, 10, 10, 11, 10, 9, 9, 9])

    def test_quadrature_verification_reaches_the_satellite(self):
        def oscillatory(t, side):
            vals = columns(4e-11*(1-np.polynomial.chebyshev.chebval(2*t-5, [0]*32+[1])))
            return vals, np.zeros_like(vals)
        result, _ = solve(oscillatory, (2., 3.), weight_length_A=100., primary=DeltaPrimary([1000], [1]),
                          output_edges_eV=[997, 998], features={"min_variation_scale_eV": .5})
        self.assertGreater(result.diagnostics["stages"], 0)
        self.assertAlmostEqual(result.single_loss_total[0], 4e-9*(1+1/1023), delta=1e-16)

    def test_verification_uncertainty_and_snapshot(self):
        def lorentz(t, side):
            v = columns((1e-10/.7753994963)/(1+((t-2.4)/.5)**2))
            return v, np.zeros_like(v)
        def inject(role, t, v, e):
            return (np.zeros_like(v), v) if role == "interpolation" else (v, v)
        with self.assertRaises(SingleLossConvergenceError) as err:
            solve(lorentz, (2., 3.), inject=inject, panel_limit=1)
        self.assertEqual((err.exception.reason, err.exception.dominant_part),
                         ("loss_budget_not_met", "verification"))
        d = err.exception.diagnostics
        record = d["error_budget"]["outputs"]["window_area.total"]
        worst = d["error_budget"]["max_ratio"]
        row = d["error_budget"]["outputs"][worst["key"]]
        for field in (f"output={worst['key']}", f"estimate={row['estimate']}", f"budget={row['budget']}",
                      f"ratio={worst['value']}"):
            self.assertIn(field, str(err.exception))
        self.assertGreater(record["verification"], record["kernel_point_error"])
        self.assertAlmostEqual(record["estimate"], 3.5585e-10, delta=3.5e-12)
        self.assertGreater(d["verification"]["rho_max"], 0)
        self.assertGreater(d["verification"]["epsilon_max"], 0)

    def test_point_floor_can_fall_after_refinement(self):
        def spike(t, side):
            return columns(np.ones(t.size)), columns(np.exp(-((t-2.5)/.001)**2))
        result, _ = solve(spike, (2., 3.), area_epsabs=.01)
        self.assertGreater(result.diagnostics["stages"], 0)
        self.assertAlmostEqual(result.window_area["total"], 1.)
        def floor(role, t, v, e):
            if role == "interpolation":
                e[:] = 1e-4
            return v, e
        kw = dict(inject=floor, panel_limit=2, initial_edges=[2, 4, 5])
        with self.assertRaises(SingleLossConvergenceError) as err:
            solve(constant, area_epsabs=2e-4, **kw)
        self.assertEqual(err.exception.reason, "kernel_error_exceeds_budget")
        result, _ = solve(constant, area_epsabs=3.5e-4, **kw)
        record = result.diagnostics["error_budget"]["outputs"]["window_area.total"]
        self.assertAlmostEqual(record["kernel_point_error"], 3e-4, delta=3e-16)

    def test_nonfinite_kernel_and_refined_unresolved_panel(self):
        def invalid(t, side):
            return np.full((11, t.size), np.inf), np.zeros((11, t.size))
        with self.assertRaises(SingleLossConvergenceError) as err:
            solve(invalid)
        self.assertEqual(err.exception.reason, "nonfinite_output")
        def noise(t, side):
            return np.zeros((11, t.size)), np.ones((11, t.size))
        with self.assertRaises(SingleLossConvergenceError) as err:
            solve(noise, (2., 2.+256*math.ulp(2.)), area_epsabs=1e-30)
        self.assertEqual(err.exception.reason, "panel_width_unresolved")
        self.assertEqual(err.exception.panel_width_ulp, 128)


class DetectorCellTests(unittest.TestCase):
    def test_affine_mass_times_signed_kernel(self):
        def triangular(t, side):
            inside = ((t > 2) | ((t == 2) & (side == "right"))) & ((t < 5) | ((t == 5) & (side == "left")))
            electron = .1*np.maximum(0., np.minimum(t-2, (5-t)/2))
            hole = .02*(t-3)*inside
            v = columns(electron)
            v[2] = v[9] = hole
            v[10] += .5*hole
            return v, np.zeros_like(v)
        r, _ = solve(triangular, (1., 6.), features=dict(jumps_eV=[2, 5], kinks_eV=[3]),
                     primary=BinnedPrimary([998, 999, 1002], [2, 1]), output_edges_eV=[995, 996, 998, 1000],
                     hole_charge=.5)
        self.assertAlmostEqual(r.window_area["total"], .165)
        self.assertAlmostEqual(r.single_loss_hole[-1], -2/225)
        for item in r.diagnostics["balance"].values():
            self.assertLessEqual(abs(item["difference"]), item["rounding_bound"])

    def test_total_is_integrated_as_an_independent_column(self):
        values = np.zeros((11, 17))
        values[8], values[9], values[10] = .1, .02, .5
        panel = Panel.from_lobatto(2., 5., values, np.zeros_like(values))
        rows = fold([panel], DeltaPrimary([1000], [1]), [990, 1001], multiplier=2., hole_charge=.5).outputs
        self.assertAlmostEqual(rows["bin[0].single_loss_total"]["value"], 3.)
        self.assertAlmostEqual(rows["bin[0].single_loss_hole"]["value"], .12)

    def test_rectangle_and_independent_multiplier(self):
        primary = DeltaPrimary([1000], [1])
        edges = [990, 994, 995, 996.5, 998, 999, 1000, 1001]
        coeff = np.array([.04,.03,.02,.02,-.007,-.01,.005,.002,.085,.015,.1])
        def rectangle(t, side):
            inside = ((t > 2) | ((t == 2) & (side == "right"))) & ((t < 5) | ((t == 5) & (side == "left")))
            v = coeff[:, None]*inside
            return v, np.zeros_like(v)
        result, panels = solve(rectangle, (1., 8.), features=dict(jumps_eV=[2, 5]),
                               primary=primary, output_edges_eV=edges)
        np.testing.assert_allclose(list(result.window_area.values()), 3*coeff, atol=2e-16, rtol=0)
        np.testing.assert_array_equal(result.zero_loss, [0, 0, 0, 0, 0, 0, 1])
        self.assertAlmostEqual(sum(result.single_loss_total), .6)
        folded = fold(panels, primary, edges, multiplier=3., hole_charge=1.)
        for name in ("single_loss_total", "single_loss_no_hole", "single_loss_hole"):
            bins = getattr(result, name)
            original = {f"bin[{i}]": v for i, v in enumerate(bins)}
            original.update({"outside."+side: row[name] for side, row in result.outside.items()})
            for region, val in original.items():
                self.assertAlmostEqual(folded.outputs[f"{region}.{name}"]["value"], 1.5*val,
                                       delta=4*math.ulp(1.)*abs(1.5*val))

    def test_exact_crossing_and_two_one_sided_limits(self):
        p = DeltaPrimary([1000], [1])
        r, _ = solve(constant, primary=p, output_edges_eV=[994, 995, 998, 1001])
        np.testing.assert_allclose(r.cutoff_sensitivity_bins["plus"]["total"], [0, .4, 0])
        np.testing.assert_allclose(r.cutoff_sensitivity_bins["minus"]["total"], [0, 0, .4])
        r, _ = solve(constant, (999.9, 1000.5), primary=p, output_edges_eV=[0, .1, .2])
        for side in ("plus", "minus"):
            np.testing.assert_array_equal(r.cutoff_sensitivity_bins[side]["total"], [0., 2*999.9*.1])

    def test_negative_weight_propagates_absolute_point_error(self):
        vals = np.full((11, 17), .1)
        err = np.zeros_like(vals)
        err[:, 2] = 1e-6
        panel = Panel.from_lobatto(2., 5., vals, err)
        row = fold([panel], DeltaPrimary([1000], [1]), [997.98, 998], multiplier=2.,
                   hole_charge=1.).outputs["bin[0].single_loss_total"]
        self.assertAlmostEqual(row["kernel_point_error"], 2e-6*.0016024037678, delta=1e-18)
        self.assertNotIn("budget", row)

    def test_exact_products_before_density_and_local_overflow_bins(self):
        p = BinnedPrimary([0, 2**900], [2**900])
        result, _ = solve(constant, primary=p, output_edges_eV=[0., 2.**-200])
        self.assertEqual(result.zero_loss[0], 2.**-200)
        self.assertEqual(result.zero_loss_per_eV[0], 1.)
        p = DeltaPrimary([980, 1000, 1030], [1, 2.**60, 2])
        r, _ = solve(constant, primary=p, output_edges_eV=[990, 1002])
        self.assertAlmostEqual(r.outside["low"]["single_loss_total"], .6)
        self.assertAlmostEqual(r.outside["high"]["single_loss_total"], 1.2)

    def test_internal_signed_sum_of_bins_and_lines(self):
        panel = Panel.from_lobatto(2., 3., np.ones((11, 17)), np.zeros((11, 17)))
        mixture = [(1, DeltaPrimary([10], [1])), (-1, BinnedPrimary([9, 11], [3]))]
        rows = fold([panel], mixture, [0, 20], multiplier=2., hole_charge=1.).outputs
        self.assertEqual(rows["bin[0].zero_loss"]["value"], -2.)
        self.assertAlmostEqual(rows["bin[0].single_loss_total"]["value"], -4.)

    def test_reciprocal_bounds_are_tightly_outward(self):
        a, e = .6519427672317586, .06645324737206787
        _, (lo, hi), why = inverse_area(a, e)
        self.assertIsNone(why)
        for x, exact, other in [(lo, 1/(Fraction(a)+Fraction(e)), math.inf),
                                (hi, 1/(Fraction(a)-Fraction(e)), -math.inf)]:
            interior = math.nextafter(x, other)
            self.assertTrue(min(x, interior) <= exact <= max(x, interior))
        self.assertLessEqual(Fraction(lo), 1/(Fraction(a)+Fraction(e)))
        self.assertGreaterEqual(Fraction(hi), 1/(Fraction(a)-Fraction(e)))
        self.assertEqual(inverse_area(1., 1.)[2], "nonpositive_or_unresolved_area")
        self.assertEqual(inverse_area(math.ulp(0.), 0.)[2], "inverse_not_representable")

