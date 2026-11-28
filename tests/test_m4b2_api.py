"""Support and display orientation through the production loss-order path."""

from fractions import Fraction as Q
from itertools import product
import unittest

import numpy as np

import xpsloss as x
from xpsloss import _loss_orders
from xpsloss._binding_axis import convert
from test_m4b2_axis import axis


def box(t, side):
    values = np.zeros((11, len(t)))
    values[[0, 8, 10]] = .2
    return values, np.zeros_like(values)


def spectrum(p, edges, ax, *, window=(1., 2.5), **changes):
    options = dict(features={}, weight_length_A=2., primary=p, output_edges_eV=edges,
                   model=x.MULTIPLE_LOSS_MODEL_ID, max_order=1, energy_axis=ax)
    return _loss_orders.integrate_orders(box, window, **(options | changes))[0]


class BindingSpectrum(unittest.TestCase):
    def test_delta_edges_and_both_outside_regions(self):
        p = x.BindingDeltaPrimary([0., 1., 3., -2., 5.], [1., 2., 4., 8., 16.])
        for edges, mass in (([0., 1., 3.], [3., 4.]), ([3., 1., 0.], [4., 3.])):
            r = spectrum(p, edges, axis(), max_order=2)
            self.assertEqual(r.zero_loss.tolist(), mass)
            for side, value, ke in (("high_binding", 16., "low"), ("low_binding", 8., "high")):
                self.assertEqual(set(r.outside[side]), {"zero_loss", "zero_loss_error", "spectrum", "spectrum_error"})
                self.assertEqual(r.outside[side]["zero_loss"], value)
                self.assertEqual(r.orders[0].outside[side]["total"], r.single_loss.outside[ke]["single_loss_total"])
            self.assertEqual(r.single_loss.zero_loss.tolist(), [4., 3.])
            self.assertEqual(r.single_loss.output_edges_eV.tolist(), [29., 31., 32.])
            self.assertIsNone(r.energy_axis["electron_energy_minus_reference_eV"])
        inc = spectrum(p, [0., 1., 3.], axis(), max_order=2)
        np.testing.assert_array_equal(inc.spectrum[::-1], r.spectrum)
        for a, b in zip(inc.orders, r.orders):
            for c in ("total", "no_hole", "hole"):
                for f in (c, c + "_error"):
                    np.testing.assert_array_equal(getattr(a, f)[::-1], getattr(b, f))
                for k, part in getattr(a, c + "_error_parts").items():
                    np.testing.assert_array_equal(part[::-1], getattr(b, c + "_error_parts")[k])

    def test_binned_primary_uses_converted_widths(self):
        u = 2.**-40
        p = x.BindingBinnedPrimary([2*u, 0., -2*u], [2., 4.])
        r = spectrum(p, [-2*u, 0., .75*u, 2*u], axis(8192.))
        self.assertEqual(r.zero_loss.tolist(), [4., 1., 1.])
        self.assertEqual(r.energy_axis["rounding_max_eV"], u/4)
        self.assertEqual(r.energy_axis["region_loss_upper_eV"], [2*u, 3*u, 4*u])

    def test_both_support_orders_and_bounds_even_without_disagreement(self):
        p = x.BindingDeltaPrimary([0.], [1.])
        a = 1. + 2.**-45
        for d, window, support in ((2.+2.**-43, (a, 2.5), (2, 1)),
                                   (2.-2.**-43, (1., 2.5), (1, 2)),
                                   (2.25+2.**-43, (a, 2.5), (2, 2))):
            r = spectrum(p, [0., 1., d], axis(8192.), window=window)
            info = r.energy_axis
            self.assertEqual((info["support_order_binding"], info["support_order_kinetic"]), support)
            self.assertEqual(r.truncation["state"], "truncated")
            self.assertEqual(r.truncation["bound_reason"], ["complete", None])
            width = max(d, 8192. - (8192. - d))
            self.assertEqual(info["region_loss_upper_eV"], [1., width])
            exact = spectrum(p, [0., 1., width], axis(16.), window=window)
            for k in ("neglected_at_or_above_region_bound", "bound_reason", "tau_per_eV"):
                self.assertEqual(r.truncation[k], exact.truncation[k])
            self.assertGreater(r.truncation["neglected_at_or_above_region_bound"][1], 0.)
            codes = [w["code"] for w in r.diagnostics["warnings"]]
            self.assertEqual(codes.count("axis_conversion_support_disagreement"), int(support[0] != support[1]))
        d = 2.+2.**-43
        accepted = spectrum(p, [0., 1., d], axis(8192.), window=(a, 2.5), neglect_tolerance_fraction=.01)
        self.assertEqual(accepted.truncation["state"], "bounded_on_display")
        with self.assertRaises(x.MultipleLossConvergenceError) as cm:
            spectrum(p, [0., 1., d], axis(8192.), window=(a, 2.5), neglect_tolerance_fraction=1e-6)
        self.assertEqual(cm.exception.required_order_estimate, 2)
        self.assertEqual(cm.exception.diagnostics["truncation"]["bound_reason"], ["complete", None])

    def test_exact_support_comparison_before_upward_rounding(self):
        a = 1.+2.**-52
        edges = [-1., 2*a]
        # At nonzero p the exact width lies just below 2a; its ceiling is 2a.
        for pos, reason in ((2.**-55, "complete"), (0., None)):
            r = spectrum(x.BindingDeltaPrimary([pos], [1.]), edges, axis(16.), window=(a, 2.5))
            self.assertEqual(r.energy_axis["region_loss_upper_eV"], [2*a])
            self.assertEqual(r.truncation["bound_reason"], [reason])

    def test_bisection_counts_and_eight_examples_against_threshold_relations(self):
        positions = [float(8*i) for i in range(9, -1, -1)]
        edges = sorted((p + k + 2.**-43 for p in positions for k in range(1, 5)), reverse=True)
        p = x.BindingDeltaPrimary(positions, [1.]*len(positions))
        display = convert(axis(8192.), p, edges)
        display.support_order(1.)
        display.compare((1., 2.), 2)
        obs = display.info["support_comparison"]
        expected = []
        for n, (bound, t), point, edge in product((1, 2), (("lower", 1), ("upper", 2)), sorted(positions), sorted(edges)):
            rb = Q(edge)-Q(point)-n*t
            rk = Q(8192.-point)-Q(8192.-edge)-n*t
            if (rb > 0)-(rb < 0) != (rk > 0)-(rk < 0):
                expected.append((n, bound, positions.index(point), edges.index(edge)))
        self.assertGreater(len(expected), 8)
        self.assertEqual(obs["computed_relation_mismatch_count"], len(expected))
        self.assertEqual([(v["order"], v["bound"], v["primary_index"], v["edge_index"])
                          for v in obs["first_mismatches"]], expected[:8])
        self.assertEqual(obs["next_relation_mismatch_count"], 1)
        self.assertEqual(obs["next_reachability_mismatch_bins"], [])  # equal vs greater, both reached

    def test_next_order_example_and_warning_sequence(self):
        r = spectrum(x.BindingDeltaPrimary([0.], [1.]), [9., 2.+2.**-43, 0.],
                     axis(8192., emission_angle_deg=30.), window=(1.+2.**-45, 2.5))
        diag = r.energy_axis["support_comparison"]
        self.assertEqual(diag["next_reachability_mismatch_bins"], [1])
        ex, = diag["first_mismatches"]
        self.assertEqual((ex["comparison"], ex["edge_index"], ex["bin_index"], ex["binding_relation"],
                          ex["kinetic_relation"]), ("next_order", 1, 1, ">", "<"))
        warnings = r.diagnostics["warnings"]
        self.assertEqual(tuple(w["message"] for w in warnings), r.warnings)
        self.assertEqual([w["code"] for w in warnings[-4:]], ["binding_axis_assumed_calibrated",
                         "energy_axis_field_not_recorded", "recorded_emission_angle_differs_from_theta",
                         "axis_conversion_support_disagreement"])

    def test_display_support_and_validation_order(self):
        p = x.BindingDeltaPrimary([0.], [1.])
        kwargs = dict(primary=p, output_edges_eV=[0., 2.+2.**-43], max_order="display_support",
                      features={}, weight_length_A=1., energy_axis=axis(8192.))
        for model, limit, expected in (("unknown", 1, "invalid_model"), (x.MULTIPLE_LOSS_MODEL_ID, 1, "order_limit_exceeded")):
            with self.assertRaises(x.MultipleLossInputError) as cm:
                _loss_orders.integrate_orders(self.fail, (1.+2.**-45, 2.5), model=model, order_limit=limit, **kwargs)
            self.assertEqual(cm.exception.reason, expected)
            self.assertNotIn("detail", vars(cm.exception))
        r = spectrum(x.BindingDeltaPrimary([0.], [0.]), [0., 4.], axis(), max_order="display_support")
        self.assertEqual(r.orders, ())
        self.assertEqual(r.energy_axis["region_loss_upper_eV"], [None])
        self.assertEqual(r.truncation["bound_reason"], ["complete"])
