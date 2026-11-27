"""First-order compatibility and structural zero cases at the public entry points."""

from fractions import Fraction as Q
import unittest
from unittest.mock import patch

import numpy as np

import xpsloss as x
from xpsloss import xps_single_loss as legacy
from xpsloss import _loss_orders as layers


class PublicOrders(unittest.TestCase):
    def setUp(self):
        self.zero = x.Material("zero", 0.0, (x.Oscillator(20.0, 0.0, 8.0, 0.0),))
        self.mat = x.Material("synthetic", 0.0, (x.Oscillator(20.0, 200.0, 8.0, 0.0),))
        self.args = dict(loss_window_eV=(2.0, 5.0), weight_length_A=20.0, epsilon_momentum="k_total")
        self.p = x.DeltaPrimary([1000.0], [1.0])

    def test_zero_needs_no_high_order_work(self):
        with patch.object(legacy, "_call_kernel", side_effect=AssertionError):
            r = x.multiple_loss_spectrum(
                self.zero,
                1000,
                self.p,
                [990.0, 999.0, 1001.0],
                **self.args,
                model=x.MULTIPLE_LOSS_MODEL_ID,
                max_order=4,
                max_convolution_points=1,
            )
        self.assertEqual(r.diagnostics["resources"]["quadrature_points_total"], 0)
        np.testing.assert_array_equal(r.spectrum, r.zero_loss)
        for order in r.orders:
            self.assertTrue(all(kind == "exact_zero" for kind in order.area_zero_kind.values()))
            self.assertFalse(any(order.area_error.values()))
            self.assertFalse(order.total_error.any())

    def test_first_order_delivery_and_calls(self):
        with patch.object(legacy, "_call_kernel", wraps=legacy._call_kernel) as calls:
            old = x.single_loss_spectrum(self.mat, 1000, self.p, [992.0, 998.0, 1001.0], **self.args)
            points = [c.args[2].tobytes() for c in calls.call_args_list]
            calls.reset_mock()
            with patch.object(legacy, "_run", wraps=legacy._run) as driver:
                r = x.multiple_loss_spectrum(
                    self.mat,
                    1000,
                    self.p,
                    [992.0, 998.0, 1001.0],
                    **self.args,
                    model=x.MULTIPLE_LOSS_MODEL_ID,
                    max_order=1,
                )
            self.assertEqual([c.args[2].tobytes() for c in calls.call_args_list], points)
            self.assertEqual(driver.call_count, 1)
            self.assertTrue(callable(driver.call_args.kwargs["on_final_panels"]))
        self.assertEqual(r.spectrum_error.tobytes(), old.spectrum_error.tobytes())
        self.assertEqual(r.orders[0].total.tobytes(), old.single_loss_total.tobytes())
        self.assertEqual(r.single_loss.warnings, old.warnings)
        self.assertEqual(r.diagnostics["first_order"]["beta_deliveries"], 1)

    def test_area_conventions(self):
        for c, m, kind in (
            (0.0, "weight_length", "exactly_one"),
            (0.1, "weight_length", "below_one"),
            (0.1, "external", "not_determined_by_convention"),
        ):
            r = x.multiple_loss_order_areas(
                self.mat,
                1000,
                loss_window_eV=(2.0, 5.0),
                epsilon_momentum="k_total",
                weight_length_A=None,
                weight_length_kind="inverse_13a_window_area",
                extra_inverse_length_per_A=c,
                model=x.MULTIPLE_LOSS_MODEL_ID,
                max_order=1,
                multiplier=m,
                multiplier_A=2.0 if m == "external" else None,
            )
            self.assertEqual(r.convention_ratio, kind)
            self.assertFalse(hasattr(r, "truncation"))
            single = r.single_loss
            self.assertEqual(
                r.order_area_ratio["value"], float(Q(single.multiplier_A) * Q(single.window_area["term_13a"]))
            )
            self.assertEqual(r.orders[0].mesh["cells"], 0)

    def test_display_order_zero(self):
        r = x.multiple_loss_spectrum(
            self.mat,
            1000,
            self.p,
            [1000.0, 1001.0],
            **self.args,
            model=x.MULTIPLE_LOSS_MODEL_ID,
            max_order="display_support",
        )
        self.assertEqual(r.orders, ())
        self.assertEqual(r.truncation["state"], "complete_on_display")
        self.assertEqual(r.truncation["computed"], 0)
        self.assertGreater(r.truncation["neglected_total_area"]["value"], 0.0)

    def test_production_panels_pass_unchanged_to_order_two(self):
        delivered, used = [], []
        original_run, original_envelopes = legacy._run, layers.first_envelopes

        def run(*args, on_final_panels, **kwargs):
            def capture(single, panels):
                delivered.append(panels)
                on_final_panels(single, panels)

            return original_run(*args, **kwargs, on_final_panels=capture)

        def envelopes(panels):
            used.append(panels)
            return original_envelopes(panels)

        with patch.object(legacy, "_run", run), patch.object(layers, "first_envelopes", envelopes):
            r = x.multiple_loss_order_areas(
                self.mat, 1000, **self.args, model=x.MULTIPLE_LOSS_MODEL_ID, max_order=2
            )
        self.assertEqual(len(delivered), 1)
        self.assertIs(delivered[0], used[0])
        self.assertEqual(r.orders[1].mesh["formal_folds"], 0)

    def test_old_exception_survives(self):
        with self.assertRaises(x.SingleLossInputError) as why:
            x.multiple_loss_order_areas(
                self.zero,
                1000,
                **self.args,
                model=x.MULTIPLE_LOSS_MODEL_ID,
                max_order=1,
                multiplier="inverse_window_area",
            )
        self.assertEqual(why.exception.reason, "multiplier_undefined")

    def test_first_order_validation_precedes_order_settings(self):
        bad_edges = [999., 980.]
        for invalid in ({"theta_deg": 5.}, {"weight_length_A": -1.}, {"loss_window_eV": (0., 5.)},
                        {"epsrel": -1.}, {"multiplier": "external"}, {}):
            kw = self.args | invalid
            with self.subTest(invalid=invalid), patch.object(legacy, "_call_kernel", side_effect=AssertionError):
                with self.assertRaises(ValueError) as old:
                    x.single_loss_spectrum(self.mat, 1000, self.p, bad_edges, **kw)
                with self.assertRaises(ValueError) as current:
                    x.multiple_loss_spectrum(self.mat, 1000, self.p, bad_edges,
                                             **kw, model="invalid", max_order=0)
                self.assertIs(type(current.exception), type(old.exception))
                self.assertEqual(str(current.exception), str(old.exception))
                self.assertEqual(vars(current.exception), vars(old.exception))
