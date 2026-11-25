"""Public contract, moving-endpoint examples and the M3a boundary."""

from unittest import TestCase
from unittest.mock import patch

import numpy as np

from xpsloss import (BinnedPrimary, DeltaPrimary, DepthAverageConvergenceError, Material, Oscillator,
                     SingleLossConvergenceError, SingleLossInputError, loss_window_area, single_loss_spectrum)
from xpsloss import xps_single_loss as api
from xpsloss._loss_panels import _initial_mesh, _diagnostics
from xpsloss._loss_roots import candidates
from xpsloss.config import load_config
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
WIDE = Material("synthetic_wide", 0, [Oscillator(20, 200, 8)])
ZERO = Material("synthetic_zero", 2, [Oscillator(3, 0, .02)])
KW = dict(epsilon_momentum="k_total", weight_length_A=20, loss_window_eV=(1, 50))
EDGES = [960, 980, 984, 988, 995, 999, 1000, 1001]


class ResolutionTests(TestCase):
    def mesh(self, m, w, mode):
        roots, basis, _ = candidates(m, 1000., mode, w)
        f = dict(singularities=[(r["z_re"], r["z_im"]) for r in roots],
                 min_half_width_eV=min(o.width_eV for o in m.oscillators)/2)
        mesh = _initial_mesh(w, f, 1000, _diagnostics({}), energy=1000., gap=m.gap_eV)
        return len(mesh), roots, basis

    def test_moving_endpoint_roots_and_kinematic_endpoint_grading(self):
        for osc, mode, window, count, family, distance in [
                ((3, 3450, .1, 96), "k_total", (49.5, 50.5), 14, "surface", .0222),
                ((1, .5, .1, 3037.62221442588), "k_parallel", (49.5, 50.5), 12, "bulk", .02536),
                ((5, .5, .1, 1.21604938271605), "k_total", (989.5, 990.5), 18, "bulk", .00503)]:
            with self.subTest(mode=mode, osc=osc):
                n, roots, basis = self.mesh(Material("synthetic", 0., [Oscillator(*osc)]), window, mode)
                self.assertEqual(n, count)
                if count != 18:
                    self.assertEqual(basis, "endpoint_roots_enumerated")
                self.assertTrue(any(r["end"] == "lower" and r["family"] == family and
                                    abs(r["distance_to_window_eV"]/distance-1) < .01 for r in roots))
        n, _, basis = self.mesh(Material("synthetic", 0, [Oscillator(15.6, 243.3598, 1.5)]),
                                (900., 999.9), "k_total")
        self.assertEqual((n, basis), (74, "fixed_poles_enumerated"))

    def test_near_degenerate_surface_roots_are_preserved(self):
        for energy, n, expected, confirmed in [(10.05, 5, "fixed_poles_enumerated", True),
                                                (10.0000001, 4, "computed_candidates_only", False)]:
            mat = Material("synthetic_doublet", 0., [Oscillator(10, 30, 1), Oscillator(energy, 30, 1)])
            count, roots, basis = self.mesh(mat, (7., 11.), "k_total")
            self.assertEqual((count, basis), (n, expected))
            positive = [r for r in roots if r["family"] == "surface" and 7 < r["z_re"] < 11]
            self.assertEqual(len(positive), 2)
            self.assertEqual(max(positive, key=lambda r: r["z_re"])["confirmed"], confirmed)

    def test_gap_grading_for_both_momenta(self):
        for mode in ("k_total", "k_parallel"):
            for name, window, count, expected_basis in [
                    ("al", (1., 50.), 35, "endpoint_roots_enumerated"),
                    ("si", (5., 50.), 15, "computed_candidates_only")]:
                m = load_config(ROOT/"benchmarks"/"inputs"/f"{name}_elfdb_v1.0.bulk.json").material
                n, _, basis = self.mesh(m, window, mode)
                self.assertEqual((n, basis), (count, expected_basis))
            for width, expected in [(.02, 401), (.2, 45)]:
                m = Material("synthetic", 2., [Oscillator(3., 6., width, .5)])
                n, _, _ = self.mesh(m, (2.1, 10.), mode)
                self.assertEqual(n, expected)


class PublicTests(TestCase):
    def test_common_validation_precedes_kernel_even_for_zero_material(self):
        with patch.object(api, "_call_kernel", side_effect=AssertionError) as calls:
            for kw, reason in [({"loss_window_eV": [1, float("nan")]}, "invalid_window"),
                               ({"epsrel": 1e-17}, "invalid_settings"), ({"panel_limit": True}, "invalid_settings")]:
                with self.assertRaises(SingleLossInputError) as exc:
                    loss_window_area(ZERO, 1000, **{**KW, **kw})
                self.assertEqual(exc.exception.reason, reason)
            with self.assertRaises(ValueError):
                loss_window_area(ZERO, 1000, **KW, theta_deg=1.)
            for primary, edges, reason in [(None, EDGES, "invalid_primary"),
                                           (DeltaPrimary([1000], [1]), [1, 1], "invalid_output_edges")]:
                with self.assertRaises(SingleLossInputError) as exc:
                    single_loss_spectrum(WIDE, 1000, primary, edges, **KW)
                self.assertEqual(exc.exception.reason, reason)
            with self.assertRaises(SingleLossInputError) as exc:
                single_loss_spectrum(WIDE, 1000, DeltaPrimary([1000], [1]), EDGES, **KW, max_cells=88)
            self.assertEqual(exc.exception.reason, "cell_limit_exceeded")
            calls.assert_not_called()

    def test_primary_input_ranges_and_types(self):
        for constructor, positions, weights in [(DeltaPrimary, [1, 2], [1]), (BinnedPrimary, [1, 1], [1]),
                                                (DeltaPrimary, [1j], [0]), (DeltaPrimary, [1], [-1]),
                                                (DeltaPrimary, [1], [2.**-901]),
                                                (DeltaPrimary, [1, 2], [2.**900, 2.**-900])]:
            with self.assertRaises(SingleLossInputError):
                constructor(positions, weights)

    def test_zero_response_and_empty_primary_support(self):
        with patch.object(api, "_call_kernel", side_effect=AssertionError):
            primary = DeltaPrimary([999, 1000, 1001], [1, 0, 2])
            result = single_loss_spectrum(ZERO, 1000, primary, [999, 1000], **KW)
            np.testing.assert_array_equal(result.zero_loss, [1.])
            self.assertEqual(result.outside["high"]["spectrum"], 2.)
            np.testing.assert_array_equal(result.single_loss_total, 0.)
            self.assertEqual(result.inverse_reason, "nonpositive_or_unresolved_area")
            self.assertEqual(result.diagnostics["resolution_basis_reason"], "zero_response")
            self.assertEqual(result.diagnostics["kernel_calls"]["calls"], 0)
            empty = single_loss_spectrum(ZERO, 1000, DeltaPrimary([1008], [0]), [998, 1009], **KW)
            self.assertEqual(set(empty.diagnostics["primary_support"].values()), {None})
            row = empty.diagnostics["error_budget"]["outputs"]["bin[0].zero_loss"]
            self.assertEqual((row["ratio"], row["flag"]), (None, "zero_budget_zero_error"))

    def test_real_kernel_limit_and_nonconvergence_record(self):
        with self.assertRaisesRegex(ValueError, "integration limit too small"):
            loss_window_area(WIDE, 1000, **KW, kernel_limit=1)
        with self.assertRaises(SingleLossConvergenceError) as exc:
            loss_window_area(WIDE, 1000, **{**KW, "loss_window_eV": [5., 6.]},
                             kernel_epsabs=1e-30, kernel_epsrel=1e-15)
        ex = exc.exception
        self.assertEqual((ex.reason, ex.inner_reason, ex.node_kind, ex.index, ex.stage),
                         ("kernel_not_converged", "panel_floor_exceeds_budget", "sentinel", 0, 0))
        self.assertEqual(ex.diagnostics["error_budget"], {"outputs": {}, "max_ratio": None})

    def test_cells_after_refinement_and_inherited_warning_provenance(self):
        with self.assertRaises(SingleLossConvergenceError) as exc:
            single_loss_spectrum(WIDE, 1000, DeltaPrimary([1000], [1]), EDGES, **KW, max_cells=89)
        self.assertEqual(exc.exception.reason, "cell_limit_reached")
        d = exc.exception.diagnostics
        self.assertEqual(d["initial_cells"], 89)
        self.assertGreater(d["error_budget"]["max_ratio"]["value"], 1.)
        self.assertEqual(d["kernel_calls"]["batch_sizes"], [32*9+1])
        for row in d["warnings"]:
            if row["source"] == "xps_depth_averaged":
                self.assertEqual(len(row["occurrences"][0]["loss_eV"]), 32*9+1)

    def test_each_batch_has_both_sentinels(self):
        def fail(mat, energy, grid, **kw):
            self.assertEqual((grid[0], grid[-1]), (1, 50))
            self.assertEqual(len(grid), 3)
            self.assertEqual(kw["weight_length_A"], 20)
            raise DepthAverageConvergenceError("budget_not_met", grid[1], 1, column="total", estimate=2., budget=1.)
        with patch.object(api, "_call_kernel", side_effect=fail):
            with self.assertRaises(SingleLossConvergenceError) as caught:
                api._single_loss_spectrum(WIDE, 1000, DeltaPrimary([1000], [1]), EDGES, **KW, kernel_batch_size=1)
        ex = caught.exception
        self.assertEqual((ex.node_kind, ex.index, ex.column, ex.stage), ("interior", 1, "total", 0))
        self.assertEqual(ex.inner_attributes["estimate"], 2.)

