"""Public contract, moving-endpoint examples and the M3a boundary."""

import math
from unittest import TestCase
from unittest.mock import Mock, patch

import numpy as np

from xpsloss import (BinnedPrimary, DeltaPrimary, DepthAverageConvergenceError, Material, Oscillator,
                     SingleLossConvergenceError, SingleLossInputError, loss_window_area, single_loss_spectrum)
from xpsloss import xps_single_loss as api
from xpsloss import _loss_panels, _loss_roots
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
        enumeration = candidates(m, 1000., mode, w)
        roots = enumeration["candidates"]
        f = dict(singularities=[(r["z_re"], r["z_im"]) for r in roots],
                 min_half_width_eV=min(o.width_eV for o in m.oscillators)/2)
        mesh = _initial_mesh(w, f, 1000, _diagnostics({}), energy=1000., gap=m.gap_eV)
        return len(mesh), roots, enumeration["resolution_basis"]

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
        # CE5 needs 4 panels: its surface roots lie on the clearance boundary (prereg 11).
        for energy, n, expected, confirmed in [(10.05, 4, "fixed_poles_enumerated", True),
                                                (10.0000001, 4, "computed_candidates_only", False)]:
            mat = Material("synthetic_doublet", 0., [Oscillator(10, 30, 1), Oscillator(energy, 30, 1)])
            count, roots, basis = self.mesh(mat, (7., 11.), "k_total")
            self.assertEqual((count, basis), (n, expected))
            positive = [r for r in roots if r["family"] == "surface" and 7 < r["z_re"] < 11]
            self.assertEqual(len(positive), 2)
            self.assertEqual(max(positive, key=lambda r: r["z_re"])["confirmed"], confirmed)

    def test_clearance_boundary_does_not_depend_on_the_last_bits(self):
        # prereg 11: a pole at distance |P|/2 from a panel of the cap width, or
        # within a relative 2**-48 inside that boundary, does not split it.
        def count(imag):
            f = dict(singularities=[(3., imag)], min_half_width_eV=.5)
            return len(_initial_mesh((2., 4.), f, 1000, _diagnostics({})))
        self.assertEqual(count(-.5), 2)
        self.assertEqual(count(math.nextafter(-.5, 0.)), 2)
        self.assertEqual(count(-.5*(1-2.**-49)), 2)
        self.assertEqual(count(-.5*(1-2.**-47)), 4)

    def test_remote_rational_check_is_visible_outside_mesh_candidates(self):
        material = Material("synthetic", 0, [Oscillator(1, .5, .1, 3037.62221442588)])
        window = (49.5, 50.5)
        d = candidates(material, 1000., "k_parallel", window)
        self.assertEqual(d["resolution_basis"], "endpoint_roots_enumerated")
        remote = d["remote_root_checks"]
        self.assertTrue(any(r["confirmed"] and r["family"] == "surface" for r in remote))
        self.assertTrue(all(r["method"] == "rational_remote_check" and r["distance_to_window_eV"] >= 1
                            for r in remote))
        self.assertTrue(d["candidates"])
        self.assertTrue(all(r["method"] == "enumerated" and r["distance_to_window_eV"] < 1
                            for r in d["candidates"]))

    def test_numeric_failure_is_local_to_one_family_and_endpoint(self):
        material = Material("synthetic", 0, [Oscillator(20, 3, 3, .5)])
        args = (material, 1000., "k_total", (1., 50.))
        def unaffected(r):
            return (r["family"], r["end"]) != ("surface", "lower")
        expected = [r for r in candidates(*args)["candidates"] if unaffected(r)]
        self.assertTrue(any(r["end"] == "lower" for r in expected))
        build = _loss_roots._family_polynomial
        bad_roots = Mock(coef=np.array([1., 1.]))
        bad_roots.trim.return_value.roots.return_value = np.array([complex(float("nan"), 0)])
        for failure in [TypeError("polynomial arithmetic"), np.linalg.LinAlgError("root solver"),
                        _loss_roots.Poly([1., np.inf]), bad_roots]:
            with self.subTest(failure=failure):
                def substitute(lines, energy, mode, end, family, index):
                    if (family, end) == ("surface", "lower"):
                        if isinstance(failure, Exception):
                            raise failure
                        return failure, _loss_roots.Poly([0., 1.])
                    return build(lines, energy, mode, end, family, index)
                with patch.object(_loss_roots, "_family_polynomial", side_effect=substitute):
                    d = candidates(*args)
                self.assertEqual(d["candidates"], expected)
                self.assertEqual((d["resolution_basis"], d["resolution_basis_reason"]),
                                 ("computed_candidates_only", "enumeration_failed"))
                warning, = d["warnings"]
                self.assertEqual((warning["family"], warning["end"], warning["code"]),
                                 ("surface", "lower", "enumeration_failed"))
                self.assertIn(warning["reason"], warning["message"])

    def test_many_oscillators_defer_input_rejection_to_m3a(self):
        lines = [Oscillator(c, a, w, (1., .5)[i % 2]) for i, (c, a, w) in
                 enumerate(zip(np.linspace(3, 50, 60), np.linspace(5, 80, 60), np.linspace(1, 6, 60)))]
        material = Material("synthetic_many", 0., lines)
        with self.assertRaisesRegex(ValueError, "static sign") as upstream:
            api.effective_cross_section_depth_averaged(material, 1000., [1., 50.],
                                                       epsilon_momentum="k_total", weight_length_A=20)
        with self.assertRaises(type(upstream.exception)) as caught:
            loss_window_area(material, 1000., **KW)
        self.assertEqual(str(caught.exception), str(upstream.exception))

    def test_many_oscillators_return_area_after_coupled_polynomial_failure(self):
        for n, energy in [(60, 10000.), (30, 100000.)]:
            with self.subTest(n=n, energy=energy):
                material = Material("synthetic_many", 0, [Oscillator(c, 3., 3., 1.) for c in np.linspace(6, 50, n)])
                wide = candidates(material, energy, "k_total", (1., 50.))
                self.assertEqual(wide["resolution_basis_reason"], "enumeration_failed")
                self.assertTrue(any(r["family"] == "bulk" for r in wide["candidates"]))
                # A narrow window keeps the real M3a integration inexpensive in the routine suite.
                window = (1., 1.01)
                kw = dict(epsilon_momentum="k_total", weight_length_A=20.)
                kernel = api.effective_cross_section_depth_averaged(material, energy, window, limit=3000, **kw)
                result = loss_window_area(material, energy, loss_window_eV=window, kernel_limit=3000, **kw)
                mean = result.window_area["total"]/(window[1]-window[0])
                self.assertLess(min(kernel.total), mean)
                self.assertLess(mean, max(kernel.total))
                d = result.diagnostics
                self.assertEqual((d["resolution_basis"], d["resolution_basis_reason"]),
                                 ("computed_candidates_only", "enumeration_failed"))
                failed = {(w["family"], w["end"]) for w in d["warnings"] if w["code"] == "enumeration_failed"}
                self.assertEqual(failed, {(f, e) for f in ("surface", "epsilon") for e in ("lower", "upper")})
                self.assertTrue(all(w["message"] in result.warnings for w in d["warnings"]))

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
    def test_integer_settings_use_the_index_protocol(self):
        integers = dict(panel_limit=np.int64(400), kernel_limit=np.int32(200), max_cells=np.uint64(1000))
        primary = DeltaPrimary([1000.], [1.])
        r = single_loss_spectrum(WIDE, 1000., primary, EDGES, **{**KW, "loss_window_eV": (5., 5.01)}, **integers)
        for name, value in integers.items():
            actual = r.diagnostics["effective_settings"][name]
            self.assertIs(type(actual), int)
            self.assertEqual(actual, value)
            for invalid in (True, np.bool_(True), np.int64(0), 1.5):
                with self.subTest(name=name, invalid=invalid), self.assertRaises(SingleLossInputError) as err:
                    single_loss_spectrum(ZERO, 1000., primary, EDGES, **KW, **{name: invalid})
                self.assertEqual(err.exception.reason, "invalid_settings")
                self.assertIsNone(err.exception.diagnostics)

    def test_m3a_inherits_callers_float_settings(self):
        settings = dict(over="ignore", divide="warn", invalid="ignore", under="ignore")
        kernel, areas = api._call_kernel, _loss_panels._areas
        def sample(*args, **kwargs):
            self.assertEqual(np.geterr(), settings)
            return kernel(*args, **kwargs)
        def integrate(panels):
            self.assertEqual(np.geterr(), {**settings, "over": "raise", "divide": "raise", "invalid": "raise"})
            return areas(panels)
        with np.errstate(**settings), patch.object(api, "_call_kernel", side_effect=sample) as m3a, \
                patch.object(_loss_panels, "_areas", side_effect=integrate) as m4a:
            loss_window_area(WIDE, 1000., **{**KW, "loss_window_eV": (5., 5.01)})
            self.assertEqual(np.geterr(), settings)
        m3a.assert_called()
        m4a.assert_called()

    def test_common_validation_precedes_kernel_even_for_zero_material(self):
        with patch.object(api, "_call_kernel", side_effect=AssertionError) as calls:
            for kw, reason in [({"loss_window_eV": [1, float("nan")]}, "invalid_window"),
                               ({"loss_window_eV": [1, 2**20_000]}, "invalid_window"),
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
        for field in (f"loss_eV={ex.loss_eV}", "index=1", "estimate=2.0", "budget=1.0"):
            self.assertIn(field, str(ex))
        self.assertTrue(str(ex).isascii())

