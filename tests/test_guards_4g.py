"""Public input and convergence witnesses from the phase 4g audits."""
from dataclasses import replace
import math
import unittest
from unittest.mock import patch

import numpy as np

import test_numerics_followup3 as frozen_followup
import test_numerics_4a as frozen_4a
import test_numerics_4f as frozen_4f
from xpsloss import (BulkConvergenceError, DielectricConvergenceError, FixedDepthConvergenceError,
                     Material, NumericalConvergenceError, Oscillator, bulk_diimfp, dielectric_response,
                     effective_cross_section_fixed_depth as fixed_depth, loss_function,
                     negative_inverse_dielectric_function as inverse)
from xpsloss.config import parse_config
from xpsloss.dielectric import (negative_inverse_dielectric_closed_form as closed,
                                optical_elf_first_moment, optical_sum_rule_lhs, tabulated_elf_kk)
from xpsloss.xps_fixed_depth import _coefficients, _include_node_roundoff, _integrate_loss, _log_breakpoints


class InputGuardTests(unittest.TestCase):
    def test_complex_coordinates_and_parameters(self):
        mat = Material("synthetic", 0, (Oscillator(15, 180, 2),))
        for z in (2+9j, np.complex128(2+9j), np.array([2+9j]), np.array([np.complex128(2+9j)], dtype=object)):
            calls = [lambda: bulk_diimfp(mat, z, 1), lambda: bulk_diimfp(mat, 1000, z),
                     lambda: fixed_depth(mat, z, 1, 10, epsilon_momentum="k_total"),
                     lambda: fixed_depth(mat, 1000, z, 10, epsilon_momentum="k_total"),
                     lambda: fixed_depth(mat, 1000, 1, z, epsilon_momentum="k_total")]
            for function in (loss_function, inverse, dielectric_response, closed):
                calls.extend((lambda f=function: f(mat, z, 0), lambda f=function: f(mat, 1, z)))
            for key in ("energy_eV", "strength_eV2", "width_eV", "dispersion"):
                calls.append(lambda k=key: replace(mat.oscillators[0], **{k: z}))
            for index, call in enumerate(calls):
                with self.subTest(z=z, boundary=index), self.assertRaises(ValueError):
                    call()

        # With the imaginary part discarded, these tables would have valid
        # shapes, increasing nodes and a finite transform at a point outside the table.
        for dtype in (complex, object):
            values = np.array([0, np.complex128(1+2j)], dtype=dtype)
            for nodes, elf in ((values, [0, 1]), ([0, 1], values)):
                with self.subTest(dtype=dtype, nodes=nodes), self.assertRaisesRegex(ValueError, "real"):
                    tabulated_elf_kk(nodes, elf, 2, low_energy="zero", high_energy="hard_cut")

    def test_grid_preserves_requested_endpoints(self):
        config = {"material": {"name": "synthetic", "gap_eV": 0, "oscillators": [
            {"energy_eV": 15, "strength_eV2": 180, "width_eV": 2, "dispersion": 0}]},
            "electron_energy_eV": 1000}
        for hi, step, expected in ((1, 1e12, [1]), (2, 1, [1, 2]), (2, 1e12, None)):
            with self.subTest(hi=hi, step=step):
                config["loss_grid_eV"] = {"min": 1, "max": hi, "step": step}
                if expected is None:
                    with self.assertRaises(ValueError):
                        parse_config(config)
                else:
                    np.testing.assert_array_equal(parse_config(config).grid(), expected)

    def test_positive_damping_underflow(self):
        t = 5e-101
        material = Material("synthetic", 0, (Oscillator(1e-100, 1e-30, 1e-240),))
        # Exact binary64 inputs give ELF 8.888888888888888e29. A zero is spurious.
        with self.assertRaises(NumericalConvergenceError):
            loss_function(material, t)
        with self.assertRaises(BulkConvergenceError):
            bulk_diimfp(material, 1000, t)
        with self.assertRaises(DielectricConvergenceError):
            inverse(material, t)
        # Negative static F and resolvable momentum limits isolate damping in fixed depth.
        mat = Material("synthetic", 0, (Oscillator(15, 1, 1e-300),))
        for loss in (1e-30, [1e-30, 2e-30]):
            with self.subTest(loss=loss), self.assertRaises(FixedDepthConvergenceError):
                fixed_depth(mat, 1000, loss, 10, epsilon_momentum="k_parallel")

    def test_fixed_depth_static_sign(self):
        material = Material("synthetic", 0, (Oscillator(10, 110, 1, 0.5),))
        with self.assertRaisesRegex(ValueError, "static sign"):
            fixed_depth(material, 1000, 3, 10, epsilon_momentum="k_parallel")

    def test_fixed_depth_only_needs_static_inverse(self):
        mat = Material("synthetic", 0, (Oscillator(5, 16, 1e-309),))
        # F(0)=-0.36; at T=3 direct epsilon overflows, but S and S/(2-S) are finite.
        r = fixed_depth(mat, 1000, 3, 10, epsilon_momentum="k_parallel")
        self.assertTrue(np.isfinite(r.total[0]))
        self.assertGreater(r.term_13a[0], 0)

    def test_gapless_optical_moments_at_tight_quadrature_tolerance(self):
        mat = Material("synthetic", 0, (Oscillator(3, 1, 1),))
        self.assertEqual(optical_sum_rule_lhs(mat, epsabs=1e-30, epsrel=1e-30), 1 / 9)
        self.assertEqual(optical_elf_first_moment(mat, epsabs=1e-30, epsrel=1e-30), math.pi / 2)


class BudgetGuardTests(unittest.TestCase):
    def test_bulk_line_and_combined_budgets(self):
        mat = Material("synthetic", 0, tuple(Oscillator(15, a, 2) for a in (289.4511329707871, 77.6947777974218)))
        # Compare errors before the 1/(pi*a0*E) normalization applied to both returned fields.
        for ea, er, message in ((1e-30, 1.2e-15, "oscillator 1"), (1.78e-15, 1e-15, "combined")):
            with self.subTest(epsabs=ea), self.assertRaisesRegex(BulkConvergenceError, message):
                bulk_diimfp(mat, 1000, 20, epsabs=ea, epsrel=er)

    def test_small_surface_term(self):
        mat = Material("synthetic", 0, (Oscillator(15, 180, 3e-7),))
        r = fixed_depth(mat, 1000, math.sqrt(105.), 10, hole_charge=0, epsilon_momentum="k_parallel",
                        epsabs=1e-23, epsrel=1e-7)
        # Exact rational Re C times an independent integral of the positive geometry weight.
        self.assertAlmostEqual(r.term_15_T1[0], 6.102412759438351e-17, delta=1e-23)

    def test_cancelling_fixed_depth_total(self):
        mat = Material("synthetic", 0, (Oscillator(15, 180, 3e-7),))
        with self.assertRaises(FixedDepthConvergenceError):
            fixed_depth(mat, 1000, 14, 10, hole_charge=1.0000077546411514,
                        epsilon_momentum="k_parallel", epsabs=1e-22, epsrel=1e-7)

    def test_term_16h_contributes_to_reported_error(self):
        # One accepted panel. Only the unit-hole term 16h carries uncertainty.
        # The public common bound must also cover its weight at charge -2.
        terms = np.zeros((1, 8))
        errors = terms.copy()
        errors[0, 7] = 1e-15
        mat = Material("synthetic", 0, (Oscillator(15, 1, 1),))
        _, bound = _integrate_loss(mat, 15, 1000, 10, 1, 0, 1, 2, -2, 1e-12, 1e-7, 1,
                                   initial=(np.array([0., 1.]), (terms, errors, errors)))
        self.assertGreaterEqual(bound, 2e-15)

    def test_bulk_quadpack_roundoff_message(self):
        # The panel's halved relative tolerance is below the roundoff floor.
        # The returned error would pass both subsequent (looser) budgets.
        mat = Material("synthetic", 0, (Oscillator(15, 180, 2, 0.01),))
        with self.assertRaisesRegex(BulkConvergenceError, "roundoff"):
            bulk_diimfp(mat, 1000, 1, epsabs=1e-30, epsrel=1.2e-14)

    def test_fixed_depth_panel_and_node_resolution(self):
        # A line inside the panel whose width is below the spacing of its center.
        mat = Material("synthetic", 0, (Oscillator(15, 1, 1e-24, 0.5),))
        with self.assertRaisesRegex(FixedDepthConvergenceError, "linewidth unresolved"):
            _log_breakpoints(mat, 20., 0.5, 3., 0.)
        # The allowed coordinate perturbation straddles the center of this line.
        q = np.ones(1)
        mat = Material("synthetic", 0, (Oscillator(15, 1, 1e-10, 1),))
        loss = 15 + 3.80998212
        coef = _coefficients(mat, loss, q)
        with self.assertRaisesRegex(FixedDepthConvergenceError, "oscillator denominator"):
            _include_node_roundoff(mat, loss, q, np.array([1e-9]), coef)

    def test_optical_moment_rejects_an_impossible_quadrature_value(self):
        # Fault injection at the integration boundary: a retained positive band
        # cannot contribute a negative first moment, even if the estimate is small.
        mat = Material("synthetic", 40, (Oscillator(15, 1, 1),))
        with patch("xpsloss.dielectric._quadrature", return_value=(-1., 0.)):
            with self.assertRaisesRegex(DielectricConvergenceError, "outside range"):
                optical_elf_first_moment(mat)

    def test_gapless_inverse_cancelling_lines(self):
        mat = Material("synthetic", 0, tuple(Oscillator(15 + offset, 100, 1e-9) for offset in (-5e-8, 5e-8)))
        with self.assertRaises(DielectricConvergenceError):
            inverse(mat, 15)

    def test_tabulated_cancellation(self):
        for nodes, values, t, tail in (([0, 0.5, 1], [0, 5e11, 1e12], math.nextafter(1., 2.), "power_law_3"),
                                        ([8.9, 9., 9.1, 10.9, 11., 11.1], [0, 1e9, 0, 0, 904458441.9762187, 0],
                                         10., "hard_cut")):
            with self.subTest(t=t), self.assertRaisesRegex(DielectricConvergenceError, "cancellation"):
                tabulated_elf_kk(nodes, values, t, low_energy="zero", high_energy=tail)

    def test_quadpack_subdivision_failure(self):
        mat = Material("synthetic", 15.01, (Oscillator(15, 1, 1e-6),))
        with self.assertRaisesRegex(DielectricConvergenceError, "subdivisions"):
            inverse(mat, 0, epsabs=1e-6, epsrel=1e-5, limit=5)

    def test_narrow_closed_kk_at_the_gap(self):
        # Independent retained-band references frozen in phase 4a. The closed
        # form actually resolves both narrow lines even at the gap threshold for the center.
        cases = [(w, gap, t, ref) for w, gap, t, ref in frozen_4a.KK if w <= 1e-12]
        self.assertEqual(len(cases), 2)
        for width, gap, t, ref in cases:
            mat = frozen_4a.material(width, 100, gap)
            with self.subTest(width=width):
                self.assertAlmostEqual(closed(mat, t).real, ref, delta=max(1e-9, 1e-7 * abs(ref)))

    def test_unresolved_frozen_points(self):
        for strength, ref in frozen_4a.FIXED:
            mat = frozen_4a.material(1e-24, strength, mixed=True)
            with self.subTest(strength=strength):
                with self.assertRaises(FixedDepthConvergenceError):
                    fixed_depth(mat, 1000, 20, 10, epsilon_momentum="k_total")
                # k_total 13a equals the bulk integral; its independent literal is testable there.
                self.assertAlmostEqual(bulk_diimfp(mat, 1000, 20)[0], ref, delta=max(1e-12, 1e-7*abs(ref)))
        mat = Material("synthetic", 1, (Oscillator(*frozen_4f.LINE),))
        with self.assertRaises(DielectricConvergenceError):
            dielectric_response(mat, math.sqrt(45.))
        for w1, w2, width, alpha, q, gap, broad, _, _ in frozen_followup.CASES:
            lines = tuple(Oscillator(w, 100, width, alpha) for w in (w1, w2))
            if broad:
                lines += (Oscillator(10, 50, 5),)
            with self.subTest(q=q, gap=gap), self.assertRaisesRegex(DielectricConvergenceError, "loss_eV=5"):
                closed(Material("synthetic", gap, lines), [4, 5], q)
