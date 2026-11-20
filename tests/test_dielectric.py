import math
from fractions import Fraction
from pathlib import Path
import unittest

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

from xpsloss import (
    Material, Oscillator, dielectric_response, loss_function, negative_inverse_dielectric_function,
)
from xpsloss.dielectric import (
    DielectricConvergenceError, DielectricGapEdgeError, dielectric_diagnostics,
    negative_inverse_dielectric_closed_form, optical_elf_first_moment, optical_sum_rule_lhs,
    optical_sum_rule_residual, tabulated_elf_kk,
)
from xpsloss.model import HBAR2_OVER_2ME_EV_A2


SOURCE = Path(__file__).resolve().parents[1] / "benchmarks/sources/zenodo6024064_v1.0"


def retained_kk(material, energy, q):
    """Independent reference using Cauchy-weighted integration and an infinite tail."""
    gap = material.gap_eV
    upper = 100.0
    ungapped = Material(material.name, 0, material.oscillators)

    def elf(u):
        return loss_function(ungapped, u, q)

    def ordinary(u):
        return 2 * u * elf(u) / ((u - energy) * (u + energy))

    if energy < gap:
        integral = quad(ordinary, gap, upper, epsabs=1e-12, epsrel=1e-12)[0]
    else:
        integral = quad(elf, gap, upper, weight="cauchy", wvar=energy,
                        epsabs=1e-12, epsrel=1e-12)[0]
        integral += quad(lambda u: elf(u) / (u + energy), gap, upper,
                         epsabs=1e-12, epsrel=1e-12)[0]
    tail = quad(ordinary, upper, math.inf, epsabs=1e-12, epsrel=1e-12)[0]
    return -1 + (integral + tail) / math.pi


def tabulated_kk_quadrature(nodes, values, energy, low, high):
    """Independently integrate the linear interpolant and its chosen tail."""
    upper = nodes[-1]
    terminal = values[-1]
    if low == "linear_to_origin":
        nodes, values = np.r_[0.0, nodes], np.r_[0.0, values]
    lower = nodes[0]
    cutoff = max(2 * upper, 2 * energy) if high == "power_law_3" else upper

    def elf(u):
        if u <= upper:
            return np.interp(u, nodes, values)
        return terminal * (upper / u)**3

    points = nodes[(lower < nodes) & (nodes < cutoff)].tolist()
    if lower < energy < cutoff:
        # Subtract on the entire finite interval, so a sample knot can be the PV point.
        level = elf(energy)
        points.append(energy)
        cauchy_log = level * math.log((cutoff - energy) / (energy - lower))
    else:
        level, cauchy_log = 0.0, 0.0

    def regular(u):
        return (elf(u) - level) / (u - energy) + elf(u) / (u + energy)

    integral = quad(regular, lower, cutoff, points=sorted(set(points)), epsabs=1e-12, epsrel=1e-12)[0]
    if cutoff > upper:
        integral += quad(lambda u: 2 * u * elf(u) / ((u - energy) * (u + energy)),
                         cutoff, math.inf, epsabs=1e-12, epsrel=1e-12)[0]
    return -1 + (integral + cauchy_log) / math.pi


class DielectricTests(unittest.TestCase):
    def test_gapless_lorentz_reciprocal_and_broadcasting(self):
        mat = Material("synthetic", 0, (Oscillator(10, 60, 2, 0.5),))
        t = np.array([0.0, 7.0, 10.0, 16.0, 1e5])
        q = np.array([[0.0], [0.5], [1.0]])
        response = dielectric_response(mat, t, q)
        center = 10 + 0.5 * HBAR2_OVER_2ME_EV_A2 * q**2
        epsilon = 1 + 60 / (center**2 - 60 - t**2 - 2j * t)
        np.testing.assert_allclose(response.epsilon, epsilon, rtol=2e-14, atol=1e-14)
        np.testing.assert_allclose(response.negative_inverse_epsilon, -1 / epsilon,
                                   rtol=2e-14, atol=1e-14)
        np.testing.assert_array_equal(response.negative_inverse_epsilon.imag, loss_function(mat, t, q))
        np.testing.assert_array_equal(response.real_quadrature_error, np.zeros((3, 5)))

    def test_hard_gap_against_retained_band_principal_value(self):
        mat = Material("synthetic", 8, (Oscillator(6, 10, 0.5, 0.1), Oscillator(11, 30, 0.5, 0.3)))
        energies = [0.0, 3.0, 7.9, 8.1, 18.0]
        for q in (0, 0.8):
            with self.subTest(q=q):
                response = dielectric_response(mat, energies, q)
                reference = [retained_kk(mat, t, q) for t in energies]
                inverse = response.negative_inverse_epsilon
                np.testing.assert_allclose(inverse.real, reference, rtol=1e-9, atol=1e-10)
                closed = negative_inverse_dielectric_closed_form(mat, energies, q)
                np.testing.assert_allclose(closed.real, reference, rtol=1e-9, atol=1e-10)
                np.testing.assert_array_equal(closed.imag, loss_function(mat, energies, q))
                np.testing.assert_allclose(inverse * response.epsilon, -1)
                np.testing.assert_array_equal(inverse.imag, loss_function(mat, energies, q))
                self.assertTrue(np.all(response.epsilon.imag >= 0))

    def test_closed_form_across_critical_damping(self):
        energies = [0, 4.0, 7.99, 8.01, 10.0, 20.0]
        for width in (2.0, 20.0 - 1e-8, 20.0, 20.0 + 1e-8, 50.0):
            for gap in (8, 14):
                with self.subTest(width=width, gap=gap):
                    mat = Material("synthetic", gap, (Oscillator(10, 50, width),))
                    reference = [retained_kk(mat, energy, 0) for energy in energies]
                    closed = negative_inverse_dielectric_closed_form(mat, energies)
                    numerical = negative_inverse_dielectric_function(mat, energies)
                    np.testing.assert_allclose(closed.real, reference, rtol=1e-9, atol=1e-10)
                    np.testing.assert_allclose(numerical.real, reference, rtol=1e-9, atol=1e-10)

    def test_narrow_resonance_below_gap(self):
        mat = Material("synthetic", 10, (Oscillator(7, 20, 1e-9),))
        energies = [0.0, 7.0, 9.0, 10.1]
        reference = [retained_kk(mat, t, 0) for t in energies]
        result = negative_inverse_dielectric_function(mat, energies)
        np.testing.assert_allclose(result.real, reference, rtol=0, atol=1e-13)

    def test_overdamped_oscillator_recovers_static_normalization(self):
        mat = Material("synthetic", 1, (Oscillator(2, 10, 1e9),))
        reference = retained_kk(mat, 0, 0)
        result = dielectric_response(mat, 0, epsabs=1e-12, epsrel=1e-11)
        self.assertAlmostEqual(result.negative_inverse_epsilon.real, reference, delta=1e-11)

    def test_hard_edge_logarithm_and_high_energy_normalization(self):
        oscillator = Oscillator(12, 50, 3)
        mat = Material("synthetic", 8, (oscillator,))
        ungapped = Material("synthetic", 0, (oscillator,))
        prefactor = loss_function(ungapped, 8) / math.pi
        for sign in (-1, 1):
            far, near = negative_inverse_dielectric_function(mat, [8 + sign * 1e-7, 8 + sign * 1e-8])
            self.assertAlmostEqual((near - far).real, prefactor * math.log(10), delta=2e-6)
        with self.assertRaisesRegex(ValueError, "hard gap edge"):
            dielectric_response(mat, 8)
        self.assertEqual(loss_function(mat, 8), 0.0)
        result = dielectric_response(mat, 1e6)
        self.assertAlmostEqual(result.negative_inverse_epsilon.real, -1, places=9)
        self.assertAlmostEqual(result.epsilon.real, 1, places=9)

    def test_hard_gap_excludes_adjacent_floats_in_all_apis(self):
        gap = 3.3
        mat = Material("synthetic", gap, (Oscillator(12, 50, 3),))
        offenders = [math.nextafter(gap, 0), gap, math.nextafter(gap, math.inf)]
        for evaluator in (negative_inverse_dielectric_function, negative_inverse_dielectric_closed_form,
                          dielectric_response):
            with self.subTest(evaluator=evaluator), self.assertRaises(DielectricGapEdgeError) as context:
                evaluator(mat, offenders)
            self.assertEqual(context.exception.offending_loss_eV, offenders)
            self.assertEqual(context.exception.exclusion_eV, 1e-12 * gap)
        results = negative_inverse_dielectric_function(mat, [gap - 1e-9, gap + 1e-9])
        self.assertTrue(np.isfinite(results).all())

    def test_near_gap_warning_reports_grid_dependent_log_term(self):
        gap = 3.3
        mat = Material("synthetic", gap, (Oscillator(15, 180, 2, 0.5),))
        losses = [gap - 1e-4, gap + 1e-9, gap + 1e-4, gap + 0.01]
        result = dielectric_diagnostics(mat, loss_eV=losses)
        self.assertAlmostEqual(result.gap_edge_warning_band_eV, 0.0033)
        self.assertEqual([point.loss_eV for point in result.near_gap_edge_points], losses[:3])
        elf_at_edge = 180 * 2 * gap / ((15**2 - gap**2)**2 + (2 * gap)**2)
        for point in result.near_gap_edge_points:
            distance = abs(point.loss_eV - gap)
            self.assertEqual(point.distance_eV, distance)
            self.assertAlmostEqual(point.edge_log_contribution, elf_at_edge / math.pi * math.log(gap / distance))
        self.assertEqual(len(result.warnings), 3)
        self.assertIn("grid offset", result.warnings[1])

    def test_smallest_positive_gap_keeps_a_nonzero_exclusion(self):
        gap = math.ulp(0.0)
        mat = Material("synthetic", gap, (Oscillator(15, 180, 2, 0.5),))
        result = dielectric_diagnostics(mat)
        self.assertEqual(result.gap_edge_exclusion_eV, gap)
        self.assertEqual(result.sub_gap_pole_status, "unresolved_at_gap_edge")
        self.assertAlmostEqual(dielectric_response(mat, 0).epsilon.real, 5)
        with self.assertRaises(DielectricGapEdgeError):
            dielectric_response(mat, [gap, math.nextafter(gap, math.inf)])

    def test_sub_gap_pole_diagnostic_matches_direct_quadrature_root(self):
        gap = 8
        mat = Material("synthetic", gap, (Oscillator(10, 50, 2, 0.5),))
        q = 0.4
        expected = brentq(lambda energy: retained_kk(mat, energy, q), 0, gap - 0.001)
        diagnostics = dielectric_diagnostics(mat, q)
        self.assertEqual(diagnostics.static_sign, "negative")
        self.assertEqual(diagnostics.sub_gap_pole_status, "located")
        self.assertAlmostEqual(diagnostics.sub_gap_pole_eV, expected, delta=1e-9)
        near_edge = Material("synthetic", 3.3, (Oscillator(100, 1, 1),))
        unresolved = dielectric_diagnostics(near_edge)
        self.assertEqual(unresolved.sub_gap_pole_status, "unresolved_at_gap_edge")
        self.assertIsNone(unresolved.sub_gap_pole_eV)

    def test_narrow_on_resonance_closed_form_reports_ill_conditioning(self):
        for gap, width in ((5, 1e-3), (9.9, 1e-4)):
            with self.subTest(gap=gap, width=width):
                mat = Material("synthetic", gap, (Oscillator(10, 50, width),))
                removed, _ = quad(lambda u: u**2 / ((u**2 - 100) * ((100 - u**2)**2 + (width*u)**2)),
                                  0, gap, epsabs=1e-15, epsrel=1e-12)
                expected = -1 - (2 / math.pi) * 50 * width * removed
                self.assertAlmostEqual(negative_inverse_dielectric_function(mat, 10).real,
                                       expected, delta=1e-10)
                with self.assertRaisesRegex(DielectricConvergenceError, "ill-conditioned"):
                    negative_inverse_dielectric_closed_form(mat, 10)

    def test_vacuum_static_pole_and_positive_static_sign(self):
        vacuum = Material("synthetic", 8, (Oscillator(10, 0, 1),))
        response = dielectric_response(vacuum, 8)
        self.assertEqual(response.negative_inverse_epsilon, -1 + 0j)
        self.assertEqual(response.epsilon, 1 + 0j)
        self.assertEqual(response.real_quadrature_error, 0.0)
        pole = Material("synthetic", 0, (Oscillator(10, 100, 1),))
        self.assertEqual(negative_inverse_dielectric_function(pole, 0), 0j)
        with self.assertRaisesRegex(ValueError, "epsilon has a pole"):
            dielectric_response(pole, 0)
        dielectric_response(pole, 1)
        positive = Material("synthetic", 0, (Oscillator(10, 110, 1),))
        self.assertAlmostEqual(negative_inverse_dielectric_function(positive, 0).real, 0.1)
        self.assertEqual(dielectric_diagnostics(positive).static_sign, "positive")
        with self.assertRaisesRegex(ValueError, "upper half-plane"):
            dielectric_response(positive, 1)

    def test_metal_closed_form_and_rounded_al_input(self):
        energies = np.array([1e-3, 0.1, 1.0, 15.6])
        width = 1.5
        # The second pair is from the Al [4] entry of Zenodo 6024064 v1.0.
        for name, center, strength in (("synthetic_drude", 10.0, 100.0), ("Al_Zenodo_v1", 15.6, 243.3598)):
            with self.subTest(center=center):
                mat = Material(name, 0, (Oscillator(center, strength, width),))
                # omega_transverse**2 = 0 for the exactly normalized Drude metal.
                restoring_force = center**2 - strength
                expected = 1 - strength / (energies**2 + 1j * width * energies - restoring_force)
                result = dielectric_response(mat, energies)
                np.testing.assert_allclose(result.epsilon, expected, rtol=1e-11, atol=0)
                self.assertTrue(np.all(result.epsilon.imag > 0))
                if restoring_force:
                    static = dielectric_response(mat, 0)
                    self.assertAlmostEqual(static.epsilon.real / (center**2 / restoring_force), 1,
                                           places=9)

    def test_marginal_static_metals_are_independent_of_last_bit(self):
        cases = (((10, 30), (20, 280)), ((3, 2.7), (7, 34.3)), ((15.6, 243.36),))
        for pairs in cases:
            with self.subTest(pairs=pairs):
                self.assertEqual(sum(Fraction(str(a)) / Fraction(str(e))**2 for e, a in pairs), 1)
                mat = Material("synthetic_drude", 0, tuple(Oscillator(e, a, 1, 0.5) for e, a in pairs))
                response = dielectric_response(mat, [0.1, 5], q_inv_A=0)
                self.assertTrue(np.all(response.epsilon.imag > 0))
                self.assertEqual(dielectric_diagnostics(mat).static_sign, "marginal")
                with self.assertRaisesRegex(ValueError, "pole at zero energy"):
                    dielectric_response(mat, [0, 5], q_inv_A=0)
                # A finite-q static point is allowed on a broadcast grid with marginal q=0 at T>0.
                dielectric_response(mat, [0, 5], q_inv_A=[0.5, 0])

    def test_zenodo_al_and_au_retained_as_rounded_metal_inputs(self):
        al_lines = (SOURCE / "ELF.txt").read_text(encoding="latin-1").splitlines()
        al_row = al_lines[al_lines.index("Al [4]") + 5]
        au_lines = (SOURCE / "ELF_Au.txt").read_text(encoding="latin-1").splitlines()
        au_rows = au_lines[au_lines.index("@") + 4:]
        for name, lines in (("Al", [al_row]), ("Au", au_rows)):
            with self.subTest(material=name):
                rows = [line.split() for line in lines if line.strip()]
                inverse_static = 1 - sum(Fraction(row[1]) / Fraction(row[0])**2 for row in rows)
                self.assertGreater(inverse_static, 0)
                mat = Material(name, 0, tuple(Oscillator(*map(float, row)) for row in rows))
                energies = np.array([0, 0.1, 1, 15.6])
                expected = 1 - sum(float(a) / (float(e)**2 - energies**2 - 1j * float(w) * energies)
                                   for e, a, w, _ in rows)
                expected[0] = float(inverse_static)
                result = dielectric_response(mat, energies)
                np.testing.assert_allclose(result.epsilon, 1 / expected, rtol=1e-8, atol=0)
                self.assertLess(1e-8, float(inverse_static))
                self.assertEqual(dielectric_diagnostics(mat).static_sign, "negative")

    def test_optical_sum_rules_gapless_closed_forms(self):
        mat = Material("synthetic", 0, (Oscillator(10, 20, 0.01), Oscillator(5, 5, 30, 0.8)))
        lhs = 20 / 10**2 + 5 / 5**2
        def integrate(f):
            near = quad(f, 0, 100, points=[5, 10], epsabs=1e-11, epsrel=1e-11)[0]
            return near + quad(f, 100, math.inf, epsabs=1e-11, epsrel=1e-11)[0]

        static = 2 / math.pi * integrate(lambda u: loss_function(mat, u) / u)
        first = integrate(lambda u: u * loss_function(mat, u))
        self.assertAlmostEqual(static, lhs, places=10)
        self.assertAlmostEqual(first, 25 * math.pi / 2, places=9)
        self.assertAlmostEqual(optical_sum_rule_lhs(mat), static, places=10)
        self.assertAlmostEqual(optical_elf_first_moment(mat), first, places=9)
        self.assertAlmostEqual(optical_sum_rule_residual(mat, 0.5), lhs - 0.5)
        with self.assertRaises(ValueError):
            optical_sum_rule_residual(mat, math.inf)

    def test_optical_sum_rules_gapped_integrals(self):
        gap = 8.0
        mat = Material("synthetic", gap, (Oscillator(10, 35, 1, 0.4), Oscillator(7, 20, 2)))
        static, _ = quad(lambda u: loss_function(mat, u) / u, gap, math.inf, epsabs=1e-12)
        first, _ = quad(lambda u: u * loss_function(mat, u), gap, math.inf, epsabs=1e-12)
        self.assertAlmostEqual(optical_sum_rule_lhs(mat), 2 * static / math.pi, places=10)
        self.assertAlmostEqual(1 + negative_inverse_dielectric_function(mat, 0).real,
                               2 * static / math.pi, places=10)
        self.assertAlmostEqual(optical_elf_first_moment(mat), first, places=9)

    def test_tabulated_linear_hard_cut_and_low_energy_choices(self):
        bounds = [0.25, 0.5, 1.0]
        energies = np.array([0, 0.1, 0.5, 0.75, 1.25, 3])
        settings = dict(high_energy="hard_cut")
        full = tabulated_elf_kk(bounds, bounds, energies, low_energy="linear_to_origin", **settings)
        expected = -1 + (2 + energies * np.log(np.abs((1 - energies) / (1 + energies)))) / math.pi
        np.testing.assert_allclose(full.real_negative_inverse_epsilon, expected, atol=1e-13, rtol=0)
        self.assertEqual(full.sample_bounds_eV, (0.25, 1.0))
        self.assertEqual(full.low_energy, "linear_to_origin")
        self.assertEqual(full.interpolation, "piecewise_linear")
        self.assertIsNone(full.tail_coefficient_eV3)
        cut = tabulated_elf_kk(bounds, bounds, 0, low_energy="zero", **settings)
        self.assertAlmostEqual(cut.real_negative_inverse_epsilon, -1 + 2 * 0.75 / math.pi)
        with self.assertRaisesRegex(ValueError, "lower cutoff"):
            tabulated_elf_kk(bounds, bounds, 0.25, low_energy="zero", **settings)
        with self.assertRaisesRegex(ValueError, "upper cutoff"):
            tabulated_elf_kk(bounds, bounds, 1, low_energy="zero", **settings)

    def test_tabulated_power_tail_through_continuous_join(self):
        queries = [0, 0.1, 0.9, np.nextafter(1.0, 0.0), 1, np.nextafter(1.0, 2.0), 1.1, 3]
        result = tabulated_elf_kk([0, 0.5, 1], [0, 0.5, 1], queries,
                                  low_energy="zero", high_energy="power_law_3")
        expected = []
        for t in queries:
            if t == 0:
                integral = 8 / 3
            elif abs(t - 1) < 1e-14:
                integral = 0
            else:
                integral = 2 * (1 - t**-2) + (t**-3 - t) * math.log((1 + t) / abs(1 - t))
            expected.append(-1 + integral / math.pi)
        np.testing.assert_allclose(result.real_negative_inverse_epsilon, expected, atol=1e-13, rtol=0)
        self.assertEqual(result.tail_coefficient_eV3, 1)
        self.assertEqual(result.terminal_log_slope, 1)
        self.assertEqual(result.high_energy, "power_law_3")
        self.assertIn("differs from -3", result.warnings[0])
        np.testing.assert_array_equal(result.real_quadrature_error, np.zeros(len(queries)))
        fitted = tabulated_elf_kk([2, 4], [1, 0.125], 0, low_energy="zero",
                                  high_energy="power_law_3")
        self.assertAlmostEqual(fitted.terminal_log_slope, -3)
        self.assertEqual(fitted.warnings, ())

    def test_tabulated_exact_segments_agree_with_quadrature(self):
        nodes = np.array([0.3, 0.8, 1.2, 2.5, 5.0, 7.0])
        values = np.array([0.1, 0.7, 1.5, 0.5, 0.09, 0.05])
        for low in ("zero", "linear_to_origin"):
            for high in ("hard_cut", "power_law_3"):
                with self.subTest(low=low, high=high):
                    queries = [0, 0.15, 0.31, 0.8, 0.95, 2.5, 6.99, 7.01, 9, 1000]
                    if high == "power_law_3":
                        queries.append(7)
                    expected = [tabulated_kk_quadrature(nodes, values, t, low, high) for t in queries]
                    actual = tabulated_elf_kk(nodes, values, queries, low_energy=low, high_energy=high)
                    np.testing.assert_allclose(actual.real_negative_inverse_epsilon, expected,
                                               rtol=0, atol=2e-12)
                    self.assertEqual(actual.method, "analytic_piecewise_linear_kk")

    def test_tabulated_invalid_inputs(self):
        for energies, elf in (([0, 1], [1, 1]), ([0, 0], [0, 1]), ([0, 1], [0, -1]),
                              ([0, math.nan], [0, 1]), ([0, 1], [0]), ([[0, 1]], [[0, 1]]),
                              ([0, 1], [0, 1j])):
            with self.subTest(energies=energies, elf=elf), self.assertRaises(ValueError):
                tabulated_elf_kk(energies, elf, 0, low_energy="zero", high_energy="hard_cut")
        with self.assertRaises(ValueError):
            tabulated_elf_kk([0, 1], [0, 1], 0, low_energy="zero", high_energy="unknown")

    def test_invalid_coordinates_and_integration_settings(self):
        mat = Material("synthetic", 0, (Oscillator(10, 60, 2),))
        for energy, q in ((-1, 0), (0, -1), (math.nan, 0), (0, math.inf), (1j, 0)):
            with self.subTest(energy=energy, q=q), self.assertRaises(ValueError):
                dielectric_response(mat, energy, q)
        for settings in ({"epsabs": 0}, {"epsrel": math.nan}, {"limit": True}):
            with self.subTest(settings=settings), self.assertRaises(ValueError):
                dielectric_response(mat, 1, **settings)
