import math
import unittest

import numpy as np
from scipy.integrate import quad

from xpsloss import (
    Material, Oscillator, dielectric_response, loss_function, negative_inverse_dielectric_function,
)
from xpsloss.dielectric import (
    negative_inverse_dielectric_closed_form, optical_elf_first_moment, optical_sum_rule_lhs,
    optical_sum_rule_residual, tabulated_elf_kk,
)
from xpsloss.model import HBAR2_OVER_2ME_EV_A2


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

    def test_vacuum_static_pole_and_unstable_epsilon(self):
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
        unstable = Material("synthetic", 0, (Oscillator(10, 110, 1),))
        self.assertAlmostEqual(negative_inverse_dielectric_function(unstable, 0).real, 0.1)
        with self.assertRaisesRegex(ValueError, "upper half-plane"):
            dielectric_response(unstable, 1)

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

    def test_optical_sum_rules_gapless_closed_forms(self):
        mat = Material("synthetic", 0, (Oscillator(10, 20, 0.01), Oscillator(5, 5, 30, 0.8)))
        lhs = 20 / 10**2 + 5 / 5**2
        self.assertAlmostEqual(optical_sum_rule_lhs(mat), lhs)
        self.assertAlmostEqual(optical_elf_first_moment(mat), 25 * math.pi / 2)
        self.assertAlmostEqual(optical_sum_rule_residual(mat, 0.5), lhs - 0.5)
        with self.assertRaises(ValueError):
            optical_sum_rule_residual(mat, math.inf)

    def test_optical_sum_rules_gapped_integrals(self):
        gap = 8.0
        mat = Material("synthetic", gap, (Oscillator(10, 35, 1, 0.4), Oscillator(7, 20, 2)))
        static, _ = quad(lambda u: loss_function(mat, u) / u, gap, math.inf, epsabs=1e-12)
        first, _ = quad(lambda u: u * loss_function(mat, u), gap, math.inf, epsabs=1e-12)
        self.assertAlmostEqual(optical_sum_rule_lhs(mat), 2 * static / math.pi, places=10)
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
