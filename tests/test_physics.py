import math
import unittest

import numpy as np
from scipy.integrate import quad

from xpsloss.model import BOHR_RADIUS_A, HBAR2_OVER_2ME_EV_A2, Material, Oscillator, bulk_diimfp, loss_function, momentum_limits


class PhysicsTests(unittest.TestCase):
    def test_energy_weighted_elf_sum_rule(self):
        # Independent moment sum rule; do not copy the implemented rational ELF.
        # For gap=0, a damped oscillator has integral T*L(T) dT = pi*A/2.
        for width in [0.4, 2.0, 12.0]:
            mat = Material("synthetic", 0, (Oscillator(15, 180, width, 0),))
            moment, err = quad(lambda t: t*loss_function(mat, t, 0), 0, np.inf,
                               epsabs=1e-8, epsrel=1e-10, limit=400)
            expected = math.pi * 180 / 2
            self.assertLess(abs(moment-expected)/expected, 1e-9)

    def test_optical_no_dispersion_matches_analytic_bulk(self):
        # Independent closed-form integral: the q-independent ELF factors out.
        mat = Material("synthetic", 0, (Oscillator(15, 180, 2, 0),))
        for energy in [200, 1000, 8000]:
            for loss in [0.01, 10.0, 15.0, 60.0]:
                k, error = bulk_diimfp(mat, energy, loss)
                log_ratio = 2 * math.log(math.sqrt(energy) + math.sqrt(energy-loss)) - math.log(loss)
                expected = 180 * 2 * loss / ((15**2-loss**2)**2 + (2*loss)**2)
                expected *= log_ratio / (math.pi * BOHR_RADIUS_A * energy)
                self.assertAlmostEqual(k / expected, 1, places=12)
                self.assertGreaterEqual(error, 0)

    def test_zero_strength_and_gap(self):
        zero = Material("zero", 0, (Oscillator(15, 0, 2, 0.5),))
        gap = Material("gapped", 3, (Oscillator(15, 180, 2, 0.5),))
        self.assertEqual(bulk_diimfp(zero, 1000, 15), (0, 0))
        np.testing.assert_array_equal(loss_function(gap, [0, 0.1, 2.9, 3], 0), np.zeros(4))
        self.assertEqual(bulk_diimfp(gap, 1000, 2), (0, 0))
        self.assertGreater(loss_function(gap, 3.1, 0), 0)

    def test_dispersion_convergence(self):
        mat = Material("synthetic", 0, (Oscillator(15, 180, 0.2, 0.5), Oscillator(30, 20, 3, 0.2)))
        for loss in [14.9, 15.0, 20.0, 60.0]:
            loose, _ = bulk_diimfp(mat, 1000, loss, epsabs=1e-8, epsrel=1e-6)
            tight, _ = bulk_diimfp(mat, 1000, loss, epsabs=1e-12, epsrel=1e-11, limit=400)
            self.assertLess(abs(loose-tight), max(1e-12, abs(tight)*1e-6))

    def test_units_and_stable_q_limits(self):
        # At E=H, the initial wavevector is exactly 1 inverse angstrom.
        e = HBAR2_OVER_2ME_EV_A2
        qmin, qmax = momentum_limits(e, e*0.75)
        self.assertAlmostEqual(qmin, 0.5)
        self.assertAlmostEqual(qmax, 1.5)
        qmin, qmax = momentum_limits(1000, 1e-12)
        self.assertGreater(qmin, 0)
        self.assertAlmostEqual(qmin*qmax, 1e-12/HBAR2_OVER_2ME_EV_A2, places=24)

    def test_invalid_physics_inputs(self):
        for args in [(0, 1, 1, 0), (1, -1, 1, 0), (1, 1, 0, 0), (1, 1, 1, -1), (math.nan, 1, 1, 0)]:
            with self.assertRaises(ValueError):
                Oscillator(*args)
        mat = Material("synthetic", 0, (Oscillator(15, 180, 2, 0.5),))
        for loss in [0, -1, 1000, 1001]:
            with self.assertRaises(ValueError):
                bulk_diimfp(mat, 1000, loss)


if __name__ == "__main__":
    unittest.main()
