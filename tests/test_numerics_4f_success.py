"""The default-setting cases of the frozen dispersive file must succeed, not be rejected.

test_numerics_4f_fu2 accepts either a resolved value or FixedDepthConvergenceError.
Its first two cases use the default accuracy (epsrel 1e-9) and are resolved by the
phase 4f implementation; this guards them against future over-rejection.
"""
import math
import unittest

from test_numerics_4f_fu2 import CASES, FIELDS
from xpsloss import Material, Oscillator, effective_cross_section_fixed_depth


class DefaultSettingSuccessTests(unittest.TestCase):
    def test_default_accuracy_cases_succeed(self):
        for loss, lines, momentum, charge, epsabs, epsrel, refs in CASES[:2]:
            self.assertEqual(epsrel, 1e-9)
            with self.subTest(loss=loss, momentum=momentum):
                mat = Material("synthetic", 0, tuple(Oscillator(*line) for line in lines))
                result = effective_cross_section_fixed_depth(mat, 1000, loss, 10, epsilon_momentum=momentum,
                                                             hole_charge=charge, epsabs=epsabs, epsrel=epsrel)
                for field, ref in zip(FIELDS, refs):
                    delta = abs(getattr(result, field)[0] - ref)
                    self.assertLessEqual(delta, max(epsabs, epsrel * abs(ref)), field)
                    self.assertLessEqual(delta, result.quadrature_error[0] + math.ulp(ref), field)


if __name__ == "__main__":
    unittest.main()
