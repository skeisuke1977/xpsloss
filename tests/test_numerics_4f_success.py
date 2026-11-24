"""All six resolved fu2 cases must keep their component accuracy and error bounds."""
import math
import unittest

import numpy as np

import test_numerics_4f as frozen
from test_numerics_4f_fu2 import CASES, FIELDS
from xpsloss import Material, Oscillator, effective_cross_section_fixed_depth


# Third point of the repeated-loss batch. Independently rounded by the shipped
# generate_references_4f.fixed(15., [LINE], True, rule, resolve_poles=True)
# at (50 digits, order 32) and (80 digits, order 48).
BATCH_TAIL = (5.6157776338929216e-05, -0.2443070844169099, -0.11146820897126988,
              0.48938435704832256, -0.13282647803664147, 0.0007059860565844847,
              -0.0015349907587858716, 8.573729759505383e-05, 0.2443044257055502,
              -0.24420894971031631, 9.547599523388455e-05)


class CoordinateSuccessTests(unittest.TestCase):
    def test_frozen_reference_dimensions(self):
        self.assertEqual(len(frozen.CASES), 10)
        self.assertEqual(set(frozen.REFERENCES), {case[0] for case in frozen.CASES})
        self.assertEqual(len(frozen.FIELDS), 11)
        for row in frozen.REFERENCES.values():
            self.assertEqual(len(row), 11)

    def test_all_six_dispersive_cases_succeed(self):
        self.assertEqual(len(CASES), 6)
        self.assertEqual(len(FIELDS), 11)
        for loss, lines, momentum, charge, epsabs, epsrel, refs in CASES:
            with self.subTest(loss=loss, momentum=momentum, charge=charge):
                mat = Material("synthetic", 0, tuple(Oscillator(*line) for line in lines))
                result = effective_cross_section_fixed_depth(mat, 1000, loss, 10, epsilon_momentum=momentum,
                                                             hole_charge=charge, epsabs=epsabs, epsrel=epsrel)
                expected = (refs,)
                if isinstance(loss, list):
                    self.assertEqual(loss, [30., 30., 15.])
                    expected = (refs, refs, BATCH_TAIL)
                np.testing.assert_array_equal(result.loss_eV, np.atleast_1d(loss))
                for i, (row, error) in enumerate(zip(expected, result.quadrature_error, strict=True)):
                    for field, ref in zip(FIELDS, row, strict=True):
                        delta = abs(getattr(result, field)[i] - ref)
                        self.assertLessEqual(delta, max(epsabs, epsrel * abs(ref)), (i, field))
                        self.assertLessEqual(delta, error + math.ulp(ref), (i, field))
