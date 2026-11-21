"""Independent references for dispersive integration nodes and the stopping edge.

The review's ref_fixed.py integrates Eqs. (13)-(16) at high precision, finding
surface lines by bisection. Two panel layouts and Gauss orders (40, 48) agree.
The bulk reference is a separate 60-digit direct integral in q, using binary64
model inputs. No reference calculation or third-party data is needed at runtime.
"""
import math
import unittest

from xpsloss import Material, Oscillator, bulk_diimfp, effective_cross_section_fixed_depth
from xpsloss.xps_fixed_depth import FixedDepthConvergenceError


FIELDS = ("term_13a", "term_13b", "term_14", "term_15_T1", "term_15_T2", "term_15_T3",
          "term_16_e", "term_16_h", "no_hole", "hole", "total")
LINE = (15, 180, 1e-6, 0.05)
# Loss, oscillator lines, momentum, charge, epsabs, epsrel, reference at the first loss.
CASES = (
    (52.24548795323831, (LINE,), "k_parallel", 1, 1e-12, 1e-9, (
        4.3660759431107286e-05, 1.7992673361208e-05, 0.00018895023858745753,
        -2.8221074410045335e-05, 7.040707046264865e-08, 6.571186190965159e-06,
        -1.5184810626087405e-05, -7.94147301609275e-08, 2.4818733947147705e-05,
        0.00018894123092775926, 0.00021375996487490697)),
    ([30., 30., 15.], (LINE,), "k_total", 1, 1e-12, 1e-9, (
        0.00018897263643722077, 0.0001250432769615234, 0.0015384928637470487,
        -0.00020089457467429824, 7.805514105419669e-06, 4.1722542565877424e-05,
        -9.494757136432818e-05, -5.082288385865093e-07, 5.9896309925995184e-05,
        0.0015457901490138818, 0.001605686458939877)),
    (30., (LINE,), "k_parallel", 0, 1e-12, 1e-11, (
        0.00018883101317496745, 0.00012499518443012228, 0.0015379134960756685,
        -0.00020080966656094905, 7.801476542203376e-06, 4.1708282088569945e-05,
        -9.491547276954845e-05, -5.080045519998742e-07, 5.980934036316216e-05,
        0.001545206968065872, 5.980934036316216e-05)),
    (11.717158173273155, ((15, 180, 5e-7, 1e-4),), "k_total", -1, 1e-12, 1e-12, (
        4.702409742485946e-10, -0.0035738063436218078, -0.00321064709496215,
        0.06739656981701202, -0.0038384808474837447, 0.05463222957083817,
        -0.11768608618758047, 0.006803719518402016, 0.0007689073268888947,
        -0.0002454084240438791, 0.0010143157509327739)),
    (15.5, ((15, 120, 4e-7, 0.05), (25, 200, 6, 0)), "k_parallel", -1, 1e-12, 1e-12, (
        0.007667383346231233, 0.0031053143226844808, 0.04777954155249117,
        -0.007639511943886983, 0.0008975729925511861, 0.0004207879825845607,
        -0.0009122537219707719, 6.044862538115752e-05, 0.0026417199856425193,
        0.04873756317042351, -0.046095843184780994)),
    (30., ((15, 180, 1e-3, 0.05),), "k_total", 0.5, 1e-20, 1e-11, (
        0.0001889965866605501, 0.00012503650781143754, 0.0015385069755897384,
        -0.00020089143419568675, 7.805525575304148e-06, 4.172273940399138e-05,
        -9.4945145404992e-05, -5.082958282319145e-07, 5.991925427530029e-05,
        0.0015458042053368107, 0.0008328213569437056)),
)


class CoordinateAccuracyTests(unittest.TestCase):
    def test_dispersive_results_resolved_or_rejected(self):
        for loss, lines, momentum, charge, epsabs, epsrel, refs in CASES:
            with self.subTest(loss=loss, lines=lines, momentum=momentum):
                mat = Material("synthetic", 0, tuple(Oscillator(*line) for line in lines))
                try:
                    result = effective_cross_section_fixed_depth(mat, 1000, loss, 10, epsilon_momentum=momentum,
                                                                 hole_charge=charge, epsabs=epsabs, epsrel=epsrel)
                except FixedDepthConvergenceError:
                    continue
                for field, ref in zip(FIELDS, refs):
                    with self.subTest(field=field):
                        delta = abs(getattr(result, field)[0] - ref)
                        self.assertLessEqual(delta, max(epsabs, epsrel * abs(ref)))
                        self.assertLessEqual(delta, result.quadrature_error[0] + math.ulp(ref))

    def test_bulk_above_former_stop_threshold(self):
        mat = Material("synthetic", 0, (Oscillator(15, 180, 0.1, 1e-20),))
        ref = 7.185536774334315e-14
        value, error = bulk_diimfp(mat, 1000, 1000 * (1 - 1.1e-5), epsabs=1e-60, epsrel=1e-13)
        self.assertLessEqual(abs(value - ref), ref * 1e-13)
        self.assertLessEqual(abs(value - ref), error + math.ulp(ref))
