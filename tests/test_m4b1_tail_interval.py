"""The all-axis remainder interval is rounded outward once (M4b-1 prereg 4.2 and 7, row DS0).

GitHub CI run 46062347681 found DS0's ``neglected_total_area`` interval one or two ulps wider than
the rebuild from the reported values on Linux and macOS: the w interval was rounded outward
before it was combined with the ratio interval. This test rebuilds the interval with exact
rationals and a single outward rounding, independently of the package's rounding helpers.
"""
from fractions import Fraction
import math
import random
from types import SimpleNamespace
import unittest

from xpsloss._order_bounds import tail_area


def down(x):
    f = float(x)
    return math.nextafter(f, -math.inf) if Fraction(f) > x else f


def up(x):
    f = float(x)
    return math.nextafter(f, math.inf) if Fraction(f) < x else f


class TailIntervalRoundsOnce(unittest.TestCase):
    def test_single_outward_rounding(self):
        rng = random.Random(20261128)
        for case in range(400):
            m = rng.uniform(1.0, 100.0)
            m_interval = (m * (1 - rng.uniform(0, 1e-7)), m * (1 + rng.uniform(0, 1e-7))) if case % 2 else None
            area = rng.uniform(0.01, 2.0)
            area_error = area * rng.uniform(0, 1e-6)
            r = rng.uniform(0.1, 0.9)
            r_interval = (r * (1 - rng.uniform(0, 1e-7)), r * (1 + rng.uniform(0, 1e-7)))
            mass = rng.choice([1.0, rng.uniform(0.5, 3.0)])
            n = rng.choice([0, 1, 3])
            single = SimpleNamespace(primary_integral=mass, single_loss_weight=float(Fraction(m) * Fraction(area)),
                                     multiplier_A=m, multiplier_interval_A=m_interval,
                                     window_area={"total": area}, window_area_error={"total": area_error})
            ratio = {"value": r, "interval": r_interval}
            got = tail_area(single, ratio, n, False, False, False)["interval"]
            lo_m, hi_m = (Fraction(m), Fraction(m)) if m_interval is None else map(Fraction, m_interval)
            corners = [x * y for x in (lo_m, hi_m)
                       for y in (Fraction(area) - Fraction(area_error), Fraction(area) + Fraction(area_error))]
            w_lo, w_hi = min(corners), max(corners)
            r_lo, r_hi = map(Fraction, r_interval)
            S = Fraction(mass)
            want = (down(S * w_lo * r_lo**n / (1 - r_lo)), up(S * w_hi * r_hi**n / (1 - r_hi)))
            self.assertEqual([x.hex() for x in got], [x.hex() for x in want], f"case {case}")


if __name__ == "__main__":
    unittest.main()
