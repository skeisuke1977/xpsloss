"""V12(g) and scaled-depth regression controls, with exact binary64 inputs."""

from fractions import Fraction
from itertools import product
import math
import unittest

from xpsloss._depth_kernel import depth_averaged_weights


def rational_weights(k, b, p):
    # Evaluate the complex factors of the Laplace transform with rational
    # Cartesian pairs. The inverse is exact on these binary64 inputs.
    k, b, p = map(Fraction, (k, b, p))
    def inv(a, c):
        return a/(a*a+c*c), -c/(a*a+c*c)
    def mul(x, y):
        return x[0]*y[0]-x[1]*y[1], x[0]*y[1]+x[1]*y[0]
    def scale(x, v):
        return v*x[0], v*x[1]
    s, shifted = inv(k, b), inv(p+k, b)
    prod = mul(s, shifted)
    j = p + 2*k
    D = k*k+b*b
    w13 = scale(prod, -1/p)
    raw = {"W13": w13, "N": (2*shifted[0]/(p*j), 0),
           "W14": scale(mul((0, 1), shifted), 2*k/(p*j)),
           "W15b": scale(mul((0, -1), shifted), 1/j),
           "W15c": (w13[0], -w13[1]), "W16e": scale(shifted, -j/(D*p)),
           "W16h": mul((0, 1), prod)}
    return {name: scale(pair, p*p/(2*b)) for name, pair in raw.items()}


class KernelTests(unittest.TestCase):
    def test_zero_crossing_absolute_error(self):
        p = math.nextafter(1.5, math.inf)
        exact = rational_weights(.5, 1., p)
        weights = depth_averaged_weights(.5, 1., p)
        for key, side in (("W13", 0), ("W15c", 0), ("W16h", 1)):
            self.assertNotEqual(exact[key][side], 0)
            self.assertGreaterEqual(Fraction(weights[key][side+2]), abs(exact[key][side]))

    def test_values_and_eight_corner_error_bounds(self):
        for point in ((.001, .5, 1.), (.05, .03, .05), (5., .03, .05), (.01, 2., .05),
                      (.5, 1., math.nextafter(1.5, math.inf))):
            widths = [4*math.ulp(x) for x in point]
            center = depth_averaged_weights(*point)
            box = depth_averaged_weights(*point, e_k=widths[0], e_b=widths[1], e_p=widths[2])
            for signs in ((0, 0, 0), *product((-1, 1), repeat=3)):
                exact = rational_weights(*(Fraction(x)+sign*Fraction(d) for x, d, sign in zip(point, widths, signs)))
                result = box if any(signs) else center
                for name, (x, y) in exact.items():
                    for j, value in enumerate((x, y)):
                        diff = abs(value-Fraction(result[name][j]))
                        self.assertLessEqual(diff, Fraction(result[name][j+2]), (point, signs, name, j))

    def test_scaled_underflow_and_range_failure(self):
        small = depth_averaged_weights(1., 1., 1e82)["W14"][0]
        self.assertTrue(math.isclose(small, 1.0000000000000002e-164, rel_tol=1e-13))
        with self.assertRaises(FloatingPointError):
            depth_averaged_weights(1., 1., 1e200)
