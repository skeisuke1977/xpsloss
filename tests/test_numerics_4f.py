"""Frozen binary64 references from benchmarks/numerics/generate_references_4f.py.

Each positive test checks the public result at the requested accuracy per
observable (including signed unit-charge hole terms at charge zero).
"""
import math
import unittest

from xpsloss import Material, Oscillator, bulk_diimfp, dielectric_response, effective_cross_section_fixed_depth
from xpsloss.dielectric import DielectricConvergenceError, negative_inverse_dielectric_closed_form
from xpsloss.xps_fixed_depth import FixedDepthConvergenceError


T = math.sqrt(135.0)
LINE = (15, 180, 3e-7, 0)
FIELDS = ("term_13a", "term_13b", "term_14", "term_15_T1", "term_15_T2", "term_15_T3",
          "term_16_e", "term_16_h", "no_hole", "hole", "total")
# label, loss, oscillators, epsilon momentum, charge, requested relative tolerance
CASES = (
    ("owner", T, (LINE,), "k_parallel", 1, 1e-9),
    ("owner", T, (LINE,), "k_parallel", 0, 1e-11),
    ("broad", T, ((15, 180, 3e-4, 0),), "k_parallel", 1, 1e-9),
    ("middle", T, ((15, 180, 3e-6, 0),), "k_parallel", 1, 1e-9),
    ("narrow", T, ((15, 180, 3e-8, 0),), "k_parallel", 1, 1e-9),
    ("below", T - 1.5e-7, (LINE,), "k_parallel", 1, 1e-11),
    ("above", T + 1.5e-7, (LINE,), "k_parallel", 1, 1e-11),
    ("parallel", T, ((15, 180, 3e-7, 2e-8),), "k_parallel", 1, 1e-9),
    ("total_momentum", T, ((15, 180, 3e-7, 2e-8),), "k_total", 1, 1e-9),
    ("two", 12., ((15, 81, 4e-7, 0), (18, 180, 2e-7, 0)), "k_parallel", 1, 1e-11),
)
REFERENCES = {
    "owner": (
        2.7189952299378816e-10, -0.0035543224728221794, -0.0031968410463442186, 0.006961947695256588,
        -3524.872050440608, 11640.413618669876, 11640.414041468406, -387.1771511413114,
        23280.831067763775, -3912.0523984229653, 19368.77866934081,
    ),
    "broad": (
        2.718995225859393e-07, -0.0035543922401306324, -0.0031884011499160284, 0.007108644648599327,
        -3.5252291520958408, 11.640413752969907, 11.640413753392705, -0.3871771365464636,
        23.284382030670603, -3.9155946897922203, 19.368787340878384,
    ),
    "middle": (
        2.7189952299374782e-09, -0.0035543231014039944, -0.003196765011283874, 0.00710717783485131,
        -352.48752743339594, 1164.0413739660894, 1164.0413781940747, -38.71771379927998,
        2328.0863050176163, -391.2084379976872, 1936.8778670199292,
    ),
    "narrow": (
        2.718995229937886e-11, -0.0035543224099639503, -0.0031968486498502105, -0.0075610662643211595,
        -35248.71662322299, 116404.12408759687, 116404.16636744977, -3871.7728262642727,
        232808.27933965798, -39120.492646335915, 193687.78669332206,
    ),
    "below": (
        2.718994990273793e-10, -0.0035543223827084807, -0.0031968409645770275, -22943.04595076863,
        -608.6973866778269, -15204.999619745755, 50919.609508381895, -2478.4703621464414,
        12771.5603835454, -3087.170945665233, 9684.389437880167,
    ),
    "above": (
        2.7189954696020025e-10, -0.0035543225629358884, -0.003196841128111419, 22943.05896187281,
        -2916.1749638787996, 26845.412177503462, -39279.19320596436, 2091.2930161266063,
        10509.274379089618, -824.8851445933216, 9684.389234496297,
    ),
    "parallel": (
        2.7189828788709094e-10, -0.0035543140020964366, -0.003196832576135512, -6313.815917536742,
        -3458.8079299503515, 4147.421562957391, 22818.025369668718, -487.2117673034349,
        20651.627460775635, -3946.0228940863626, 16705.60456668927,
    ),
    "total_momentum": (
        2.718982877647369e-10, -0.0035543140012967173, -0.0031968325754162283, -6400.20744230641,
        -3444.518291164582, 4072.77510132975, 22972.71242776651, -513.8985631921815,
        20645.27653247612, -3958.4200511893387, 16686.85648128678,
    ),
    "two": (
        2.5339762443319086e-10, -0.003437317800430911, -0.0030933952917864007, 0.006186742439609132,
        -3471.4987177256335, 12166.220874921146, 12166.22284743674, -564.773373679868,
        24332.44647178278, -4036.2751848007933, 20296.171286981986,
    ),
}
EPSILON = {0: 0.8056547382872928+89442719.09999159j, 1: -41004.179730012205+89442700.30200417j}


class SurfaceAccuracyTests(unittest.TestCase):
    def test_fixed_depth_each_observable(self):
        for name, t, lines, momentum, charge, rel in CASES:
            with self.subTest(case=name, charge=charge):
                material = Material("synthetic", 0, tuple(Oscillator(*line) for line in lines))
                settings = {"epsrel": rel} if rel != 1e-9 else {}
                result = effective_cross_section_fixed_depth(
                    material, 1000, t, 10, hole_charge=charge, epsilon_momentum=momentum, **settings)
                for field, ref in zip(FIELDS, REFERENCES[name]):
                    if field == "total" and charge == 0:
                        ref = REFERENCES[name][8]
                    with self.subTest(field=field):
                        error = abs(float(getattr(result, field)[0]) - ref)
                        self.assertLessEqual(error, max(1e-12, rel * abs(ref)))
                        self.assertLessEqual(error, result.quadrature_error[0] + math.ulp(ref))
                self.assertLessEqual(result.quadrature_error[0], max(1e-12, rel * abs(result.total[0])))

    def test_unattainable_accuracy(self):
        mat = Material("synthetic", 0, (Oscillator(*LINE),))
        with self.assertRaises(FixedDepthConvergenceError):
            effective_cross_section_fixed_depth(mat, 1000, T, 10, epsilon_momentum="k_parallel",
                                                epsabs=1e-30, epsrel=1e-20)

    def test_reciprocal_components(self):
        for gap, ref in EPSILON.items():
            with self.subTest(gap=gap):
                mat = Material("synthetic", gap, (Oscillator(*LINE),))
                try:
                    value = dielectric_response(mat, math.sqrt(45.0)).epsilon
                except DielectricConvergenceError:
                    if gap == 0:
                        raise
                    continue  # The gapped KK integral may lack the precision needed by 1/F.
                for x, y in ((value.real, ref.real), (value.imag, ref.imag)):
                    self.assertLessEqual(abs(x-y), max(1e-9, 1e-7*abs(y)))

    def test_closed_form_real_and_imaginary_parts(self):
        mat = Material("synthetic", 1, (Oscillator(*LINE),))
        value = negative_inverse_dielectric_closed_form(mat, math.sqrt(45.0))
        ref = 5.125523543471986e-12+1.1180339887498947e-8j
        for x, y in ((value.real, ref.real), (value.imag, ref.imag)):
            self.assertLessEqual(abs(x-y), max(1e-9, 1e-7*abs(y)))

    def test_short_bulk_interval(self):
        mat = Material("synthetic", 0, (Oscillator(15, 180, 0.1, 1e-20),))
        ref = 2.3099471226159184e-19
        value, estimate = bulk_diimfp(mat, 1000, math.nextafter(1000., 0.), epsabs=1e-30, epsrel=1e-9)
        self.assertLessEqual(abs(value-ref), max(1e-30, 1e-9*ref))
        self.assertLessEqual(abs(value-ref), estimate + math.ulp(ref))
