"""Finite-width regression data from benchmarks/numerics/generate_references_4a.py.

Bulk and k_total 13a use an 80-digit partial-fraction integral. Hard-gap KK
and optical moments use real-axis quadrature with resolved flanks and PV
subtraction. Triangles integrate the actual binary64 piecewise-linear ELF.
The adjacent-point references evaluate the complex oscillator at 80 digits.
Only standard package dependencies are needed to run these tests.
"""

import contextlib
import io
import json
import math
from pathlib import Path
import tempfile
import unittest

from xpsloss import (Material, Oscillator, bulk_diimfp, dielectric_response,
                     effective_cross_section_fixed_depth, loss_function)
from xpsloss.cli import main
from xpsloss.config import load_config
from xpsloss.dielectric import (DielectricConvergenceError, negative_inverse_dielectric_closed_form,
                                negative_inverse_dielectric_function, optical_elf_first_moment,
                                optical_sum_rule_lhs, optical_sum_rule_residual, tabulated_elf_kk)
from xpsloss.xps_fixed_depth import FixedDepthConvergenceError


ROOT = Path(__file__).resolve().parents[1]
BROAD = Oscillator(15, 10, 20, 0.5)
K_ABS = 1e-12  # Å^-1 eV^-1, including term 13a
F_ABS = 1e-9   # dimensionless
MOMENT_ABS = 1e-8  # eV²
REL = 1e-7

# Width (eV), single K, mixed K. The 1e-5 mixed row is the owner's example.
WIDTHS = (
    (1, 0.000973872813344618, 0.0010121170447319936),
    (0.1, 0.0008639681199304697, 0.0009022123513178453),
    (0.01, 0.0008517472671590883, 0.0008899914985464639),
    (0.001, 0.0008505139197337995, 0.000888758151121175),
    (0.0001, 0.0008503904735710441, 0.0008886347049584197),
    (1e-05, 0.0008503781278417825, 0.000888622359229158),
    (1e-06, 0.0008503768932577276, 0.0008886211246451032),
    (1e-07, 0.0008503767697992108, 0.0008886210011865863),
    (1e-09, 0.0008503767562187727, 0.0008886209876061484),
    (1e-12, 0.0008503767560817338, 0.0008886209874691093),
)
# A_narrow (eV²), mixed K. A_broad = 10 eV².
STRENGTHS = (
    (0.001, 3.824895571030801e-05),
    (0.01, 3.82914746167001e-05),
    (0.1, 3.871666368062099e-05),
    (1, 4.2968554319829905e-05),
    (10, 8.548746071191902e-05),
    (50, 0.00027446037801009287),
    (100, 0.0005106765246328102),
    (180, 0.000888622359229158),
    (1000, 0.004762567163841722),
)
# Binary64 optical centers: shifts of -5, -0.5, 0, +0.5, +5 and +100 widths
# relative to the lower (alpha=0.5) and upper (alpha=0.001) momentum limit.
EDGES = (
    (19.949493161166536, 0.5, 1e-07, 0.08155146483107366),
    (19.949493611166535, 0.5, 1e-07, 0.06317668667430791),
    (19.949493661166535, 0.5, 1e-07, 0.04213045169364175),
    (19.949493711166536, 0.5, 1e-07, 0.021084194010287823),
    (19.949494161166534, 0.5, 1e-07, 0.002708781830272134),
    (19.949503661166535, 0.5, 1e-07, 0.00017202760401988426),
    (19.949493656166535, 0.5, 1e-09, 0.08155256542567167),
    (19.949493660666537, 0.5, 1e-09, 0.06317702749791938),
    (19.949493661166535, 0.5, 1e-09, 0.042130758335852656),
    (19.949493661666533, 0.5, 1e-09, 0.021084545625680472),
    (19.949493666166536, 0.5, 1e-09, 0.002709056914568307),
    (19.949493761166536, 0.5, 1e-09, 0.00017222471431963476),
    (16.040100512677668, 0.001, 1e-07, 7.230917276212717e-05),
    (16.040100962677666, 0.001, 1e-07, 0.0003066782235611053),
    (16.040101012677667, 0.001, 1e-07, 0.0005751120841746262),
    (16.040101062677667, 0.001, 1e-07, 0.0008435459483473125),
    (16.040101512677666, 0.001, 1e-07, 0.0010779150994213776),
    (16.040111012677666, 0.001, 1e-07, 0.0011102735828476486),
    (16.040101007677666, 0.001, 1e-09, 7.230904179259365e-05),
    (16.04010101217767, 0.001, 1e-09, 0.00030667868279943766),
    (16.040101012677667, 0.001, 1e-09, 0.0005751119466462401),
    (16.040101013177665, 0.001, 1e-09, 0.000843545204197913),
    (16.040101017677667, 0.001, 1e-09, 0.0010779148400021042),
    (16.040101112677668, 0.001, 1e-09, 0.0011102707856247077),
)
KK = (
    (1e-05, 5, 20, -1.5714285698589214),
    (1e-05, 16, 20, -1.0000012256306672),
    (1e-05, 14.99995, 14.9999500001, 91493.11585868032),
    (1e-12, 15, 20, -1.2857142857149322),
    (1e-24, 15, 20, -1.2857142857142858),
)
OPTICAL = (
    (1e-05, 5, 0.44444437638992373, 157.07963173194267),
    (1e-10, 15.00000000005, 0.1111123617724055, 39.270350198984595),
    (1e-24, 15, 0.2222222222222222, 78.53981633974483),
)
TRIANGLES = (
    (0.0001, -1.054567409064183),
    (1e-06, -1.0545674090192405),
    (1e-08, -1.054567413575006),
    (1e-12, -1.0545722601345207),
)
FIXED = (
    (1e-08, 3.824423143461871e-05),
    (180, 0.0008886209874689722),
)


def material(width, strength=180, gap=0, center=15, alpha=0.5, mixed=False):
    lines = (Oscillator(center, strength, width, alpha),)
    return Material("synthetic", gap, lines + ((BROAD,) if mixed else ()))


class FiniteWidthTests(unittest.TestCase):
    def close(self, value, reference, absolute):
        self.assertAlmostEqual(value, reference, delta=max(absolute, REL * abs(reference)))

    def check_bulk(self, mat, ref):
        value, error = bulk_diimfp(mat, 1000, 20)
        self.close(value, ref, K_ABS)
        self.assertGreaterEqual(error, 0)

    def test_bulk_width_series(self):
        for width, single, mixed in WIDTHS:
            for background, ref in ((False, single), (True, mixed)):
                with self.subTest(width=width, mixed=background):
                    self.check_bulk(material(width, mixed=background), ref)

    def test_bulk_strength_series(self):
        for strength, ref in STRENGTHS:
            with self.subTest(strength=strength):
                self.check_bulk(material(1e-5, strength, mixed=True), ref)

    def test_bulk_at_both_limits_and_on_either_side(self):
        for center, alpha, width, ref in EDGES:
            with self.subTest(center=center, alpha=alpha, width=width):
                self.check_bulk(material(width, center=center, alpha=alpha, mixed=True), ref)

    def test_shipped_inputs(self):
        for name, ref in (("al", 0.001516988937396985), ("au", 0.0007527858133448525),
                          ("si", 0.0021632071306761785)):
            with self.subTest(name=name):
                cfg = load_config(ROOT / f"benchmarks/inputs/{name}_elfdb_v1.0.bulk.json")
                self.check_bulk(cfg.material, ref)

    def test_adjacent_to_gapless_resonance(self):
        t = math.nextafter(15., math.inf)
        for width, ref in ((1e-14, -210303540032860.28+591951840329175.8j),
                           (1e-16, -1875014305356693.8+52776960793242.195j)):
            with self.subTest(width=width):
                mat = material(width, 100)
                f = dielectric_response(mat, t).negative_inverse_epsilon
                self.close(f.real, ref.real, F_ABS)
                self.close(f.imag, ref.imag, F_ABS)
                self.close(loss_function(mat, t), ref.imag, F_ABS)

    def test_gapped_kk_and_closed_form(self):
        for width, gap, t, ref in KK:
            for evaluator in (negative_inverse_dielectric_function, negative_inverse_dielectric_closed_form):
                with self.subTest(width=width, gap=gap, t=t, method=evaluator.__name__):
                    try:
                        value = evaluator(material(width, 100, gap), t).real
                    except DielectricConvergenceError as exc:
                        # Only the exceptionally narrow lines at the integration edge
                        # may be explicitly unresolved. The near-gap 1e-5 case must pass.
                        self.assertLessEqual(width, 1e-12)
                        self.assertIn("oscillator", str(exc))
                        self.assertIn("loss_eV", str(exc))
                    else:
                        self.close(value, ref, F_ABS)

    def test_optical_moments_and_static_response(self):
        for width, gap, lhs, first in OPTICAL:
            with self.subTest(width=width, gap=gap):
                mat = material(width, 100, gap)
                l = optical_sum_rule_lhs(mat)
                m = optical_elf_first_moment(mat)
                self.close(l, lhs, F_ABS)
                self.close(m, first, MOMENT_ABS)
                self.close(optical_sum_rule_residual(mat, 0.5), lhs - 0.5, F_ABS)
                self.assertTrue(0 <= l <= 100 / 225)
                self.assertTrue(0 <= m <= 50 * math.pi)
                try:
                    static = negative_inverse_dielectric_function(mat, 0).real
                except DielectricConvergenceError:
                    self.assertLessEqual(width, 1e-10)
                else:
                    self.close(1 + static, lhs, F_ABS)

    def test_tabulated_triangles(self):
        for h, ref in TRIANGLES:
            with self.subTest(h=h):
                r = tabulated_elf_kk([14, 15-h, 15, 15+h, 16], [0, 0, 1/h, 0, 0], 20,
                                     low_energy="zero", high_energy="hard_cut")
                self.close(r.real_negative_inverse_epsilon, ref, F_ABS)

    def test_fixed_depth_13a(self):
        cases = [(1e-5, 180, 0.000888622359229158)]
        cases.extend((1e-24, a, ref) for a, ref in FIXED)
        for width, strength, ref in cases:
            with self.subTest(width=width, strength=strength):
                try:
                    result = effective_cross_section_fixed_depth(material(width, strength, mixed=True),
                        1000, 20, 10, epsilon_momentum="k_total")
                except FixedDepthConvergenceError as exc:
                    self.assertEqual(width, 1e-24)
                    self.assertIn("loss_eV", str(exc))
                    self.assertIn("oscillator", str(exc))
                else:
                    self.close(result.term_13a[0], ref, K_ABS)

    def test_cli_unresolved_and_invalid(self):
        # Exercise the actual CLI dispatcher, including an already existing directory.
        base = json.loads((ROOT / "examples/synthetic_bulk.json").read_text())
        base["loss_grid_eV"] = {"min": 20, "max": 20, "step": 1}
        base["material"]["oscillators"][0]["width_eV"] = 1e-5
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for eps, exists, status in ((1e-30, False, "not_converged"),
                                        (1e-30, True, "error"), (-1, False, "error")):
                with self.subTest(eps=eps, exists=exists):
                    base["integration"].update(epsrel=eps, epsabs=eps)
                    path = root / "config.json"
                    path.write_text(json.dumps(base))
                    out = root / f"result-{eps}-{exists}"
                    if exists:
                        out.mkdir()
                        (out / "sentinel").write_text("keep")
                    stderr = io.StringIO()
                    with contextlib.redirect_stderr(stderr):
                        rc = main(["bulk", str(path), "--out", str(out)])
                    self.assertNotEqual(rc, 0)
                    failure = json.loads(stderr.getvalue())
                    self.assertEqual(failure["status"], status)
                    if status == "not_converged":
                        self.assertIn("oscillator", failure["message"])
                        self.assertIn("loss_eV", failure["message"])
                    self.assertFalse((out / "run.json").exists())
                    if exists:
                        self.assertEqual(failure["error_type"], "FileExistsError")
                        self.assertEqual(failure["message"], f"output path already exists: {out}")
                        self.assertEqual((out / "sentinel").read_text(), "keep")
