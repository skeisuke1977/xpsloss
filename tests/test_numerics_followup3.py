"""R4 closed-form regressions. See benchmarks/numerics/generate_references_4a_followup3.py."""
import unittest

from xpsloss import DielectricConvergenceError, Material, Oscillator, negative_inverse_dielectric_function
from xpsloss.dielectric import negative_inverse_dielectric_closed_form


# Pair centers, common width and dispersion, q, gap, broad third oscillator,
# and independent 80-digit references rounded to binary64 at T=4 and 5 eV.
# First two are the worst q=0 and q>0 pairs. Third and fourth are the
# worst gap-3-plus-broad cases for reviewer seeds 99 and 2026.
CASES = (
    (4.999982983456365, 5.000017016543635, 4.100498394012326e-5, 0., 0., 1., False,
     21.222222412071833+0.00040498749597868346j, -0.8497089794156738+577604.188358128j),
    (4.999944554138466, 4.999974612266964, 2.223177927651792e-5, .01, .03257012904474573, 1., False,
     21.222222328705158+0.00021957312880819635j, -1.378507381166783+636219.3252422858j),
    (4.99997247124135, 5.000019552585014, 5.992149831747236e-5, .01, .010231064750821165, 3., True,
     21.80133292630203+0.13471198894187486j, -0.21116452087616444+412736.7829882366j),
    (3.1529704517641695, 3.153037025764051, .00011977626260052459, .03494293247578475,
     3.7247048538669136, 3., True,
     21.801352290278608+0.13530314710812238j, -0.14255296171211293+255135.7604712584j),
)


class ClosedFormCancellationTests(unittest.TestCase):
    def check_components(self, value, ref):
        self.assertAlmostEqual(value.real, ref.real, delta=max(1e-9, 1e-7 * abs(ref.real)))
        self.assertAlmostEqual(value.imag, ref.imag, delta=max(1e-9, 1e-7 * abs(ref.imag)))

    def test_closed_form_cancelling_lines_and_off_resonance_control(self):
        for w1, w2, width, alpha, q, gap, broad, off, peak in CASES:
            with self.subTest(q=q, gap=gap):
                lines = tuple(Oscillator(w, 100, width, alpha) for w in (w1, w2))
                if broad:
                    lines += (Oscillator(10, 50, 5),)
                mat = Material("synthetic", gap, lines)
                self.check_components(negative_inverse_dielectric_closed_form(mat, 4, q), off)
                try:
                    values = negative_inverse_dielectric_closed_form(mat, [4, 5], q)
                except DielectricConvergenceError as exc:
                    self.assertIn("closed-form KK oscillator sum", str(exc))
                    self.assertIn("loss_eV=5, index=(1,)", str(exc))
                    self.check_components(negative_inverse_dielectric_function(mat, 5, q), peak)
                else:
                    self.check_components(values[0], off)
                    self.check_components(values[1], peak)
