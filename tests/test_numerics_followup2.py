"""R1 pairs, from benchmarks/numerics/generate_references_4a_followup2.py."""
import unittest

from xpsloss import Material, Oscillator, negative_inverse_dielectric_function


# Trial, two optical centers, common width, q, Re F (gaps 0 and 1), Im F.
# References use a single common denominator at 80 digits and direct KK
# quadrature for the removed band. Trial 151 was the largest reported error.
PAIRS = (
    (2, 14.96637967944886, 14.966386222902322, 4.039911867566605e-6, 0.9393302828278689,
     -1.0721654346955802, -1.0721654346804064, 910849.0502199461),
    (3, 14.99981341940762, 14.999813441536304, 1.2983443910253402e-8, 0.06997757488524128,
     -16.404648264624118, -16.40464826462407, 262989714.82088336),
    (50, 14.99998411304997, 14.999984115862214, 1.2113577830901672e-8, 0.02041924010161186,
     -96.78251635383407, -96.78251635383403, 1044403525.1893809),
    (151, 14.998881600906143, 14.998881686550286, 2.00386798800049e-7, 0.17132810843264237,
     -0.6349596891941269, -0.6349596891933743, 56261021.188386776),
)


class DetuningCancellationTests(unittest.TestCase):
    def test_dispersed_pairs_meet_each_component_tolerance(self):
        for trial, w1, w2, width, q, gapless, gapped, imag in PAIRS:
            lines = tuple(Oscillator(w, 100, width, 0.01) for w in (w1, w2))
            for gap, real in ((0, gapless), (1, gapped)):
                with self.subTest(trial=trial, gap=gap):
                    material = Material("synthetic", gap, lines)
                    response = negative_inverse_dielectric_function(material, [14.25, 15.0], q)[1]
                    self.assertAlmostEqual(response.real, real, delta=max(1e-9, 1e-7 * abs(real)))
                    self.assertAlmostEqual(response.imag, imag, delta=max(1e-9, 1e-7 * abs(imag)))
