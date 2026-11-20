"""Reproduce the cancellation-pair test literals with Python and mpmath 1.3.0.

Run python benchmarks/numerics/generate_references_4a_followup2.py. Combine the
two oscillators over one common complex denominator at 80 decimal digits. For
a 1 eV gap subtract the removed band's direct KK integral. Inputs are the
binary64 coordinates of the indicated zero-based fu1_cancel_pair.py trials.
"""
import math
import random

import mpmath as mp


def generate():
    mp.mp.dps = 80
    rng = random.Random(1)
    for trial in range(152):
        width = 10 ** rng.uniform(-8, -4)
        detuning = rng.uniform(0.2, 2) * width / 2
        target = detuning / 3.55e-5 * rng.uniform(0.3, 0.99)
        q = math.sqrt(target / (0.01 * 3.80998212))
        if trial not in (2, 3, 50, 151):
            continue
        shift = 0.01 * 3.80998212 * q * q
        w1, w2 = 15 - shift - detuning, 15 - shift + detuning
        exact_shift = mp.mpf(0.01) * mp.mpf(3.80998212) * mp.mpf(q)**2
        c1, c2 = mp.mpf(w1) + exact_shift, mp.mpf(w2) + exact_shift

        def pair(u):
            z = u**2 + mp.j * mp.mpf(width) * u
            d1, d2 = c1**2 - z, c2**2 - z
            return 100 * (d1 + d2) / (d1 * d2)

        full = pair(mp.mpf(15)) - 1
        removed = mp.quad(lambda u: 2*u*pair(u).imag / (mp.pi*(u*u - 15**2)), [0, 1])
        print((trial, w1, w2, width, q, float(full.real), float(full.real - removed), float(full.imag)))


if __name__ == "__main__":
    generate()
