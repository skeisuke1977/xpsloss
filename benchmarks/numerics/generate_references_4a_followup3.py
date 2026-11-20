"""Print R4 literals with python benchmarks/numerics/generate_references_4a_followup3.py.

Using mpmath at 80 digits, subtract the gap's direct real-axis KK integral
from the analytic gapless response. Both losses exceed the gap, so the
integral on [0, gap] is smooth. All parameters originate as binary64.
"""
import mpmath as mp


# Centers, common width and dispersion of the pair, q, gap, broad third line.
CASES = (
    (4.999982983456365, 5.000017016543635, 4.100498394012326e-5, 0., 0., 1., False),
    (4.999944554138466, 4.999974612266964, 2.223177927651792e-5, .01, .03257012904474573, 1., False),
    (4.99997247124135, 5.000019552585014, 5.992149831747236e-5, .01, .010231064750821165, 3., True),
    (3.1529704517641695, 3.153037025764051, .00011977626260052459,
     .03494293247578475, 3.7247048538669136, 3., True),
)


def references(case):
    w1, w2, width, alpha, q, gap, broad = case
    shift = mp.mpf(alpha) * mp.mpf(3.80998212) * mp.mpf(q)**2
    lines = [(mp.mpf(w) + shift, 100, mp.mpf(width)) for w in (w1, w2)]
    if broad:
        lines.append((10, 50, 5))

    def response(u):
        return mp.fsum(a / (c*c - u*u - mp.j*g*u) for c, a, g in lines)

    for t in (4, 5):
        bare = response(mp.mpf(t))
        removed = mp.quad(lambda u: 2*u*response(u).imag / (mp.pi*(u*u-t*t)), [0, gap])
        yield complex(float(bare.real - 1 - removed), float(bare.imag))


if __name__ == "__main__":
    mp.mp.dps = 80
    for case in CASES:
        print(tuple(references(case)))
