"""Generate the 4f literals offline: python benchmarks/numerics/generate_references_4f.py.

mpmath only. Form epsilon from a common oscillator denominator, then integrate
the uncontracted braces in Simonsen et al. (1997), Eqs. (13)-(16), directly in
k_parallel. Eq. (17) is evaluated by subtracting squared momenta. Two different
orders and working precisions must give identical rounded binary64 literals.
No xpsloss or SciPy functions, package partitions, or error estimates are used.
"""
import math

import mpmath as mp


T = math.sqrt(135.0)  # The actual binary64 API argument, not an exact radical.
LINE = (15, 180, 3e-7, 0)
CASES = {
    "owner": (T, [LINE], False),
    "broad": (T, [(15, 180, 3e-4, 0)], False),
    "middle": (T, [(15, 180, 3e-6, 0)], False),
    "narrow": (T, [(15, 180, 3e-8, 0)], False),
    "below": (T - 1.5e-7, [LINE], False),
    "above": (T + 1.5e-7, [LINE], False),
    "parallel": (T, [(15, 180, 3e-7, 2e-8)], False),
    "total_momentum": (T, [(15, 180, 3e-7, 2e-8)], True),
    # 81/(15^2-12^2) + 180/(18^2-12^2) = 2 at zero damping.
    "two": (12., [(15, 81, 4e-7, 0), (18, 180, 2e-7, 0)], False),
}


def fixed(loss, lines, total_momentum, rule):
    t, h, a0 = map(mp.mpf, (loss, 3.80998212, 0.529177210903))
    energy, depth = 1000, 10
    omega = t / (2 * mp.sqrt(h * energy))
    low, high = (mp.sqrt((mp.sqrt(energy) + sign * mp.sqrt(energy - t))**2 / h - omega**2)
                 for sign in (-1, 1))
    phase = mp.exp(mp.j * omega * depth)
    prefactor = 1 / (2 * mp.pi * a0 * energy * depth * omega)
    lines = [tuple(map(mp.mpf, line)) for line in lines]

    def densities(k):
        q2 = k**2 + (omega**2 if total_momentum else 0)
        den = [(w + alpha * h * q2)**2 - t**2 - mp.j * g * t for w, _, g, alpha in lines]
        product = mp.fprod(den)
        numerator = mp.fsum(a * mp.fprod(den[:i] + den[i+1:]) for i, (_, a, _, _) in enumerate(lines))
        epsilon = product / (product - numerator)
        bulk = 1 / epsilon - 1
        reflection = (epsilon - 1) / (epsilon + 1)
        surface = (epsilon - 1) * (epsilon + 2) / (epsilon * (epsilon + 1))
        p, m = omega + mp.j*k, omega - mp.j*k
        n, decay = p*m, mp.exp(-k*depth)
        braces = (p**2 * (decay/phase - 1) - m**2 * (decay*phase - 1)
                  + (1/phase - decay) * (phase - decay) * n)
        pre15 = -m * (decay - 1/phase) / n
        pre16 = -reflection * p / (phase*n)
        terms = (
            2*depth*omega * mp.j * bulk/n,
            bulk * braces/n**2,
            bulk * (p*(1-decay/phase) - m*(decay**2-decay/phase)) / (omega*n),
            pre15 * surface * (phase-decay)*p/n,
            pre15 * surface * decay/omega,
            pre15 * reflection * phase*m/n,
            pre16 * (2*mp.j*k*phase - p*decay)/n,
            pre16 * decay/omega,
        )
        return [prefactor * k * z.real for z in terms]

    terms = [mp.mpf(0)] * 8
    left = low
    while left < high:
        right = min(2 * left, high)
        for node, weight in zip(*rule):
            k = ((right-left)*node + right+left)/2
            for i, value in enumerate(densities(k)):
                terms[i] += (right-left)*weight*value/2
        left = right
    electron = mp.fsum(terms[i] for i in (0, 1, 3, 5, 6))
    hole = mp.fsum(terms[i] for i in (2, 4, 7))
    return tuple(map(float, [*terms, electron, hole, electron+hole]))


def other_references():
    for loss, gamma in ((T, 3e-7), (math.sqrt(105.), 3e-7), (math.nextafter(15., 0.), 3e-18)):
        u = mp.mpf(loss)
        den = 225 - u**2 - mp.j * mp.mpf(gamma) * u
        eps = den / (den - 180)
        reflection = (eps - 1) / (eps + 1)
        surface = (eps - 1) * (eps + 2) / (eps * (eps + 1))
        print("R and C", loss, gamma, complex(reflection), complex(surface))

    t, width = mp.mpf(math.sqrt(45.0)), mp.mpf(3e-7)

    def band(u):
        return 180 / (225 - u*u - mp.j*width*u)

    for gap in (0, 1):
        removed = mp.quad(lambda u: 2*u*band(u).imag/(mp.pi*(u*u-t*t)), [0, gap])
        f = band(t) - 1 - removed
        print("F and epsilon, gap", gap, complex(f), complex(-1/f))
    e, t = 1000, mp.mpf(math.nextafter(1000., 0.))
    h, a0 = map(mp.mpf, (3.80998212, 0.529177210903))
    q0, q1 = ((mp.sqrt(e) + sign*mp.sqrt(e-t))/mp.sqrt(h) for sign in (-1, 1))

    def integrand(q):
        center = 15 + mp.mpf(1e-20)*h*q*q
        return (180/(center**2-t*t-mp.j*mp.mpf(0.1)*t)).imag/q
    print("bulk", float(mp.quad(integrand, [q0, q1])/(mp.pi*a0*e)))

    for gap, gamma in ((100., 2.), (1e6, 2.), (1e7, 30.), (1e7, 100.)):
        # Direct static KK on the retained band, x = gap/u.
        strength = mp.mpf(gap**3 / 10)
        def tail(x):
            return x*x / ((1 - (15*x/gap)**2)**2 + (gamma*x/gap)**2)

        f = -1 + 2 * strength * gamma / (mp.pi * mp.mpf(gap)**3) * mp.quad(tail, [0, 1])
        print("closed static", gap, gamma, float(f))


if __name__ == "__main__":
    previous = {}
    for digits, order in ((50, 32), (80, 48)):
        mp.mp.dps = digits
        rule = mp.gauss_quadrature(order, "legendre")
        for name, args in CASES.items():
            result = fixed(*args, rule)
            if name in previous and result != previous[name]:
                raise ArithmeticError(f"{name}: precision or order has not converged")
            previous[name] = result
            if digits == 80:
                print(name, result)
    other_references()
