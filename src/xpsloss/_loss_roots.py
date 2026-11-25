"""Resolution candidates at fixed and moving momentum endpoints (M4a design 3.2)."""

import cmath
from fractions import Fraction
import math

import numpy as np
from numpy.polynomial import Polynomial as Poly


def distance(z, left, right):
    return math.hypot(max(left-z.real, 0., z.real-right), z.imag)


def _families(oscillators, energy_polynomials, loss):
    denominators = [w*w - loss*loss - 1j*osc.width_eV*loss
                    for osc, w in zip(oscillators, energy_polynomials)]
    for i, denominator in enumerate(denominators):
        yield "bulk", i, denominator
    product = Poly([1])
    numerator = Poly([0])
    for osc, denominator in zip(oscillators, denominators):
        numerator = numerator*denominator + osc.strength_eV2*product
        product = product*denominator
    for name, constant in (("surface", 2), ("epsilon", 1)):
        yield name, constant, constant*product - numerator


def _response(z, energy, mode, end, lines, family, index, *, factored=False):
    """Original rational equation and its z derivative on the physical branch."""
    shift = derivative = 0j
    if end != "both":
        s = cmath.sqrt(1-z/energy)
        if s == 0:
            raise ZeroDivisionError
        t = z/energy/(1+s)
        if end == "lower":
            if mode == "k_parallel":
                shift = energy*t**3*(1-t/4)
                derivative = t*t*(3-t)/(2*s)
            else:
                shift, derivative = energy*t*t, t/s
        else:
            shift, derivative = energy*(1+s)**2, -(1+s)/s
            if mode == "k_parallel":
                shift -= z*z/(4*energy)
                derivative -= z/(2*energy)
    terms, slopes = [], []
    cancelled = False
    for i, osc in enumerate(lines):
        w = osc.energy_eV + osc.dispersion*shift
        real_part = (w-z)*(w+z) if factored else w*w-z*z
        d = real_part - 1j*osc.width_eV*z
        dp = 2*w*osc.dispersion*derivative - 2*z - 1j*osc.width_eV
        scale = abs(w)**2 + abs(z)**2 + osc.width_eV*abs(z)
        if family == "bulk":
            if i == index:
                return d, dp, abs(d)/scale, False
        else:
            cancelled |= abs(d) <= 1e-10*scale
            if d == 0:
                return 0j, 0j, math.inf, True
            terms.append(osc.strength_eV2/d)
            slopes.append(terms[-1]*dp/d)
    f = index - sum(terms)
    return f, sum(slopes), abs(f)/(1+sum(map(abs, terms))), cancelled


def _polish(root, neighbors, equation):
    # Keep each iterate in the original root's Voronoi cell. A Newton
    # trajectory attracted to a different simple root cannot replace this one.
    z = root
    for iteration in range(12):
        value, slope, residual, cancelled = equation(z)
        if cancelled or not slope or (iteration and residual < 8*np.finfo(float).eps):
            break
        trial = z - value/slope
        if (abs(trial-root) > (1+abs(root))/2 or
                any(abs(trial-other) < abs(trial-root) for other in neighbors)):
            z = root
            break
        if trial == z:
            break
        z = trial
    _, _, residual, cancelled = equation(z)
    return z, residual, cancelled


def _remote_single_root(t0, energy, mode, end, osc, family, index):
    """Confirm a remote, ill-conditioned one-oscillator root in local t.

    A correction can be smaller than an ulp of z while changing A/D by
    order one. Carry the real and imaginary parts of t as exact rationals.
    The 2-by-2 multiplication matrices implement complex arithmetic.
    """
    def matrix(x, y=0):
        x, y = Fraction(x), Fraction(y)
        return np.array([[x, -y], [y, x]], dtype=object)

    def scalar(m):
        return complex(float(m[0, 0]), float(m[1, 0]))

    t = matrix(t0.real, t0.imag)
    identity = matrix(1)
    e = Fraction(energy)
    target = 0 if family == "bulk" else Fraction(osc.strength_eV2)/index
    for _ in range(8):
        z = e*(t @ (2*identity-t))
        zp = 2*e*(identity-t)
        u = t if end == "lower" else 2*identity-t
        q = e*(u @ u)
        qp = (2 if end == "lower" else -2)*e*u
        if mode == "k_parallel":
            q -= (z @ z)/(4*e)
            qp -= (z @ zp)/(2*e)
        w = matrix(osc.energy_eV)+Fraction(osc.dispersion)*q
        wp = Fraction(osc.dispersion)*qp
        damping = matrix(0, osc.width_eV)
        d = (w @ w)-(z @ z)-(damping @ z)
        dp = 2*(w @ wp)-2*(z @ zp)-(damping @ zp)
        v, slope = scalar(d-matrix(target)), scalar(dp)
        scale = abs(scalar(w))**2+abs(scalar(z))**2+osc.width_eV*abs(scalar(z))
        if family == "bulk":
            residual, cancelled = abs(v)/scale, False
        else:
            s = osc.strength_eV2/scalar(d)
            residual, cancelled = abs(index-s)/(1+abs(s)), abs(scalar(d)) <= 1e-10*scale
        if residual <= 1e-8 or not slope:
            break
        step = v/slope
        t -= matrix(step.real, step.imag)
    return scalar(z), residual, cancelled


def candidates(material, energy, mode, window):
    lines = [o for o in material.oscillators if o.strength_eV2]
    fixed = all(o.dispersion == 0 for o in lines)
    left, right = window
    records = []
    all_confirmed = True
    t = Poly([0, 1])
    for end in (("both",) if fixed else ("lower", "upper")):
        # At either end use the local expansion about z=0. For the upper
        # branch its momentum variable is 2-t, where t=1-sqrt(1-z/E0).
        zpoly = t if fixed else energy*t*(2-t)
        momentum = t if end == "lower" else 2-t
        shift = energy*momentum**2
        if mode == "k_parallel":
            shift -= zpoly**2/(4*energy)
        ws = [o.energy_eV + o.dispersion*shift for o in lines]
        for family, index, polynomial in _families(lines, ws, zpoly):
            variables = [complex(v) for v in polynomial.trim().roots() if fixed or v.real < 1]
            roots = [complex(zpoly(v)) for v in variables]
            for i, root in enumerate(roots):
                def equation(z):
                    return _response(z, energy, mode, end, lines, family, index)
                z, residual, cancelled = _polish(root, roots[:i]+roots[i+1:], equation)
                if not cancelled:
                    # Check the rational equation in both algebraic forms at
                    # nearly coincident poles, where D itself loses digits.
                    residual = max(residual, _response(z, energy, mode, end, lines, family, index,
                                                       factored=True)[2])
                if (not fixed and len(lines) == 1 and residual > 1e-8
                        and distance(root, left, right) >= right-left):
                    trial, err, cancel = _remote_single_root(variables[i], energy, mode, end,
                                                             lines[0], family, index)
                    if all(abs(trial-root) <= abs(trial-other) for other in roots[:i]+roots[i+1:]):
                        z, residual, cancelled = trial, err, cancel
                if cancelled:
                    continue
                confirmed = residual <= 1e-8
                all_confirmed &= confirmed
                dist = distance(z, left, right)
                if dist >= right-left:
                    continue
                if any(row["family"] == family and row["end"] == end
                       and abs(z-complex(row["z_re"], row["z_im"])) <= 1e-9*(1+abs(z))
                       for row in records):
                    continue
                records.append(dict(z_re=z.real, z_im=z.imag, family=family, end=end,
                                    residual=residual, confirmed=confirmed, near_origin=abs(z)<left/2,
                                    method="enumerated", distance_to_window_eV=dist))
    if material.gap_eV or not all_confirmed:
        return records, "computed_candidates_only", ("hard_gap" if material.gap_eV else "unconfirmed_root")
    return records, ("fixed_poles_enumerated" if fixed else "endpoint_roots_enumerated"), None

