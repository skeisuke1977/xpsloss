"""Depth factors for the eight normal-emission densities (M3A_DESIGN section 3.1)."""

from fractions import Fraction
import math

import numpy as np


_EPS = np.finfo(float).eps


def fixed_depth_weights(k, omega, depth):
    denom = omega**2 + k**2
    plus = (omega + 1j*k) / denom
    minus = plus.conj()
    # Contract conjugate phases first. Every remaining phase carries exp(-ka).
    attenuation = -k * depth
    decay = np.exp(attenuation)
    phase = -1j * (omega * depth)
    u = decay * np.exp(phase)
    diff = -np.expm1(attenuation) - decay * np.expm1(phase)
    norm = abs(diff)**2 / denom
    mixed = u * diff.conj()                     # u − exp(-2ka)

    w13b = 2 * (plus**2 * diff).imag
    weight_14 = plus*diff + minus*mixed
    w15b = minus * mixed
    w15c = minus**2 * diff.conj()
    w16e = plus * (2j*k/denom - plus*u)
    w16h = plus * u
    return denom, plus, u, diff, norm, w13b, weight_14, w15b, w15c, w16e, w16h


def assemble_terms(response, reflection, surface, pref, omega, w13a, w13b,
                   norm, weight_14, w15b, w15c, w16e, w16h):
    """Project complex response coefficients, preserving the fixed-depth reduction order."""
    bulk = -response
    term_13a = w13a * response.imag
    # The A and B terms of 13b form a conjugate pair. C supplies |1-u|².
    term_13b = pref * (bulk.imag * w13b + bulk.real * norm)
    term_14 = pref * (bulk * weight_14).real / omega
    term_15_T1 = pref * surface.real * norm
    term_15_T2 = pref * (surface * w15b).real / omega
    term_15_T3 = pref * (reflection * w15c).real
    term_16_e = -pref * (reflection * w16e).real
    term_16_h = -pref * (reflection * w16h).real / omega
    terms = np.stack((term_13a, term_13b, term_14, term_15_T1, term_15_T2,
                      term_15_T3, term_16_e, term_16_h), axis=-2)
    return terms


class _Bound:
    """A Cartesian center and absolute radius. Positive divisors stay away from zero."""

    def __init__(self, value, error, rounding):
        self.v = value
        self.e = error
        self.u = rounding

    def add(self, other, sign=1):
        value = self.v + sign * other.v
        error = self.e + other.e + self.u * abs(value)
        return _Bound(value, (1 + 8*self.u)*error, self.u)

    def mul(self, other):
        value = self.v * other.v
        error = abs(self.v)*other.e + (abs(other.v) + other.e)*self.e
        return _Bound(value, (1 + 8*self.u)*(error + self.u*abs(value)), self.u)

    def div(self, other):
        value = self.v / other.v
        # |x/y - X/Y| <= (e_x + |x/y| e_y) / (|y| - e_y).
        lower = abs(other.v) - other.e
        if self.u:
            lower = np.nextafter(lower, -math.inf)
        if np.any(lower <= 0):
            return _Bound(value, math.inf, self.u)
        error = (self.e + (1 + self.u)*abs(value)*other.e) / lower
        return _Bound(value, (1 + 8*self.u)*(error + self.u*abs(value)), self.u)


def _weights(k, b, p):
    # p is already inside every numerator. Only n is a signed difference.
    L = p.add(k)
    J = L.add(k)
    bb = b.mul(b)
    D = k.mul(k).add(bb)
    H = L.mul(L).add(bb)
    n = k.mul(L).add(bb, -1)
    two_b = b.add(b)
    common = p.div(two_b.mul(D).mul(H))
    base = p.div(J.mul(H))
    zero = _Bound(0., 0., p.u)
    w13r = common.mul(n)
    w13i = common.mul(b).mul(J)
    norm = base.mul(L).div(b)
    w14r = base.mul(k)
    w15br = base.mul(p).div(_Bound(2, 0, p.u))
    w15bi = w15br.mul(L).div(b)

    def parts(x, y, sign_x=1, sign_y=1):
        return sign_x*x.v, sign_y*y.v, x.e, y.e

    return {
        "W13": parts(w13r, w13i, -1), "N": parts(norm, zero),
        "W14": parts(w14r, w14r.mul(L).div(b)),
        "W15b": parts(w15br, w15bi, -1, -1),
        "W15c": parts(w13r, w13i, -1, -1),
        "W16e": parts(common.mul(J).mul(L), w13i, -1),
        "W16h": parts(p.mul(w13i), p.mul(w13r)),
    }


def _extended_weights(k, b, p, e_k, e_b, e_p):
    """Exact rational intermediates when a product may leave the normal exponent range."""
    args = [_Bound(Fraction(float(x)), Fraction(float(e)), 0)
            for x, e in ((k, e_k), (b, e_b), (p, e_p))]
    result = {}
    for key, (x, y, ex, ey) in _weights(*args).items():
        try:
            a, c = float(x), float(y)
        except OverflowError as exc:
            raise FloatingPointError("scaled depth weight overflow") from exc
        if (x and not a) or (y and not c):
            raise FloatingPointError("scaled depth weight underflow")
        def radius(v, rounded, error):
            r = error + abs(v - Fraction(rounded))
            if not r or r == math.inf:
                return float(r)
            try:
                return math.nextafter(float(r), math.inf)
            except OverflowError:
                return math.inf
        result[key] = a, c, radius(x, a, ex), radius(y, c, ey)
    return result


def depth_averaged_weights(k, b, p, *, e_k=0.0, e_b=0.0, e_p=0.0):
    """Scaled p**2 W/(2b), as (real, imaginary, absolute real/imaginary errors).

    The same operations accept a scalar or broadcast arrays of momentum nodes.
    In the normal range reserve one epsilon for center rounding, widen each
    computed radius for its positive arithmetic, and round divisor bounds down.
    """
    if not np.isfinite(p):
        raise FloatingPointError("reciprocal length outside numerical range")
    if p <= 0 or np.min(b) <= 0 or np.min(k) < 0:
        raise ValueError("depth kernel requires k >= 0, b > 0 and p > 0")
    if min(np.min(k), np.min(b), p) >= 2.0**-40 and max(np.max(k), np.max(b), p) <= 2.0**40:
        args = [_Bound(x, e, _EPS) for x, e in ((k, e_k), (b, e_b), (p, e_p))]
        return {key: (x, y, ex*(1 + 32*_EPS), ey*(1 + 32*_EPS))
                for key, (x, y, ex, ey) in _weights(*args).items()}
    arrays = np.broadcast_arrays(k, b, p, e_k, e_b, e_p)
    if not arrays[0].ndim:
        return _extended_weights(*arrays)
    rows = [_extended_weights(*point) for point in zip(*(a.flat for a in arrays))]
    return {name: tuple(np.array(v).reshape(arrays[0].shape) for v in zip(*(r[name] for r in rows)))
            for name in rows[0]}
