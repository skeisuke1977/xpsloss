"""Retained-band D3 at fixed-depth nodes and resolution of above-gap surface lines."""

import math

import numpy as np

from .dielectric import CLOSED_ABS, CLOSED_CENTER_REL, CLOSED_REL, CRITICAL_SPLITTING_FRACTION
from .model import (
    COORDINATE_ULPS,
    FLOAT_EPSILON,
    HBAR2_OVER_2ME_EV_A2,
    NumericalConvergenceError,
    RECIPROCAL_ROUNDOFF,
    SUM_ROUNDOFF,
    _MIN_DAMPING,
    _accuracy,
)


def _closed_band(t, gap, center, strength, width):
    """Array form of D3; coalescing poles use the analytic divided difference."""
    half = width / 2
    split = np.sqrt(abs(half - center)) * np.sqrt(half + center)
    under = half < center
    p = np.asarray(half + np.where(under, 1j * split, split))
    jp = np.empty(p.shape, complex)
    jb = np.empty_like(jp)
    h, w = split[under], center[under]
    near = (gap - w) + half * (half / (w + h))
    far = gap + h
    norm = np.hypot(near, half)
    angle = 0.5 * (np.arctan2(half, near) + np.arctan2(half, far))
    log_ratio = np.empty_like(h)
    far_edge = gap > 2 * w
    log_ratio[far_edge] = 0.5 * np.log1p(4 * (gap / norm[far_edge]) * (h[far_edge] / norm[far_edge]))
    close = ~far_edge
    log_ratio[close] = np.log(np.hypot(far[close], half) / norm[close])
    jp[under] = (angle + 0.5j * log_ratio) / p[under]
    jb[under] = jp[under].conj()

    over = ~under
    b = center[over] * (center[over] / p.real[over])
    jp[over] = np.arctan(p.real[over] / gap) / p.real[over]
    jb[over] = np.arctan(b / gap) / b
    diff = 2 * np.where(under, 1j * split, split)
    critical = abs(diff) < width * CRITICAL_SPLITTING_FRACTION
    jd = np.empty_like(p)
    jd_scale = np.empty_like(center)
    h = gap / (gap**2 + center[critical] ** 2)
    d = diff[critical]
    ratio = h.astype(complex)
    nonzero = d != 0
    ratio[nonzero] = np.arctan(h[nonzero] * d[nonzero]) / d[nonzero]
    jd[critical] = (jb[critical] - ratio) / (p[critical] * width)
    jd_scale[critical] = (abs(jb[critical]) + abs(ratio)) / (abs(p[critical]) * width)
    other = ~critical
    jd[other] = (jb[other] - jp[other]) / (width * diff[other])
    jd_scale[other] = (abs(jb[other]) + abs(jp[other])) / (width * abs(diff[other]))

    detuning = (center - t) * (center + t)
    damping = width * t
    denominator = detuning**2 + damping**2
    scale = (2 / math.pi) * strength * width / denominator
    coefficient = center**4 + (t * p) ** 2
    regular = coefficient * jd
    log_term = (t / 2) * np.log1p(2 * (np.minimum(t, gap) / abs(gap - t)))
    value = scale * (regular + (log_term - t**2 * jb)).real
    magnitude = scale * (abs(regular) + log_term + t**2 * abs(jb))
    error = RECIPROCAL_ROUNDOFF * (
        np.maximum(1.0, center**2 / np.sqrt(denominator)) * magnitude + scale * abs(coefficient) * jd_scale
    )
    if not np.all(
        np.isfinite(value) & np.isfinite(error) & (error <= _accuracy(value, CLOSED_ABS, CLOSED_REL))
    ):
        raise NumericalConvergenceError("closed-form D3 is ill-conditioned")
    return value, np.array(error), denominator


def _removed_moment(g, t):
    """Upper bounds for integral_0^g u^2/(t^2-u^2) du, valid only for t>g."""
    return np.minimum(g**3 / (3 * (t - g) * (t + g)), (t / 2) * np.log1p(2 * g / (t - g)))


def gapped_sum(material, loss, q, relative_z):
    """S and separate Re/Im estimates, including the finite interval around a node."""
    t, q, rz = np.broadcast_arrays(loss, q, relative_z)
    g = material.gap_eV
    s = np.zeros(t.shape, complex)
    sr, si = np.zeros(t.shape), np.zeros(t.shape)
    d3_error = np.zeros(t.shape)
    above = t > g
    below = ~above
    active = sum(bool(o.strength_eV2) for o in material.oscillators)
    for o in material.oscillators:
        a, width = o.strength_eV2, o.width_eV
        if not a:
            continue
        if np.any(width * t[above] < _MIN_DAMPING):
            raise NumericalConvergenceError("hard-gap oscillator damping unresolved")
        shift = o.dispersion * HBAR2_OVER_2ME_EV_A2 * q**2
        center = o.energy_eV + shift
        if np.any(
            (shift > 0)
            & (COORDINATE_ULPS * np.spacing(center) > CLOSED_CENTER_REL * np.maximum(width, abs(g - center)))
        ):
            raise NumericalConvergenceError("closed-form D3 dispersed center unresolved at the gap")
        value, err, den = _closed_band(t, g, center, a, width)
        d3_error += err + active * SUM_ROUNDOFF * abs(value)
        # Undispersed centers are exact inputs. Otherwise include the addition
        # to the optical center as well as the log node, offset and recoil.
        dw = np.where(shift > 0, (rz + 8 * FLOAT_EPSILON) * shift + FLOAT_EPSILON * center, 0.0)
        dx = dw * (2 * center + dw) * (1 + 8 * FLOAT_EPSILON)
        if np.any(below):
            if np.any(value[below] + err[below] < 0):
                raise NumericalConvergenceError("positive retained band unresolved")
            err[below] += (value[below] + err[below]) * np.expm1(dx[below] / (width * g))
        if np.any(above):
            ta, xa, da = t[above], center[above], dx[above]
            m = np.nextafter((xa - g) * (xa + g) * (1 - 8 * FLOAT_EPSILON) - da, -np.inf)
            if np.any(m <= 0):
                raise NumericalConvergenceError("removed-band derivative interval reaches the gap")
            # Minimum |D| throughout the interval; the damping is constant.
            dist = np.maximum(width * ta, np.sqrt(den[above]) * (1 - 8 * FLOAT_EPSILON) - da)
            imag = a * (width * ta) / den[above]
            bare_change = (a / dist) * (da / dist)
            # d(Im S)/dx = -2 Im S * detuning/|D|^2.
            imag_change = (2 * a * (width * ta) / dist**2) * (da / dist)
            removed_change = (4 / math.pi) * a * width * _removed_moment(g, ta) * (da / m) / m**2
            err[above] += bare_change + removed_change
            s.imag[above] += imag
            si[above] += imag_change + (RECIPROCAL_ROUNDOFF + active * SUM_ROUNDOFF) * imag
        s.real += value
        sr += err + active * SUM_ROUNDOFF * abs(value)
    if np.any(d3_error + active * SUM_ROUNDOFF > _accuracy(s.real - 1.0)):
        raise NumericalConvergenceError("closed-form D3 oscillator sum is ill-conditioned")
    if not np.all(np.isfinite(s) & np.isfinite(sr) & np.isfinite(si)):
        raise NumericalConvergenceError("hard-gap response outside numerical range")
    return s, sr, si


def _surface_intervals(material, t, left, right, offset):
    """Enclose Re S and total complex variation on each log panel (z=k^2).

    The real Lorentzian extrema occur at detuning +/- damping. The
    removed-band correction is positive above g. Its derivative augments
    the full Lorentzian derivative in the variation bound.
    """
    g = material.gap_eV
    zlo, zhi = np.exp(2 * left) + offset**2, np.exp(2 * right) + offset**2
    lo, hi, left_lo, left_hi, variation, im_min = (np.zeros(left.shape) for _ in range(6))
    guard = 32 * FLOAT_EPSILON * (1 + np.maximum(abs(left), abs(right)))
    moment = _removed_moment(g, t)
    for o in material.oscillators:
        if not o.strength_eV2:
            continue
        slope = o.dispersion * HBAR2_OVER_2ME_EV_A2
        a, width = o.strength_eV2, o.width_eV
        wlo = (o.energy_eV + slope * zlo) * (1 - guard)
        whi = (o.energy_eV + slope * zhi) * (1 + guard)
        # Enclose the factored detuning, including the product rounding.
        dl = (wlo - t) * (wlo + t)
        dh = (whi - t) * (whi + t)
        dl -= guard * abs(dl)
        dh += guard * abs(dh)
        b = width * t

        def real(d):
            return a * d / (d * d + b * b)

        fl, fh = real(dl), real(dh)
        minimum = np.where((dl <= -b) & (-b <= dh), -a / (2 * b), np.minimum(fl, fh))
        maximum = np.where((dl <= b) & (b <= dh), a / (2 * b), np.maximum(fl, fh))
        m = (wlo - g) * (wlo + g)
        if np.any(m <= 0):
            raise NumericalConvergenceError("surface search cannot bound the removed band")
        correction = (2 / math.pi) * a * width * moment / m**2
        lo = np.nextafter(lo + (minimum - guard * (abs(minimum) + 1)), -np.inf)
        hi = np.nextafter(hi + (maximum + correction + guard * (abs(maximum) + correction + 1)), np.inf)
        # Also enclose the sum at the lower detunings and propagate over
        # the whole panel with the complex derivative. The second enclosure
        # accounts for the joint variation of the oscillator sum.
        left_lo = np.nextafter(left_lo + (fl - guard * (abs(fl) + 1)), -np.inf)
        left_hi = np.nextafter(left_hi + (fl + correction + guard * (abs(fl) + correction + 1)), np.inf)
        dx = np.nextafter(np.maximum((whi - wlo) * (whi + wlo), dh - dl), np.inf)
        distance = np.maximum(0.0, np.maximum(dl, -dh))
        change = dx * (a / (distance**2 + b**2) + 2 * correction / m) * (1 + guard)
        variation = np.nextafter(variation + change, np.inf)
        lower_imag = (a * b / (np.maximum(abs(dl), abs(dh)) ** 2 + b * b)) * (1 - guard)
        im_min = np.nextafter(im_min + lower_imag, 0.0)
    lo = np.maximum(lo, np.nextafter(left_lo - variation, -np.inf))
    hi = np.minimum(hi, np.nextafter(left_hi + variation, np.inf))
    if not np.all(np.isfinite(lo) & np.isfinite(hi) & np.isfinite(variation) & (im_min > 0)):
        raise NumericalConvergenceError("surface search bounds unresolved")
    # For a remaining candidate, change in |2-S| is at most the complex
    # variation. This also resolves a tangency or two crossings in a panel.
    return (lo <= 2) & (2 <= hi) & (variation > 0.5 * im_min)


def resolve_surface(material, loss, edges, offset, limit):
    if loss <= material.gap_eV or not any(o.strength_eV2 and o.dispersion for o in material.oscillators):
        return edges
    left, right = edges[:-1], edges[1:]
    parts = [edges]
    count = len(left)
    while left.size:
        candidate = _surface_intervals(material, loss, left, right, offset)
        left, right = left[candidate], right[candidate]
        mids = (left + right) / 2
        count += len(mids)
        if count > limit or np.any((mids <= left) | (mids >= right)):
            raise NumericalConvergenceError("surface resonance unresolved within panel limit")
        parts.append(mids)
        left, right = np.concatenate((left, mids)), np.concatenate((mids, right))
    return np.sort(np.concatenate(parts))
