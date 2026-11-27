"""Cumulative truncation bounds and fixed-multiplier interval arithmetic."""

from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
import math

import numpy as np
from numpy.polynomial import chebyshev as cheb

from ._fold import _up
from ._order_quadrature import area_terms, tau_grid
from ._order_types import MultipleLossConvergenceError
from ._single_loss_types import COLUMNS, EPS


def outward(lo, hi):
    return -_up(-lo), _up(hi)


def times(x, y):
    corners = [a * b for a, b in product(x, y)]
    return min(corners), max(corners)


def around(x, e):
    return Q(x) - Q(e), Q(x) + Q(e)


def area_interval(single, name):
    bounds = single.multiplier_interval_A or (single.multiplier_A,) * 2
    return outward(*times(tuple(map(Q, bounds)), around(single.window_area[name], single.window_area_error[name])))


def ratio_record(single):
    v = float(Q(single.multiplier_A) * Q(single.window_area["term_13a"]))
    interval = area_interval(single, "term_13a")
    return dict(value=v, interval=interval, error=_up(max(Q(v) - Q(interval[0]), Q(interval[1]) - Q(v))))


def multiplier_rho(single, n):
    if single.multiplier_interval_A is None:
        return Q(0)
    lo, hi = (Q(x) / Q(single.multiplier_A) for x in single.multiplier_interval_A)
    return Q(_up(max(1 - lo**n, hi**n - 1)))


def corrected(value, fixed, rho):
    extra = rho * (abs(Q(value)) + Q(fixed))
    return _up(Q(fixed) + extra), _up(extra)


def error_from_parts(value, parts, rho):
    """The higher-order total encloses the sum of the reported parts."""
    fixed = sum(map(Q, parts))
    extra = _up(rho * (abs(Q(value)) + Q(_up(fixed))))
    return _up(fixed + Q(extra)), extra


def sensitivity(first, panels, seeds, single, n, rho):
    center, error = {}, {}
    h = single.cutoff_sensitivity
    m, a = Q(single.multiplier_A), Q(single.loss_window_eV[0])
    he = {c: _up(m * a * Q(e) + 2 * Q(EPS) * abs(Q(h[c]))) for c, e in zip(COLUMNS, panels[0].errors[:, 0])}
    r, hr = Q(first.area["term_13a"]), Q(h["term_13a"])
    ri = around(r, seeds["term_13a"])

    def power(interval, k):
        lo, hi = interval
        ends = [lo**k, hi**k]
        if k and k % 2 == 0 and lo < 0 < hi:
            ends.append(Q(0))
        return min(ends), max(ends)

    rp, rq = power(ri, n - 1), power(ri, n - 2)
    for c in COLUMNS:
        w = Q(first.area[c])
        value = float(Q(h[c]) * r ** (n - 1) + (n - 1) * w * hr * r ** (n - 2))
        base = times(around(h[c], he[c]), rp)
        term = times(times(around(w, seeds[c]), around(hr, he["term_13a"])), rq)
        fixed = _up(max(Q(value) - base[0] - (n - 1) * term[0], base[1] + (n - 1) * term[1] - Q(value)))
        center[c], error[c] = value, corrected(value, fixed, rho)[0]
    return center, error


def log_up(value):
    """Round toward +infinity, preserving log(0) for an absent term."""
    return np.where(np.isneginf(value), -np.inf, np.nextafter(value, np.inf))


def log_geometric_factor(logr):
    """Upper estimate of -log(1-R) for an upper log(R) strictly below zero."""
    denominator = np.nextafter(-np.expm1(logr), 0.0)
    return log_up(-np.log(denominator))


def log_power(log_value, n):
    if not n:
        return np.zeros_like(log_value)
    # The support-based tolerance search can supply an integer beyond 2**53.
    limit = np.finfo(float).max
    lo, hi = outward(Q(n), Q(n)) if n <= int(limit) else (limit, math.inf)
    scaled = np.where(log_value > 0, hi, lo) * log_value
    return np.where(np.isneginf(log_value), -np.inf, np.nextafter(scaled, np.inf))


def exp_bound(log_value):
    """Convert a log upper estimate, rounding a positive underflow to the first subnormal."""
    if log_value == -math.inf:
        return 0.0
    if not log_value <= math.log(np.finfo(float).max):
        return None
    bound = math.nextafter(math.exp(log_value), math.inf)
    return bound if math.isfinite(bound) else None


def log_laplace(panels, col, taus, m):
    """Integrate an upper envelope against the exponential, separately from the L1 error."""
    total = np.full(len(taus), -np.inf)
    for p in panels:
        a, b = p.left_eV, p.right_eV
        values = p.values[col]
        if np.all(values == values[0]):
            edges, heights = [a, b], [abs(values[0])]
        else:
            # Bound |P| on each small interval by its linear interpolant plus
            # h^2 max|P''|/8. Integrate the exponential analytically there.
            coeff = p.coefficients[col]
            x = np.linspace(-1, 1, 65)
            y = abs(cheb.chebval(x, coeff))
            remainder = abs(cheb.chebder(coeff, 2)).sum() * (2 / 64) ** 2 / 8
            heights = np.maximum(y[:-1], y[1:]) + remainder + 64 * EPS * abs(coeff).sum()
            edges = np.linspace(a, b, 65)
        for left, right, height in zip(edges[:-1], edges[1:], heights):
            w = _up(Q(right) - Q(left))
            if height == 0 or w == 0:
                continue
            log_weights = np.full(len(taus), math.nextafter(math.log(w), math.inf))
            nz = taus > 0
            z = np.nextafter(taus[nz] * w, np.inf)
            small = z < EPS
            # For very small tau*w the width itself is an upper estimate.
            use = np.flatnonzero(nz)[~small]
            numerator = np.nextafter(-np.expm1(-z[~small]), np.inf)
            log_weights[use] = log_up(log_up(np.log(numerator)) - np.nextafter(np.log(taus[use]), -np.inf))
            log_mass = log_up(math.nextafter(math.log(height), math.inf) + log_weights)
            attenuation = np.maximum(0.0, np.nextafter(taus * left, -np.inf))
            total = log_up(np.logaddexp(total, log_up(log_mass - attenuation)))
        error = area_terms(p)[col]
        if error:
            attenuation = np.maximum(0.0, np.nextafter(taus * a, -np.inf))
            log_error = log_up(math.nextafter(math.log(error), math.inf) - attenuation)
            total = log_up(np.logaddexp(total, log_error))
    return log_up(log_up(total + math.nextafter(math.log(m), math.inf)) + math.log1p(64 * EPS))


def tail_area(single, ratio, n, negative, first_zero, repeat_zero):
    mass = Q(single.primary_integral)
    w = Q(single.single_loss_weight)
    wi = tuple(map(Q, area_interval(single, "total")))

    def result(value=None, interval=None, reason=None):
        return dict(value=value, interval=interval, reason=reason)

    if not mass or first_zero or (repeat_zero and n):
        return result(0.0, (0.0, 0.0))
    if repeat_zero:
        return result(float(mass * w), outward(mass * wi[0], mass * wi[1]))
    if negative:
        return result(reason="repeated_kernel_not_nonnegative")
    if wi == (0, 0):
        return result(0.0, (0.0, 0.0))
    lo, hi = map(Q, ratio["interval"])
    if hi >= 1:
        diverges = lo >= 1 and not wi[0] <= 0 <= wi[1]
        return result(reason="diverges" if diverges else "undetermined")
    r = Q(ratio["value"])
    # Nonnegative repeat kernel restricts the physical lower endpoint to zero.
    lo = max(0, lo)
    bounds = times(wi, (lo**n / (1 - lo), hi**n / (1 - hi)))
    return result(float(mass * w * r**n / (1 - r)), outward(mass * bounds[0], mass * bounds[1]))


def truncation(
    single, panels, geometry, requested, n, support, ratio, tol, negative, diagnostics, tail_model=None
):
    a = single.loss_window_eV[0]
    taus = np.r_[0.0, tau_grid(a)]
    m = (single.multiplier_interval_A or (single.multiplier_A,) * 2)[1]
    logk = log_laplace(panels, COLUMNS.index("total"), taus, m)
    logr = log_laplace(panels, 0, taus, m)
    mass = single.primary_abs_integral
    pmax = Q(geometry.points[-1]) if geometry.points else Q(0)
    available = logr < 0

    @lru_cache(maxsize=2)
    def tail_at_order(count):
        with np.errstate(divide="ignore", invalid="ignore", over="ignore", under="ignore"):
            if tail_model is not None:
                return tail_model(taus, count, logk, logr)
            log_tail = np.full(len(taus), np.inf)
            power = log_power(logr[available], count)
            log_tail[available] = log_up(log_up(logk[available] + power) + log_geometric_factor(logr[available]))
            return log_tail

    # Each grid value is an independent upper estimate. Work in log space
    # through the minimization, also when the final positive bound underflows.
    def bound(index, count):
        u = pmax - Q(geometry.edges[index])
        if not mass or u < (count + 1) * Q(a):
            return 0.0, "complete", None
        log_tail = tail_at_order(count)
        with np.errstate(invalid="ignore", over="ignore"):
            log_scale = log_up(math.nextafter(math.log(mass), math.inf) + log_up(taus * _up(u)))
            candidates = log_up(log_scale + log_tail)
        if np.all(np.isposinf(log_tail)):
            return None, "no_usable_tau", None
        # An overflowing exponent combined with a vanishing tail carries no
        # finite binary64 estimate; another grid point may still suffice.
        candidates = np.where(np.isnan(candidates), np.inf, candidates)
        i = int(np.argmin(candidates))
        if np.isneginf(log_tail[0]) and (np.isneginf(logk[0]) or (count and np.isneginf(logr[0]))):
            return 0.0, None, float(taus[i])
        value = exp_bound(float(candidates[i]))
        if value is None:
            return None, "not_representable", None
        return max(value, math.ulp(0.0)), None, float(taus[i])

    rows = [bound(i, n) for i in range(len(geometry.edges) - 1)]
    values, reasons, best = map(list, zip(*rows))
    estimate = None
    if tol is not None:
        # Cumulative bound at the lowest display edge controls every bin.
        lo, hi = 0, support
        while lo < hi:
            mid = (lo + hi) // 2
            value = bound(0, mid)[0]
            if value is not None and value <= tol * mass:
                hi = mid
            else:
                lo = mid + 1
        estimate = lo
    state = "complete_on_display" if n >= support else "truncated"
    ans = dict(
        requested=requested,
        computed=n,
        support_order=support,
        state=state,
        neglect_tolerance_fraction=tol,
        neglected_at_or_above_region_bound=values,
        bound_reason=reasons,
        tau_per_eV=best,
        required_order_estimate=estimate,
    )
    if tail_model is None:
        ans["neglected_total_area"] = tail_area(
            single, ratio, n, negative, np.isneginf(logk[0]), np.isneginf(logr[0])
        )
    if state == "truncated" and tol is not None:
        reason = (
            "truncation_bound_unavailable"
            if None in values
            else "truncation_bound_exceeds_tolerance"
            if max(values) > tol * mass
            else None
        )
        if reason:
            raise MultipleLossConvergenceError(
                reason,
                stage="truncation",
                required_order_estimate=estimate,
                diagnostics={**diagnostics, "truncation": ans},
            )
        ans["state"] = "bounded_on_display"
    return ans
