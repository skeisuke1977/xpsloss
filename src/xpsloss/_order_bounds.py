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
            w = right - left
            if height == 0 or w == 0:
                continue
            log_weights = np.full(len(taus), math.log(w))
            nz = taus > 0
            z = taus[nz] * w
            small = z < EPS
            log_weights[nz] = np.log(-np.expm1(-np.where(small, 1.0, z))) - np.log(taus[nz])
            log_weights[np.flatnonzero(nz)[small]] = math.log(w)
            total = np.logaddexp(total, math.log(height) + log_weights - taus * left)
        error = area_terms(p)[col]
        if error:
            total = np.logaddexp(total, math.log(error) - taus * a)
    return total + math.log(m) + math.log1p(64 * EPS)


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
            power = float(count) * logr[available] if count else 0.0
            log_tail[available] = logk[available] + power - np.log(-np.expm1(logr[available]))
            return log_tail

    # Each grid value is an independent upper estimate. Work in log space
    # through the minimization, also when the final positive bound underflows.
    def bound(index, count):
        u = pmax - Q(geometry.edges[index])
        if not mass or u < (count + 1) * Q(a):
            return 0.0, "complete", None
        log_tail = tail_at_order(count)
        with np.errstate(invalid="ignore", over="ignore"):
            candidates = math.log(mass) + taus * float(u) + log_tail
        if np.all(np.isposinf(log_tail)):
            return None, "no_usable_tau", None
        # An overflowing exponent combined with a vanishing tail carries no
        # finite binary64 estimate; another grid point may still suffice.
        candidates = np.where(np.isnan(candidates), np.inf, candidates)
        i = int(np.argmin(candidates))
        if candidates[i] > math.log(np.finfo(float).max):
            return None, "not_representable", None
        if np.isneginf(log_tail[0]) and (np.isneginf(logk[0]) or (count and np.isneginf(logr[0]))):
            return 0.0, None, float(taus[i])
        value = math.nextafter(math.exp(float(candidates[i])), math.inf)
        if not math.isfinite(value):
            return None, "not_representable", None
        return value, None, float(taus[i])

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
