"""Quadrature on intersections of two Chebyshev panels (M4b design 6.1)."""

from bisect import bisect_left, bisect_right
from fractions import Fraction as Q
from itertools import islice, pairwise
import math

import numpy as np
from numpy.polynomial import chebyshev as cheb
from numpy.polynomial.legendre import leggauss

from ._loss_panels import CC, Panel, _nodes
from ._order_types import MultipleLossConvergenceError
from ._single_loss_types import EPS

GAUSS, WEIGHT = leggauss(17)


def mandatory_knots(features, window, order, limit):
    a, b = map(Q, window)
    atoms = {a, b}
    for key in ("jumps_eV", "kinks_eV"):
        atoms.update(Q(x) for x in features.get(key, ()) if a < x < b)
    sums = {Q(0)}
    for _ in range(order):
        new = set()
        for x in sums:
            new.update(x + y for y in atoms)
            if len(new) > limit + 1:
                raise MultipleLossConvergenceError(
                    "order_budget_not_met", order=order, stage="order_refinement", panel_limit=limit
                )
        sums = new
    edges = list(map(float, sorted(sums)))
    for pair in pairwise(edges):
        if not _nodes(*pair)[2]:
            raise MultipleLossConvergenceError(
                "order_mesh_unresolved", order=order, stage="order_refinement", panel_eV=pair
            )
    return edges


def edges_of(panels):
    return [panels[0].left_eV, *(p.right_eV for p in panels)]


def area_terms(panel):
    width = panel.right_eV - panel.left_eV
    weights = CC * (width / 2)
    return (
        panel.interpolation
        + width * panel.verification
        + panel.errors @ weights
        + 12 * EPS * (abs(panel.values) @ weights)
    )


def first_envelopes(panels):
    e = np.array([np.maximum(p.errors.max(axis=1), p.epsilon) + p.rho + p.coeff_diff for p in panels])
    return np.stack((e, np.zeros_like(e), np.zeros_like(e)), axis=1)


def affine_abs_mean(left, right):
    """Exact mean of the absolute value of the line between two endpoints."""
    a, b = abs(left), abs(right)
    result = (a + b) * 0.5
    crossing = (left < 0) != (right < 0)
    # Two triangles separated at the affine zero.
    return result - np.where(crossing, np.minimum(a, b) * (np.maximum(a, b) / np.where(a + b, a + b, 1)), 0.0)


class Polynomials:
    def __init__(self, panels, envelopes, columns):
        self.edges = edges_of(panels)
        self.ratios = [x.as_integer_ratio() for x in self.edges]
        self.coeff = np.stack([p.coefficients[columns] for p in panels])
        self.envelope = envelopes[:, :, columns]
        self.remainder = np.sum(abs(self.coeff[:, :, 2:]), axis=2)
        self.amplitude = np.sum(abs(self.coeff), axis=2)
        self.slope = abs(self.coeff) @ (np.arange(17) ** 2)
        self.sign_known = abs(self.coeff[:, :, 0]) >= abs(self.coeff[:, :, 1]) + self.remainder

    def abs_mean(self, indices, center, half, signed):
        coeff = self.coeff[indices]
        v = coeff[:, :, 0] + coeff[:, :, 1] * center[:, None]
        delta = coeff[:, :, 1] * half[:, None]
        upper = affine_abs_mean(v - delta, v + delta) + self.remainder[indices]
        # If the polynomial has one sign, its signed Gauss mean suffices.
        mean = np.where(self.sign_known[indices], abs(signed), upper)
        return mean + 64 * EPS * self.amplitude[indices]

    def norms(self):
        mean = np.einsum("pck,k->pc", self.coeff, np.array([0.0 if k % 2 else 1 / (1 - k * k) for k in range(17)]))
        n = len(self.coeff)
        absolute = self.abs_mean(np.arange(n), np.zeros(n), np.ones(n), mean)
        return np.sum(absolute * np.diff(self.edges)[:, None], axis=0)


def integer_edges(ratios, denominator):
    return [a * (denominator // b) for a, b in ratios]


def intersections(left, right, nodes):
    """Intersect in integer units of a common dyadic denominator, including tiny slivers."""
    ratios = [x.as_integer_ratio() for x in nodes]
    den = max(b for _, b in (*ratios, *left.ratios, *right.ratios))
    a, b = integer_edges(left.ratios, den), integer_edges(right.ratios, den)
    for k, s in enumerate(integer_edges(ratios, den)):
        t = max(b[0], s - a[-1])
        end = min(b[-1], s - a[0])
        i, j = bisect_left(a, s - t) - 1, bisect_right(b, t) - 1
        while t < end:
            stop = min(s - a[i], b[j + 1], end)
            width = stop - t
            da, db = a[i + 1] - a[i], b[j + 1] - b[j]
            # Normalized centers and half-widths rounded after exact subtraction.
            yield (
                k,
                i,
                j,
                width / den,
                (2 * s - t - stop - a[i] - a[i + 1]) / da,
                width / da,
                (t + stop - b[j] - b[j + 1]) / db,
                width / db,
            )
            t = stop
            if t == s - a[i]:
                i -= 1
            if t == b[j + 1]:
                j += 1


def _quadrature(left, right, nodes, charge):
    cells = intersections(left, right, nodes)
    out = np.zeros((4, len(nodes), left.coeff.shape[1]))
    # Batch independent 17-point rules. Every positive intersection is charged.
    while batch := list(islice(cells, min(256, charge(0) // 17 + 1))):
        k, i, j, width, x, hx, y, hy = map(np.asarray, zip(*batch))
        u = cheb.chebvander(x[:, None] - hx[:, None] * GAUSS, 16)
        v = cheb.chebvander(y[:, None] + hy[:, None] * GAUSS, 16)
        lv = left.coeff[i] @ u.transpose(0, 2, 1)
        rv = right.coeff[j] @ v.transpose(0, 2, 1)
        product = lv * rv
        mean = (product @ WEIGHT) * 0.5
        al = left.abs_mean(i, x, hx, (lv @ WEIGHT) * 0.5)
        ar = right.abs_mean(j, y, hy, (rv @ WEIGHT) * 0.5)
        el, er = left.envelope[i], right.envelope[j, 0]
        inh = (al + el.sum(axis=1)) * er + el[:, 0] * ar
        prop = (el[:, 1] + el[:, 2]) * ar
        # Coordinate, polynomial evaluation and quadrature roundoff, including
        # signed cancellation in an individual convolution.
        rounding = EPS * (
            128 * left.amplitude[i] * right.amplitude[j]
            + 4 * (left.slope[i] * right.amplitude[j] + left.amplitude[i] * right.slope[j])
        )
        for result, values in zip(out, (mean, inh, prop, rounding)):
            np.add.at(result, k, width[:, None] * values)
        charge(17 * len(batch))
    return out.transpose(0, 2, 1)


def build_panel(left, right, bounds, order, charge, inject, transform=None):
    nodes, check, separated = _nodes(*bounds)
    if not separated:
        raise MultipleLossConvergenceError(
            "order_mesh_unresolved", order=order, stage="order_refinement", panel_eV=bounds
        )
    result = _quadrature(left, right, np.concatenate((nodes, check)), charge)
    if transform is not None:
        transform(result, order)
    inter, ver = result[:, :, :17], result[:, :, 17:]
    if inject is not None:
        for block, role, positions in ((inter, "interpolation", nodes), (ver, "verification", check)):
            block[0], block[3] = inject(order, role, bounds, positions, block[0], block[3])
    p = Panel.from_lobatto(*bounds, inter[0], inter[3])
    p.check(ver[0], ver[3])
    env = result[1:].max(axis=2)
    env[2] += p.rho + p.coeff_diff
    return p, env, dict(zip(("inherited", "propagated", "local"), inter[1:]))


def minimum(panel, column):
    values = panel.values[column]
    if np.all(values == values[0]):
        return float(values[0])
    coeff = panel.coefficients[column]
    roots = cheb.chebroots(cheb.chebder(coeff))
    x = roots.real[(abs(roots.imag) < 64 * EPS) & (abs(roots.real) < 1)]
    return min(float(values[0]), float(values[-1]), *(cheb.chebval(x, coeff).tolist()))


def tau_grid(a):
    top = min(math.log10(np.finfo(float).max), math.log10(40.0) - math.log10(a))
    return np.array(
        [10 ** (-3 + i / 20) for i in range(max(0, math.floor((top + 3) * 20) + 1)) if -3 + i / 20 <= top]
    )
