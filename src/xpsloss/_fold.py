"""Integrate each detector region directly against loss-panel polynomials."""

from bisect import bisect_left, bisect_right
from dataclasses import dataclass
from fractions import Fraction as Q
from heapq import merge
from itertools import chain, pairwise
import math

import numpy as np
from numpy.polynomial.chebyshev import chebvander
from numpy.polynomial.legendre import leggauss

from ._single_loss_types import BinnedPrimary, COLUMNS, COMPONENTS, EPS, PARTS, QUANTITIES

GX, GW = leggauss(9)
LOBATTO = -np.cos(np.pi*np.arange(17)/16)
# Column q is the Chebyshev expansion of cardinal polynomial q.
CARDINAL = np.cos(np.pi*np.outer(np.arange(17), np.arange(16, -1, -1))/16)/8
CARDINAL[[0, -1], :] *= .5
CARDINAL[:, [0, -1]] *= .5
COMP_INDEX = [COLUMNS.index(c) for c in COMPONENTS]


def _difference(x, d):
    """TwoSum(x,-d), retained exactly even inside a sub-ulp cell."""
    hi = x-d
    v = hi-x
    lo = (x-(hi-v)) + (-d-v)
    return Q(hi) + Q(lo)


class Geometry:
    def __init__(self, primary, edges):
        # The public API validates a single positive measure. Internally this
        # also represents sums of signed delta and uniform-bin measures.
        items = primary if isinstance(primary, (tuple, list)) else (primary,)
        self.lines, self.bins = [], []
        points = set()
        for item in items:
            coefficient, p = item if isinstance(item, tuple) else (1, item)
            if coefficient == 0:
                continue
            if isinstance(p, BinnedPrimary):
                for a, b, f in zip(p.edges_eV[:-1], p.edges_eV[1:], p.intensities):
                    if f:
                        self.bins.append((Q(a), Q(b), Q(coefficient)*Q(f)/(Q(b)-Q(a))))
                        points.update((float(a), float(b)))
            else:
                for a, f in zip(p.energies_eV, p.intensities):
                    if f:
                        self.lines.append((Q(a), Q(coefficient)*Q(f)))
                        points.add(float(a))
        self.points = sorted(points)
        self.edges = tuple(map(float, edges))
        self.regions = [(None, Q(edges[0])), *zip(map(Q, edges[:-1]), map(Q, edges[1:])),
                        (Q(edges[-1]), None)]
        self.names = ["outside.low", *(f"bin[{i}]" for i in range(len(edges)-1)), "outside.high"]

    def breaks(self, region, a, b):
        def differences(d):
            lo = bisect_left(self.points, d+a)
            hi = bisect_right(self.points, d+b)
            for i in range(lo, hi):
                t = _difference(self.points[i], float(d))
                if a < t < b:
                    yield t
        # At most two ordered streams. Count coincident edges once without
        # allocating the Cartesian product of primary and detector edges.
        previous = None
        for t in merge(*(differences(d) for d in self.regions[region] if d is not None)):
            if t != previous:
                yield t
                previous = t

    def cell_count(self, bounds, limit):
        panel_edges = set(map(Q, bounds))
        base = (len(bounds)-1)*len(self.regions)
        if base > limit:
            return base
        for region in range(len(self.regions)):
            for t in self.breaks(region, Q(bounds[0]), Q(bounds[-1])):
                if t not in panel_edges:
                    base += 1
                    if base > limit:
                        return base
        return base

    def mass(self, region, t, side):
        """Exact G(t), including its one-sided limit at a delta crossing."""
        lo, hi = self.regions[region]
        lower = None if lo is None else lo+t
        upper = None if hi is None else hi+t
        # At T+0 the detected energy moves to the left.
        include_lower = side == "left" or (side == "point" and hi is not None)
        include_upper = side == "right" or (side == "point" and region == len(self.regions)-2)
        mass = Q(0)
        for p, f in self.lines:
            if (lower is None or p > lower or (include_lower and p == lower)) and (
                    upper is None or p < upper or (include_upper and p == upper)):
                mass += f
        for p, q, density in self.bins:
            start = p if lower is None else max(p, lower)
            stop = q if upper is None else min(q, upper)
            if start < stop:
                mass += (stop-start)*density
        return mass


def _scaled(a, b):
    """Scale a binary64 value by an exact rational, rounding only the result."""
    return float(Q(float(a))*b)


def _entry(value, parts):
    return dict(value=float(value), **dict(zip(PARTS, map(float, parts))), estimate=math.fsum(parts))


@dataclass
class Folded:
    outputs: dict
    panel_errors: dict
    support_widths: dict


def _cell(panel, left, right, g0, g1, length):
    """Nine Gauss points integrate a degree-16 polynomial times affine G."""
    width = right-left
    p0, p1 = Q(panel.left_eV), Q(panel.right_eV)
    middle = float((left+right-p0-p1)/(p1-p0))
    span = float(width/(p1-p0))
    basis = chebvander(middle + span*GX, 16) @ CARDINAL
    norm = max(abs(g0), abs(g1))
    start, stop = float(g0/norm), float(g1/norm)
    gauss_mass = (start*(1-GX) + stop*(1+GX))*(GW/4)
    products = gauss_mass[:, None]*basis
    weights = products.sum(axis=0)
    abs_products = np.abs(products).sum(axis=0)
    scale = Q(length)*width*norm
    values, errors = panel.values[COMP_INDEX], panel.errors[COMP_INDEX]
    point_error = errors @ abs(weights)
    # Rounding in the weights acts on the uncertain nodal values as well.
    uncertain_nodes = np.abs(values)+errors
    roundoff = 48*EPS*(uncertain_nodes @ abs_products)
    # Exact cell endpoints and exact affine overlaps leave only the rounding
    # of the normalized endpoint masses: integral |delta G| <= width *
    # (|delta G0| + |delta G1|)/2.
    envelope, slope = panel.component_bounds
    mass_error = (abs(g0-Q(start)*norm)+abs(g1-Q(stop)*norm))/2
    edge_scale = Q(length)*width*mass_error
    # |T_k'(x)| <= k^2 bounds the effect of rounding local Gauss coordinates.
    coordinate_error = (abs(Q(middle)-(left+right-p0-p1)/(p1-p0))
                        + abs(Q(span)-width/(p1-p0)) + Q(3*EPS)*(abs(Q(middle))+Q(span)))
    roundoff += float(coordinate_error)*slope
    vals = [_scaled(v, scale) for v in values @ weights]
    parts = []
    for j, i in enumerate(COMP_INDEX):
        interpolation = min(width*Q(panel.coeff_diff[i]), Q(panel.interpolation[i]))
        parts.append([float(Q(length)*norm*interpolation),
                      _scaled(panel.verification[i], scale), _scaled(point_error[j], scale),
                      _scaled(envelope[j], edge_scale),
                      _scaled(roundoff[j], scale)])
    return vals, parts


def fold(panels, primary, output_edges_eV, *, multiplier, hole_charge):
    """Fold a primary or a sequence of (signed coefficient, primary) pairs."""
    return _fold(panels, Geometry(primary, output_edges_eV), multiplier)


def _fold(panels, geometry, multiplier, lower_minus=None):
    outputs, panel_errors, support_widths = {}, {}, {}
    first = panels[0]
    endpoints = {"plus": (first.values[:, 0], first.errors[:, 0]),
                 "minus": lower_minus if lower_minus is not None else
                          (first.values[:, 0], first.errors[:, 0])}
    for region, name in enumerate(geometry.names):
        mass = geometry.mass(region, Q(0), "point")
        zero = float(mass)
        outputs[name+".zero_loss"] = _entry(zero, [0, 0, 0, 0, float(abs(mass-Q(zero)))])
        sums = []
        per_panel, support = [], []
        for p in panels:
            a, b = Q(p.left_eV), Q(p.right_eV)
            cuts = chain((a,), geometry.breaks(region, a, b), (b,))
            cells, errors = [], []
            active_width = Q(0)
            for left, right in pairwise(cuts):
                g0 = geometry.mass(region, left, "right")
                g1 = geometry.mass(region, right, "left")
                if g0 or g1:
                    active_width += right-left
                    v, e = _cell(p, left, right, g0, g1, multiplier)
                    cells.append(v)
                    errors.append(e)
            sums.append([math.fsum(c[j] for c in cells) for j in range(3)])
            per_panel.append(np.sum(errors, axis=0) if errors else np.zeros((3, 5)))
            support.append(float(active_width))
        errors = np.array(per_panel)
        widths = np.asarray(support)
        for j, component in enumerate(COMPONENTS):
            key = name+".single_loss_"+component
            outputs[key] = _entry(math.fsum(s[j] for s in sums), errors[:, j].sum(axis=0))
            panel_errors[key] = errors[:, j].sum(axis=1)
            support_widths[key] = widths
        satellite = outputs[name+".single_loss_total"]
        spectrum = zero+satellite["value"]
        parts = [satellite[k] for k in PARTS]
        parts[-1] += outputs[name+".zero_loss"]["estimate"] + EPS*abs(spectrum)
        outputs[name+".spectrum"] = _entry(spectrum, parts)

        for side, direction in (("plus", "right"), ("minus", "left")):
            g = geometry.mass(region, Q(first.left_eV), direction)
            v, e = endpoints[side]
            factor = Q(multiplier)*Q(first.left_eV)
            # Exact overlap and product, including uncertainty in the rounded G.
            eg = abs(g-Q(float(g)))
            for c, i in zip(COMPONENTS, COMP_INDEX):
                value = _scaled(v[i], factor*g)
                err = float(factor*(abs(g)*Q(e[i]) + eg*(abs(Q(v[i]))+Q(e[i]))))
                outputs[f"{name}.cutoff_sensitivity.{side}.{c}"] = _entry(
                    value, [0, 0, err, 0, 3*EPS*abs(value)])
    return Folded(outputs, panel_errors, support_widths)

