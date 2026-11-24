"""Simonsen, Yubero & Tougaard (1997), Eqs. (13)–(17), at normal emission."""

from dataclasses import dataclass
from decimal import Context, Decimal
from fractions import Fraction
import math

import numpy as np
from scipy.optimize import brentq
from scipy.special import roots_legendre

from .dielectric import _checked_static_response
from .model import (
    BOHR_RADIUS_A, DECIMAL_DIGITS, DEFAULT_REL, DETUNING_REL, FLANK_RATIO, FLOAT_EPSILON, HBAR2_OVER_2ME_EV_A2, NEARBY_FACTOR,
    RECIPROCAL_ROUNDOFF,
    Material, NumericalConvergenceError, _PanelLimitError, _accuracy, _exact_oscillator_sum, _first_loss, _loss_point,
    _oscillator_response, finite_number,
)


REFLECTION_ROUNDOFF = 8 * FLOAT_EPSILON  # Complex oscillator sum followed by subtraction from 2.
PANEL_ROUNDOFF = 32 * FLOAT_EPSILON  # Floor on the embedded 32-node quadrature rule.
LOG_K_PANEL_SPAN = 1  # Bound the ratio of k endpoints to e before adaptive refinement.
MOMENTUM_ULPS = 4  # exp(log k), squaring and the offset in k_parallel².
FIXED_ABS = 1e-12  # Default absolute quadrature error in inverse angstrom per eV.
FIXED_REL = 1e-9  # Relative target for each returned cross section.

FIXED_DEPTH_MODEL_ID = "simonsen1997_fixed_depth_v0"
EPSILON_MOMENTA = ("k_parallel", "k_total")
TERM_NAMES = ("term_13a", "term_13b", "term_14", "term_15_T1", "term_15_T2",
              "term_15_T3", "term_16_e", "term_16_h")
_LOW_RULE = roots_legendre(16)
_HIGH_RULE = roots_legendre(32)
_NODES = np.concatenate((_LOW_RULE[0], _HIGH_RULE[0]))
_ELECTRON = [0, 1, 3, 5, 6]
_HOLE = [2, 4, 7]
LOSS_BATCH = 8  # Share array evaluations of the initial rule across nearby loss points.
_MIN_NORMAL = np.finfo(float).tiny


class FixedDepthConvergenceError(NumericalConvergenceError):
    """Quadrature could not resolve the requested absolute or relative accuracy."""


@dataclass(frozen=True)
class FixedDepthResult:
    """All K arrays are in eV^-1 Å^-1 per solid path, including the vacuum terms."""

    electron_energy_eV: float
    depth_A: float
    path_length_A: float
    theta_deg: float
    hole_charge: float
    epsilon_momentum: str
    loss_eV: np.ndarray
    term_13a: np.ndarray
    term_13b: np.ndarray
    term_14: np.ndarray
    term_15_T1: np.ndarray
    term_15_T2: np.ndarray
    term_15_T3: np.ndarray
    term_16_e: np.ndarray
    term_16_h: np.ndarray
    no_hole: np.ndarray
    hole: np.ndarray
    total: np.ndarray
    quadrature_error: np.ndarray


def _parallel_limits(energy, loss):
    """Eq. (17), factoring k±² − Ω0² before evaluating either square root."""
    s, r = math.sqrt(energy), np.sqrt(energy - loss)
    plus = s + r
    minus = loss / plus
    scale = 2 * s * math.sqrt(HBAR2_OVER_2ME_EV_A2)
    return (minus * np.sqrt(minus * (3*s + r)) / scale,
            plus * np.sqrt(plus * (3*s - r)) / scale)


def _oscillator_sum(material, loss, q):
    """E4b: compute S = 1 − 1/ε directly, retaining weak responses near ε = 1."""
    response = np.zeros_like(q, dtype=complex)
    real_error = np.zeros_like(q, dtype=float)
    imag_error = np.zeros_like(q, dtype=float)
    rounding = REFLECTION_ROUNDOFF * len(material.oscillators)
    t = np.broadcast_to(loss, q.shape)
    for osc in material.oscillators:
        if osc.strength_eV2:
            term, detuning_error = _oscillator_response(osc, t, q, with_detuning_error=True)
            response += term
            # For a perturbed detuning, |dS| <= |S| |dD/D| to first
            # order, whereas Im S is proportional to 1/|D|^2. The line
            # evaluator reconstructs D before its relative error reaches 1e-10.
            real_error += rounding * abs(term.real) + 2 * detuning_error * abs(term)
            imag_error += (rounding + 3 * detuning_error) * term.imag
    return response, real_error, imag_error


def _bounded_sum(material, loss, q):
    """Fast E4b with a uniform bound on the detuning error over these panels."""
    lo, hi = float(np.min(loss)), float(np.max(loss))
    qmax = float(np.max(q))
    q2 = q*q
    response = np.zeros(q.shape, dtype=complex)
    magnitude = np.zeros(q.shape)
    largest = 0.0
    for osc in material.oscillators:
        if not osc.strength_eV2:
            continue
        w, width = osc.energy_eV, osc.width_eV
        slope = osc.dispersion * HBAR2_OVER_2ME_EV_A2
        # p=T-W>0, shift<=2p: (W+T+shift)<=3T and |D|>=width*T.
        # For shift>2p use |D|>=(shift-p)*(W+T+shift). For T<=W
        # the detuning has no subtraction of the dispersion from p.
        relative = (RECIPROCAL_ROUNDOFF * max(1., 3*(hi-w)/width)
                    if slope and hi > w else 0.0)
        center_max = w + slope * qmax * qmax
        scale = max(center_max*center_max, hi*hi, width*hi)
        upper = 2 * scale * scale
        lower = (width * lo) * (width * lo)
        # Certify the range of the squared denominator and its reciprocal
        # before using real arithmetic instead of scaled complex division.
        if (relative >= DETUNING_REL or lower < _MIN_NORMAL or not math.isfinite(upper)
                or osc.strength_eV2 / upper < _MIN_NORMAL
                or not math.isfinite(osc.strength_eV2 / lower)):
            return _oscillator_sum(material, loss, q)
        shift = slope * q2
        delta = ((w - loss) + shift) * ((w + shift) + loss)
        damping = width * loss
        weight = osc.strength_eV2 / (delta*delta + damping*damping)
        real = weight * delta
        response.real += real
        response.imag += weight * damping
        magnitude += abs(real)
        largest = max(largest, relative)
    rounding = REFLECTION_ROUNDOFF * len(material.oscillators)
    # Sum of magnitudes is bounded by sum |Re S_i| + Im S (positive damping).
    sr = rounding * magnitude + (2*largest) * (magnitude + response.imag)
    si = (rounding + 3*largest) * response.imag
    return response, sr, si


def _node_squares(left, right, points, loss, energy, total_momentum):
    """Reconstruct exp(2 log k) at the selected GL nodes of the actual panels."""
    offset = 0
    if total_momentum:
        offset = Fraction(float(loss))**2 / (4 * Fraction(HBAR2_OVER_2ME_EV_A2) * Fraction(energy))
    ctx = Context(prec=DECIMAL_DIGITS)
    for i in np.flatnonzero(points.any(axis=1)):
        lo, hi = Decimal(float(left[i])), Decimal(float(right[i]))
        center, span = ctx.add(lo, hi), ctx.subtract(hi, lo)
        for j in np.flatnonzero(points[i]):
            log_z = ctx.fma(span, Decimal(float(_NODES[j])), center)
            yield Fraction(ctx.exp(log_z)) + offset


def _exact_points(material, loss, q, points, coefficients, *, node_squares=None):
    """Replace the selected Cartesian coefficients and their rounding bounds."""
    (s, r, c), errors = coefficients
    losses = np.broadcast_to(loss, q.shape)
    for row in np.argwhere(points):
        i = tuple(row)
        z = next(node_squares) if node_squares is not None else None
        x, y = _exact_oscillator_sum(material, losses[i], q[i], q_squared=z)
        # Put Re S and Im S over one integer denominator. R and C are
        # needed only as floats: leave their numerator/denominator pairs
        # unreduced, and round once after the final subtraction in C.
        den = math.lcm(x.denominator, y.denominator)
        a = x.numerator * (den // x.denominator)
        b = y.numerator * (den // y.denominator)
        complement = 2 * den - a
        norm = complement**2 + b**2
        # With positive loss and width, Im S > 0 for any active line.
        # With no active lines S = 0. Retain a defensive pole check.
        if not norm:
            raise FixedDepthConvergenceError(f"surface pole at loss_eV={losses[i]:g}")
        rx, ry = a*complement - b**2, 2*b*den
        common = den * norm
        s[i], r[i] = complex(x, y), complex(rx / norm, ry / norm)
        c[i] = complex((2*a*norm - rx*den) / common, (2*b*norm - ry*den) / common)
    for value, (er, ei) in zip((s, r, c), errors):
        er[points] = FLOAT_EPSILON * abs(value.real[points])
        ei[points] = FLOAT_EPSILON * abs(value.imag[points])


def _coefficients(material, loss, q, epsrel=None, *, bounded=False):
    """Return S, R and C and absolute errors of their Cartesian components."""
    if bounded:
        s, sr, si = _bounded_sum(material, loss, q)
    else:
        s, sr, si = _oscillator_sum(material, loss, q)
    rounding = REFLECTION_ROUNDOFF * len(material.oscillators)
    d = abs(2 - s)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        r = s / (2 - s)
        # Im R = 2 Im S / |2-S|^2, including at a large bulk response.
        r.imag = (s.imag / d) * (2 / d)
        c = 2 * s - r
        relative = (sr + si) / d
        # R = 2/(2-S) - 1. Bound the input perturbation and rounding
        # of the real projection of the complex division separately.
        rr = (2 / d) * relative / (1 - relative) + rounding * abs(r)
        ri = ((2 / d) * (si / d) + abs(r.imag) * relative * (2 + relative)) / (1 - relative)**2
        ri += rounding * abs(r.imag)
        cr = 2*sr + rr + rounding*abs(c.real)
        ci = 2*si + ri + rounding*abs(c.imag)
        uncertain = ~((relative < 1) & np.isfinite(cr) & np.isfinite(ci))
        if epsrel is not None:
            # Coefficient-only tests supply epsrel here. Integration uses its weighted budgets.
            for value, error in ((s.real, sr), (s.imag, si), (r.real, rr), (r.imag, ri), (c.real, cr), (c.imag, ci)):
                uncertain |= ~(error <= epsrel * abs(value))
    coefficients = (s, r, c), ((sr, si), (rr, ri), (cr, ci))
    if np.any(uncertain):
        if bounded:
            return _coefficients(material, loss, q, epsrel)
        _exact_points(material, loss, q, uncertain, coefficients)
    return coefficients


def _include_node_roundoff(material, loss, q, relative_z, coefficients):
    """Propagate uncertainty in z=q², even after exact evaluation at rounded q."""
    q2 = q*q
    ds = np.zeros_like(q)
    for osc in material.oscillators:
        if not (osc.dispersion and osc.strength_eV2):
            continue
        shift = osc.dispersion * HBAR2_OVER_2ME_EV_A2 * q2
        center = osc.energy_eV + shift
        dshift = relative_z * shift
        delta = ((osc.energy_eV - loss) + shift) * (center + loss)
        # For D=(W+shift)^2-T^2-i*width*T, this bounds the finite
        # change in D. Allow the same shift uncertainty when estimating
        # the distance to the pole with a binary64 detuning.
        change = dshift * (2*center + dshift)
        distance = np.maximum(osc.width_eV*loss, abs(delta) - change) * (1 - REFLECTION_ROUNDOFF)
        fraction = change / distance
        if np.any(fraction >= 1):
            raise FixedDepthConvergenceError("oscillator denominator unresolved at quadrature node")
        # Bound the derivative over the entire coordinate interval. This
        # same bound can be rescaled for a more accurately located node.
        ds += (osc.strength_eV2 / (distance - change)) * (fraction / (1 - fraction))

    (s, _, _), (se, re, ce) = coefficients
    # R=2/(2-S)-1. Allow for the uncertainty already present at fixed q
    # in both denominators of the finite reciprocal difference.
    distance = abs(2 - s) - (se[0] + se[1])
    if np.any(ds >= distance):
        raise FixedDepthConvergenceError("surface denominator unresolved at quadrature node")
    dr = (2 / (distance - ds)) * (ds / (distance - ds))
    for (er, ei), change in zip((se, re, ce), (ds, dr, 2*ds + dr)):
        er += change
        ei += change
    return ds, dr


def _term_integrands(k, coefficients, energy, depth, omega, *, select_bounds=False):
    """Eight real densities in d ln(k_parallel), in the order of M1 §4.2."""
    (response, reflection, surface), (ds, dr, dc) = coefficients
    bulk = -response

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

    measure = k**2 / (math.pi * BOHR_RADIUS_A * energy)
    pref = measure / (2 * depth * omega)
    w13a = measure / denom
    term_13a = w13a * response.imag
    # The A and B terms of 13b form a conjugate pair. C supplies |1-u|².
    w13b = 2 * (plus**2 * diff).imag
    term_13b = pref * (bulk.imag * w13b + bulk.real * norm)
    weight_14 = plus*diff + minus*mixed
    term_14 = pref * (bulk * weight_14).real / omega
    term_15_T1 = pref * surface.real * norm
    w15b = minus * mixed
    w15c = minus**2 * diff.conj()
    w16e = plus * (2j*k/denom - plus*u)
    w16h = plus * u
    term_15_T2 = pref * (surface * w15b).real / omega
    term_15_T3 = pref * (reflection * w15c).real
    term_16_e = -pref * (reflection * w16e).real
    term_16_h = -pref * (reflection * w16h).real / omega
    terms = np.stack((term_13a, term_13b, term_14, term_15_T1, term_15_T2,
                      term_15_T3, term_16_e, term_16_h), axis=-2)

    def projected(weight, coefficient_size, uncertainty):
        bound = abs(weight.real)*uncertainty[0] + abs(weight.imag)*uncertainty[1]
        if select_bounds:
            return bound
        return bound + abs(weight)*coefficient_size

    bulk_size = PANEL_ROUNDOFF * abs(bulk)
    reflection_size = PANEL_ROUNDOFF * abs(reflection)
    surface_size = PANEL_ROUNDOFF * abs(surface)
    error13b = ds[1] * abs(w13b) + ds[0] * norm
    if not select_bounds:
        error13b += PANEL_ROUNDOFF * (2 * abs(bulk.imag) * abs(diff) / denom + abs(bulk.real) * norm)

    errors = np.stack((
        w13a * ds[1],
        pref * error13b,
        pref * projected(weight_14, bulk_size, ds) / omega,
        pref * dc[0] * norm,
        pref * projected(w15b, surface_size, dc) / omega,
        pref * projected(w15c, reflection_size, dr),
        pref * projected(w16e, reflection_size, dr),
        pref * projected(w16h, reflection_size, dr) / omega,
    ), axis=-2)
    return terms, errors


def _line_positions(material, loss, q_offset, upper_squared):
    """Centers and half-width scales in k_parallel², used only for panel layout."""
    dispersive = []
    fixed = 0j
    for index, osc in enumerate(material.oscillators):
        if not osc.strength_eV2:
            continue
        slope = osc.dispersion * HBAR2_OVER_2ME_EV_A2
        if slope:
            pole = (loss - osc.energy_eV) / slope - q_offset**2
            dispersive.append((osc, slope, pole))
            yield pole, osc.width_eV / (2 * slope), f"bulk oscillator {index}", loss / slope
        else:
            fixed += osc.strength_eV2 / complex((osc.energy_eV-loss)*(osc.energy_eV+loss),
                                                -osc.width_eV*loss)
    if not dispersive:
        return
    if len(dispersive) == 1 and fixed == 0:
        # For one oscillator the surface denominator is 2(Ω²-T²-iΓT) − A.
        osc, slope, _ = dispersive[0]
        center = math.hypot(loss, math.sqrt(osc.strength_eV2 / 2))
        yield ((center - osc.energy_eV) / slope - q_offset**2,
               osc.width_eV * loss / (2 * center * slope), "surface oscillator sum", center / slope)
        return

    # Surface lines belong to the oscillator *sum*. In the narrow-line limit
    # Re(S) decreases between successive bulk poles in z=k_parallel². Locate
    # its crossings of 2 and linearize the damping there to seed both flanks.
    poles = {pole for _, _, pole in dispersive}
    stop = NEARBY_FACTOR * upper_squared
    edges = sorted({0.0, stop, *(p for p in poles if 0 < p < stop)})
    for left, right in zip(edges, edges[1:]):
        def secular(z):
            if z == left and z in poles:
                return math.inf
            if z == right and z in poles:
                return -math.inf
            value = fixed.real - 2
            for osc, slope, pole in dispersive:
                center = osc.energy_eV + slope * (z + q_offset**2)
                value += osc.strength_eV2 / (slope * (z - pole) * (center + loss))
            return value

        if secular(left) < 0 or secular(right) > 0:
            continue
        try:
            z = brentq(secular, left, right, xtol=np.finfo(float).tiny)
        except RuntimeError as exc:
            raise FixedDepthConvergenceError(f"surface oscillator sum unresolved at loss_eV={loss:g}") from exc
        if z in poles:
            raise FixedDepthConvergenceError(f"unresolved oscillator and surface lines at loss_eV={loss:g}")
        damping, gradient = fixed.imag, 0.0
        for osc, slope, pole in dispersive:
            center = osc.energy_eV + slope * (z + q_offset**2)
            detuning = slope * (z - pole) * (center + loss)
            weight = osc.strength_eV2 / detuning**2
            damping += weight * osc.width_eV * loss
            gradient += weight * 2 * center * slope
        yield z, damping / gradient, "surface oscillator sum", z + q_offset**2


def _log_breakpoints(material, loss, lower, upper, q_offset):
    """Log panels, split at the dispersive oscillator and surface lines."""
    lo, hi = math.log(lower), math.log(upper)
    points = list(np.linspace(lo, hi, math.ceil((hi - lo) / LOG_K_PANEL_SPAN) + 1))
    lower_squared, upper_squared = lower**2, upper**2
    for center, width, line, scale in _line_positions(material, loss, q_offset, upper_squared):
        distance = max(lower_squared - center, center - upper_squared, 0.0)
        # Propagate coordinate roundoff through exp(log k), k² and the
        # dispersed energy. Coincident quadrature nodes from both rules can
        # otherwise give a spuriously tiny difference for an unresolved line.
        spacing = math.ulp(max(abs(center), q_offset**2, scale))
        if MOMENTUM_ULPS * spacing > DEFAULT_REL * max(width, distance):
            raise FixedDepthConvergenceError(
                f"fixed-depth {line}, loss_eV={loss:g}: linewidth unresolved on momentum panels")
        if lower_squared < center < upper_squared:
            points.append(math.log(center) / 2)
        # Geometric flanks resolve the tails even on a much wider log panel.
        offset = width
        while 0 < offset < FLANK_RATIO * max(abs(center - lower_squared), abs(upper_squared - center)):
            for flank in (center - offset, center + offset):
                if lower_squared < flank < upper_squared:
                    points.append(math.log(flank) / 2)
            offset *= FLANK_RATIO
    return np.unique(points)


def _panel_integrals(edges_left, edges_right, integrand):
    mid = (edges_left + edges_right) / 2
    half = (edges_right - edges_left) / 2
    # Roundoff in the affine log node becomes relative error after exp.
    # Include squaring, the optional normal momentum offset and the shift.
    log_size = np.maximum(abs(edges_left), abs(edges_right))
    relative_z = MOMENTUM_ULPS * FLOAT_EPSILON * (2 + log_size[:, None])
    values, errors = integrand(np.exp(mid[:, None] + half[:, None] * _NODES), relative_z, edges_left, edges_right)
    split = len(_LOW_RULE[0])

    def weighted(array, nodes, weights):
        return half[:, None] * np.einsum("ptn,n->pt", array[..., nodes], weights)

    low = weighted(values, slice(None, split), _LOW_RULE[1])
    high = weighted(values, slice(split, None), _HIGH_RULE[1])
    magnitude = weighted(abs(values), slice(split, None), _HIGH_RULE[1])
    high_error = weighted(errors, slice(split, None), _HIGH_RULE[1])
    low_error = weighted(errors, slice(None, split), _LOW_RULE[1])
    # Rounding affects both sides of the embedded difference as well as
    # the returned high-rule value. The weights of both rules are positive.
    quadrature = np.maximum(abs(high - low), PANEL_ROUNDOFF * magnitude)
    arithmetic = 2 * high_error + low_error
    return high, quadrature + arithmetic, arithmetic


def _component_sums(terms, charge_weight):
    """Combine signed terms with charge, or their absolute errors with a positive charge weight."""
    electron = terms[..., _ELECTRON].sum(axis=-1)
    hole = terms[..., _HOLE].sum(axis=-1)
    return np.stack((electron, hole, electron + charge_weight * hole), axis=-1)


def _initial_edges(material, loss, lower, upper, q_offset, limit):
    edges = _log_breakpoints(material, loss, lower, upper, q_offset)
    if len(edges) - 1 > limit:
        raise _PanelLimitError(f"integration limit too small for {len(edges) - 1} initial oscillator panels")
    return edges


def _initial_panels(material, losses, energy, depth, omega, offsets, lower, upper, limit):
    """Evaluate both rules for a batch of loss grids; each grid keeps its own panel boundaries."""
    edges = [_initial_edges(material, t, lo, hi, off, limit)
             for t, lo, hi, off in zip(losses, lower, upper, offsets)]
    counts = [len(e) - 1 for e in edges]
    owners = np.repeat(np.arange(len(edges)), counts)
    loss = losses[owners, None]
    frequency = omega[owners, None]
    offset = offsets[owners, None]
    dispersive = any(o.dispersion and o.strength_eV2 for o in material.oscillators)

    def integrand(k, relative_z, _left, _right):
        q = np.hypot(k, offset) if offset.any() else k
        coefficients = _coefficients(material, loss, q if dispersive else np.zeros_like(loss), bounded=True)
        if dispersive:
            _include_node_roundoff(material, loss, q, relative_z, coefficients)
        return _term_integrands(k, coefficients, energy, depth, frequency)

    panels = _panel_integrals(np.concatenate([e[:-1] for e in edges]),
                              np.concatenate([e[1:] for e in edges]), integrand)
    starts = np.cumsum([0, *counts])
    return [(e, tuple(a[start:stop] for a in panels)) for e, start, stop in zip(edges, starts, starts[1:])]


def _integrate_loss(material, loss, energy, depth, omega, q_offset, lower, upper,
                    charge, epsabs, epsrel, limit, initial=None):
    edges = initial[0] if initial is not None else _initial_edges(material, loss, lower, upper, q_offset, limit)
    # Refining one panel must not move either endpoint of its neighbor.
    left, right = edges[:-1].copy(), edges[1:].copy()

    constant = None
    constant_exact = False
    density_budget = None
    def integrand(k, relative_z, panel_left, panel_right):
        nonlocal constant, constant_exact
        if constant is None and not any(o.dispersion and o.strength_eV2 for o in material.oscillators):
            constant = _coefficients(material, loss, np.zeros(1))
        q = np.hypot(k, q_offset) if q_offset else k
        coefficients = constant or _coefficients(material, loss, q)
        if constant is None:
            ds, dr = _include_node_roundoff(material, loss, q, relative_z, coefficients)
        if density_budget is not None and not constant_exact:
            _, bounds = _term_integrands(k, coefficients, energy, depth, omega, select_bounds=True)
            points = (bounds > density_budget[:, None]).any(axis=-2)
            if np.any(points):
                if constant is None:
                    squares = _node_squares(panel_left, panel_right, points, loss, energy, bool(q_offset))
                    _exact_points(material, loss, q, points, coefficients, node_squares=squares)
                    # Rounding the GL abscissa now affects log(q²) only
                    # through the panel width. Retain this and the residual
                    # of the decimal exponential after replacing S, R and C.
                    residual = 4*FLOAT_EPSILON*(panel_right-panel_left) + 10.0**(4-DECIMAL_DIGITS)
                    scale = (residual / relative_z[:, 0])[np.nonzero(points)[0]]
                    for (er, ei), bound in zip(coefficients[1], (ds, dr, 2*ds + dr)):
                        er[points] += scale * bound[points]
                        ei[points] += scale * bound[points]
                else:
                    _exact_points(material, loss, np.zeros(1), np.ones(1, dtype=bool), constant)
                    constant_exact = True
        return _term_integrands(k, coefficients, energy, depth, omega)

    # The same estimate covers individual terms, both unit-charge components
    # and their user-weighted sum, even when signed contributions cancel.
    charge_weight = max(1, abs(charge))
    values, errors, arithmetic = initial[1] if initial is not None else _panel_integrals(left, right, integrand)
    while True:
        terms = values.sum(axis=0)
        term_error = errors.sum(axis=0)
        term_budget = _accuracy(terms, epsabs, epsrel)
        component_budget = _accuracy(_component_sums(terms, charge), epsabs, epsrel)
        component_error = _component_sums(term_error, charge_weight)
        if np.all(term_error <= term_budget) and np.all(component_error <= component_budget):
            return terms, component_error[-1]
        if initial is not None:
            # Tighten the oscillator-specific bounds before choosing any
            # expensive point recomputations or subdividing the batch's mesh.
            initial = None
            values, errors, arithmetic = _panel_integrals(left, right, integrand)
            continue
        if density_budget is None:
            rounding = arithmetic.sum(axis=0)
            budgets = np.concatenate((term_budget, component_budget))
            estimates = np.concatenate((term_error, component_error))
            rounding = np.concatenate((rounding, _component_sums(rounding, charge_weight)))
            quadrature = np.maximum(estimates - rounding, 0.0)
            affected = (estimates > budgets) & ((rounding >= budgets) | (quadrature < budgets))
            if np.any(affected):
                # Allocate the remaining budget after the embedded estimate to
                # terms of affected outputs only. If quadrature also exceeds
                # its budget, its panels will need refinement independently.
                residual = np.where(quadrature < budgets, budgets - quadrature, budgets)
                allowance = np.where(affected, residual, np.inf)
                remaining, component_left = allowance[:8], allowance[8:]
                remaining[_ELECTRON] = np.minimum(remaining[_ELECTRON], component_left[0] / len(_ELECTRON))
                remaining[_HOLE] = np.minimum(remaining[_HOLE], component_left[1] / len(_HOLE))
                total_share = component_left[2] / (len(_ELECTRON) + charge_weight * len(_HOLE))
                density_budget = np.minimum(remaining, total_share) / (3 * (edges[-1] - edges[0]))
                values, errors, arithmetic = _panel_integrals(left, right, integrand)
                continue
        if len(left) >= limit:
            raise FixedDepthConvergenceError(
                f"fixed-depth oscillator sum quadrature failed to converge at loss_eV={loss:g} "
                f"with {limit} panels (estimated error {component_error[-1]:.3g})")
        # Refine where the largest fraction of a requested budget is spent.
        priority = np.max(errors / term_budget, axis=1)
        priority = np.maximum(priority, np.max(_component_sums(errors, charge_weight) / component_budget, axis=1))
        i = priority.argmax()
        midpoint = (left[i] + right[i]) / 2
        children_left = np.array([left[i], midpoint])
        children_right = np.array([midpoint, right[i]])
        child_values, child_errors, child_arithmetic = _panel_integrals(children_left, children_right, integrand)
        right[i], values[i], errors[i] = midpoint, child_values[0], child_errors[0]
        left = np.append(left, midpoint)
        right = np.append(right, children_right[1])
        values = np.vstack((values, child_values[1]))
        errors = np.vstack((errors, child_errors[1]))
        arithmetic[i] = child_arithmetic[0]
        arithmetic = np.vstack((arithmetic, child_arithmetic[1]))


def effective_cross_section_fixed_depth(
    material: Material, electron_energy_eV, loss_eV, depth_A, *, theta_deg=0.0,
    hole_charge=1.0, epsilon_momentum, epsabs=FIXED_ABS, epsrel=FIXED_REL, limit=200,
) -> FixedDepthResult:
    """Fixed-depth K for a scalar or 1-D loss grid, returned as 1-D arrays.

    θ is measured from the surface normal. All hole terms refer to unit charge.
    epsabs is in inverse angstrom per eV, epsrel applies to each term and
    component, and limit caps panels per loss. quadrature_error includes
    coefficient rounding and excludes material and model uncertainties.
    Every nonzero line must be resolved on the momentum grid. Strength alone
    cannot bound its peak or its effect on the coupled surface denominator.
    """
    if not isinstance(material, Material):
        raise ValueError("material must be a Material")
    if material.gap_eV > 0:
        raise ValueError("fixed-depth gapped materials are deferred to M2-2")
    if not isinstance(epsilon_momentum, str) or epsilon_momentum not in EPSILON_MOMENTA:
        raise ValueError("epsilon_momentum must be 'k_parallel' or 'k_total'")
    theta = finite_number(theta_deg, "theta_deg")
    if theta != 0:
        raise ValueError("only normal emission (theta_deg=0) is supported; see Q3/Q15")
    energy = finite_number(electron_energy_eV, "electron_energy_eV")
    depth = finite_number(depth_A, "depth_A")
    charge = finite_number(hole_charge, "hole_charge")
    if depth <= 0 or energy <= 0:
        raise ValueError("depth_A and electron_energy_eV must be positive")
    grid = np.atleast_1d(loss_eV)
    if grid.dtype.kind not in "fiu" or grid.ndim != 1 or not grid.size:
        raise ValueError("loss_eV must be a real scalar or a nonempty one-dimensional grid")
    grid = np.array(grid, dtype=float)
    if not np.all((grid > 0) & (grid < energy)):
        raise ValueError("require 0 < loss_eV < electron_energy_eV; the zero-loss delta term is excluded")
    epsabs, epsrel = finite_number(epsabs, "epsabs"), finite_number(epsrel, "epsrel")
    if epsabs <= 0 or epsrel <= 0:
        raise ValueError("epsabs and epsrel must be positive")
    if isinstance(limit, bool) or not isinstance(limit, int) or limit < 1:
        raise ValueError("limit must be a positive integer")

    with np.errstate(over="ignore", divide="ignore", invalid="ignore", under="ignore"):
        lower, upper = _parallel_limits(energy, grid)
        omega = grid / (2 * math.sqrt(HBAR2_OVER_2ME_EV_A2 * energy))
        bad = ~((lower > 0) & (lower < upper) & np.isfinite(upper) & (omega > 0) & np.isfinite(omega))
        if np.any(bad):
            raise FixedDepthConvergenceError(f"parallel momentum bounds unresolved at {_first_loss(grid, bad)}")
        offsets = omega if epsilon_momentum == "k_total" else np.zeros_like(omega)
        # For nonnegative dispersions static F is largest at the smallest q.
        i = int(np.hypot(lower, offsets).argmin())
        try:
            _checked_static_response(material, np.hypot(lower[i], offsets[i]))
        except NumericalConvergenceError as exc:
            raise FixedDepthConvergenceError(f"fixed-depth at {_loss_point(grid, (i,))}: {exc}") from exc
        terms = np.empty((len(TERM_NAMES), grid.size))
        error = np.empty(grid.size)
        batch = None
        for i, t in enumerate(grid):
            try:
                with np.errstate(over="raise", divide="raise", invalid="raise"):
                    if i % LOSS_BATCH == 0 and grid.size > 1:
                        block = slice(i, i + LOSS_BATCH)
                        try:
                            batch = _initial_panels(material, grid[block], energy, depth, omega[block],
                                                    offsets[block], lower[block], upper[block], limit)
                        except (ArithmeticError, NumericalConvergenceError, _PanelLimitError):
                            # Replay individually to report the original loss and grid index.
                            batch = None
                    terms[:, i], error[i] = _integrate_loss(
                        material, t, energy, depth, omega[i], offsets[i], lower[i], upper[i],
                        charge, epsabs, epsrel, limit, batch[i % LOSS_BATCH] if batch is not None else None)
            except _PanelLimitError as exc:
                raise ValueError(f"fixed-depth at {_loss_point(grid, (i,))}: {exc}") from exc
            except (ArithmeticError, NumericalConvergenceError) as exc:
                raise FixedDepthConvergenceError(f"fixed-depth oscillator sum at {_loss_point(grid, (i,))}: {exc}") from exc
        no_hole, hole, total = _component_sums(terms.T, charge).T
        bad = ~(np.isfinite(total) & np.isfinite(hole) & np.isfinite(no_hole))
        if np.any(bad):
            raise FixedDepthConvergenceError(f"fixed-depth oscillator sum outside numerical range at {_first_loss(grid, bad)}")
    return FixedDepthResult(energy, depth, depth, theta, charge, epsilon_momentum,
                            grid, *terms, no_hole, hole, total, error)
