"""Simonsen, Yubero & Tougaard (1997), Eqs. (13)–(17), at normal emission."""

from dataclasses import dataclass
import math

import numpy as np
from scipy.optimize import brentq
from scipy.special import roots_legendre

from .dielectric import dielectric_response
from .model import (
    BOHR_RADIUS_A, HBAR2_OVER_2ME_EV_A2, Material, _oscillator_energy, finite_number,
)


FIXED_DEPTH_MODEL_ID = "simonsen1997_fixed_depth_v0"
EPSILON_MOMENTA = ("k_parallel", "k_total")
TERM_NAMES = ("term_13a", "term_13b", "term_14", "term_15_T1", "term_15_T2",
              "term_15_T3", "term_16_e", "term_16_h")
_LOW_RULE = roots_legendre(16)
_HIGH_RULE = roots_legendre(32)


class FixedDepthConvergenceError(RuntimeError):
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
    for osc in material.oscillators:
        if osc.strength_eV2:
            center = _oscillator_energy(osc, q)
            detuning = (center - loss) * (center + loss)
            response += osc.strength_eV2 / (detuning - 1j * osc.width_eV * loss)
    return response


def _term_integrands(k, material, loss, energy, depth, omega, q_offset):
    """Eight real densities in d ln(k_parallel), in the order of M1 §4.2."""
    response = _oscillator_sum(material, loss, np.hypot(k, q_offset))
    bulk = -response                            # 1/ε − 1
    reflection = response / (2 - response)      # (ε − 1)/(ε + 1)
    surface = 2 * response - reflection         # (ε − 1)(ε + 2)/[ε(ε + 1)]

    denom = omega**2 + k**2
    plus = (omega + 1j*k) / denom
    minus = plus.conj()
    # Contract conjugate phases first. Every remaining phase carries exp(-ka).
    exponent = -k*depth - 1j*(omega*depth)
    u = np.exp(exponent)
    diff = -np.expm1(exponent)                  # 1 − exp(-ka − iΩa)
    norm = abs(diff)**2 / denom
    mixed = u * diff.conj()                     # u − exp(-2ka)

    measure = k**2 / (math.pi * BOHR_RADIUS_A * energy)
    pref = measure / (2 * depth * omega)
    term_13a = measure * response.imag / denom
    # The A and B terms of 13b form a conjugate pair. C supplies |1-u|².
    term_13b = pref * (2 * bulk.imag * (plus**2 * diff).imag + bulk.real * norm)
    term_14 = pref * (bulk * (plus*diff + minus*mixed)).real / omega
    term_15_T1 = pref * surface.real * norm
    term_15_T2 = pref * (surface * minus * mixed).real / omega
    term_15_T3 = pref * (reflection * minus**2 * diff.conj()).real
    term_16_e = -pref * (reflection * plus * (2j*k/denom - plus*u)).real
    term_16_h = -pref * (reflection * plus * u).real / omega
    return np.stack((term_13a, term_13b, term_14, term_15_T1, term_15_T2,
                     term_15_T3, term_16_e, term_16_h), axis=-1)


def _line_positions(material, loss, q_offset, upper_squared):
    """Centers and half-width scales in k_parallel², used only for panel layout."""
    dispersive = []
    fixed = 0j
    for osc in material.oscillators:
        if not osc.strength_eV2:
            continue
        slope = osc.dispersion * HBAR2_OVER_2ME_EV_A2
        if slope:
            pole = (loss - osc.energy_eV) / slope - q_offset**2
            dispersive.append((osc, slope, pole))
            yield pole, osc.width_eV / (2 * slope)
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
               osc.width_eV * loss / (2 * center * slope))
        return

    # Surface lines belong to the oscillator *sum*. In the narrow-line limit
    # Re(S) decreases between successive bulk poles in z=k_parallel². Locate
    # its crossings of 2 and linearize the damping there to seed both flanks.
    poles = {pole for _, _, pole in dispersive}
    stop = 2 * upper_squared
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
            raise FixedDepthConvergenceError(f"surface line unresolved at loss_eV={loss:g}") from exc
        if z in poles:
            raise FixedDepthConvergenceError(f"unresolved surface and bulk lines at loss_eV={loss:g}")
        damping, gradient = fixed.imag, 0.0
        for osc, slope, pole in dispersive:
            center = osc.energy_eV + slope * (z + q_offset**2)
            detuning = slope * (z - pole) * (center + loss)
            weight = osc.strength_eV2 / detuning**2
            damping += weight * osc.width_eV * loss
            gradient += weight * 2 * center * slope
        yield z, damping / gradient


def _log_breakpoints(material, loss, lower, upper, q_offset):
    """Log panels, split at the dispersive oscillator and surface lines."""
    lo, hi = math.log(lower), math.log(upper)
    points = list(np.linspace(lo, hi, math.ceil(hi - lo) + 1))
    lower_squared, upper_squared = lower**2, upper**2
    for center, width in _line_positions(material, loss, q_offset, upper_squared):
        if lower_squared < center < upper_squared:
            points.append(math.log(center) / 2)
        # Geometric flanks resolve the tails even on a much wider log panel.
        offset = width
        while 0 < offset < 8 * max(abs(center - lower_squared), abs(upper_squared - center)):
            for flank in (center - offset, center + offset):
                if lower_squared < flank < upper_squared:
                    points.append(math.log(flank) / 2)
            offset *= 8
    return np.unique(points)


def _panel_integrals(edges_left, edges_right, integrand, error_weights):
    mid = (edges_left + edges_right) / 2
    half = (edges_right - edges_left) / 2

    def integrate(rule):
        nodes, weights = rule
        values = integrand(np.exp(mid[:, None] + half[:, None] * nodes))
        weighted = (half[:, None] * weights)[..., None] * values
        return weighted.sum(axis=1), abs(weighted).sum(axis=1)

    low, _ = integrate(_LOW_RULE)
    high, magnitude = integrate(_HIGH_RULE)
    uncertainty = np.maximum(abs(high - low), 32 * np.finfo(float).eps * magnitude)
    return high, uncertainty @ error_weights


def _integrate_loss(material, loss, energy, depth, omega, q_offset, lower, upper,
                    charge, epsabs, epsrel, limit):
    edges = _log_breakpoints(material, loss, lower, upper, q_offset)
    # Refining one panel must not move either endpoint of its neighbor.
    left, right = edges[:-1].copy(), edges[1:].copy()
    if len(left) > limit:
        raise ValueError(f"integration limit too small for initial panels at loss_eV={loss:g}")

    def integrand(k):
        return _term_integrands(k, material, loss, energy, depth, omega, q_offset)

    # The same estimate covers individual terms, both unit-charge components
    # and their user-weighted sum, even when signed contributions cancel.
    error_weights = np.ones(len(TERM_NAMES))
    hole_terms = [2, 4, 7]
    error_weights[hole_terms] = max(1, abs(charge))
    values, errors = _panel_integrals(left, right, integrand, error_weights)
    while True:
        terms = values.sum(axis=0)
        error = errors.sum()
        no_hole = terms[[0, 1, 3, 5, 6]].sum()
        total = no_hole + charge * terms[hole_terms].sum()
        if error <= max(epsabs, epsrel * abs(total)):
            return terms, error
        if len(left) >= limit:
            raise FixedDepthConvergenceError(
                f"fixed-depth quadrature failed to converge at loss_eV={loss:g} "
                f"with {limit} panels (estimated error {error:.3g})")
        i = errors.argmax()
        midpoint = (left[i] + right[i]) / 2
        children_left = np.array([left[i], midpoint])
        children_right = np.array([midpoint, right[i]])
        child_values, child_errors = _panel_integrals(children_left, children_right, integrand, error_weights)
        right[i], values[i], errors[i] = midpoint, child_values[0], child_errors[0]
        left = np.append(left, midpoint)
        right = np.append(right, children_right[1])
        values = np.vstack((values, child_values[1]))
        errors = np.append(errors, child_errors[1])


def _finite_real(value, label):
    if np.iscomplexobj(value):
        raise ValueError(f"{label} must be real")
    return finite_number(value, label)


def effective_cross_section_fixed_depth(
    material: Material, electron_energy_eV, loss_eV, depth_A, *, theta_deg=0.0,
    hole_charge=1.0, epsilon_momentum, epsabs=1e-12, epsrel=1e-9, limit=200,
) -> FixedDepthResult:
    """Fixed-depth K for a scalar or 1-D loss grid, returned as 1-D arrays.

    θ is measured from the surface normal. All hole terms refer to unit charge.
    epsabs is in eV^-1 Å^-1, epsrel is relative to total, and limit caps panels
    per loss. The quadrature estimate excludes material and model uncertainties.
    """
    if not isinstance(material, Material):
        raise ValueError("material must be a Material")
    if material.gap_eV > 0:
        raise ValueError("fixed-depth gapped materials are deferred to M2-2")
    if not isinstance(epsilon_momentum, str) or epsilon_momentum not in EPSILON_MOMENTA:
        raise ValueError("epsilon_momentum must be 'k_parallel' or 'k_total'")
    theta = _finite_real(theta_deg, "theta_deg")
    if theta != 0:
        raise ValueError("only normal emission (theta_deg=0) is supported; see Q3/Q15")
    energy = _finite_real(electron_energy_eV, "electron_energy_eV")
    depth = _finite_real(depth_A, "depth_A")
    charge = _finite_real(hole_charge, "hole_charge")
    if depth <= 0 or energy <= 0:
        raise ValueError("depth_A and electron_energy_eV must be positive")
    grid = np.atleast_1d(loss_eV)
    if grid.dtype.kind not in "fiu" or grid.ndim != 1 or not grid.size:
        raise ValueError("loss_eV must be a real scalar or a nonempty one-dimensional grid")
    grid = np.array(grid, dtype=float)
    if not np.all((grid > 0) & (grid < energy)):
        raise ValueError("require 0 < loss_eV < electron_energy_eV; the zero-loss delta term is excluded")
    epsabs, epsrel = _finite_real(epsabs, "epsabs"), _finite_real(epsrel, "epsrel")
    if epsabs <= 0 or epsrel <= 0:
        raise ValueError("epsabs and epsrel must be positive")
    if isinstance(limit, bool) or not isinstance(limit, int) or limit < 1:
        raise ValueError("limit must be a positive integer")

    with np.errstate(over="raise", divide="raise", invalid="raise", under="ignore"):
        try:
            lower, upper = _parallel_limits(energy, grid)
            if not np.all((lower > 0) & (lower < upper)):
                raise ValueError("parallel momentum bounds outside numerical range")
            omega = grid / (2 * math.sqrt(HBAR2_OVER_2ME_EV_A2 * energy))
            offsets = omega if epsilon_momentum == "k_total" else np.zeros_like(omega)
            # For nonnegative dispersions, static F is largest at the smallest
            # sampled q. Apply the same upper-half-plane pole check as M1 there.
            i = np.hypot(lower, offsets).argmin()
            dielectric_response(material, grid[i], np.hypot(lower[i], offsets[i]))
            terms = np.empty((len(TERM_NAMES), grid.size))
            error = np.empty(grid.size)
            for i, t in enumerate(grid):
                terms[:, i], error[i] = _integrate_loss(
                    material, t, energy, depth, omega[i], offsets[i], lower[i], upper[i],
                    charge, epsabs, epsrel, limit)
            no_hole = terms[0] + terms[1] + terms[3] + terms[5] + terms[6]
            hole = terms[2] + terms[4] + terms[7]
            total = no_hole + charge * hole
        except (FloatingPointError, OverflowError, ZeroDivisionError) as exc:
            raise ValueError("fixed-depth calculation outside numerical range") from exc
    return FixedDepthResult(energy, depth, depth, theta, charge, epsilon_momentum,
                            grid, *terms, no_hole, hole, total, error)
