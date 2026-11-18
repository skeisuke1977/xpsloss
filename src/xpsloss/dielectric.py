"""Causal reconstruction of the oscillator ELF, with the exp(-i omega t) convention."""

import cmath
from dataclasses import dataclass
import math

import numpy as np
from scipy.integrate import quad

from .model import Material, _loss_coordinates, _oscillator_energy, finite_number, loss_function


DIELECTRIC_MODEL_ID = "causal_oscillator_elf_v1"


class DielectricConvergenceError(RuntimeError):
    """The numerical dielectric response could not be resolved."""


@dataclass(frozen=True)
class DielectricResponse:
    """Dimensionless responses and the quadrature error estimate for Re[-1/epsilon]."""

    negative_inverse_epsilon: complex | np.ndarray
    epsilon: complex | np.ndarray
    real_quadrature_error: float | np.ndarray


@dataclass(frozen=True)
class TabulatedKK:
    """Real part of F and the settings for reconstructing it from one ELF table."""

    real_negative_inverse_epsilon: float | np.ndarray
    real_quadrature_error: float | np.ndarray
    sample_bounds_eV: tuple[float, float]
    low_energy: str
    high_energy: str
    tail_coefficient_eV3: float | None
    terminal_log_slope: float | None
    epsabs: float
    epsrel: float
    limit: int
    interpolation: str = "piecewise_linear"


def _quadrature(integrand, upper, points, settings, *, lower=0):
    epsabs, epsrel, limit = settings
    if len(points) >= limit:
        raise DielectricConvergenceError("dielectric integration limit is too small for oscillator breakpoints")
    result = quad(integrand, lower, upper, points=points, epsabs=epsabs,
                  epsrel=epsrel, limit=limit, full_output=True)
    if len(result) > 3:
        raise DielectricConvergenceError(f"dielectric integration failed to converge: {result[3]}")
    return result[0], result[1]


def _breakpoints(start, stop, scale):
    """Space quadrature breakpoints geometrically away from a narrow feature."""
    points = []
    span = abs(stop - start)
    direction = 1 if stop > start else -1
    offset = scale
    while 0 < offset < span:
        point = start + direction * offset
        if point != start and point != stop:
            points.append(point)
        offset *= 8
    return points


def _removed_kk_band(energy, gap, resonance, strength, width, settings):
    """KK contribution from [0, gap], with analytic principal-value subtraction."""
    def denominator(u):
        return ((resonance - u) * (resonance + u))**2 + (width * u)**2

    prefactor = 2 * strength * width / math.pi
    if 0 < energy < gap or 0 < energy - gap <= min(gap / 2, max(width, resonance - gap)):
        anchor = min(energy, gap)
        at_anchor = denominator(anchor)
        resonance_squared = resonance**2

        def integrand(u):
            difference = resonance * (resonance - anchor) + anchor * (resonance - u)
            divided_difference = difference * (resonance_squared + u * anchor)
            regular = (prefactor / at_anchor) * (divided_difference / denominator(u))
            if anchor == energy:
                return regular
            return regular * ((u - anchor) / (u - energy)) * ((u + anchor) / (u + energy))

        log_ratio = math.log(abs(gap - energy)) - math.log(gap + energy)
        analytic = (prefactor / at_anchor) * (anchor / 2) * (anchor / energy) * log_ratio
    elif energy == 0:
        def integrand(u):
            return prefactor / denominator(u)

        analytic = 0.0
    else:
        def integrand(u):
            return (prefactor / denominator(u)) * (u / (u - energy)) * (u / (u + energy))

        analytic = 0.0

    points = [p for p in _breakpoints(resonance, 0, max(width, resonance - gap)) if p < gap]
    overdamped_scale = resonance * (resonance / width)
    points.extend(_breakpoints(0, gap, overdamped_scale))
    integral, error = _quadrature(integrand, gap, sorted(set(points)), settings)
    return integral + analytic, error


def _retained_kk_band(energy, gap, resonance, strength, width, settings):
    """Integrate [gap, infinity) in x=gap/u when the resonance is below the gap."""
    r, w, s = resonance / gap, width / gap, energy / gap
    relative_spacing = (gap - resonance) / gap
    prefactor = (2 / math.pi) * (strength / gap / gap) * w

    def detuning(x):
        return relative_spacing + r * (1 - x)

    def denominator(x):
        return (detuning(x) * (1 + r * x))**2 + (w * x)**2

    if energy <= gap / 2:
        def integrand(x):
            return prefactor * x**2 / (denominator(x) * (1 - (s * x)**2))

        analytic = 0.0
    else:
        anchor = min(1.0, gap / energy)
        at_anchor = denominator(anchor)

        def integrand(x):
            product = (r * x) * (r * anchor)
            difference = detuning(x) + r * x * detuning(anchor)
            regular = (prefactor / at_anchor) * (difference * (1 + product) / denominator(x))
            if energy > gap:
                return -regular * anchor**2
            gap_distance = (gap - energy) / gap + s * (1 - x)
            kernel = ((x - 1) * (x + 1)) / (gap_distance * (1 + s * x))
            return regular * kernel

        log_ratio = math.log(gap + energy) - math.log(abs(gap - energy))
        analytic = (prefactor / at_anchor) * anchor**2 * (gap / energy / 2) * log_ratio

    points = _breakpoints(1.0, 0, max(relative_spacing, w))
    points.extend(_breakpoints(0, 1.0, 1 / w))
    integral, error = _quadrature(integrand, 1.0, sorted(set(points)), settings)
    return integral + analytic, error


def _inverse_response(material, energies, q, epsabs, epsrel, limit):
    real = np.full(energies.shape, -1.0)
    error = np.zeros(energies.shape)
    oscillators = tuple(o for o in material.oscillators if o.strength_eV2 > 0)
    if not oscillators or not energies.size:
        return real.astype(complex), error
    gap = material.gap_eV
    if gap > 0 and np.any(energies == gap):
        raise ValueError("dielectric response diverges at the exact hard gap edge")

    try:
        with np.errstate(over="raise", invalid="raise", divide="raise"):
            for osc in oscillators:
                resonance = _oscillator_energy(osc, q)
                if gap == 0:
                    delta = (resonance - energies) * (resonance + energies)
                    real += (osc.strength_eV2 / (delta - 1j * osc.width_eV * energies)).real
                    continue
                settings = (epsabs / len(oscillators), epsrel, limit)
                for index in np.ndindex(energies.shape):
                    energy, center = float(energies[index]), float(resonance[index])
                    parameters = (energy, gap, center, osc.strength_eV2, osc.width_eV, settings)
                    if center < gap:
                        contribution, estimate = _retained_kk_band(*parameters)
                    else:
                        correction, estimate = _removed_kk_band(*parameters)
                        delta = (center - energy) * (center + energy)
                        bare = osc.strength_eV2 / complex(delta, -osc.width_eV * energy)
                        contribution = bare.real - correction
                    real[index] += contribution
                    error[index] += estimate
            result = real + 1j * loss_function(material, energies, q)
    except (FloatingPointError, OverflowError, ZeroDivisionError) as exc:
        raise ValueError("dielectric response is outside the numerical range") from exc
    if not np.all(np.isfinite(result)) or not np.all(np.isfinite(error)):
        raise ValueError("nonfinite dielectric response")
    return result, error


def _response_inputs(loss_eV, q_inv_A, epsabs, epsrel, limit):
    if np.iscomplexobj(loss_eV) or np.iscomplexobj(q_inv_A):
        raise ValueError("loss and q must be real")
    t, q = _loss_coordinates(loss_eV, q_inv_A)
    epsabs = finite_number(epsabs, "epsabs")
    epsrel = finite_number(epsrel, "epsrel")
    if epsabs <= 0 or epsrel <= 0:
        raise ValueError("epsabs and epsrel must be positive")
    if isinstance(limit, bool) or not isinstance(limit, int) or limit < 1:
        raise ValueError("limit must be a positive integer")
    return t, q, epsabs, epsrel


def negative_inverse_dielectric_function(
    material: Material, loss_eV, q_inv_A=0.0, *, epsabs=1e-9, epsrel=1e-7, limit=200,
):
    """Return F=-1/epsilon, with the existing ELF as Im(F) and F(infinity)=-1.

    Loss and q are nonnegative with NumPy broadcasting. At a positive hard
    gap edge, F diverges for any material with nonzero oscillator strength.
    Quadrature settings apply when gap_eV > 0.
    """
    t, q, epsabs, epsrel = _response_inputs(loss_eV, q_inv_A, epsabs, epsrel, limit)
    response, _ = _inverse_response(material, t, q, epsabs, epsrel, limit)
    return complex(response) if response.ndim == 0 else response


def dielectric_response(
    material: Material, loss_eV, q_inv_A=0.0, *, epsabs=1e-9, epsrel=1e-7, limit=200,
):
    """Return F, its causal reciprocal epsilon=-1/F, and the Re(F) error estimate.

    Reject a positive static F at any sampled q (epsilon would have an
    upper-half-plane pole) or a requested point where F vanishes. A sharp gap
    can produce a real-axis pole in epsilon below the absorption onset.
    """
    t, q, epsabs, epsrel = _response_inputs(loss_eV, q_inv_A, epsabs, epsrel, limit)
    response, error = _inverse_response(material, t, q, epsabs, epsrel, limit)
    unique_q = np.unique(q)
    static, static_error = _inverse_response(
        material, np.zeros_like(unique_q), unique_q, epsabs, epsrel, limit,
    )
    if np.any(static.real > static_error):
        raise ValueError("positive static -1/epsilon gives epsilon a pole in the upper half-plane")
    if np.any((static_error > 0) & (abs(static.real) <= static_error)):
        raise DielectricConvergenceError("the sign of the static -1/epsilon is unresolved")
    if np.any(response == 0):
        raise ValueError("epsilon has a pole where -1/epsilon vanishes")
    if np.any(abs(response) <= error):
        raise DielectricConvergenceError("epsilon is unresolved near a zero of -1/epsilon")
    try:
        with np.errstate(over="raise", invalid="raise", divide="raise"):
            epsilon = -1 / response
    except FloatingPointError as exc:
        raise ValueError("epsilon is outside the numerical range") from exc
    if response.ndim == 0:
        return DielectricResponse(complex(response), complex(epsilon), float(error))
    return DielectricResponse(response, epsilon, error)


def _edge_log(energy, edge):
    """Return log((edge + energy)/abs(edge - energy)) without subtracting logs."""
    return math.log1p(2 * (min(energy, edge) / abs(edge - energy)))


def _closed_kk_band(t, gap, center, strength, width):
    """Integrate the retained rational KK kernel by partial fractions (Eqs. D2-D3)."""
    half_width = width / 2
    splitting = cmath.sqrt(half_width - center) * cmath.sqrt(half_width + center)
    p = half_width + splitting
    b = center * (center / p)  # p*b = center**2, also in the overdamped limit
    j_b = cmath.atan(b / gap) / b
    difference = 2 * splitting
    if abs(difference) < width / 4:
        h = gap / (gap**2 + center**2)
        atan_ratio = cmath.atan(h * difference) / difference if difference else h
        j_divided = (j_b - atan_ratio) / (p * width)
    else:
        j_p = cmath.atan(p / gap) / p
        j_divided = (j_b - j_p) / (width * difference)

    denominator = ((center - t) * (center + t))**2 + (width * t)**2
    scale = (2 / math.pi) * strength * width / denominator
    regular = (center**4 + (t * p)**2) * j_divided
    log_term = (t / 2) * _edge_log(t, gap)
    pole = log_term - t**2 * j_b
    value = scale * (regular + pole).real
    roundoff_scale = scale * (abs(regular) + log_term + t**2 * abs(j_b))
    if np.finfo(float).eps * roundoff_scale > max(1e-10, 1e-9 * abs(value)):
        raise DielectricConvergenceError("closed-form KK is ill-conditioned; use numerical KK")
    return value


def negative_inverse_dielectric_closed_form(material: Material, loss_eV, q_inv_A=0.0):
    """Return F=-1/epsilon using an independent closed form for hard-gap KK.

    The retained-band rational integral uses Eqs. D2-D3 in DIELECTRIC_DESIGN.md.
    The numerical API is preferred where cancellation limits the closed form.
    """
    if material.gap_eV == 0 or all(o.strength_eV2 == 0 for o in material.oscillators):
        return negative_inverse_dielectric_function(material, loss_eV, q_inv_A)
    t, q, _, _ = _response_inputs(loss_eV, q_inv_A, 1e-9, 1e-7, 200)
    gap = material.gap_eV
    if np.any(t == gap):
        raise ValueError("dielectric response diverges at the exact hard gap edge")
    real = np.full(t.shape, -1.0)
    try:
        with np.errstate(over="raise", invalid="raise", divide="raise"):
            for osc in material.oscillators:
                if osc.strength_eV2 == 0:
                    continue
                centers = _oscillator_energy(osc, q)
                for index in np.ndindex(t.shape):
                    real[index] += _closed_kk_band(
                        float(t[index]), gap, float(centers[index]), osc.strength_eV2, osc.width_eV,
                    )
            result = real + 1j * loss_function(material, t, q)
    except (FloatingPointError, OverflowError, ZeroDivisionError) as exc:
        raise ValueError("closed-form dielectric response is outside the numerical range") from exc
    if not np.all(np.isfinite(result)):
        raise ValueError("nonfinite closed-form dielectric response")
    return complex(result) if result.ndim == 0 else result


def _optical_moment(material, first_moment, settings):
    gap = material.gap_eV
    if gap == 0:
        if first_moment:
            return (math.pi / 2) * math.fsum(osc.strength_eV2 for osc in material.oscillators)
        return math.fsum((osc.strength_eV2 / osc.energy_eV) / osc.energy_eV
                         for osc in material.oscillators)

    # Transform the positive retained integral on [gap, infinity) with x=gap/u.
    oscillators = tuple(osc for osc in material.oscillators if osc.strength_eV2 > 0)
    if not oscillators:
        return 0.0
    epsabs, epsrel, limit = settings
    settings = (epsabs / len(oscillators), epsrel, limit)
    total = 0.0
    for osc in oscillators:
        center, width = osc.energy_eV, osc.width_eV
        peak_x = gap / center
        damping = width / gap
        prefactor = osc.strength_eV2 * damping
        if not first_moment:
            prefactor = (2 / math.pi) * (prefactor / gap / gap)

        def integrand(x):
            detuning = (peak_x - x) / peak_x
            denominator = (detuning * (1 + x / peak_x))**2 + (damping * x)**2
            return prefactor * (1 if first_moment else x**2) / denominator

        span = (width / center) * peak_x
        points = _breakpoints(0, 1, gap / width)
        points.extend(_breakpoints(0, 1, span))
        if peak_x < 1:
            points.append(peak_x)
            points.extend(_breakpoints(peak_x, 0, span))
            points.extend(_breakpoints(peak_x, 1, span))
        else:
            points.extend(_breakpoints(1, 0, max(peak_x - 1, span)))
        value, _ = _quadrature(integrand, 1, sorted(set(points)), settings)
        total += value
    return total


def optical_sum_rule_lhs(material: Material, *, epsabs=1e-9, epsrel=1e-7, limit=200):
    """Return the dimensionless LHS of Eq. D4, using quadrature only for a positive gap."""
    _, _, epsabs, epsrel = _response_inputs(0, 0, epsabs, epsrel, limit)
    return _checked_optical_moment(material, False, (epsabs, epsrel, limit))


def optical_sum_rule_residual(
    material: Material, reference_inverse_epsilon, *, epsabs=1e-9, epsrel=1e-7, limit=200,
):
    """Return the sum-rule LHS minus (1 - reference Re[1/epsilon(0)])."""
    reference = finite_number(reference_inverse_epsilon, "reference_inverse_epsilon")
    return optical_sum_rule_lhs(material, epsabs=epsabs, epsrel=epsrel, limit=limit) - (1 - reference)


def optical_elf_first_moment(material: Material, *, epsabs=1e-9, epsrel=1e-7, limit=200):
    """Return the first moment in eV**2 (Eq. D4), using quadrature for a positive gap."""
    _, _, epsabs, epsrel = _response_inputs(0, 0, epsabs, epsrel, limit)
    return _checked_optical_moment(material, True, (epsabs, epsrel, limit))


def _checked_optical_moment(material, first_moment, settings):
    try:
        value = _optical_moment(material, first_moment, settings)
    except (OverflowError, ZeroDivisionError) as exc:
        raise ValueError("optical integral is outside the numerical range") from exc
    if not math.isfinite(value):
        raise ValueError("nonfinite optical integral")
    return value


def _power_law_tail(t, edge, level):
    """Return pi times the KK contribution of level*(edge/u)**3 above edge."""
    ratio = t / edge
    if ratio < 0.03:
        return 2 * level * sum(ratio**(2 * n) / (2 * n + 3) for n in range(5))
    inverse_ratio = edge / t
    return level * inverse_ratio**2 * (inverse_ratio * _edge_log(t, edge) - 2)


def _tabulated_kk_point(t, nodes, values, tail_level, settings):
    """Integrate the interpolant, combining the two endpoint logs at a matched tail."""
    lower, upper = float(nodes[0]), float(nodes[-1])
    points = nodes[1:-1].tolist()
    tail = 0.0
    if t == 0 or t < lower / 2 or t > 2 * upper:
        def integrand(u):
            return 2 * (u / (u + t)) * (float(np.interp(u, nodes, values)) / (u - t))

        if tail_level:
            tail = _power_law_tail(t, upper, tail_level)
    else:
        anchor = min(max(t, lower), upper)
        level = float(np.interp(anchor, nodes, values))
        segment = min(np.searchsorted(nodes, anchor, side="right") - 1, len(nodes) - 2)
        slope = float((values[segment + 1] - values[segment]) / (nodes[segment + 1] - nodes[segment]))
        if lower < t < upper and t not in points:
            points.append(t)

        def integrand(u):
            sample = float(np.interp(u, nodes, values))
            if nodes[segment] <= u <= nodes[segment + 1] or u == t:
                regular = slope if u == t else slope * ((u - anchor) / (u - t))
            else:
                regular = (sample - level) / (u - t)
            return regular + sample / (u + t)

        if t == lower and level:
            raise ValueError("tabulated KK diverges at the discontinuous lower cutoff")
        if t != lower and level:
            tail -= level * math.log(abs(lower - t))
        upper_weight = level
        if tail_level and t >= upper / 2:
            ratio = upper / t
            tail_weight = tail_level * ratio**3
            upper_weight -= tail_weight
            tail += tail_weight * math.log(upper + t) - 2 * tail_level * ratio**2
        elif tail_level:
            tail += _power_law_tail(t, upper, tail_level)
        if t == upper and upper_weight:
            raise ValueError("tabulated KK diverges at the discontinuous upper cutoff")
        if t != upper and upper_weight:
            tail += upper_weight * math.log(abs(upper - t))

    value, error = _quadrature(integrand, upper, points, settings, lower=lower)
    return -1 + (value + tail) / math.pi, error / math.pi


def tabulated_elf_kk(
    sample_eV, elf, loss_eV, *, low_energy, high_energy, epsabs=1e-9, epsrel=1e-7,
):
    """Reference numerical KK for one sampled ELF, with recorded extrapolations.

    Between samples use linear interpolation. Below the first sample choose
    "zero" or "linear_to_origin". Above the last sample choose "hard_cut"
    or "power_law_3" (continuous u**-3 tail). Losses can be scalars or arrays.
    """
    t, _, epsabs, epsrel = _response_inputs(loss_eV, 0, epsabs, epsrel, 200)
    if low_energy not in ("zero", "linear_to_origin"):
        raise ValueError("low_energy must be 'zero' or 'linear_to_origin'")
    if high_energy not in ("hard_cut", "power_law_3"):
        raise ValueError("high_energy must be 'hard_cut' or 'power_law_3'")
    if np.iscomplexobj(sample_eV) or np.iscomplexobj(elf):
        raise ValueError("sample energies and ELF must be real")
    nodes, values = np.asarray(sample_eV, dtype=float), np.asarray(elf, dtype=float)
    if nodes.ndim != 1 or nodes.shape != values.shape or nodes.size < 2:
        raise ValueError("sample energies and ELF must be matching 1-D arrays with at least two points")
    if not np.all(np.isfinite(nodes)) or not np.all(np.isfinite(values)):
        raise ValueError("sample energies and ELF must be finite")
    if nodes[0] < 0 or np.any(nodes[1:] <= nodes[:-1]) or np.any(values < 0):
        raise ValueError("sample energies must be nonnegative and strictly increasing, and ELF nonnegative")
    lower, upper = float(nodes[0]), float(nodes[-1])
    if lower == 0 and values[0] != 0:
        raise ValueError("the ELF sample at zero energy must be zero")

    terminal_slope = None
    if nodes[-2] > 0 and values[-2] > 0 and values[-1] > 0:
        log_spacing = math.log1p((upper - float(nodes[-2])) / float(nodes[-2]))
        terminal_slope = (math.log(float(values[-1])) - math.log(float(values[-2]))) / log_spacing
    tail_level = float(values[-1]) if high_energy == "power_law_3" else 0.0
    if low_energy == "linear_to_origin" and lower > 0:
        nodes, values = np.insert(nodes, 0, 0.0), np.insert(values, 0, 0.0)
    limit = max(200, 2 * len(nodes))
    settings = (epsabs, epsrel, limit)
    real, error = np.empty(t.shape), np.empty(t.shape)
    try:
        tail_coefficient = tail_level * upper**3 if tail_level else 0.0
        for index in np.ndindex(t.shape):
            real[index], error[index] = _tabulated_kk_point(
                float(t[index]), nodes, values, tail_level, settings,
            )
    except (OverflowError, ZeroDivisionError, FloatingPointError) as exc:
        raise ValueError("tabulated KK is outside the numerical range") from exc
    if not math.isfinite(tail_coefficient) or not np.all(np.isfinite(real)) or not np.all(np.isfinite(error)):
        raise ValueError("nonfinite tabulated KK result")
    if t.ndim == 0:
        real, error = float(real), float(error)
    return TabulatedKK(
        real, error, (lower, upper), low_energy, high_energy,
        tail_coefficient if high_energy == "power_law_3" else None,
        terminal_slope, epsabs, epsrel, limit,
    )
