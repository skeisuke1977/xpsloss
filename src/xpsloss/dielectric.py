"""Causal reconstruction of the oscillator ELF, with the exp(-i omega t) convention."""

import cmath
from dataclasses import dataclass
import math

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

from .model import Material, _loss_coordinates, _oscillator_energy, finite_number, loss_function


DIELECTRIC_MODEL_ID = "causal_oscillator_elf_v1"
FLOAT_EPSILON = np.finfo(float).eps


class DielectricConvergenceError(RuntimeError):
    """The numerical dielectric response could not be resolved."""


class DielectricGapEdgeError(ValueError):
    """One or more evaluation energies are within the hard-gap exclusion."""

    remedy = "Shift the loss grid so all points are farther from gap_eV than gap_edge_exclusion_eV."

    def __init__(self, gap_eV, exclusion_eV, energies):
        super().__init__("dielectric response is singular at the hard gap edge; nearby points are excluded")
        self.gap_eV = gap_eV
        self.exclusion_eV = exclusion_eV
        self.offending_loss_eV = energies.tolist()


@dataclass(frozen=True)
class DielectricResponse:
    """Dimensionless responses and the quadrature error estimate for Re[-1/epsilon]."""

    negative_inverse_epsilon: complex | np.ndarray
    epsilon: complex | np.ndarray
    real_quadrature_error: float | np.ndarray


@dataclass(frozen=True)
class DielectricDiagnostics:
    """Static classification and the real-axis pole diagnostic at one q."""

    static_negative_inverse_epsilon: float
    static_sign_tolerance: float
    static_status: str
    gap_edge_exclusion_eV: float
    sub_gap_pole_eV: float | None
    sub_gap_pole_status: str


@dataclass(frozen=True)
class TabulatedKK:
    """Real part of F and the interpolation/extrapolation diagnostics.

    real_quadrature_error is zero because the transform is analytic. It
    excludes rounding, sampling, measurement and extrapolation errors.
    """

    real_negative_inverse_epsilon: float | np.ndarray
    real_quadrature_error: float | np.ndarray
    sample_bounds_eV: tuple[float, float]
    low_energy: str
    high_energy: str
    tail_coefficient_eV3: float | None
    terminal_log_slope: float | None
    warnings: tuple[str, ...]
    interpolation: str = "piecewise_linear"
    method: str = "analytic_piecewise_linear_kk"


def _quadrature(integrand, upper, points, settings):
    epsabs, epsrel, limit = settings
    if len(points) >= limit:
        raise ValueError("dielectric integration limit is too small for oscillator breakpoints")
    result = quad(integrand, 0, upper, points=points, epsabs=epsabs,
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


def _gap_edge_exclusion_eV(material):
    gap = material.gap_eV
    if gap == 0 or all(o.strength_eV2 == 0 for o in material.oscillators):
        return 0.0
    return min(gap / 2, max(32 * math.ulp(gap), 1e-12 * gap))


def _check_gap_edge(material, energies):
    gap = material.gap_eV
    if gap == 0 or all(o.strength_eV2 == 0 for o in material.oscillators):
        return
    exclusion = _gap_edge_exclusion_eV(material)
    offending = energies[abs(energies - gap) <= exclusion]
    if offending.size:
        raise DielectricGapEdgeError(gap, exclusion, np.unique(offending))


def _static_response(material, q, epsabs, epsrel, limit):
    response, error = _inverse_response(material, np.zeros_like(q), q, epsabs, epsrel, limit)
    real = response.real
    count = sum(o.strength_eV2 > 0 for o in material.oscillators)
    # Scale from 1 and the positive static KK contribution F(0)+1.
    tolerance = error + (8 * (count + 1) * FLOAT_EPSILON) * (1 + abs(real + 1))
    return real, tolerance


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
    gap edge, F diverges for any material with nonzero oscillator strength;
    points within the documented exclusion interval raise a gap edge error.
    Quadrature settings apply for a positive gap.
    """
    t, q, epsabs, epsrel = _response_inputs(loss_eV, q_inv_A, epsabs, epsrel, limit)
    _check_gap_edge(material, t)
    response, _ = _inverse_response(material, t, q, epsabs, epsrel, limit)
    return complex(response) if response.ndim == 0 else response


def dielectric_response(
    material: Material, loss_eV, q_inv_A=0.0, *, epsabs=1e-9, epsrel=1e-7, limit=200,
):
    """Return F, epsilon=-1/F, and the Re(F) error estimate.

    Reject a clearly positive static F at any sampled q. A static F within
    the roundoff and quadrature tolerance is marginal: allow T > 0 and
    report an epsilon pole at T = 0.
    """
    t, q, epsabs, epsrel = _response_inputs(loss_eV, q_inv_A, epsabs, epsrel, limit)
    _check_gap_edge(material, t)
    unique_q, q_indices = np.unique(q, return_inverse=True)
    static, tolerance = _static_response(material, unique_q, epsabs, epsrel, limit)
    if np.any(static > tolerance):
        raise ValueError("positive static -1/epsilon gives epsilon a pole in the upper half-plane")
    marginal = (abs(static) <= tolerance)[q_indices].reshape(t.shape)
    if np.any((t == 0) & marginal):
        raise ValueError("epsilon has a pole at zero energy for a marginal static response")
    response, error = _inverse_response(material, t, q, epsabs, epsrel, limit)
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


def dielectric_diagnostics(
    material: Material, q_inv_A=0.0, *, epsabs=1e-9, epsrel=1e-7, limit=200,
):
    """Classify static F and locate a sub-gap epsilon pole at one wavevector."""
    _, q, epsabs, epsrel = _response_inputs(0, q_inv_A, epsabs, epsrel, limit)
    if q.ndim:
        raise ValueError("dielectric diagnostics require a single q_inv_A value")
    static, tolerance = _static_response(material, q, epsabs, epsrel, limit)
    static, tolerance = float(static), float(tolerance)
    exclusion = _gap_edge_exclusion_eV(material)
    if abs(static) <= tolerance:
        return DielectricDiagnostics(static, tolerance, "marginal", exclusion, None, "marginal_static")
    if static > tolerance:
        return DielectricDiagnostics(static, tolerance, "unstable", exclusion, None, "absent")
    if exclusion == 0:
        return DielectricDiagnostics(static, tolerance, "stable", exclusion, None, "absent")

    def sub_gap_real(energy):
        value, _ = _inverse_response(material, np.asarray(energy), q, epsabs, epsrel, limit)
        return float(value.real)

    gap = material.gap_eV
    upper = math.nextafter(gap, 0.0)
    if sub_gap_real(upper) < 0:
        # A zero is guaranteed by the logarithmic divergence, but the interval
        # between the nearest float below the gap and the gap cannot be resolved.
        return DielectricDiagnostics(static, tolerance, "stable", exclusion, None,
                                     "unresolved_at_gap_edge")
    pole = brentq(sub_gap_real, 0, upper, xtol=math.ulp(gap), rtol=4 * FLOAT_EPSILON)
    pole_status = "within_gap_edge_exclusion" if gap - pole <= exclusion else "located"
    return DielectricDiagnostics(static, tolerance, "stable", exclusion, pole, pole_status)


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
    resonance_condition = max(1.0, center**2 / math.sqrt(denominator))
    if 16 * FLOAT_EPSILON * resonance_condition * roundoff_scale > max(1e-10, 1e-9 * abs(value)):
        raise DielectricConvergenceError("closed-form KK is ill-conditioned; use numerical KK")
    return value


def negative_inverse_dielectric_closed_form(material: Material, loss_eV, q_inv_A=0.0):
    """Return F=-1/epsilon using an independent closed form for hard-gap KK.

    The retained-band rational integral uses Eqs. D2-D3 in DIELECTRIC_DESIGN.md.
    A zero gap uses the same analytic E4b as the standard API. Both paths
    take Im(F) from the same oscillator ELF.
    The numerical API is preferred where cancellation limits the closed form.
    """
    if material.gap_eV == 0 or all(o.strength_eV2 == 0 for o in material.oscillators):
        return negative_inverse_dielectric_function(material, loss_eV, q_inv_A)
    t, q, _, _ = _response_inputs(loss_eV, q_inv_A, 1e-9, 1e-7, 200)
    gap = material.gap_eV
    _check_gap_edge(material, t)
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


def _linear_cauchy_outside(pole, nodes, values):
    """Integrate each linear segment over 1/(u-pole), away from a cutoff."""
    ratio = np.diff(nodes) / (nodes[:-1] - pole)
    logarithm = np.log1p(ratio)
    small = abs(ratio) < 1e-3
    # 1 - log(1+r)/r loses precision when |r| is small.
    shape = np.empty_like(ratio)
    r = ratio[small]
    shape[small] = r * (0.5 + r * (-1 / 3 + r * (0.25 - r / 5)))
    shape[~small] = 1 - logarithm[~small] / ratio[~small]
    return float(values[:-1] @ logarithm + np.diff(values) @ shape)


def _linear_cauchy(pole, nodes, values, slopes, tail_level=0.0):
    """Exact segment transform with internal endpoint logarithms collected."""
    lower, upper = float(nodes[0]), float(nodes[-1])
    if pole < -upper or pole > 2 * upper or 0 <= pole < lower / 2:
        result = _linear_cauchy_outside(pole, nodes, values)
        if tail_level:
            result += _power_law_tail(pole, upper, tail_level)
        return result

    dist = nodes[1:-1] - pole
    nonzero = dist != 0
    dist = dist[nonzero]
    log_scale = math.log(upper)
    # The inner log coefficient is (slope_right - slope_left)*(u_i-pole).
    # Its product with log|u_i-pole| has a zero limit at a sample node.
    inner = np.diff(slopes)[nonzero] @ (dist * (np.log(abs(dist)) - log_scale))
    left = float(slopes[0]) * (lower - pole) - float(values[0])
    right = float(values[-1]) - float(slopes[-1]) * (upper - pole)
    tail = 0.0
    if tail_level and pole >= upper / 2:
        ratio = upper / pole
        tail_weight = tail_level * ratio**3
        right -= tail_weight
        tail = tail_weight * math.log1p(pole / upper) - 2 * tail_level * ratio**2
    elif tail_level:
        tail = _power_law_tail(pole, upper, tail_level)
    if pole == lower and left:
        raise ValueError("tabulated KK diverges at the discontinuous lower cutoff")
    if pole == upper and right:
        raise ValueError("tabulated KK diverges at the discontinuous upper cutoff")
    left_log = left * (math.log(abs(lower - pole)) - log_scale) if left else 0.0
    right_log = right * (math.log(abs(upper - pole)) - log_scale) if right else 0.0
    return math.fsum((float(values[-1]) - float(values[0]), inner, left_log, right_log, tail))


def tabulated_elf_kk(sample_eV, elf, loss_eV, *, low_energy, high_energy):
    """Analytic KK for the piecewise-linear interpolant of one sampled ELF.

    Between samples use linear interpolation. Below the first sample choose
    "zero" or "linear_to_origin". Above the last sample choose "hard_cut"
    or "power_law_3" (continuous u**-3 tail). Losses can be scalars or arrays.
    """
    if np.iscomplexobj(loss_eV):
        raise ValueError("loss must be real")
    t, _ = _loss_coordinates(loss_eV, 0)
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
    warnings = ()
    if high_energy == "power_law_3":
        if terminal_slope is None:
            warnings = ("The terminal log-slope is undefined; the u^-3 tail cannot be checked.",)
        elif abs(terminal_slope + 3) > 0.5:
            warnings = ("The terminal log-slope differs from -3 by more than 0.5; check the u^-3 extrapolation.",)
    if low_energy == "linear_to_origin" and lower > 0:
        nodes, values = np.insert(nodes, 0, 0.0), np.insert(values, 0, 0.0)
    real = np.empty(t.shape)
    try:
        with np.errstate(over="raise", invalid="raise", divide="raise"):
            tail_coefficient = tail_level * upper**3 if tail_level else 0.0
            slopes = np.diff(values) / np.diff(nodes)
            for index in np.ndindex(t.shape):
                energy = float(t[index])
                positive = _linear_cauchy(energy, nodes, values, slopes, tail_level)
                negative = _linear_cauchy(-energy, nodes, values, slopes)
                real[index] = -1 + (positive + negative) / math.pi
    except (OverflowError, ZeroDivisionError, FloatingPointError) as exc:
        raise ValueError("tabulated KK is outside the numerical range") from exc
    if not math.isfinite(tail_coefficient) or not np.all(np.isfinite(real)):
        raise ValueError("nonfinite tabulated KK result")
    if t.ndim == 0:
        real, error = float(real), 0.0
    else:
        error = np.zeros(t.shape)
    return TabulatedKK(
        real, error, (lower, upper), low_energy, high_energy,
        tail_coefficient if high_energy == "power_law_3" else None,
        terminal_slope, warnings,
    )
