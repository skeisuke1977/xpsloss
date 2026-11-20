"""Causal reconstruction of the oscillator ELF, with the exp(-i omega t) convention."""

import cmath
from dataclasses import dataclass, replace
import math

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

from .model import (COORDINATE_ULPS, DEFAULT_ABS, DEFAULT_REL, FLANK_RATIO, FLOAT_EPSILON, NEARBY_FACTOR,
                    RECIPROCAL_ROUNDOFF, SUM_ROUNDOFF, Material, NumericalConvergenceError,
                    _PanelLimitError, _accuracy, _first_loss, _loss_coordinates, _loss_point, _oscillator_energy,
                    _oscillator_response, finite_number, loss_function)


DIELECTRIC_MODEL_ID = "causal_oscillator_elf_v1"
GAP_EDGE_ULPS = 32  # Leave several representable neighbors outside the singular endpoint.
GAP_EDGE_REL = 1e-12  # Relative hard-edge exclusion specified in DIELECTRIC_DESIGN.
KK_REFINEMENT = 8  # Reduce the preceding absolute error per line by three binary digits per pass.
TAB_ROUNDOFF = 32 * FLOAT_EPSILON  # Local weights, endpoint logarithms and compensated reduction.
SMALL_SEGMENT_RATIO = 0.01  # Eighth-order log1p weight has relative remainder below 1e-16 here.
TAIL_SERIES_RATIO = 0.03  # Five even powers suffice below this loss/cutoff ratio.
TAIL_SLOPE_WARNING = 0.5  # Warn when the measured terminal exponent differs from -3 by half a unit.
CLOSED_ABS = 1e-10  # Per-line cancellation budget for the auxiliary partial-fraction check.
CLOSED_REL = 1e-9  # Leave two relative guard digits before combining closed-form terms.
CLOSED_CENTER_REL = 1e-8  # Reserve one relative digit for the gap-to-pole distance at finite q.
CRITICAL_SPLITTING_FRACTION = 0.25  # Use a divided difference when the two damped poles coalesce.

GAP_EDGE_WARNING_FRACTION = 1e-3  # Diagnostic band expressed as a fraction of the gap.
BRENT_MIN_REL = 4 * FLOAT_EPSILON  # Smallest relative root tolerance supported by brentq.


class DielectricConvergenceError(NumericalConvergenceError):
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
    """Dimensionless responses and estimated absolute error of the real KK integral.

    For a hard gap the estimate includes quadrature and arithmetic used to
    reconstruct and sum the bands. Gapless analytic responses have zero
    quadrature error; their component roundoff is checked separately.
    """

    negative_inverse_epsilon: complex | np.ndarray
    epsilon: complex | np.ndarray
    real_quadrature_error: float | np.ndarray


@dataclass(frozen=True)
class GapEdgePoint:
    """A sampled loss in the warning band outside the gap exclusion."""

    loss_eV: float
    distance_eV: float
    edge_log_contribution: float


@dataclass(frozen=True)
class DielectricDiagnostics:
    """Static sign and real-axis pole diagnostics for the oscillator representation.

    At the given q, static_sign compares Re[-1/epsilon(q, 0)] to the tolerance
    τ0 (static_sign_tolerance): below -τ0 is negative, within ±τ0 is marginal,
    and above τ0 is positive. A positive sign means epsilon of this representation
    has an upper-half-plane pole, so dielectric_response rejects the response.
    This diagnostic does not test material stability.
    """

    static_negative_inverse_epsilon: float
    static_sign_tolerance: float
    static_sign: str
    gap_edge_exclusion_eV: float
    sub_gap_pole_eV: float | None
    sub_gap_pole_status: str
    gap_edge_warning_band_eV: float = 0.0
    near_gap_edge_points: tuple[GapEdgePoint, ...] = ()

    @property
    def warnings(self):
        return tuple(
            f"Near the hard gap edge, Re(F) depends logarithmically on the grid offset "
            f"(loss_eV={point.loss_eV:.17g}, distance_eV={point.distance_eV:.17g}, "
            f"edge_log_contribution={point.edge_log_contribution:.9g})."
            for point in self.near_gap_edge_points
        )


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
        raise _PanelLimitError("dielectric integration limit is too small for oscillator breakpoints")
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
        offset *= FLANK_RATIO
    return points


def _removed_kk_band(energy, gap, resonance, strength, width, settings):
    """KK contribution from [0, gap], with analytic principal-value subtraction."""
    def denominator(u):
        return ((resonance - u) * (resonance + u))**2 + (width * u)**2

    prefactor = 2 * strength * width / math.pi
    if 0 < energy < gap or 0 < energy - gap <= min(gap / NEARBY_FACTOR, max(width, resonance - gap)):
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
    # Above the gap the anchor is the endpoint. The regularized integrand
    # still turns over at gap-u ~ energy-gap, independently of the linewidth.
    if energy > gap:
        points.extend(_breakpoints(gap, 0, energy - gap))
    overdamped_scale = resonance * (resonance / width)
    points.extend(_breakpoints(0, gap, overdamped_scale))
    integral, error = _quadrature(integrand, gap, sorted(set(points)), settings)
    return integral + analytic, error + SUM_ROUNDOFF * (abs(integral) + abs(analytic))


def _retained_kk_band(energy, gap, resonance, strength, width, settings):
    """Integrate [gap, infinity) in x=gap/u when the resonance is below the gap."""
    r, w, s = resonance / gap, width / gap, energy / gap
    relative_spacing = (gap - resonance) / gap
    prefactor = (2 / math.pi) * (strength / gap / gap) * w

    def detuning(x):
        return (1 - r * x) if r > NEARBY_FACTOR else relative_spacing + r * (1 - x)

    def denominator(x):
        return (detuning(x) * (1 + r * x))**2 + (w * x)**2

    if energy <= gap / NEARBY_FACTOR:
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
    if gap / NEARBY_FACTOR < energy < gap:
        points.extend(_breakpoints(1.0, 0, (gap - energy) / gap))
    points.extend(_breakpoints(0, 1.0, 1 / w))
    integral, error = _quadrature(integrand, 1.0, sorted(set(points)), settings)
    return integral + analytic, error + SUM_ROUNDOFF * (abs(integral) + abs(analytic))


def _gapped_real(lines, energy, gap, epsabs, epsrel, limit, location):
    """Sum KK bands, refining to a shared absolute budget when terms cancel."""
    count = len(lines)
    # Use the same pilot and refinement sequence at every requested tolerance.
    # A looser tolerance accepts an earlier pass or the same pass.
    abs_goal, rel_goal = DEFAULT_ABS / count, DEFAULT_REL
    while True:
        values, estimates = [-1.0], []
        for i, center, osc, bare_real in lines:
            try:
                flank = max(abs(gap - center), osc.width_eV)
                if COORDINATE_ULPS * math.ulp(center) > epsrel * flank:
                    raise DielectricConvergenceError("failed to converge at the resonance flank")
                params = (energy, gap, center, osc.strength_eV2, osc.width_eV, (abs_goal, rel_goal, limit))
                overdamped = osc.width_eV >= 2 * center and center * (center / osc.width_eV) <= gap
                if center < gap or overdamped:
                    value, estimate = _retained_kk_band(*params)
                else:
                    correction, estimate = _removed_kk_band(*params)
                    value = bare_real - correction
                    estimate += SUM_ROUNDOFF * (abs(bare_real) + abs(correction))
            except _PanelLimitError as exc:
                raise ValueError(f"gapped KK oscillator {i}, {location}: {exc}") from exc
            except (ArithmeticError, DielectricConvergenceError) as exc:
                raise DielectricConvergenceError(f"gapped KK oscillator {i}, {location}: {exc}") from exc
            values.append(value)
            estimates.append(estimate)
        try:
            real = math.fsum(values)
            error = math.fsum(estimates) + SUM_ROUNDOFF * math.fsum(map(abs, values))
        except (OverflowError, ValueError) as exc:
            raise DielectricConvergenceError(f"gapped KK oscillator sum outside numerical range at {location}") from exc
        if math.isfinite(real) and math.isfinite(error) and error <= _accuracy(real, epsabs, epsrel):
            return real, error
        # Divide the observed error between all oscillators, using only an
        # absolute target on the next pass. Stop if rounding prevents progress.
        next_goal = error / (KK_REFINEMENT * count)
        if not math.isfinite(next_goal) or next_goal <= 0 or (rel_goal == 0 and next_goal >= abs_goal):
            raise DielectricConvergenceError(f"gapped KK oscillator sum failed to converge at {location}")
        abs_goal, rel_goal = next_goal, 0.0


def _inverse_response(material, energies, q, epsabs, epsrel, limit):
    real = np.full(energies.shape, -1.0)
    error = np.zeros(energies.shape)
    imag = np.zeros(energies.shape)
    gap = material.gap_eV
    bands = []
    magnitude = np.ones(energies.shape)
    active = 0
    with np.errstate(over="ignore", invalid="ignore"):
        for i, osc in enumerate(material.oscillators):
            if not osc.strength_eV2:
                continue
            active += 1
            try:
                bare = _oscillator_response(osc, energies, q, for_real_sum=True)
            except NumericalConvergenceError as exc:
                raise DielectricConvergenceError(f"dielectric oscillator {i}: {exc}") from exc
            imag += bare.imag
            if gap:
                bands.append((i, _oscillator_energy(osc, q), osc, bare.real))
            else:
                real += bare.real
                magnitude += abs(bare.real)
        if bands:
            for index in np.ndindex(energies.shape):
                lines = [(i, float(center[index]), osc, float(bare[index])) for i, center, osc, bare in bands]
                real[index], error[index] = _gapped_real(
                    lines, float(energies[index]), gap, epsabs, epsrel, limit, _loss_point(energies, index))
        result = real + 1j * np.where(energies > gap, imag, 0.0)
    bad = ~np.isfinite(result)
    if not gap:
        bad |= SUM_ROUNDOFF * active * magnitude > _accuracy(real, epsabs, epsrel)
    if np.any(bad):
        raise DielectricConvergenceError(f"dielectric oscillator sum failed to converge at {_first_loss(energies, bad)}")
    return result, error


def _gap_edge_exclusion_eV(material):
    gap = material.gap_eV
    if gap == 0 or all(o.strength_eV2 == 0 for o in material.oscillators):
        return 0.0
    return max(math.ulp(0.0), min(gap / 2, max(GAP_EDGE_ULPS * math.ulp(gap), GAP_EDGE_REL * gap)))


def _check_gap_edge(material, energies):
    gap = material.gap_eV
    exclusion = _gap_edge_exclusion_eV(material)
    if exclusion == 0:
        return
    # The minimum positive gap can itself be one ulp from the regular static point.
    offending = energies[(energies > 0) & (abs(energies - gap) <= exclusion)]
    if offending.size:
        raise DielectricGapEdgeError(gap, exclusion, np.unique(offending))


def _static_response(material, q, epsabs, epsrel, limit):
    response, error = _inverse_response(material, np.zeros_like(q), q, epsabs, epsrel, limit)
    real = response.real
    count = sum(o.strength_eV2 > 0 for o in material.oscillators)
    # Scale from 1 and the positive static KK contribution F(0)+1.
    tolerance = error + (COORDINATE_ULPS * (count + 1) * FLOAT_EPSILON) * (1 + abs(real + 1))
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
    material: Material, loss_eV, q_inv_A=0.0, *, epsabs=DEFAULT_ABS, epsrel=DEFAULT_REL, limit=200,
):
    """Return F=-1/epsilon, with the existing ELF as Im(F) and F(infinity)=-1.

    Loss and q are nonnegative with NumPy broadcasting. At a positive hard
    gap edge, F diverges for any material with nonzero oscillator strength;
    points within the documented exclusion interval raise a gap edge error.
    At the supplied binary64 q, accuracy is checked on Re F itself, including
    close to a dispersed resonance. Quadrature is needed only for a positive gap.
    """
    t, q, epsabs, epsrel = _response_inputs(loss_eV, q_inv_A, epsabs, epsrel, limit)
    _check_gap_edge(material, t)
    response, _ = _inverse_response(material, t, q, epsabs, epsrel, limit)
    return complex(response) if response.ndim == 0 else response


def dielectric_response(
    material: Material, loss_eV, q_inv_A=0.0, *, epsabs=DEFAULT_ABS, epsrel=DEFAULT_REL, limit=200,
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
        raise ValueError("Re[-1/epsilon(q, 0)] exceeds the static sign tolerance: "
                         "epsilon of the oscillator representation has a pole in the upper half-plane")
    marginal = (abs(static) <= tolerance)[q_indices].reshape(t.shape)
    if np.any((t == 0) & marginal):
        raise ValueError("epsilon has a pole at zero energy for a marginal static response")
    response, error = _inverse_response(material, t, q, epsabs, epsrel, limit)
    if np.any(response == 0):
        raise ValueError("epsilon has a pole where -1/epsilon vanishes")
    unresolved = abs(response) <= error
    if np.any(unresolved):
        raise DielectricConvergenceError(f"epsilon unresolved near a zero of -1/epsilon at {_first_loss(t, unresolved)}")
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        epsilon = -1 / response
    bad = ~np.isfinite(epsilon)
    if np.any(bad):
        raise DielectricConvergenceError(f"epsilon oscillator sum outside numerical range at {_first_loss(t, bad)}")
    if response.ndim == 0:
        return DielectricResponse(complex(response), complex(epsilon), float(error))
    return DielectricResponse(response, epsilon, error)


def _near_gap_edge_points(material, q, energies, exclusion, band):
    gap = material.gap_eV
    distances = abs(energies - gap)
    near = (distances > exclusion) & (distances <= band)
    if not np.any(near):
        return ()
    coefficient = loss_function(replace(material, gap_eV=0.0), gap, q) / math.pi
    log_gap = math.log(gap)
    return tuple(
        GapEdgePoint(float(energy), float(distance), coefficient * (log_gap - math.log(distance)))
        for energy, distance in zip(energies[near], distances[near])
    )


def _sub_gap_pole(material, q, exclusion, settings):
    def sub_gap_real(energy):
        value, _ = _inverse_response(material, np.asarray(energy), q, *settings)
        return float(value.real)

    gap = material.gap_eV
    upper = math.nextafter(gap, 0.0)
    if upper == 0 or sub_gap_real(upper) < 0:
        # A zero is guaranteed by the logarithmic divergence, but the interval
        # between the nearest float below the gap and the gap cannot be resolved.
        return None, "unresolved_at_gap_edge"
    pole = brentq(sub_gap_real, 0, upper, xtol=math.ulp(gap), rtol=BRENT_MIN_REL)
    pole_status = "within_gap_edge_exclusion" if gap - pole <= exclusion else "located"
    return pole, pole_status


def dielectric_diagnostics(
    material: Material, q_inv_A=0.0, *, loss_eV=(), epsabs=DEFAULT_ABS, epsrel=DEFAULT_REL, limit=200,
):
    """Classify static F, locate a sub-gap pole and flag optional near-edge samples."""
    _, q, epsabs, epsrel = _response_inputs(0, q_inv_A, epsabs, epsrel, limit)
    if q.ndim:
        raise ValueError("dielectric diagnostics require a single q_inv_A value")
    if np.iscomplexobj(loss_eV):
        raise ValueError("loss must be real")
    energies, _ = _loss_coordinates(loss_eV, 0)
    exclusion = _gap_edge_exclusion_eV(material)
    band = GAP_EDGE_WARNING_FRACTION * material.gap_eV if exclusion else 0.0
    near_edge_points = _near_gap_edge_points(material, q, energies, exclusion, band)
    settings = (epsabs, epsrel, limit)
    static, tolerance = _static_response(material, q, *settings)
    static, tolerance = float(static), float(tolerance)
    pole, pole_status = None, "absent"
    if abs(static) <= tolerance:
        static_sign, pole_status = "marginal", "marginal_static"
    elif static > tolerance:
        static_sign = "positive"
    else:
        static_sign = "negative"
        if exclusion:
            pole, pole_status = _sub_gap_pole(material, q, exclusion, settings)
    return DielectricDiagnostics(static, tolerance, static_sign, exclusion, pole, pole_status,
                                 band, near_edge_points)


def _edge_log(energy, edge):
    """Return log((edge + energy)/abs(edge - energy)) without subtracting logs."""
    return math.log1p(2 * (min(energy, edge) / abs(edge - energy)))


def _closed_kk_band(t, gap, center, strength, width):
    """Return the retained-band integral (Eqs. D2-D3) and its absolute roundoff estimate."""
    half_width = width / 2
    splitting = cmath.sqrt(half_width - center) * cmath.sqrt(half_width + center)
    p = half_width + splitting
    if half_width < center:
        # atan(p/gap) has a logarithmic branch point at p=i*gap.
        # Compute gap-Im(p) without rounding away the width² correction.
        h = splitting.imag
        near = (gap - center) + half_width * (half_width / (center + h))
        far = gap + h
        angle = 0.5 * (math.atan2(half_width, near) + math.atan2(half_width, far))
        # atan(z) = [log(1+iz) - log(1-iz)]/(2i).
        log_ratio = math.log(math.hypot(far, half_width) / math.hypot(near, half_width))
        j_p = complex(angle, 0.5 * log_ratio) / p
        j_b = j_p.conjugate()
    else:
        b = center * (center / p)  # p*b = center**2, also in the overdamped limit
        j_b = cmath.atan(b / gap) / b
        j_p = cmath.atan(p / gap) / p
    difference = 2 * splitting
    if abs(difference) < width * CRITICAL_SPLITTING_FRACTION:
        h = gap / (gap**2 + center**2)
        atan_ratio = cmath.atan(h * difference) / difference if difference else h
        j_divided = (j_b - atan_ratio) / (p * width)
    else:
        j_divided = (j_b - j_p) / (width * difference)

    denominator = ((center - t) * (center + t))**2 + (width * t)**2
    scale = (2 / math.pi) * strength * width / denominator
    regular = (center**4 + (t * p)**2) * j_divided
    log_term = (t / 2) * _edge_log(t, gap)
    pole = log_term - t**2 * j_b
    value = scale * (regular + pole).real
    roundoff_scale = scale * (abs(regular) + log_term + t**2 * abs(j_b))
    resonance_condition = max(1.0, center**2 / math.sqrt(denominator))
    roundoff = RECIPROCAL_ROUNDOFF * resonance_condition * roundoff_scale
    if not math.isfinite(roundoff) or roundoff > _accuracy(value, CLOSED_ABS, CLOSED_REL):
        raise DielectricConvergenceError("closed-form KK is ill-conditioned; use numerical KK")
    return value, roundoff


def negative_inverse_dielectric_closed_form(material: Material, loss_eV, q_inv_A=0.0):
    """Return F=-1/epsilon using an independent closed form for hard-gap KK.

    The retained-band rational integral uses Eqs. D2-D3 in DIELECTRIC_DESIGN.md.
    A zero gap uses the same analytic E4b as the standard API. Both paths
    take Im(F) from the same oscillator ELF.
    The numerical API is preferred where cancellation limits the closed form.
    """
    if material.gap_eV == 0 or all(o.strength_eV2 == 0 for o in material.oscillators):
        return negative_inverse_dielectric_function(material, loss_eV, q_inv_A)
    t, q, _, _ = _response_inputs(loss_eV, q_inv_A, DEFAULT_ABS, DEFAULT_REL, 200)
    gap = material.gap_eV
    _check_gap_edge(material, t)
    real = np.full(t.shape, -1.0)
    magnitude = np.ones(t.shape)
    band_roundoff = np.zeros(t.shape)
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        for i, osc in enumerate(material.oscillators):
            if osc.strength_eV2 == 0:
                continue
            centers = _oscillator_energy(osc, q)
            for index in np.ndindex(t.shape):
                energy, center = float(t[index]), float(centers[index])
                try:
                    if (osc.dispersion and q[index] and COORDINATE_ULPS * math.ulp(center)
                            > CLOSED_CENTER_REL * max(osc.width_eV, abs(gap - center))):
                        raise DielectricConvergenceError("dispersed center unresolved at the gap")
                    value, estimate = _closed_kk_band(energy, gap, center, osc.strength_eV2, osc.width_eV)
                    real[index] += value
                    magnitude[index] += abs(value)
                    band_roundoff[index] += estimate
                except (ArithmeticError, DielectricConvergenceError) as exc:
                    raise DielectricConvergenceError(f"closed-form KK oscillator {i}, {_loss_point(t, index)}: {exc}") from exc
        try:
            result = real + 1j * loss_function(material, t, q)
        except NumericalConvergenceError as exc:
            raise DielectricConvergenceError(f"closed-form KK: {exc}") from exc
    # Opposite bands can cancel in Re F; their absolute error estimates still add.
    rounding = band_roundoff + SUM_ROUNDOFF * len(material.oscillators) * magnitude
    bad = ~(np.isfinite(result) & (rounding <= _accuracy(real)))
    if np.any(bad):
        raise DielectricConvergenceError(
            f"closed-form KK oscillator sum unresolved at {_first_loss(t, bad)}; use numerical KK")
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
    for i, osc in enumerate(material.oscillators):
        if not osc.strength_eV2:
            continue
        center, width = osc.energy_eV, osc.width_eV
        if width < center and gap <= NEARBY_FACTOR * center:
            # z = u - center²/u, theta = atan(z/width) maps the Lorentzian
            # to d(theta). The two weights are (1 ± z/hypot(z,2*center))/2.
            # Integrating the second term in cos(theta) gives asinh.
            z = (gap - center) * (1 + center / gap)
            phi = math.atan2(width, z)
            splitting = math.sqrt(2 * center - width) * math.sqrt(2 * center + width)
            odd = (width / splitting) * math.asinh(splitting / math.hypot(width, z))
            if first_moment:
                value = (osc.strength_eV2 / 2) * (phi + odd)
            else:
                value = (osc.strength_eV2 / center / center / math.pi) * (phi - odd)
            estimate = RECIPROCAL_ROUNDOFF * abs(value)
            if not math.isfinite(value) or estimate > _accuracy(value, settings[0], epsrel):
                raise DielectricConvergenceError(
                    f"optical moment oscillator {i}, loss_eV=[{gap:g}, infinity): roundoff")
            total += value
            continue
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
        try:
            value, error = _quadrature(integrand, 1, sorted(set(points)), settings)
            full = (math.pi / 2 * osc.strength_eV2 if first_moment
                    else osc.strength_eV2 / center / center)
            if (not 0 <= value <= full + error or not math.isfinite(error)
                    or error > _accuracy(value, settings[0], epsrel)):
                raise DielectricConvergenceError("optical moment outside range or error budget")
        except DielectricConvergenceError as exc:
            raise DielectricConvergenceError(
                f"optical moment oscillator {i}, loss_eV=[{gap:g}, infinity): {exc}") from exc
        total += value
    return total


def optical_sum_rule_lhs(material: Material, *, epsabs=DEFAULT_ABS, epsrel=DEFAULT_REL, limit=200):
    """Return the dimensionless retained optical sum-rule integral (Eq. D4)."""
    _, _, epsabs, epsrel = _response_inputs(0, 0, epsabs, epsrel, limit)
    return _checked_optical_moment(material, False, (epsabs, epsrel, limit))


def optical_sum_rule_residual(
    material: Material, reference_inverse_epsilon, *, epsabs=DEFAULT_ABS, epsrel=DEFAULT_REL, limit=200,
):
    """Return the sum-rule LHS minus (1 - reference Re[1/epsilon(0)])."""
    reference = finite_number(reference_inverse_epsilon, "reference_inverse_epsilon")
    return optical_sum_rule_lhs(material, epsabs=epsabs, epsrel=epsrel, limit=limit) - (1 - reference)


def optical_elf_first_moment(material: Material, *, epsabs=DEFAULT_ABS, epsrel=DEFAULT_REL, limit=200):
    """Return the retained ELF first moment in eV**2 (Eq. D4)."""
    _, _, epsabs, epsrel = _response_inputs(0, 0, epsabs, epsrel, limit)
    return _checked_optical_moment(material, True, (epsabs, epsrel, limit))


def _checked_optical_moment(material, first_moment, settings):
    try:
        value = _optical_moment(material, first_moment, settings)
    except ArithmeticError as exc:
        raise DielectricConvergenceError(
            f"optical oscillator integral, loss_eV=[{material.gap_eV:g}, infinity): {exc}") from exc
    if not math.isfinite(value):
        raise DielectricConvergenceError(
            f"nonfinite optical oscillator integral, loss_eV=[{material.gap_eV:g}, infinity)")
    return value


def _power_law_tail(t, edge, level):
    """Return pi times the KK contribution of level*(edge/u)**3 above edge."""
    ratio = t / edge
    if ratio < TAIL_SERIES_RATIO:
        return 2 * level * sum(ratio**(2 * n) / (2 * n + 3) for n in range(5))
    inverse_ratio = edge / t
    return level * inverse_ratio**2 * (inverse_ratio * _edge_log(t, edge) - 2)


def _outside_segments(pole, nodes, levels):
    """Regular segments, oriented from the endpoint nearest the pole."""
    if len(nodes) < 2:
        return 0.0, 0.0
    sign = 1
    if pole > nodes[-1]:
        nodes, levels = nodes[::-1], levels[::-1]
        sign = -1
    ratio = np.diff(nodes) / (nodes[:-1] - pole)
    log = np.log1p(ratio)
    small = ratio < SMALL_SEGMENT_RATIO
    # Weight of the far endpoint: 1 - log1p(r)/r. Use its convergent
    # series when subtracting from one would lose significant digits.
    weight = np.empty_like(ratio)
    r = ratio[small]
    series = np.zeros_like(r)
    for k in range(8, 0, -1):
        series = 1 / (k + 1) - r * series
    weight[small] = r * series
    weight[~small] = 1 - log[~small] / ratio[~small]
    near = (log - weight) * levels[:-1]
    far = weight * levels[1:]
    size = np.sum(abs(near) + abs(far)) + np.sum(abs(np.diff(levels)[~small]))
    return sign * math.fsum(near + far), size


def _linear_cauchy(pole, nodes, values, tail_level=0.0):
    """PV around the pole; return the transform and its roundoff scale."""
    lower, upper = float(nodes[0]), float(nodes[-1])
    if pole == lower and values[0]:
        raise ValueError("tabulated KK diverges at the discontinuous lower cutoff")
    if pole == upper and values[-1] != tail_level:
        raise ValueError("tabulated KK diverges at the discontinuous upper cutoff")
    near_tail = tail_level and upper / NEARBY_FACTOR <= pole <= NEARBY_FACTOR * upper
    if (pole < lower or pole > upper) and not near_tail:
        regular, size = _outside_segments(pole, nodes, values)
        tail = _power_law_tail(pole, upper, tail_level) if tail_level else 0.0
        return regular + tail, size + abs(tail)

    # Subtract the interpolated level throughout the finite band. This makes
    # the absolute integral of a smooth residual bounded under mesh refinement.
    # On each segment touching the pole, the regular integrand is just its slope.
    left = max(0, nodes.searchsorted(pole, "left") - 1)
    right = min(len(nodes) - 1, nodes.searchsorted(pole, "right"))
    if left == right:  # just outside the continuous tail join
        level = values[left]
    else:
        j = nodes.searchsorted(pole, "left")
        if nodes[j] == pole:
            level = values[j]
        else:
            near, far = (left, right) if pole - nodes[left] <= nodes[right] - pole else (right, left)
            fraction = (pole - nodes[near]) / (nodes[far] - nodes[near])
            level = values[near] + fraction * (values[far] - values[near])
    residual = values - level
    a, size_a = _outside_segments(pole, nodes[:left + 1], residual[:left + 1])
    b, size_b = _outside_segments(pole, nodes[right:], residual[right:])
    terms = [a, b, values[right] - values[left]]
    # Restore the integral of the subtracted constant. Collect the log at
    # the cutoff with the tail's log before evaluating a continuous join.
    if level and pole != lower:
        terms.append(level * (math.log(upper) - math.log(abs(lower - pole))))
    coefficient = level
    if near_tail:
        ratio = upper / pole
        tail_weight = tail_level * ratio**3
        coefficient -= tail_weight
        terms.extend((tail_weight * math.log1p(pole / upper), -2 * tail_level * ratio**2))
    elif tail_level:
        terms.append(_power_law_tail(pole, upper, tail_level))
    if pole != upper:
        logarithm = math.log(abs(upper - pole)) - math.log(upper)
        terms.append(coefficient * logarithm)
        if near_tail:
            size_b += (abs(level) + abs(tail_weight)) * abs(logarithm)
    size = size_a + size_b + abs(values[left]) + abs(values[right]) + math.fsum(map(abs, terms[2:]))
    return math.fsum(terms), size


def tabulated_elf_kk(sample_eV, elf, loss_eV, *, low_energy, high_energy):
    """Analytic KK for the piecewise-linear interpolant of one sampled ELF.

    Between samples use linear interpolation. Below the first sample choose
    "zero" or "linear_to_origin". Above the last sample choose "hard_cut"
    or "power_law_3" (continuous u**-3 tail). Losses can be scalars or arrays.
    The analytic quadrature error is zero. An independent roundoff check
    applies the dimensionless component tolerance to the final Re F.
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
        elif abs(terminal_slope + 3) > TAIL_SLOPE_WARNING:
            warnings = ("The terminal log-slope differs from -3 by more than 0.5; check the u^-3 extrapolation.",)
    if low_energy == "linear_to_origin" and lower > 0:
        nodes, values = np.insert(nodes, 0, 0.0), np.insert(values, 0, 0.0)
    real = np.empty(t.shape)
    try:
        tail_coefficient = float(tail_level) * upper**3 if tail_level else 0.0
    except OverflowError as exc:
        raise DielectricConvergenceError(f"tabulated KK tail outside numerical range at cutoff_eV={upper:g}") from exc
    with np.errstate(over="raise", invalid="raise", divide="raise"):
        for index in np.ndindex(t.shape):
            energy = float(t[index])
            try:
                positive, p_size = _linear_cauchy(energy, nodes, values, tail_level)
                negative, n_size = _linear_cauchy(-energy, nodes, values)
                real[index] = -1 + math.fsum((positive, negative)) / math.pi
                # Analytic integration gives zero quadrature error. Test its
                # rounding separately, once, against the final Re F budget.
                rounding = TAB_ROUNDOFF * (1 + (p_size + n_size) / math.pi)
                if not math.isfinite(real[index]) or not rounding <= _accuracy(real[index]):
                    raise ArithmeticError("unresolved segment cancellation")
            except ArithmeticError as exc:
                raise DielectricConvergenceError(f"tabulated KK at {_loss_point(t, index)}: {exc}") from exc
    if not math.isfinite(tail_coefficient):
        raise DielectricConvergenceError(f"tabulated KK tail outside numerical range at cutoff_eV={upper:g}")
    if t.ndim == 0:
        real, error = float(real), 0.0
    else:
        error = np.zeros(t.shape)
    return TabulatedKK(
        real, error, (lower, upper), low_energy, high_energy,
        tail_coefficient if high_energy == "power_law_3" else None,
        terminal_slope, warnings,
    )
