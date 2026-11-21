"""Elementary dielectric bulk losses, with all energies in eV and lengths in Å.

ELF: Pauly, Yubero & Tougaard (2023), doi:10.3390/nano13020339, eqs. 3a/3b.
Bulk DIIMFP: Tougaard et al. (2022), doi:10.1002/sia.7095, eq. 12;
also QUEELS-XPS v2.1 manual, p. 10, eq. 2 (infinite medium).
No surface, core-hole, interference, depth averaging, or multiple scattering.
The finite-q oscillator dispersion is a model assumption, not measured data.
"""

from dataclasses import dataclass
from decimal import Decimal, localcontext
from fractions import Fraction
import math

import numpy as np
from scipy.integrate import quad


# Rounded physical constants (angstrom/eV convention), adequate to this
# starter's numerical precision; not a claim of metrological exactness.
BOHR_RADIUS_A = 0.529177210903
HBAR2_OVER_2ME_EV_A2 = 3.80998212
MODEL_ID = "bulk_drude_lindhard_v0"

FLOAT_EPSILON = np.finfo(float).eps
_LOG_MIN_NORMAL = math.log(np.finfo(float).tiny)
DEFAULT_ABS = 1e-9  # Absolute floor for dimensionless response components.
DEFAULT_REL = 1e-7  # Relative accuracy of the published numerical examples.
SUM_ROUNDOFF = 4 * FLOAT_EPSILON  # Allow for additions and the normalization of a sum.
COORDINATE_ULPS = 8  # Products, squares and endpoint arithmetic in a dispersed coordinate.
RECIPROCAL_ROUNDOFF = 16 * FLOAT_EPSILON  # Products and complex division in an oscillator response.
DETUNING_REL = 1e-10  # Relative detuning target for an individual line in ELF and fixed-depth integrands.
DECIMAL_DIGITS = 80  # Guard digits for differences of nearly equal binary64 energies.
NEARBY_FACTOR = 2  # Within a factor of two, switch coordinates to resolve nearby features.
BULK_EDGE_FRACTION = 1e-5  # Reconstruct the endpoint when subtraction loses five decimal digits.
LOG_PANEL_SPAN = 4  # Sample at least once per four units of log(s) or asinh(detuning/halfwidth).
FLANK_RATIO = 8  # Geometric shoulders resolve a linewidth and its tails with bounded panel count.


def _accuracy(value, epsabs=DEFAULT_ABS, epsrel=DEFAULT_REL):
    """Absolute error budget of the reported quantity, including near zero."""
    return np.maximum(epsabs, epsrel * np.abs(value))


def _loss_point(loss, index):
    index = tuple(map(int, index))
    return f"loss_eV={float(loss[index]):.17g}, index={index}"


def _first_loss(loss, mask):
    return _loss_point(loss, np.unravel_index(np.flatnonzero(mask)[0], loss.shape))


# Model-specific scopes. Existing bulk runs retain their original exclusions.
UNSUPPORTED_PHYSICS = (
    "surface_losses", "core_hole", "intrinsic_extrinsic_interference",
    "depth_averaging", "elastic_scattering", "multiple_inelastic_scattering",
    "primary_spectrum_convolution", "full_imfp",
)
FIXED_DEPTH_UNSUPPORTED_PHYSICS = (
    "oblique_emission", "gapped_materials", "depth_averaging", "elastic_scattering",
    "multiple_inelastic_scattering", "primary_spectrum_convolution", "full_imfp",
)


class _PanelLimitError(ValueError):
    """The specified integration limit cannot hold the initial partition."""


class NumericalConvergenceError(RuntimeError):
    """A valid numerical calculation could not be resolved."""


class BulkConvergenceError(NumericalConvergenceError):
    """Momentum integration could not resolve a finite-width oscillator."""


def finite_number(value, label):
    if isinstance(value, bool) or not isinstance(value, (int, float, np.number)):
        raise ValueError(f"{label} must be a finite number")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{label} must be finite")
    return result


@dataclass(frozen=True)
class Oscillator:
    energy_eV: float
    strength_eV2: float
    width_eV: float
    dispersion: float = 0.0

    def __post_init__(self):
        for key in ("energy_eV", "strength_eV2", "width_eV", "dispersion"):
            object.__setattr__(self, key, finite_number(getattr(self, key), key))
        if self.energy_eV <= 0 or self.width_eV <= 0:
            raise ValueError("oscillator energy_eV and width_eV must be positive")
        if self.strength_eV2 < 0 or self.dispersion < 0:
            raise ValueError("oscillator strength_eV2 and dispersion must be nonnegative")


@dataclass(frozen=True)
class Material:
    name: str
    gap_eV: float
    oscillators: tuple[Oscillator, ...]

    def __post_init__(self):
        if not isinstance(self.name, str) or not self.name.strip():
            raise ValueError("material name must be a nonempty string")
        object.__setattr__(self, "gap_eV", finite_number(self.gap_eV, "gap_eV"))
        if self.gap_eV < 0:
            raise ValueError("gap_eV must be nonnegative")
        if not isinstance(self.oscillators, (tuple, list)) or not self.oscillators:
            raise ValueError("at least one oscillator is required")
        if not all(isinstance(o, Oscillator) for o in self.oscillators):
            raise ValueError("oscillators must contain Oscillator objects")
        object.__setattr__(self, "oscillators", tuple(self.oscillators))


def _loss_coordinates(loss_eV, q_inv_A):
    t, q = np.broadcast_arrays(np.asarray(loss_eV, dtype=float), np.asarray(q_inv_A, dtype=float))
    if not np.all(np.isfinite(t)) or not np.all(np.isfinite(q)):
        raise ValueError("loss and q must be finite")
    if np.any(t < 0) or np.any(q < 0):
        raise ValueError("loss and q must be nonnegative")
    return t, q


def _oscillator_energy(oscillator, q):
    return oscillator.energy_eV + oscillator.dispersion * HBAR2_OVER_2ME_EV_A2 * q**2


def _exact_oscillator_sum(material, loss, q, *, q_squared=None):
    """Exact rational S at binary64 q or at a reconstructed quadrature node."""
    t = Fraction(float(loss))
    recoil = Fraction(HBAR2_OVER_2ME_EV_A2) * (Fraction(float(q))**2 if q_squared is None else q_squared)
    real = imag = Fraction(0)
    for osc in material.oscillators:
        if osc.strength_eV2:
            center = Fraction(osc.energy_eV) + Fraction(osc.dispersion) * recoil
            delta = (center - t) * (center + t)
            damping = Fraction(osc.width_eV) * t
            weight = Fraction(osc.strength_eV2) / (delta**2 + damping**2)
            real += weight * delta
            imag += weight * damping
    return real, imag


def _oscillator_response(osc, t, q, *, for_real_sum=False, with_detuning_error=False):
    """E4b, using the same factored detuning for both parts of the reciprocal."""
    with np.errstate(over="ignore", divide="ignore", invalid="ignore", under="ignore"):
        shift = osc.dispersion * HBAR2_OVER_2ME_EV_A2 * q**2 if osc.dispersion else 0.0
        distance = (osc.energy_eV - t) + shift
        if osc.dispersion:
            scale, rtol = shift, DETUNING_REL
            if for_real_sum:
                # In shift - (T-W) only the cancelled part amplifies roundoff.
                # Other lines can cancel Re F, so allow only the sum's roundoff.
                scale = np.minimum(shift, np.maximum(t - osc.energy_eV, 0.0))
                rtol = SUM_ROUNDOFF
            uncertain = ((scale > 0) & np.isfinite(shift)
                         & (RECIPROCAL_ROUNDOFF * scale >= rtol * abs(distance)))
            if np.any(uncertain):
                distance = np.array(distance, copy=True)
                with localcontext() as ctx:
                    ctx.prec = DECIMAL_DIGITS
                    slope = Decimal(osc.dispersion) * Decimal(HBAR2_OVER_2ME_EV_A2)
                    for row in np.argwhere(uncertain):
                        i = tuple(row)
                        distance[i] = float(Decimal(osc.energy_eV) - Decimal(float(t[i]))
                                            + slope * Decimal(float(q[i]))**2)
        delta = distance * ((osc.energy_eV + shift) + t)
        damping = osc.width_eV * t
        response = osc.strength_eV2 / (delta - 1j * damping)
        if with_detuning_error and osc.dispersion:
            # Only the cancelled part of (W-T)+shift needs an extra
            # allowance. Decimal-reconstructed distances are rounded once.
            cancelled = np.minimum(shift, np.maximum(t - osc.energy_eV, 0.0))
            delta_error = (RECIPROCAL_ROUNDOFF * cancelled) * ((osc.energy_eV + shift) + t)
            relative = np.where(uncertain, 0.0, delta_error / np.hypot(delta, damping))
    bad = ~(np.isfinite(delta) & np.isfinite(damping) & np.isfinite(response))
    if np.any(bad):
        raise NumericalConvergenceError(f"oscillator response outside numerical range at {_first_loss(t, bad)}")
    if with_detuning_error:
        return response, relative if osc.dispersion else 0.0
    return response


def loss_function(material: Material, loss_eV, q_inv_A=0.0):
    """Return dimensionless Im[-1/epsilon(q,T)] with NumPy broadcasting.

    The gap convention is L=0 for T <= Eg, including the exact threshold.
    This function does not reconstruct the real part of epsilon or enforce
    optical sum rules. Values at T=0 are set to zero.
    """
    t, q = _loss_coordinates(loss_eV, q_inv_A)
    result = np.zeros(t.shape, dtype=float)
    with np.errstate(over="ignore", invalid="ignore"):
        for i, osc in enumerate(material.oscillators):
            if osc.strength_eV2:
                try:
                    result += _oscillator_response(osc, t, q).imag
                except NumericalConvergenceError as exc:
                    raise NumericalConvergenceError(f"ELF oscillator {i}: {exc}") from exc
    result = np.where(t > material.gap_eV, result, 0.0)
    bad = ~np.isfinite(result)
    if np.any(bad):
        raise NumericalConvergenceError(f"nonfinite ELF oscillator sum at {_first_loss(t, bad)}")
    return float(result) if result.ndim == 0 else result


def momentum_limits(electron_energy_eV: float, loss_eV: float):
    """Kinematic q bounds in inverse Å; require 0 < T < E.

    q_minus is rationalized to avoid cancellation when T/E is small.
    These are direct dielectric-scattering bounds; exchange corrections and
    full-IMFP integration limits are outside this starter's scope.
    """
    e = finite_number(electron_energy_eV, "electron_energy_eV")
    t = finite_number(loss_eV, "loss_eV")
    if not 0 < t < e:
        raise ValueError("require 0 < loss_eV < electron_energy_eV")
    s = math.sqrt(e) + math.sqrt(e - t)
    root_h = math.sqrt(HBAR2_OVER_2ME_EV_A2)
    qmin, qmax = t / (s * root_h), s / root_h
    if not 0 < qmin < qmax or not math.isfinite(qmax):
        raise ValueError("momentum limits outside numerical range")
    return qmin, qmax


def _edge_detuning(e, t, osc, upper, shift, peak):
    distance = shift - peak
    if abs(distance) > BULK_EDGE_FRACTION * shift:
        return distance
    # The rounded kinematic limit can otherwise move by a substantial fraction
    # of a width. Form the difference from the actual binary64 model inputs.
    with localcontext() as ctx:
        ctx.prec = DECIMAL_DIGITS
        energy, loss = Decimal(e), Decimal(t)
        plus = energy.sqrt() + (energy - loss).sqrt()
        bound = plus if upper else loss / plus
        return float(Decimal(osc.energy_eV) - loss + Decimal(osc.dispersion) * bound**2)


def _bulk_line(osc, e, t, epsabs, epsrel, limit):
    plus = math.sqrt(e) + math.sqrt(e - t)
    damping = osc.width_eV * t
    if not 0 < damping < math.inf:
        raise ArithmeticError("damping outside numerical range")
    peak = t - osc.energy_eV

    def log_momentum_ratio():
        ratio = math.sqrt(e - t) / math.sqrt(t)
        if math.isfinite(ratio):
            return 2 * math.asinh(ratio)
        # The log difference is well-conditioned for such a large ratio.
        return 2 * math.log(plus) - math.log(t)

    if not osc.dispersion:
        delta = -peak * (osc.energy_eV + t)
        value = (osc.strength_eV2 / complex(delta, -damping)).imag * log_momentum_ratio()
        return value, COORDINATE_ULPS * math.ulp(value)

    # s = alpha H q²; dq/q = ds/(2s). Away from the peak log(s) resolves
    # the kinematic lower limit. In [peak/2, 2*peak] the denominator s is
    # smooth, and d = s-peak = (Gamma/2) sinh(v) resolves both flanks even
    # when adding a half-width to the optical energy would round back to it.
    lo = osc.dispersion * (t / plus)**2
    hi = osc.dispersion * plus**2
    if not 0 < lo < hi < math.inf:
        raise ArithmeticError("unresolved kinematic bounds")
    width = osc.width_eV / 2
    parts = []

    def elf(detuning, shift):
        delta = detuning * ((osc.energy_eV + t) + shift)
        return (osc.strength_eV2 / complex(delta, -damping)).imag

    def log_panel(left, right):
        if left >= right:
            return
        if left == lo and right == hi:
            # log(s_hi/s_lo) = 4 asinh(sqrt((E-T)/T)). Keep this span
            # independent of rounding of alpha and of the two endpoints.
            span = 2 * log_momentum_ratio()
        else:
            ratio = (right - left) / left
            # If the ratio overflows, the log endpoints are at least 709
            # apart and their difference is well-conditioned.
            span = math.log1p(ratio) if math.isfinite(ratio) else math.log(right) - math.log(left)

        # Anchor wide intervals at the upper end: exp(u) could overflow
        # even when the shifted energy itself is representable.
        def from_upper(u):
            shift = right * math.exp(u) if u >= _LOG_MIN_NORMAL else math.exp(math.log(right) + u)
            return 0.5 * elf(shift - peak, shift)

        parts.append((from_upper, -span, 0.0))

    left, right = max(lo, peak / NEARBY_FACTOR), min(hi, NEARBY_FACTOR * peak)
    if left < right:
        log_panel(lo, left)
        dl = _edge_detuning(e, t, osc, False, lo, peak) if left == lo else left - peak
        dr = _edge_detuning(e, t, osc, True, hi, peak) if right == hi else right - peak

        def centered(v):
            d = width * math.sinh(v)
            shift = peak + d
            return 0.5 * elf(d, shift) * (width * math.cosh(v) / shift)

        parts.append((centered, math.asinh(dl / width), math.asinh(dr / width)))
        log_panel(right, hi)
    else:
        log_panel(lo, hi)

    values, errors = [], []
    # Reserve the initial log panels for every region before adaptive subdivision.
    counts = [max(1, math.ceil((upper - lower) / LOG_PANEL_SPAN)) for _, lower, upper in parts]
    reserved = sum(counts)
    if reserved > limit:
        raise _PanelLimitError(f"integration limit too small for {reserved} initial oscillator panels")
    remaining = limit
    for (function, lower, upper), panels in zip(parts, counts):
        reserved -= panels
        if upper <= lower:
            raise ArithmeticError("unresolved logarithmic interval")
        points = np.linspace(lower, upper, panels + 1)[1:-1]
        # Summing max(abs_i, rel_i) for separate positive integrals may
        # exceed max(abs_total, rel_total). Reserve half for each budget.
        integral = quad(function, lower, upper, points=points, epsabs=epsabs / (2 * len(parts)),
                        epsrel=epsrel / 2, limit=remaining - reserved, full_output=True)
        if len(integral) > 3:
            raise ArithmeticError(integral[3])
        values.append(integral[0])
        errors.append(integral[1])
        remaining -= integral[2]["last"]
    return math.fsum(values), math.fsum(errors)


def bulk_diimfp(material: Material, electron_energy_eV: float, loss_eV: float,
                *, epsabs: float = DEFAULT_ABS, epsrel: float = DEFAULT_REL, limit: int = 200):
    """Return (DIIMFP, estimated quadrature error), both in Å^-1 eV^-1.

    Integrate each oscillator in L(q,T) d(log q)/(pi*a0*E) separately.
    The numerical error estimate excludes physical model/input uncertainty.
    """
    e = finite_number(electron_energy_eV, "electron_energy_eV")
    t = finite_number(loss_eV, "loss_eV")
    if not 0 < t < e:
        raise ValueError("require 0 < loss_eV < electron_energy_eV")
    epsabs = finite_number(epsabs, "epsabs")
    epsrel = finite_number(epsrel, "epsrel")
    if epsabs <= 0 or epsrel <= 0:
        raise ValueError("epsabs and epsrel must be positive")
    if isinstance(limit, bool) or not isinstance(limit, int) or limit < 1:
        raise ValueError("limit must be a positive integer")
    if t <= material.gap_eV or all(o.strength_eV2 == 0 for o in material.oscillators):
        return 0.0, 0.0
    active = sum(o.strength_eV2 > 0 for o in material.oscillators)
    values, errors = [], []
    for i, osc in enumerate(material.oscillators):
        if not osc.strength_eV2:
            continue
        try:
            val, err = _bulk_line(osc, e, t, epsabs / active, epsrel, limit)
            if not (math.isfinite(val) and val > 0 and math.isfinite(err) and err >= 0
                    and err <= _accuracy(val, epsabs / active, epsrel)):
                raise ArithmeticError("integral or error estimate outside tolerance")
        except _PanelLimitError as exc:
            raise ValueError(f"bulk oscillator {i}, loss_eV={t:g}: {exc}") from exc
        except (ArithmeticError, ValueError) as exc:
            raise BulkConvergenceError(f"bulk oscillator {i}, loss_eV={t:g}: {exc}") from exc
        values.append(val)
        errors.append(err)
    try:
        val, err = math.fsum(values), math.fsum(errors)
        if err > _accuracy(val, epsabs, epsrel):
            raise ArithmeticError("combined quadrature error exceeds tolerance")
        result = (val / e / (math.pi * BOHR_RADIUS_A), err / e / (math.pi * BOHR_RADIUS_A))
        if result[0] <= 0 or not all(map(math.isfinite, result)):
            raise ArithmeticError("normalization outside numerical range")
    except ArithmeticError as exc:
        raise BulkConvergenceError(f"bulk oscillator sum, loss_eV={t:g}: {exc}") from exc
    return result
