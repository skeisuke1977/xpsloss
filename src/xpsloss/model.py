"""Elementary dielectric bulk losses, with all energies in eV and lengths in Å.

ELF: Pauly, Yubero & Tougaard (2023), doi:10.3390/nano13020339, eqs. 3a/3b.
Bulk DIIMFP: Tougaard et al. (2022), doi:10.1002/sia.7095, eq. 12;
also QUEELS-XPS v2.1 manual, p. 10, eq. 2 (infinite medium).
No surface, core-hole, interference, depth averaging, or multiple scattering.
The finite-q oscillator dispersion is a model assumption, not measured data.
"""

from dataclasses import dataclass
import math
import warnings

import numpy as np
from scipy.integrate import IntegrationWarning, quad


# Rounded physical constants (angstrom/eV convention), adequate to this
# starter's numerical precision; not a claim of metrological exactness.
BOHR_RADIUS_A = 0.529177210903
HBAR2_OVER_2ME_EV_A2 = 3.80998212
MODEL_ID = "bulk_drude_lindhard_v0"
UNSUPPORTED_PHYSICS = (
    "surface_losses", "core_hole", "intrinsic_extrinsic_interference",
    "depth_averaging", "elastic_scattering", "multiple_inelastic_scattering",
    "primary_spectrum_convolution", "full_imfp", "dielectric_real_part",
)


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


def loss_function(material: Material, loss_eV, q_inv_A=0.0):
    """Return dimensionless Im[-1/epsilon(q,T)] with NumPy broadcasting.

    The gap convention is L=0 for T <= Eg, including the exact threshold.
    This function does not reconstruct the real part of epsilon or enforce
    optical sum rules. Values at T=0 are set to zero.
    """
    t, q = np.broadcast_arrays(np.asarray(loss_eV, dtype=float), np.asarray(q_inv_A, dtype=float))
    if not np.all(np.isfinite(t)) or not np.all(np.isfinite(q)):
        raise ValueError("loss and q must be finite")
    if np.any(t < 0) or np.any(q < 0):
        raise ValueError("loss and q must be nonnegative")
    result = np.zeros(t.shape, dtype=float)
    try:
        with np.errstate(over="raise", invalid="raise", divide="raise"):
            for osc in material.oscillators:
                resonance = osc.energy_eV + osc.dispersion * HBAR2_OVER_2ME_EV_A2 * q**2
                denominator = (resonance**2 - t**2)**2 + (osc.width_eV * t)**2
                result += osc.strength_eV2 * osc.width_eV * t / denominator
    except FloatingPointError as exc:
        raise ValueError("numeric overflow in ELF; check parameter scale") from exc
    result = np.where(t > material.gap_eV, result, 0.0)
    if not np.all(np.isfinite(result)):
        raise ValueError("nonfinite ELF result")
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


def bulk_diimfp(material: Material, electron_energy_eV: float, loss_eV: float,
                *, epsabs: float = 1e-9, epsrel: float = 1e-7, limit: int = 200):
    """Return (DIIMFP, estimated quadrature error), both in Å^-1 eV^-1.

    Integrate L(q,T) d(log q)/(pi*a0*E) with adaptive quadrature.
    The numerical error estimate excludes physical model/input uncertainty.
    """
    e = finite_number(electron_energy_eV, "electron_energy_eV")
    t = finite_number(loss_eV, "loss_eV")
    qmin, qmax = momentum_limits(e, t)
    epsabs = finite_number(epsabs, "epsabs")
    epsrel = finite_number(epsrel, "epsrel")
    if epsabs <= 0 or epsrel <= 0:
        raise ValueError("epsabs and epsrel must be positive")
    if isinstance(limit, bool) or not isinstance(limit, int) or limit < 1:
        raise ValueError("limit must be a positive integer")
    if t <= material.gap_eV or all(o.strength_eV2 == 0 for o in material.oscillators):
        return 0.0, 0.0
    # Resolve any oscillator resonance that lies within the allowed q range.
    points = []
    for osc in material.oscillators:
        if osc.dispersion > 0 and t > osc.energy_eV:
            qr = math.sqrt((t - osc.energy_eV) / (osc.dispersion * HBAR2_OVER_2ME_EV_A2))
            if qmin < qr < qmax:
                points.append(math.log(qr))
    with warnings.catch_warnings():
        warnings.simplefilter("error", IntegrationWarning)
        try:
            val, err = quad(
                lambda u: loss_function(material, t, math.exp(u)),
                math.log(qmin), math.log(qmax), points=sorted(set(points)) or None,
                epsabs=epsabs, epsrel=epsrel, limit=limit,
            )
        except IntegrationWarning as exc:
            raise ValueError(f"bulk integration failed to converge: {exc}") from exc
    prefactor = 1.0 / (math.pi * BOHR_RADIUS_A * e)
    val, err = val * prefactor, err * prefactor
    if not math.isfinite(val) or not math.isfinite(err) or val < 0:
        raise ValueError("invalid bulk integration result")
    return val, err
