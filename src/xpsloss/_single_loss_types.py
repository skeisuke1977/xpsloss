"""Inputs, results and failure records shared by both loss-window layers."""

from dataclasses import dataclass
from fractions import Fraction
import math

import numpy as np

from .model import NumericalConvergenceError, finite_number
from .xps_depth_averaged import COLUMNS

COMPONENTS = ("total", "no_hole", "hole")
QUANTITIES = ("zero_loss", "single_loss_total", "single_loss_no_hole", "single_loss_hole", "spectrum")
PARTS = ("interpolation", "verification", "kernel_point_error", "edge", "rounding")
EPS = math.ulp(1.)
_ATTRIBUTES = ("required_panels", "panel_eV", "panel_width_ulp", "dominant_part", "node_kind",
               "loss_eV", "index", "column", "stage", "inner_exception", "inner_reason",
               "inner", "inner_attributes")


class _Failure:
    def __init__(self, reason, **attributes):
        self.reason = reason
        for name in _ATTRIBUTES:
            setattr(self, name, None)
        self.__dict__.update(attributes)
        super().__init__(reason)


class SingleLossInputError(_Failure, ValueError):
    pass


class SingleLossConvergenceError(_Failure, NumericalConvergenceError):
    pass


def vector(value, reason, *, increasing=False):
    try:
        if isinstance(value, np.ndarray) and value.dtype.kind in "fiu":
            arr = value
        else:
            # Avoid coercing mixed booleans to integers, and accept large
            # JSON integers such as 2**900 when representable in binary64.
            items = np.asarray(value, dtype=object)
            if items.ndim != 1:
                raise ValueError
            arr = np.asarray([finite_number(x, reason) for x in items])
        if arr.ndim != 1 or not arr.size:
            raise ValueError
        arr = np.array(arr, dtype=np.float64)
        if not np.isfinite(arr).all() or np.any(arr < 0):
            raise ValueError
        if increasing and (len(arr) < 2 or np.any(arr[:-1] >= arr[1:])):
            raise ValueError
    except (TypeError, ValueError, OverflowError) as exc:
        raise SingleLossInputError(reason) from exc
    arr.setflags(write=False)
    return arr


def _primary(obj, field, binned):
    reason = "invalid_primary"
    energies = vector(getattr(obj, field), reason, increasing=binned)
    intensities = vector(obj.intensities, reason)
    if (len(energies) != len(intensities) + binned
            or np.any((intensities > 0) & (intensities < 2.**-900))
            or sum(map(Fraction, intensities)) > 2**900):
        raise SingleLossInputError(reason)
    for key in ("intensity_unit", "description", "included_broadening", "truncation", "intrinsic_satellites"):
        value = getattr(obj, key)
        if value is None and key != "intensity_unit":
            continue
        if not isinstance(value, str) or (key == "intensity_unit" and not value.strip()):
            raise SingleLossInputError(reason)
        try:
            value.encode("utf-8")
        except UnicodeError as exc:
            raise SingleLossInputError(reason) from exc
    object.__setattr__(obj, field, energies)
    object.__setattr__(obj, "intensities", intensities)


@dataclass(frozen=True)
class DeltaPrimary:
    energies_eV: np.ndarray
    intensities: np.ndarray
    intensity_unit: str = "arb. unit"
    description: str | None = None
    included_broadening: str | None = None
    truncation: str | None = None
    intrinsic_satellites: str | None = None

    def __post_init__(self):
        _primary(self, "energies_eV", False)


@dataclass(frozen=True)
class BinnedPrimary:
    edges_eV: np.ndarray
    intensities: np.ndarray
    intensity_unit: str = "arb. unit"
    description: str | None = None
    included_broadening: str | None = None
    truncation: str | None = None
    intrinsic_satellites: str | None = None

    def __post_init__(self):
        _primary(self, "edges_eV", True)


@dataclass
class LossWindowAreaResult:
    electron_energy_eV: float | None
    loss_window_eV: tuple
    weight_length_A: float
    epsilon_momentum: str | None
    hole_charge: float
    theta_deg: float
    window_area: dict
    window_area_error: dict
    window_area_budget: dict
    single_loss_weight: float
    single_loss_weight_error: float
    inverse_window_area_A: float | None
    inverse_window_area_interval_A: tuple | None
    inverse_reason: str | None
    cutoff_sensitivity: dict
    cutoff_sensitivity_error: dict
    diagnostics: dict
    warnings: tuple


@dataclass
class SingleLossSpectrumResult(LossWindowAreaResult):
    output_edges_eV: np.ndarray
    primary: DeltaPrimary | BinnedPrimary
    primary_integral: float
    primary_abs_integral: float
    zero_loss: np.ndarray
    zero_loss_error: np.ndarray
    single_loss_total: np.ndarray
    single_loss_total_error: np.ndarray
    single_loss_no_hole: np.ndarray
    single_loss_no_hole_error: np.ndarray
    single_loss_hole: np.ndarray
    single_loss_hole_error: np.ndarray
    spectrum: np.ndarray
    spectrum_error: np.ndarray
    zero_loss_per_eV: np.ndarray
    zero_loss_per_eV_error: np.ndarray
    single_loss_total_per_eV: np.ndarray
    single_loss_total_per_eV_error: np.ndarray
    single_loss_no_hole_per_eV: np.ndarray
    single_loss_no_hole_per_eV_error: np.ndarray
    single_loss_hole_per_eV: np.ndarray
    single_loss_hole_per_eV_error: np.ndarray
    spectrum_per_eV: np.ndarray
    spectrum_per_eV_error: np.ndarray
    cutoff_sensitivity_bins: dict
    cutoff_sensitivity_bins_error: dict
    outside: dict


def inverse_area(area, error):
    a, e = Fraction(area), Fraction(error)
    if a <= e:
        return None, None, "nonpositive_or_unresolved_area"
    try:
        lower, upper = 1/(a+e), 1/(a-e)
        lo, hi = float(lower), float(upper)
        if lo > lower:
            lo = math.nextafter(lo, -math.inf)
        if hi < upper:
            hi = math.nextafter(hi, math.inf)
        if math.isfinite(hi):
            return float(1/a), (lo, hi), None
    except OverflowError:
        pass
    return None, None, "inverse_not_representable"

