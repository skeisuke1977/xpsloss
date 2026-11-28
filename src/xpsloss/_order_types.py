"""Records and preflight checks for the finite-window loss-order expansion."""

from dataclasses import dataclass
from fractions import Fraction as Q
from numbers import Integral

import numpy as np

from ._fold import Geometry
from ._single_loss_types import BinnedPrimary, DeltaPrimary, LossWindowAreaResult, SingleLossInputError, vector
from .model import NumericalConvergenceError, finite_number

MULTIPLE_LOSS_MODEL_ID = "first_loss_13a_chain_v0"
PARTS = ("inherited", "propagated", "local")
ORDER_DEFAULTS = dict(
    order_epsrel=1e-6,
    order_epsabs_fraction=1e-9,
    order_limit=64,
    max_convolution_points=200_000_000,
    neglect_tolerance_fraction=None,
)


class _OrderFailure:
    def __init__(self, reason, **details):
        self.reason = reason
        for field in (
            "argument requested_order order_limit resource lower_bound limit order stage inner estimate "
            "measured exceeded key side column suggestion required_order_estimate panel_eV panel_limit "
            "diagnostics"
        ).split():
            setattr(self, field, details.get(field))
        if reason == "invalid_energy_axis":
            self.field, self.detail = details.get("field"), details.get("detail")
        elif reason == "axis_conversion_not_injective":
            self.binding_eV, self.kinetic_eV, self.roles = (details.get(k) for k in ("binding_eV", "kinetic_eV", "roles"))
        message = reason
        if self.suggestion:
            message += (
                f": check the first-order mesh and absolute kernel accuracy for {self.key}, "
                "and try lowering epsrel (for example, epsrel=1e-7)"
            )
            if self.suggestion != "epsrel":
                message += f". Also check the first-order absolute tolerance {self.suggestion}"
        elif details.get("detail"):
            message += f" ({details['detail']})"
        super().__init__(message)


class MultipleLossInputError(_OrderFailure, ValueError):
    pass


class MultipleLossConvergenceError(_OrderFailure, NumericalConvergenceError):
    pass


@dataclass
class OrderAreas:
    order: int
    area: dict
    area_error: dict
    area_error_parts: dict
    area_zero_kind: dict
    cutoff_sensitivity: dict
    cutoff_sensitivity_error: dict
    mesh: dict


@dataclass
class SpectrumOrder(OrderAreas):
    total: np.ndarray
    total_error: np.ndarray
    total_error_parts: dict
    no_hole: np.ndarray
    no_hole_error: np.ndarray
    no_hole_error_parts: dict
    hole: np.ndarray
    hole_error: np.ndarray
    hole_error_parts: dict
    outside: dict


@dataclass
class MultipleLossOrderAreasResult:
    single_loss: LossWindowAreaResult
    orders: tuple
    order_area_ratio: dict
    convention_ratio: str
    diagnostics: dict
    warnings: tuple
    model_id: str
    depends_on_model_id: str
    kernel_energy_policy: str = "fixed_E0"
    repeated_kernel: str = "term_13a"
    hole_treatment: str = "first_event_only"


@dataclass(kw_only=True)
class MultipleLossSpectrumResult(MultipleLossOrderAreasResult):
    """Displayed bins follow the input axis. single_loss always retains the ascending KE bins."""

    zero_loss: np.ndarray
    zero_loss_error: np.ndarray
    spectrum: np.ndarray
    spectrum_error: np.ndarray
    outside: dict
    truncation: dict
    energy_axis: dict | None = None


def positive_integer(value):
    return isinstance(value, Integral) and not isinstance(value, (bool, np.bool_)) and value >= 1


def positive_real(value):
    try:
        return finite_number(value, "setting") > 0
    except (ValueError, OverflowError):
        return False


def preflight(model, requested, window, primary, edges, settings, zero=False, display=None):
    geometry = None
    if primary is not None:
        if not isinstance(primary, (BinnedPrimary, DeltaPrimary)):
            raise SingleLossInputError("invalid_primary")
        geometry = Geometry(primary, vector(edges, "invalid_output_edges", increasing=True))
    return plan_orders(model, requested, window, geometry, settings, zero, display)


def plan_orders(model, requested, window, geometry, settings, zero, display=None):
    if not isinstance(model, str) or model != MULTIPLE_LOSS_MODEL_ID:
        raise MultipleLossInputError("invalid_model")
    by_support = isinstance(requested, str) and requested == "display_support" and geometry is not None
    if not by_support and not positive_integer(requested):
        raise MultipleLossInputError("invalid_max_order")
    tol = settings["neglect_tolerance_fraction"]
    if tol is not None and not positive_real(tol):
        raise MultipleLossInputError("invalid_neglect_tolerance")
    for key in ("order_limit", "max_convolution_points", "order_epsrel", "order_epsabs_fraction"):
        check = positive_integer if key in ("order_limit", "max_convolution_points") else positive_real
        if not check(settings[key]):
            raise MultipleLossInputError("invalid_settings", argument=key)
    # Fraction accepts Python binary64, but not every NumPy floating scalar.
    # Preserve the accepted numerical settings when forming rational budgets.
    for key, value in settings.items():
        if isinstance(value, np.floating):
            settings[key] = float(value)
        elif isinstance(value, np.integer):
            settings[key] = int(value)
    support = None
    if display is not None:
        support = display.support_order(window[0])
    elif geometry is not None:
        support = max(0, (Q(geometry.points[-1]) - Q(geometry.edges[0])) // Q(window[0])) if geometry.points else 0
    n = support if by_support else int(requested)
    if n > settings["order_limit"]:
        raise MultipleLossInputError(
            "order_limit_exceeded", requested_order=n, order_limit=settings["order_limit"]
        )
    lower = 0 if zero or n < 2 else 527 * (n * (n + 1) // 2 - 1)
    if lower > settings["max_convolution_points"]:
        raise MultipleLossInputError(
            "work_lower_bound_exceeded",
            lower_bound=lower,
            resource="max_convolution_points",
            limit=settings["max_convolution_points"],
        )
    if geometry and not zero and positive_integer(settings["max_cells"]):
        regions = len(geometry.names)
        first = settings["max_cells"] // regions + 1
        if first <= n:
            raise MultipleLossInputError(
                "work_lower_bound_exceeded",
                order=first,
                lower_bound=regions * first,
                resource="max_cells",
                limit=settings["max_cells"],
            )
    return n, support, geometry, lower
