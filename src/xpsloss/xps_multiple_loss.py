"""Approximate normal-emission XPS loss orders on a finite single-event window.

The first event uses the complete depth-averaged kernel at E0. Each later event
uses only term 13a at the same E0. The multiplier is constant across events.
"""

from . import _loss_orders, xps_single_loss as single
from ._loss_panels import checked_window, observe_mesh
from ._order_types import (
    MULTIPLE_LOSS_MODEL_ID,
    ORDER_DEFAULTS,
    MultipleLossConvergenceError,
    MultipleLossInputError,
    MultipleLossOrderAreasResult,
    MultipleLossSpectrumResult,
    preflight,
)
from ._single_loss_types import BinnedPrimary, DeltaPrimary, SingleLossInputError
from .model import Material
from .xps_fixed_depth import FIXED_ABS, FIXED_REL

__all__ = [
    "MULTIPLE_LOSS_MODEL_ID",
    "MultipleLossInputError",
    "MultipleLossConvergenceError",
    "MultipleLossSpectrumResult",
    "MultipleLossOrderAreasResult",
    "multiple_loss_spectrum",
    "multiple_loss_order_areas",
]


def _run_orders(material, electron_energy_eV, primary, edges, args, defaults):
    first_settings = {k: args[k] for k in defaults}
    settings = {**first_settings, **{k: args.get(k, v) for k, v in ORDER_DEFAULTS.items()}}
    window = checked_window(args["loss_window_eV"], electron_energy_eV)
    if not isinstance(material, Material):
        raise ValueError("material must be a Material")
    active = [o for o in material.oscillators if o.strength_eV2]
    plan = preflight(args["model"], args["max_order"], window, primary, edges, settings, zero=not active)
    conventions = {
        k: args[k]
        for k in (
            "weight_length_kind",
            "extra_inverse_length_per_A",
            "weight_length_seed_A",
            "multiplier",
            "multiplier_A",
        )
    }
    deliveries = []

    def receive(result, panels):
        deliveries.append(panels)

    with observe_mesh() as trace:
        result = single._run(
            material,
            electron_energy_eV,
            primary,
            edges,
            window,
            args["weight_length_A"],
            args["epsilon_momentum"],
            args["theta_deg"],
            args["hole_charge"],
            first_settings,
            on_final_panels=receive,
            **conventions,
        )
    dispersive = args["epsilon_momentum"] == "k_parallel" and any(o.dispersion for o in active)
    return _loss_orders.finish(
        result,
        deliveries[0],
        plan,
        args["max_order"],
        settings,
        {},
        deliveries=len(deliveries),
        zero=not active,
        dispersive=dispersive,
        trace=trace,
    )[0]


def multiple_loss_spectrum(
    material,
    electron_energy_eV,
    primary,
    output_edges_eV,
    *,
    loss_window_eV,
    weight_length_A,
    epsilon_momentum,
    model,
    max_order,
    neglect_tolerance_fraction=None,
    order_limit=64,
    max_convolution_points=200_000_000,
    order_epsrel=1e-6,
    order_epsabs_fraction=1e-9,
    theta_deg=0.0,
    hole_charge=1.0,
    epsrel=1e-6,
    area_epsabs=1e-10,
    spectrum_epsabs_fraction=1e-9,
    panel_limit=400,
    max_cells=2_000_000,
    kernel_epsabs=FIXED_ABS,
    kernel_epsrel=FIXED_REL,
    kernel_limit=200,
    weight_length_kind="external_constant",
    extra_inverse_length_per_A=0.0,
    weight_length_seed_A=None,
    multiplier="weight_length",
    multiplier_A=None,
) -> MultipleLossSpectrumResult:
    """Fold each loss order directly into detector bins and both outside regions (KE in eV).

    max_order is a positive integer or 'display_support'. A positive optional
    tolerance limits the estimated missing cumulative display intensity as a
    fraction of the primary intensity. Order areas are dimensionless weights.
    """
    if not isinstance(primary, (BinnedPrimary, DeltaPrimary)):
        raise SingleLossInputError("invalid_primary")
    return _run_orders(material, electron_energy_eV, primary, output_edges_eV, locals(), single.SPECTRUM_DEFAULTS)


def multiple_loss_order_areas(
    material,
    electron_energy_eV,
    *,
    loss_window_eV,
    weight_length_A,
    epsilon_momentum,
    model,
    max_order,
    order_limit=64,
    max_convolution_points=200_000_000,
    order_epsrel=1e-6,
    order_epsabs_fraction=1e-9,
    theta_deg=0.0,
    hole_charge=1.0,
    epsrel=1e-6,
    area_epsabs=1e-10,
    panel_limit=400,
    kernel_epsabs=FIXED_ABS,
    kernel_epsrel=FIXED_REL,
    kernel_limit=200,
    weight_length_kind="external_constant",
    extra_inverse_length_per_A=0.0,
    weight_length_seed_A=None,
    multiplier="weight_length",
    multiplier_A=None,
) -> MultipleLossOrderAreasResult:
    """Return eleven dimensionless areas and cutoff sensitivities per requested integer order."""
    return _run_orders(material, electron_energy_eV, None, None, locals(), single.AREA_DEFAULTS)
