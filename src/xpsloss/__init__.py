"""Dielectric response, XPS cross sections and finite-window single-loss spectra."""

__version__ = "0.3.0.dev0"

from .dielectric import (
    DielectricConvergenceError, DielectricResponse, dielectric_response,
    negative_inverse_dielectric_function,
)
from .model import BulkConvergenceError, Material, NumericalConvergenceError, Oscillator, bulk_diimfp, loss_function
from .xps_fixed_depth import (
    FixedDepthConvergenceError, FixedDepthGapDomainError, FixedDepthResult, effective_cross_section_fixed_depth,
)

from .xps_depth_averaged import (
    DepthAverageConvergenceError, DepthAveragedResult, effective_cross_section_depth_averaged,
)
from .xps_single_loss import (
    SINGLE_LOSS_MODEL_ID, BinnedPrimary, DeltaPrimary, LossWindowAreaResult, SingleLossConvergenceError,
    SingleLossInputError, SingleLossSpectrumResult, loss_window_area, single_loss_spectrum,
)

__all__ = [
    "DepthAverageConvergenceError", "DepthAveragedResult", "effective_cross_section_depth_averaged",
    "BulkConvergenceError", "NumericalConvergenceError",
    "DielectricConvergenceError", "DielectricResponse", "Material", "Oscillator",
    "bulk_diimfp", "dielectric_response", "loss_function", "negative_inverse_dielectric_function",
    "FixedDepthConvergenceError", "FixedDepthGapDomainError", "FixedDepthResult", "effective_cross_section_fixed_depth",
    "SINGLE_LOSS_MODEL_ID", "BinnedPrimary", "DeltaPrimary", "LossWindowAreaResult", "SingleLossSpectrumResult",
    "SingleLossInputError", "SingleLossConvergenceError", "loss_window_area", "single_loss_spectrum",
]
