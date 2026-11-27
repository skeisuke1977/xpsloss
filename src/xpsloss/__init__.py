"""Independent dielectric response, single-event XPS losses and approximate
finite-window loss-order spectra at normal emission."""

__version__ = "0.4.0.dev0"

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
    SINGLE_LOSS_MODEL_ID, SINGLE_LOSS_CONVENTIONS_MODEL_ID, BinnedPrimary, DeltaPrimary, LossWindowAreaResult, SingleLossConvergenceError,
    SingleLossInputError, SingleLossSpectrumResult, loss_window_area, single_loss_spectrum,
)

from .xps_multiple_loss import (
    MULTIPLE_LOSS_MODEL_ID, MultipleLossConvergenceError, MultipleLossInputError,
    MultipleLossOrderAreasResult, MultipleLossSpectrumResult, multiple_loss_order_areas, multiple_loss_spectrum,
)

__all__ = [
    "DepthAverageConvergenceError", "DepthAveragedResult", "effective_cross_section_depth_averaged",
    "BulkConvergenceError", "NumericalConvergenceError",
    "DielectricConvergenceError", "DielectricResponse", "Material", "Oscillator",
    "bulk_diimfp", "dielectric_response", "loss_function", "negative_inverse_dielectric_function",
    "FixedDepthConvergenceError", "FixedDepthGapDomainError", "FixedDepthResult", "effective_cross_section_fixed_depth",
    "SINGLE_LOSS_MODEL_ID", "SINGLE_LOSS_CONVENTIONS_MODEL_ID", "BinnedPrimary", "DeltaPrimary", "LossWindowAreaResult", "SingleLossSpectrumResult",
    "SingleLossInputError", "SingleLossConvergenceError", "loss_window_area", "single_loss_spectrum",
    "MULTIPLE_LOSS_MODEL_ID", "MultipleLossConvergenceError", "MultipleLossInputError",
    "MultipleLossOrderAreasResult", "MultipleLossSpectrumResult", "multiple_loss_order_areas", "multiple_loss_spectrum",
]
