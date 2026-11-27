"""Independent dielectric response, bulk, fixed-depth and depth-averaged XPS losses,
and the single-loss spectrum with finite loss-window areas at normal emission.

Not a full QUEELS-XPS implementation: oblique emission, trajectory-dependent
weight lengths, multiple scattering and full spectra beyond single loss are
not implemented.
"""

__version__ = "0.3.0"

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

__all__ = [
    "DepthAverageConvergenceError", "DepthAveragedResult", "effective_cross_section_depth_averaged",
    "BulkConvergenceError", "NumericalConvergenceError",
    "DielectricConvergenceError", "DielectricResponse", "Material", "Oscillator",
    "bulk_diimfp", "dielectric_response", "loss_function", "negative_inverse_dielectric_function",
    "FixedDepthConvergenceError", "FixedDepthGapDomainError", "FixedDepthResult", "effective_cross_section_fixed_depth",
    "SINGLE_LOSS_MODEL_ID", "SINGLE_LOSS_CONVENTIONS_MODEL_ID", "BinnedPrimary", "DeltaPrimary", "LossWindowAreaResult", "SingleLossSpectrumResult",
    "SingleLossInputError", "SingleLossConvergenceError", "loss_window_area", "single_loss_spectrum",
]
