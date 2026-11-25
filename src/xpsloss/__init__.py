"""Independent dielectric response, bulk, fixed-depth and depth-averaged XPS losses.

Not a full QUEELS-XPS implementation: oblique emission, trajectory-dependent
weight lengths, loss-window areas, spectra and multiple scattering are not
implemented.
"""

__version__ = "0.2.0.dev0"

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

__all__ = [
    "DepthAverageConvergenceError", "DepthAveragedResult", "effective_cross_section_depth_averaged",
    "BulkConvergenceError", "NumericalConvergenceError",
    "DielectricConvergenceError", "DielectricResponse", "Material", "Oscillator",
    "bulk_diimfp", "dielectric_response", "loss_function", "negative_inverse_dielectric_function",
    "FixedDepthConvergenceError", "FixedDepthGapDomainError", "FixedDepthResult", "effective_cross_section_fixed_depth",
]
