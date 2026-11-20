"""Independent dielectric response, bulk losses and fixed-depth XPS cross sections.

Not a full QUEELS-XPS implementation: depth averaging, oblique emission, gapped
materials in the fixed-depth model and spectra are not implemented.
"""

__version__ = "0.1.0"

from .dielectric import (
    DielectricConvergenceError, DielectricResponse, dielectric_response,
    negative_inverse_dielectric_function,
)
from .model import BulkConvergenceError, Material, NumericalConvergenceError, Oscillator, bulk_diimfp, loss_function
from .xps_fixed_depth import (
    FixedDepthConvergenceError, FixedDepthResult, effective_cross_section_fixed_depth,
)

__all__ = [
    "BulkConvergenceError", "NumericalConvergenceError",
    "DielectricConvergenceError", "DielectricResponse", "Material", "Oscillator",
    "bulk_diimfp", "dielectric_response", "loss_function", "negative_inverse_dielectric_function",
    "FixedDepthConvergenceError", "FixedDepthResult", "effective_cross_section_fixed_depth",
]
