"""Independent dielectric response and bulk losses, not a full QUEELS-XPS implementation."""

__version__ = "0.0.1"

from .dielectric import (
    DielectricConvergenceError, DielectricResponse, dielectric_response,
    negative_inverse_dielectric_function,
)
from .model import Material, Oscillator, bulk_diimfp, loss_function

__all__ = [
    "DielectricConvergenceError", "DielectricResponse", "Material", "Oscillator",
    "bulk_diimfp", "dielectric_response", "loss_function", "negative_inverse_dielectric_function",
]
