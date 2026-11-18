"""Independent bulk-loss starter; not a full QUEELS-XPS implementation."""

__version__ = "0.0.1"

from .model import Material, Oscillator, bulk_diimfp, loss_function

__all__ = ["Material", "Oscillator", "bulk_diimfp", "loss_function"]
