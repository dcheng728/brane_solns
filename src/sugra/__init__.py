"""
sugra — Symbolic supergravity solution verifier.

Define a metric, form fields, and scalars, then verify whether
the equations of motion are satisfied (Einstein, dilaton, Bianchi,
self-duality) using both symbolic simplification and numerical
spot-checking.
"""

from .geometry import Metric, HarmonicFunction, warped_product
from .forms import (
    FormField,
    exterior_derivative,
    hodge_star,
    form_norm_squared,
    form_contraction,
    form_stress_energy,
)
from .verifier import Verifier
__all__ = [
    # Geometry
    "Metric",
    "HarmonicFunction",
    "warped_product",
    # Forms
    "FormField",
    "exterior_derivative",
    "hodge_star",
    "form_norm_squared",
    "form_contraction",
    "form_stress_energy",
    # Verifier
    "Verifier",
]
