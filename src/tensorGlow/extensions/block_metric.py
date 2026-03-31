"""Scalar field for dimensional reduction ansatze.

Provides ScalarField: an abstract scalar field (e.g. dilaton) whose
partial derivatives on each child space are represented as TensorHeads.
"""

import sympy as sp
from ..core.tensor_head import TensorHead
from ..core.symmetry import TensorSymmetry


class ScalarField:
    """An abstract scalar field (e.g. dilaton).

    The field is represented as a SymPy Symbol so that expressions like
    ``exp(-2*Phi/3)`` are standard SymPy expressions supporting ``diff``.

    Parameters
    ----------
    name : str
    depends_on : set of IndexType
        Child spaces the field depends on. partial_m(Phi) is nonzero
        only if m's IndexType is in depends_on.
    """

    def __init__(self, name, depends_on=None):
        self.name = name
        self.symbol = sp.Symbol(name, real=True)
        self.depends_on = set(depends_on) if depends_on else set()
        self._partial_heads = {}

    def partial_head(self, index_type):
        """Get the TensorHead for partial_m(Phi) on a given child space.

        Returns None if the field doesn't depend on that space.
        """
        if index_type not in self.depends_on:
            return None
        if index_type not in self._partial_heads:
            self._partial_heads[index_type] = TensorHead(
                f'd{self.name}', [index_type],
                TensorSymmetry.no_symmetry(1)
            )
        return self._partial_heads[index_type]

    def __repr__(self):
        return f"ScalarField({self.name!r})"
