"""Tensor declarations (TensorHead)."""

from .index import Index
from .symmetry import TensorSymmetry


class TensorHead:
    """Declaration of a named tensor with specified index types and symmetry.

    This is the "type" — not a specific expression with indices, but the
    declaration that e.g. "R is a rank-4 tensor on Lorentz with Riemann symmetry."

    Parameters
    ----------
    name : str
    index_types : list of IndexType
        One per slot.
    symmetry : TensorSymmetry, optional
        Defaults to no symmetry.

    Usage
    -----
    >>> R = TensorHead('R', [L, L, L, L], TensorSymmetry.riemann())
    >>> expr = R(a, b, -c, -d)   # creates a TensorProduct
    """

    def __init__(self, name, index_types, symmetry=None):
        self.name = name
        self.index_types = list(index_types)
        self.symmetry = symmetry or TensorSymmetry.no_symmetry(len(index_types))
        if self.symmetry.rank != len(index_types):
            raise ValueError(
                f"Symmetry rank {self.symmetry.rank} != number of index types "
                f"{len(index_types)}"
            )

    @property
    def rank(self):
        return len(self.index_types)

    def __call__(self, *indices):
        """Create a tensor expression with specific indices.

        Parameters
        ----------
        *indices : Index
            Must match the declared index_types.

        Returns
        -------
        TensorProduct with coefficient 1 and one TensorAtom.
        """
        if len(indices) != self.rank:
            raise ValueError(
                f"{self.name} expects {self.rank} indices, got {len(indices)}"
            )
        for i, (idx, itype) in enumerate(zip(indices, self.index_types)):
            if not isinstance(idx, Index):
                raise TypeError(
                    f"Slot {i} of {self.name}: expected Index, got {type(idx).__name__}"
                )
            if idx.index_type is not itype:
                raise TypeError(
                    f"Slot {i} of {self.name}: expected IndexType "
                    f"{itype.name!r}, got {idx.index_type.name!r}"
                )

        # Import here to avoid circular dependency
        from .expr import TensorAtom, TensorProduct
        atom = TensorAtom(self, tuple(indices))
        return TensorProduct.from_atom(atom)

    def __repr__(self):
        return f"TensorHead({self.name!r}, rank={self.rank})"

    def __eq__(self, other):
        if not isinstance(other, TensorHead):
            return NotImplemented
        return self is other

    def __hash__(self):
        return id(self)
