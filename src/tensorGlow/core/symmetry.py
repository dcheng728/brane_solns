"""Tensor symmetry declarations via BSGS (Base and Strong Generating Set).

Wraps SymPy's combinatorics for use with the tensorGlow canonicalization engine.
"""

from sympy.combinatorics import Permutation
from sympy.combinatorics.tensor_can import (
    get_symmetric_group_sgs,
    bsgs_direct_product,
    riemann_bsgs,
)


class TensorSymmetry:
    """Symmetry of a tensor's index slots, stored as BSGS.

    The BSGS (base, strong generating set) is the format required by
    sympy.combinatorics.tensor_can.canonicalize.

    Parameters
    ----------
    base : list of int
        BSGS base.
    generators : list of Permutation
        BSGS strong generating set.
    rank : int
        Number of index slots.
    """

    def __init__(self, base, generators, rank):
        self.base = list(base)
        self.generators = list(generators)
        self.rank = rank

    def __repr__(self):
        return f"TensorSymmetry(rank={self.rank})"

    @classmethod
    def no_symmetry(cls, rank):
        """No index symmetry."""
        # Identity permutation on rank + 2 elements (tensor_can convention:
        # rank index slots + 2 sign indicators)
        size = rank + 2 if rank > 0 else 2
        base = []
        generators = [Permutation(size - 1)]  # identity
        return cls(base, generators, rank)

    @classmethod
    def fully_symmetric(cls, rank):
        """Fully symmetric in all slots: T_{(a1 a2 ... an)}."""
        if rank <= 1:
            return cls.no_symmetry(rank)
        base, gens = get_symmetric_group_sgs(rank, antisym=False)
        return cls(base, gens, rank)

    @classmethod
    def fully_antisymmetric(cls, rank):
        """Fully antisymmetric in all slots: T_{[a1 a2 ... an]}."""
        if rank <= 1:
            return cls.no_symmetry(rank)
        base, gens = get_symmetric_group_sgs(rank, antisym=True)
        return cls(base, gens, rank)

    @classmethod
    def riemann(cls):
        """Riemann tensor symmetry: R_{abcd}.

        - Antisymmetric in (a,b)
        - Antisymmetric in (c,d)
        - Symmetric under (a,b) <-> (c,d) pair exchange
        - First Bianchi identity (monoterm part encoded in BSGS)

        Uses SymPy's built-in riemann_bsgs which produces generators
        on 2*4+2 = 6 elements (4 index slots + 2 sign slots).
        """
        base, gens = riemann_bsgs
        return cls(list(base), list(gens), 4)

    @classmethod
    def symmetric_pair(cls, rank, pair_indices):
        """Symmetric under exchange of a specific pair of slots.

        Parameters
        ----------
        rank : int
        pair_indices : tuple of (int, int)
            The two slot positions that are symmetric.
        """
        i, j = pair_indices
        size = 2 * rank + 2
        arr = list(range(size))
        # Swap slots 2*i,2*i+1 with 2*j,2*j+1
        arr[2 * i], arr[2 * j] = 2 * j, 2 * i
        arr[2 * i + 1], arr[2 * j + 1] = 2 * j + 1, 2 * i + 1
        gen = Permutation(arr)
        return cls([min(2 * i, 2 * j)], [gen, Permutation(size - 1)], rank)

    @classmethod
    def block(cls, *symmetries):
        """Direct product of symmetries for different slot groups.

        Parameters
        ----------
        *symmetries : TensorSymmetry
            Symmetries for consecutive slot blocks.
        """
        if len(symmetries) == 0:
            return cls.no_symmetry(0)
        if len(symmetries) == 1:
            return symmetries[0]

        base = symmetries[0].base
        gens = symmetries[0].generators
        for s in symmetries[1:]:
            base, gens = bsgs_direct_product(base, gens, s.base, s.generators)

        total_rank = sum(s.rank for s in symmetries)
        return cls(base, gens, total_rank)
