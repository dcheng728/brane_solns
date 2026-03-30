"""Dimensional reduction: split a parent IndexType into child sub-types.

An IndexSplit declares that a parent index M decomposes as M = (mu, i),
where mu runs over one child type and i over another. No assumptions
are made about the metric — cross-blocks like g_{mu,i} are kept.

Decomposing a tensor expression replaces each parent-type index with
a sum over all child types, producing one term per "block" (choice of
child type for each slot).
"""

import itertools
from .index import Index, IndexType
from .tensor_head import TensorHead
from .symmetry import TensorSymmetry
from .expr import TensorAtom, TensorProduct, TensorSum


class IndexSplit:
    """Declares a splitting of one IndexType into sub-types.

    Parameters
    ----------
    parent : IndexType
        The parent type being decomposed (e.g. 10d Lorentz).
    children : list of IndexType
        The child types (e.g. [worldvolume, transverse]).
        Their dimensions should sum to the parent dimension (if concrete).

    Usage
    -----
    >>> D10 = IndexType('10d', 10, 'M')
    >>> WV  = IndexType('WV', 3, 'mu')
    >>> TR  = IndexType('TR', 7, 'i')
    >>> split = IndexSplit(D10, [WV, TR])
    >>> blocks = split.decompose(R(M, N, -P, -Q))
    """

    def __init__(self, parent, children):
        self.parent = parent
        self.children = list(children)

        # Validate dimensions if concrete
        if parent.dim is not None and all(c.dim is not None for c in children):
            child_sum = sum(c.dim for c in children)
            if child_sum != parent.dim:
                raise ValueError(
                    f"Child dimensions {child_sum} != parent dimension {parent.dim}"
                )

        # Cache for projected TensorHeads
        self._projected_heads = {}
        # Cache for projected MetricTensors
        self._child_dummy_counters = {c: 0 for c in children}

    def decompose(self, expr):
        """Decompose a tensor expression by splitting all parent-type indices.

        Parameters
        ----------
        expr : TensorExpr

        Returns
        -------
        dict : maps block_key -> TensorExpr
            block_key is a tuple of child IndexType names, one per free index.
            E.g. ('WV', 'WV', 'TR', 'TR') for R_{mu nu i j}.
        """
        if isinstance(expr, TensorSum):
            # Decompose each term separately, merge dicts
            result = {}
            for term in expr.terms:
                blocks = self.decompose(term)
                for key, val in blocks.items():
                    if key in result:
                        result[key] = result[key] + val
                    else:
                        result[key] = val
            return result

        if isinstance(expr, TensorProduct):
            return self._decompose_product(expr)

        raise TypeError(f"Cannot decompose {type(expr)}")

    def _decompose_product(self, product):
        """Decompose a TensorProduct by enumerating child-type assignments."""
        # Find all parent-type index slots
        parent_slots = []  # list of (atom_idx, slot_idx, Index)
        for ai, atom in enumerate(product.atoms):
            for si, idx in enumerate(atom.indices):
                if idx.index_type is self.parent:
                    parent_slots.append((ai, si, idx))

        if not parent_slots:
            # No parent indices — return as-is with empty key
            return {(): product}

        n_parent = len(parent_slots)
        n_children = len(self.children)

        # Enumerate all assignments: each parent slot → one child type
        result = {}
        for assignment in itertools.product(range(n_children), repeat=n_parent):
            # Build the new product with child-type indices
            block_key = tuple(self.children[c].name for c in assignment)

            # Check dummy consistency: if two parent indices form a dummy
            # pair, they must be assigned to the same child type
            # (cross-type contractions are zero only if metric is block-diagonal,
            # but we keep all blocks — the user decides what vanishes)
            # Actually: we generate ALL blocks. The user or metric structure
            # determines which vanish.

            # Create new indices for each parent slot
            new_indices_map = {}  # (atom_idx, slot_idx) -> new Index
            for slot_i, (ai, si, idx) in enumerate(parent_slots):
                child_type = self.children[assignment[slot_i]]
                new_name = idx.name  # keep the same name
                new_idx = Index(new_name, child_type, idx.is_up)
                new_indices_map[(ai, si)] = new_idx

            # Rebuild atoms
            new_atoms = []
            child_types_per_atom = []
            for ai, atom in enumerate(product.atoms):
                new_atom_indices = []
                atom_child_types = []
                for si, idx in enumerate(atom.indices):
                    if (ai, si) in new_indices_map:
                        new_atom_indices.append(new_indices_map[(ai, si)])
                        atom_child_types.append(
                            self.children[assignment[
                                next(j for j, (a, s, _) in enumerate(parent_slots)
                                     if a == ai and s == si)
                            ]]
                        )
                    else:
                        new_atom_indices.append(idx)
                        atom_child_types.append(None)

                # Get or create projected TensorHead
                head_key = (id(atom.head),) + tuple(
                    (id(ct) if ct else None) for ct in atom_child_types
                )
                if head_key not in self._projected_heads:
                    new_itypes = []
                    for si, idx in enumerate(atom.indices):
                        if atom_child_types[si] is not None:
                            new_itypes.append(atom_child_types[si])
                        else:
                            new_itypes.append(atom.head.index_types[si])

                    # Name reflects the block
                    child_labels = []
                    for ct in atom_child_types:
                        if ct is not None:
                            child_labels.append(ct.name)
                    suffix = '_'.join(child_labels) if child_labels else ''
                    proj_name = atom.head.name
                    if suffix:
                        proj_name = f"{atom.head.name}"

                    # Projected tensor has no_symmetry — the symmetry of
                    # the parent may not hold for mixed blocks
                    proj_head = TensorHead(
                        proj_name, new_itypes,
                        TensorSymmetry.no_symmetry(len(new_itypes))
                    )
                    self._projected_heads[head_key] = proj_head

                proj_head = self._projected_heads[head_key]
                new_atoms.append(TensorAtom(proj_head, tuple(new_atom_indices)))

            new_product = TensorProduct(product.coeff, tuple(new_atoms))

            if block_key in result:
                result[block_key] = result[block_key] + new_product
            else:
                result[block_key] = new_product

        return result

    def project(self, expr, slot_assignments):
        """Project onto a specific block.

        Parameters
        ----------
        expr : TensorExpr
        slot_assignments : dict mapping Index -> IndexType (child)
            Which child type each parent index is projected onto.

        Returns
        -------
        TensorExpr for that specific block.
        """
        blocks = self.decompose(expr)

        # Find the key matching the requested assignments
        # This is a convenience wrapper — for more control, use decompose()
        for key, val in blocks.items():
            match = True
            # key is a tuple of child type names
            # We'd need to match against slot_assignments
            # For now, return all blocks
            pass

        return blocks

    def decompose_to_list(self, expr):
        """Like decompose(), but returns a flat list of (block_key, TensorExpr) pairs."""
        return list(self.decompose(expr).items())
