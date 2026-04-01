"""Metric tensor, index raising/lowering, and metric contraction."""

import sympy as sp
from .index import Index, IndexType
from .tensor_head import TensorHead
from .symmetry import TensorSymmetry
from .expr import TensorAtom, TensorProduct, TensorSum, ScalarExpr


class MetricTensor:
    """A metric tensor on a given IndexType.

    Provides raising/lowering of indices and the Kronecker delta.

    Parameters
    ----------
    index_type : IndexType
    name : str
        Display name (default 'g').
    """

    def __init__(self, index_type, name='g'):
        self.index_type = index_type
        self.name = name

        self.head = TensorHead(
            name, [index_type, index_type],
            TensorSymmetry.fully_symmetric(2)
        )
        self.inv_head = TensorHead(
            name, [index_type, index_type],
            TensorSymmetry.fully_symmetric(2)
        )
        self.delta_head = TensorHead(
            'delta', [index_type, index_type],
            TensorSymmetry.no_symmetry(2)
        )

    def __call__(self, i, j):
        """Dispatch on index variance:
        - g(-a, -b)  → covariant metric $g_{\mu\nu}$
        - g(a, b)    → inverse metric $g^{\mu\nu}$
        - g(a, -b) or g(-a, b) → Kronecker delta $\delta^\mu{}_\nu$
        """
        both_down = not i.is_up and not j.is_up
        both_up = i.is_up and j.is_up
        if both_down:
            return self.head(i, j)
        if both_up:
            return self.inv_head(i, j)
        # Mixed: one up, one down → Kronecker delta
        if i.matches(j):
            dim = self.index_type.dim
            return ScalarExpr(sp.sympify(dim) if dim is not None else sp.Symbol('D'))
        return self.delta_head(i, j)

    def raise_index(self, expr, idx_to_raise, dummy_name=None):
        """Raise a covariant index using the inverse metric.

        T_{...a...} → g^{ab} T_{...b...}

        Parameters
        ----------
        expr : TensorProduct
        idx_to_raise : Index (must be covariant, i.e., is_up=False)
        dummy_name : str, optional
            Name for the new dummy index. Auto-generated if None.
        """
        if idx_to_raise.is_up:
            raise ValueError(f"Index {idx_to_raise} is already contravariant")

        if dummy_name is None:
            dummy_name = _fresh_name(expr, idx_to_raise.index_type.dummy_prefix)

        dummy_up = Index(dummy_name, idx_to_raise.index_type, is_up=True)
        dummy_down = Index(dummy_name, idx_to_raise.index_type, is_up=False)

        # Replace idx_to_raise with dummy_down in expr
        new_expr = _replace_index(expr, idx_to_raise, dummy_down)

        # Multiply by g^{raised_name, dummy_up}
        raised = Index(idx_to_raise.name, idx_to_raise.index_type, is_up=True)
        metric_term = self.inv_head(raised, dummy_up)

        return metric_term * new_expr

    def lower_index(self, expr, idx_to_lower, dummy_name=None):
        """Lower a contravariant index using the metric.

        T^{...a...} → g_{ab} T^{...b...}
        """
        if not idx_to_lower.is_up:
            raise ValueError(f"Index {idx_to_lower} is already covariant")

        if dummy_name is None:
            dummy_name = _fresh_name(expr, idx_to_lower.index_type.dummy_prefix)

        dummy_up = Index(dummy_name, idx_to_lower.index_type, is_up=True)
        dummy_down = Index(dummy_name, idx_to_lower.index_type, is_up=False)

        new_expr = _replace_index(expr, idx_to_lower, dummy_up)
        lowered = Index(idx_to_lower.name, idx_to_lower.index_type, is_up=False)
        metric_term = self.head(lowered, dummy_down)

        return metric_term * new_expr

    def contract_metrics(self, expr):
        """Simplify metric contractions in an expression.

        Applies rules:
        - g_{ab} g^{bc} → delta_a^c
        - g_{ab} T^{a...} → T_{b...} (absorb metric)
        - delta_a^a → dim
        """
        if isinstance(expr, TensorSum):
            return TensorSum(tuple(
                self.contract_metrics(t) for t in expr.terms
            ))
        if not isinstance(expr, TensorProduct):
            return expr

        return _contract_metrics_product(expr, self)


def _fresh_name(expr, prefix):
    """Generate a dummy name not already used in the expression."""
    used = set()
    if isinstance(expr, TensorProduct):
        for idx in expr.all_indices:
            used.add(idx.name)
    counter = 0
    name = f"{prefix}_{counter}"
    while name in used:
        counter += 1
        name = f"{prefix}_{counter}"
    return name


def _replace_index(product, old_idx, new_idx):
    """Replace an index in a TensorProduct."""
    new_atoms = []
    for atom in product.atoms:
        new_indices = tuple(
            new_idx if (idx.name == old_idx.name
                        and idx.index_type is old_idx.index_type
                        and idx.is_up == old_idx.is_up)
            else idx
            for idx in atom.indices
        )
        new_atoms.append(TensorAtom(atom.head, new_indices))
    return TensorProduct(product.coeff, tuple(new_atoms))


def _contract_metrics_product(product, metric):
    """Attempt to eliminate metric tensors from a product by absorbing them."""
    atoms = list(product.atoms)
    coeff = product.coeff
    changed = True

    while changed:
        changed = False
        for i, atom in enumerate(atoms):
            is_metric_like = (atom.head is metric.head
                              or atom.head is metric.inv_head
                              or atom.head is metric.delta_head)
            if not is_metric_like:
                continue

            # This atom is a metric, inverse metric, or delta
            m_idx0, m_idx1 = atom.indices
            is_inverse = (atom.head is metric.inv_head)
            is_delta = (atom.head is metric.delta_head)

            # Try to find another atom that contracts with one of the metric's indices
            for j, other in enumerate(atoms):
                if i == j:
                    continue
                for s, o_idx in enumerate(other.indices):
                    # Check if metric index 0 contracts with this
                    if m_idx0.matches(o_idx):
                        if is_delta:
                            # delta^a_b T^{b...} = T^{a...}: surviving index
                            # keeps m_idx1's name with o_idx's variance
                            new_idx = Index(m_idx1.name, m_idx1.index_type,
                                            o_idx.is_up)
                        else:
                            new_idx = _compute_replacement(m_idx1, o_idx)
                        new_other_indices = list(other.indices)
                        new_other_indices[s] = new_idx
                        atoms[j] = TensorAtom(other.head, tuple(new_other_indices))
                        atoms.pop(i)
                        changed = True
                        break
                    if m_idx1.matches(o_idx):
                        if is_delta:
                            new_idx = Index(m_idx0.name, m_idx0.index_type,
                                            o_idx.is_up)
                        else:
                            new_idx = _compute_replacement(m_idx0, o_idx)
                        new_other_indices = list(other.indices)
                        new_other_indices[s] = new_idx
                        atoms[j] = TensorAtom(other.head, tuple(new_other_indices))
                        atoms.pop(i)
                        changed = True
                        break
                if changed:
                    break
            if changed:
                break

    # Handle delta traces: delta_a^a → dim
    new_atoms = []
    for atom in atoms:
        if atom.head is metric.delta_head:
            i0, i1 = atom.indices
            if i0.matches(i1):
                dim = metric.index_type.dim
                if dim is not None:
                    coeff *= dim
                    continue
        new_atoms.append(atom)

    if not new_atoms:
        return TensorProduct(coeff, ())
    return TensorProduct(coeff, tuple(new_atoms))


def _compute_replacement(metric_idx, contracted_idx):
    """Compute the replacement index when absorbing a metric contraction.

    When g_{ab} contracts with T^{a...}, index 'a' is consumed and
    the surviving metric index 'b' replaces 'a' in T. The replacement
    has the name of the surviving metric index and the variance that
    was on the contracted slot of T (but flipped to match absorption).
    """
    # The replacement takes the name of metric_idx but the variance
    # of contracted_idx (since raising/lowering flips it)
    return Index(metric_idx.name, metric_idx.index_type,
                 not metric_idx.is_up)
