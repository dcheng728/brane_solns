"""Covariant derivative operator and derivative algebra.

Provides the CovDerivative operator, DerivativeAtom for unevaluated
derivatives, Leibniz rule expansion, and covariant derivative commutators.
"""

import sympy as sp
from .index import Index
from .tensor_head import TensorHead
from .symmetry import TensorSymmetry
from .expr import TensorAtom, TensorProduct, TensorSum, ScalarExpr


class DerivativeAtom:
    """Represents D_a(T^{bc...}) as an unevaluated derivative.

    This is treated as a tensor with indices [a] + indices_of_operand.
    It participates in tensor products and canonicalization.

    Attributes
    ----------
    deriv_index : Index
        The derivative index (covariant).
    operand_head : TensorHead
        The tensor being differentiated.
    operand_indices : tuple of Index
        The indices on the operand.
    """

    def __init__(self, deriv_index, operand_head, operand_indices):
        self.deriv_index = deriv_index
        self.operand_head = operand_head
        self.operand_indices = tuple(operand_indices)

        # Create a synthetic TensorHead for this derivative object
        # D_a T_{bc} is a rank-(1 + rank_T) tensor
        all_itypes = [deriv_index.index_type] + list(operand_head.index_types)
        self._head = TensorHead(
            f'D({operand_head.name})',
            all_itypes,
            TensorSymmetry.no_symmetry(len(all_itypes))
        )

    @property
    def all_indices(self):
        return (self.deriv_index,) + self.operand_indices

    @property
    def head(self):
        return self._head

    def to_atom(self):
        """Convert to a TensorAtom for use in expressions."""
        return TensorAtom(self._head, self.all_indices)

    def __repr__(self):
        op_str = ', '.join(repr(i) for i in self.operand_indices)
        return f"D({self.deriv_index})({self.operand_head.name}({op_str}))"


class CovDerivative:
    """Covariant derivative operator.

    Can be associated with a Levi-Civita connection (from a metric)
    and/or a gauge connection.

    Parameters
    ----------
    index_type : IndexType
        The manifold's tangent bundle.
    metric : MetricTensor, optional
        For Levi-Civita connection (D_a g_{bc} = 0).
    name : str
        Display name.
    """

    def __init__(self, index_type, metric=None, name='D'):
        self.index_type = index_type
        self.metric = metric
        self.name = name
        # Cache derivative TensorHeads
        self._deriv_heads = {}

    def __call__(self, deriv_idx, expr):
        """Apply D_{deriv_idx} to a tensor expression.

        Parameters
        ----------
        deriv_idx : Index
            The derivative index (will be made covariant).
        expr : TensorProduct or TensorSum

        Returns
        -------
        TensorExpr
        """
        # Ensure derivative index is covariant
        if deriv_idx.is_up:
            deriv_idx = -deriv_idx

        if isinstance(expr, TensorSum):
            # Linearity
            return TensorSum(tuple(
                self(deriv_idx, term) for term in expr.terms
            ))

        if isinstance(expr, TensorProduct):
            return self._apply_to_product(deriv_idx, expr)

        if isinstance(expr, ScalarExpr):
            # D_a(scalar) = partial_a(scalar)
            head = self._get_partial_head()
            atom = TensorAtom(head, (deriv_idx,))
            return TensorProduct(expr.expr, (atom,))

        raise TypeError(f"Cannot differentiate {type(expr)}")

    def _apply_to_product(self, deriv_idx, product):
        """Apply derivative to a TensorProduct.

        For a single-atom product (coeff * T_{bc}), creates D_a(T_{bc}).
        For multi-atom products, returns an unevaluated derivative.
        Use expand_leibniz() to expand.
        """
        if len(product.atoms) == 0:
            # Derivative of a scalar coefficient
            return TensorProduct(0, ())

        if len(product.atoms) == 1:
            atom = product.atoms[0]

            # Metric compatibility: D_a(g_{bc}) = 0
            if self.metric and atom.head is self.metric.head:
                return TensorProduct(0, ())
            if self.metric and atom.head is self.metric.inv_head:
                return TensorProduct(0, ())

            # Create derivative tensor
            d_head = self._get_deriv_head(atom.head)
            d_atom = TensorAtom(d_head, (deriv_idx,) + atom.indices)
            return TensorProduct(product.coeff, (d_atom,))

        # Multi-atom: create unevaluated derivative.
        # Store as a special "D_a(product)" TensorProduct.
        # For now, just return Leibniz-expanded form.
        return self._leibniz(deriv_idx, product)

    def _leibniz(self, deriv_idx, product):
        """Expand D_a(T1 * T2 * ... * Tn) via Leibniz rule.

        D_a(coeff * T1 * T2 * ... * Tn) =
            coeff * (D_a(T1)*T2*...*Tn + T1*D_a(T2)*...*Tn + ...)
        """
        terms = []
        n = len(product.atoms)
        for i in range(n):
            # Differentiate atom i, leave others unchanged
            atom_i = product.atoms[i]

            # Metric compatibility check
            if self.metric and (atom_i.head is self.metric.head
                                or atom_i.head is self.metric.inv_head):
                continue  # D_a(g) = 0, this term vanishes

            d_head = self._get_deriv_head(atom_i.head)
            d_atom = TensorAtom(d_head, (deriv_idx,) + atom_i.indices)

            other_atoms = product.atoms[:i] + (d_atom,) + product.atoms[i+1:]
            terms.append(TensorProduct(product.coeff, other_atoms))

        if not terms:
            return TensorProduct(0, ())
        if len(terms) == 1:
            return terms[0]
        return TensorSum(tuple(terms))

    def leibniz(self, expr):
        """Expand all multi-atom derivatives via Leibniz rule."""
        # For now, this is done automatically in __call__
        return expr

    def commutator(self, idx_a, idx_b, expr):
        """Compute [D_a, D_b] applied to expr.

        For a vector V^c: [D_a, D_b] V^c = R^c_{dab} V^d
        For a covector W_c: [D_a, D_b] W_c = -R^d_{cab} W_d
        For general tensor: sum over each index.

        Parameters
        ----------
        idx_a, idx_b : Index
        expr : TensorProduct (single-atom)

        Returns
        -------
        TensorExpr with Riemann tensor terms.
        """
        if not isinstance(expr, TensorProduct) or len(expr.atoms) != 1:
            raise ValueError("commutator currently supports single-atom expressions")

        atom = expr.atoms[0]
        head = atom.head
        tensor_indices = atom.indices

        # Ensure derivative indices are covariant
        if idx_a.is_up:
            idx_a = -idx_a
        if idx_b.is_up:
            idx_b = -idx_b

        if self._riemann is None:
            raise ValueError("Need a Riemann tensor to compute commutator. "
                             "Set via set_riemann().")

        terms = []
        used_names = set(i.name for i in tensor_indices)
        used_names.add(idx_a.name)
        used_names.add(idx_b.name)

        for k, idx_k in enumerate(tensor_indices):
            # Generate a fresh dummy index
            dummy_name = _fresh_dummy(idx_k.index_type, used_names)
            used_names.add(dummy_name)

            if idx_k.is_up:
                # Upper index: +R^{idx_k}_{dummy, a, b} * T^{...dummy...}
                dummy_down = Index(dummy_name, idx_k.index_type, is_up=False)
                dummy_up = Index(dummy_name, idx_k.index_type, is_up=True)

                r_atom = TensorAtom(self._riemann,
                                    (idx_k, dummy_down, idx_a, idx_b))
                new_t_indices = list(tensor_indices)
                new_t_indices[k] = dummy_up
                t_atom = TensorAtom(head, tuple(new_t_indices))
                terms.append(TensorProduct(expr.coeff, (r_atom, t_atom)))
            else:
                # Lower index: -R^{dummy}_{idx_k, a, b} * T_{...dummy...}
                dummy_up = Index(dummy_name, idx_k.index_type, is_up=True)
                dummy_down = Index(dummy_name, idx_k.index_type, is_up=False)

                r_atom = TensorAtom(self._riemann,
                                    (dummy_up, idx_k, idx_a, idx_b))
                new_t_indices = list(tensor_indices)
                new_t_indices[k] = dummy_down
                t_atom = TensorAtom(head, tuple(new_t_indices))
                terms.append(TensorProduct(-expr.coeff, (r_atom, t_atom)))

        if not terms:
            return TensorProduct(0, ())
        if len(terms) == 1:
            return terms[0]
        return TensorSum(tuple(terms))

    def set_riemann(self, riemann_head):
        """Set the Riemann tensor head for commutator computations.

        Parameters
        ----------
        riemann_head : TensorHead
            Must be rank 4 with Riemann symmetry. Convention: R^a_{bcd}.
        """
        self._riemann = riemann_head

    _riemann = None

    def _get_deriv_head(self, tensor_head):
        """Get or create the TensorHead for D(T)."""
        if tensor_head in self._deriv_heads:
            return self._deriv_heads[tensor_head]

        itypes = [self.index_type] + list(tensor_head.index_types)
        name = f'D{tensor_head.name}'
        # D_a T_{bc...} has no special symmetry in general
        d_head = TensorHead(name, itypes, TensorSymmetry.no_symmetry(len(itypes)))
        self._deriv_heads[tensor_head] = d_head
        return d_head

    def _get_partial_head(self):
        """Get the partial derivative head for scalars."""
        return TensorHead('partial', [self.index_type],
                          TensorSymmetry.no_symmetry(1))


class PartialDerivative:
    """Bare partial derivative operator — no connection, no metric compatibility.

    Unlike CovDerivative, partial_m(g_{np}) is NOT automatically zero.
    This is needed for computing Christoffel symbols from the defining formula.

    Parameters
    ----------
    index_type : IndexType
    name : str
    """

    def __init__(self, index_type, name='partial'):
        self.index_type = index_type
        self.name = name
        self._deriv_heads = {}

    def __call__(self, deriv_idx, expr):
        """Apply partial_{deriv_idx} to expr."""
        if deriv_idx.is_up:
            deriv_idx = -deriv_idx

        if isinstance(expr, TensorSum):
            return TensorSum(tuple(self(deriv_idx, t) for t in expr.terms))

        if isinstance(expr, TensorProduct):
            return self._leibniz(deriv_idx, expr)

        if isinstance(expr, ScalarExpr):
            head = self._get_deriv_head_scalar()
            atom = TensorAtom(head, (deriv_idx,))
            return TensorProduct(expr.expr, (atom,))

        raise TypeError(f"Cannot differentiate {type(expr)}")

    def _leibniz(self, deriv_idx, product):
        """Expand via Leibniz rule. No metric compatibility — all terms kept."""
        terms = []
        n = len(product.atoms)
        for i in range(n):
            atom_i = product.atoms[i]
            d_head = self._get_deriv_head(atom_i.head)
            d_atom = TensorAtom(d_head, (deriv_idx,) + atom_i.indices)
            other_atoms = product.atoms[:i] + (d_atom,) + product.atoms[i+1:]
            terms.append(TensorProduct(product.coeff, other_atoms))

        if not terms:
            return TensorProduct(0, ())
        if len(terms) == 1:
            return terms[0]
        return TensorSum(tuple(terms))

    def _get_deriv_head(self, tensor_head):
        if tensor_head in self._deriv_heads:
            return self._deriv_heads[tensor_head]
        itypes = [self.index_type] + list(tensor_head.index_types)
        d_head = TensorHead(
            f'd{tensor_head.name}', itypes,
            TensorSymmetry.no_symmetry(len(itypes))
        )
        self._deriv_heads[tensor_head] = d_head
        return d_head

    def _get_deriv_head_scalar(self):
        return TensorHead('d', [self.index_type],
                          TensorSymmetry.no_symmetry(1))


def _fresh_dummy(index_type, used):
    """Generate a fresh dummy name not in used."""
    prefix = index_type.dummy_prefix
    counter = 0
    while True:
        name = f"{prefix}_{counter}"
        if name not in used:
            return name
        counter += 1
