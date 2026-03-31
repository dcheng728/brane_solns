"""Tensor operators: partial derivative, covariant derivative.

An operator acts on a tensor expression and produces a new tensor
expression with one additional index. Operator applications are
first-class nodes in the expression tree, preserving the structural
distinction between the derivative index and the operand indices.

    partial_a(Gamma(c, -d, -b))   is an OperatorExpr, not a flat TensorAtom
    D_a(V(b))                      likewise

Operators support:
- Nesting: partial_a(partial_b(g(-c, -d)))
- Leibniz rule: partial_a(V(b) * W(c)) = partial_a(V(b))*W(c) + V(b)*partial_a(W(c))
- Metric compatibility: D_a(g(-b, -c)) = 0
"""

from .index import Index
from .expr import TensorExpr, TensorAtom, TensorProduct, TensorSum, ScalarExpr


class OperatorExpr(TensorExpr):
    """Application of a derivative operator to a tensor expression.

    Attributes
    ----------
    operator : Operator
        The operator being applied (partial or D).
    deriv_index : Index
        The derivative index (covariant).
    operand : TensorExpr
        The expression being differentiated.
    """

    def __init__(self, operator, deriv_index, operand):
        self.operator = operator
        self.deriv_index = deriv_index if not deriv_index.is_up else -deriv_index
        self.operand = operand

    @property
    def free_indices(self):
        return (self.deriv_index,) + self.operand.free_indices

    def __repr__(self):
        return f"{self.operator.name}_{{{self.deriv_index}}}({self.operand})"

    def __eq__(self, other):
        if not isinstance(other, OperatorExpr):
            return NotImplemented
        return (self.operator is other.operator
                and self.deriv_index == other.deriv_index
                and self.operand == other.operand)

    def __hash__(self):
        return hash((id(self.operator), self.deriv_index, self.operand))

    def __neg__(self):
        return _ScaledOperator(-1, self)

    def __mul__(self, other):
        import sympy as sp
        if isinstance(other, (int, float, sp.Basic)):
            return _ScaledOperator(sp.sympify(other), self)
        if isinstance(other, TensorExpr):
            return _ProductWithOperator(self, other)
        return NotImplemented

    def __rmul__(self, other):
        import sympy as sp
        if isinstance(other, (int, float, sp.Basic)):
            return _ScaledOperator(sp.sympify(other), self)
        if isinstance(other, TensorExpr):
            return _ProductWithOperator(self, other)
        return NotImplemented

    def __add__(self, other):
        left = self._to_sum()
        if isinstance(other, OperatorExpr):
            right = other._to_sum()
        elif isinstance(other, (_ScaledOperator, _ProductWithOperator)):
            right = other._to_sum()
        elif isinstance(other, TensorExpr):
            right = other._to_sum()
        else:
            return NotImplemented
        return TensorSum(left.terms + right.terms)

    def __radd__(self, other):
        if isinstance(other, TensorExpr):
            return other.__add__(self)
        return NotImplemented

    def __sub__(self, other):
        return self + (-other)

    def expand_leibniz(self):
        """Expand the Leibniz rule if the operand is a product.

        partial_a(f * T1 * T2) = (partial_a f) * T1 * T2
                                + f * partial_a(T1) * T2
                                + f * T1 * partial_a(T2)

        Returns self unchanged if the operand is a single atom.
        """
        op = self.operand
        if isinstance(op, TensorSum):
            # Linearity: distribute over sum
            return TensorSum(tuple(
                OperatorExpr(self.operator, self.deriv_index, t).expand_leibniz()
                for t in op.terms
            ))
        if isinstance(op, TensorProduct) and len(op.atoms) > 1:
            return self.operator._leibniz(self.deriv_index, op)
        return self

    def _to_sum(self):
        # Wrap in a TensorProduct with coeff=1 so it can participate in sums
        # This is a workaround — OperatorExpr in sums needs more thought
        return TensorSum((self._as_product(),))

    def _as_product(self):
        """Wrap this OperatorExpr into the product framework for algebra."""
        # For now, OperatorExpr participates in sums by being wrapped
        return _OperatorProduct(self)


class _OperatorProduct(TensorProduct):
    """A TensorProduct wrapper around an OperatorExpr for algebra compatibility."""

    def __init__(self, op_expr):
        self._op_expr = op_expr
        super().__init__(1, ())

    @property
    def free_indices(self):
        return self._op_expr.free_indices

    def __repr__(self):
        return repr(self._op_expr)


class Operator:
    """Base class for tensor derivative operators.

    Subclasses: Partial, CovariantDerivative.
    """

    def __init__(self, name, index_type):
        self.name = name
        self.index_type = index_type

    def __call__(self, deriv_index, expr):
        """Apply the operator to an expression.

        Parameters
        ----------
        deriv_index : Index
        expr : TensorExpr

        Returns
        -------
        OperatorExpr, TensorExpr, or 0
        """
        if deriv_index.is_up:
            deriv_index = -deriv_index

        if isinstance(expr, TensorSum):
            # Linearity
            return TensorSum(tuple(self(deriv_index, t) for t in expr.terms))

        if isinstance(expr, ScalarExpr):
            return self._on_scalar(deriv_index, expr)

        if isinstance(expr, TensorProduct):
            if len(expr.atoms) == 0:
                return TensorProduct(0, ())
            if len(expr.atoms) == 1 and expr.coeff == 1:
                result = self._on_atom(deriv_index, expr.atoms[0])
                if result is not None:
                    return result
                return OperatorExpr(self, deriv_index, expr)
            # Multi-atom or nontrivial coefficient: Leibniz
            return self._leibniz(deriv_index, expr)

        if isinstance(expr, OperatorExpr):
            # Nesting: op1(op2(T))
            return OperatorExpr(self, deriv_index, expr)

        return OperatorExpr(self, deriv_index, expr)

    def _on_atom(self, deriv_index, atom):
        """Override in subclasses for special rules (e.g. metric compatibility).

        Return None to fall through to generic OperatorExpr.
        """
        return None

    def _on_scalar(self, deriv_index, scalar_expr):
        """Derivative of a pure scalar."""
        return OperatorExpr(self, deriv_index, scalar_expr)

    def _leibniz(self, deriv_index, product):
        """Expand via Leibniz rule."""
        terms = []
        coeff = product.coeff
        atoms = product.atoms

        # Derivative of coefficient (if nontrivial)
        if coeff != 1 and coeff != 0:
            # For now, treat coefficient as constant (no chain rule on SymPy exprs)
            # The BlockMetric handles scalar field chain rule separately
            pass

        # Leibniz on tensor factors
        for i in range(len(atoms)):
            atom_i = atoms[i]
            d_atom_i = self._on_atom(deriv_index, atom_i)
            if d_atom_i is not None and _is_zero(d_atom_i):
                continue  # e.g. D(g) = 0
            if d_atom_i is None:
                d_atom_i = OperatorExpr(self, deriv_index,
                                        TensorProduct(1, (atom_i,)))
            other = atoms[:i] + atoms[i+1:]
            if isinstance(d_atom_i, OperatorExpr):
                # Wrap: d_atom_i * other_atoms
                if other:
                    other_prod = TensorProduct(coeff, other)
                    terms.append(_ProductWithOperator(d_atom_i, other_prod))
                else:
                    terms.append(d_atom_i if coeff == 1
                                 else _ScaledOperator(coeff, d_atom_i))
            else:
                # d_atom_i is already a TensorProduct or similar
                new_atoms = (d_atom_i.atoms if isinstance(d_atom_i, TensorProduct)
                             else (d_atom_i,))
                terms.append(TensorProduct(coeff, new_atoms + other))

        if not terms:
            return TensorProduct(0, ())
        if len(terms) == 1:
            return terms[0]
        return TensorSum(tuple(terms))


class Partial(Operator):
    """Partial derivative operator.

    No connection, no metric compatibility.
    Partials commute: partial_a(partial_b(T)) = partial_b(partial_a(T)).

    Parameters
    ----------
    index_type : IndexType
    """

    def __init__(self, index_type):
        super().__init__('partial', index_type)


class CovariantD(Operator):
    """Covariant derivative (Levi-Civita connection).

    Satisfies metric compatibility: D_a(g_{bc}) = 0.
    Commutator produces Riemann: [D_a, D_b] V^c = R^c_{dab} V^d.

    Parameters
    ----------
    index_type : IndexType
    metric : MetricTensor
    """

    def __init__(self, index_type, metric=None):
        super().__init__('D', index_type)
        self.metric = metric
        self._riemann = None

    def _on_atom(self, deriv_index, atom):
        """Metric compatibility: D_a(g_{bc}) = 0."""
        if self.metric is not None:
            if (atom.head is self.metric.head
                    or atom.head is self.metric.inv_head):
                return TensorProduct(0, ())
        return None  # fall through to generic

    def set_riemann(self, riemann_head):
        self._riemann = riemann_head

    def commutator(self, idx_a, idx_b, expr):
        """Compute [D_a, D_b] applied to a single-atom expression.

        For V^c: [D_a, D_b] V^c = R^c_{dab} V^d
        For W_c: [D_a, D_b] W_c = -R^d_{cab} W_d
        """
        if not isinstance(expr, TensorProduct) or len(expr.atoms) != 1:
            raise ValueError("commutator currently supports single-atom expressions")

        atom = expr.atoms[0]
        tensor_indices = atom.indices

        if idx_a.is_up:
            idx_a = -idx_a
        if idx_b.is_up:
            idx_b = -idx_b

        if self._riemann is None:
            raise ValueError("Need a Riemann tensor. Call set_riemann().")

        terms = []
        used_names = set(i.name for i in tensor_indices)
        used_names.add(idx_a.name)
        used_names.add(idx_b.name)

        for k, idx_k in enumerate(tensor_indices):
            dummy_name = _fresh_dummy(idx_k.index_type, used_names)
            used_names.add(dummy_name)

            if idx_k.is_up:
                dummy_down = Index(dummy_name, idx_k.index_type, is_up=False)
                dummy_up = Index(dummy_name, idx_k.index_type, is_up=True)
                r_atom = TensorAtom(self._riemann,
                                    (idx_k, dummy_down, idx_a, idx_b))
                new_t_indices = list(tensor_indices)
                new_t_indices[k] = dummy_up
                t_atom = TensorAtom(atom.head, tuple(new_t_indices))
                terms.append(TensorProduct(expr.coeff, (r_atom, t_atom)))
            else:
                dummy_up = Index(dummy_name, idx_k.index_type, is_up=True)
                dummy_down = Index(dummy_name, idx_k.index_type, is_up=False)
                r_atom = TensorAtom(self._riemann,
                                    (dummy_up, idx_k, idx_a, idx_b))
                new_t_indices = list(tensor_indices)
                new_t_indices[k] = dummy_down
                t_atom = TensorAtom(atom.head, tuple(new_t_indices))
                terms.append(TensorProduct(-expr.coeff, (r_atom, t_atom)))

        if not terms:
            return TensorProduct(0, ())
        if len(terms) == 1:
            return terms[0]
        return TensorSum(tuple(terms))


# ── Display helpers for products involving OperatorExpr ───────────────

class _ProductWithOperator(TensorExpr):
    """Represents op_expr * other_product for display in sums."""

    def __init__(self, op_expr, other):
        self.op_expr = op_expr
        self.other = other

    @property
    def free_indices(self):
        return self.op_expr.free_indices + (self.other.free_indices if hasattr(self.other, 'free_indices') else ())

    def __repr__(self):
        return f"{self.op_expr}*{self.other}"

    def __neg__(self):
        return _ProductWithOperator(self.op_expr, -self.other)

    def __add__(self, other):
        left = self._to_sum()
        if hasattr(other, '_to_sum'):
            right = other._to_sum()
        else:
            return NotImplemented
        return TensorSum(left.terms + right.terms)

    def __radd__(self, other):
        if hasattr(other, '_to_sum'):
            return other.__add__(self)
        return NotImplemented

    def __sub__(self, other):
        return self + (-other)

    @property
    def coeff(self):
        return getattr(self.other, 'coeff', 1)

    @property
    def atoms(self):
        return getattr(self.other, 'atoms', ())

    def _to_sum(self):
        return TensorSum((self,))


class _ScaledOperator(TensorExpr):
    """Represents coeff * OperatorExpr."""

    def __init__(self, coeff, op_expr):
        import sympy as sp
        self._coeff = sp.sympify(coeff)
        self.op_expr = op_expr

    @property
    def free_indices(self):
        return self.op_expr.free_indices

    def __repr__(self):
        if self._coeff == 1:
            return repr(self.op_expr)
        if self._coeff == -1:
            return f"-{self.op_expr}"
        return f"{self._coeff}*{self.op_expr}"

    def __neg__(self):
        return _ScaledOperator(-self._coeff, self.op_expr)

    def __add__(self, other):
        left = self._to_sum()
        if hasattr(other, '_to_sum'):
            right = other._to_sum()
        else:
            return NotImplemented
        return TensorSum(left.terms + right.terms)

    def __radd__(self, other):
        if hasattr(other, '_to_sum'):
            return other.__add__(self)
        return NotImplemented

    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, other):
        import sympy as sp
        if isinstance(other, (int, float, sp.Basic)):
            return _ScaledOperator(self._coeff * sp.sympify(other), self.op_expr)
        if isinstance(other, TensorExpr):
            return _ProductWithOperator(self.op_expr, other) if self._coeff == 1 else _ProductWithOperator(self.op_expr, self._coeff * other)
        return NotImplemented

    def __rmul__(self, other):
        import sympy as sp
        if isinstance(other, (int, float, sp.Basic)):
            return _ScaledOperator(self._coeff * sp.sympify(other), self.op_expr)
        return NotImplemented

    @property
    def coeff(self):
        return self._coeff

    @property
    def atoms(self):
        return ()

    def _to_sum(self):
        return TensorSum((self,))


# ── Helpers ──────────────────────────────────────────────────────────

def _is_zero(expr):
    if isinstance(expr, TensorProduct) and expr.coeff == 0:
        return True
    if isinstance(expr, ScalarExpr) and expr.expr == 0:
        return True
    return False


def _fresh_dummy(index_type, used):
    prefix = index_type.dummy_prefix
    counter = 0
    while True:
        name = f"{prefix}_{counter}"
        if name not in used:
            return name
        counter += 1
