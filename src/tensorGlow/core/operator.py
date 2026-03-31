"""Operator algebra for tensorGlow.

Operators act on expressions via multiplication: op * expr = Apply(op, expr).
An operator with a fixed index is a BoundOperator.

    partial(-a) * g(-b, -c)  →  Apply(partial, -a, g(-b,-c))
    D(-a) * V(b)             →  Apply(D, -a, V(b))

Operators compose:
    partial(-a) * partial(-b) * g(-c,-d)
    = partial_a(partial_b(g_{cd}))

Operator definitions:
    partial  — primitive, commuting
    D        — defined as partial + connection (computed, not hardcoded)
"""

import sympy as sp
from .index import Index
from .expr import Expr, Tensor, Apply, Prod, Sum, Scalar


# ═════════════════════════════════════════════════════════════════════
# Operator base
# ═════════════════════════════════════════════════════════════════════

class Operator:
    """A derivative operator on a given index type.

    Calling with an index returns a BoundOperator:
        partial(-a)  →  BoundOperator(partial, -a)

    The BoundOperator applies to expressions via *:
        partial(-a) * expr  →  Apply(partial, -a, expr)
    """

    def __init__(self, name, index_type):
        self.name = name
        self.index_type = index_type

    def __call__(self, index):
        """Bind the operator to a specific index."""
        if not isinstance(index, Index):
            raise TypeError(f"Operator expects an Index, got {type(index)}")
        return BoundOperator(self, index)

    def _on_apply(self, deriv_index, operand):
        """Hook for subclasses to intercept application.

        Return None to use the default Apply node.
        Return an Expr to override (e.g. D(g) = 0).
        """
        return None

    def _leibniz(self, deriv_index, prod):
        """Expand Leibniz rule on a Prod.

        D_a(coeff * f1 * f2 * ...) = partial_a(coeff * f1 * f2 * ...)
                                   + connection terms for free indices of the whole product.

        The partial_a part uses Leibniz: partial_a(f1*f2) = partial_a(f1)*f2 + f1*partial_a(f2).
        Connection terms are added by _on_apply for the free indices of the whole expression.

        For a plain Partial operator, there are no connection terms, so this just does Leibniz.
        """
        terms = []
        coeff = prod.coeff
        factors = prod.factors

        # Leibniz with partial only (not the full operator)
        for i in range(len(factors)):
            fi = factors[i]
            # Apply the PARTIAL part only to each factor
            partial_applied = Apply(Partial(self.index_type) if not isinstance(self, Partial) else self,
                                    deriv_index, Prod(1, (fi,)))
            other = factors[:i] + factors[i+1:]
            if other:
                other_prod = Prod(coeff, other)
                terms.append(partial_applied * other_prod)
            else:
                terms.append(coeff * partial_applied if coeff != 1 else partial_applied)

        # Connection terms for the FREE indices of the whole product
        if not isinstance(self, Partial) and hasattr(self, 'christoffel') and self.christoffel is not None:
            from .expr import replace_index
            free = prod.free_indices
            used = set(i.name for i in free)
            used.add(deriv_index.name)
            if hasattr(prod, 'all_indices'):
                used.update(i.name for i in prod.all_indices)

            Gamma = self.christoffel
            for idx_k in free:
                sigma_name = _fresh_dummy(idx_k.index_type, used)
                used.add(sigma_name)
                sigma_up = Index(sigma_name, idx_k.index_type, is_up=True)
                sigma_down = Index(sigma_name, idx_k.index_type, is_up=False)

                replaced = replace_index(prod, idx_k,
                                         sigma_up if idx_k.is_up else sigma_down)

                if idx_k.is_up:
                    gamma = Tensor(Gamma, (idx_k, deriv_index, sigma_down))
                    terms.append(gamma * replaced)
                else:
                    gamma = Tensor(Gamma, (sigma_up, deriv_index, idx_k))
                    terms.append(sp.S.NegativeOne * gamma * replaced)

        if not terms:
            return Prod(0, ())
        if len(terms) == 1:
            return terms[0]
        return Sum(tuple(terms))

    def __repr__(self):
        return self.name


# ═════════════════════════════════════════════════════════════════════
# BoundOperator
# ═════════════════════════════════════════════════════════════════════

class BoundOperator(Expr):
    """An operator with its index fixed, ready to apply via *.

    BoundOperator * Expr → Apply(op, index, Expr)
    BoundOperator * BoundOperator * Expr → nested Apply
    """

    def __init__(self, op, index):
        self.op = op
        self.index = index if not index.is_up else -index

    @property
    def free_indices(self):
        # A bound operator by itself has one free index (the derivative index)
        return (self.index,)

    def __mul__(self, other):
        """Apply the operator to an expression."""
        import sympy as sp
        if isinstance(other, (int, float, sp.Basic)):
            # op * scalar: treat as op applied to scalar
            return _apply_op(self.op, self.index, Scalar(sp.sympify(other)))
        if isinstance(other, BoundOperator):
            # op1 * op2 → ComposedBound(op1, op2)
            return ComposedBound((self, other))
        if isinstance(other, ComposedBound):
            return ComposedBound((self,) + other.bound_ops)
        if isinstance(other, Expr):
            return _apply_op(self.op, self.index, other)
        return NotImplemented

    def __rmul__(self, other):
        """scalar * BoundOperator or Expr * BoundOperator.

        Expr * BoundOperator = composition: "first apply op, then multiply by Expr".
        So (Expr * BoundOp) * target = Expr * (BoundOp * target).
        """
        import sympy as sp
        if isinstance(other, (int, float, sp.Basic)):
            return ScaledBound(sp.sympify(other), self)
        if isinstance(other, Expr):
            return ComposedAction(other, self)
        return NotImplemented

    def __repr__(self):
        return f"{self.op.name}({self.index})"


class ScaledBound:
    """coeff * BoundOperator — applies as coeff * op(expr)."""

    def __init__(self, coeff, bound_op):
        import sympy as sp
        self.coeff = sp.sympify(coeff)
        self.bound_op = bound_op

    def __mul__(self, other):
        if isinstance(other, Expr):
            result = _apply_op(self.bound_op.op, self.bound_op.index, other)
            return self.coeff * result
        return NotImplemented

    def __repr__(self):
        return f"{self.coeff}*{self.bound_op}"


class ComposedAction:
    """An Expr composed with a BoundOperator.

    (Expr * BoundOp) * target = Expr * (BoundOp * target)

    This means "first apply the operator, then multiply by the expression."
    """

    def __init__(self, expr, bound_op):
        self.expr = expr
        self.bound_op = bound_op

    def __mul__(self, other):
        """Apply: first the bound operator, then multiply by self.expr."""
        if isinstance(other, BoundOperator):
            # Compose further: (E * op1) * op2 = E * (op1 composed with op2)
            return ComposedAction(self.expr, ComposedBound((self.bound_op, other))
                                  if isinstance(self.bound_op, BoundOperator)
                                  else ComposedBound(self.bound_op.bound_ops + (other,)))
        if isinstance(other, ComposedAction):
            # (E1 * op1) * (E2 * op2) — chain
            return ComposedAction(self.expr * other.expr,
                                  ComposedBound((self.bound_op, other.bound_op))
                                  if isinstance(self.bound_op, BoundOperator)
                                  else self.bound_op)
        if isinstance(other, Expr):
            # Apply: first bound_op to other, then multiply by self.expr
            if isinstance(self.bound_op, BoundOperator):
                applied = _apply_op(self.bound_op.op, self.bound_op.index, other)
            elif isinstance(self.bound_op, ComposedBound):
                applied = self.bound_op * other
            else:
                applied = other
            return self.expr * applied
        return NotImplemented

    def __repr__(self):
        return f"{self.expr}*{self.bound_op}"


class ComposedBound:
    """A chain of bound operators waiting for an operand.

    partial(-a) * partial(-b)  →  ComposedBound([bound_a, bound_b])
    ComposedBound * expr       →  Apply(a, Apply(b, expr))
    """

    def __init__(self, bound_ops):
        self.bound_ops = tuple(bound_ops)

    def __mul__(self, other):
        if isinstance(other, BoundOperator):
            return ComposedBound(self.bound_ops + (other,))
        if isinstance(other, ComposedBound):
            return ComposedBound(self.bound_ops + other.bound_ops)
        if isinstance(other, Expr):
            # Apply from right to left: op1 * op2 * expr = op1(op2(expr))
            result = other
            for bound in reversed(self.bound_ops):
                result = _apply_op(bound.op, bound.index, result)
            return result
        return NotImplemented

    def __repr__(self):
        return ' * '.join(repr(b) for b in self.bound_ops)


# ═════════════════════════════════════════════════════════════════════
# Concrete operators
# ═════════════════════════════════════════════════════════════════════

class Partial(Operator):
    """Partial derivative. Commuting, no connection."""

    def __init__(self, index_type):
        super().__init__('partial', index_type)


class CovariantD(Operator):
    """Covariant derivative (Levi-Civita), defined from partial + Christoffel.

    D_c T^{a1...}_{b1...} = partial_c T^{a1...}_{b1...}
                           + Gamma^{ak}_{c sigma} T^{...sigma...}_{...}  (for each upper index)
                           - Gamma^{sigma}_{c bk} T^{...}_{...sigma...}  (for each lower index)

    Metric compatibility D_a(g_{bc}) = 0 is NOT hardcoded — it is a consequence
    of the Christoffel formula. To verify it, expand D(g) and substitute Gamma.
    """

    def __init__(self, index_type, metric=None, christoffel=None, partial_op=None):
        super().__init__('D', index_type)
        self.metric = metric
        self.christoffel = christoffel  # TensorHead for Gamma^a_{bc}
        self.partial_op = partial_op or Partial(index_type)
        self._riemann = None

    def set_christoffel(self, christoffel_head):
        self.christoffel = christoffel_head

    def set_riemann(self, riemann_head):
        self._riemann = riemann_head

    def _on_apply(self, deriv_index, operand):
        """Expand D_c(expr) = partial_c(expr) + Gamma connection terms.

        Works on any expression by inspecting its FREE indices only.
        For each free upper index: +Gamma^{idx}_{c sigma} * expr[idx -> sigma]
        For each free lower index: -Gamma^{sigma}_{c idx} * expr[idx -> sigma]

        Dummy (contracted) indices are NOT touched — this is the correct
        definition of the covariant derivative on a tensor of any rank.
        """
        from .expr import replace_index

        if self.christoffel is None:
            return None

        free = operand.free_indices
        if not free:
            return self.partial_op(deriv_index) * operand

        used = set(i.name for i in free)
        used.add(deriv_index.name)
        # Also collect dummy names to avoid clashes
        if hasattr(operand, 'all_indices'):
            used.update(i.name for i in operand.all_indices)

        terms = []

        # Term 1: partial_c(expr)
        terms.append(self.partial_op(deriv_index) * operand)

        # Connection terms for each FREE index only
        Gamma = self.christoffel
        for idx_k in free:
            sigma_name = _fresh_dummy(idx_k.index_type, used)
            used.add(sigma_name)
            sigma_up = Index(sigma_name, idx_k.index_type, is_up=True)
            sigma_down = Index(sigma_name, idx_k.index_type, is_up=False)

            replaced = replace_index(operand, idx_k,
                                     sigma_up if idx_k.is_up else sigma_down)

            if idx_k.is_up:
                gamma = Tensor(Gamma, (idx_k, deriv_index, sigma_down))
                terms.append(gamma * replaced)
            else:
                gamma = Tensor(Gamma, (sigma_up, deriv_index, idx_k))
                terms.append(sp.S.NegativeOne * gamma * replaced)

        if len(terms) == 1:
            return terms[0]
        return Sum(tuple(terms))


# ═════════════════════════════════════════════════════════════════════
# Application logic
# ═════════════════════════════════════════════════════════════════════

def _apply_op(op, deriv_index, expr):
    """Apply operator to expression, respecting special rules."""
    # Linearity over Sum
    if isinstance(expr, Sum):
        return Sum(tuple(_apply_op(op, deriv_index, t) for t in expr.terms))

    # Leibniz over multi-factor Prod
    if isinstance(expr, Prod) and len(expr.factors) > 1:
        return op._leibniz(deriv_index, expr)

    # Check for operator-specific expansion (e.g. D = partial + Gamma)
    result = op._on_apply(deriv_index, expr)
    if result is not None:
        return result

    # Default: create Apply node
    return Apply(op, deriv_index, expr)


# ═════════════════════════════════════════════════════════════════════
# Helpers
# ═════════════════════════════════════════════════════════════════════

def _is_zero(expr):
    if isinstance(expr, Prod) and expr.coeff == 0:
        return True
    if isinstance(expr, Scalar) and expr.expr == 0:
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
