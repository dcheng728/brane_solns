"""Expression tree for tensorGlow.

Five node types, all subclasses of Expr:

    Scalar(sympy_expr)            — a SymPy expression with no indices
    Tensor(head, indices)         — a named tensor with specific indices
    Apply(operator, index, expr)  — operator application: partial_a(g(-b,-c))
    Prod(coeff, factors)          — coeff * expr1 * expr2 * ...
    Sum(terms)                    — expr1 + expr2 + ...

Every node has .free_indices (possibly empty for scalars).
A tensor is an expression with free index structure; a scalar is
the special case with no free indices.
"""

import sympy as sp
from .index import Index


# ═════════════════════════════════════════════════════════════════════
# Base class
# ═════════════════════════════════════════════════════════════════════

class Expr:
    """Base class for all tensor expressions."""

    @property
    def free_indices(self):
        raise NotImplementedError

    @property
    def rank(self):
        return len(self.free_indices)

    # ── Arithmetic ──────────────────────────────────────────────────

    def __add__(self, other):
        if isinstance(other, (int, float, sp.Basic)):
            other = Scalar(sp.sympify(other))
        if not isinstance(other, Expr):
            return NotImplemented
        return Sum(_flatten_sum(self) + _flatten_sum(other))

    def __radd__(self, other):
        if isinstance(other, (int, float, sp.Basic)):
            return Scalar(sp.sympify(other)) + self
        return NotImplemented

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return (-self) + other

    def __neg__(self):
        return self * sp.S.NegativeOne

    def __mul__(self, other):
        if isinstance(other, (int, float, sp.Basic)):
            other_s = sp.sympify(other)
            # Absorb scalar into Prod or create one
            if isinstance(self, Prod):
                return Prod(self.coeff * other_s, self.factors)
            if isinstance(self, Scalar):
                return Scalar(self.expr * other_s)
            if isinstance(self, Sum):
                return Sum(tuple(t * other_s for t in self.terms))
            return Prod(other_s, (self,))
        if isinstance(other, Scalar):
            return self * other.expr
        if not isinstance(other, Expr):
            return NotImplemented
        # Check for operator types — they handle Expr * Op via __rmul__
        from .operator import BoundOperator, ComposedBound, ComposedAction
        if isinstance(other, (BoundOperator, ComposedBound, ComposedAction)):
            return other.__rmul__(self)
        # Distribute over Sum: (A) * (B + C) = A*B + A*C
        if isinstance(other, Sum):
            return Sum(tuple(self * t for t in other.terms))
        if isinstance(self, Sum):
            return Sum(tuple(t * other for t in self.terms))
        # Expr * Expr → Prod
        left_coeff, left_factors = _as_prod(self)
        right_coeff, right_factors = _as_prod(other)
        return Prod(left_coeff * right_coeff, left_factors + right_factors)

    def __rmul__(self, other):
        if isinstance(other, (int, float, sp.Basic)):
            return self.__mul__(other)
        if isinstance(other, Expr):
            return other.__mul__(self)
        return NotImplemented

    def __truediv__(self, other):
        if isinstance(other, (int, float, sp.Basic)):
            return self * (sp.S.One / sp.sympify(other))
        return NotImplemented

    # ── Simplification ──────────────────────────────────────────────

    def canon_bp(self):
        """Canonicalize using Butler-Portugal."""
        from .canonicalize import canonicalize_expr
        return canonicalize_expr(self)

    def simplify(self, coeff_func=sp.cancel):
        """Canonicalize + collect + simplify coefficients."""
        expr = self.canon_bp()
        if isinstance(expr, Sum):
            expr = expr.collect()
            if coeff_func:
                expr = expr.simplify_coefficients(coeff_func)
        return expr

    def __repr__(self):
        from ..utils.io import to_console
        return to_console(self)

    def to_latex(self):
        from ..utils.io import to_latex
        return to_latex(self)


# ═════════════════════════════════════════════════════════════════════
# Scalar
# ═════════════════════════════════════════════════════════════════════

class Scalar(Expr):
    """A SymPy scalar expression with no free indices."""

    def __init__(self, expr):
        self.expr = sp.sympify(expr)

    @property
    def free_indices(self):
        return ()

    def __eq__(self, other):
        if isinstance(other, Scalar):
            return self.expr == other.expr
        if isinstance(other, (int, float, sp.Basic)):
            return self.expr == sp.sympify(other)
        return NotImplemented

    def __hash__(self):
        return hash(self.expr)


# ═════════════════════════════════════════════════════════════════════
# Tensor
# ═════════════════════════════════════════════════════════════════════

class Tensor(Expr):
    """A named tensor with specific indices.

    Attributes
    ----------
    head : TensorHead
    indices : tuple of Index
    """

    def __init__(self, head, indices):
        self.head = head
        self.indices = tuple(indices)

    @property
    def free_indices(self):
        # In a standalone Tensor, all indices are free
        # (dummy detection happens at the Prod level)
        return self.indices

    def __eq__(self, other):
        if not isinstance(other, Tensor):
            return NotImplemented
        return self.head is other.head and self.indices == other.indices

    def __hash__(self):
        return hash((id(self.head), self.indices))


# ═════════════════════════════════════════════════════════════════════
# Apply
# ═════════════════════════════════════════════════════════════════════

class Apply(Expr):
    """Application of an operator to an expression.

    Attributes
    ----------
    op : object
        The operator (e.g. Partial, CovariantD).
    deriv_index : Index
        The operator's index (covariant).
    operand : Expr
        What the operator acts on.
    """

    def __init__(self, op, deriv_index, operand):
        self.op = op
        self.deriv_index = deriv_index if not deriv_index.is_up else -deriv_index
        self.operand = operand

    @property
    def free_indices(self):
        return (self.deriv_index,) + self.operand.free_indices

    def __eq__(self, other):
        if not isinstance(other, Apply):
            return NotImplemented
        return (self.op is other.op
                and self.deriv_index == other.deriv_index
                and self.operand == other.operand)

    def __hash__(self):
        return hash((id(self.op), self.deriv_index, self.operand))

    def expand(self):
        """Expand the operator application using the operator's defining formula.

        For D: D_a(T) → partial_a(T) + Gamma connection terms.
        For partial: stays as Apply (already primitive).
        """
        if hasattr(self.op, 'expand_apply'):
            return self.op.expand_apply(self.deriv_index, self.operand)
        return self

    def expand_leibniz(self):
        """Expand Leibniz rule if operand is a Prod."""
        if isinstance(self.operand, Sum):
            return Sum(tuple(
                Apply(self.op, self.deriv_index, t).expand_leibniz()
                for t in self.operand.terms
            ))
        if isinstance(self.operand, Prod) and len(self.operand.factors) > 1:
            return self.op._leibniz(self.deriv_index, self.operand)
        return self


# ═════════════════════════════════════════════════════════════════════
# Prod
# ═════════════════════════════════════════════════════════════════════

class Prod(Expr):
    """Product of expressions with a scalar coefficient.

    Attributes
    ----------
    coeff : sp.Expr
        Scalar coefficient.
    factors : tuple of Expr
        The factors (Tensor, Apply, Scalar, etc.).
    """

    def __init__(self, coeff, factors):
        self.coeff = sp.sympify(coeff)
        self.factors = tuple(factors)

    @property
    def free_indices(self):
        """Free indices = all indices minus contracted (dummy) pairs."""
        all_idx = self.all_indices
        dummy_names = set()
        for (_, _), (_, _) in self.dummy_pairs:
            pass  # just collecting names
        seen = {}
        for idx in all_idx:
            key = (idx.name, id(idx.index_type))
            seen.setdefault(key, []).append(idx)
        dummy_keys = {k for k, v in seen.items()
                      if len(v) == 2 and v[0].is_up != v[1].is_up}
        return tuple(idx for idx in all_idx
                     if (idx.name, id(idx.index_type)) not in dummy_keys)

    @property
    def all_indices(self):
        """All indices across all factors, in order."""
        result = []
        for f in self.factors:
            if isinstance(f, Tensor):
                result.extend(f.indices)
            elif isinstance(f, Apply):
                result.extend(f.free_indices)
            elif isinstance(f, Prod):
                result.extend(f.all_indices)
        return tuple(result)

    @property
    def dummy_pairs(self):
        """Find contracted index pairs across factors."""
        all_idx = []
        for fi, f in enumerate(self.factors):
            if isinstance(f, Tensor):
                for si, idx in enumerate(f.indices):
                    all_idx.append((fi, si, idx))
            elif isinstance(f, (Apply, Prod)):
                for si, idx in enumerate(f.free_indices):
                    all_idx.append((fi, si, idx))

        pairs = []
        used = set()
        for i, (fi, si, idx_i) in enumerate(all_idx):
            if i in used:
                continue
            for j, (fj, sj, idx_j) in enumerate(all_idx):
                if j <= i or j in used:
                    continue
                if idx_i.matches(idx_j):
                    pairs.append(((fi, si), (fj, sj)))
                    used.add(i)
                    used.add(j)
                    break
        return pairs

    @property
    def tensors(self):
        """Return only the Tensor factors (for canonicalization)."""
        return tuple(f for f in self.factors if isinstance(f, Tensor))

    @property
    def atoms(self):
        """Backward compatibility: same as tensors."""
        return self.tensors

    def __eq__(self, other):
        if not isinstance(other, Prod):
            return NotImplemented
        return self.coeff == other.coeff and self.factors == other.factors

    def __hash__(self):
        return hash((self.coeff, self.factors))


# ═════════════════════════════════════════════════════════════════════
# Sum
# ═════════════════════════════════════════════════════════════════════

class Sum(Expr):
    """Sum of expressions.

    Attributes
    ----------
    terms : tuple of Expr
    """

    def __init__(self, terms):
        flat = []
        for t in terms:
            if isinstance(t, Sum):
                flat.extend(t.terms)
            elif isinstance(t, Expr):
                # Filter zeros
                if isinstance(t, Prod) and t.coeff == 0:
                    continue
                if isinstance(t, Scalar) and t.expr == 0:
                    continue
                flat.append(t)
            else:
                raise TypeError(f"Sum expects Expr terms, got {type(t)}")
        self.terms = tuple(flat)

    @property
    def free_indices(self):
        if not self.terms:
            return ()
        return self.terms[0].free_indices

    def collect(self):
        """Combine terms with identical structure, summing coefficients."""
        groups = {}
        for term in self.terms:
            if isinstance(term, Prod):
                key = term.factors
                if key in groups:
                    groups[key] = (groups[key][0] + term.coeff, term.factors)
                else:
                    groups[key] = (term.coeff, term.factors)
            else:
                # Non-Prod terms: use repr as key
                key = repr(term)
                if key in groups:
                    # Can't easily sum non-Prod terms
                    groups[key + f'_{id(term)}'] = (1, term)
                else:
                    groups[key] = (1, term)

        new_terms = []
        for key, (coeff, factors_or_expr) in groups.items():
            coeff = sp.cancel(coeff)
            if coeff == 0:
                continue
            if isinstance(factors_or_expr, tuple):
                new_terms.append(Prod(coeff, factors_or_expr))
            else:
                new_terms.append(factors_or_expr)
        return Sum(tuple(new_terms))

    def simplify_coefficients(self, func=sp.cancel):
        """Apply simplification to all Prod coefficients."""
        return Sum(tuple(
            Prod(func(t.coeff), t.factors) if isinstance(t, Prod) else t
            for t in self.terms
        ))

    def __len__(self):
        return len(self.terms)

    def __iter__(self):
        return iter(self.terms)

    def __getitem__(self, idx):
        return self.terms[idx]


# ═════════════════════════════════════════════════════════════════════
# Helpers
# ═════════════════════════════════════════════════════════════════════

def _format_indexed(name, indices):
    """Format a tensor name with upper/lower indices: R^{ab}_{cd}."""
    if not indices:
        return name
    groups = []  # list of (is_up, [names])
    for idx in indices:
        idx_name = idx.name
        if groups and groups[-1][0] == idx.is_up:
            groups[-1][1].append(idx_name)
        else:
            groups.append((idx.is_up, [idx_name]))
    result = name
    for is_up, names in groups:
        joined = ' '.join(names)
        if is_up:
            result += f'^{{{joined}}}'
        else:
            result += f'_{{{joined}}}'
    return result


def _as_prod(expr):
    """Extract (coeff, factors) from any Expr."""
    if isinstance(expr, Prod):
        return (expr.coeff, expr.factors)
    if isinstance(expr, Scalar):
        return (expr.expr, ())
    # Tensor, Apply, etc. — wrap as a single factor
    return (sp.S.One, (expr,))


def _flatten_sum(expr):
    """Return a tuple of terms for Sum construction."""
    if isinstance(expr, Sum):
        return expr.terms
    if isinstance(expr, Prod):
        return (expr,)
    if isinstance(expr, Scalar):
        return (Prod(expr.expr, ()),) if expr.expr != 0 else ()
    # Tensor, Apply — wrap in Prod
    return (Prod(sp.S.One, (expr,)),)


# ═════════════════════════════════════════════════════════════════════
# Index substitution
# ═════════════════════════════════════════════════════════════════════

def replace_index(expr, old_idx, new_idx):
    """Replace all occurrences of old_idx with new_idx in an expression.

    Matches by name + index_type + variance (exact match).
    """
    if isinstance(expr, Scalar):
        return expr

    if isinstance(expr, Tensor):
        new_indices = tuple(
            new_idx if (i.name == old_idx.name
                        and i.index_type is old_idx.index_type
                        and i.is_up == old_idx.is_up)
            else i
            for i in expr.indices
        )
        return Tensor(expr.head, new_indices)

    if isinstance(expr, Apply):
        new_deriv = (new_idx if (expr.deriv_index.name == old_idx.name
                                  and expr.deriv_index.index_type is old_idx.index_type
                                  and expr.deriv_index.is_up == old_idx.is_up)
                     else expr.deriv_index)
        new_operand = replace_index(expr.operand, old_idx, new_idx)
        return Apply(expr.op, new_deriv, new_operand)

    if isinstance(expr, Prod):
        new_factors = tuple(replace_index(f, old_idx, new_idx) for f in expr.factors)
        return Prod(expr.coeff, new_factors)

    if isinstance(expr, Sum):
        new_terms = tuple(replace_index(t, old_idx, new_idx) for t in expr.terms)
        return Sum(new_terms)

    return expr


# ═════════════════════════════════════════════════════════════════════
# Backward compatibility aliases
# ═════════════════════════════════════════════════════════════════════

TensorExpr = Expr
TensorAtom = Tensor
TensorProduct = Prod
TensorSum = Sum
ScalarExpr = Scalar
