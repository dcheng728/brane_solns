"""Tensor expression tree: TensorAtom, TensorProduct, TensorSum."""

import sympy as sp
from .index import Index
from .utils import fresh_dummy_name


class TensorExpr:
    """Abstract base class for tensor expressions."""

    def __add__(self, other):
        if isinstance(other, (int, float, sp.Basic)):
            other = ScalarExpr(sp.sympify(other))
        if not isinstance(other, TensorExpr):
            return NotImplemented
        left = self._to_sum()
        right = other._to_sum()
        return TensorSum(left.terms + right.terms)

    def __radd__(self, other):
        if isinstance(other, (int, float, sp.Basic)):
            return ScalarExpr(sp.sympify(other)) + self
        return NotImplemented

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return (-self) + other

    def __neg__(self):
        return self.__mul__(sp.S.NegativeOne)

    def __mul__(self, other):
        if isinstance(other, (int, float, sp.Basic)):
            other_s = sp.sympify(other)
            if isinstance(self, TensorProduct):
                return TensorProduct(self.coeff * other_s, self.atoms)
            elif isinstance(self, TensorSum):
                return TensorSum(tuple(t * other_s for t in self.terms))
            elif isinstance(self, ScalarExpr):
                return ScalarExpr(self.expr * other_s)
        if isinstance(other, ScalarExpr):
            return self * other.expr
        if isinstance(other, TensorProduct):
            if isinstance(self, TensorProduct):
                return _mul_products(self, other)
            elif isinstance(self, TensorSum):
                return TensorSum(tuple(t * other for t in self.terms))
            elif isinstance(self, ScalarExpr):
                return TensorProduct(self.expr * other.coeff, other.atoms)
        if isinstance(other, TensorSum):
            if isinstance(self, (TensorProduct, TensorSum)):
                return TensorSum(tuple(
                    self_t * other_t
                    for self_t in self._to_sum().terms
                    for other_t in other.terms
                ))
            elif isinstance(self, ScalarExpr):
                return TensorSum(tuple(
                    TensorProduct(self.expr * t.coeff, t.atoms)
                    for t in other.terms
                ))
        return NotImplemented

    def __rmul__(self, other):
        if isinstance(other, (int, float, sp.Basic)):
            return self.__mul__(other)
        if isinstance(other, TensorExpr):
            return other.__mul__(self)
        return NotImplemented

    def __truediv__(self, other):
        if isinstance(other, (int, float, sp.Basic)):
            return self * (sp.S.One / sp.sympify(other))
        return NotImplemented

    def _to_sum(self):
        """Convert to TensorSum for uniform handling."""
        if isinstance(self, TensorSum):
            return self
        if isinstance(self, TensorProduct):
            return TensorSum((self,))
        if isinstance(self, ScalarExpr):
            return TensorSum((TensorProduct(self.expr, ()),))
        raise TypeError(f"Cannot convert {type(self)} to TensorSum")

    @property
    def free_indices(self):
        """Free indices (those not contracted). Returns a tuple of Index."""
        raise NotImplementedError

    @property
    def rank(self):
        return len(self.free_indices)

    def canon_bp(self):
        """Canonicalize using Butler-Portugal algorithm."""
        from .canonicalize import canonicalize_expr
        return canonicalize_expr(self)

    def simplify(self, coeff_func=sp.cancel):
        """Full simplification: canonicalize + collect + simplify coefficients."""
        expr = self.canon_bp()
        if isinstance(expr, TensorSum):
            expr = expr.collect()
            if coeff_func:
                expr = expr.simplify_coefficients(coeff_func)
        return expr

    def equals(self, other):
        """Check equality by canonicalizing the difference."""
        diff = (self - other).simplify()
        if isinstance(diff, TensorSum) and len(diff.terms) == 0:
            return True
        if isinstance(diff, TensorProduct) and diff.coeff == 0:
            return True
        if isinstance(diff, ScalarExpr) and diff.expr == 0:
            return True
        return False

    def to_latex(self):
        from .latex import to_latex
        return to_latex(self)


class TensorAtom:
    """A single tensor head applied to specific indices (not an expression by itself).

    Attributes
    ----------
    head : TensorHead
    indices : tuple of Index
    """

    def __init__(self, head, indices):
        self.head = head
        self.indices = tuple(indices)

    def __repr__(self):
        # Derivative heads: show as op_name(deriv_indices, operand(operand_indices))
        info = getattr(self.head, '_deriv_info', None)
        if info is not None:
            op_name, operand_name, n_deriv = info
            deriv_idx = ', '.join(repr(i) for i in self.indices[:n_deriv])
            if operand_name:
                operand_idx = ', '.join(repr(i) for i in self.indices[n_deriv:])
                return f"{op_name}_{{{deriv_idx}}}({operand_name}({operand_idx}))"
            else:
                return f"{op_name}({deriv_idx})"
        idx_str = ', '.join(repr(i) for i in self.indices)
        return f"{self.head.name}({idx_str})"

    def __eq__(self, other):
        if not isinstance(other, TensorAtom):
            return NotImplemented
        return self.head is other.head and self.indices == other.indices

    def __hash__(self):
        return hash((id(self.head), self.indices))


class TensorProduct(TensorExpr):
    """A monomial: coeff * atom_1 * atom_2 * ...

    Attributes
    ----------
    coeff : sp.Expr
        Scalar coefficient.
    atoms : tuple of TensorAtom
        The tensor factors.
    """

    def __init__(self, coeff, atoms):
        self.coeff = sp.sympify(coeff)
        self.atoms = tuple(atoms)

    @classmethod
    def from_atom(cls, atom):
        """Create a TensorProduct from a single TensorAtom with coeff=1."""
        return cls(sp.S.One, (atom,))

    @property
    def all_indices(self):
        """All indices across all atoms, in order."""
        result = []
        for atom in self.atoms:
            result.extend(atom.indices)
        return tuple(result)

    @property
    def dummy_pairs(self):
        """Find contracted index pairs.

        Returns list of ((atom_i, slot_i), (atom_j, slot_j)) pairs
        where the two indices have the same name+type but opposite variance.
        """
        all_idx = []
        for ai, atom in enumerate(self.atoms):
            for si, idx in enumerate(atom.indices):
                all_idx.append((ai, si, idx))

        pairs = []
        used = set()
        for i, (ai, si, idx_i) in enumerate(all_idx):
            if i in used:
                continue
            for j, (aj, sj, idx_j) in enumerate(all_idx):
                if j <= i or j in used:
                    continue
                if idx_i.matches(idx_j):
                    pairs.append(((ai, si), (aj, sj)))
                    used.add(i)
                    used.add(j)
                    break
        return pairs

    @property
    def free_indices(self):
        """Indices that are not contracted (appear only once)."""
        all_idx = self.all_indices
        dummy_names = set()
        for (ai, si), (aj, sj) in self.dummy_pairs:
            dummy_names.add(self.atoms[ai].indices[si].name)

        return tuple(idx for idx in all_idx
                     if idx.name not in dummy_names)

    def __repr__(self):
        if not self.atoms:
            return repr(self.coeff)
        atoms_str = '*'.join(repr(a) for a in self.atoms)
        if self.coeff == 1:
            return atoms_str
        if self.coeff == -1:
            return f"-{atoms_str}"
        return f"{self.coeff}*{atoms_str}"

    def __eq__(self, other):
        if not isinstance(other, TensorProduct):
            return NotImplemented
        return (self.coeff == other.coeff
                and self.atoms == other.atoms)

    def __hash__(self):
        return hash((self.coeff, self.atoms))


class TensorSum(TensorExpr):
    """Sum of TensorProduct terms.

    Attributes
    ----------
    terms : tuple of TensorProduct
    """

    def __init__(self, terms):
        # Flatten nested sums and filter zeros
        flat = []
        for t in terms:
            if isinstance(t, TensorSum):
                flat.extend(t.terms)
            elif isinstance(t, TensorProduct):
                if t.coeff != 0:
                    flat.append(t)
            elif isinstance(t, TensorExpr):
                # Accept any TensorExpr (OperatorExpr, _ScaledOperator, etc.)
                if not (hasattr(t, 'coeff') and t.coeff == 0):
                    flat.append(t)
            else:
                raise TypeError(f"TensorSum expects TensorExpr terms, got {type(t)}")
        self.terms = tuple(flat)

    @property
    def free_indices(self):
        if not self.terms:
            return ()
        return self.terms[0].free_indices

    def collect(self):
        """Combine terms with identical tensor structure, summing coefficients."""
        groups = {}
        for term in self.terms:
            # Key = the atoms with coeff factored out
            key = term.atoms
            # Also need to track the index pattern
            idx_key = tuple(
                (a.head, a.indices) for a in term.atoms
            )
            full_key = idx_key
            if full_key in groups:
                groups[full_key] = (groups[full_key][0] + term.coeff, term.atoms)
            else:
                groups[full_key] = (term.coeff, term.atoms)

        new_terms = []
        for key, (coeff, atoms) in groups.items():
            coeff = sp.cancel(coeff)
            if coeff != 0:
                new_terms.append(TensorProduct(coeff, atoms))
        return TensorSum(tuple(new_terms))

    def simplify_coefficients(self, func=sp.cancel):
        """Apply a simplification function to all scalar coefficients."""
        return TensorSum(tuple(
            TensorProduct(func(t.coeff), t.atoms) for t in self.terms
        ))

    def expand(self):
        """Distribute products over sums."""
        return self  # Already in sum-of-products form

    def __repr__(self):
        if not self.terms:
            return "0"
        return ' + '.join(repr(t) for t in self.terms)

    def __len__(self):
        return len(self.terms)

    def __iter__(self):
        return iter(self.terms)

    def __getitem__(self, idx):
        return self.terms[idx]


class ScalarExpr(TensorExpr):
    """A scalar (rank-0) tensor expression wrapping a SymPy Expr."""

    def __init__(self, expr):
        self.expr = sp.sympify(expr)

    @property
    def free_indices(self):
        return ()

    def __repr__(self):
        return repr(self.expr)

    def __eq__(self, other):
        if isinstance(other, ScalarExpr):
            return self.expr == other.expr
        if isinstance(other, (int, float, sp.Basic)):
            return self.expr == sp.sympify(other)
        return NotImplemented

    def __hash__(self):
        return hash(self.expr)


# ── Helpers ──────────────────────────────────────────────────────────────

def _mul_products(left, right):
    """Multiply two TensorProducts."""
    new_coeff = left.coeff * right.coeff
    new_atoms = left.atoms + right.atoms
    return TensorProduct(new_coeff, new_atoms)
