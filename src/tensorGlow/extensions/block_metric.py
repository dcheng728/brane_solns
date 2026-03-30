"""Block metric computation engine for dimensional reduction.

Given a metric expressed as blocks in a parent→children index split,
computes Christoffel symbols, Riemann tensor, Ricci tensor, and Ricci
scalar from the defining formulas of Riemannian geometry.

No assumptions about block-diagonality — all cross-blocks are kept.
"""

import sympy as sp
import itertools
from ..core.index import Index, IndexType, indices
from ..core.tensor_head import TensorHead
from ..core.symmetry import TensorSymmetry
from ..core.expr import TensorAtom, TensorProduct, TensorSum, ScalarExpr
from ..core.derivative import PartialDerivative


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


class BlockMetric:
    """A metric expressed as blocks in a dimensional reduction ansatz.

    Parameters
    ----------
    split : IndexSplit
        The parent → children decomposition.
    blocks : dict of (IndexType, IndexType) → TensorExpr or sp.Expr
        Metric blocks. Keys are pairs of child types.
        Values are tensorGlow expressions in child-space tensors,
        with SymPy scalar coefficients that may involve ScalarFields.
        For a scalar block (e.g. g_{zz}), use a plain SymPy expression.
    inv_blocks : dict, optional
        Inverse metric blocks. Same format. If not provided, the user
        must supply them (automatic inversion is not yet implemented).
    scalar_fields : list of ScalarField
        All scalar fields appearing in block coefficients.
    truncation : set of IndexType
        Child types where partial derivatives vanish (e.g. the circle
        direction in 11d → IIA where fields are z-independent).
    """

    def __init__(self, split, blocks, inv_blocks=None, scalar_fields=None,
                 truncation=None):
        self.split = split
        self.children = split.children
        self.blocks = {}
        for key, val in blocks.items():
            if isinstance(val, (int, float, sp.Basic)):
                val = ScalarExpr(sp.sympify(val))
            self.blocks[key] = val

        if inv_blocks is not None:
            self.inv_blocks = {}
            for key, val in inv_blocks.items():
                if isinstance(val, (int, float, sp.Basic)):
                    val = ScalarExpr(sp.sympify(val))
                self.inv_blocks[key] = val
        else:
            self.inv_blocks = None

        self.scalar_fields = list(scalar_fields or [])
        self.truncation = set(truncation or [])

        # Create partial derivative operators for each child type
        self._partials = {}
        for child in self.children:
            self._partials[child] = PartialDerivative(child, f'partial_{child.name}')

    def get_block(self, type_row, type_col):
        """Get metric block g_{type_row, type_col}."""
        return self.blocks.get((type_row, type_col), ScalarExpr(sp.S.Zero))

    def get_inv_block(self, type_row, type_col):
        """Get inverse metric block g^{type_row, type_col}."""
        if self.inv_blocks is None:
            raise ValueError("Inverse metric blocks not provided. "
                             "Pass inv_blocks to BlockMetric constructor.")
        return self.inv_blocks.get((type_row, type_col), ScalarExpr(sp.S.Zero))

    def _differentiate(self, deriv_idx, expr):
        """Differentiate a block expression w.r.t. deriv_idx.

        Handles:
        - Truncation (d/dz = 0)
        - Scalar field chain rule
        - Leibniz rule on tensor factors
        """
        deriv_type = deriv_idx.index_type

        # Truncation
        if deriv_type in self.truncation:
            return TensorProduct(0, ())

        if isinstance(expr, ScalarExpr):
            # Pure scalar: differentiate via chain rule on scalar fields
            return self._diff_scalar(deriv_idx, expr.expr)

        if isinstance(expr, TensorSum):
            terms = [self._differentiate(deriv_idx, t) for t in expr.terms]
            return TensorSum(tuple(terms))

        if isinstance(expr, TensorProduct):
            return self._diff_product(deriv_idx, expr)

        raise TypeError(f"Cannot differentiate {type(expr)}")

    def _diff_scalar(self, deriv_idx, sympy_expr):
        """Differentiate a SymPy scalar expression via chain rule on scalar fields."""
        terms = []
        for sf in self.scalar_fields:
            d_coeff = sp.diff(sympy_expr, sf.symbol)
            if d_coeff != 0:
                ph = sf.partial_head(deriv_idx.index_type)
                if ph is not None:
                    atom = TensorAtom(ph, (deriv_idx,))
                    terms.append(TensorProduct(d_coeff, (atom,)))
        if not terms:
            return TensorProduct(0, ())
        if len(terms) == 1:
            return terms[0]
        return TensorSum(tuple(terms))

    def _diff_product(self, deriv_idx, product):
        """Differentiate a TensorProduct via Leibniz + scalar chain rule."""
        terms = []
        coeff = product.coeff
        atoms = product.atoms

        # 1. Differentiate the scalar coefficient (chain rule)
        if coeff != 1 and coeff != 0:
            d_coeff_expr = self._diff_scalar(deriv_idx, coeff)
            if not _is_zero(d_coeff_expr):
                # d(coeff) * atoms
                if isinstance(d_coeff_expr, TensorProduct):
                    new_atoms = d_coeff_expr.atoms + atoms
                    terms.append(TensorProduct(d_coeff_expr.coeff, new_atoms))
                elif isinstance(d_coeff_expr, TensorSum):
                    for t in d_coeff_expr.terms:
                        new_atoms = t.atoms + atoms
                        terms.append(TensorProduct(t.coeff, new_atoms))

        # 2. Differentiate each tensor factor (Leibniz)
        partial = self._partials.get(deriv_idx.index_type)
        if partial is not None:
            for i, atom in enumerate(atoms):
                d_head = partial._get_deriv_head(atom.head)
                d_atom = TensorAtom(d_head, (deriv_idx,) + atom.indices)
                other = atoms[:i] + (d_atom,) + atoms[i+1:]
                terms.append(TensorProduct(coeff, other))

        if not terms:
            return TensorProduct(0, ())
        if len(terms) == 1:
            return terms[0]
        return TensorSum(tuple(terms))

    def christoffel(self):
        r"""Compute all Christoffel symbol blocks.

        .. math::
            \Gamma^A_{BC} = \frac{1}{2} g^{AD}
            (\partial_B g_{DC} + \partial_C g_{DB} - \partial_D g_{BC})

        Returns
        -------
        dict : (type_A, type_B, type_C) → TensorExpr
        """
        result = {}
        children = self.children

        for type_A, type_B, type_C in itertools.product(children, repeat=3):
            total = None

            # Fresh indices for this block
            a = _fresh_idx(type_A, 'a', is_up=True)
            b = _fresh_idx(type_B, 'b', is_up=False)
            c = _fresh_idx(type_C, 'c', is_up=False)

            for type_D in children:
                d_down = _fresh_idx(type_D, 'd', is_up=False)

                # g^{AD}
                g_inv_AD = self.get_inv_block(type_A, type_D)
                if _is_zero(g_inv_AD):
                    continue

                # partial_B(g_{DC})
                g_DC = self.get_block(type_D, type_C)
                term1 = self._differentiate(b, g_DC) if not _is_zero(g_DC) else None

                # partial_C(g_{DB})
                g_DB = self.get_block(type_D, type_B)
                term2 = self._differentiate(c, g_DB) if not _is_zero(g_DB) else None

                # partial_D(g_{BC})
                g_BC = self.get_block(type_B, type_C)
                term3 = self._differentiate(d_down, g_BC) if not _is_zero(g_BC) else None

                # Combine: 1/2 * g^{AD} * (term1 + term2 - term3)
                bracket = _add_exprs(term1, term2, _neg(term3))
                if bracket is None or _is_zero(bracket):
                    continue

                contribution = sp.Rational(1, 2) * g_inv_AD * bracket
                total = contribution if total is None else total + contribution

            result[(type_A, type_B, type_C)] = total if total is not None else TensorProduct(0, ())

        return result

    def ricci_scalar(self, christoffel_blocks=None):
        r"""Compute the Ricci scalar by contracting Ricci tensor with inverse metric.

        This is a convenience that computes Christoffel → Riemann → Ricci → scalar.
        For now, returns the Christoffel blocks as a starting point.
        Full Riemann computation will be added.
        """
        if christoffel_blocks is None:
            christoffel_blocks = self.christoffel()
        return christoffel_blocks


# ── Helpers ──────────────────────────────────────────────────────────────

def _fresh_idx(index_type, base_name, is_up=True):
    """Create a fresh Index on a given type."""
    return Index(base_name, index_type, is_up)


def _is_zero(expr):
    """Check if an expression is zero."""
    if expr is None:
        return True
    if isinstance(expr, TensorProduct) and expr.coeff == 0:
        return True
    if isinstance(expr, ScalarExpr) and expr.expr == 0:
        return True
    if isinstance(expr, TensorSum) and len(expr.terms) == 0:
        return True
    return False


def _neg(expr):
    """Negate an expression, handling None."""
    if expr is None or _is_zero(expr):
        return None
    return -expr


def _add_exprs(*exprs):
    """Add expressions, skipping None/zero."""
    non_zero = [e for e in exprs if e is not None and not _is_zero(e)]
    if not non_zero:
        return None
    result = non_zero[0]
    for e in non_zero[1:]:
        result = result + e
    return result
