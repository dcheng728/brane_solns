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
from ..core.expr import (TensorAtom, TensorProduct, TensorSum, ScalarExpr,
                         replace_index)
from ..core.derivative import PartialDerivative
from ..core.metric import MetricTensor


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
    blocks : dict of (IndexType, IndexType) → callable, MetricTensor, or sp.Expr
        Metric blocks. Keys are pairs of child types.
        Values are either:
        - Callables ``(idx1, idx2) → TensorExpr`` (e.g. MetricTensor objects)
        - SymPy scalars (for constant scalar blocks)
    inv_blocks : dict, optional
        Inverse metric blocks. Same format. For MetricTensor blocks,
        use ``metric.inv`` as the value.
    scalar_fields : list of ScalarField
        All scalar fields appearing in block coefficients.
    truncation : set of IndexType
        Child types where partial derivatives vanish (e.g. the torus
        direction in F-theory where fields are T²-independent).
    """

    def __init__(self, split, blocks, inv_blocks=None, scalar_fields=None,
                 truncation=None):
        self.split = split
        self.children = split.children

        self.blocks = {}
        for key, val in blocks.items():
            self.blocks[key] = _as_block_callable(val)

        if inv_blocks is not None:
            self.inv_blocks = {}
            for key, val in inv_blocks.items():
                self.inv_blocks[key] = _as_block_callable(val)
        else:
            self.inv_blocks = None

        self.scalar_fields = list(scalar_fields or [])
        self.truncation = set(truncation or [])

        # Create partial derivative operators for each child type
        self._partials = {}
        for child in self.children:
            self._partials[child] = PartialDerivative(child, f'partial_{child.name}')

        # Counter for unique dummy names in _get_gamma
        self._dummy_counter = 0

    def get_block(self, type_row, type_col, idx1, idx2):
        """Get metric block g_{type_row, type_col}(idx1, idx2)."""
        block_fn = self.blocks.get((type_row, type_col))
        if block_fn is None:
            return TensorProduct(0, ())
        return block_fn(idx1, idx2)

    def get_inv_block(self, type_row, type_col, idx1, idx2):
        """Get inverse metric block g^{type_row, type_col}(idx1, idx2)."""
        if self.inv_blocks is None:
            raise ValueError("Inverse metric blocks not provided. "
                             "Pass inv_blocks to BlockMetric constructor.")
        block_fn = self.inv_blocks.get((type_row, type_col))
        if block_fn is None:
            return TensorProduct(0, ())
        return block_fn(idx1, idx2)

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

        Each block has canonical free indices: ``_ca`` (up), ``_cb`` (down),
        ``_cc`` (down).  Internal dummies use ``_cd``.

        Returns
        -------
        dict : (type_A, type_B, type_C) → TensorExpr
        """
        result = {}
        children = self.children

        for type_A, type_B, type_C in itertools.product(children, repeat=3):
            total = None

            # Canonical free indices for this block
            a = Index('_ca', type_A, is_up=True)
            b = Index('_cb', type_B, is_up=False)
            c = Index('_cc', type_C, is_up=False)

            for type_D in children:
                d_up = Index('_cd', type_D, is_up=True)
                d_down = Index('_cd', type_D, is_up=False)

                # g^{AD}(a, d_up)
                g_inv_AD = self.get_inv_block(type_A, type_D, a, d_up)
                if _is_zero(g_inv_AD):
                    continue

                # partial_B g_{DC}
                g_DC = self.get_block(type_D, type_C, d_down, c)
                term1 = self._differentiate(b, g_DC) if not _is_zero(g_DC) else None

                # partial_C g_{DB}
                g_DB = self.get_block(type_D, type_B, d_down, b)
                term2 = self._differentiate(c, g_DB) if not _is_zero(g_DB) else None

                # partial_D g_{BC}
                g_BC = self.get_block(type_B, type_C, b, c)
                term3 = self._differentiate(d_down, g_BC) if not _is_zero(g_BC) else None

                # Combine: 1/2 * g^{AD} * (term1 + term2 - term3)
                bracket = _add_exprs(term1, term2, _neg(term3))
                if bracket is None or _is_zero(bracket):
                    continue

                contribution = sp.Rational(1, 2) * g_inv_AD * bracket
                total = contribution if total is None else total + contribution

            result[(type_A, type_B, type_C)] = (
                total if total is not None else TensorProduct(0, ())
            )

        return result

    def riemann(self, christoffel_blocks=None):
        r"""Compute all Riemann tensor blocks :math:`R^A{}_{BCD}`.

        .. math::
            R^A{}_{BCD} = \partial_C \Gamma^A_{DB}
                        - \partial_D \Gamma^A_{CB}
                        + \Gamma^A_{CE} \Gamma^E_{DB}
                        - \Gamma^A_{DE} \Gamma^E_{CB}

        Parameters
        ----------
        christoffel_blocks : dict, optional
            Pre-computed Christoffel blocks from :meth:`christoffel`.

        Returns
        -------
        dict : (type_A, type_B, type_C, type_D) → TensorExpr
            Free indices: ``_ra`` (up/A), ``_rb`` (down/B),
            ``_rc`` (down/C), ``_rd`` (down/D).
        """
        if christoffel_blocks is None:
            christoffel_blocks = self.christoffel()

        self._dummy_counter = 0
        result = {}
        children = self.children

        for type_A, type_B, type_C, type_D in itertools.product(children, repeat=4):
            # Free indices for this Riemann block
            a = Index('_ra', type_A, is_up=True)
            b = Index('_rb', type_B, is_up=False)
            c = Index('_rc', type_C, is_up=False)
            d = Index('_rd', type_D, is_up=False)

            total = None

            # --- Term 1: +partial_c Gamma^a_{db} ---
            gamma_ADB = self._get_gamma(christoffel_blocks,
                                        type_A, type_D, type_B, a, d, b)
            if not _is_zero(gamma_ADB):
                term1 = self._differentiate(c, gamma_ADB)
                if not _is_zero(term1):
                    total = _add_exprs(total, term1)

            # --- Term 2: -partial_d Gamma^a_{cb} ---
            gamma_ACB = self._get_gamma(christoffel_blocks,
                                        type_A, type_C, type_B, a, c, b)
            if not _is_zero(gamma_ACB):
                term2 = self._differentiate(d, gamma_ACB)
                if not _is_zero(term2):
                    total = _add_exprs(total, _neg(term2))

            # --- Terms 3,4: Gamma*Gamma contractions ---
            for type_E in children:
                e_up = Index('_re', type_E, is_up=True)
                e_down = Index('_re', type_E, is_up=False)

                # +Gamma^a_{ce} Gamma^e_{db}
                gamma_ACE = self._get_gamma(christoffel_blocks,
                                            type_A, type_C, type_E,
                                            a, c, e_down)
                gamma_EDB = self._get_gamma(christoffel_blocks,
                                            type_E, type_D, type_B,
                                            e_up, d, b)
                if not _is_zero(gamma_ACE) and not _is_zero(gamma_EDB):
                    term3 = gamma_ACE * gamma_EDB
                    total = _add_exprs(total, term3)

                # -Gamma^a_{de} Gamma^e_{cb}
                gamma_ADE = self._get_gamma(christoffel_blocks,
                                            type_A, type_D, type_E,
                                            a, d, e_down)
                gamma_ECB = self._get_gamma(christoffel_blocks,
                                            type_E, type_C, type_B,
                                            e_up, c, b)
                if not _is_zero(gamma_ADE) and not _is_zero(gamma_ECB):
                    term4 = gamma_ADE * gamma_ECB
                    total = _add_exprs(total, _neg(term4))

            result[(type_A, type_B, type_C, type_D)] = (
                total if total is not None else TensorProduct(0, ())
            )

        return result

    def riemann_lower(self, riemann_blocks=None, christoffel_blocks=None):
        r"""Compute all-lower Riemann tensor blocks :math:`R_{ABCD}`.

        .. math::
            R_{ABCD} = g_{AE} \, R^E{}_{BCD}

        Parameters
        ----------
        riemann_blocks : dict, optional
        christoffel_blocks : dict, optional

        Returns
        -------
        dict : (type_A, type_B, type_C, type_D) → TensorExpr
        """
        if riemann_blocks is None:
            riemann_blocks = self.riemann(christoffel_blocks)

        result = {}
        children = self.children

        for type_A, type_B, type_C, type_D in itertools.product(children, repeat=4):
            total = None

            # Fresh index for the lowered slot
            a_down = Index('_rla', type_A, is_up=False)

            for type_E in children:
                e_up = Index('_rle', type_E, is_up=True)
                e_down = Index('_rle', type_E, is_up=False)

                # g_{AE}(a_down, e_down)
                g_AE = self.get_block(type_A, type_E, a_down, e_down)
                if _is_zero(g_AE):
                    continue

                # R^E_{BCD}
                r_block = riemann_blocks.get((type_E, type_B, type_C, type_D))
                if r_block is None or _is_zero(r_block):
                    continue

                # Rename R's up-index _ra → e_up for contraction with g_{AE}
                r_renamed = replace_index(
                    r_block, Index('_ra', type_E, True), e_up)

                contribution = g_AE * r_renamed
                total = contribution if total is None else total + contribution

            result[(type_A, type_B, type_C, type_D)] = (
                total if total is not None else TensorProduct(0, ())
            )

        return result

    def _get_gamma(self, blocks, tA, tB, tC, a, b, c):
        """Retrieve a Christoffel block with specific free indices.

        Renames the canonical free indices (_ca, _cb, _cc) to the requested
        indices (a, b, c) and gives internal dummies unique names to prevent
        clashes when blocks are multiplied in the Riemann formula.
        """
        block = blocks.get((tA, tB, tC))
        if block is None or _is_zero(block):
            return TensorProduct(0, ())

        # Rename canonical free indices → requested indices
        result = replace_index(block, Index('_ca', tA, True), a)
        result = replace_index(result, Index('_cb', tB, False), b)
        result = replace_index(result, Index('_cc', tC, False), c)

        # Rename internal dummies to unique names
        count = self._dummy_counter
        self._dummy_counter += 1
        new_name = f'_g{count}'
        for child in self.children:
            result = replace_index(result,
                                   Index('_cd', child, True),
                                   Index(new_name, child, True))
            result = replace_index(result,
                                   Index('_cd', child, False),
                                   Index(new_name, child, False))

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

def _as_block_callable(val):
    """Wrap a block value as a callable ``(idx1, idx2) → TensorExpr``."""
    if callable(val):
        return val
    if isinstance(val, ScalarExpr):
        s = val.expr
        return lambda i, j: ScalarExpr(s)
    scalar = sp.sympify(val)
    return lambda i, j: ScalarExpr(scalar)


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
