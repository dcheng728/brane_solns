"""Curvature tensors derived from a metric.

Provides RiemannGeometry: declares Christoffel, Riemann, Ricci, Weyl,
scalar curvature, Einstein — both as named TensorHeads and as defining
formulas expressed in tensorGlow's abstract index notation.

Supports dimensional reduction via index splitting, metric truncation,
and partial derivative truncation.
"""

import itertools
import sympy as sp
from .index import Index, indices
from .tensor_head import TensorHead
from .symmetry import TensorSymmetry
from .expr import Tensor, Apply, Prod, Sum, Scalar, replace_index
from .metric import MetricTensor
from .operator import Partial, CovariantD
from .derivative import PartialDerivative


class RiemannGeometry:
    """Curvature tensor hierarchy for a given metric.

    Declares all standard curvature objects with correct symmetries,
    and provides the defining formulas as abstract tensor expressions.

    Supports dimensional reduction: split a parent index type into
    child types, set metric cross-blocks to zero, and truncate
    partial derivatives along certain directions.

    Parameters
    ----------
    metric : MetricTensor
    index_type : IndexType
    """

    def __init__(self, metric, index_type, coordinates=None, metric_components=None):
        self.metric = metric
        self.index_type = index_type

        # Concrete coordinate evaluation (optional)
        self._concrete = None
        if coordinates is not None and metric_components is not None:
            from .coordinates import Geometry
            self._concrete = Geometry(coordinates, metric_components)

        # Operators (first-class)
        self.partial = Partial(index_type)

        # Christoffel symbol: Gamma^a_{bc}, symmetric in (b,c)
        self.Christoffel = TensorHead(
            'Gamma', [index_type] * 3,
            TensorSymmetry.no_symmetry(3)
        )

        # Covariant derivative: defined from partial + Christoffel
        self.D = CovariantD(index_type, metric,
                            christoffel=self.Christoffel,
                            partial_op=self.partial)

        # Riemann tensor: R_{abcd} with Riemann symmetry
        self.Riemann = TensorHead(
            'R', [index_type] * 4,
            TensorSymmetry.riemann()
        )

        # Ricci tensor: R_{ab}, symmetric
        self.Ricci = TensorHead(
            'Ric', [index_type] * 2,
            TensorSymmetry.fully_symmetric(2)
        )

        # Ricci scalar: R (rank 0)
        self.RicciScalar = TensorHead(
            'R', [],
            TensorSymmetry.no_symmetry(0)
        )

        # Weyl tensor: C_{abcd}, Riemann symmetry
        self.Weyl = TensorHead(
            'C', [index_type] * 4,
            TensorSymmetry.riemann()
        )

        # Einstein tensor: G_{ab}, symmetric
        self.Einstein = TensorHead(
            'G', [index_type] * 2,
            TensorSymmetry.fully_symmetric(2)
        )

        # Wire up Riemann for commutator
        self.D.set_riemann(self.Riemann)

        # Dimensional reduction settings
        self._split = None
        self._metric_zeros = set()   # set of (IndexType, IndexType)
        self._partial_zeros = set()  # set of IndexType

    # ──────────────────────────────────────────────────────────────
    # Dimensional reduction: index splitting and truncation
    # ──────────────────────────────────────────────────────────────

    def split_indices(self, split):
        """Declare an index split for dimensional reduction.

        Parameters
        ----------
        split : IndexSplit
            Decomposition of the parent index type into child types.
        """
        self._split = split
        return self

    def set_metric_zero(self, type1, type2):
        """Declare that the metric block g_{type1, type2} vanishes.

        Applies to both g and g^{-1} with those child types.

        Parameters
        ----------
        type1, type2 : IndexType
            Child index types for which the metric block is zero.
        """
        self._metric_zeros.add((type1, type2))
        return self

    def set_partial_zero(self, child_type):
        """Declare that partial derivatives along child_type vanish.

        E.g. for truncation d_a = 0 on a torus direction.

        Parameters
        ----------
        child_type : IndexType
        """
        self._partial_zeros.add(child_type)
        return self

    def decompose(self, expr):
        """Decompose an expression into blocks by splitting parent-type indices.

        Each parent-type dummy pair is expanded into a sum over child types.
        Each parent-type free index is assigned to a child type (the block key).
        Zero rules (metric cross-blocks, truncated partials) eliminate terms.

        Parameters
        ----------
        expr : Expr
            An abstract tensor expression with parent-type indices.

        Returns
        -------
        dict : block_key → Expr
            block_key is a tuple of child IndexTypes, one per free
            parent-type index (in order of first appearance).
            Only non-zero blocks are included.
        """
        if self._split is None:
            raise ValueError("No index split defined. Call split_indices() first.")
        return _decompose_expr(
            expr, self._split, self.metric,
            self._metric_zeros, self._partial_zeros
        )

    # ──────────────────────────────────────────────────────────────
    # Defining formulas as abstract tensor expressions
    # ──────────────────────────────────────────────────────────────

    def christoffel_formula(self, c, a, b):
        r"""The Christoffel symbol from the metric:

        .. math::
            \Gamma^c_{ab} = \frac{1}{2} g^{cd}
            (\partial_a g_{db} + \partial_b g_{ad} - \partial_d g_{ab})

        Parameters
        ----------
        c : Index (contravariant, free)
        a, b : Index (covariant, free)

        Returns
        -------
        TensorExpr
        """
        d_name = _fresh_name('d', c, a, b)
        d_up = Index(d_name, self.index_type, is_up=True)
        d_down = Index(d_name, self.index_type, is_up=False)

        g = self.metric
        p = self.partial

        # partial_a g_{db}
        term1 = p(a) * g(d_down, b)
        # partial_b g_{ad}
        term2 = p(b) * g(a, d_down)
        # partial_d g_{ab}
        term3 = p(d_down) * g(a, b)

        return sp.Rational(1, 2) * g.inv(c, d_up) * (term1 + term2 - term3)

    def geodesic_equation(self, a, U):
        r"""The geodesic equation:

        .. math::
            \frac{D U^a}{D\tau} = U^b \nabla_b U^a = 0

        Equivalently:
        .. math::
            \ddot{x}^a + \Gamma^a_{bc} \dot{x}^b \dot{x}^c = 0

        Parameters
        ----------
        a : Index (contravariant, free)
        U : TensorHead (rank-1, the tangent vector / 4-velocity)

        Returns
        -------
        TensorExpr : the acceleration (should vanish for geodesics)
        """
        b_name = _fresh_name('b', a)
        c_name = _fresh_name('c', a, Index(b_name, self.index_type))
        b_up = Index(b_name, self.index_type, is_up=True)
        b_down = Index(b_name, self.index_type, is_up=False)
        c_up = Index(c_name, self.index_type, is_up=True)
        c_down = Index(c_name, self.index_type, is_up=False)

        Gamma = self.Christoffel

        # Gamma^a_{bc} U^b U^c
        return Gamma(a, -b_down, -c_down) * U(b_up) * U(c_up)

    def riemann_formula(self, a, b, c, d):
        r"""The Riemann tensor from Christoffel symbols.

        Always computes the canonical form first:

        .. math::
            R^{\alpha}{}_{\beta\gamma\delta}
              = \partial_\gamma \Gamma^\alpha_{\delta\beta}
              - \partial_\delta \Gamma^\alpha_{\gamma\beta}
              + \Gamma^\alpha_{\gamma\epsilon} \Gamma^\epsilon_{\delta\beta}
              - \Gamma^\alpha_{\delta\epsilon} \Gamma^\epsilon_{\gamma\beta}

        then raises/lowers indices with the metric to match the
        requested index structure.  This is necessary because Christoffel
        symbols are not tensors — their slot positions in the definition
        are fixed and cannot be raised/lowered.

        Parameters
        ----------
        a, b, c, d : Index
            Free indices with any variance.  The result has the index
            structure implied by the ``is_up`` flags.

        Returns
        -------
        TensorExpr
        """
        user = [a, b, c, d]
        canonical_up = [True, False, False, False]  # R^._{ . . . }

        # Fresh internal names that won't collide with user names
        used = {idx.name for idx in user}
        internal_names = []
        for base in ['_r', '_s', '_t', '_u']:
            name = base
            cnt = 0
            while name in used:
                name = f"{base}{cnt}"
                cnt += 1
            used.add(name)
            internal_names.append(name)

        i0_up   = Index(internal_names[0], self.index_type, is_up=True)
        i1_down = Index(internal_names[1], self.index_type, is_up=False)
        i2_down = Index(internal_names[2], self.index_type, is_up=False)
        i3_down = Index(internal_names[3], self.index_type, is_up=False)

        e_name = _fresh_name('e', *user, i0_up, i1_down, i2_down, i3_down)
        e_up   = Index(e_name, self.index_type, is_up=True)
        e_down = Index(e_name, self.index_type, is_up=False)

        Gamma = self.Christoffel
        p = self.partial

        # R^{i0}_{i1 i2 i3}  (canonical form)
        term1 = p(i2_down) * Gamma(i0_up, i3_down, i1_down)
        term2 = p(i3_down) * Gamma(i0_up, i2_down, i1_down)
        term3 = Gamma(i0_up, i2_down, e_down) * Gamma(e_up, i3_down, i1_down)
        term4 = Gamma(i0_up, i3_down, e_down) * Gamma(e_up, i2_down, i1_down)

        result = term1 - term2 + term3 - term4

        # Adjust each slot to match the user's requested variance
        g = self.metric
        for slot in range(4):
            iname = internal_names[slot]
            want_up = user[slot].is_up
            have_up = canonical_up[slot]

            if want_up == have_up:
                # Variance matches — rename internal index to user's name
                result = replace_index(
                    result,
                    Index(iname, self.index_type, is_up=True),
                    Index(user[slot].name, self.index_type, is_up=True))
                result = replace_index(
                    result,
                    Index(iname, self.index_type, is_up=False),
                    Index(user[slot].name, self.index_type, is_up=False))
            elif have_up and not want_up:
                # Canonical is up, user wants down → lower with g_{user, internal}
                i_down = Index(iname, self.index_type, is_up=False)
                result = g(user[slot], i_down) * result
            else:
                # Canonical is down, user wants up → raise with g^{user, internal}
                i_up = Index(iname, self.index_type, is_up=True)
                result = g.inv(user[slot], i_up) * result

        return result

    def ricci_from_riemann(self, a, b):
        r"""Ricci tensor as contraction of Riemann:

        .. math::
            \text{Ric}_{ab} = R^c{}_{acb}

        Parameters
        ----------
        a, b : Index (covariant)
        """
        dummy = _fresh_name('c', a, b)
        c_up = Index(dummy, self.index_type, is_up=True)
        c_down = Index(dummy, self.index_type, is_up=False)
        return Prod(1, (Tensor(self.Riemann, (c_up, a, c_down, b)),))

    def scalar_from_ricci(self):
        r"""Ricci scalar:

        .. math::
            R = g^{ab} \text{Ric}_{ab}
        """
        a_up = Index('a', self.index_type, is_up=True)
        a_down = Index('a', self.index_type, is_up=False)
        b_up = Index('b', self.index_type, is_up=True)
        b_down = Index('b', self.index_type, is_up=False)
        return (self.metric.inv_head(a_up, b_up)
                * self.Ricci(a_down, b_down))

    def einstein_from_ricci(self, a, b):
        r"""Einstein tensor:

        .. math::
            G_{ab} = \text{Ric}_{ab} - \frac{1}{2} R \, g_{ab}
        """
        ric_term = self.Ricci(a, b)
        scalar_term = sp.Rational(-1, 2) * self.RicciScalar() * self.metric(a, b)
        return ric_term + scalar_term

    def weyl_decomposition(self, a, b, c, d):
        r"""Weyl tensor in terms of Riemann, Ricci, scalar curvature:

        .. math::
            C_{abcd} = R_{abcd}
                     - \frac{2}{D-2} (g_{a[c} \text{Ric}_{d]b}
                       - g_{b[c} \text{Ric}_{d]a})
                     + \frac{2}{(D-1)(D-2)} R \, g_{a[c} g_{d]b}
        """
        D = self.index_type.dim
        if D is None:
            D = sp.Symbol('D')

        R = self.Riemann
        Ric = self.Ricci
        g = self.metric
        Rsc = self.RicciScalar

        result = R(a, b, c, d)
        factor1 = sp.Rational(-1, 1) / (D - 2)
        result = result + factor1 * (
            g(a, c) * Ric(d, b) - g(a, d) * Ric(c, b)
            - g(b, c) * Ric(d, a) + g(b, d) * Ric(c, a)
        )
        factor2 = sp.Rational(2, 1) / ((D - 1) * (D - 2))
        result = result + factor2 * Rsc() * (
            g(a, c) * g(b, d) - g(a, d) * g(b, c)
        )
        return result


    # ──────────────────────────────────────────────────────────────
    # Concrete evaluation (when coordinates + metric are provided)
    # ──────────────────────────────────────────────────────────────

    def evaluate_christoffel(self):
        """Evaluate Christoffel symbols with concrete coordinates.

        Requires coordinates and metric_components passed to __init__.
        """
        if self._concrete is None:
            raise ValueError("No coordinates defined. Pass coordinates= and "
                             "metric_components= to RiemannGeometry.")
        return self._concrete.christoffel()

    def evaluate_riemann(self):
        """Evaluate Riemann tensor with concrete coordinates."""
        if self._concrete is None:
            raise ValueError("No coordinates defined.")
        return self._concrete.riemann()

    def evaluate_ricci(self):
        """Evaluate Ricci tensor with concrete coordinates."""
        if self._concrete is None:
            raise ValueError("No coordinates defined.")
        return self._concrete.ricci()

    def evaluate_ricci_scalar(self):
        """Evaluate Ricci scalar with concrete coordinates."""
        if self._concrete is None:
            raise ValueError("No coordinates defined.")
        return self._concrete.ricci_scalar()

    def evaluate_geodesic(self):
        """Evaluate geodesic equations with concrete coordinates."""
        if self._concrete is None:
            raise ValueError("No coordinates defined.")
        return self._concrete.geodesic_equations()


def _fresh_name(base, *existing_indices):
    """Pick a name not clashing with existing indices."""
    used = set()
    for i in existing_indices:
        if isinstance(i, Index):
            used.add(i.name)
    name = base
    counter = 0
    while name in used:
        name = f"{base}{counter}"
        counter += 1
    return name


# ═════════════════════════════════════════════════════════════════════
# Decomposition helpers
# ═════════════════════════════════════════════════════════════════════

def _decompose_expr(expr, split, metric, metric_zeros, partial_zeros):
    """Decompose an expression by splitting parent-type indices into child types."""
    if isinstance(expr, Sum):
        # Determine canonical free index names from the first term
        # (all terms in a valid Sum should have the same free indices)
        free_names = _get_free_names(expr, split.parent)
        result = {}
        for term in expr.terms:
            blocks = _decompose_expr_inner(
                term, split, metric, metric_zeros, partial_zeros, free_names)
            for key, val in blocks.items():
                if key in result:
                    result[key] = result[key] + val
                else:
                    result[key] = val
        return result

    return _decompose_expr_inner(expr, split, metric, metric_zeros, partial_zeros)


def _decompose_expr_inner(expr, split, metric, metric_zeros, partial_zeros,
                          free_names=None):
    """Inner decomposition, with optional canonical free_names from parent Sum."""
    if isinstance(expr, Sum):
        if free_names is None:
            free_names = _get_free_names(expr, split.parent)
        result = {}
        for term in expr.terms:
            blocks = _decompose_expr_inner(
                term, split, metric, metric_zeros, partial_zeros, free_names)
            for key, val in blocks.items():
                if key in result:
                    result[key] = result[key] + val
                else:
                    result[key] = val
        return result

    if isinstance(expr, Prod):
        return _decompose_prod(expr, split, metric, metric_zeros, partial_zeros,
                               free_names)

    if isinstance(expr, (Tensor, Apply)):
        return _decompose_prod(Prod(1, (expr,)), split, metric, metric_zeros,
                               partial_zeros, free_names)

    if isinstance(expr, Scalar):
        return {(): expr}

    raise TypeError(f"Cannot decompose {type(expr)}")


def _get_free_names(expr, parent):
    """Get sorted free index names of the parent type from an expression."""
    free = expr.free_indices
    return sorted(set(idx.name for idx in free if idx.index_type is parent))


def _decompose_prod(prod, split, metric, metric_zeros, partial_zeros,
                    free_names=None):
    """Decompose a Prod by enumerating all child-type assignments.

    Parameters
    ----------
    free_names : list of str, optional
        Canonical ordering of free index names. If provided, ensures
        consistent block keys across terms in a Sum.
    """
    parent = split.parent
    children = split.children

    # Collect all parent-type indices from the full expression tree
    all_idx = _collect_all_indices(prod)
    parent_idx = [idx for idx in all_idx if idx.index_type is parent]

    if not parent_idx:
        # No parent indices — return as-is
        return {(): prod}

    # Group by name, classify as free or dummy
    by_name = {}
    for idx in parent_idx:
        by_name.setdefault(idx.name, []).append(idx)

    if free_names is None:
        free_names = []
    else:
        free_names = list(free_names)
    dummy_names = []
    for name in sorted(by_name.keys()):
        occurrences = by_name[name]
        ups = [i for i in occurrences if i.is_up]
        downs = [i for i in occurrences if not i.is_up]
        if len(ups) == 1 and len(downs) == 1:
            dummy_names.append(name)
        elif name not in free_names:
            free_names.append(name)
    free_names.sort()

    n_children = len(children)
    result = {}

    for free_assign in itertools.product(range(n_children), repeat=len(free_names)):
        for dummy_assign in itertools.product(range(n_children), repeat=len(dummy_names)):
            # Build name → child_type mapping
            name_to_child = {}
            for i, name in enumerate(free_names):
                name_to_child[name] = children[free_assign[i]]
            for i, name in enumerate(dummy_names):
                name_to_child[name] = children[dummy_assign[i]]

            # Replace all parent-type indices with child-type versions
            new_prod = _replace_parent_type(prod, parent, name_to_child)

            # Check zero rules — if any factor is zero, skip
            if _has_zero(new_prod, metric, metric_zeros, partial_zeros):
                continue

            # Block key: tuple of child types for each free index
            block_key = tuple(children[free_assign[i]] for i in range(len(free_names)))

            if block_key in result:
                result[block_key] = result[block_key] + new_prod
            else:
                result[block_key] = new_prod

    return result


def _collect_all_indices(expr):
    """Collect all Index objects from an expression tree, in order."""
    if isinstance(expr, Tensor):
        return list(expr.indices)
    if isinstance(expr, Apply):
        return [expr.deriv_index] + _collect_all_indices(expr.operand)
    if isinstance(expr, Prod):
        result = []
        for f in expr.factors:
            result.extend(_collect_all_indices(f))
        return result
    if isinstance(expr, Scalar):
        return []
    return []


def _replace_parent_type(expr, parent, name_to_child):
    """Replace all parent-type indices with child-type versions.

    Parameters
    ----------
    expr : Expr
    parent : IndexType
        The parent index type to replace.
    name_to_child : dict
        Maps index name → child IndexType.
    """
    if isinstance(expr, Tensor):
        new_indices = tuple(
            Index(idx.name, name_to_child[idx.name], idx.is_up)
            if idx.index_type is parent and idx.name in name_to_child
            else idx
            for idx in expr.indices
        )
        return Tensor(expr.head, new_indices)

    if isinstance(expr, Apply):
        di = expr.deriv_index
        if di.index_type is parent and di.name in name_to_child:
            new_di = Index(di.name, name_to_child[di.name], di.is_up)
        else:
            new_di = di
        new_operand = _replace_parent_type(expr.operand, parent, name_to_child)
        return Apply(expr.op, new_di, new_operand)

    if isinstance(expr, Prod):
        new_factors = tuple(
            _replace_parent_type(f, parent, name_to_child) for f in expr.factors
        )
        return Prod(expr.coeff, new_factors)

    if isinstance(expr, Sum):
        new_terms = tuple(
            _replace_parent_type(t, parent, name_to_child) for t in expr.terms
        )
        return Sum(new_terms)

    return expr


def _has_zero(expr, metric, metric_zeros, partial_zeros):
    """Check if any zero rule makes this expression vanish.

    Returns True if a metric cross-block is zero or a partial
    derivative is along a truncated direction.
    """
    if isinstance(expr, Tensor):
        # Check metric and inverse metric for zero cross-blocks
        if expr.head is metric.head or expr.head is metric.inv_head:
            if len(expr.indices) == 2:
                t1 = expr.indices[0].index_type
                t2 = expr.indices[1].index_type
                if (t1, t2) in metric_zeros:
                    return True
        return False

    if isinstance(expr, Apply):
        # Partial derivative along a truncated direction → zero
        if expr.deriv_index.index_type in partial_zeros:
            return True
        # Also check if the operand itself is zero
        return _has_zero(expr.operand, metric, metric_zeros, partial_zeros)

    if isinstance(expr, Prod):
        # If any factor is zero, the whole product is zero
        return any(_has_zero(f, metric, metric_zeros, partial_zeros) for f in expr.factors)

    return False
