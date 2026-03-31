"""Curvature tensors derived from a metric.

Provides RiemannGeometry: declares Christoffel, Riemann, Ricci, Weyl,
scalar curvature, Einstein — both as named TensorHeads and as defining
formulas expressed in tensorGlow's abstract index notation.
"""

import sympy as sp
from .index import Index, indices
from .tensor_head import TensorHead
from .symmetry import TensorSymmetry
from .expr import Tensor, Prod, Sum, Scalar
from .metric import MetricTensor
from .operator import Partial, CovariantD
from .derivative import PartialDerivative


class RiemannGeometry:
    """Curvature tensor hierarchy for a given metric.

    Declares all standard curvature objects with correct symmetries,
    and provides the defining formulas as abstract tensor expressions.

    Parameters
    ----------
    metric : MetricTensor
    index_type : IndexType
    """

    def __init__(self, metric, index_type):
        self.metric = metric
        self.index_type = index_type

        # Operators (new, first-class)
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
        r"""The Riemann tensor from Christoffel symbols:

        .. math::
            R^a{}_{bcd} = \partial_c \Gamma^a_{db}
                        - \partial_d \Gamma^a_{cb}
                        + \Gamma^a_{ce} \Gamma^e_{db}
                        - \Gamma^a_{de} \Gamma^e_{cb}

        Parameters
        ----------
        a : Index (contravariant, free)
        b, c, d : Index (covariant, free)

        Returns
        -------
        TensorExpr
        """
        e_name = _fresh_name('e', a, b, c, d)
        e_up = Index(e_name, self.index_type, is_up=True)
        e_down = Index(e_name, self.index_type, is_up=False)

        Gamma = self.Christoffel
        p = self.partial

        # partial_c Gamma^a_{db}
        term1 = p(-c) * Gamma(a, -d, -b)
        # partial_d Gamma^a_{cb}
        term2 = p(-d) * Gamma(a, -c, -b)
        # Gamma^a_{ce} Gamma^e_{db}
        term3 = Gamma(a, -c, -e_down) * Gamma(e_up, -d, -b)
        # Gamma^a_{de} Gamma^e_{cb}
        term4 = Gamma(a, -d, -e_down) * Gamma(e_up, -c, -b)

        return term1 - term2 + term3 - term4

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
