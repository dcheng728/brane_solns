"""Curvature tensors derived from a metric.

Provides RiemannGeometry: a convenient container that declares Christoffel
symbols, Riemann, Ricci, scalar curvature, Weyl, and Einstein tensors with
their correct symmetries and inter-relations.
"""

import sympy as sp
from .index import Index, indices
from .tensor_head import TensorHead
from .symmetry import TensorSymmetry
from .expr import TensorAtom, TensorProduct, TensorSum
from .metric import MetricTensor
from .derivative import CovDerivative


class RiemannGeometry:
    """Curvature tensor hierarchy for a given metric.

    Creates all standard curvature objects with correct symmetries:
    Christoffel, Riemann, Ricci, Weyl, scalar curvature, Einstein.

    Parameters
    ----------
    metric : MetricTensor
    index_type : IndexType

    Attributes
    ----------
    Christoffel : TensorHead  — Gamma^a_{bc}, symmetric in (b,c)
    Riemann : TensorHead      — R^a_{bcd}, Riemann symmetry
    Ricci : TensorHead         — R_{ab}, symmetric
    RicciScalar : TensorHead   — R (scalar, rank 0)
    Weyl : TensorHead          — C_{abcd}, Riemann symmetry + trace-free
    Einstein : TensorHead      — G_{ab}, symmetric
    D : CovDerivative          — Levi-Civita covariant derivative
    """

    def __init__(self, metric, index_type):
        self.metric = metric
        self.index_type = index_type

        # Christoffel symbol: Gamma^a_{bc}, symmetric in (b,c)
        # No full symmetry — only the lower pair is symmetric.
        # For now, declare no_symmetry(3) and rely on user convention.
        self.Christoffel = TensorHead(
            'Gamma', [index_type] * 3,
            TensorSymmetry.no_symmetry(3)
        )

        # Riemann tensor: R^a_{bcd} or R_{abcd}
        # We declare it with all-covariant Riemann symmetry
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

        # Weyl tensor: C_{abcd}, same symmetries as Riemann (+ trace-free)
        self.Weyl = TensorHead(
            'C', [index_type] * 4,
            TensorSymmetry.riemann()
        )

        # Einstein tensor: G_{ab}, symmetric
        self.Einstein = TensorHead(
            'G', [index_type] * 2,
            TensorSymmetry.fully_symmetric(2)
        )

        # Covariant derivative (Levi-Civita)
        self.D = CovDerivative(index_type, metric, name='D')
        self.D.set_riemann(self.Riemann)

    def ricci_from_riemann(self, a, b):
        """Return Ric_{ab} = R^c_{acb} as a substitution expression.

        Parameters
        ----------
        a, b : Index (covariant)

        Returns
        -------
        TensorProduct: contraction of Riemann
        """
        dummy = _fresh('c', a, b, self.index_type)
        c_up = Index(dummy, self.index_type, is_up=True)
        c_down = Index(dummy, self.index_type, is_up=False)
        return TensorProduct.from_atom(
            TensorAtom(self.Riemann, (c_up, a, c_down, b))
        )

    def scalar_from_ricci(self):
        """Return R = g^{ab} Ric_{ab}.

        Returns
        -------
        TensorProduct: contraction of metric inverse with Ricci
        """
        a_up = Index('a', self.index_type, is_up=True)
        a_down = Index('a', self.index_type, is_up=False)
        b_up = Index('b', self.index_type, is_up=True)
        b_down = Index('b', self.index_type, is_up=False)

        return (self.metric.inv_head(a_up, b_up)
                * self.Ricci(a_down, b_down))

    def einstein_from_ricci(self, a, b):
        """Return G_{ab} = Ric_{ab} - (1/2) R g_{ab}.

        Parameters
        ----------
        a, b : Index (covariant)
        """
        ric_term = self.Ricci(a, b)
        scalar_term = sp.Rational(-1, 2) * self.RicciScalar() * self.metric(a, b)
        return ric_term + scalar_term

    def weyl_decomposition(self, a, b, c, d):
        """Express Weyl in terms of Riemann, Ricci, scalar curvature.

        C_{abcd} = R_{abcd}
                 - 2/(D-2) (g_{a[c} Ric_{d]b} - g_{b[c} Ric_{d]a})
                 + 2/((D-1)(D-2)) R g_{a[c} g_{d]b}

        Parameters
        ----------
        a, b, c, d : Index (all covariant)
        """
        D = self.index_type.dim
        if D is None:
            D = sp.Symbol('D')

        R = self.Riemann
        Ric = self.Ricci
        g = self.metric
        Rsc = self.RicciScalar

        result = R(a, b, c, d)

        # -2/(D-2) * (g_{ac} Ric_{db} - g_{ad} Ric_{cb}
        #            - g_{bc} Ric_{da} + g_{bd} Ric_{ca})
        factor1 = sp.Rational(-1, 1) / (D - 2)
        result = result + factor1 * (
            g(a, c) * Ric(d, b) - g(a, d) * Ric(c, b)
            - g(b, c) * Ric(d, a) + g(b, d) * Ric(c, a)
        )

        # + 2/((D-1)(D-2)) * R * (g_{ac} g_{bd} - g_{ad} g_{bc})
        factor2 = sp.Rational(2, 1) / ((D - 1) * (D - 2))
        result = result + factor2 * Rsc() * (
            g(a, c) * g(b, d) - g(a, d) * g(b, c)
        )

        return result


def _fresh(base, idx_a, idx_b, index_type):
    """Pick a name not clashing with existing indices."""
    used = set()
    for i in [idx_a, idx_b]:
        if isinstance(i, Index):
            used.add(i.name)
    name = base
    counter = 0
    while name in used:
        name = f"{base}{counter}"
        counter += 1
    return name
