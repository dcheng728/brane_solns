"""
Standard identities of Riemannian geometry, derived by tensorGlow.

Christoffel symbol, geodesic equation, Riemann tensor,
covariant derivative, commutator, curvature hierarchy,
metric compatibility, Leibniz rule.
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

import sympy as sp
from tensorGlow import *


D = sp.Symbol('D')
L = IndexType('Lorentz', dim=D, dummy_prefix='L')
a, b, c, d, e = indices('a b c d e', L)

g = MetricTensor(L, name='g')
geom = RiemannGeometry(g, L)
R = geom.Riemann
U = TensorHead('U', [L])
V = TensorHead('V', [L])
W = TensorHead('W', [L])


# Q3(a): Christoffel from metric
with Problem("Q3(a)", "Christoffel symbol") as p:
    p.answer("Gamma^c_{ab}", geom.christoffel_formula(c, -a, -b))


# Q3(b): Geodesic equation
with Problem("Q3(b)", "Geodesic equation") as p:
    p.step("Covariant derivative of U",
           "D_b U^a", geom.D(-b) * U(a))
    p.answer("U^b D_b U^a", U(b) * geom.D(-b) * U(a))


# Q4: 4-acceleration
with Problem("Q4", "4-acceleration") as p:
    p.answer("a^a = U^b D_b U^a", U(b) * geom.D(-b) * U(a))


# Riemann from Christoffels
with Problem("Riemann tensor", "From Christoffels") as p:
    p.answer("R^a_{bcd}", geom.riemann_formula(a, -b, -c, -d))


# Riemann symmetries
with Problem("Riemann symmetries", "") as p:
    p.step("Antisymmetry in last pair",
           "R_{abdc}", R(-a, -b, -d, -c).canon_bp())
    p.step("Pair exchange",
           "R_{cdab}", R(-c, -d, -a, -b).canon_bp())
    p.step("Vanishing sum",
           "R_{abcd} + R_{abdc}", (R(-a,-b,-c,-d) + R(-a,-b,-d,-c)).canon_bp())


# Commutator: computed by expanding D twice
with Problem("[D_a, D_b]", "Computed from D = partial + Gamma") as p:
    p.step("D_a D_b V^c",
           "D_a D_b V^c", geom.D(-a) * geom.D(-b) * V(c))
    p.step("D_b D_a V^c",
           "D_b D_a V^c", geom.D(-b) * geom.D(-a) * V(c))


# Curvature hierarchy
with Problem("Curvature hierarchy", "") as p:
    p.step("Ricci from Riemann",
           "Ric_{ab}", geom.ricci_from_riemann(-a, -b))
    p.step("Ricci scalar",
           "R", geom.scalar_from_ricci())
    p.step("Einstein tensor",
           "G_{ab}", geom.einstein_from_ricci(-a, -b))


# Metric compatibility: derived, not hardcoded
with Problem("Metric compatibility", "Derived from D = partial + Gamma") as p:
    p.answer("D_a g_{bc}", geom.D(-a) * g(-b, -c))


# Leibniz rule
with Problem("Leibniz rule", "") as p:
    p.answer("D_c(V^a V^b)", geom.D(-c) * (V(a) * V(b)))


# Symmetric x antisymmetric
with Problem("Symmetric x antisymmetric", "") as p:
    S = TensorHead('S', [L, L], TensorSymmetry.fully_symmetric(2))
    A = TensorHead('A', [L, L], TensorSymmetry.fully_antisymmetric(2))
    p.answer("S_{ab} A^{ab}", (S(-a, -b) * A(a, b)).canon_bp())
