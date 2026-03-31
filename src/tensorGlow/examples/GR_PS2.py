"""
GR Problem Sheet 2 -- Tensors
Claudia de Rham, Imperial College, 2024

Q3: Christoffel symbols and geodesic equation
Q4: 4-acceleration, covariant derivative, commutator
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


with Problem("Q3(a)", "Christoffel symbol") as p:
    p.answer("Gamma^c_{ab}", geom.christoffel_formula(c, -a, -b))

with Problem("Q3(b)", "Geodesic equation") as p:
    p.step("Covariant derivative of tangent vector",
           "D_b U^a", geom.D(-b, U(a)))
    p.answer("U^b D_b U^a", U(b) * geom.D(-b, U(a)))

with Problem("Q4", "4-acceleration") as p:
    p.answer("a^a = U^b D_b U^a", U(b) * geom.D(-b, U(a)))

with Problem("Riemann tensor", "From Christoffels") as p:
    p.answer("R^a_{bcd}", geom.riemann_formula(a, -b, -c, -d))

with Problem("Riemann symmetries", "") as p:
    p.step("Antisymmetry in last pair",
           "R_{abdc}", R(-a, -b, -d, -c).canon_bp())
    p.step("Pair exchange",
           "R_{cdab}", R(-c, -d, -a, -b).canon_bp())
    p.step("Vanishing sum",
           "R_{abcd} + R_{abdc}", (R(-a,-b,-c,-d) + R(-a,-b,-d,-c)).canon_bp())

with Problem("Commutator [D_a, D_b]", "") as p:
    p.step("On a contravariant vector",
           "[D_a, D_b] V^c", geom.D.commutator(a, b, V(c)))
    p.step("On a covariant vector",
           "[D_a, D_b] W_c", geom.D.commutator(a, b, W(-c)))

with Problem("Curvature hierarchy", "") as p:
    p.step("Ricci from Riemann",
           "Ric_{ab}", geom.ricci_from_riemann(-a, -b))
    p.step("Ricci scalar",
           "R", geom.scalar_from_ricci())
    p.step("Einstein tensor",
           "G_{ab}", geom.einstein_from_ricci(-a, -b))

with Problem("Identities", "") as p:
    p.step("Metric compatibility",
           "D_a g_{bc}", geom.D(-a, g(-b, -c)))
    p.step("Leibniz rule",
           "D_c(V^a V^b)", geom.D(-c, V(a) * V(b)))
    S = TensorHead('S', [L, L], TensorSymmetry.fully_symmetric(2))
    A = TensorHead('A', [L, L], TensorSymmetry.fully_antisymmetric(2))
    p.step("Symmetric x antisymmetric",
           "S_{ab} A^{ab}", (S(-a, -b) * A(a, b)).canon_bp())
