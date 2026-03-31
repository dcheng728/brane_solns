"""
GR Problem Sheet 2 -- Tensors
Claudia de Rham, Imperial College, 2024

Q3: Christoffel symbols and geodesic equation in polar coordinates
Q4: Uniformly accelerated particle in Minkowski space
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

import sympy as sp
from tensorGlow import *


# =====================================================================
# Q3: Life in Polar Coordinates
# =====================================================================

r, theta = sp.symbols('r theta')

L = IndexType('polar', dim=2, dummy_prefix='L')
g = MetricTensor(L, name='g')

geom = RiemannGeometry(g, L,
    coordinates=(r, theta),
    metric_components={
        ('r', 'theta'): 0,
        ('r', 'r'): 1,
        ('theta', 'theta'): r**2,
    }
)

with Problem("Q3(a)", "Christoffel symbols in polar coordinates") as p:
    chris = geom.evaluate_christoffel()
    for (c, a, b), val in sorted(chris.items()):
        p.step(f"Gamma^{c}_{{{a},{b}}}", rhs=val)

with Problem("Q3(b)", "Geodesic equation in polar coordinates") as p:
    for name, eq in geom.evaluate_geodesic():
        p.step(f"{name} component", rhs=f"{eq} = 0")


# =====================================================================
# Q4: Particle in Acceleration
# =====================================================================

tau = sp.Symbol('tau')
a_const = sp.Symbol('a', positive=True)

t_traj = sp.sinh(a_const * tau) / a_const
x_traj = sp.cosh(a_const * tau) / a_const

with Problem("Q4(a)", "Hyperbola") as p:
    p.answer("x^2 - t^2", sp.simplify(x_traj**2 - t_traj**2))

with Problem("Q4(b)", "4-velocity and norm") as p:
    Ut = sp.diff(t_traj, tau)
    Ux = sp.diff(x_traj, tau)
    p.step("U^t", rhs=sp.simplify(Ut))
    p.step("U^x", rhs=sp.simplify(Ux))
    p.answer("U.U = eta_{mu nu} U^mu U^nu", sp.simplify(-Ut**2 + Ux**2))

with Problem("Q4(c)", "4-acceleration and norm") as p:
    at = sp.diff(t_traj, tau, 2)
    ax = sp.diff(x_traj, tau, 2)
    p.step("a^t", rhs=sp.simplify(at))
    p.step("a^x", rhs=sp.simplify(ax))
    p.answer("a.a = eta_{mu nu} a^mu a^nu", sp.simplify(-at**2 + ax**2))
    p.step("U.a", rhs=sp.simplify(-Ut * at + Ux * ax))
