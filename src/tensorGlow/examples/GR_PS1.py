"""
GR Problem Sheet 1 — Tensor Manipulations (Section A)
Claudia de Rham, Imperial College, 2024

Solved symbolically using tensorGlow.
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

import sympy as sp
from tensorGlow import *


# ═══════════════════════════════════════════════════════════════════════
# Setup: 4d Minkowski spacetime
# ═══════════════════════════════════════════════════════════════════════

L = IndexType('Lorentz', dim=4, dummy_prefix='L')
mu, nu, rho, sigma = indices('mu nu rho sigma', L)
a, b, c, d = indices('a b c d', L)
g = MetricTensor(L, name='eta')


# ═══════════════════════════════════════════════════════════════════════
# Q3. Matrix multiplication in index notation
# ═══════════════════════════════════════════════════════════════════════

print("=" * 60)
print("Q3. Matrix multiplication in index notation")
print("=" * 60)

# Use a general IndexType for n-dimensional matrices
Mat = IndexType('Mat', dim=sp.Symbol('n'), dummy_prefix='k')
i, j, k, l = indices('i j k l', Mat)

A = TensorHead('A', [Mat, Mat])
B = TensorHead('B', [Mat, Mat])
C = TensorHead('C', [Mat, Mat])
D = TensorHead('D', [Mat, Mat])

# (1) A = B(C + D)  →  A_{ij} = B_{ik}(C_{kj} + D_{kj})
print("\n(1) A = B(C + D)")
expr1 = B(-i, -k) * (C(k, -j) + D(k, -j))
print(f"    A_{{ij}} = {expr1}")

# (2) A = BCD  →  A_{ij} = B_{ik} C_{kl} D_{lj}
print("\n(2) A = BCD")
expr2 = B(-i, -k) * C(k, -l) * D(l, -j)
print(f"    A_{{ij}} = {expr2}")


# ═══════════════════════════════════════════════════════════════════════
# Q5. Kronecker delta contractions
# ═══════════════════════════════════════════════════════════════════════

print("\n" + "=" * 60)
print("Q5. Kronecker delta contractions")
print("=" * 60)

N = IndexType('N', dim=sp.Symbol('n'), dummy_prefix='L')
a, b, c, d = indices('a b c d', N)
gN = MetricTensor(N, name='delta')

# delta^a_a = n  (self-contraction resolves automatically)
print(f"\n  delta^a_a = {gN.delta(a, -a)}")

# delta^a_b delta^b_a: contract deltas, then the trace resolves
expr_dd = gN.delta(a, -b) * gN.delta(b, -a)
print(f"  delta^a_b delta^b_a = {gN.contract_metrics(expr_dd)}")

# delta^a_b delta^b_c delta^c_d = delta^a_d
expr_ddd = gN.delta(a, -b) * gN.delta(b, -c) * gN.delta(c, -d)
print(f"  delta^a_b delta^b_c delta^c_d = {gN.contract_metrics(expr_ddd)}")


# ═══════════════════════════════════════════════════════════════════════
# Q6. Dummy index relabeling
# (Z_abc + Z_cab + Z_bca) X^a X^b X^c = 3 Z_abc X^a X^b X^c
# ═══════════════════════════════════════════════════════════════════════

print("\n" + "=" * 60)
print("Q6. Dummy index relabeling")
print("=" * 60)

L6 = IndexType('L', dim=sp.Symbol('n'), dummy_prefix='L')
a, b, c = indices('a b c', L6)
Z = TensorHead('Z', [L6, L6, L6])
X = TensorHead('X', [L6])

t1 = Z(-a, -b, -c) * X(a) * X(b) * X(c)
t2 = Z(-c, -a, -b) * X(a) * X(b) * X(c)
t3 = Z(-b, -c, -a) * X(a) * X(b) * X(c)

# Canonicalize each term — they should all become identical
c1 = t1.canon_bp()
c2 = t2.canon_bp()
c3 = t3.canon_bp()

print(f"\n  Z_abc X^a X^b X^c  canonicalizes to:  {c1}")
print(f"  Z_cab X^a X^b X^c  canonicalizes to:  {c2}")
print(f"  Z_bca X^a X^b X^c  canonicalizes to:  {c3}")
print(f"\n  All three identical? {c1.atoms == c2.atoms == c3.atoms and c1.coeff == c2.coeff == c3.coeff}")

# Sum collects to 3x
total = (t1 + t2 + t3).canon_bp()
print(f"  Sum = {total}")
print("  => Coefficient is 3, confirming the identity.")


# ═══════════════════════════════════════════════════════════════════════
# Q7. S_{ab} A^{ab} = 0  (symmetric × antisymmetric contraction)
# ═══════════════════════════════════════════════════════════════════════

print("\n" + "=" * 60)
print("Q7. Symmetric × antisymmetric contraction vanishes")
print("=" * 60)

L7 = IndexType('L', dim=sp.Symbol('n'), dummy_prefix='L')
a, b = indices('a b', L7)

S = TensorHead('S', [L7, L7], TensorSymmetry.fully_symmetric(2))
A = TensorHead('A', [L7, L7], TensorSymmetry.fully_antisymmetric(2))

result = (S(-a, -b) * A(a, b)).canon_bp()
print(f"\n  S_{{ab}} A^{{ab}} = {result}")
print("  => Vanishes identically, as expected.")

# Also verify: A^{ba} = g^{bc} g^{ad} A_{cd}
# Since g is symmetric and A is antisymmetric, A^{ab} = -A^{ba}
expr_swap = A(b, a).canon_bp()
print(f"\n  A^{{ba}} canonicalizes to: {expr_swap}")
print("  => Picks up a minus sign (antisymmetry).")
