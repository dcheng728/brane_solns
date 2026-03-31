"""
GR Problem Sheet 1 -- Section A: Tensor Manipulations
Claudia de Rham, Imperial College, 2024
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

import sympy as sp
from tensorGlow import *


n = sp.Symbol('n')
N = IndexType('N', dim=n, dummy_prefix='L')
a, b, c, d = indices('a b c d', N)


with Problem("Q3", "Matrix multiplication in index notation") as p:
    Mat = IndexType('Mat', dim=n, dummy_prefix='k')
    i, j, k, l = indices('i j k l', Mat)
    B = TensorHead('B', [Mat, Mat])
    C = TensorHead('C', [Mat, Mat])
    D = TensorHead('D', [Mat, Mat])

    p.step("A = B(C + D)",
           "A_{ij}", B(-i, -k) * (C(k, -j) + D(k, -j)))
    p.step("A = BCD",
           "A_{ij}", B(-i, -k) * C(k, -l) * D(l, -j))


with Problem("Q5", "Kronecker delta contractions") as p:
    g = MetricTensor(N, name='delta')

    p.step("Trace",
           "delta^a_a", g.delta(a, -a))
    p.step("Double contraction",
           "delta^a_b delta^b_a",
           g.contract_metrics(g.delta(a, -b) * g.delta(b, -a)))
    p.step("Triple contraction",
           "delta^a_b delta^b_c delta^c_d",
           g.contract_metrics(g.delta(a, -b) * g.delta(b, -c) * g.delta(c, -d)))


with Problem("Q6", "Dummy index relabeling") as p:
    Z = TensorHead('Z', [N, N, N])
    X = TensorHead('X', [N])

    t1 = (Z(-a, -b, -c) * X(a) * X(b) * X(c)).canon_bp()
    t2 = (Z(-c, -a, -b) * X(a) * X(b) * X(c)).canon_bp()
    t3 = (Z(-b, -c, -a) * X(a) * X(b) * X(c)).canon_bp()

    p.step("Canonicalize Z_{abc} X^a X^b X^c", rhs=t1)
    p.step("Canonicalize Z_{cab} X^a X^b X^c", rhs=t2)
    p.step("Canonicalize Z_{bca} X^a X^b X^c", rhs=t3)
    p.answer("(Z_{abc} + Z_{cab} + Z_{bca}) X^a X^b X^c",
             (t1 + t2 + t3).canon_bp())


with Problem("Q7", "Symmetric x antisymmetric contraction") as p:
    S = TensorHead('S', [N, N], TensorSymmetry.fully_symmetric(2))
    A = TensorHead('A', [N, N], TensorSymmetry.fully_antisymmetric(2))

    p.step("Antisymmetry", "A^{ba}", A(b, a).canon_bp())
    p.answer("S_{ab} A^{ab}", (S(-a, -b) * A(a, b)).canon_bp())
