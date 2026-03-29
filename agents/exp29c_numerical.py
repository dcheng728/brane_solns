"""Exp29c: Numerical verification of KK formula for (1,1)-string."""
import sys; sys.path.insert(0, 'src')
import sympy as sp
from sympy import Rational as R, cancel
from sugra import Metric

y0, y1 = sp.symbols('y0 y1', real=True)
t, x1 = sp.symbols('t x1', real=True)
z1s, z2s = sp.symbols('z1 z2', real=True)

# Explicit H — no H(r), just a rational function
H = 1 + 1/(y0**2 + y1**2)**3

M0 = sp.Matrix([[H**R(1,2), 0], [0, H**R(-1,2)]])
Lambda = sp.Matrix([[1, 0], [1, 1]])
M = Lambda * M0 * Lambda.T
M_inv = M.inv()

coords = [t, x1, z1s, z2s, y0, y1]

# F1 6d metric
g_f1 = sp.zeros(6, 6)
g_f1[0,0] = -1/H; g_f1[1,1] = 1/H
g_f1[2,2] = M0[0,0]; g_f1[3,3] = M0[1,1]
g_f1[4,4] = 1; g_f1[5,5] = 1

# (1,1) 6d metric
g_11 = sp.zeros(6, 6)
g_11[0,0] = -1/H; g_11[1,1] = 1/H
g_11[2,2] = M[0,0]; g_11[2,3] = M[0,1]
g_11[3,2] = M[1,0]; g_11[3,3] = M[1,1]
g_11[4,4] = 1; g_11[5,5] = 1

print("Computing F1 6d Ricci...")
Ric_f1 = Metric(g_f1, coords).ricci_tensor(simplify_func=cancel)
print("Computing (1,1) 6d Ricci...")
Ric_11 = Metric(g_11, coords).ricci_tensor(simplify_func=cancel)

# 4d base Ricci
g4 = sp.zeros(4,4)
g4[0,0] = -1/H; g4[1,1] = 1/H; g4[2,2] = 1; g4[3,3] = 1
print("Computing 4d Ricci...")
Ric_4 = Metric(g4, [t, x1, y0, y1]).ricci_tensor(simplify_func=cancel)

# KK term
dMinv_y0 = sp.Matrix([[sp.diff(M_inv[i,j], y0) for j in range(2)] for i in range(2)])
dMinv_y1 = sp.Matrix([[sp.diff(M_inv[i,j], y1) for j in range(2)] for i in range(2)])
dM_y0 = sp.Matrix([[sp.diff(M[i,j], y0) for j in range(2)] for i in range(2)])
dM_y1 = sp.Matrix([[sp.diff(M[i,j], y1) for j in range(2)] for i in range(2)])

kk_01 = R(1,4) * (dMinv_y0 * dM_y1).trace()
kk_00 = R(1,4) * (dMinv_y0 * dM_y0).trace()

pt = {y0: R(3,2), y1: R(1,2)}

print("\nNumerical check at y0=3/2, y1=1/2:")
for label, i6, j6, i4, j4, kk in [
    ("y0,y0", 4, 4, 2, 2, kk_00),
    ("y0,y1", 4, 5, 2, 3, kk_01)
]:
    v_f1 = float(Ric_f1[i6, j6].subs(pt))
    v_11 = float(Ric_11[i6, j6].subs(pt))
    v_4 = float(Ric_4[i4, j4].subs(pt))
    v_kk = float(kk.subs(pt))
    print(f"\n[{label}]:")
    print(f"  F1 Ric6   = {v_f1:.10f}")
    print(f"  (1,1) Ric6 = {v_11:.10f}")
    print(f"  F1 vs (1,1) diff = {v_11-v_f1:.2e}")
    print(f"  R4        = {v_4:.10f}")
    print(f"  KK        = {v_kk:.10f}")
    print(f"  R4+KK     = {v_4+v_kk:.10f}")
    print(f"  KK formula diff = {v_11 - v_4 - v_kk:.2e}")
