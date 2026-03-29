"""Quick float check: verify KK formula for (1,1)-string numerically using SymPy."""
import sys; sys.path.insert(0, 'src')
import sympy as sp
from sympy import Rational as R, sqrt, cancel, Float
from sugra import Metric

# Use a SIMPLE explicit H to keep expressions manageable
y0, y1 = sp.symbols('y0 y1', real=True)
t, x1 = sp.symbols('t x1', real=True)
z1s, z2s = sp.symbols('z1 z2', real=True)

# H = 1 + a/(y0²+y1²) — simpler than 1/r^6
a = sp.Symbol('a', positive=True)
rsq = y0**2 + y1**2
H = 1 + a/rsq

M0 = sp.Matrix([[sp.sqrt(H), 0], [0, 1/sp.sqrt(H)]])
Lambda = sp.Matrix([[1, 0], [1, 1]])
M = Lambda * M0 * Lambda.T
M_inv = M.inv()

coords = [t, x1, z1s, z2s, y0, y1]

# (1,1) 6d metric
g_11 = sp.zeros(6, 6)
g_11[0,0] = -1/H; g_11[1,1] = 1/H
g_11[2,2] = M[0,0]; g_11[2,3] = M[0,1]; g_11[3,2] = M[1,0]; g_11[3,3] = M[1,1]
g_11[4,4] = 1; g_11[5,5] = 1

# F1 6d metric
g_f1 = sp.zeros(6, 6)
g_f1[0,0] = -1/H; g_f1[1,1] = 1/H
g_f1[2,2] = M0[0,0]; g_f1[3,3] = M0[1,1]
g_f1[4,4] = 1; g_f1[5,5] = 1

# 4d base
g4 = sp.zeros(4, 4)
g4[0,0] = -1/H; g4[1,1] = 1/H; g4[2,2] = 1; g4[3,3] = 1

print("Computing Ricci tensors with H = 1 + a/(y0²+y1²)...")
print("(1,1) 6d...")
Ric_11 = Metric(g_11, coords).ricci_tensor(simplify_func=cancel)
print("F1 6d...")
Ric_f1 = Metric(g_f1, coords).ricci_tensor(simplify_func=cancel)
print("4d base...")
Ric_4 = Metric(g4, [t, x1, y0, y1]).ricci_tensor(simplify_func=cancel)

# KK term
dMinv_y0 = sp.Matrix([[sp.diff(M_inv[i,j], y0) for j in range(2)] for i in range(2)])
dMinv_y1 = sp.Matrix([[sp.diff(M_inv[i,j], y1) for j in range(2)] for i in range(2)])
dM_y0 = sp.Matrix([[sp.diff(M[i,j], y0) for j in range(2)] for i in range(2)])
dM_y1 = sp.Matrix([[sp.diff(M[i,j], y1) for j in range(2)] for i in range(2)])

kk_00 = R(1,4) * (dMinv_y0 * dM_y0).trace()
kk_01 = R(1,4) * (dMinv_y0 * dM_y1).trace()
kk_11 = R(1,4) * (dMinv_y1 * dM_y1).trace()

# Evaluate at specific point
pt = {y0: R(3,2), y1: R(1,2), a: R(1)}

print("\nEvaluating at y0=3/2, y1=1/2, a=1:")
print(f"  H = {H.subs(pt)}, r² = {rsq.subs(pt)}")

print("\n--- SL(2,Z) invariance: ℛ(1,1) vs ℛ(F1) ---")
for label, i, j in [("y0,y0", 4, 4), ("y0,y1", 4, 5), ("y1,y1", 5, 5),
                     ("t,t", 0, 0), ("z1,z1", 2, 2), ("z1,z2", 2, 3)]:
    v11 = float(Ric_11[i,j].subs(pt))
    vf1 = float(Ric_f1[i,j].subs(pt))
    diff = v11 - vf1
    status = "✓" if abs(diff) < 1e-10 else "✗"
    print(f"  [{label}]: (1,1)={v11:.8f}, F1={vf1:.8f}, diff={diff:.2e} {status}")

print("\n--- KK formula: ℛ₆(1,1) vs R₄ + KK ---")
for label, i6, j6, i4, j4, kk_val in [
    ("y0,y0", 4, 4, 2, 2, kk_00),
    ("y0,y1", 4, 5, 2, 3, kk_01),
    ("y1,y1", 5, 5, 3, 3, kk_11),
]:
    v6 = float(Ric_11[i6,j6].subs(pt))
    v4 = float(Ric_4[i4,j4].subs(pt))
    vkk = float(kk_val.subs(pt))
    diff = v6 - v4 - vkk
    status = "✓" if abs(diff) < 1e-10 else "✗"
    print(f"  [{label}]: ℛ₆={v6:.8f}, R₄={v4:.8f}, KK={vkk:.8f}, diff={diff:.2e} {status}")

print("\n--- F1 KK formula check ---")
dMinv_y0_f1 = sp.Matrix([[sp.diff((M0.inv())[i,j], y0) for j in range(2)] for i in range(2)])
dMinv_y1_f1 = sp.Matrix([[sp.diff((M0.inv())[i,j], y1) for j in range(2)] for i in range(2)])
dM_y0_f1 = sp.Matrix([[sp.diff(M0[i,j], y0) for j in range(2)] for i in range(2)])
dM_y1_f1 = sp.Matrix([[sp.diff(M0[i,j], y1) for j in range(2)] for i in range(2)])

kk_f1_00 = R(1,4) * (dMinv_y0_f1 * dM_y0_f1).trace()
kk_f1_01 = R(1,4) * (dMinv_y0_f1 * dM_y1_f1).trace()
kk_f1_11 = R(1,4) * (dMinv_y1_f1 * dM_y1_f1).trace()

for label, i6, j6, i4, j4, kk_val in [
    ("y0,y0", 4, 4, 2, 2, kk_f1_00),
    ("y0,y1", 4, 5, 2, 3, kk_f1_01),
    ("y1,y1", 5, 5, 3, 3, kk_f1_11),
]:
    v6 = float(Ric_f1[i6,j6].subs(pt))
    v4 = float(Ric_4[i4,j4].subs(pt))
    vkk = float(kk_val.subs(pt))
    diff = v6 - v4 - vkk
    status = "✓" if abs(diff) < 1e-10 else "✗"
    print(f"  [{label}]: ℛ₆(F1)={v6:.8f}, R₄+KK={v4+vkk:.8f}, diff={diff:.2e} {status}")

print("\n--- KK trace SL(2,Z) invariance ---")
for label, kk11, kkf1 in [("y0,y0", kk_00, kk_f1_00), ("y0,y1", kk_01, kk_f1_01), ("y1,y1", kk_11, kk_f1_11)]:
    v11 = float(kk11.subs(pt))
    vf1 = float(kkf1.subs(pt))
    print(f"  KK[{label}]: (1,1)={v11:.8f}, F1={vf1:.8f}, diff={v11-vf1:.2e}")
