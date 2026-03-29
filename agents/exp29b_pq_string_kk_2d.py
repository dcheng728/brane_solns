"""Experiment 29b: (1,1)-string KK decomposition — 2 transverse dimensions.

Same as exp29 but with 2 transverse coordinates instead of 8,
to avoid SymPy simplification issues with r = sqrt(y0²+...+y7²).

Tests the KK formula with non-diagonal M + form flux in a clean setting.
"""
import sys; sys.path.insert(0, 'src')
import sympy as sp
from sympy import Rational as R, cancel, simplify, sqrt, Function

print("=" * 70)
print("EXP 29b: (1,1)-string KK — 2 transverse coords")
print("=" * 70)

# ===================================================================
# Setup with 2 transverse coords
# ===================================================================
t, x1 = sp.symbols('t x1', real=True)
z1s, z2s = sp.symbols('z1 z2', real=True)
y0, y1 = sp.symbols('y0 y1', real=True)
r = sp.Symbol('r', positive=True)

coords_6 = [t, x1, z1s, z2s, y0, y1]
coords_4 = [t, x1, y0, y1]

# H = H(r), r = sqrt(y0² + y1²) — use explicit H(r) for clean derivatives
H = Function('H')(sqrt(y0**2 + y1**2))

# (1,1)-string torus
M0 = sp.Matrix([[H**R(1, 2), 0], [0, H**R(-1, 2)]])
Lambda = sp.Matrix([[1, 0], [1, 1]])
M = cancel(Lambda * M0 * Lambda.T)
M_inv = cancel(M.inv())
print(f"det(M) = {cancel(M.det())}")

# String-frame 10d (4d representative) metric
g_sf_4 = sp.zeros(4, 4)
g_sf_4[0, 0] = -1/H  # t
g_sf_4[1, 1] = 1/H   # x1
g_sf_4[2, 2] = 1      # y0
g_sf_4[3, 3] = 1      # y1

# 6d metric (4d + torus)
g_sf_6 = sp.zeros(6, 6)
g_sf_6[0, 0] = -1/H
g_sf_6[1, 1] = 1/H
g_sf_6[2, 2] = M[0, 0]
g_sf_6[2, 3] = M[0, 1]
g_sf_6[3, 2] = M[1, 0]
g_sf_6[3, 3] = M[1, 1]
g_sf_6[4, 4] = 1
g_sf_6[5, 5] = 1

# Also F1 6d metric for comparison
g_f1_6 = sp.zeros(6, 6)
g_f1_6[0, 0] = -1/H
g_f1_6[1, 1] = 1/H
g_f1_6[2, 2] = M0[0, 0]
g_f1_6[3, 3] = M0[1, 1]
g_f1_6[4, 4] = 1
g_f1_6[5, 5] = 1

from sugra import Metric

# ===================================================================
# Part A: Compute all Ricci tensors
# ===================================================================
print("\nComputing Ricci tensors...")
metric_6 = Metric(g_sf_6, coords_6)
Ric_6 = metric_6.ricci_tensor(simplify_func=cancel)
print("  (1,1) 6d done")

metric_4 = Metric(g_sf_4, coords_4)
Ric_4 = metric_4.ricci_tensor(simplify_func=cancel)
print("  10d (4d rep) done")

metric_f1 = Metric(g_f1_6, coords_6)
Ric_f1 = metric_f1.ricci_tensor(simplify_func=cancel)
print("  F1 6d done")

# ===================================================================
# Part B: KK formula for 10d directions
# ===================================================================
print("\n" + "=" * 70)
print("Part B: 10d KK formula: ℛ₆ = R₄ + (1/4)Tr(∂M⁻¹·∂M)")
print("=" * 70)

def dM(coord):
    return sp.Matrix([[cancel(sp.diff(M[i,j], coord)) for j in range(2)] for i in range(2)])

def dMinv(coord):
    return sp.Matrix([[cancel(sp.diff(M_inv[i,j], coord)) for j in range(2)] for i in range(2)])

# Map: 4d index → 6d index
idx_map = {0: 0, 1: 1, 2: 4, 3: 5}  # t→0, x1→1, y0→4, y1→5

all_10d_pass = True
for m in range(4):
    for n in range(m, 4):
        cm = coords_4[m]
        cn = coords_4[n]

        # KK term
        trace = cancel((dMinv(cm) * dM(cn)).trace())
        kk = cancel(R(1, 4) * trace)

        # 6d Ricci
        i6, j6 = idx_map[m], idx_map[n]
        ric6 = cancel(Ric_6[i6, j6])

        # 4d Ricci
        ric4 = cancel(Ric_4[m, n])

        rhs = cancel(ric4 + kk)
        diff = cancel(ric6 - rhs)

        if diff != 0:
            diff = simplify(diff)

        status = "✓" if diff == 0 else "✗"
        if status == "✗":
            all_10d_pass = False
            print(f"  [{cm},{cn}]: MISMATCH")
            print(f"    ℛ₆  = {ric6}")
            print(f"    R₄  = {ric4}")
            print(f"    KK  = {kk}")
            print(f"    diff = {diff}")
        elif ric6 != 0:
            print(f"  [{cm},{cn}]: {status}  (ℛ₆ = {ric6})")
        else:
            print(f"  [{cm},{cn}]: {status}  (= 0)")

if all_10d_pass:
    print("\n  ★ 10d KK formula VERIFIED ✓")
else:
    print("\n  ✗ Some components failed — checking numerically...")

    # Numerical check: substitute H = 1 + 1/r, evaluate at y0=1, y1=2
    from sympy import Subs
    Hnum = 1 + 1/sqrt(y0**2 + y1**2)
    num_subs = {H: Hnum, y0: sp.Rational(1), y1: sp.Rational(2)}

    print("\n  Numerical verification at y0=1, y1=2, H=1+1/r:")
    for m in range(4):
        for n in range(m, 4):
            cm = coords_4[m]
            cn = coords_4[n]
            trace = cancel((dMinv(cm) * dM(cn)).trace())
            kk = R(1, 4) * trace
            i6, j6 = idx_map[m], idx_map[n]
            ric6 = Ric_6[i6, j6]
            ric4 = Ric_4[m, n]

            # Evaluate numerically — need to handle H(r) substitution
            r_val = sqrt(sp.Rational(1)**2 + sp.Rational(2)**2)
            H_val = 1 + 1/r_val
            Hp_val = sp.diff(1 + 1/sqrt(y0**2+y1**2), y0).subs({y0: 1, y1: 2})
            # Just try direct numerical eval
            try:
                v6 = float(ric6.subs({y0: 1, y1: 2}).doit().rewrite(sp.Piecewise).simplify())
                v4 = float(ric4.subs({y0: 1, y1: 2}).doit().rewrite(sp.Piecewise).simplify())
                vkk = float(kk.subs({y0: 1, y1: 2}).doit().rewrite(sp.Piecewise).simplify())
                print(f"    [{cm},{cn}]: ℛ₆={v6:.6f}, R₄+KK={v4+vkk:.6f}, diff={v6-v4-vkk:.2e}")
            except Exception as e:
                print(f"    [{cm},{cn}]: numerical eval failed: {e}")

# ===================================================================
# Part C: SL(2,Z) invariance of 6d Ricci (10d directions)
# ===================================================================
print("\n" + "=" * 70)
print("Part C: ℛ^{SF}(1,1) vs ℛ^{SF}(F1) for 10d directions")
print("=" * 70)

ricci_match = True
for m in range(4):
    for n in range(m, 4):
        i6, j6 = idx_map[m], idx_map[n]
        v11 = cancel(Ric_6[i6, j6])
        vf1 = cancel(Ric_f1[i6, j6])
        diff = cancel(v11 - vf1)
        if diff != 0:
            diff = simplify(diff)
        if diff != 0:
            ricci_match = False
            print(f"  [{coords_4[m]},{coords_4[n]}]: DIFFER, diff = {diff}")
        elif v11 != 0:
            print(f"  [{coords_4[m]},{coords_4[n]}]: match ✓")

if ricci_match:
    print("\n  ★ ℛ^{SF}(1,1) = ℛ^{SF}(F1) for all 10d directions ✓")

# ===================================================================
# Part D: Torus KK formula
# ===================================================================
print("\n" + "=" * 70)
print("Part D: Torus KK formula")
print("=" * 70)

sqrt_g_4 = 1/H  # √|det g^S_4d|
transverse = [(y0, 4), (y1, 5)]  # (coord, 6d index)

div_MinvdM = sp.zeros(2, 2)
for c in range(2):
    for b in range(2):
        total = sp.Integer(0)
        for (yk, _) in transverse:
            MinvdM_k = cancel(M_inv * dM(yk))
            flux = cancel(sqrt_g_4 * 1 * MinvdM_k[c, b])  # g^{yk}=1
            total += sp.diff(flux, yk)
        div_MinvdM[c, b] = cancel(total / sqrt_g_4)

all_torus_pass = True
labels = ['z1', 'z2']
for a in range(2):
    for b in range(a, 2):
        ric_ab = cancel(Ric_6[2+a, 2+b])
        rhs = sp.Integer(0)
        for c in range(2):
            rhs += M[a, c] * div_MinvdM[c, b]
        rhs = cancel(R(-1, 2) * rhs)
        diff = cancel(ric_ab - rhs)
        if diff != 0:
            diff = simplify(diff)
        if diff == 0:
            print(f"  ℛ[{labels[a]},{labels[b]}] = {ric_ab}  ✓")
        else:
            print(f"  ℛ[{labels[a]},{labels[b]}]: MISMATCH! diff={diff}")
            all_torus_pass = False

if all_torus_pass:
    print("\n  ★ Torus KK formula VERIFIED ✓")

# ===================================================================
# Part E: Cross terms
# ===================================================================
print("\n" + "=" * 70)
print("Part E: Cross terms")
print("=" * 70)

cross_pass = True
for m_6d in [0, 1, 4, 5]:
    for a in range(2):
        val = cancel(Ric_6[m_6d, 2+a])
        if val != 0:
            val = simplify(val)
            if val != 0:
                cross_pass = False
                nm = {0:'t', 1:'x1', 4:'y0', 5:'y1'}[m_6d]
                print(f"  ℛ[{nm},{labels[a]}] = {val}  ✗")
if cross_pass:
    print("  ALL zero ✓")

# ===================================================================
# Part F: Torus covariance ℛ_{ab}(1,1) = Λ ℛ_{cd}(F1) Λ^T
# ===================================================================
print("\n" + "=" * 70)
print("Part F: Torus SL(2,Z) covariance")
print("=" * 70)

R_torus_11 = sp.Matrix([[cancel(Ric_6[2,2]), cancel(Ric_6[2,3])],
                         [cancel(Ric_6[3,2]), cancel(Ric_6[3,3])]])
R_torus_f1 = sp.Matrix([[cancel(Ric_f1[2,2]), cancel(Ric_f1[2,3])],
                         [cancel(Ric_f1[3,2]), cancel(Ric_f1[3,3])]])

expected = cancel(Lambda * R_torus_f1 * Lambda.T)
torus_cov = True
for i in range(2):
    for j in range(2):
        diff = cancel(R_torus_11[i,j] - expected[i,j])
        if diff != 0:
            diff = simplify(diff)
        if diff != 0:
            torus_cov = False
            print(f"  [{labels[i]},{labels[j]}]: MISMATCH! diff={diff}")

if torus_cov:
    print("  ★ ℛ_{ab}(1,1) = Λ·ℛ(F1)·Λ^T ✓")

# ===================================================================
# CONCLUSIONS
# ===================================================================
print("\n" + "=" * 70)
print("CONCLUSIONS")
print("=" * 70)
all_pass = all_10d_pass and all_torus_pass and cross_pass and ricci_match and torus_cov
print(f"""
  1. 10d KK:     ℛ = R + (1/4)Tr(∂M⁻¹·∂M)      {'✓' if all_10d_pass else '✗'}
  2. Torus KK:   ℛ = -(1/2)M·div(M⁻¹∂M)         {'✓' if all_torus_pass else '✗'}
  3. Cross:      ℛ_{ma} = 0                        {'✓' if cross_pass else '✗'}
  4. SL(2,Z):    ℛ(1,1) = ℛ(F1) (10d)             {'✓' if ricci_match else '✗'}
  5. Torus cov:  ℛ → ΛℛΛ^T                        {'✓' if torus_cov else '✗'}
""")
