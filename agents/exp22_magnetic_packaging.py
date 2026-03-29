"""Experiment 22: Magnetic packaging theorem verification.

Exp21 verified the G₄ packaging theorem for ELECTRIC branes:
  G₄ = F₃ ∧ dz_a  →  (G₄²)_{mn}/3! = (c^T M⁻¹ c)·(F₃²)_{mn}/2!

For magnetic branes (NS5, D5), G₄ = (*₄F₃) ∧ dz_a has TRANSVERSE legs
(y_i, y_j, y_k, z_a) instead of worldvolume legs (t, x_i, z_a, y_k).

Question: Does the same packaging identity hold for magnetic G₄?
  (G₄^{mag})²_{mn}/3! = (c^T M⁻¹ c)·((*₄F₃)²)_{mn}/2!  ?

The Hodge dual introduces √g factors which could break the simple
factorization.

Plan:
  A. NS5: Build G₄ = H₃^{mag} ∧ dz₂, verify packaging against 10d H₃^{mag}
  B. D5:  Build G₄ = F₃^{mag} ∧ dz₁, verify packaging
  C. S-dual magnetic pair: verify D5 = S(NS5) at packaging level
  D. Summary and comparison with electric sector
"""
import sys; sys.path.insert(0, 'src')
import sympy as sp
from sympy import Rational as R, ln, cancel
from sugra import (HarmonicFunction, Metric,
                   FormField, form_contraction, form_norm_squared)
from itertools import combinations
from sympy.combinatorics import Permutation

print("=" * 70)
print("EXP 22: Magnetic packaging theorem")
print("=" * 70)

# ===================================================================
# Part A: NS5 magnetic packaging
# ===================================================================
print("\nPART A: NS5 — magnetic G₄ packaging")
print("-" * 50)

# NS5: 6d worldvolume + 4d transverse
wv_coords = list(sp.symbols('t x1 x2 x3 x4 x5', real=True))
z1s, z2s = sp.symbols('z1 z2', real=True)
trans_coords = list(sp.symbols('y0:4', real=True))
coords = wv_coords + [z1s, z2s] + trans_coords
D = 12

hf = HarmonicFunction(transverse_coords=trans_coords)
H = sp.Function('H')(hf.r_expr)

# --- NS5 12d Einstein-frame metric ---
g_ef_ns5 = sp.zeros(D, D)
g_ef_ns5[0, 0] = -H**R(-1, 4)   # t
for k in range(1, 6):
    g_ef_ns5[k, k] = H**R(-1, 4) # x1,...,x5
g_ef_ns5[6, 6] = H**R(-1, 2)    # z1 (= e^{-Φ} = H^{-1/2})
g_ef_ns5[7, 7] = H**R(1, 2)     # z2 (= e^Φ = H^{1/2})
for k in range(4):
    g_ef_ns5[8 + k, 8 + k] = H**R(3, 4)  # transverse

M_ns5 = sp.Matrix([[H**R(-1, 2), 0], [0, H**R(1, 2)]])
M_ns5_inv = sp.Matrix([[H**R(1, 2), 0], [0, H**R(-1, 2)]])

metric_ns5 = Metric(g_ef_ns5, coords)

# --- Build magnetic H₃ = *₄^{flat} dH ---
def levi_civita_4(i, j, k, l):
    perm = [i, j, k, l]
    if len(set(perm)) < 4:
        return 0
    sign = 1
    for a in range(4):
        for b in range(a+1, 4):
            if perm[a] > perm[b]:
                sign *= -1
    return sign

# G₄(NS5) = H₃^{mag} ∧ dz₂
# Components: G₄_{y_i, y_j, y_k, z₂} = ε_{ijkl} ∂_{y_l} H
G4_ns5 = FormField(rank=4, dim=D)
for i, j, k in combinations(range(4), 3):
    l = [x for x in range(4) if x not in (i, j, k)][0]
    eps = levi_civita_4(i, j, k, l)
    if eps == 0:
        continue
    y_l = trans_coords[l]
    dH = sp.diff(H, y_l)
    # 12d indices: y_i=8+i, z₂=7
    original = [8+i, 8+j, 8+k, 7]
    sorted_idx = sorted(original)
    perm_map = [sorted_idx.index(x) for x in original]
    sort_sign = Permutation(perm_map).signature()
    G4_ns5[tuple(sorted_idx)] = sp.Rational(eps) * sort_sign * dH

print("G₄(NS5) sample components:")
for idx in sorted(G4_ns5.nonzero_components.keys())[:4]:
    names = [coords[i].name if hasattr(coords[i], 'name') else str(coords[i]) for i in idx]
    print(f"  G₄[{','.join(names)}] = {cancel(G4_ns5[idx])}")

# --- 12d form contractions ---
print("\nComputing 12d form contraction for NS5...")
FF_ns5_12d = form_contraction(G4_ns5, metric_ns5)
S_ns5_12d = form_norm_squared(G4_ns5, metric_ns5)
S_ns5_12d_val = hf.substitute(cancel(S_ns5_12d))
print(f"|G₄|²(NS5, 12d) = {S_ns5_12d_val}")

# --- 10d magnetic H₃ and its contraction ---
# 10d NS5 Einstein-frame metric
coords_10 = wv_coords + trans_coords
g_10d_ns5 = sp.zeros(10, 10)
g_10d_ns5[0, 0] = -H**R(-1, 4)
for k in range(1, 6):
    g_10d_ns5[k, k] = H**R(-1, 4)
for k in range(4):
    g_10d_ns5[6 + k, 6 + k] = H**R(3, 4)

metric_10d_ns5 = Metric(g_10d_ns5, coords_10)

H3_mag = FormField(rank=3, dim=10)
for i, j, k in combinations(range(4), 3):
    l = [x for x in range(4) if x not in (i, j, k)][0]
    eps = levi_civita_4(i, j, k, l)
    if eps == 0:
        continue
    y_l = trans_coords[l]
    dH = sp.diff(H, y_l)
    # 10d indices: y_i = 6+i
    original = [6+i, 6+j, 6+k]
    sorted_idx = sorted(original)
    perm_map = [sorted_idx.index(x) for x in original]
    sort_sign = Permutation(perm_map).signature()
    H3_mag[tuple(sorted_idx)] = sp.Rational(eps) * sort_sign * dH

FF_H3_10d = form_contraction(H3_mag, metric_10d_ns5)
S_H3_10d = form_norm_squared(H3_mag, metric_10d_ns5)
S_H3_10d_val = hf.substitute(cancel(S_H3_10d))
print(f"|H₃^{{mag}}|²(10d) = {S_H3_10d_val}")

# --- Packaging test ---
# c = (0, 1) for NS5 (NSNS), c^T M⁻¹ c = (M⁻¹)^{22} = H^{-1/2}
c_ns5 = sp.Matrix([0, 1])
cTMinvc_ns5 = cancel((c_ns5.T * M_ns5_inv * c_ns5)[0, 0])
cTMinvc_ns5_val = hf.substitute(cancel(cTMinvc_ns5))
print(f"\nc^T M⁻¹ c (NS5) = {cTMinvc_ns5_val}")

print("\nPackaging test: (G₄²)_{mn}/3! = (c^T M⁻¹ c)·(H₃^{mag})²_{mn}/2!")
pack_ns5 = True
test_blocks = [
    (0, 0, 't'), (1, 1, 'x1'),      # worldvolume
    (8, 6, 'y0'), (9, 7, 'y1'),      # transverse
]

for i_12, i_10, name in test_blocks:
    ff_12 = hf.substitute(cancel(FF_ns5_12d[i_12, i_12]))
    ff_10 = hf.substitute(cancel(FF_H3_10d[i_10, i_10]))
    lhs = cancel(ff_12 / 6)   # G₄²/3!
    rhs = hf.substitute(cancel(cTMinvc_ns5 * ff_10 / 2))  # c^T M⁻¹ c · H₃²/2!
    diff = cancel(lhs - rhs)
    status = '✓' if diff == 0 else f'✗ diff={diff}'
    if diff != 0:
        pack_ns5 = False
    print(f"  [{name}]: G₄²/3!={lhs}, pred={rhs}, {status}")

# Norm check
S_pred_ns5 = hf.substitute(cancel(cTMinvc_ns5 * S_H3_10d))
norm_diff_ns5 = cancel(S_ns5_12d_val - S_pred_ns5)
print(f"\n  |G₄|² = {S_ns5_12d_val}")
print(f"  c^T M⁻¹ c · |H₃^{{mag}}|² = {S_pred_ns5}")
print(f"  Norm match: {'✓' if norm_diff_ns5 == 0 else '✗ diff=' + str(norm_diff_ns5)}")

print(f"\n  NS5 Packaging: {'★ VERIFIED ★' if pack_ns5 and norm_diff_ns5 == 0 else 'FAILED'}")

# --- Torus-direction packaging ---
print("\nTorus packaging (NS5):")
# FF[z2,z2]: G₄ has z₂ leg, so this should be nonzero
# FF[z1,z1]: no z₁ leg in G₄, should be zero
ff_z1z1_ns5 = hf.substitute(cancel(FF_ns5_12d[6, 6]))
ff_z2z2_ns5 = hf.substitute(cancel(FF_ns5_12d[7, 7]))
ff_z1z2_ns5 = hf.substitute(cancel(FF_ns5_12d[6, 7]))
print(f"  FF[z1,z1] = {ff_z1z1_ns5}  (expected 0)")
print(f"  FF[z2,z2] = {ff_z2z2_ns5}")
print(f"  FF[z1,z2] = {ff_z1z2_ns5}  (expected 0)")

# Torus packaging: FF_{za,zb} = c_a·(M⁻¹c)_b × normalization × |H₃|²
# For NS5 c=(0,1): FF_{z2,z2} = 1·(M⁻¹c)₂·... = (M⁻¹)^{22}·...
# FF_{z1,*} = 0 since c₁=0.
torus_z1_ok = ff_z1z1_ns5 == 0 and ff_z1z2_ns5 == 0
print(f"\n  z1 components zero (c₁=0): {'✓' if torus_z1_ok else '✗'}")

# ===================================================================
# Part B: D5 magnetic packaging
# ===================================================================
print("\n" + "=" * 70)
print("PART B: D5 — magnetic G₄ packaging")
print("-" * 50)

# D5: same worldvolume/transverse structure as NS5 but z₁↔z₂
# M = diag(H^{1/2}, H^{-1/2}) — S-dual of NS5
# G₄ = F₃^{mag} ∧ dz₁

g_ef_d5 = sp.zeros(D, D)
g_ef_d5[0, 0] = -H**R(-1, 4)
for k in range(1, 6):
    g_ef_d5[k, k] = H**R(-1, 4)
g_ef_d5[6, 6] = H**R(1, 2)     # z1 (= e^Φ = H^{1/2})
g_ef_d5[7, 7] = H**R(-1, 2)    # z2 (= e^{-Φ} = H^{-1/2})
for k in range(4):
    g_ef_d5[8 + k, 8 + k] = H**R(3, 4)

M_d5 = sp.Matrix([[H**R(1, 2), 0], [0, H**R(-1, 2)]])
M_d5_inv = sp.Matrix([[H**R(-1, 2), 0], [0, H**R(1, 2)]])

metric_d5 = Metric(g_ef_d5, coords)

# G₄(D5) = F₃^{mag} ∧ dz₁  (same magnetic 3-form, but ∧ dz₁)
G4_d5 = FormField(rank=4, dim=D)
for i, j, k in combinations(range(4), 3):
    l = [x for x in range(4) if x not in (i, j, k)][0]
    eps = levi_civita_4(i, j, k, l)
    if eps == 0:
        continue
    y_l = trans_coords[l]
    dH = sp.diff(H, y_l)
    # 12d indices: y_i=8+i, z₁=6
    original = [8+i, 8+j, 8+k, 6]
    sorted_idx = sorted(original)
    perm_map = [sorted_idx.index(x) for x in original]
    sort_sign = Permutation(perm_map).signature()
    G4_d5[tuple(sorted_idx)] = sp.Rational(eps) * sort_sign * dH

print("G₄(D5) sample components:")
for idx in sorted(G4_d5.nonzero_components.keys())[:4]:
    names = [coords[i].name if hasattr(coords[i], 'name') else str(coords[i]) for i in idx]
    print(f"  G₄[{','.join(names)}] = {cancel(G4_d5[idx])}")

# --- 12d form contractions ---
print("\nComputing 12d form contraction for D5...")
FF_d5_12d = form_contraction(G4_d5, metric_d5)
S_d5_12d = form_norm_squared(G4_d5, metric_d5)
S_d5_12d_val = hf.substitute(cancel(S_d5_12d))
print(f"|G₄|²(D5, 12d) = {S_d5_12d_val}")

# --- D5 10d metric ---
g_10d_d5 = sp.zeros(10, 10)
g_10d_d5[0, 0] = -H**R(-1, 4)
for k in range(1, 6):
    g_10d_d5[k, k] = H**R(-1, 4)
for k in range(4):
    g_10d_d5[6 + k, 6 + k] = H**R(3, 4)

metric_10d_d5 = Metric(g_10d_d5, coords_10)

# F₃^{mag} in 10d — same form structure as H₃^{mag}
F3_mag = FormField(rank=3, dim=10)
for i, j, k in combinations(range(4), 3):
    l = [x for x in range(4) if x not in (i, j, k)][0]
    eps = levi_civita_4(i, j, k, l)
    if eps == 0:
        continue
    y_l = trans_coords[l]
    dH = sp.diff(H, y_l)
    original = [6+i, 6+j, 6+k]
    sorted_idx = sorted(original)
    perm_map = [sorted_idx.index(x) for x in original]
    sort_sign = Permutation(perm_map).signature()
    F3_mag[tuple(sorted_idx)] = sp.Rational(eps) * sort_sign * dH

FF_F3_10d = form_contraction(F3_mag, metric_10d_d5)
S_F3_10d = form_norm_squared(F3_mag, metric_10d_d5)
S_F3_10d_val = hf.substitute(cancel(S_F3_10d))
print(f"|F₃^{{mag}}|²(10d) = {S_F3_10d_val}")

# c = (1, 0) for D5 (RR), c^T M⁻¹ c = (M⁻¹)^{11} = H^{-1/2}
c_d5 = sp.Matrix([1, 0])
cTMinvc_d5 = cancel((c_d5.T * M_d5_inv * c_d5)[0, 0])
cTMinvc_d5_val = hf.substitute(cancel(cTMinvc_d5))
print(f"\nc^T M⁻¹ c (D5) = {cTMinvc_d5_val}")

print("\nPackaging test: (G₄²)_{mn}/3! = (c^T M⁻¹ c)·(F₃^{mag})²_{mn}/2!")
pack_d5 = True
for i_12, i_10, name in test_blocks:
    ff_12 = hf.substitute(cancel(FF_d5_12d[i_12, i_12]))
    ff_10 = hf.substitute(cancel(FF_F3_10d[i_10, i_10]))
    lhs = cancel(ff_12 / 6)
    rhs = hf.substitute(cancel(cTMinvc_d5 * ff_10 / 2))
    diff = cancel(lhs - rhs)
    status = '✓' if diff == 0 else f'✗ diff={diff}'
    if diff != 0:
        pack_d5 = False
    print(f"  [{name}]: G₄²/3!={lhs}, pred={rhs}, {status}")

# Norm check
S_pred_d5 = hf.substitute(cancel(cTMinvc_d5 * S_F3_10d))
norm_diff_d5 = cancel(S_d5_12d_val - S_pred_d5)
print(f"\n  |G₄|² = {S_d5_12d_val}")
print(f"  c^T M⁻¹ c · |F₃^{{mag}}|² = {S_pred_d5}")
print(f"  Norm match: {'✓' if norm_diff_d5 == 0 else '✗ diff=' + str(norm_diff_d5)}")

print(f"\n  D5 Packaging: {'★ VERIFIED ★' if pack_d5 and norm_diff_d5 == 0 else 'FAILED'}")

# --- Torus-direction packaging ---
print("\nTorus packaging (D5):")
ff_z1z1_d5 = hf.substitute(cancel(FF_d5_12d[6, 6]))
ff_z2z2_d5 = hf.substitute(cancel(FF_d5_12d[7, 7]))
ff_z1z2_d5 = hf.substitute(cancel(FF_d5_12d[6, 7]))
print(f"  FF[z1,z1] = {ff_z1z1_d5}")
print(f"  FF[z2,z2] = {ff_z2z2_d5}  (expected 0)")
print(f"  FF[z1,z2] = {ff_z1z2_d5}  (expected 0)")

torus_z2_ok = ff_z2z2_d5 == 0 and ff_z1z2_d5 == 0
print(f"\n  z2 components zero (c₂=0): {'✓' if torus_z2_ok else '✗'}")

# ===================================================================
# Part C: S-duality of magnetic packaging
# ===================================================================
print("\n" + "=" * 70)
print("PART C: S-duality — NS5 ↔ D5 comparison")
print("-" * 50)

# Under S-duality: z₁↔z₂, Φ→-Φ
# NS5 c^T M⁻¹ c = H^{-1/2} (e^{-Φ}), D5 c^T M⁻¹ c = H^{-1/2} (e^{+Φ})
# But H^{-1/2} is the SAME since e^{-Φ}(NS5) = H^{-1/2} and e^{Φ}(D5) = H^{-1/2}
# (D5 has Φ = -(1/2)lnH, so e^Φ = H^{-1/2})

print(f"c^T M⁻¹ c comparison:")
print(f"  NS5: {cTMinvc_ns5_val}")
print(f"  D5:  {cTMinvc_d5_val}")
s_dual_coupling = cancel(cTMinvc_ns5_val - cTMinvc_d5_val) == 0
print(f"  Equal: {'✓' if s_dual_coupling else '✗'}")

# The 10d metrics are identical (NS5 and D5 have same warp factors)
# so |H₃^{mag}|²(NS5) = |F₃^{mag}|²(D5) with same metric
print(f"\n|3-form|² comparison (10d):")
print(f"  NS5 |H₃^{{mag}}|² = {S_H3_10d_val}")
print(f"  D5  |F₃^{{mag}}|² = {S_F3_10d_val}")
s_dual_norm = cancel(S_H3_10d_val - S_F3_10d_val) == 0
print(f"  Equal: {'✓' if s_dual_norm else '✗'}")

# The 12d norms should therefore also be equal
print(f"\n12d |G₄|² comparison:")
print(f"  NS5: {S_ns5_12d_val}")
print(f"  D5:  {S_d5_12d_val}")
s_dual_g4 = cancel(S_ns5_12d_val - S_d5_12d_val) == 0
print(f"  Equal: {'✓' if s_dual_g4 else '✗'}")

# ===================================================================
# Part D: Electric vs Magnetic comparison
# ===================================================================
print("\n" + "=" * 70)
print("PART D: Electric vs magnetic packaging — structural comparison")
print("-" * 50)

print("""
PACKAGING THEOREM STRUCTURE:

For ANY brane with G₄ = c_a F₃^a ∧ dz_a:

  (G₄²)_{mn}/3! = (c^T M⁻¹ c) · (F₃²)_{mn}/2!
  |G₄|² = (c^T M⁻¹ c) · |F₃|²

ELECTRIC sector (F1, D1):
  F₃^a has worldvolume legs (t, x1, y_k)
  G₄ = F₃ ∧ dz_a has legs (t, x1, y_k, z_a)
  FF_{mn} ≠ 0 for worldvolume AND transverse (m ∈ {t,x1} or {y_k})
  |G₄|² < 0 (timelike leg)

MAGNETIC sector (NS5, D5):
  F₃^a has transverse legs (y_i, y_j, y_k)
  G₄ = (*₄F₃) ∧ dz_a has legs (y_i, y_j, y_k, z_a)
  FF_{mn} = 0 for worldvolume (m ∈ {t,x1,...,x5})
  FF_{mn} ≠ 0 for transverse (m ∈ {y_k})
  |G₄|² > 0 (all spacelike legs)

The packaging identity is PURELY KINEMATIC — it follows from the
factored structure G₄ = F₃ ∧ dz and the block-diagonal metric
g = g^{10d} ⊕ M. The electric/magnetic distinction only affects
which 10d components are nonzero, not the validity of the identity.
""")

# ===================================================================
# Part E: Verification summary
# ===================================================================
print("=" * 70)
print("VERIFICATION SUMMARY")
print("=" * 70)

all_pass = pack_ns5 and (norm_diff_ns5 == 0) and pack_d5 and (norm_diff_d5 == 0)
print(f"""
MAGNETIC PACKAGING THEOREM:
  NS5 (NSNS magnetic, c=(0,1)): {'★ VERIFIED ★' if pack_ns5 and norm_diff_ns5 == 0 else 'FAILED'}
  D5  (RR magnetic, c=(1,0)):   {'★ VERIFIED ★' if pack_d5 and norm_diff_d5 == 0 else 'FAILED'}

S-DUALITY CONSISTENCY:
  c^T M⁻¹ c equal:  {'✓' if s_dual_coupling else '✗'}
  |F₃|² equal:      {'✓' if s_dual_norm else '✗'}
  |G₄|² equal:      {'✓' if s_dual_g4 else '✗'}

TORUS COMPONENTS:
  NS5 FF[z1]=0 (c₁=0): {'✓' if torus_z1_ok else '✗'}
  D5  FF[z2]=0 (c₂=0): {'✓' if torus_z2_ok else '✗'}

COMPLETE PACKAGING THEOREM (electric + magnetic):
  F1  (NSNS electric):  exp21 ✓
  D1  (RR electric):    exp21 ✓
  (1,1) (non-diagonal): exp21 ✓
  NS5 (NSNS magnetic):  {'✓' if pack_ns5 else '✗'}
  D5  (RR magnetic):    {'✓' if pack_d5 else '✗'}
  D3  (self-dual F₅):   exp14 (separate sector)

{'★★★ ALL BRANES PASS — PACKAGING THEOREM FULLY VERIFIED ★★★' if all_pass else 'SOME TESTS FAILED'}
""")
