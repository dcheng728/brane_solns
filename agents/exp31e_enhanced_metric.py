"""Experiment 31e: Verify the UNIVERSAL 12d formula with enhanced metric.

FORMULA:
  ℛ_{MN} = (1/2)[FF_{MN}/3! - (1/4)|G₄|² · ĝ_{MN}]

where ĝ_{MN} = g_{MN} + M̃_{MN}  ("enhanced metric")
  g_{MN} = full 12d metric
  M̃_{MN} = torus metric embedded in 12d (0 for non-torus indices)

Effect:
  10d directions: ĝ = g → standard T with λ=1/4
  Torus directions: ĝ = 2g → effective λ=1/2

This is equivalent to: ℛ = (1/2)[FF/3! - (1/4)|G₄|²g] - (1/8)|G₄|²M̃
"""
import sys; sys.path.insert(0, 'src')
import sympy as sp
from sympy import Rational as R, cancel
from sympy.combinatorics import Permutation
from sugra import (HarmonicFunction, warped_product, Metric,
                   FormField, exterior_derivative,
                   form_contraction, form_norm_squared)
from itertools import combinations


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


def verify_enhanced_formula(name, Ric, metric, FF, norm_val, M, hf,
                            torus_indices, labels):
    """ℛ = (1/2)[FF/3! - (1/4)|G₄|²·ĝ] where ĝ = g + M̃."""
    print(f"\n{name}: ℛ = (1/2)[FF/3! - (1/4)|G₄|²·ĝ]")
    print("-"*60)
    all_pass = True
    for i12, label in labels:
        ric_val = hf.substitute(cancel(Ric[i12, i12]))
        ff_val = hf.substitute(cancel(FF[i12, i12]))
        g_val = hf.substitute(cancel(metric.matrix[i12, i12]))

        # Enhanced metric: ĝ = g + M̃
        if i12 in torus_indices:
            a_idx = torus_indices.index(i12)
            m_tilde = hf.substitute(cancel(M[a_idx, a_idx]))
        else:
            m_tilde = sp.Integer(0)
        g_hat = cancel(g_val + m_tilde)

        T = cancel(R(1,2) * (ff_val/6 - R(1,4)*norm_val*g_hat))
        diff = cancel(ric_val - T)
        status = "✓" if diff == 0 else f"✗ residual={diff}"
        if diff != 0:
            all_pass = False
        print(f"  [{label}]: ĝ={g_hat}, diff={diff}  {status}")

    # Off-diagonal torus
    i_z1, i_z2 = torus_indices
    ric_cross = hf.substitute(cancel(Ric[i_z1, i_z2]))
    ff_cross = hf.substitute(cancel(FF[i_z1, i_z2]))
    g_cross = hf.substitute(cancel(metric.matrix[i_z1, i_z2]))
    m_cross = hf.substitute(cancel(M[0, 1]))
    g_hat_cross = cancel(g_cross + m_cross)
    T_cross = cancel(R(1,2)*(ff_cross/6 - R(1,4)*norm_val*g_hat_cross))
    diff_cross = cancel(ric_cross - T_cross)
    status = "✓" if diff_cross == 0 else f"✗ residual={diff_cross}"
    if diff_cross != 0:
        all_pass = False
    print(f"  [z1z2]: diff={diff_cross}  {status}")
    print(f"  {'ALL PASS ✓' if all_pass else 'SOME FAIL ✗'}")
    return all_pass


z1s, z2s = sp.symbols('z1 z2', real=True)


# ===================================================================
# F1-STRING
# ===================================================================
print("="*70)
print("F1-STRING (NSNS electric)")
print("="*70)

wv_f1 = list(sp.symbols('t x1', real=True))
trans_f1 = list(sp.symbols('y0:8', real=True))
coords_f1 = wv_f1 + [z1s, z2s] + trans_f1

hf_f1 = HarmonicFunction(transverse_coords=trans_f1)
H_f1 = sp.Function('H')(hf_f1.r_expr)

m_f1 = warped_product(
    [H_f1**R(-3,4), H_f1**R(1,2), H_f1**R(-1,2), H_f1**R(1,4)],
    [2, 1, 1, 8], ['lorentzian','euclidean','euclidean','euclidean'],
    coords_f1)
M_f1 = sp.Matrix([[H_f1**R(1,2), 0], [0, H_f1**R(-1,2)]])

C3 = FormField(rank=3, dim=12)
C3[(0, 1, 3)] = 1/H_f1
G4_f1 = exterior_derivative(C3, coords_f1)

print("Computing Ricci and form data...")
Ric_f1 = m_f1.ricci_tensor(simplify_func=cancel)
FF_f1 = form_contraction(G4_f1, m_f1)
norm_f1 = hf_f1.substitute(cancel(form_norm_squared(G4_f1, m_f1)))
print(f"|G₄|²(F1) = {norm_f1}")

verify_enhanced_formula("F1", Ric_f1, m_f1, FF_f1, norm_f1, M_f1, hf_f1,
                        [2, 3], [(0,'t'),(1,'x1'),(2,'z1'),(3,'z2'),(4,'y0'),(5,'y1')])


# ===================================================================
# D1-STRING (RR electric, S-dual of F1)
# ===================================================================
print("\n\n" + "="*70)
print("D1-STRING (RR electric)")
print("="*70)

# D1: z1↔z2 swapped vs F1 → H^{-3/4}ds²_{1,1} + H^{-1/2}dz₁² + H^{1/2}dz₂² + H^{1/4}ds²_8
m_d1 = warped_product(
    [H_f1**R(-3,4), H_f1**R(-1,2), H_f1**R(1,2), H_f1**R(1,4)],
    [2, 1, 1, 8], ['lorentzian','euclidean','euclidean','euclidean'],
    coords_f1)
M_d1 = sp.Matrix([[H_f1**R(-1,2), 0], [0, H_f1**R(1,2)]])

# G₄ = F₃ ∧ dz₁: C_{t,x1,z1} = H^{-1}
C3_d1 = FormField(rank=3, dim=12)
C3_d1[(0, 1, 2)] = 1/H_f1  # t=0, x1=1, z1=2
G4_d1 = exterior_derivative(C3_d1, coords_f1)

print("Computing Ricci and form data...")
Ric_d1 = m_d1.ricci_tensor(simplify_func=cancel)
FF_d1 = form_contraction(G4_d1, m_d1)
norm_d1 = hf_f1.substitute(cancel(form_norm_squared(G4_d1, m_d1)))
print(f"|G₄|²(D1) = {norm_d1}")

verify_enhanced_formula("D1", Ric_d1, m_d1, FF_d1, norm_d1, M_d1, hf_f1,
                        [2, 3], [(0,'t'),(1,'x1'),(2,'z1'),(3,'z2'),(4,'y0'),(5,'y1')])


# ===================================================================
# NS5-BRANE (NSNS magnetic)
# ===================================================================
print("\n\n" + "="*70)
print("NS5-BRANE (NSNS magnetic)")
print("="*70)

wv_ns5 = list(sp.symbols('t x1 x2 x3 x4 x5', real=True))
trans_ns5 = list(sp.symbols('y0:4', real=True))
coords_ns5 = wv_ns5 + [z1s, z2s] + trans_ns5

hf_ns5 = HarmonicFunction(transverse_coords=trans_ns5)
H_ns5 = sp.Function('H')(hf_ns5.r_expr)

m_ns5 = warped_product(
    [H_ns5**R(-1,4), H_ns5**R(-1,2), H_ns5**R(1,2), H_ns5**R(3,4)],
    [6, 1, 1, 4], ['lorentzian','euclidean','euclidean','euclidean'],
    coords_ns5)
M_ns5 = sp.Matrix([[H_ns5**R(-1,2), 0], [0, H_ns5**R(1,2)]])

# G₄ = H₃^{flat} ∧ dz₂ with H₃_{ijk} = ε_{ijkl} ∂_{yl}H
G4_ns5 = FormField(rank=4, dim=12)
y_base_ns5 = 8
z2_idx_ns5 = 7

for i, j, k in combinations(range(4), 3):
    l = [x for x in range(4) if x not in (i, j, k)][0]
    eps = levi_civita_4(i, j, k, l)
    if eps == 0:
        continue
    dH = sp.diff(H_ns5, trans_ns5[l])
    original = [y_base_ns5+i, y_base_ns5+j, y_base_ns5+k, z2_idx_ns5]
    sorted_idx = tuple(sorted(original))
    perm_map = [sorted_idx.index(x) for x in original]
    sort_sign = Permutation(perm_map).signature()
    G4_ns5[sorted_idx] = sp.Rational(eps) * sort_sign * dH

print("Computing Ricci and form data...")
Ric_ns5 = m_ns5.ricci_tensor(simplify_func=cancel)
FF_ns5 = form_contraction(G4_ns5, m_ns5)
norm_ns5 = hf_ns5.substitute(cancel(form_norm_squared(G4_ns5, m_ns5)))
print(f"|G₄|²(NS5) = {norm_ns5}")

verify_enhanced_formula("NS5", Ric_ns5, m_ns5, FF_ns5, norm_ns5, M_ns5, hf_ns5,
                        [6, 7], [(0,'t'),(1,'x1'),(6,'z1'),(7,'z2'),(8,'y0'),(9,'y1')])


# ===================================================================
# D5-BRANE (RR magnetic, S-dual of NS5)
# ===================================================================
print("\n\n" + "="*70)
print("D5-BRANE (RR magnetic)")
print("="*70)

# D5: z1↔z2 swapped vs NS5
m_d5 = warped_product(
    [H_ns5**R(-1,4), H_ns5**R(1,2), H_ns5**R(-1,2), H_ns5**R(3,4)],
    [6, 1, 1, 4], ['lorentzian','euclidean','euclidean','euclidean'],
    coords_ns5)
M_d5 = sp.Matrix([[H_ns5**R(1,2), 0], [0, H_ns5**R(-1,2)]])

# G₄ = F₃^{mag} ∧ dz₁ (S-dual: uses z1 instead of z2)
G4_d5 = FormField(rank=4, dim=12)
z1_idx_ns5 = 6

for i, j, k in combinations(range(4), 3):
    l = [x for x in range(4) if x not in (i, j, k)][0]
    eps = levi_civita_4(i, j, k, l)
    if eps == 0:
        continue
    dH = sp.diff(H_ns5, trans_ns5[l])
    original = [y_base_ns5+i, y_base_ns5+j, y_base_ns5+k, z1_idx_ns5]
    sorted_idx = tuple(sorted(original))
    perm_map = [sorted_idx.index(x) for x in original]
    sort_sign = Permutation(perm_map).signature()
    G4_d5[sorted_idx] = sp.Rational(eps) * sort_sign * dH

print("Computing Ricci and form data...")
Ric_d5 = m_d5.ricci_tensor(simplify_func=cancel)
FF_d5 = form_contraction(G4_d5, m_d5)
norm_d5 = hf_ns5.substitute(cancel(form_norm_squared(G4_d5, m_d5)))
print(f"|G₄|²(D5) = {norm_d5}")

verify_enhanced_formula("D5", Ric_d5, m_d5, FF_d5, norm_d5, M_d5, hf_ns5,
                        [6, 7], [(0,'t'),(1,'x1'),(6,'z1'),(7,'z2'),(8,'y0'),(9,'y1')])


# ===================================================================
# SUMMARY
# ===================================================================
print("\n\n" + "="*70)
print("★★★ UNIVERSAL 12d FORMULA ★★★")
print("="*70)
print("""
  ℛ_{MN} = (1/2) [(G₄²)_{MN}/3! - (1/4)|G₄|² ĝ_{MN}]

where:
  ĝ_{MN} = g_{MN} + M̃_{MN}    (enhanced metric)
  M̃_{MN} = M_{ab} δ^a_M δ^b_N  (torus metric in 12d, zero outside torus)

Equivalently:
  ĝ_{mn} = g_{mn}  for 10d directions (λ_eff = 1/4)
  ĝ_{ab} = 2g_{ab}  for torus (λ_eff = 1/2)

The formula uses:
  λ = 1/4 = (p-1)/(D-2)|_{p=3, D=10}  (10d value, NOT 12d value 3/10)
  but with the enhanced metric ĝ = g + M̃ that doubles the torus trace.

Physical interpretation:
  The torus metric M̃ contributes an EXTRA trace term proportional to |G₄|²,
  reflecting the coupling between form flux and torus moduli.
""")
