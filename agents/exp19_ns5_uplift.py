"""Experiment 19: NS5-brane 12d uplift.

The NS5-brane tests the magnetic sector of the 12d equation.

10d NS5:
  String frame: ds² = ds²_{1,5} + H ds²_4  (H harmonic in R⁴)
  Dilaton: e^{2Φ} = H, i.e., Φ = (1/2)ln H
  H₃ = *_4 dH  (magnetic 3-form in transverse R⁴)

The NS5 H₃ is MAGNETIC: H_{ijk} = ε_{ijkl} ∂_l H / H²...
Actually for a p-form with H₃ = *_4 dH, the components are:
  H_{y_i y_j y_k} = ε_{ijkl} g^{ll}_{transverse} ∂_l H × √(det g_4)
In string frame g_4 = H δ, so √det = H², g^{ll} = 1/H.
  H_{ijk} = ε_{ijkl} (1/H) ∂_l H × H² = H ε_{ijkl} ∂_l H

Wait, let me be more careful. *_4 dH means:
  (dH)_l = ∂_l H
  (*_4 dH)_{ijk} = √g_4 ε_{ijkl} g^{ll} (dH)_l
For diagonal g_4 = H δ_{ll} (string frame):
  √g_4 = H², g^{ll} = 1/H
  (*_4 dH)_{ijk} = H² ε_{ijkl} (1/H) ∂_l H = H ε_{ijkl} ∂_l H

Hmm, but this should satisfy dH₃ = 0 and d(*_{10}H₃) ~ δ-function.
The Bianchi identity dH₃ = 0 means H₃ is closed. Let me verify.

Actually: for the NS5, the equation is:
  d(e^{-2Φ} *_{10} H₃) = 0 (Bianchi)  and  dH₃ = 0 (EOM from B-field)

With Φ = (1/2)ln H: e^{-2Φ} = 1/H.

12d uplift:
  M = diag(τ₂, 1/τ₂) = diag(e^{-Φ}, e^Φ) = diag(H^{-1/2}, H^{1/2})
  ds²_{12} = g^S + M dz² = ds²_{1,5} + H ds²_4 + H^{-1/2}dz₁² + H^{1/2}dz₂²

Note: this is the OPPOSITE torus assignment compared to F1!
  F1: M = diag(H^{1/2}, H^{-1/2})
  NS5: M = diag(H^{-1/2}, H^{1/2})
These are related by z₁↔z₂ (S-duality), as expected since NS5 and D1 are
the same type (opposite sector), and F1 and NS5 are electromagnetic duals.

G₄ = H₃ ∧ dz₂  (H₃ magnetic in R⁴)
  G₄_{y_i,y_j,y_k,z₂} = H ε_{ijkl} ∂_l H

Plan:
  A. Build the NS5 12d metric (string frame + Einstein frame)
  B. Compute Ricci tensor
  C. Build G₄ = H₃ ∧ dz₂ and compute form contraction
  D. Test the exp13 Einstein-frame equation
  E. Test the KK decomposition formulas from exp18
  F. Compare with F1 (electromagnetic duality)
"""
import sys; sys.path.insert(0, 'src')
import sympy as sp
from sympy import Rational as R, ln
from sugra import (HarmonicFunction, Metric,
                   FormField, exterior_derivative,
                   form_contraction, form_norm_squared)

print("=" * 70)
print("EXP 19: NS5-brane 12d uplift")
print("=" * 70)

# NS5: transverse R⁴ with harmonic function H = 1 + Q/r²
# We use 4 transverse coordinates y0,...,y3.
wv_coords = list(sp.symbols('t x1 x2 x3 x4 x5', real=True))
z1s, z2s = sp.symbols('z1 z2', real=True)
trans_coords = list(sp.symbols('y0:4', real=True))
coords = wv_coords + [z1s, z2s] + trans_coords
D = 12

hf = HarmonicFunction(transverse_coords=trans_coords)
H = sp.Function('H')(hf.r_expr)

# ===================================================================
# Part A: NS5 12d Einstein-frame metric
# ===================================================================
print("\nPART A: NS5 12d metric")
print("-" * 50)

# Einstein frame: ds² = H^{-1/4}ds²_{1,5} + H^{3/4}ds²_4
# Dilaton: Φ = (1/2)ln H → e^{-Φ} = H^{-1/2}
# Torus: M = diag(e^{-Φ}, e^Φ) = diag(H^{-1/2}, H^{1/2})

g_ef = sp.zeros(D, D)
# Worldvolume (6d)
g_ef[0, 0] = -H**R(-1, 4)  # t
for k in range(1, 6):
    g_ef[k, k] = H**R(-1, 4)  # x1,...,x5
# Torus
g_ef[6, 6] = H**R(-1, 2)  # z1 (= e^{-Φ})
g_ef[7, 7] = H**R(1, 2)   # z2 (= e^Φ)
# Transverse (4d)
for k in range(4):
    g_ef[8 + k, 8 + k] = H**R(3, 4)

print("NS5 Einstein-frame 12d metric:")
for i, name in [(0, 't'), (1, 'x1'), (6, 'z1'), (7, 'z2'), (8, 'y0')]:
    print(f"  g[{name}] = {g_ef[i,i]}")

metric_ef = Metric(g_ef, coords)
print("\nComputing 12d Ricci tensor (Einstein frame)...")
Ric_ef = metric_ef.ricci_tensor(simplify_func=sp.cancel)

print("\nRicci tensor per block (with ℛ/g):")
for i, name in [(0, 't'), (1, 'x1'), (6, 'z1'), (7, 'z2'), (8, 'y0'), (9, 'y1')]:
    ric_val = hf.substitute(sp.cancel(Ric_ef[i, i]))
    g_val = hf.substitute(sp.cancel(g_ef[i, i]))
    ratio = sp.cancel(ric_val / g_val) if g_val != 0 else None
    print(f"  ℛ[{name}] = {ric_val}  (ℛ/g = {ratio})")

# ===================================================================
# Part B: Build G₄ = H₃ ∧ dz₂
# ===================================================================
print("\n" + "=" * 70)
print("PART B: G₄ form field for NS5")
print("-" * 50)

# For NS5, H₃ is magnetic: H₃ = *_4 dH (Hodge dual in transverse R⁴)
# In the Einstein-frame transverse metric g_trans = H^{3/4} δ:
# *_{EF} is with respect to the EF metric.
#
# Actually, let me build H₃ directly. The components of H₃ = *_4^{flat} dH
# in flat transverse space (before warping) are:
# H_{ijk} = ε_{ijkl} ∂_l H  (flat-space Hodge dual)
#
# But in the warped metric, the proper Hodge dual is:
# (*_4 dH)_{ijk} = √g_4 ε_{ijkl} g^{ll} ∂_l H
#
# For Einstein frame: g_4 = H^{3/4}δ, √g_4 = H³, g^{ll} = H^{-3/4}
# (*_4 dH)_{ijk} = H³ × H^{-3/4} × ε_{ijkl} ∂_l H = H^{9/4} ε_{ijkl} ∂_l H
#
# Hmm, but what we really want is the POTENTIAL whose exterior derivative
# is H₃. For the magnetic case, H₃ is closed (dH₃ = 0 when H is harmonic),
# and locally H₃ = dB₂.

# Actually, the simplest approach: build the 3-form H₃ by its components
# and then form G₄ = H₃ ∧ dz₂.

# For the NS5, the standard expression is:
# H₃ = *_4^{flat} dH where *_4^{flat} is w.r.t. flat R⁴ metric
# So H_{y_i,y_j,y_k} = ε_{ijkl} ∂_{y_l} H

# But this needs to be "covariantized" for the curved metric?
# No — H₃ is a 3-FORM (tensor), defined with flat components.
# It's the same tensor regardless of metric; the metric only enters
# when computing norms and contractions.

# The Levi-Civita SYMBOL ε_{0123} = +1 (not the tensor).
# H₃ = dB₂ is the field strength; its components in coord basis are:
# H_{y_i y_j y_k} = ε_{ijkl} ∂_{y_l} H

# For G₄ = H₃ ∧ dz₂:
# G₄_{y_i,y_j,y_k,z₂} = H₃_{y_i,y_j,y_k} = ε_{ijkl} ∂_{y_l} H

# Note: ε_{ijkl} is the Levi-Civita SYMBOL (purely combinatorial +1/-1/0).
# In 4d: ε_{0123} = +1, and other components by antisymmetry.

from itertools import combinations

G4 = FormField(rank=4, dim=D)

# Populate G₄ components. Indices: 8,9,10,11 are y0,y1,y2,y3; 7 is z2.
# For each triple (i,j,k) of transverse indices, and the remaining l:
# G₄_{y_i,y_j,y_k,z₂} = ε_{ijkl} ∂_{y_l} H
trans_idx = list(range(4))  # internal indices 0,1,2,3 for y

def levi_civita_4(i, j, k, l):
    """Levi-Civita symbol in 4d."""
    perm = [i, j, k, l]
    if len(set(perm)) < 4:
        return 0
    # Count inversions
    sign = 1
    for a in range(4):
        for b in range(a+1, 4):
            if perm[a] > perm[b]:
                sign *= -1
    return sign

for i, j, k in combinations(trans_idx, 3):
    # Find the remaining index l
    l = [x for x in trans_idx if x not in (i, j, k)][0]
    eps = levi_civita_4(i, j, k, l)
    if eps == 0:
        continue
    # G₄_{8+i, 8+j, 8+k, 7} = ε_{ijkl} ∂_{y_l} H
    y_l = trans_coords[l]
    dH = sp.diff(H, y_l)
    # The form indices in 12d are: y_i=8+i, y_j=8+j, y_k=8+k, z₂=7
    idx = tuple(sorted([8+i, 8+j, 8+k, 7]))
    # Need to determine sign from sorting
    original = [8+i, 8+j, 8+k, 7]
    sorted_idx = sorted(original)
    # Count sort permutation sign
    from sympy.combinatorics import Permutation
    perm_map = [sorted_idx.index(x) for x in original]
    sort_sign = Permutation(perm_map).signature()
    G4[tuple(sorted_idx)] = sp.Rational(eps) * sort_sign * dH

# Print some components
print("G₄ components (nonzero):")
count = 0
for idx in sorted(G4.nonzero_components.keys()):
    val = sp.cancel(G4[idx])
    if val != 0 and count < 6:
        idx_names = [coords[i].name if hasattr(coords[i], 'name') else str(coords[i]) for i in idx]
        print(f"  G₄[{','.join(idx_names)}] = {val}")
        count += 1

# ===================================================================
# Part C: Form contractions
# ===================================================================
print("\n" + "=" * 70)
print("PART C: Form contractions and norm")
print("-" * 50)

FF = form_contraction(G4, metric_ef)
S = form_norm_squared(G4, metric_ef)

S_val = hf.substitute(sp.cancel(S))
print(f"|G₄|² = {S_val}")

print("\nFF_{MN}/3! per block:")
for i, name in [(0, 't'), (1, 'x1'), (6, 'z1'), (7, 'z2'), (8, 'y0'), (9, 'y1')]:
    ff_val = hf.substitute(sp.cancel(FF[i, i]))
    print(f"  FF[{name}]/3! = {sp.cancel(ff_val / 6)}")

# ===================================================================
# Part D: Test Einstein-frame equation (exp13 formula)
# ===================================================================
print("\n" + "=" * 70)
print("PART D: Test exp13 formula (Einstein frame, λ=1/4)")
print("-" * 50)

# exp13 formula: ℛ_{MN} = (1/2)[FF/3! - (1/4)|G₄|²g] + (1/2)|∂Ψ|²P
# Let's first test the simpler: ℛ_{MN} = (1/2)[FF_{MN}/3! - λ|G₄|²g_{MN}]
# with λ = 1/4 and see what's left.

lambda_val = R(1, 4)
print(f"\nTesting T_{'{MN}'} = (1/2)[FF/3! - ({lambda_val})|G₄|²g]:")

residuals = {}
for i, name in [(0, 't'), (1, 'x1'), (6, 'z1'), (7, 'z2'), (8, 'y0'), (9, 'y1')]:
    ric_val = hf.substitute(sp.cancel(Ric_ef[i, i]))
    ff_val = hf.substitute(sp.cancel(FF[i, i]))
    g_val = hf.substitute(sp.cancel(g_ef[i, i]))
    T_val = sp.cancel(R(1, 2) * (ff_val / 6 - lambda_val * S_val * g_val))
    diff = sp.cancel(ric_val - T_val)
    residuals[name] = diff
    status = '✓' if diff == 0 else f'residual={diff}'
    print(f"  [{name}]: ℛ={ric_val}, T={T_val}, {status}")

# Check if residuals have a pattern (like the torus projection from exp13)
print("\nResidual analysis (R - T with λ=1/4):")
for name in residuals:
    r = residuals[name]
    if r != 0:
        print(f"  D[{name}] = {r}")

# ===================================================================
# Part E: KK decomposition check
# ===================================================================
print("\n" + "=" * 70)
print("PART E: KK decomposition (string frame)")
print("-" * 50)

# String-frame NS5 metric
g_sf_12 = sp.zeros(D, D)
g_sf_12[0, 0] = -1          # t (flat worldvolume)
for k in range(1, 6):
    g_sf_12[k, k] = 1       # x1,...,x5
g_sf_12[6, 6] = H**R(-1, 2) # z1
g_sf_12[7, 7] = H**R(1, 2)  # z2
for k in range(4):
    g_sf_12[8 + k, 8 + k] = H  # y_k (transverse, warped)

M_ns5 = sp.Matrix([[H**R(-1, 2), 0], [0, H**R(1, 2)]])

print("String-frame NS5 12d metric:")
for i, name in [(0, 't'), (6, 'z1'), (7, 'z2'), (8, 'y0')]:
    print(f"  g[{name}] = {g_sf_12[i,i]}")

metric_sf_12 = Metric(g_sf_12, coords)
print("\nComputing string-frame 12d Ricci...")
Ric_sf_12 = metric_sf_12.ricci_tensor(simplify_func=sp.cancel)

# 10d string-frame base metric (no torus)
coords_10 = wv_coords + trans_coords
g_sf_10 = sp.zeros(10, 10)
g_sf_10[0, 0] = -1
for k in range(1, 6):
    g_sf_10[k, k] = 1
for k in range(4):
    g_sf_10[6 + k, 6 + k] = H

metric_sf_10 = Metric(g_sf_10, coords_10)
Ric_sf_10 = metric_sf_10.ricci_tensor(simplify_func=sp.cancel)

# KK formula: ℛ_{mn}(12d) = R_{mn}(10d) + (1/4)Tr(∂_m M⁻¹·∂_n M)
M_ns5_inv = sp.cancel(M_ns5.inv())

print("\nVerify KK formula: ℛ_{mn} = R_{mn}^{SF} + (1/4)Tr(∂_m M⁻¹·∂_n M)")
kk_pass = True
for i_12, i_10, name in [(0, 0, 't'), (1, 1, 'x1'), (8, 6, 'y0'), (9, 7, 'y1')]:
    R12 = hf.substitute(sp.cancel(Ric_sf_12[i_12, i_12]))
    R10 = hf.substitute(sp.cancel(Ric_sf_10[i_10, i_10]))

    tr_val = sp.Integer(0)
    for a in range(2):
        for b in range(2):
            dMinv = sp.diff(M_ns5_inv[a, b], coords_10[i_10])
            dM = sp.diff(M_ns5[b, a], coords_10[i_10])
            tr_val += dMinv * dM
    tr_val = hf.substitute(sp.cancel(R(1, 4) * tr_val))

    predicted = sp.cancel(R10 + tr_val)
    diff = sp.cancel(R12 - predicted)
    status = '✓' if diff == 0 else f'✗ diff={diff}'
    if diff != 0:
        kk_pass = False
    print(f"  [{name}]: R12={R12}, R10+KK={predicted}, {status}")

print(f"\n  10d KK: {'★ VERIFIED ★' if kk_pass else 'FAILED'}")

# Torus KK (corrected Pope formula from exp18c)
sqrt_g_ns5 = sp.Integer(1)
# det(g^S_10) for NS5: g_tt×...×g_55×g_y0...g_y3 = (-1)^1×(1)^5×H^4 = -H^4
# √|det| = H²
sqrt_g_sf = H**2

print("\nVerify torus KK: ℛ_{ab} = -(1/2) M_{ac} div(M⁻¹∂M)^c_b")

def pope_mixed_ns5(c, b):
    result = sp.Integer(0)
    for nu in range(10):
        g_inv = sp.cancel(1 / g_sf_10[nu, nu])
        MinvdM = sp.Integer(0)
        for e in range(2):
            MinvdM += M_ns5_inv[c, e] * sp.diff(M_ns5[e, b], coords_10[nu])
        MinvdM = sp.cancel(MinvdM)
        if MinvdM == 0:
            continue
        flux = sp.cancel(sqrt_g_sf * g_inv * MinvdM)
        result += sp.diff(flux, coords_10[nu])
    return sp.cancel(result / sqrt_g_sf)

torus_pass = True
for a, b, name in [(0, 0, 'z1'), (1, 1, 'z2'), (0, 1, 'z1z2')]:
    predicted = sp.Integer(0)
    for c in range(2):
        div_cb = pope_mixed_ns5(c, b)
        predicted += M_ns5[a, c] * div_cb
    predicted = hf.substitute(sp.cancel(-R(1, 2) * predicted))
    actual = hf.substitute(sp.cancel(Ric_sf_12[a+6, b+6]))
    diff = sp.cancel(actual - predicted)
    status = '✓' if diff == 0 else f'✗ diff={diff}'
    if diff != 0:
        torus_pass = False
    print(f"  [{name}]: actual={actual}, predicted={predicted}, {status}")

print(f"\n  Torus KK: {'★ VERIFIED ★' if torus_pass else 'FAILED'}")

# ===================================================================
# Part F: 10d string-frame equation for NS5
# ===================================================================
print("\n" + "=" * 70)
print("PART F: 10d string-frame equation for NS5")
print("-" * 50)

# NS5 in string frame: R^S + 2∇∂Φ = (1/4)(H₃²)
# where Φ = (1/2)ln H and H₃ = *_4^{flat} dH

# Build H₃ form in 10d
H3 = FormField(rank=3, dim=10)
for i, j, k in combinations(range(4), 3):
    l = [x for x in range(4) if x not in (i, j, k)][0]
    eps = levi_civita_4(i, j, k, l)
    if eps == 0:
        continue
    y_l = trans_coords[l]
    dH = sp.diff(H, y_l)
    idx_10 = tuple(sorted([6+i, 6+j, 6+k]))
    original = [6+i, 6+j, 6+k]
    sorted_idx = sorted(original)
    from sympy.combinatorics import Permutation
    perm_map = [sorted_idx.index(x) for x in original]
    sort_sign = Permutation(perm_map).signature()
    H3[tuple(sorted_idx)] = sp.Rational(eps) * sort_sign * dH

FF_H3 = form_contraction(H3, metric_sf_10)
S_H3 = form_norm_squared(H3, metric_sf_10)
S_H3_val = hf.substitute(sp.cancel(S_H3))
print(f"|H₃|²_{'{SF}'} = {S_H3_val}")

print("\n(H₃²)_{mn}/2! per block:")
for i, name in [(0, 't'), (1, 'x1'), (6, 'y0'), (7, 'y1')]:
    ff_val = hf.substitute(sp.cancel(FF_H3[i, i]))
    print(f"  (H₃²)[{name}]/2 = {sp.cancel(ff_val / 2)}")

# ∇∂Φ in string frame
Phi_ns5 = R(1, 2) * ln(H)
chris_sf = metric_sf_10.christoffel(simplify_func=sp.cancel)

print("\n10d equation check: R^S + 2∇∂Φ = (1/4)(H₃²)_{mn}")
eq_pass = True
for i_10, name in [(0, 't'), (1, 'x1'), (6, 'y0'), (7, 'y1')]:
    R_val = hf.substitute(sp.cancel(Ric_sf_10[i_10, i_10]))

    # ∇_i∂_iΦ = ∂²Φ/∂x_i² - Σ_ρ Γ^ρ_{ii} ∂_ρΦ
    d2Phi = sp.diff(sp.diff(Phi_ns5, coords_10[i_10]), coords_10[i_10])
    conn_term = sp.Integer(0)
    for rho in range(10):
        key = (rho, i_10, i_10)
        if key in chris_sf:
            conn_term += chris_sf[key] * sp.diff(Phi_ns5, coords_10[rho])
    nabla_Phi = hf.substitute(sp.cancel(d2Phi - conn_term))

    ff_val = hf.substitute(sp.cancel(FF_H3[i_10, i_10]))

    lhs = sp.cancel(R_val + 2 * nabla_Phi)
    rhs = sp.cancel(R(1, 4) * ff_val)
    diff = sp.cancel(lhs - rhs)
    status = '✓' if diff == 0 else f'✗ diff={diff}'
    if diff != 0:
        eq_pass = False
    print(f"  [{name}]: LHS={lhs}, RHS={rhs}, {status}")

print(f"\n  10d SF equation: {'★ VERIFIED ★' if eq_pass else 'FAILED'}")

# ===================================================================
# Part G: Summary
# ===================================================================
print("\n" + "=" * 70)
print("PART G: Summary")
print("=" * 70)
print("""
NS5-brane 12d uplift results:

1. 12d metric (Einstein frame):
   ds² = H^{-1/4}ds²_{1,5} + H^{3/4}ds²_4 + H^{-1/2}dz₁² + H^{1/2}dz₂²
   Torus: M = diag(H^{-1/2}, H^{1/2})  [opposite to F1]

2. KK decomposition (string frame):
   ℛ_{mn}(12d) = R_{mn}(10d,SF) + (1/4)Tr(∂M⁻¹·∂M)  [same formula as F1]
   ℛ_{ab}(12d) = corrected Pope formula  [same formula as F1]

3. 10d equation: R^S + 2∇∂Φ = (1/4)(H₃²)  [same form as F1, with Φ=(1/2)lnH]

4. The NS5 uses MAGNETIC H₃ = *_4 dH (legs in transverse R⁴)
   vs F1's ELECTRIC H₃ = d(H⁻¹) ∧ vol_{wv} (legs in worldvolume)
""")
