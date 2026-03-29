"""Experiment 16: (p,q)-string uplift — non-diagonal torus metric.

Tests the SL(2,Z)-covariant 12d equation with BOTH doublet terms active.

Strategy:
  Apply SL(2,Z) transformation Λ = [[1,0],[1,1]] to the F1 solution.
  This maps F1 → (1,1)-string with:
    - Non-diagonal torus metric M = Λ M₀ Λ^T
    - Both F₃ components active: F₃'^1 = -H₃, F₃'^2 = H₃
    - G₄ = F₃'^a ∧ dz'_a = -H₃∧dz₁ + H₃∧dz₂  (both terms active)

  Since the 12d metric ds²₁₂ = ds²₁₀ + M_{ab}dz^a dz^b has a CROSS TERM
  dz₁dz₂, the metric is non-diagonal. This tests:
    1. The general Christoffel/Ricci code for non-diagonal metrics
    2. Form contractions with non-diagonal metrics
    3. Tr(∂M ∂M^{-1}) for non-diagonal M
    4. The full SL(2,Z)-covariant equation

Verify:
  ℛ_{MN} = (1/2)[G₄²_{MN}/3! − (1/4)|G₄|² 𝒢_{MN}]
          − (1/4) Tr(∂_M ℳ · ∂_N ℳ⁻¹) · Π_{MN}

  with the (1,1)-string data.
"""
import sys; sys.path.insert(0, 'src')
import sympy as sp
from sympy import Rational as R
from sugra import (HarmonicFunction, Metric,
                   FormField, exterior_derivative,
                   form_contraction, form_norm_squared)

# ===================================================================
# Setup: (1,1)-string 12d metric
# ===================================================================
print("=" * 70)
print("EXP 16: (1,1)-string 12d uplift (non-diagonal torus)")
print("=" * 70)

wv_coords = list(sp.symbols('t x1', real=True))
z1s, z2s = sp.symbols('z1 z2', real=True)
harmonic_coords = list(sp.symbols('y0:8', real=True))
coords = wv_coords + [z1s, z2s] + harmonic_coords  # 12 coords

hf = HarmonicFunction(transverse_coords=harmonic_coords)
H = sp.Function('H')(hf.r_expr)
y0 = sp.Symbol('y0', real=True)

# F1 torus metric: M₀ = diag(H^{1/2}, H^{-1/2})
# (z1 = v has g_{z1}=H^{1/2}, z2 = u has g_{z2}=H^{-1/2})
M0 = sp.Matrix([[H**R(1, 2), 0],
                 [0, H**R(-1, 2)]])

# SL(2,Z) transformation: Λ = [[1,0],[1,1]]
# Maps F1 (1,0)-string → (1,1)-string
Lam = sp.Matrix([[1, 0], [1, 1]])
print(f"\nSL(2,Z) matrix Λ = {Lam.tolist()}")
print(f"det(Λ) = {Lam.det()}")

# New torus metric: M = Λ M₀ Λ^T
M_pq = sp.cancel(Lam * M0 * Lam.T)
print(f"\nTorus metric M = Λ M₀ Λ^T =")
for i in range(2):
    print(f"  [{sp.cancel(M_pq[i, 0])}, {sp.cancel(M_pq[i, 1])}]")

# Build the full 12d metric matrix
# ds²₁₂ = H^{-3/4} ds²_{1,1} + M_{ab} dz^a dz^b + H^{1/4} ds²_8
D = 12
g_matrix = sp.zeros(D, D)

# Worldvolume: t, x1 (indices 0, 1)
g_matrix[0, 0] = -H**R(-3, 4)  # -dt²
g_matrix[1, 1] = H**R(-3, 4)   # +dx1²

# Torus: z1, z2 (indices 2, 3) — from M_{ab}
g_matrix[2, 2] = M_pq[0, 0]
g_matrix[2, 3] = M_pq[0, 1]
g_matrix[3, 2] = M_pq[1, 0]
g_matrix[3, 3] = M_pq[1, 1]

# Transverse: y0...y7 (indices 4..11)
for k in range(8):
    g_matrix[4 + k, 4 + k] = H**R(1, 4)

print("\nConstructing 12d Metric (non-diagonal)...")
metric_pq = Metric(g_matrix, coords)
print(f"  Diagonal: {metric_pq.is_diagonal}")

print("Computing 12d Ricci tensor...")
Ric_pq = metric_pq.ricci_tensor(simplify_func=sp.cancel)
print("Done.")

# ===================================================================
# Part A: G₄ for the (1,1)-string
# ===================================================================
print("\n" + "=" * 70)
print("PART A: 4-form G₄ for (1,1)-string")
print("=" * 70)

# Under Λ = [[1,0],[1,1]], the 3-form doublet transforms as:
# F₃'^a = (Λ^{-T})^a_b F₃^b
# Λ^{-1} = [[1,0],[-1,1]], Λ^{-T} = [[1,-1],[0,1]]
#
# For F1: F₃^1 = 0, F₃^2 = H₃ (NSNS 3-form, B_{tx}=H^{-1})
# → F₃'^1 = (Λ^{-T})^1_2 H₃ = -H₃
#   F₃'^2 = (Λ^{-T})^2_2 H₃ = H₃
#
# G₄ = F₃'^a ∧ dz_a = -H₃ ∧ dz₁ + H₃ ∧ dz₂
#
# Since H₃ = d(H^{-1}) ∧ dt ∧ dx₁, the potential is:
#   C_{t,x1,z1} = -H^{-1}  (from -H₃ ∧ dz₁)
#   C_{t,x1,z2} = +H^{-1}  (from +H₃ ∧ dz₂)

print("\nF₃ doublet after SL(2,Z):")
print("  F₃'^1 = -H₃  (RR)")
print("  F₃'^2 = +H₃  (NSNS)")
print("  G₄ = -H₃∧dz₁ + H₃∧dz₂ = H₃∧(dz₂ - dz₁)")

C_pq = FormField(rank=3, dim=12)
C_pq[(0, 1, 2)] = -1 / H  # C_{t,x1,z1} = -H^{-1}
C_pq[(0, 1, 3)] = 1 / H   # C_{t,x1,z2} = +H^{-1}
F_pq = exterior_derivative(C_pq, coords)

print("\nNon-zero G₄ components:")
for idx, val in F_pq.nonzero_components.items():
    names = [str(coords[i]) for i in idx]
    print(f"  G₄[{','.join(names)}] = {val}")

# ===================================================================
# Part B: Compute form contractions and verify equation
# ===================================================================
print("\n" + "=" * 70)
print("PART B: Form contractions with non-diagonal metric")
print("=" * 70)

FF_pq = form_contraction(F_pq, metric_pq)
S_pq = form_norm_squared(F_pq, metric_pq)
S_pq_val = hf.substitute(sp.cancel(S_pq))
print(f"\n|G₄|² = {S_pq_val}")

# Compare with F1 value: |G₄|² should be the SAME (G₄ is SL(2,Z) invariant)
# F1 had |G₄|² = -H'²/H^{9/4}
S_f1_expected = -hf.Hp**2 / hf.H**R(9, 4)
print(f"|G₄|²(F1) = {S_f1_expected}")
print(f"Match: {sp.cancel(S_pq_val - S_f1_expected) == 0}")

# ===================================================================
# Part C: Torus moduli kinetic term Tr(∂M ∂M^{-1})
# ===================================================================
print("\n" + "=" * 70)
print("PART C: Torus moduli kinetic term (non-diagonal M)")
print("=" * 70)

M_inv_pq = sp.cancel(M_pq.inv())
print(f"\nM^{{-1}} =")
for i in range(2):
    print(f"  [{sp.cancel(M_inv_pq[i, 0])}, {sp.cancel(M_inv_pq[i, 1])}]")

# det(M) = det(Λ)² det(M₀) = 1·1 = 1
det_M = sp.cancel(M_pq.det())
print(f"\ndet(M) = {det_M}")

# Tr(∂_M ℳ · ∂_N ℳ^{-1}) for each pair of 12d directions
# Only non-zero when M,N are transverse (y-directions), since M depends on H(r)
# and r depends only on y-coords.

# Compute ∂_k M and ∂_k M^{-1} for a generic transverse direction
# M = M(H), H = H(r), r = r(y_k)
# ∂_k M = (dM/dH)(dH/dr)(∂r/∂y_k) = M'(H) · H'(r) · (y_k/r)

# For the scalar kinetic tensor contracted over torus:
# (1/4) Tr(∂_M M · ∂_N M^{-1}) for M,N both y-indices:
# = (1/4) Σ_{a,b} (∂_M M)_{ab} (∂_N M^{-1})_{ba}
# = (1/4) (y_M/r)(y_N/r) (H'/H)² × Σ_{a,b} [H d(M_{ab})/dH] [H d(M^{-1}_{ba})/dH] / H²

# Let's compute Tr(dM/dH · dM^{-1}/dH)
dM_dH = sp.cancel(sp.diff(M_pq, H))
dMinv_dH = sp.cancel(sp.diff(M_inv_pq, H))

print(f"\ndM/dH =")
for i in range(2):
    print(f"  [{sp.cancel(dM_dH[i, 0])}, {sp.cancel(dM_dH[i, 1])}]")

print(f"\ndM^(-1)/dH =")
for i in range(2):
    print(f"  [{sp.cancel(dMinv_dH[i, 0])}, {sp.cancel(dMinv_dH[i, 1])}]")

trace_product = sp.cancel(sp.trace(dM_dH * dMinv_dH))
print(f"\nTr(dM/dH · dM^(-1)/dH) = {trace_product}")

# For F1 (diagonal M₀ = diag(H^{1/2}, H^{-1/2})):
dM0_dH = sp.cancel(sp.diff(M0, H))
dM0inv_dH = sp.cancel(sp.diff(M0.inv(), H))
trace_f1 = sp.cancel(sp.trace(dM0_dH * dM0inv_dH))
print(f"Tr(dM₀/dH · dM₀^(-1)/dH) [F1] = {trace_f1}")
print(f"Match: {sp.cancel(trace_product - trace_f1) == 0}")

# The full kinetic term for transverse direction y_k:
# (1/4) Tr(∂_k M · ∂_k M^{-1}) = (1/4)(H')²(y_k²/r²) × Tr(dM/dH · dM^{-1}/dH)
# For the invariant contraction |∂Ψ|²:
# (1/4) g^{yk,yk} Σ_k (H')²(y_k²/r²) × trace_product
# = (1/4)(1/H^{1/4})(H')² × trace_product
# This should equal (1/2)(∂σ)² = H'²/(4H^{9/4}) as for F1

kinetic_scalar = R(1, 4) * hf.Hp**2 / hf.H**R(1, 4) * trace_product
kinetic_scalar = sp.cancel(kinetic_scalar)
print(f"\n(1/4)|∂|² × Tr = {kinetic_scalar}")

kinetic_f1 = hf.Hp**2 / (4 * hf.H**R(9, 4))
# We expect: (1/4) × Tr × (H')²/H^{1/4} should simplify to same thing
# trace_product should be -1/(2H²) for both F1 and (1,1)
# then (1/4)(H'²/H^{1/4})(-1/(2H²)) = -H'²/(8H^{9/4})
# Hmm, need to check sign conventions

print(f"Expected (F1 value) = -(1/2)(∂σ)² = {sp.cancel(-kinetic_f1/2)}")

# ===================================================================
# Part D: The NAIVE T^{G₄}(λ=1/4) formula FAILS for non-diagonal M
# ===================================================================
print("\n" + "=" * 70)
print("PART D: T^{G₄}(λ=1/4) fails — because |G₄|² is NOT SL(2,Z) invariant")
print("=" * 70)

lam = R(1, 4)

print("""
KEY INSIGHT: |G₄|² = G^{MNPQ}G_{MNPQ} uses the 12d metric to raise indices.
Since the torus metric changed (M → ΛMΛ^T), raising z-indices uses M^{-1},
and |G₄|² changes even though G₄ as a form is SL(2,Z) invariant.

The T^{G₄}(λ=1/4) formula contracts z-indices with M^{-1}_{ab},
but the correct 10d IIB doublet equation contracts with M_{ab}.
These differ for non-diagonal M.
""")

print("Quick check — T^{G₄}(λ=1/4) for t-component:")
R_t = hf.substitute(sp.cancel(Ric_pq[0, 0]))
ff_t = hf.substitute(sp.cancel(FF_pq[0, 0]))
g_t = hf.substitute(sp.cancel(metric_pq.matrix[0, 0]))
T_G4_t = sp.cancel(R(1, 2) * (ff_t / 6 - lam * S_pq_val * g_t))
diff_t = sp.cancel(R_t - T_G4_t)
print(f"  diff[t] = {diff_t}  (should be 0 if formula worked)")

# ===================================================================
# Part E: CORRECT formulation — doublet contraction with M_{ab}
# ===================================================================
print("\n" + "=" * 70)
print("PART E: Correct SL(2,Z)-covariant equation using doublet F₃^a")
print("=" * 70)

print("""
The correct 10d IIB Einstein equation (SL(2,Z) covariant) is:
  R_{mn}(10d) = (1/4)Tr(∂_m M ∂_n M^{-1})
              + (1/4) M_{ab} [F₃^a·F₃^b_{mn}/2! − (1/4)(F₃^a·F₃^b)g_{mn}]

where M_{ab} = torus metric contracts the doublet (NOT M^{-1}!).

In 12d, KK gives ℛ_{mn} = R_{mn}(10d) − (1/4)Tr(∂_m M ∂_n M^{-1})
→ scalar kinetic terms cancel for 10d directions, leaving:

  ★ ℛ_{mn}(12d) = (1/4) M_{ab} [(F₃^a·F₃^b)_{mn}/2! − (1/4)(F₃^a·F₃^b)g_{mn}]

where F₃^a·F₃^b means the 10d contraction of the rank-3 forms.
""")

# Build the two F₃ forms in 10d-style (but embedded in 12d coordinates)
# F₃^1 = -H₃, F₃^2 = +H₃
# H₃ = d(H^{-1} dt ∧ dx₁) → H₃_{t,x1,yk} = -(H'/H²)(yk/r)

# We need (F₃^a·F₃^b)_{mn} contracted using the 10D metric
# (i.e., raising only the 10d indices, NOT the z-indices)
# For our case, F₃^1 = -H₃ and F₃^2 = +H₃, so:
#   F₃^a · F₃^b = c_a c_b (H₃·H₃)
# where c = (-1, +1)

# First, compute H₃ contractions using the 10d metric embedded in 12d
# H₃ lives in the 10d subspace: indices {0,1,4,...,11}
# H₃_{t,x1,yk} for k=0,...,7

# The 10d metric (for the 10d directions):
# g_{tt} = -H^{-3/4}, g_{x1} = H^{-3/4}, g_{yk} = H^{1/4}
# These are the SAME as in the F1 metric — they don't change under SL(2,Z)!

# Build H₃ as a 10d form
C_h3 = FormField(rank=2, dim=12)
C_h3[(0, 1)] = 1 / H  # B_{t,x1} = H^{-1}
H3_form = exterior_derivative(C_h3, coords)

# Build a "10d metric" for contracting H₃ — use the 12d metric but only
# the 10d block. Since we only contract (t,x1,y0,...,y7), the off-diagonal
# z-terms don't matter.
print("Computing H₃ contractions in 10d subspace...")
FF_H3_10d = form_contraction(H3_form, metric_pq)
S_H3_10d = form_norm_squared(H3_form, metric_pq)
S_H3_val = hf.substitute(sp.cancel(S_H3_10d))
print(f"  |H₃|²(using 12d metric) = {S_H3_val}")

# For the F1 case: |H₃|² should be -H'²/H^{11/4} (10d value)
S_H3_f1 = -hf.Hp**2 / hf.H**R(11, 4)
print(f"  |H₃|²(F1, 10d) = {S_H3_f1}")
print(f"  Match: {sp.cancel(S_H3_val - S_H3_f1) == 0}")

# Now the doublet-contracted stress energy:
# T^{doublet}_{mn} = (1/4) M_{ab} c_a c_b [(H₃²)_{mn}/2 − (1/4)|H₃|²g_{mn}]
#
# With c = (-1, +1) and M = [[√H, √H], [√H, (H+1)/√H]]:
# c^T M c = (-1)²M₁₁ + 2(-1)(+1)M₁₂ + (+1)²M₂₂
#         = M₁₁ - 2M₁₂ + M₂₂

cTMc = sp.cancel(M_pq[0, 0] - 2 * M_pq[0, 1] + M_pq[1, 1])
print(f"\nc^T M c = M₁₁ − 2M₁₂ + M₂₂ = {cTMc}")

# For F1: c=(0,1), c^T M₀ c = M₀_{22} = H^{-1/2}
cTM0c = M0[1, 1]
print(f"c^T M₀ c (F1) = {sp.cancel(cTM0c)}")
print(f"Match: {sp.cancel(cTMc - cTM0c) == 0}")
print("→ c^T M c is SL(2,Z) INVARIANT (as expected: charge bilinear)")

# Now verify the equation for ALL 10d directions
print("\n--- Verifying: ℛ_{mn} = (1/4)(c^TMc)[(H₃²)_{mn}/2 − (1/4)|H₃|²g] ---\n")

all_10d_pass = True
for i, name in [(0, 't'), (1, 'x1'), (4, 'y0'), (5, 'y1')]:
    R_val = hf.substitute(sp.cancel(Ric_pq[i, i]))
    ff_H3_val = hf.substitute(sp.cancel(FF_H3_10d[i, i]))
    g_val = hf.substitute(sp.cancel(metric_pq.matrix[i, i]))
    cTMc_val = hf.substitute(sp.cancel(cTMc))

    T_doublet = sp.cancel(R(1, 4) * cTMc_val * (ff_H3_val / 2 - R(1, 4) * S_H3_val * g_val))
    diff = sp.cancel(R_val - T_doublet)

    status = '✓' if diff == 0 else f'✗ diff={diff}'
    if diff != 0:
        all_10d_pass = False
    print(f"  [{name}] (i={i}): {status}")
    if diff != 0 and i <= 5:
        print(f"    ℛ = {R_val}")
        print(f"    T = {T_doublet}")

# Check remaining y-coords
for k in range(2, 8):
    i = 4 + k
    R_val = hf.substitute(sp.cancel(Ric_pq[i, i]))
    ff_H3_val = hf.substitute(sp.cancel(FF_H3_10d[i, i]))
    g_val = hf.substitute(sp.cancel(metric_pq.matrix[i, i]))
    cTMc_val = hf.substitute(sp.cancel(cTMc))
    T_doublet = sp.cancel(R(1, 4) * cTMc_val * (ff_H3_val / 2 - R(1, 4) * S_H3_val * g_val))
    diff = sp.cancel(R_val - T_doublet)
    if diff != 0:
        all_10d_pass = False
        print(f"  [y{k}] (i={i}): ✗ diff={diff}")

print(f"\n  10d directions: {'★ ALL PASS ★' if all_10d_pass else 'SOME FAILED'}")

# ===================================================================
# Part F: Torus directions
# ===================================================================
print("\n" + "=" * 70)
print("PART F: Torus directions — doublet formulation + scalar kinetic")
print("=" * 70)

# For torus directions, the KK reduction gives:
# ℛ_{ab}(12d) = KK_{ab} (no 10d Ricci for internal directions)
# The KK term includes the scalar kinetic energy from M_{ab}.
#
# Also, for torus indices, the form G₄ has specific FF_{ab} structure.
# Let's compute the residual ℛ_{ab} - T^{doublet}_{ab} and identify the pattern.

# First: does the doublet form stress-energy contribute to torus indices?
# F₃^a lives in 10d, so (F₃^a·F₃^b)_{z_c, z_d} = 0 (no z-indices in F₃)
# The only contribution is from the trace term: -(1/4)|F₃^a·F₃^b|g_{z_c z_d}
# But wait — the FULL G₄ stress-energy does have torus components.

# Let's compute T^{G₄} at the torus indices as a reference
print("\nG₄ stress-energy at torus indices:")
for i in range(2, 4):
    for j in range(i, 4):
        ff_val = hf.substitute(sp.cancel(FF_pq[i, j]))
        ni = 'z1' if i == 2 else 'z2'
        nj = 'z1' if j == 2 else 'z2'
        print(f"  FF(G₄)[{ni},{nj}] = {ff_val}")

# For the doublet formulation at torus indices:
# The 10d equation for the scalar sector is:
#   0 = □M_{ab} + (form terms in M)
# This KK-lifts to the torus components of ℛ.
#
# The full torus ℛ should be sourced by:
#   ℛ_{ab} = T^{G₄}_{ab}(λ) + T^{scalar}_{ab}
# where T^{scalar} encodes the M equation of motion.

# Let's compute ℛ and form residuals at torus block
print("\nTorus block residuals:")
torus_pass = True
for i in range(2, 4):
    for j in range(i, 4):
        R_val = hf.substitute(sp.cancel(Ric_pq[i, j]))
        ff_val = hf.substitute(sp.cancel(FF_pq[i, j]))
        g_val = hf.substitute(sp.cancel(metric_pq.matrix[i, j]))

        ni = 'z1' if i == 2 else 'z2'
        nj = 'z1' if j == 2 else 'z2'
        print(f"\n  [{ni},{nj}]:")
        print(f"    ℛ = {R_val}")

        # Try: ℛ_{ab} = (1/2)[FF(G₄)_{ab}/3! − (1/4)|G₄|²g_{ab}] + (1/2)|∂Ψ|²g_{ab}
        # where |∂Ψ|² = H'²/(4H^{9/4}) as for F1 (since Tr is SL(2,Z) invariant)
        dPsi_sq = hf.Hp**2 / (4 * hf.H**R(9, 4))

        T_form = sp.cancel(R(1, 2) * (ff_val / 6 - lam * S_pq_val * g_val))
        T_scalar = sp.cancel(R(1, 2) * dPsi_sq * g_val)
        T_total = sp.cancel(T_form + T_scalar)
        diff = sp.cancel(R_val - T_total)
        print(f"    T^G₄(λ=1/4) + (1/2)|∂Ψ|²g = {T_total}")
        print(f"    diff = {diff}")
        if diff != 0:
            torus_pass = False

        # Alternative: use doublet contraction
        # For torus, F₃ has no z-legs, so (F₃^a·F₃^b)_{z,z} = 0
        # Only the trace term contributes:
        # T^{doublet}_{ab} = (1/4)(c^TMc)(0 - (1/4)|H₃|²g_{ab})
        #                  = -(1/16)(c^TMc)|H₃|² g_{ab}
        cTMc_val = hf.substitute(sp.cancel(cTMc))
        T_dbl_torus = sp.cancel(-R(1, 16) * cTMc_val * S_H3_val * g_val)
        resid_dbl = sp.cancel(R_val - T_dbl_torus)
        print(f"    T^doublet(trace only) = {T_dbl_torus}")
        print(f"    ℛ − T^doublet = {resid_dbl}")

        # The residual should be the torus scalar kinetic from KK
        # Try: residual = -(1/4)(∂M·∂M^{-1})_{ab} contracted over y-directions
        # For ab torus indices: this is the KK formula for internal Ricci
        # ℛ_{ab} = -(1/2) M_{ac} □_{10d} M^{-1}_{cb} + ...

print(f"\n  Torus with exp15 formula: {'★ PASS ★' if torus_pass else 'FAIL'}")

# ===================================================================
# Part G: Cross-check — off-diagonal 10d×torus Ricci
# ===================================================================
print("\n" + "=" * 70)
print("PART G: Off-diagonal components (10d × torus)")
print("=" * 70)

print("\nOff-diagonal components ℛ_{m,a} (should be zero):")
off_diag_ok = True
for m in [0, 1, 4]:
    for a in [2, 3]:
        R_val = sp.cancel(Ric_pq[m, a])
        name_m = str(coords[m])
        name_a = str(coords[a])
        if R_val != 0:
            off_diag_ok = False
            print(f"  ℛ[{name_m},{name_a}] = {R_val} ✗")
        else:
            print(f"  ℛ[{name_m},{name_a}] = 0 ✓")

print(f"\n  Off-diagonal: {'★ PASS ★' if off_diag_ok else 'FAIL'}")

# ===================================================================
# Part H: Relation between T^{G₄} and doublet T — derive general formula
# ===================================================================
print("\n" + "=" * 70)
print("PART H: Rewrite doublet equation in terms of G₄")
print("=" * 70)

print("""
For DIAGONAL M (F1/D1), T^{G₄}(λ=1/4) directly gives the correct answer
because raising z-indices with M^{-1} gives a factor M^{-1}_{aa}, and
the 10d doublet contraction uses M_{aa}, and for a single charge vector:
  M^{-1}_{aa} × 1/2 = M_{aa} × (something)
These happen to match because g_{z1}·g_{z2} = det(M) = 1.

For NON-DIAGONAL M, the relationship between M^{-1}_{ab} (used in G₄
contractions) and M_{ab} (used in doublet IIB equation) is:
  M_{ab} = (cofactor of M^{-1})_{ab}
Since det(M) = 1: M_{ab} = adj(M^{-1})_{ab}.

Can we express the doublet equation purely in terms of G₄?
""")

# The doublet equation for 10d directions is:
# ℛ_{mn} = (1/4)(c^T M c)[(H₃²)_{mn}/2 − (1/4)|H₃|²g]
#
# And (G₄²)_{mn}/3! = (c^T M^{-1} c)(H₃²)_{mn}/2
#     |G₄|² = (c^T M^{-1} c)|H₃|²
#
# So: (H₃²)/2 = (G₄²)_{mn}/(3! · c^T M^{-1} c)
#     |H₃|² = |G₄|²/(c^T M^{-1} c)

# Compute c^T M^{-1} c
cTMinvc = sp.cancel(M_inv_pq[0, 0] - 2 * M_inv_pq[0, 1] + M_inv_pq[1, 1])
print(f"c^T M^{{-1}} c = {hf.substitute(sp.cancel(cTMinvc))}")
print(f"c^T M c = {hf.substitute(sp.cancel(cTMc))}")

# Product
product = sp.cancel(cTMc * cTMinvc)
print(f"(c^T M c)(c^T M^{{-1}} c) = {hf.substitute(sp.cancel(product))}")

# For F1: c^T M₀^{-1} c = M₀^{-1}_{22} = H^{1/2}, c^T M₀ c = H^{-1/2}
# Product = 1.
cTM0invc = sp.cancel(M0.inv()[1, 1])
print(f"\nF1: c^T M₀^{{-1}} c = {sp.cancel(cTM0invc)}")
print(f"    c^T M₀ c = {sp.cancel(cTM0c)}")
print(f"    product = {sp.cancel(cTM0c * cTM0invc)}")

# So for general M:
# ℛ_{mn} = (1/4)(c^T M c) / (c^T M^{-1} c) × [(G₄²)_{mn}/3! − (1/4)|G₄|²g]
# = (c^T M c) / (2 c^T M^{-1} c) × T^{G₄}(λ=1/4)

ratio = sp.cancel(cTMc / (2 * cTMinvc))
print(f"\nRescaling factor (c^TMc)/(2c^TM^{{-1}}c) = {hf.substitute(sp.cancel(ratio))}")

# Verify this rescaling works:
print("\n--- Verifying: ℛ = ratio × T^{G₄}(λ=1/4) for 10d directions ---\n")
ratio_val = hf.substitute(sp.cancel(ratio))

all_rescaled_pass = True
for i, name in [(0, 't'), (1, 'x1'), (4, 'y0'), (5, 'y1')]:
    R_val = hf.substitute(sp.cancel(Ric_pq[i, i]))
    ff_val = hf.substitute(sp.cancel(FF_pq[i, i]))
    g_val = hf.substitute(sp.cancel(metric_pq.matrix[i, i]))

    T_G4 = sp.cancel(R(1, 2) * (ff_val / 6 - lam * S_pq_val * g_val))
    T_rescaled = sp.cancel(ratio_val * T_G4)
    diff = sp.cancel(R_val - T_rescaled)

    status = '✓' if diff == 0 else f'✗ diff={diff}'
    if diff != 0:
        all_rescaled_pass = False
    print(f"  [{name}]: {status}")

# y-coords
for k in range(2, 8):
    i = 4 + k
    R_val = hf.substitute(sp.cancel(Ric_pq[i, i]))
    ff_val = hf.substitute(sp.cancel(FF_pq[i, i]))
    g_val = hf.substitute(sp.cancel(metric_pq.matrix[i, i]))
    T_G4 = sp.cancel(R(1, 2) * (ff_val / 6 - lam * S_pq_val * g_val))
    T_rescaled = sp.cancel(ratio_val * T_G4)
    diff = sp.cancel(R_val - T_rescaled)
    if diff != 0:
        all_rescaled_pass = False
        print(f"  [y{k}]: ✗")

print(f"\n  Rescaled T^G₄: {'★ ALL PASS ★' if all_rescaled_pass else 'FAIL'}")

# ===================================================================
# Part I: ROOT CAUSE — Einstein frame changes under SL(2,Z)
# ===================================================================
print("\n" + "=" * 70)
print("PART I: ROOT CAUSE — Einstein frame changes under SL(2,Z)!")
print("=" * 70)

print("""
★ KEY INSIGHT: Under SL(2,Z), the dilaton changes. Since the 10d Einstein
frame metric is g^E = e^{-Φ/2} g^S (where g^S is string frame), the
10d Einstein frame metric is NOT SL(2,Z) invariant.

For the (1,1)-string obtained from F1 by Λ = [[1,0],[1,1]]:
  τ → τ' = τ/(τ+1)  with τ(F1) = iH^{1/2}
  → τ' = iH^{1/2}/(1+iH^{1/2}) → τ₂' = H^{1/2}/(1+H)

  e^{-Φ'} = τ₂' = H^{1/2}/(1+H)  ≠  H^{1/2} = e^{-Φ}(F1)

So the (1,1) Einstein frame has an extra factor (1+H)^{-1/2}.
The 12d metric I computed above uses the WRONG 10d block.
""")

# The 12d metric ds²₁₂ = g^E + M dz² uses Einstein frame for the 10d part.
# This metric is NOT manifestly SL(2,Z) covariant because g^E changes.
# The exp15 "SL(2,Z) covariance" was verified only for F1 ↔ D1 (S-duality),
# which swaps z₁↔z₂ AND flips the dilaton sign, preserving the Einstein frame.

# STRING FRAME is SL(2,Z) invariant: g^S = H^{-1}ds²_{wv} + ds²_8
# In string frame, the 12d metric ds²₁₂ = g^S + M dz² IS covariant.

print("Building string-frame 12d F1 metric for comparison...")

g_sf = sp.zeros(D, D)
g_sf[0, 0] = -1 / H
g_sf[1, 1] = 1 / H
g_sf[2, 2] = M0[0, 0]  # H^{1/2}
g_sf[3, 3] = M0[1, 1]  # H^{-1/2}
for k in range(8):
    g_sf[4 + k, 4 + k] = sp.Integer(1)

metric_sf_f1 = Metric(g_sf, coords)
print(f"  Diagonal: {metric_sf_f1.is_diagonal}")
print("Computing string-frame Ricci (F1)...")
Ric_sf_f1 = metric_sf_f1.ricci_tensor(simplify_func=sp.cancel)
print("Done.")

# G₄ for F1 in string frame
C_f1_sf = FormField(rank=3, dim=12)
C_f1_sf[(0, 1, 3)] = 1 / H
F_f1_sf = exterior_derivative(C_f1_sf, coords)
FF_f1_sf = form_contraction(F_f1_sf, metric_sf_f1)
S_f1_sf = form_norm_squared(F_f1_sf, metric_sf_f1)
S_f1_sf_val = hf.substitute(sp.cancel(S_f1_sf))

print(f"\n  |G₄|²(string frame, F1) = {S_f1_sf_val}")

# Check if ℛ^{SF} = T^{G₄}(λ) works for some λ
print("\nString-frame F1: ℛ/g per block:")
for i, name in [(0, 't'), (2, 'z1'), (3, 'z2'), (4, 'y0')]:
    R_val = hf.substitute(sp.cancel(Ric_sf_f1[i, i]))
    g_val = hf.substitute(sp.cancel(metric_sf_f1.matrix[i, i]))
    R_over_g = sp.cancel(R_val / g_val)
    print(f"  [{name}]: ℛ/g = {R_over_g}")

# Also compute (1,1) string-frame 12d metric
print("\nBuilding string-frame 12d (1,1)-string metric...")
g_sf_pq = sp.zeros(D, D)
g_sf_pq[0, 0] = -1 / H
g_sf_pq[1, 1] = 1 / H
g_sf_pq[2, 2] = M_pq[0, 0]
g_sf_pq[2, 3] = M_pq[0, 1]
g_sf_pq[3, 2] = M_pq[1, 0]
g_sf_pq[3, 3] = M_pq[1, 1]
for k in range(8):
    g_sf_pq[4 + k, 4 + k] = sp.Integer(1)

metric_sf_pq = Metric(g_sf_pq, coords)
print(f"  Diagonal: {metric_sf_pq.is_diagonal}")
print("Computing string-frame Ricci (1,1)...")
Ric_sf_pq = metric_sf_pq.ricci_tensor(simplify_func=sp.cancel)
print("Done.")

# G₄ in string frame for (1,1): same G₄ form, different metric
FF_pq_sf = form_contraction(F_pq, metric_sf_pq)
S_pq_sf = form_norm_squared(F_pq, metric_sf_pq)
S_pq_sf_val = hf.substitute(sp.cancel(S_pq_sf))

print(f"\n  |G₄|²(string frame, (1,1)) = {S_pq_sf_val}")
print(f"  |G₄|²(string frame, F1)    = {S_f1_sf_val}")
print(f"  Match: {sp.cancel(S_pq_sf_val - S_f1_sf_val) == 0}")

# Check ℛ^{SF}(1,1) vs ℛ^{SF}(F1) for 10d directions
print("\n--- String frame: ℛ(1,1) vs ℛ(F1) for 10d directions ---\n")
sf_match = True
for i, name in [(0, 't'), (1, 'x1'), (4, 'y0')]:
    R_f1 = hf.substitute(sp.cancel(Ric_sf_f1[i, i]))
    R_pq = hf.substitute(sp.cancel(Ric_sf_pq[i, i]))
    diff = sp.cancel(R_f1 - R_pq)
    status = '✓' if diff == 0 else f'✗ diff={diff}'
    if diff != 0:
        sf_match = False
    print(f"  [{name}]: ℛ_F1={R_f1}, ℛ_(1,1)={R_pq}, match={status}")

print(f"\n  10d Ricci (string frame) same for F1 and (1,1): {'✓' if sf_match else '✗'}")

# Now check if the G₄ equation holds in string frame
print("\n--- String frame: ℛ = T^{G₄}(λ) for F1? ---\n")

# Try λ = 1/4 first
for test_lam, lam_name in [(R(1, 4), "1/4"), (R(3, 10), "3/10"), (R(1, 2), "1/2")]:
    all_ok = True
    for i in [0, 1, 2, 3, 4]:
        R_val = hf.substitute(sp.cancel(Ric_sf_f1[i, i]))
        ff_val = hf.substitute(sp.cancel(FF_f1_sf[i, i]))
        g_val = hf.substitute(sp.cancel(metric_sf_f1.matrix[i, i]))
        T = sp.cancel(R(1, 2) * (ff_val / 6 - test_lam * S_f1_sf_val * g_val))
        diff = sp.cancel(R_val - T)
        if diff != 0:
            all_ok = False
    print(f"  λ={lam_name}: {'PASS' if all_ok else 'FAIL'}")

# ===================================================================
# Summary
# ===================================================================
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print(f"""
★★★ EXP 16: (p,q)-string uplift — CRITICAL FINDING ★★★

1. The 12d metric ds² = g^{{Einstein}}_{{mn}}dx^m dx^n + M_{{ab}}dz^a dz^b
   is NOT manifestly SL(2,Z) covariant because g^E changes under SL(2,Z)
   (the dilaton changes: e^{{-Φ'}} = H^{{1/2}}/(1+H) for (1,1)-string).

2. The exp15 "SL(2,Z) covariance" was verified only for S-duality (F1↔D1),
   which is a Z₂ subgroup swapping z₁↔z₂. This works because S-duality
   flips the dilaton sign AND swaps torus coordinates, preserving g^E.

3. For a general (p,q)-string, the correct 12d metric in Einstein frame has
   DIFFERENT 10d warp factors. The equation T^{{G₄}}(λ=1/4) does NOT hold.

4. SL(2,Z) INVARIANTS confirmed:
   - Tr(dM/dH·dM^{{-1}}/dH) is SL(2,Z) invariant ✓
   - c^T M c (charge bilinear) is SL(2,Z) invariant ✓
   - The 10d STRING frame metric is SL(2,Z) invariant ✓

5. A truly covariant 12d formulation should use STRING frame for the 10d
   block, where g^S = H^{{-1}}ds²_{{wv}} + ds²_8 is invariant.
   The string-frame equation will involve additional dilaton factors.

6. 10d Ricci in string frame: F1 vs (1,1) same? {sf_match}
""")
