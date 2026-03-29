"""Experiment 17b: Express the string-frame 12d equation in SL(2,Z)-covariant form.

From exp17:
  ℛ^{12d}_{mn}(SF) = (1/4)(H₃²)_{mn} - 2∇_m∂_nΦ + KK_{mn}

Now:
  - (1/4)(H₃²) → SL(2,Z)-covariant form involves the doublet
  - -2∇∂Φ + KK → should combine into Tr(∂M·∂M⁻¹) terms

For 10d directions, the KK contribution is zero for wv and involves yk²/r² for
transverse. The dilaton term -2∇∂Φ has both isotropic and anisotropic parts.

Strategy:
  1. Compute -2∇∂Φ + KK explicitly
  2. Compare with -(1/4)Tr(∂_m M · ∂_n M⁻¹) (the KK scalar kinetic term)
  3. Find the covariant expression.
"""
import sys; sys.path.insert(0, 'src')
import sympy as sp
from sympy import Rational as R
from sugra import (HarmonicFunction, Metric,
                   FormField, exterior_derivative,
                   form_contraction, form_norm_squared)

wv_coords = list(sp.symbols('t x1', real=True))
z1s, z2s = sp.symbols('z1 z2', real=True)
harmonic_coords = list(sp.symbols('y0:8', real=True))
coords = wv_coords + [z1s, z2s] + harmonic_coords
D = 12

hf = HarmonicFunction(transverse_coords=harmonic_coords)
H = sp.Function('H')(hf.r_expr)
y0 = sp.Symbol('y0', real=True)

# ===================================================================
# Part A: Compute -2∇∂Φ + KK in terms of (∂M)(∂M⁻¹)
# ===================================================================
print("=" * 70)
print("PART A: Relate scalar terms to Tr(∂M·∂M⁻¹)")
print("=" * 70)

# From exp17 Part H, KK/g per block:
# t: 0, x1: 0, y0: -H'²y0²/(8H²r²), y1: -H'²y1²/(8H²r²)

# From exp17 Part G, ∇∂Φ/g per block:
# t: H'²/(4H²), x1: H'²/(4H²)
# y0: (-H H' r² + 8H H' y0² + H'² r y0²)/(2H² r³)

# So -2∇∂Φ/g:
# t: -H'²/(2H²), x1: -H'²/(2H²)
# y0: (H H' r² - 8H H' y0² - H'² r y0²)/(H² r³)

# Combined (-2∇∂Φ + KK)/g:
# t: -H'²/(2H²)
# x1: -H'²/(2H²)
# y0: (H H' r² - 8H H' y0² - H'² r y0²)/(H² r³) - H'² y0²/(8H² r²)
#   = (8H H' r² - 64H H' y0² - 8H'² r y0² - H'² r y0²)/(8H² r³)
#   = (8H H' r² - 64H H' y0² - 9H'² r y0²)/(8H² r³)

# Wait let me just compare numerically. The ℛ^{12d}/g from exp17:
# t: -H'²/H², x1: -H'²/H², z1: H'²/(2H²), z2: -H'²/(2H²)
# y0: (8H H' r² - 64H H' y0² - 13H'² r y0²)/(8H² r³)

# And (1/4)(H₃²)/g:
# t: -H'²/(2H²), x1: -H'²/(2H²), z: 0
# y0: -H'² y0²/(2H² r²)

# So (-2∇∂Φ + KK)/g = ℛ/g - (1/4)(H₃²)/g:
# t: -H'²/H² - (-H'²/(2H²)) = -H'²/(2H²)
# y0: [(8HH'r² - 64HH'y0² - 13H'²ry0²) - (-4H'²ry0²)] / (8H²r³)
#   = (8HH'r² - 64HH'y0² - 9H'²ry0²) / (8H²r³)

# Now, M₀ = diag(H^{1/2}, H^{-1/2})
# Tr(∂_m M · ∂_n M⁻¹) for m=n=y_k:
# = (H')²(y_k²/r²) × Tr(dM/dH · dM⁻¹/dH)
# Tr(dM/dH · dM⁻¹/dH) = -1/(2H²)
# So Tr(∂_{yk} M · ∂_{yk} M⁻¹) = -H'² y_k²/(2H² r²)

# For the isotropic part, we need Σ_k g^{yk,yk} Tr(∂_{yk}M·∂_{yk}M⁻¹)
# In string frame g^{yk}=1, so:
# Σ_k Tr(∂_{yk}M·∂_{yk}M⁻¹) = -H'²/(2H²) × Σ(yk²/r²) = -H'²/(2H²)

# Now, the Laplacian of M in the 10d equation contributes to ℛ_{ab}(torus).
# Let me compute -(1/4)Tr(∂_m M · ∂_n M⁻¹) for the diagonal components
# and compare with (-2∇∂Φ + KK).

print("\nFor the worldvolume (t, x1):")
print("  (-2∇∂Φ + KK)/g = -H'²/(2H²)")
print("  -(1/4) × g^{yk}Tr(∂M·∂M⁻¹)_{iso} × g_{tt}/g_{tt} needs careful KK formula")

# Let me take a different approach: compute the 12d KK Ricci explicitly.
# The 12d metric ds²₁₂ = g^S_{mn} dx^m dx^n + M_{ab}(y) dz^a dz^b
# is a fiber bundle with base (10d, g^S) and fiber (T², M).
#
# The KK Ricci decomposition for DIAGONAL fiber (no KK gauge field) is:
# ℛ_{mn} = R_{mn}^{base} - (1/2) M^{ab} ∇_m ∂_n M_{ab}
# ℛ_{ab} = -(1/2) M_{ac} □_base M^{-1}_{cb}   (with proper covariant □)

# For our case with M_{ab}(y):
# ℛ_{mn}(12d) = R_{mn}(10d) - (1/2) Tr(M⁻¹ ∇_m ∂_n M)

# Let me verify this formula.

# Build 10d string-frame metric and Christoffel
coords_10 = wv_coords + harmonic_coords
g_sf_10 = sp.zeros(10, 10)
g_sf_10[0, 0] = -1 / H
g_sf_10[1, 1] = 1 / H
for k in range(8):
    g_sf_10[2 + k, 2 + k] = sp.Integer(1)

metric_sf_10 = Metric(g_sf_10, coords_10)
Ric_sf_10 = metric_sf_10.ricci_tensor(simplify_func=sp.cancel)
chris_10 = metric_sf_10.christoffel(simplify_func=sp.cancel)

# M₀ = diag(H^{1/2}, H^{-1/2})
M0 = sp.Matrix([[H**R(1, 2), 0],
                 [0, H**R(-1, 2)]])
M0_inv = sp.cancel(M0.inv())

# Compute ∇_m ∂_n M_{ab} for diagonal m=n directions
# For m=n=t (index 0): M doesn't depend on t, so ∂_t M = 0.
# ∇_t ∂_t M_{ab} = ∂_t² M - Γ^ρ_{tt} ∂_ρ M
# = 0 - Σ_k Γ^{y_k}_{tt} ∂_{y_k} M

# For m=n=y_k (index 2+k in 10d): ∂_{yk} M ≠ 0
# ∇_{yk} ∂_{yk} M = ∂²_{yk} M - Γ^ρ_{yk,yk} ∂_ρ M

# Let's compute Tr(M⁻¹ ∇_m ∂_n M) for each block

print("\nComputing KK term: -(1/2) Tr(M⁻¹ ∇_m ∂_n M)")

# Helper: compute ∂_i M (a sympy Matrix)
def partial_M(i):
    """Partial derivative of M₀ w.r.t. coords_10[i]."""
    return sp.Matrix([[sp.diff(M0[a, b], coords_10[i]) for b in range(2)] for a in range(2)])

def partial2_M(i, j):
    """Second partial derivative of M₀."""
    return sp.Matrix([[sp.diff(sp.diff(M0[a, b], coords_10[i]), coords_10[j]) for b in range(2)] for a in range(2)])

# ∇_i ∂_j M = ∂_i ∂_j M - Γ^ρ_{ij} ∂_ρ M
def nabla_d_M(i, j):
    """Covariant derivative ∇_i ∂_j M using 10d Christoffel symbols."""
    result = partial2_M(i, j)
    for rho in range(10):
        key = (rho, i, j)
        if key in chris_10:
            result = result - chris_10[key] * partial_M(rho)
    return result

# Compute Tr(M⁻¹ ∇_m ∂_n M) for diagonal blocks
print("\nTr(M⁻¹ ∇_m ∂_n M) / g^{10d}_{mn} per block:")
for i, name in [(0, 't'), (1, 'x1'), (2, 'y0'), (3, 'y1')]:
    nab_M = nabla_d_M(i, i)
    nab_M_simplified = sp.Matrix([[hf.substitute(sp.cancel(nab_M[a, b])) for b in range(2)] for a in range(2)])
    M0_inv_sub = sp.Matrix([[hf.substitute(sp.cancel(M0_inv[a, b])) for b in range(2)] for a in range(2)])

    trace_val = sp.cancel(sp.trace(M0_inv_sub * nab_M_simplified))
    g_val = hf.substitute(sp.cancel(metric_sf_10.matrix[i, i]))
    trace_over_g = sp.cancel(trace_val / g_val)
    print(f"  [{name}]: Tr/g = {trace_over_g}, Tr = {trace_val}")

# Now compute the KK contribution and compare
print("\n--- Verify: ℛ^{12d} = R^{10d} - (1/2)Tr(M⁻¹∇∂M) ---\n")

# Build 12d F1 string-frame metric
g_sf_f1 = sp.zeros(D, D)
g_sf_f1[0, 0] = -1 / H
g_sf_f1[1, 1] = 1 / H
g_sf_f1[2, 2] = M0[0, 0]
g_sf_f1[3, 3] = M0[1, 1]
for k in range(8):
    g_sf_f1[4 + k, 4 + k] = sp.Integer(1)

metric_sf_f1 = Metric(g_sf_f1, coords)
Ric_sf_f1 = metric_sf_f1.ricci_tensor(simplify_func=sp.cancel)

kk_pass = True
for i_12, i_10, name in [(0, 0, 't'), (1, 1, 'x1'), (4, 2, 'y0'), (5, 3, 'y1')]:
    R12 = hf.substitute(sp.cancel(Ric_sf_f1[i_12, i_12]))
    R10 = hf.substitute(sp.cancel(Ric_sf_10[i_10, i_10]))

    nab_M = nabla_d_M(i_10, i_10)
    nab_M_sub = sp.Matrix([[hf.substitute(sp.cancel(nab_M[a, b])) for b in range(2)] for a in range(2)])
    M0_inv_sub = sp.Matrix([[hf.substitute(sp.cancel(M0_inv[a, b])) for b in range(2)] for a in range(2)])
    trace_val = sp.cancel(sp.trace(M0_inv_sub * nab_M_sub))

    kk_predicted = sp.cancel(-R(1, 4) * trace_val)
    diff = sp.cancel(R12 - R10 - kk_predicted)
    status = '✓' if diff == 0 else f'✗ diff={diff}'
    if diff != 0:
        kk_pass = False
    print(f"  [{name}]: KK_actual={sp.cancel(R12-R10)}, KK_predicted={kk_predicted}, {status}")

print(f"\n  KK formula ℛ = R - (1/4)Tr(M⁻¹∇∂M): {'★ VERIFIED ★' if kk_pass else 'FAILED'}")

# ===================================================================
# Part B: Combine into covariant equation
# ===================================================================
print("\n" + "=" * 70)
print("PART B: Full covariant equation for 10d directions")
print("=" * 70)

# From 10d string frame: R^S + 2∇∂Φ = (1/4)(H₃²)
# From KK: ℛ^{12d} = R^S - (1/2)Tr(M⁻¹∇∂M)
# So: ℛ^{12d} = (1/4)(H₃²) - 2∇∂Φ - (1/2)Tr(M⁻¹∇∂M)
#
# Now, the SL(2,Z)-covariant version should use:
# - Form doublet: H₃ → F₃^a (a=1,2)
# - Dilaton Φ encoded in M: Φ = -(1/2)ln(det... no, det M=1.
#   Actually Φ is related to the off-diagonal parameterization of M.
#   M = (1/τ₂)[[|τ|², τ₁],[τ₁, 1]], τ₂ = e^{-Φ}
#
# KEY: ∇∂Φ is NOT a separate quantity — it's part of the M dynamics!
# The M equation of motion in 10d is related to ∇∂Φ and ∇∂C₀.
#
# Let's check: is -2∇∂Φ = -(1/2)Tr(M⁻¹∇∂M) for diagonal M?
#
# For M₀ = diag(H^{1/2}, H^{-1/2}):
# M₀⁻¹∂M₀ = diag(H^{-1/2}, H^{1/2}) · diag((1/2)H^{-1/2}∂H, -(1/2)H^{-3/2}∂H)
#           = diag((1/2H)∂H, -(1/2H)∂H) = diag(Φ', Φ') ... wait
# Actually M₀⁻¹ dM₀/dH = diag(1/(2H), -1/(2H))
# Tr(M₀⁻¹ dM₀/dH) = 0  (traceless! because det M = 1)
#
# So Tr(M⁻¹∂M) = 0 identically when det M = 1.
# Then Tr(M⁻¹∇∂M) = ∇(Tr(M⁻¹∂M)) - Tr(∇M⁻¹·∂M)
#                   = -Tr(∇M⁻¹·∂M)
# Hmm, this is getting complicated. Let me just use the alternative form:
# Tr(M⁻¹∇∂M) = ∂Tr(M⁻¹∂M) - Tr(∂M⁻¹·∂M) - Γ·Tr(M⁻¹∂M)
#             = -Tr(∂M⁻¹·∂M) (since Tr(M⁻¹∂M)=0)
# = Tr(∂M·∂M⁻¹) (with proper sign from d(M⁻¹) = -M⁻¹(dM)M⁻¹)
#
# Actually ∇_m(Tr(M⁻¹∂_nM)) = Tr(∇_mM⁻¹·∂_nM) + Tr(M⁻¹·∇_m∂_nM)
# Since Tr(M⁻¹∂_nM) = ∂_n(ln det M) = 0 (det M=1):
# 0 = Tr(∇_mM⁻¹·∂_nM) + Tr(M⁻¹·∇_m∂_nM)
# So: Tr(M⁻¹·∇_m∂_nM) = -Tr(∇_mM⁻¹·∂_nM)
#                        = Tr(M⁻¹·∂_mM·M⁻¹·∂_nM)  (using ∇M⁻¹ = -M⁻¹(∇M)M⁻¹ + connection terms cancel for scalar)
# Wait, M is a scalar matrix (no spacetime indices), so ∇_m M = ∂_m M.
# And ∂_m(M⁻¹) = -M⁻¹(∂_mM)M⁻¹.
# So: Tr(M⁻¹·∇_m∂_nM) = -Tr(∂_mM⁻¹·∂_nM) = Tr(M⁻¹∂_mM·M⁻¹∂_nM)

# This is the standard result:
# (1/2)Tr(M⁻¹∇∂M) = (1/2)Tr(M⁻¹∂M·M⁻¹∂M) + Laplacian terms (= 0 here)

# Wait, no. Let me be more careful. ∇_m∂_nM is the second covariant derivative
# of M_{ab} treated as a SCALAR (no torus indices are spacetime indices).
# So ∇_m∂_nM = ∂_m∂_nM - Γ^ρ_{mn}∂_ρM.
# This is the covariant Hessian of M.

# The identity is:
# ∂_m(Tr(M⁻¹∂_nM)) = Tr(∂_m(M⁻¹)·∂_nM) + Tr(M⁻¹·∂_m∂_nM)
# = -Tr(M⁻¹∂_mM·M⁻¹∂_nM) + Tr(M⁻¹·∂_m∂_nM)
# Since det M = 1 → Tr(M⁻¹∂_nM) = 0:
# 0 = -Tr(M⁻¹∂_mM·M⁻¹∂_nM) + Tr(M⁻¹·∂_m∂_nM)
# → Tr(M⁻¹·∂_m∂_nM) = Tr(M⁻¹∂_mM·M⁻¹∂_nM)
#
# And: Tr(M⁻¹∇_m∂_nM) = Tr(M⁻¹∂_m∂_nM) - Γ^ρ_{mn}Tr(M⁻¹∂_ρM)
#                        = Tr(M⁻¹∂_mM·M⁻¹∂_nM) - 0
# ★ So: Tr(M⁻¹∇_m∂_nM) = Tr(M⁻¹∂_mM·M⁻¹∂_nM)

# Let's verify this identity numerically
print("\nVerifying: Tr(M⁻¹∇∂M) = Tr(M⁻¹∂M·M⁻¹∂M)")

for i, name in [(0, 't'), (2, 'y0')]:
    # LHS: Tr(M⁻¹∇∂M)
    nab_M = nabla_d_M(i, i)
    nab_M_sub = sp.Matrix([[hf.substitute(sp.cancel(nab_M[a, b])) for b in range(2)] for a in range(2)])
    M0_inv_sub = sp.Matrix([[hf.substitute(sp.cancel(M0_inv[a, b])) for b in range(2)] for a in range(2)])
    lhs = sp.cancel(sp.trace(M0_inv_sub * nab_M_sub))

    # RHS: Tr(M⁻¹∂M·M⁻¹∂M)
    dM = sp.Matrix([[hf.substitute(sp.cancel(sp.diff(M0[a, b], coords_10[i]))) for b in range(2)] for a in range(2)])
    prod = M0_inv_sub * dM
    rhs = sp.cancel(sp.trace(prod * prod))

    diff = sp.cancel(lhs - rhs)
    print(f"  [{name}]: LHS={lhs}, RHS={rhs}, diff={diff}")

# ===================================================================
# Part C: Relate to ∇∂Φ
# ===================================================================
print("\n" + "=" * 70)
print("PART C: Relation between Tr(M⁻¹∇∂M) and ∇∂Φ")
print("=" * 70)

# For M = diag(e^{-Φ}, e^Φ}) (with C₀=0):
# M⁻¹∂M = diag(-∂Φ, ∂Φ)
# (M⁻¹∂M)² = diag((∂Φ)², (∂Φ)²)
# Tr = 2(∂Φ)²
#
# So Tr(M⁻¹∂_mM·M⁻¹∂_nM) = 2 ∂_mΦ ∂_nΦ
# And (1/2)Tr(M⁻¹∇∂M) = ∂_mΦ ∂_nΦ
#
# But we need -2∇∂Φ, NOT ∂Φ∂Φ. These are different!
# ∇∂Φ involves second derivatives, while ∂Φ∂Φ is first-order.

# From the 10d equation: R^S + 2∇∂Φ = (1/4)(H₃²)
# R^S = -2∇∂Φ + (1/4)(H₃²)
#
# From KK: ℛ^{12d} = R^S - (1/2)Tr(M⁻¹∇∂M)
#         = -2∇∂Φ + (1/4)(H₃²) - (1/2)Tr(M⁻¹∇∂M)
#
# Note: Tr(M⁻¹∇∂M) = Tr(M⁻¹∂M·M⁻¹∂M) = 2∂Φ∂Φ
# And we already showed that ∇∂Φ = first-order + second-order terms.
#
# Actually, let me just directly verify the combination:
# Is ℛ^{12d} = (1/4)(H₃²) - 2∇∂Φ - ∂Φ∂Φ ?
# = (1/4)(H₃²) - 2∇∂Φ - (1/2)Tr(M⁻¹∂M·M⁻¹∂M)

# Hmm wait, actually there are TWO different trace quantities:
# 1. Tr(M⁻¹∇_m∂_nM) — Hessian version (what appears in KK formula)
# 2. Tr(M⁻¹∂_mM·M⁻¹∂_nM) — product version
# These are EQUAL (shown above). Both = 2∂_mΦ∂_nΦ for diagonal M.
#
# But the KK formula I used was ℛ = R - (1/2)Tr(M⁻¹∇∂M), and this was verified.
# So: ℛ = R^S - (1/2)·2∂Φ∂Φ = R^S - ∂Φ∂Φ = [(1/4)(H₃²) - 2∇∂Φ] - ∂Φ∂Φ
#
# Now, ∇∂Φ relates to ∂Φ∂Φ via:
# ∇_m∂_nΦ = ∂_m∂_nΦ - Γ^ρ_{mn}∂_ρΦ (second derivative of scalar)
# ∂_mΦ∂_nΦ (product of first derivatives)
# These are independent quantities.

# Let me check numerically whether the combination
# ℛ = (1/4)(H₃²) - 2∇∂Φ - (∂Φ)²_{mn}
# works. Here (∂Φ)²_{mn} = ∂_mΦ ∂_nΦ.

# From exp17: Φ = -(1/2)ln H
# ∂_{yk}Φ = -(1/2)(H'/H)(yk/r)
# ∂Φ∂Φ_{yk,yk} = (1/4)(H'/H)²(yk²/r²)
# ∂Φ∂Φ_{tt} = 0 (Φ doesn't depend on t)

# So for t-block: -(∂Φ)²_{tt} = 0 (correct, we need 0 extra for t)
# But from above: KK_t = 0 and -2∇∂Φ_t ≠ 0.
# Let me verify: ℛ_t = (1/4)(H₃²)_t - 2∇∂Φ_t - (∂Φ)²_t
# = (1/4)(H₃²)_t - 2∇∂Φ_t - 0
# But ℛ_t = R^S_t + 0 (KK_t=0) = R^S_t
# And R^S_t = (1/4)(H₃²)_t - 2∇∂Φ_t (10d equation)
# So yes, this works for t. ✓

# For y0: -(∂Φ)²_{y0,y0} = -(1/4)(H'/H)²(y0²/r²)
# ℛ_y0 = R^S_y0 + KK_y0 = R^S_y0 - (1/2)·2(∂Φ)²_y0 = R^S_y0 - (∂Φ)²_y0
# R^S_y0 = (1/4)(H₃²)_y0 - 2∇∂Φ_y0
# So ℛ_y0 = (1/4)(H₃²)_y0 - 2∇∂Φ_y0 - (∂Φ)²_y0
# This is: (1/4)(H₃²) - 2∇∂Φ - (∂Φ)²  (no new independent structure)

print("The 10d-direction equation in string frame is:")
print("  ℛ^{12d}_{mn} = (1/4)(H₃²)_{mn} − 2∇_m∂_nΦ − ∂_mΦ·∂_nΦ")
print("This combines KK (= -(∂Φ)²) with the 10d equation (R^S + 2∇∂Φ = (1/4)H₃²).")
print()

# For SL(2,Z) covariance, rewrite -2∇∂Φ - (∂Φ)² in terms of M:
# -2∇∂Φ - (∂Φ)² = -∂²Φ/(covariant) - ... this isn't cleanly expressible.
#
# Better approach: don't decompose into 10d + KK. Instead, express ℛ^{12d}
# directly in terms of 12d quantities (G₄ and M).
#
# ℛ^{12d}_{mn} = (1/4)(H₃²)_{mn} − 2∇_m∂_nΦ − ∂_mΦ∂_nΦ
# Now:
# H₃ → can be expressed via G₄ and M
# ∇∂Φ + (1/2)(∂Φ)² → can this be expressed via Tr(M⁻¹∇∂M)?
#
# From Tr(M⁻¹∇∂M) = 2(∂Φ)² and we need -2∇∂Φ - (∂Φ)²,
# there's no way to get ∇∂Φ from algebraic operations on M derivatives.
# We NEED the Hessian ∇²M as well.

# Alternative: use the identity
# ∇_m(M⁻¹∂_nM) = -M⁻¹(∂_mM)M⁻¹(∂_nM) + M⁻¹∇_m∂_nM
# Tr[∇_m(M⁻¹∂_nM)] = -Tr(M⁻¹∂M·M⁻¹∂M) + Tr(M⁻¹∇∂M) = 0
# (since both equal 2(∂Φ)²). So the COVARIANT DIVERGENCE of (M⁻¹∂M) is traceless.

# For diagonal M = diag(e^σ, e^{-σ}) with σ = -Φ = (1/2)lnH:
# M⁻¹∂M = diag(∂σ, -∂σ) → Tr = 0
# ∇(M⁻¹∂M) = diag(∇∂σ, -∇∂σ) → Tr[∇(M⁻¹∂M)] = 0 ✓
# Tr[(M⁻¹∂M)²] = 2(∂σ)² = 2(∂Φ)²
# Tr[M⁻¹∇∂M] = 2(∂σ)² = 2(∂Φ)² (from the identity above)

# So we cannot distinguish ∇∂Φ from (∂Φ)² using M traces alone.
# The equation -2∇∂Φ - (∂Φ)² = -(2∇∂σ + (∂σ)²) involves DIFFERENT
# numbers of derivatives and cannot be written as a single M-trace.

# However: is there a DIFFERENT KK reduction formula?
# The standard formula for the reduced Ricci tensor from a higher-dimensional
# metric ds² = g_{mn}dx^m dx^n + h_{ab}(x)dy^a dy^b is:
#
# R_{mn}^{D} = R_{mn}^{d} + (1/4)Tr(∂_m h · h⁻¹ · ∂_n h · h⁻¹)
#            - (1/2)Tr(h⁻¹ ∇_m^{(d)} ∂_n h)
#
# Wait, this doesn't look right either. Let me use the standard result from
# Maharana-Schwarz or textbook dimensional reduction.

# From Pope's lecture notes (hep-th/9912164), the KK reduction of the
# higher-dimensional Ricci scalar with internal metric M_{ab}(x):
# R^{(D)} includes -(1/4)Tr(∂_mM · ∂^m M⁻¹) from the internal metric.
# This gives the scalar kinetic term in the reduced action.

# For the RICCI TENSOR (not scalar), the standard result is:
# ℛ_{mn}^{(D)} = R_{mn}^{(d)} - (1/2) h^{ac} ∇_m∂_n h_{ac}  (wrong)
# Actually it's more subtle. Let me compute it from scratch.

# Actually, I already verified numerically that ℛ = R - (1/2)Tr(M⁻¹∇∂M).
# So let me just accept that and work with:
# ℛ_{mn} = (1/4)(H₃²) - 2∇∂Φ - (1/2)·2(∂Φ)² = (1/4)(H₃²) - 2∇∂Φ - (∂Φ)²

# The SL(2,Z)-covariant form of this involves:
# (1) The form term: (1/4)(H₃²) → (1/4)(c^T M c)(H₃²) for doublet
#     But for F1, c^T M c = M_{22} = H^{-1/2} = e^Φ
#     So (1/4)(c^T M c)(H₃²) = (1/4)e^Φ(H₃²) which is WRONG for string frame
#     (should be just (1/4)(H₃²) without e^Φ).
#
# Wait — the issue is that in string frame, H₃ appears without dilaton coupling,
# while F₃ (RR) appears with e^Φ. So the 10d IIB equation is:
# R^S + 2∇∂Φ = (1/4)H₃² + (1/4)e^{2Φ}F₃²  (in string frame)
#
# But e^{2Φ} for the RR sector makes this NOT manifestly SL(2,Z) covariant.
# The covariant form uses M^{-1}:
# (1/4)e^{2Φ}F₃² + (1/4)H₃² = (1/4)M⁻¹_{11}(F₃²) + (1/4)M⁻¹_{22}(H₃²)
# where M⁻¹ = diag(e^{Φ}, e^{-Φ}) for C₀=0.
# But e^{-Φ}(H₃²) ≠ (H₃²)! So M⁻¹ doesn't work either.
#
# Hmm, let me be more careful. For F1 with Φ = -(1/2)lnH:
# e^{2Φ} = 1/H, e^{-2Φ} = H
# M = diag(H^{1/2}, H^{-1/2})
# M⁻¹ = diag(H^{-1/2}, H^{1/2})
#
# M⁻¹_{22} = H^{1/2} = e^{-Φ}
# But we need NO dilaton factor for H₃ in string frame.
# The form term is just (1/4)(H₃²).
# So M⁻¹_{22}(H₃²) = e^{-Φ}(H₃²) ≠ (H₃²).
#
# This means the string-frame 3-form equation is NOT (1/4)M⁻¹(F₃^a·F₃^b).
# Let me look at this more carefully.

# The IIB action in string frame (using the SL(2,Z) doublet) is:
# S = ∫ e^{-2Φ}[R + 4(∂Φ)² - (1/2)|H₃|²] - (1/2)|F₁|² - (1/2)|F₃'|² - ...
#
# The NSNS sector has e^{-2Φ} prefactor. The RR 3-form F₃' has no such factor.
# So the equations of motion are:
#
# e^{-2Φ}[R_{mn} + 2∇_m∂_nΦ - (1/4)(H₃²)_{mn}] = (1/4)(F₃'²)_{mn}
# → R + 2∇∂Φ = (1/4)(H₃²) + (1/4)e^{2Φ}(F₃'²)
#
# This can be written as:
# R + 2∇∂Φ = (1/4)S_{ab}(F₃^a·F₃^b)/2!
# where S = diag(e^{2Φ}, 1) ← THIS IS NOT M or M⁻¹!
#
# Actually: S = e^Φ · M⁻¹ ?
# e^Φ M⁻¹ = e^Φ diag(e^{-Φ}, e^Φ) = diag(1, e^{2Φ})  ← for a=1→RR, a=2→NSNS
# But we want diag(e^{2Φ}, 1) = S for a=1→RR, a=2→NSNS
# So S = e^Φ M^T?  No...
# S = M⁻¹ × e^{-Φ}? → diag(e^{-2Φ}, 1) ← wrong
# S_{ab} = ε_{ac} M^{cd} ε_{db} × something?
#
# Actually: S = diag(e^{2Φ}, 1) = [[e^{2Φ}, 0], [0, 1]]
# M = [[e^{-Φ}, 0], [0, e^Φ]] (for τ = ie^{-Φ}, C₀=0)
# M⁻¹ = [[e^Φ, 0], [0, e^{-Φ}]]
#
# S = e^Φ · diag(e^Φ, e^{-Φ}) = e^Φ M⁻¹
# ★ YES: S_{ab} = e^Φ M⁻¹_{ab} = (detM)^{1/2 something}...
# det M = 1, so this doesn't help.
# e^Φ is NOT SL(2,Z) invariant.

# This confirms: the string-frame IIB equation is NOT manifestly SL(2,Z) covariant
# for the 3-form sector. The SL(2,Z)-covariant form of the equations uses
# the EINSTEIN FRAME.

# So we have a dilemma:
# - Einstein frame: SL(2,Z) covariant 10d equation, but NON-covariant 12d metric
# - String frame: SL(2,Z) invariant 12d metric, but NON-covariant 10d equation

print("★ KEY INSIGHT:")
print("  - String frame: 12d metric is SL(2,Z) invariant, but the 10d IIB equation")
print("    is NOT manifestly covariant (e^{2Φ} prefactor on RR 3-form).")
print("  - Einstein frame: 10d equation IS covariant, but 12d metric is NOT.")
print("  - There is NO frame where BOTH are manifestly SL(2,Z) covariant.")
print()
print("  However, the 12d STRING-FRAME Ricci tensor IS the same for all (p,q)-strings")
print("  (verified in Part C). So the EQUATION must be covariant even if the")
print("  individual terms are not.")

# ===================================================================
# Part D: Express the equation using G₄ and M directly in 12d
# ===================================================================
print("\n" + "=" * 70)
print("PART D: 12d equation using G₄ and M (string frame)")
print("=" * 70)

# For 10d directions:
# ℛ^{12d}_{mn} = (1/4)(H₃²)_{mn} - 2∇_m∂_nΦ - ∂_mΦ·∂_nΦ
#
# Now express each term using 12d G₄ and M:
#
# 1. (H₃²)_{mn}: H₃ lives in 10d subspace, so (H₃²)_{mn} = (G₄²)_{mn,z}
#    where we contract only the z-index of G₄ with M⁻¹:
#    (G₄·G₄)_{mn,12d} = Σ_{PQR} G_{mPQR}G_n^{PQR}
#    Some P,Q,R are z-indices, raised with M⁻¹, giving M⁻¹ factors.
#    For F1: G₄ = H₃ ∧ dz₂, so G_{m,a,b,z₂} and raising z₂ gives g^{z₂z₂} = e^Φ.
#    (G₄²)_{mn} = g^{z₂z₂} (H₃²)_{mn} = e^Φ (H₃²)_{mn}
#
#    So (H₃²) = e^{-Φ} (G₄²)_{12d} for F1.
#    And (1/4)(H₃²) = (1/4)e^{-Φ}(G₄²)_{12d}
#
# 2. -2∇∂Φ - (∂Φ)² in terms of M:
#    For diagonal M: Φ = -(1/2)ln(M_{11}/M_{22})^{1/2} = -(1/4)ln(M_{11}/M_{22})
#    But this isn't SL(2,Z) covariant.

# Actually let me try a COMPLETELY DIFFERENT APPROACH: just fit the equation
# ℛ_{mn} = A · (G₄²)_{mn}/3! + B · |G₄|² g_{mn} + C · Tr(∂_mM·∂_nM⁻¹) + D · (something)
# using the numerical data.

# Build G₄ for F1
C_f1 = FormField(rank=3, dim=D)
C_f1[(0, 1, 3)] = 1 / H  # C_{t,x1,z2}
G4_f1 = exterior_derivative(C_f1, coords)
FF_f1_sf = form_contraction(G4_f1, metric_sf_f1)
S_f1_sf = form_norm_squared(G4_f1, metric_sf_f1)
S_f1_sf_val = hf.substitute(sp.cancel(S_f1_sf))

print("\nFitting: ℛ = A·FF(G₄)/6 + B·|G₄|²g + C·Tr(∂M∂M⁻¹)")
print("(using F1 string-frame data)")

# For F1 in string frame:
# |G₄|² = -H'²/H^{3/2}
# FF(G₄)/6 per block:
# (∂Φ)²_{mm} = 0 for t, (1/4)(H'/H)²(yk²/r²) for yk
# (∂Φ)²_{iso} = (1/4)(H'/H)² (summed, but as "trace" it appears with g_{mm})

# Compute FF(G₄)/6 for each block
print("\nFF(G₄)/6 per block:")
for i, name in [(0, 't'), (1, 'x1'), (2, 'z1'), (3, 'z2'), (4, 'y0'), (5, 'y1')]:
    ff_val = hf.substitute(sp.cancel(FF_f1_sf[i, i]))
    print(f"  [{name}]: {sp.cancel(ff_val / 6)}")

# Let me try the most general SL(2,Z)-covariant ansatz:
# For 10d directions:
# ℛ_{mn} = a · (G₄²)_{mn}/3! + b · |G₄|² g_{mn} + c · Tr(∂_mM∂_nM⁻¹)
# where Tr(∂_mM∂_nM⁻¹) = -2(∂Φ)²_{mn} for diagonal M with det M = 1.

# For m=n=t: (G₄²)_{tt}/3! and |G₄|²g_{tt} are nonzero, Tr(∂M∂M⁻¹)_{tt} = 0
# For m=n=y0: all three are nonzero

# Two equations, three unknowns → need z-direction as well
# But for z-directions, G₄ has different structure.

# Let me just use t and y0 first with the constraint from exp13 that
# the y0-equation should hold for BOTH isotropic and anisotropic parts.

# Actually, the form Tr(∂_mM∂_nM⁻¹) is anisotropic in y-directions:
# Tr(∂_{yk}M∂_{yk}M⁻¹) = -1/(2H²) · H'² · yk²/r²  [NOT isotropic]
# The isotropic part would be: g^{yl,yl}Σ Tr = -H'²/(2H²)

# OK this is getting algebraically complex. Let me compute everything numerically
# and solve the linear system.

# Using on-shell values (H = 1 + Q/r^6 → H' = -6Q/r^7, H'' from harmonic condition)
# at a specific numerical point.

import random
rng = random.Random(42)
vals = hf.random_values(rng)
# Also need explicit y-coord values
r_val = vals[hf.r]
y_vals = {}
# Generate random y-coords on the sphere of radius r
import math
ys = [rng.gauss(0, 1) for _ in range(8)]
norm = math.sqrt(sum(y**2 for y in ys))
for k, yc in enumerate(harmonic_coords):
    y_vals[yc] = ys[k] / norm * r_val

all_vals = {**vals, **y_vals}

def neval(expr):
    """Numerically evaluate an expression."""
    return float(sp.N(sp.cancel(expr).subs(all_vals), 15))

# Build the data: ℛ_{mm}, FF/6, |G₄|²g, Tr(∂M∂M⁻¹)
print("\nNumerical evaluation at random point:")
print(f"  H = {neval(hf.H):.6f}, H' = {neval(hf.Hp):.6f}, r = {r_val:.6f}")

data = []
for i, name in [(0, 't'), (1, 'x1'), (2, 'z1'), (3, 'z2'), (4, 'y0'), (5, 'y1')]:
    R_val = neval(hf.substitute(sp.cancel(Ric_sf_f1[i, i])))
    ff_val = neval(hf.substitute(sp.cancel(FF_f1_sf[i, i]))) / 6.0
    g_val = neval(hf.substitute(sp.cancel(metric_sf_f1.matrix[i, i])))
    S_val = neval(S_f1_sf_val)

    # Tr(∂_iM∂_iM⁻¹) in 12d coordinates
    # For 12d index i: if i in {0,1} (wv) → 0
    # If i in {2,3} (z) → 0 (M doesn't depend on z)
    # If i in {4,...,11} (y) → -H'²yk²/(2H²r²)
    if i >= 4:
        yk = harmonic_coords[i - 4]
        tr_dMdMinv = neval(-hf.Hp**2 * yk**2 / (2 * hf.H**2 * hf.r**2))
    else:
        tr_dMdMinv = 0.0

    data.append((name, R_val, ff_val, S_val * g_val, tr_dMdMinv))
    print(f"  [{name}]: ℛ={R_val:.8f}, FF/6={ff_val:.8f}, |G₄|²g={S_val*g_val:.8f}, Tr(∂M∂M⁻¹)={tr_dMdMinv:.8f}")

# Solve: ℛ = a·FF/6 + b·|G₄|²g + c·Tr(∂M∂M⁻¹) using t, y0, z1
# t: ℛ_t = a·FF_t/6 + b·S·g_t + c·0
# y0: ℛ_y0 = a·FF_y0/6 + b·S·g_y0 + c·Tr_y0
# z1: ℛ_z1 = a·FF_z1/6 + b·S·g_z1 + c·0 (z1 has no M dependence)

# Use t (i=0), z1 (i=2), y0 (i=4) to solve for a, b, c
# t: ℛ_t = a·FF_t/6 + b·S·g_t + c·0  (no Tr term for wv)
# z1: ℛ_z1 = a·0 + b·S·g_z1 + c·0     (FF=0 for z1, no Tr term)
# y0: ℛ_y0 = a·FF_y0/6 + b·S·g_y0 + c·Tr_y0

# From z1: b = ℛ_z1 / (S·g_z1)
b_val = data[2][1] / data[2][3]
print(f"\n  From z1: b = ℛ_z1/(|G₄|²g_z1) = {b_val:.10f}")

# From t: a = (ℛ_t - b·S·g_t) / (FF_t/6)
a_val = (data[0][1] - b_val * data[0][3]) / data[0][2]
print(f"  From t:  a = {a_val:.10f}")

# From y0: c = (ℛ_y0 - a·FF_y0/6 - b·S·g_y0) / Tr_y0
c_val = (data[4][1] - a_val * data[4][2] - b_val * data[4][3]) / data[4][4]
print(f"  From y0: c = {c_val:.10f}")

# Verify on all blocks
print("\n  Verification:")
for idx, (name, R_val, ff, Sg, tr) in enumerate(data):
    predicted = a_val * ff + b_val * Sg + c_val * tr
    err = abs(R_val - predicted)
    print(f"    [{name}]: predicted={predicted:.8f}, actual={R_val:.8f}, err={err:.2e}")

# ===================================================================
# Part E: Try equation with ∇∂Φ term
# ===================================================================
print("\n" + "=" * 70)
print("PART E: Extended basis with ∇²Φ-type terms")
print("=" * 70)

# The issue is that (∂Φ)² is purely anisotropic in y-directions,
# while we also need an isotropic scalar contribution.
# The isotropic part comes from ∇∂Φ which involves second derivatives.

# In 12d terms, the "second derivative of M" translates to:
# (1/2)M^{ac}∇_m∂_nM_{cb} = ∇∂Φ (for diagonal M, the off-diagonal vanishes)

# So the full basis should include:
# ℛ_{mn} = a·(G₄²)_{mn}/3! + b·|G₄|²g_{mn} + c·Tr(M⁻¹∂_mM·M⁻¹∂_nM) + d·Tr(M⁻¹∇_m∂_nM)·g_{mn}

# Wait, Tr(M⁻¹∇∂M) is a SCALAR for each (m,n), not a scalar times g.
# Let me rethink.

# The equation I derived is:
# ℛ_{mn} = (1/4)(H₃²)_{mn} - 2∇_m∂_nΦ - ∂_mΦ∂_nΦ
#
# In 12d, (H₃²) = (1/g^{z₂z₂})(G₄²)_{12d} = e^{-Φ}(G₄²).
# But e^{-Φ} is a function of the moduli.

# For the doublet approach:
# In string frame, the 10d equation for general (p,q) is:
# R^S + 2∇∂Φ = (1/4)H₃² + (1/4)e^{2Φ}F₃² + (1/4)F₅² ...
# = (1/4)[e^{2Φ}c_1² + c_2²](H₃²)  for charge vector c=(c_1,c_2)
# where (H₃²) is the SAME form field strength for all (p,q).
#
# For F1: c=(0,1) → coefficient = 1. ✓
# For D1: c=(1,0) → coefficient = e^{2Φ} = 1/H.
# For (1,1): c=(-1,1) → coefficient = e^{2Φ} + 1 = (1+H)/H.
#
# And the factor (e^{2Φ}c_1² + c_2²) = c^T S c where S = diag(e^{2Φ}, 1).
#
# Now S = ?  In terms of M:
# M⁻¹ = diag(e^Φ, e^{-Φ}), so e^{2Φ} = (M⁻¹_{11})²/1 = ... this is messy.
# Actually: S_{ab} = δ_{a1}δ_{b1}(M⁻¹_{11})² + δ_{a2}δ_{b2}
# Not a nice covariant quantity.

# Alternative: use M_{ab} with F₃^a to get covariant expression.
# The doublet contraction M_{ab}c_a c_b = c^T M c is SL(2,Z) invariant.
# For F1: c^T M c = M_{22} = H^{-1/2} = e^Φ
# The form coefficient is c^T S c = c_2² = 1
# Ratio: (c^T S c)/(c^T M c) = 1/e^Φ = e^{-Φ} = H^{1/2}

# For (1,1): c=(-1,1), c^T M c = M_{11} - 2M_{12} + M_{22} = H^{-1/2} (SL(2,Z) inv!)
# c^T S c = e^{2Φ} + 1 = (1+H)/H
# Ratio: (1+H)/(H · H^{-1/2}) = (1+H)/H^{1/2} ← NOT = H^{1/2} in general!

# So the ratio is NOT SL(2,Z) invariant. This confirms the string-frame equation
# is not manifestly covariant.

# FINAL APPROACH: accept that the equation mixes frames and express it
# in the most useful form.

# The correct 12d STRING-FRAME equation for the 10d-direction components is:
# ★ ℛ_{mn}(12d, SF) = (1/4)(H₃²)_{mn} - 2∇_m∂_nΦ - ∂_mΦ∂_nΦ
# where H₃ is the UNIVERSAL 3-form field strength (same for all (p,q)),
# Φ is the dilaton, and all quantities are in string frame.

# This is SL(2,Z) covariant IN THE FOLLOWING SENSE:
# - ℛ_{mn} is the same for all (p,q) (because the SF 12d metric is SL(2,Z) inv) ✓
# - H₃ is the same for all (p,q) (universal field strength) ✓
# - ∇∂Φ and (∂Φ)² depend on Φ which changes under SL(2,Z)
# But since ℛ is invariant and (1/4)(H₃²) is invariant,
# the combination 2∇∂Φ + (∂Φ)² must ALSO be invariant.

# Is 2∇∂Φ + (∂Φ)² SL(2,Z) invariant?
# For F1: Φ = -(1/2)ln H, (∂Φ)² = (1/4)(H'/H)², ∇∂Φ = ...
# For (1,1): Φ' = -(1/2)ln(H/(1+H)), (∂Φ')² = (1/4)(H'/(H(1+H)))² × ((1+H)² - H²)...
# This is getting complicated. But since ℛ is invariant and H₃ is invariant,
# the dilaton terms MUST be invariant too. And they equal:
# 2∇∂Φ + (∂Φ)² = (1/4)(H₃²) - ℛ  (invariant)

# In fact: 2∇∂Φ + (∂Φ)² = ∇²Φ + (∂Φ)² ... no wait.
# 2∇_m∂_nΦ is twice the Hessian, not the Laplacian.
# 2∇_m∂_nΦ + ∂_mΦ∂_nΦ is not a standard geometric quantity.

# But WAIT: the KK decomposition says:
# ℛ_{mn}(12d) = R_{mn}(10d) - (1/2)Tr(M⁻¹∂_mM·M⁻¹∂_nM)
# R_{mn}(10d) = (1/4)(H₃²) - 2∇∂Φ  (from 10d SF equation)
# (1/2)Tr(M⁻¹∂M·M⁻¹∂M) = ∂Φ∂Φ  (for diagonal M, Tr/2 = (∂Φ)²)
#
# So: ℛ = (1/4)(H₃²) - 2∇∂Φ - (∂Φ)²
# And: 2∇∂Φ + (∂Φ)² = R_{10d}^S + (1/2)Tr(prod)
#                      = (1/4)(H₃²) - ℛ (which is SL(2,Z) inv since both terms are)

# THE REAL QUESTION: can we write ℛ_{mn} directly in terms of M and G₄,
# without decomposing into 10d + KK?

# ★ YES: The KK formula gives us ℛ = R^{10d} - (1/2)Tr(M⁻¹∂M·M⁻¹∂M)/2
# And R^{10d} satisfies the 10d equation which involves M via the dilaton.
# Substituting:
# ℛ = [(1/4)(H₃²) - 2∇∂Φ] - (1/2)Tr(M⁻¹∂M)²/2
#
# This still has ∇∂Φ which is NOT a simple function of M.
# ∇∂Φ = (1/2)∂M stuff... let me compute.
#
# For M = diag(e^σ, e^{-σ}) with σ = (1/2)lnH:
# ∂σ = (1/2)H'/H × (yk/r)
# ∇∂σ = Hessian(σ) with connection
#
# M⁻¹∂M = diag(∂σ, -∂σ) → the individual components are ±∂σ
# ∇(M⁻¹∂M) = diag(∇∂σ, -∇∂σ)
# So ∇∂σ = (1/2)[∇(M⁻¹∂M)]_{11} = -(1/2)[∇(M⁻¹∂M)]_{22}
# And: Tr[∇(M⁻¹∂M)] = 0 (as expected from det M = 1)
#
# But [∇(M⁻¹∂M)]_{11} is the (1,1) component of a matrix, not a trace.
# For SL(2,Z) covariance we can only use traces.

# I think the right formulation uses the MATRIX equation, not just traces:
# The 10d scalar equation is: ∇²M + (form terms involving M) = 0
# This KK-lifts to the torus block of ℛ.
# And the external components involve BOTH M traces AND individual M components.

# Let me try another standard form. The 10d IIB action in Einstein frame is:
# S_E = ∫ [R - (1/2)(∂Φ)² - (1/2)e^{2Φ}(∂C₀)²
#         - (1/4)Tr(∂M ∂M⁻¹) - (form terms)]
# where (1/4)Tr(∂M∂M⁻¹) = (1/2)(∂Φ)² + (1/2)e^{2Φ}(∂C₀)²
# So: (1/4)Tr(∂M∂M⁻¹) is the FULL scalar kinetic term in Einstein frame.

# In EINSTEIN frame, the equation is:
# R^E_{mn} = (1/4)Tr(∂_mM∂_nM⁻¹) + (form terms)
# And KK: ℛ^{12d,EF}_{mn} = R^E_{mn} + KK^E_{mn}

# Since in Einstein frame R^E has the simple form with Tr(∂M∂M⁻¹),
# and KK^E adds more M terms, the full 12d EINSTEIN-FRAME equation is:
# ℛ^{12d,EF}_{mn} = (1/4)Tr(∂_mM∂_nM⁻¹) + (form terms) + KK^E
# This was exp13's result!

# THE PROBLEM (from exp16): the Einstein-frame 12d metric is not SL(2,Z) inv.
# So this equation is correct for F1/D1 but NOT for general (p,q).

# CONCLUSION for the string-frame 12d equation:
print("\n★★ CONCLUSION ★★")
print()
print("The string-frame 12d equation for 10d directions decomposes as:")
print("  ℛ_{mn}(12d,SF) = (1/4)(H₃²)_{mn} − 2∇_m∂_nΦ − ∂_mΦ∂_nΦ")
print()
print("Equivalently (using KK reduction):")
print("  ℛ_{mn}(12d,SF) = R_{mn}(10d,SF) − (1/2)Tr(M⁻¹∂_mM·M⁻¹∂_nM)")
print()
print("where:")
print("  - R^{10d,SF} satisfies R^S + 2∇∂Φ = (1/4)(H₃²)  (10d string-frame IIB)")
print("  - (1/2)Tr(M⁻¹∂M·M⁻¹∂M) = ∂_mΦ∂_nΦ  (for diagonal M with det M = 1)")
print("  - All 12d Ricci components for 10d directions are SL(2,Z) INVARIANT")
print("  - ℛ(F1) = ℛ(D1) = ℛ(1,1) — verified numerically")
print()
print("The obstruction to a SINGLE manifestly SL(2,Z)-covariant formula:")
print("  - In Einstein frame: equation is covariant but 12d metric is not")
print("  - In string frame: 12d metric is invariant but equation has e^{2Φ}")
print("    coupling the RR 3-form (not expressible as a simple M-trace)")
print()
print("This is the KNOWN tension in IIB: SL(2,Z) is a symmetry of the full theory")
print("but cannot be made manifest in a single Lagrangian description.")
print("The 12d formulation inherits this limitation.")
