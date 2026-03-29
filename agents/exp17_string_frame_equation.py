"""Experiment 17: String-frame 12d field equation — SL(2,Z) covariant.

From exp16: the Einstein-frame 12d metric is NOT SL(2,Z) covariant because
g^E = e^{-Φ/2} g^S changes under SL(2,Z). But the string-frame 10d metric
g^S = H^{-1} ds²_{wv} + ds²_8 IS invariant for all (p,q)-strings.

Strategy:
  1. Build string-frame 12d metric: ds²₁₂ = g^S + M_{ab} dz^a dz^b
     for F1 and (1,1)-string.
  2. Compute Ricci tensors and verify ℛ^S(F1) = ℛ^S(1,1) for 10d directions.
  3. The 10d IIB string-frame equation is:
       R^S_{mn} = -(1/2)∇_m∂_nΦ + (1/4) M_{ab}[F₃^a·F₃^b_{mn}/2! - (1/6)|F₃^a·F₃^b|g^S_{mn}]
     (need to verify the exact form — coefficients may differ from Einstein frame)
  4. Compute the KK decomposition: ℛ^{12d}_{MN} = R^S_{mn} + KK terms
  5. Express the full 12d string-frame equation in SL(2,Z)-covariant form.

The string-frame 10d IIB equations (from Polchinski or standard references):
  R^S_{mn} + 2∇_m∂_nΦ = (1/4)(H₃)²_{mn}   (for F1 with only NSNS 3-form)

More generally with the doublet:
  R^S_{mn} + 2∇_m∂_nΦ = (1/4) e^Φ (F₃·F₃)_{mn} + (1/4)(H₃·H₃)_{mn}
                       = (1/4) M_{ab}^{-1} (F₃^a·F₃^b)_{mn}/2!

Wait — need to be careful. In string frame, the NSNS 3-form H₃ appears with
coefficient 1, while the RR 3-form F₃ has e^Φ. In the doublet notation:
  F₃^1 = C₃ (RR), F₃^2 = H₃ (NSNS)
  M^{-1} = diag(e^Φ, e^{-Φ}) (diagonal F1 case)

So M^{-1}_{ab}(F₃^a·F₃^b) = e^Φ(F₃·F₃) + e^{-Φ}(H₃·H₃).

For F1: only H₃, so this reduces to e^{-Φ}(H₃·H₃) = H^{1/2}(H₃·H₃)... but
the standard equation has just (H₃)² with no dilaton prefactor in string frame.

Let me compute everything numerically and let the algebra tell us the answer.
"""
import sys; sys.path.insert(0, 'src')
import sympy as sp
from sympy import Rational as R, ln, sqrt
from sugra import (HarmonicFunction, Metric,
                   FormField, exterior_derivative,
                   form_contraction, form_norm_squared)

# ===================================================================
# Setup: coordinates and harmonic function
# ===================================================================
print("=" * 70)
print("EXP 17: String-frame 12d equation")
print("=" * 70)

wv_coords = list(sp.symbols('t x1', real=True))
z1s, z2s = sp.symbols('z1 z2', real=True)
harmonic_coords = list(sp.symbols('y0:8', real=True))
coords = wv_coords + [z1s, z2s] + harmonic_coords  # 12 coords
D = 12

hf = HarmonicFunction(transverse_coords=harmonic_coords)
H = sp.Function('H')(hf.r_expr)
y0 = sp.Symbol('y0', real=True)

# ===================================================================
# Part A: F1 string-frame 12d metric and Ricci
# ===================================================================
print("\n" + "=" * 70)
print("PART A: F1 string-frame 12d metric")
print("=" * 70)

# String-frame 10d F1: g^S = H^{-1} ds²_{1,1} + ds²_8
# Torus: M₀ = diag(H^{1/2}, H^{-1/2})  [z1=v, z2=u]
M0 = sp.Matrix([[H**R(1, 2), 0],
                 [0, H**R(-1, 2)]])

g_sf_f1 = sp.zeros(D, D)
g_sf_f1[0, 0] = -1 / H       # t (string frame)
g_sf_f1[1, 1] = 1 / H        # x1 (string frame)
g_sf_f1[2, 2] = M0[0, 0]     # z1 = H^{1/2}
g_sf_f1[3, 3] = M0[1, 1]     # z2 = H^{-1/2}
for k in range(8):
    g_sf_f1[4 + k, 4 + k] = sp.Integer(1)  # y_k (string frame: flat transverse)

metric_sf_f1 = Metric(g_sf_f1, coords)
print(f"  F1 string-frame 12d metric built (diagonal={metric_sf_f1.is_diagonal})")

print("Computing Ricci tensor (F1 string frame)...")
Ric_sf_f1 = metric_sf_f1.ricci_tensor(simplify_func=sp.cancel)
print("Done.")

# Print ℛ/g for each block
print("\nℛ^{SF}_{MM}/g_{MM} per block (F1):")
for i, name in [(0, 't'), (1, 'x1'), (2, 'z1'), (3, 'z2'), (4, 'y0'), (5, 'y1')]:
    R_val = hf.substitute(sp.cancel(Ric_sf_f1[i, i]))
    g_val = hf.substitute(sp.cancel(metric_sf_f1.matrix[i, i]))
    R_over_g = sp.cancel(R_val / g_val)
    print(f"  [{name}]: ℛ/g = {R_over_g}")

# ===================================================================
# Part B: (1,1)-string string-frame 12d metric and Ricci
# ===================================================================
print("\n" + "=" * 70)
print("PART B: (1,1)-string string-frame 12d metric")
print("=" * 70)

# SL(2,Z): Λ = [[1,0],[1,1]]
Lam = sp.Matrix([[1, 0], [1, 1]])
M_pq = sp.cancel(Lam * M0 * Lam.T)
print(f"Torus metric M = Λ M₀ Λ^T:")
for i in range(2):
    print(f"  [{sp.cancel(M_pq[i, 0])}, {sp.cancel(M_pq[i, 1])}]")

# String-frame 12d metric for (1,1): SAME 10d string-frame block + new M
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
print(f"  (1,1) string-frame 12d metric built (diagonal={metric_sf_pq.is_diagonal})")

print("Computing Ricci tensor ((1,1) string frame)...")
Ric_sf_pq = metric_sf_pq.ricci_tensor(simplify_func=sp.cancel)
print("Done.")

# ===================================================================
# Part C: Compare Ricci tensors — 10d directions
# ===================================================================
print("\n" + "=" * 70)
print("PART C: ℛ^{SF}(F1) vs ℛ^{SF}(1,1) for 10d directions")
print("=" * 70)

all_10d_match = True
for i, name in [(0, 't'), (1, 'x1'), (4, 'y0'), (5, 'y1')]:
    R_f1 = hf.substitute(sp.cancel(Ric_sf_f1[i, i]))
    R_pq = hf.substitute(sp.cancel(Ric_sf_pq[i, i]))
    diff = sp.cancel(R_f1 - R_pq)
    status = '✓' if diff == 0 else f'✗ diff={diff}'
    if diff != 0:
        all_10d_match = False
    print(f"  [{name}]: F1={R_f1}, (1,1)={R_pq}, match={status}")

# Check more y-coords
for k in range(2, 8):
    i = 4 + k
    R_f1 = hf.substitute(sp.cancel(Ric_sf_f1[i, i]))
    R_pq = hf.substitute(sp.cancel(Ric_sf_pq[i, i]))
    diff = sp.cancel(R_f1 - R_pq)
    if diff != 0:
        all_10d_match = False
        print(f"  [y{k}]: ✗ diff={diff}")

print(f"\n  10d Ricci (string frame) match: {'★ YES ★' if all_10d_match else 'NO'}")

# ===================================================================
# Part D: G₄ form and contractions in string frame
# ===================================================================
print("\n" + "=" * 70)
print("PART D: G₄ contractions in string frame")
print("=" * 70)

# F1: G₄ = H₃ ∧ dz₂, C_{t,x1,z2} = H^{-1}
C_f1 = FormField(rank=3, dim=D)
C_f1[(0, 1, 3)] = 1 / H  # C_{t,x1,z2}
G4_f1 = exterior_derivative(C_f1, coords)

FF_f1_sf = form_contraction(G4_f1, metric_sf_f1)
S_f1_sf = form_norm_squared(G4_f1, metric_sf_f1)
S_f1_sf_val = hf.substitute(sp.cancel(S_f1_sf))
print(f"|G₄|²(SF, F1) = {S_f1_sf_val}")

# (1,1): G₄ = -H₃∧dz₁ + H₃∧dz₂
C_pq = FormField(rank=3, dim=D)
C_pq[(0, 1, 2)] = -1 / H  # C_{t,x1,z1} = -H^{-1}
C_pq[(0, 1, 3)] = 1 / H   # C_{t,x1,z2} = +H^{-1}
G4_pq = exterior_derivative(C_pq, coords)

FF_pq_sf = form_contraction(G4_pq, metric_sf_pq)
S_pq_sf = form_norm_squared(G4_pq, metric_sf_pq)
S_pq_sf_val = hf.substitute(sp.cancel(S_pq_sf))
print(f"|G₄|²(SF, (1,1)) = {S_pq_sf_val}")
print(f"Match: {sp.cancel(S_f1_sf_val - S_pq_sf_val) == 0}")

# ===================================================================
# Part E: Build H₃ and compute doublet contractions
# ===================================================================
print("\n" + "=" * 70)
print("PART E: H₃ doublet contractions in string frame")
print("=" * 70)

# H₃ lives in the 10d subspace
# For contraction, we need a 10d-style metric.
# In string frame, the 10d metric is g^S = H^{-1}ds²_{wv} + ds²_8
# But since H₃ only has (t,x1,yk) indices, we can use the full 12d metric
# because the z-block doesn't participate.

# Build H₃ form (same for all (p,q) — it's the universal 3-form field strength)
C_h3 = FormField(rank=2, dim=D)
C_h3[(0, 1)] = 1 / H  # B_{t,x1} = H^{-1}
H3_form = exterior_derivative(C_h3, coords)

print("Non-zero H₃ components:")
for idx, val in H3_form.nonzero_components.items():
    names = [str(coords[i]) for i in idx]
    v = hf.substitute(sp.cancel(val))
    print(f"  H₃[{','.join(names)}] = {v}")

# Compute H₃ contractions using string-frame F1 metric
# (same metric for 10d directions across all (p,q))
FF_H3_sf = form_contraction(H3_form, metric_sf_f1)
S_H3_sf = form_norm_squared(H3_form, metric_sf_f1)
S_H3_sf_val = hf.substitute(sp.cancel(S_H3_sf))
print(f"\n|H₃|²(string frame) = {S_H3_sf_val}")

# (H₃²)_{MN} per block
print("\n(H₃·H₃)_{MM}/g_{MM} per block (string frame):")
for i, name in [(0, 't'), (1, 'x1'), (2, 'z1'), (3, 'z2'), (4, 'y0'), (5, 'y1')]:
    ff_val = hf.substitute(sp.cancel(FF_H3_sf[i, i]))
    g_val = hf.substitute(sp.cancel(metric_sf_f1.matrix[i, i]))
    if g_val != 0:
        ratio = sp.cancel(ff_val / g_val)
    else:
        ratio = 'N/A (g=0)'
    print(f"  [{name}]: FF/g = {ratio}, FF = {ff_val}")

# ===================================================================
# Part F: Dilaton kinetic term in string frame
# ===================================================================
print("\n" + "=" * 70)
print("PART F: Dilaton kinetic term ∂Φ∂Φ")
print("=" * 70)

# Φ = -(1/2) ln H → ∂_k Φ = -(1/2)(H'/H)(y_k/r)
# ∂_m Φ ∂_n Φ = (1/4)(H')²/H² (y_m y_n / r²) for m,n transverse
# g^{mn} ∂_m Φ ∂_n Φ = (1/4)(H')²/H² × g^{yk,yk} × Σ(yk²/r²)
#   = (1/4)(H')²/H² × 1 (string frame: g^{yk}=1) × 1
#   = (1/4)(H')²/H²  = (∂Φ)²

dPhi_sq = R(1, 4) * hf.Hp**2 / hf.H**2
print(f"|∂Φ|²(string frame) = {dPhi_sq}")

# Second covariant derivative: ∇_m∂_nΦ
# For m=n=y_k (transverse, flat in string frame except for warp factors...
# wait, string frame transverse is FLAT: g_{yk}=1.
# So ∇_m = ∂_m for the transverse y-directions? NO — the full 12d metric
# has off-diagonal Christoffel symbols connecting y to wv and z.
# But for the 10d string-frame metric, Γ^{y}_{y,y} comes from the flat transverse
# and the warped worldvolume. Need to compute this properly.

# Let me compute ∇∂Φ numerically via the Ricci identity instead.
# The 10d string-frame IIB equation for F1 (only NSNS sector) is:
#   R^S_{mn} + 2∇_m∂_nΦ - (1/4)(H₃)²_{mn} = 0
# So ∇_m∂_nΦ = (1/2)[(1/4)(H₃²)_{mn} - R^S_{mn}]
# But we need the 10d R^S, not the 12d one.
# 12d = 10d + KK torus contribution, so let's first get the 10d string-frame Ricci.

# Build 10d string-frame metric
coords_10 = wv_coords + harmonic_coords
g_sf_10 = sp.zeros(10, 10)
g_sf_10[0, 0] = -1 / H
g_sf_10[1, 1] = 1 / H
for k in range(8):
    g_sf_10[2 + k, 2 + k] = sp.Integer(1)

metric_sf_10 = Metric(g_sf_10, coords_10)
print("Computing 10d string-frame Ricci tensor...")
Ric_sf_10 = metric_sf_10.ricci_tensor(simplify_func=sp.cancel)
print("Done.")

print("\n10d string-frame ℛ^S/g per block:")
for i, name in [(0, 't'), (1, 'x1'), (2, 'y0'), (3, 'y1')]:
    R_val = hf.substitute(sp.cancel(Ric_sf_10[i, i]))
    g_val = hf.substitute(sp.cancel(metric_sf_10.matrix[i, i]))
    print(f"  [{name}]: ℛ/g = {sp.cancel(R_val / g_val)}, ℛ = {R_val}")

# ===================================================================
# Part G: Verify 10d string-frame F1 equation
# ===================================================================
print("\n" + "=" * 70)
print("PART G: Verify 10d string-frame equation R^S + 2∇∂Φ = (1/4)(H₃)²")
print("=" * 70)

# Compute ∇_m∂_nΦ for diagonal directions
# Φ = -(1/2) ln H
# ∂_k Φ = -(1/2) H'/H · y_k/r  for k transverse
# ∂_t Φ = ∂_x1 Φ = 0

# ∇_m∂_nΦ = ∂_m∂_nΦ - Γ^ρ_{mn}∂_ρΦ
# Since ∂Φ only has transverse components, Γ^{y_k}_{mn} are the relevant ones.

# For the string frame metric, compute Christoffel symbols
print("Computing 10d Christoffel symbols (string frame)...")
chris_10 = metric_sf_10.christoffel(simplify_func=sp.cancel)
print(f"  {len(chris_10)} non-zero symbols")

# Compute ∇∂Φ for diagonal components
# For m=n=t (index 0):
# ∂_t∂_tΦ = 0 (Φ independent of t)
# Γ^ρ_{tt}∂_ρΦ = Σ_k Γ^{y_k}_{tt} · ∂_{y_k}Φ
# ∂_{y_k}Φ = -(1/2)(H'/H)(y_k/r)

nabla_dPhi = sp.zeros(10, 10)

for m in range(10):
    for n in range(m, 10):
        # ∂_m∂_nΦ
        # Φ = -(1/2)ln H, H = H(r), r = r(y_k)
        # ∂_m Φ = 0 unless m is transverse (index >= 2 in 10d)
        # ∂_m∂_n Φ = 0 unless both m,n are transverse
        d2Phi = sp.Integer(0)
        if m >= 2 and n >= 2:  # both transverse
            ym = coords_10[m]
            yn = coords_10[n]
            # ∂_m∂_n Φ = -(1/2) ∂_m[H'/H · y_n/r]
            Phi_expr = -R(1, 2) * sp.log(H)
            # Use chain rule
            d2Phi = sp.diff(sp.diff(Phi_expr, yn), ym)

        # Γ^ρ_{mn} ∂_ρ Φ  (only transverse ρ contribute)
        gamma_term = sp.Integer(0)
        for rho in range(2, 10):  # transverse indices in 10d
            key = (rho, m, n)
            if key in chris_10:
                yr = coords_10[rho]
                dPhi_rho = sp.diff(-R(1, 2) * sp.log(H), yr)
                gamma_term += chris_10[key] * dPhi_rho

        nabla_dPhi[m, n] = d2Phi - gamma_term
        nabla_dPhi[n, m] = nabla_dPhi[m, n]

# Substitute and simplify
print("\n∇_m∂_nΦ / g_{mn} per block:")
for i, name in [(0, 't'), (1, 'x1'), (2, 'y0'), (3, 'y1')]:
    val = hf.substitute(sp.cancel(nabla_dPhi[i, i]))
    g_val = hf.substitute(sp.cancel(metric_sf_10.matrix[i, i]))
    ratio = sp.cancel(val / g_val)
    print(f"  [{name}]: ∇∂Φ/g = {ratio}, ∇∂Φ = {val}")

# Now verify: R^S_{mn} + 2∇_m∂_nΦ - (1/4)(H₃²)_{mn} = 0
# Need H₃ contractions in 10d string frame
# Build H₃ in 10d
C_h3_10 = FormField(rank=2, dim=10)
C_h3_10[(0, 1)] = 1 / H
H3_10 = exterior_derivative(C_h3_10, coords_10)
FF_H3_10 = form_contraction(H3_10, metric_sf_10)
S_H3_10 = form_norm_squared(H3_10, metric_sf_10)
S_H3_10_val = hf.substitute(sp.cancel(S_H3_10))
print(f"\n|H₃|²(10d string frame) = {S_H3_10_val}")

print("\nVerifying: R^S + 2∇∂Φ = (1/4)(H₃²) for diagonal components:")
all_sf_10d_pass = True
for i, name in [(0, 't'), (1, 'x1'), (2, 'y0'), (3, 'y1')]:
    R_val = hf.substitute(sp.cancel(Ric_sf_10[i, i]))
    nabla_val = hf.substitute(sp.cancel(nabla_dPhi[i, i]))
    ff_val = hf.substitute(sp.cancel(FF_H3_10[i, i]))

    lhs = sp.cancel(R_val + 2 * nabla_val)
    rhs = sp.cancel(R(1, 4) * ff_val)
    diff = sp.cancel(lhs - rhs)
    status = '✓' if diff == 0 else f'✗ diff={diff}'
    if diff != 0:
        all_sf_10d_pass = False
    print(f"  [{name}]: R^S + 2∇∂Φ = {lhs}, (1/4)H₃² = {rhs}, {status}")

# Also try: R^S + 2∇∂Φ = (1/4)(H₃²) - (1/6)|H₃|²g
# (trace-modified version)
print("\nTrying: R^S + 2∇∂Φ = (1/4)[(H₃²) - (1/3)|H₃|²g] ?")
for i, name in [(0, 't'), (1, 'x1'), (2, 'y0'), (3, 'y1')]:
    R_val = hf.substitute(sp.cancel(Ric_sf_10[i, i]))
    nabla_val = hf.substitute(sp.cancel(nabla_dPhi[i, i]))
    ff_val = hf.substitute(sp.cancel(FF_H3_10[i, i]))
    g_val = hf.substitute(sp.cancel(metric_sf_10.matrix[i, i]))

    lhs = sp.cancel(R_val + 2 * nabla_val)
    rhs = sp.cancel(R(1, 4) * (ff_val - R(1, 3) * S_H3_10_val * g_val))
    diff = sp.cancel(lhs - rhs)
    status = '✓' if diff == 0 else f'✗ diff={diff}'
    print(f"  [{name}]: {status}")

# If the basic form doesn't work, try to find the right coefficient
# R^S + 2∇∂Φ = α (H₃²)_{mn} + β |H₃|² g_{mn}
# For each block: α ff + β S g = R + 2∇∂Φ
# Two unknowns, use t and y0 blocks to solve
print("\nSolving for coefficients: R^S + 2∇∂Φ = α(H₃²) + β|H₃|²g")
alpha, beta = sp.symbols('alpha beta')

R_t = hf.substitute(sp.cancel(Ric_sf_10[0, 0]))
R_y0 = hf.substitute(sp.cancel(Ric_sf_10[2, 2]))
nd_t = hf.substitute(sp.cancel(nabla_dPhi[0, 0]))
nd_y0 = hf.substitute(sp.cancel(nabla_dPhi[2, 2]))
ff_t = hf.substitute(sp.cancel(FF_H3_10[0, 0]))
ff_y0 = hf.substitute(sp.cancel(FF_H3_10[2, 2]))
g_t = hf.substitute(sp.cancel(metric_sf_10.matrix[0, 0]))
g_y0 = hf.substitute(sp.cancel(metric_sf_10.matrix[2, 2]))

lhs_t = sp.cancel(R_t + 2 * nd_t)
lhs_y0 = sp.cancel(R_y0 + 2 * nd_y0)

eq1 = sp.Eq(alpha * ff_t + beta * S_H3_10_val * g_t, lhs_t)
eq2 = sp.Eq(alpha * ff_y0 + beta * S_H3_10_val * g_y0, lhs_y0)

sol = sp.solve([eq1, eq2], [alpha, beta])
print(f"  α = {sol.get(alpha, 'N/A')}")
print(f"  β = {sol.get(beta, 'N/A')}")

# Verify with remaining blocks
if alpha in sol and beta in sol:
    a_val, b_val = sol[alpha], sol[beta]
    print(f"\nVerifying α={a_val}, β={b_val} for all blocks:")
    for i, name in [(0, 't'), (1, 'x1'), (2, 'y0'), (3, 'y1')]:
        R_val = hf.substitute(sp.cancel(Ric_sf_10[i, i]))
        nabla_val = hf.substitute(sp.cancel(nabla_dPhi[i, i]))
        ff_val = hf.substitute(sp.cancel(FF_H3_10[i, i]))
        g_val = hf.substitute(sp.cancel(metric_sf_10.matrix[i, i]))

        lhs = sp.cancel(R_val + 2 * nabla_val)
        rhs = sp.cancel(a_val * ff_val + b_val * S_H3_10_val * g_val)
        diff = sp.cancel(lhs - rhs)
        status = '✓' if diff == 0 else f'✗ diff={diff}'
        print(f"  [{name}]: {status}")

# ===================================================================
# Part H: KK decomposition — 12d vs 10d Ricci in string frame
# ===================================================================
print("\n" + "=" * 70)
print("PART H: KK decomposition ℛ^{12d} = R^{10d} + KK")
print("=" * 70)

# The 12d Ricci for the F1 case (from Part A) vs the 10d string-frame Ricci
# KK_{mn} = ℛ^{12d}_{mn} - R^{10d}_{mn}
# This should encode the torus moduli kinetic energy

print("\nKK contribution: ℛ^{12d} − R^{10d} per block (F1):")
for i_12, i_10, name in [(0, 0, 't'), (1, 1, 'x1'), (4, 2, 'y0'), (5, 3, 'y1')]:
    R_12 = hf.substitute(sp.cancel(Ric_sf_f1[i_12, i_12]))
    R_10 = hf.substitute(sp.cancel(Ric_sf_10[i_10, i_10]))
    kk = sp.cancel(R_12 - R_10)
    g_val = hf.substitute(sp.cancel(metric_sf_f1.matrix[i_12, i_12]))
    kk_over_g = sp.cancel(kk / g_val)
    print(f"  [{name}]: KK = {kk}, KK/g = {kk_over_g}")

# Also check: is KK ∝ Tr(∂M ∂M^{-1})?
# Tr(dM₀/dH · dM₀^{-1}/dH)
M0_inv = sp.cancel(M0.inv())
dM0_dH = sp.cancel(sp.diff(M0, H))
dM0inv_dH = sp.cancel(sp.diff(M0_inv, H))
trace_M = sp.cancel(sp.trace(dM0_dH * dM0inv_dH))
print(f"\nTr(dM₀/dH · dM₀^{{-1}}/dH) = {trace_M}")

# The KK term for y_k direction should be:
# KK_{yk,yk} = -(1/4) g^{12d}_{yk,yk} × Σ_l (H')²(y_l²/r²) × Tr stuff
# In string frame g^{yk}=1, so:
# KK_{yk} = -(1/4)(y_k²/r²)(H')² × Tr_product (for aniso part)
# Plus isotropic part from the Laplacian term

# Let me compute Tr(∂_yk M · ∂_yk M^{-1}) for the F1 diagonal case
# ∂_yk M = (dM/dH)(H')(y_k/r)
# So Tr(∂M·∂M^{-1}) evaluated at y_k:
# = (H')²(y_k²/r²) × Tr(dM/dH · dM^{-1}/dH)

# For the (1,1) torus metric:
M_pq_inv = sp.cancel(M_pq.inv())
dMpq_dH = sp.cancel(sp.diff(M_pq, H))
dMpqinv_dH = sp.cancel(sp.diff(M_pq_inv, H))
trace_Mpq = sp.cancel(sp.trace(dMpq_dH * dMpqinv_dH))
print(f"Tr(dM/dH · dM^{{-1}}/dH) [(1,1)] = {trace_Mpq}")
print(f"Tr match F1 vs (1,1): {sp.cancel(trace_M - trace_Mpq) == 0}")

# ===================================================================
# Part I: Full 12d string-frame equation for F1
# ===================================================================
print("\n" + "=" * 70)
print("PART I: Derive full 12d string-frame equation (F1)")
print("=" * 70)

# The 12d equation should be:
# ℛ^{12d}_{MN} = [form terms] + [scalar terms]
#
# From 10d IIB string frame: R^S + 2∇∂Φ = (form terms)
# So: ℛ^{12d} = R^S + KK
#             = (form terms) - 2∇∂Φ + KK
#
# First, let me just directly check: does ℛ^{12d}(F1, SF) satisfy a simple
# equation with G₄ and M?

# Try: ℛ = (1/2)[FF(G₄)/3! - λ|G₄|²g] for various λ
print("\nDirect test: ℛ^{SF,12d} = (1/2)[FF/3! - λ|G₄|²g] ?")
FF_G4_f1_sf = FF_f1_sf  # already computed in Part D
for test_lam in [R(1, 4), R(3, 10), R(1, 2), R(1, 6), R(1, 8)]:
    all_ok = True
    for i in [0, 1, 2, 3, 4]:
        R_val = hf.substitute(sp.cancel(Ric_sf_f1[i, i]))
        ff_val = hf.substitute(sp.cancel(FF_G4_f1_sf[i, i]))
        g_val = hf.substitute(sp.cancel(metric_sf_f1.matrix[i, i]))
        T = sp.cancel(R(1, 2) * (ff_val / 6 - test_lam * S_f1_sf_val * g_val))
        diff = sp.cancel(R_val - T)
        if diff != 0:
            all_ok = False
            break
    print(f"  λ={test_lam}: {'PASS' if all_ok else 'FAIL'}")

# Since simple T^{G₄}(λ) won't work (the dilaton term doesn't have the right
# structure), let's decompose: what is the RESIDUAL after the best-fit form term?

# Find λ that works for wv (t block):
R_t_12 = hf.substitute(sp.cancel(Ric_sf_f1[0, 0]))
ff_t_12 = hf.substitute(sp.cancel(FF_G4_f1_sf[0, 0]))
g_t_12 = hf.substitute(sp.cancel(metric_sf_f1.matrix[0, 0]))

# ℛ_t = (1/2)(FF_t/6 - λ S g_t)
# → λ = (FF_t/6 - 2ℛ_t) / (S g_t)
lam_t = sp.cancel((ff_t_12 / 6 - 2 * R_t_12) / (S_f1_sf_val * g_t_12))
print(f"\nλ that fits t-block: {lam_t}")

# Check if this λ works for other blocks
print(f"\nUsing λ = {lam_t}:")
for i, name in [(0, 't'), (1, 'x1'), (2, 'z1'), (3, 'z2'), (4, 'y0'), (5, 'y1')]:
    R_val = hf.substitute(sp.cancel(Ric_sf_f1[i, i]))
    ff_val = hf.substitute(sp.cancel(FF_G4_f1_sf[i, i]))
    g_val = hf.substitute(sp.cancel(metric_sf_f1.matrix[i, i]))
    T = sp.cancel(R(1, 2) * (ff_val / 6 - lam_t * S_f1_sf_val * g_val))
    diff = sp.cancel(R_val - T)
    diff_over_g = sp.cancel(diff / g_val) if g_val != 0 else diff
    status = '✓' if diff == 0 else f'✗ (ℛ-T)/g = {diff_over_g}'
    print(f"  [{name}]: {status}")

# ===================================================================
# Part J: Compute residuals and identify structure
# ===================================================================
print("\n" + "=" * 70)
print("PART J: Residual structure ℛ - T^{G₄}(λ) per block")
print("=" * 70)

# Compute residuals for several λ values and tabulate
for test_lam_val, lam_name in [(R(1, 4), "1/4"), (R(1, 8), "1/8")]:
    print(f"\n--- λ = {lam_name} ---")
    for i, name in [(0, 't'), (1, 'x1'), (2, 'z1'), (3, 'z2'), (4, 'y0'), (5, 'y1')]:
        R_val = hf.substitute(sp.cancel(Ric_sf_f1[i, i]))
        ff_val = hf.substitute(sp.cancel(FF_G4_f1_sf[i, i]))
        g_val = hf.substitute(sp.cancel(metric_sf_f1.matrix[i, i]))
        T = sp.cancel(R(1, 2) * (ff_val / 6 - test_lam_val * S_f1_sf_val * g_val))
        resid = sp.cancel(R_val - T)
        resid_g = sp.cancel(resid / g_val) if g_val != 0 else resid
        print(f"  [{name}]: (ℛ-T)/g = {resid_g}")

# ===================================================================
# Part K: Try doublet form with explicit M dependence
# ===================================================================
print("\n" + "=" * 70)
print("PART K: Doublet formulation in string frame")
print("=" * 70)

# In string frame, the IIB equation involves:
# For the NSNS 3-form sector: (1/4)(H₃²)_{mn}  (no dilaton factor in SF)
# For the RR 3-form sector: (1/4)e^{2Φ}(F₃²)_{mn}
#
# In doublet notation with M^{-1}_{ab} = diag(e^{Φ}/τ₂, e^{-Φ}/τ₂) ...
# Actually, M = (1/τ₂) [[|τ|², τ₁], [τ₁, 1]]
# M^{-1} = τ₂ [[1, -τ₁], [-τ₁, |τ|²]]
#
# For F1 with τ=iH^{1/2}: τ₁=0, τ₂=H^{1/2}, |τ|²=H
# M = diag(H^{1/2}, H^{-1/2})
# M^{-1} = diag(H^{-1/2}, H^{1/2})
#
# The string-frame equation uses:
# T^{SF}_{mn} = (1/4) Σ_{a,b} S_{ab} (F₃^a·F₃^b)_{mn}/2!
# where S = diag(e^{2Φ}, 1) = diag(H^{-1}, 1) for F1
# This is NOT simply M or M^{-1}.
#
# Actually S = M^{-1} × τ₂ × diag(e^Φ, e^{-Φ}) ... hmm let me think.
#
# Standard IIB in string frame (Polchinski Vol 2, eq 12.1.30):
# R_{mn} + 2∇_m∂_nΦ = (1/12) F₃·F₃_mn  [F₃ = RR 3-form field strength]
#                    + (1/12) H₃·H₃_mn   [H₃ = NSNS]
#
# Wait — but the RR field strength actually has an e^{2Φ} prefactor in the action:
# S = ∫ e^{-2Φ}(R + 4(∂Φ)²) - (1/2)|H₃|² - (1/2)e^{2Φ}|F₁|² - (1/2)|F₃'|² ...
# where F₃' = F₃ - C₀ H₃.
#
# The equations of motion in string frame are:
# R_{mn} + 2∇_m∂_nΦ = (1/4)H₃²_{mn} + (1/4)e^{2Φ}F₁ terms + (1/4)F₃'^2_{mn} + ...
#
# Actually for the 3-form sector it's:
# (1/4)H₃²_{mn} + (1/4)F₃'²_{mn}
# where H₃ and F₃' enter with the SAME coefficient (both have 1/4).

# For F1: only H₃ ≠ 0, F₃ = 0.
# R^S + 2∇∂Φ = (1/4)(H₃²)_{mn}  →  (1/4)(H₃·H₃)_{mm} for diagonal

# We already verified this in Part G. Let's just directly compute
# what the 12d equation looks like.

# KEY APPROACH: express everything in terms of 12d quantities.
# ℛ^{12d}_{mn} = R^{10d,S}_{mn} + KK_{mn}
#              = [(1/4)(H₃²)_{mn} - 2∇∂Φ] + KK_{mn}

# The dilaton and KK terms should combine into torus moduli terms.
# Let's compute each piece numerically.

print("\nDecomposition: ℛ^{12d} = (1/4)H₃² − 2∇∂Φ + KK")
print("All divided by g_{mm}:\n")
for i_12, i_10, name in [(0, 0, 't'), (1, 1, 'x1'), (4, 2, 'y0'), (5, 3, 'y1')]:
    R12 = hf.substitute(sp.cancel(Ric_sf_f1[i_12, i_12]))
    R10 = hf.substitute(sp.cancel(Ric_sf_10[i_10, i_10]))
    kk = sp.cancel(R12 - R10)
    nd = hf.substitute(sp.cancel(nabla_dPhi[i_10, i_10]))
    ff = hf.substitute(sp.cancel(FF_H3_10[i_10, i_10]))
    g_val = hf.substitute(sp.cancel(metric_sf_f1.matrix[i_12, i_12]))

    H3_term = sp.cancel(R(1, 4) * ff)
    dilaton_term = sp.cancel(-2 * nd)
    total = sp.cancel(H3_term + dilaton_term + kk)

    print(f"  [{name}]:")
    print(f"    (1/4)H₃²  /g = {sp.cancel(H3_term / g_val)}")
    print(f"    -2∇∂Φ     /g = {sp.cancel(dilaton_term / g_val)}")
    print(f"    KK         /g = {sp.cancel(kk / g_val)}")
    print(f"    sum        /g = {sp.cancel(total / g_val)}")
    print(f"    ℛ(12d)    /g = {sp.cancel(R12 / g_val)}")
    print(f"    check: {sp.cancel(total - R12) == 0}")

# ===================================================================
# Part L: Torus directions
# ===================================================================
print("\n" + "=" * 70)
print("PART L: Torus block equation (z1, z2)")
print("=" * 70)

# For torus directions in the 12d metric:
# ℛ^{12d}_{ab} = KK formula for internal metric M_{ab}
# The 10d Ricci for internal directions doesn't exist (they're extra dimensions).
# The full equation comes from ∂M/∂Φ equation of motion.

for i, name in [(2, 'z1'), (3, 'z2')]:
    R_val = hf.substitute(sp.cancel(Ric_sf_f1[i, i]))
    ff_val = hf.substitute(sp.cancel(FF_G4_f1_sf[i, i]))
    g_val = hf.substitute(sp.cancel(metric_sf_f1.matrix[i, i]))
    print(f"  [{name}]: ℛ = {R_val}, FF(G₄) = {ff_val}, g = {g_val}")
    print(f"    ℛ/g = {sp.cancel(R_val / g_val)}")

# ===================================================================
# Part M: Verify (1,1) equation with doublet
# ===================================================================
print("\n" + "=" * 70)
print("PART M: (1,1)-string — same equation?")
print("=" * 70)

# For 10d directions: if the 10d string-frame Ricci is the same,
# and the KK term only depends on Tr(dM·dM^{-1}) which is SL(2,Z) invariant,
# then the 12d Ricci is the same for all (p,q)-strings.

# Let's verify this explicitly.
print("\n(1,1) vs F1: 12d Ricci in string frame, 10d directions:")
all_12d_match = True
for i, name in [(0, 't'), (1, 'x1'), (4, 'y0'), (5, 'y1')]:
    R_f1 = hf.substitute(sp.cancel(Ric_sf_f1[i, i]))
    R_pq = hf.substitute(sp.cancel(Ric_sf_pq[i, i]))
    diff = sp.cancel(R_f1 - R_pq)
    status = '✓' if diff == 0 else f'✗'
    if diff != 0:
        all_12d_match = False
        print(f"  [{name}]: F1={R_f1}, (1,1)={R_pq}, diff={diff}")
    else:
        print(f"  [{name}]: {status}")

for k in range(2, 8):
    i = 4 + k
    R_f1 = hf.substitute(sp.cancel(Ric_sf_f1[i, i]))
    R_pq = hf.substitute(sp.cancel(Ric_sf_pq[i, i]))
    diff = sp.cancel(R_f1 - R_pq)
    if diff != 0:
        all_12d_match = False
        print(f"  [y{k}]: ✗")

print(f"\n  10d directions (12d Ricci, SF): {'★ MATCH ★' if all_12d_match else 'MISMATCH'}")

# For torus directions: G₄ and ℛ differ between F1 and (1,1)
# but the EQUATION should be the same (SL(2,Z) covariant)
print("\nTorus block: (1,1)-string")
for i in range(2, 4):
    for j in range(i, 4):
        R_val = hf.substitute(sp.cancel(Ric_sf_pq[i, j]))
        ff_val = hf.substitute(sp.cancel(FF_pq_sf[i, j]))
        ni = 'z1' if i == 2 else 'z2'
        nj = 'z1' if j == 2 else 'z2'
        g_val = hf.substitute(sp.cancel(metric_sf_pq.matrix[i, j]))
        print(f"  [{ni},{nj}]: ℛ = {R_val}, g = {g_val}")

# ===================================================================
# Summary
# ===================================================================
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print("""
EXP 17: String-frame 12d equation

1. String-frame 12d metric: ds²₁₂ = g^S_{mn}dx^m dx^n + M_{ab}dz^a dz^b
   where g^S = H^{-1}ds²_{wv} + ds²_8 (SL(2,Z) invariant)

2. 10d string-frame Ricci: F1 vs (1,1) match for 10d directions?
   → See Part C results above

3. 12d Ricci in string frame: F1 vs (1,1) match for 10d directions?
   → See Part M results above

4. KK decomposition: ℛ^{12d} = R^{10d} + KK, where KK encodes torus moduli.
   → See Part H results above

5. 10d equation: R^S + 2∇∂Φ = α(H₃²) + β|H₃|²g
   → See Part G results (coefficients α, β)

6. The 12d equation combines: form terms + dilaton + KK → SL(2,Z) covariant?
   → See Part K decomposition
""")
