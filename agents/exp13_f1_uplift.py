"""Experiment 13: F1-string uplift — derive 12d T from 10d IIB equations.

Key insight: the 12d metric IS the F1-string (fundamental string) uplift:
  ds²₁₂ = ds²₁₀(F1, Einstein frame) + M_{ij}(Φ) dz^i dz^j

where:
  ds²₁₀ = H^{-3/4} ds²_{1,1} + H^{1/4} ds²_8
  Φ = -(1/2) ln H  →  τ₂ = e^{-Φ} = H^{1/2}
  M = diag(H^{-1/2}, H^{1/2})  [for C₀=0]

This gives: ds²₁₂ = H^{-3/4} ds²_{1,1} + H^{-1/2} dz₂² + H^{1/2} dz₁² + H^{1/4} ds²_8
matching the metric in the problem (with z₁↔z₂ labeling: z₂=u, z₁=v).

The F1 string in 10d satisfies:
  R_{mn} = (1/2)∂_mΦ∂_nΦ + (1/2)e^{-Φ}[H₃²_{mn}/2! - (1/4)|H₃|² g_{mn}]

where H₃ = dB₂ with B_{tx} = H^{-1} (→ D1 form in 12d: G₄ = H₃ ∧ dz₂).

This experiment:
Part A: Verify F1-string in 10d
Part B: Compute 12d ℛ_{MN} - 10d R_{mn} = KK torus contribution
Part C: Identify the KK terms and express the full 12d equation
"""
import sys; sys.path.insert(0, 'src')
import sympy as sp
from sympy import Rational as R
from sugra import (HarmonicFunction, warped_product,
                   FormField, exterior_derivative,
                   form_stress_energy, form_contraction,
                   form_norm_squared)

# ===================================================================
# Part A: Verify F1-string in 10d
# ===================================================================
print("="*70)
print("PART A: F1-string verification in 10d")
print("="*70)

# 10d setup
wv_coords_10 = list(sp.symbols('t x1', real=True))
harmonic_coords = list(sp.symbols('y0:8', real=True))
coords_10 = wv_coords_10 + harmonic_coords
y0 = sp.Symbol('y0', real=True)

hf = HarmonicFunction(transverse_coords=harmonic_coords)
H = sp.Function('H')(hf.r_expr)

# 10d F1-string metric (Einstein frame): H^{-3/4} ds²_{1,1} + H^{1/4} ds²_8
metric_10 = warped_product(
    warp_factors=[H**R(-3,4), H**R(1,4)],
    block_dims=[2, 8],
    block_signatures=['lorentzian', 'euclidean'],
    coordinates=coords_10,
)

print("Computing 10d Ricci tensor...")
Ric_10 = metric_10.ricci_tensor(simplify_func=sp.cancel)
print("Done.\n")

# H₃ = dB₂ with B_{tx} = H^{-1}
B2 = FormField(rank=2, dim=10)
B2[(0, 1)] = 1/H  # B_{t,x1} = H^{-1}
H3 = exterior_derivative(B2, coords_10)

print("Non-zero H₃ components:")
for idx, val in H3.nonzero_components.items():
    print(f"  H₃{list(idx)} = {val}")

# H₃ contraction and norm in 10d
FF_H3 = form_contraction(H3, metric_10)
S_H3 = form_norm_squared(H3, metric_10)

# Dilaton: Φ = -(1/2) ln H → ∂_{yk}Φ = -(1/2)(H'/H)(yk/r)
# T^Φ contribution to R: (1/2) ∂_mΦ ∂_nΦ (trace-reversed = no subtraction)
# Only non-zero for m,n both transverse

# 10d Einstein equation: R_{mn} = (1/2)∂_mΦ∂_nΦ + T^{H₃}_{mn}
# where T^{H₃}_{mn} = (1/2)e^{-Φ}[FF_{mn}/2 - (1/4)|H₃|²g_{mn}]
# and (p-1)/(D-2) = 2/8 = 1/4 for H₃ in 10d.

# e^{-Φ} = H^{1/2}
ePhiInv = H**R(1,2)

print("\n10d F1-string: R_{mn} vs T_{mn} = (1/2)∂Φ² + T^{H₃}")
print("-" * 60)

for i, name in [(0,'t'), (1,'x1'), (2,'y0'), (3,'y1')]:
    R_val = hf.substitute(sp.cancel(Ric_10[i,i]))

    # H₃ stress-energy (trace-reversed): (1/2) e^{-Φ} [FF/2 - (1/4)|H₃|²g]
    ff = hf.substitute(sp.cancel(FF_H3[i,i]))
    s = hf.substitute(sp.cancel(S_H3))
    g = hf.substitute(sp.cancel(metric_10.matrix[i,i]))
    epi = hf.substitute(sp.cancel(ePhiInv))

    T_H3 = R(1,2) * epi * (ff/2 - R(1,4)*s*g)
    T_H3 = sp.cancel(T_H3)

    # Dilaton kinetic: (1/2) ∂_mΦ ∂_nΦ
    # ∂_{yk}Φ = -(1/2)(H'/H)(yk/r), only for transverse directions (i >= 2 in 10d coords)
    if i >= 2:  # transverse
        dPhi_sq_ii = R(1,4) * hf.Hp**2 / hf.H**2 * sp.Symbol(str(coords_10[i]))**2 / hf.r**2
        T_dil = R(1,2) * dPhi_sq_ii
    else:
        T_dil = sp.Integer(0)

    T_total = sp.cancel(T_H3 + T_dil) if T_dil != 0 else T_H3
    diff = sp.cancel(R_val - T_total)

    print(f"\n  [{name}] (i={i}):")
    print(f"    R₁₀     = {R_val}")
    print(f"    T^H₃    = {T_H3}")
    print(f"    T^Φ     = {T_dil}")
    print(f"    T_total = {T_total}")
    print(f"    R-T     = {diff}")
    if diff == 0:
        print(f"    ✓ MATCH")
    else:
        print(f"    ✗ MISMATCH")

# ===================================================================
# Part B: Compare 12d Ricci with 10d Ricci (KK contribution)
# ===================================================================
print("\n" + "="*70)
print("PART B: KK decomposition — 12d Ricci vs 10d Ricci")
print("="*70)

# 12d setup
z1s, z2s = sp.symbols('z1 z2', real=True)
coords_12 = wv_coords_10 + [z1s, z2s] + harmonic_coords

metric_12 = warped_product(
    warp_factors=[H**R(-3,4), H**R(1,2), H**R(-1,2), H**R(1,4)],
    block_dims=[2, 1, 1, 8],
    block_signatures=['lorentzian', 'euclidean', 'euclidean', 'euclidean'],
    coordinates=coords_12,
)

print("Computing 12d Ricci tensor...")
Ric_12 = metric_12.ricci_tensor(simplify_func=sp.cancel)
print("Done.\n")

# Map between 10d and 12d indices
# 10d: 0=t, 1=x1, 2=y0, 3=y1, ..., 9=y7
# 12d: 0=t, 1=x1, 2=z1, 3=z2, 4=y0, 5=y1, ..., 11=y7
idx_map_10_to_12 = {0: 0, 1: 1}  # t, x1
for k in range(8):
    idx_map_10_to_12[2+k] = 4+k  # yk

print("KK contribution: ΔR = ℛ_{mn}(12d) - R_{mn}(10d) for 10d indices")
print("-" * 60)

for i10, name in [(0,'t'), (1,'x1'), (2,'y0'), (3,'y1')]:
    i12 = idx_map_10_to_12[i10]

    R_10 = hf.substitute(sp.cancel(Ric_10[i10, i10]))
    R_12 = hf.substitute(sp.cancel(Ric_12[i12, i12]))

    delta = sp.cancel(R_12 - R_10)

    # Compare with dilaton kinetic (1/2)∂Φ∂Φ
    if i10 >= 2:  # transverse
        dPhi_sq_ii = R(1,4) * hf.Hp**2 / hf.H**2 * sp.Symbol(str(coords_10[i10]))**2 / hf.r**2
        dil_kin = R(1,2) * dPhi_sq_ii
    else:
        dil_kin = sp.Integer(0)

    extra = sp.cancel(delta - dil_kin) if dil_kin != 0 else delta

    print(f"\n  [{name}] (10d idx={i10}, 12d idx={i12}):")
    print(f"    R₁₀         = {R_10}")
    print(f"    ℛ₁₂         = {R_12}")
    print(f"    ΔR = ℛ-R    = {delta}")
    print(f"    (1/2)∂Φ²    = {dil_kin}")
    print(f"    ΔR-(1/2)∂Φ² = {extra}")

# Also compute ℛ for the z-directions (no 10d counterpart)
print("\n\nTorus Ricci components:")
for i12, name in [(2, 'z1'), (3, 'z2')]:
    R_z = hf.substitute(sp.cancel(Ric_12[i12, i12]))
    g_z = hf.substitute(sp.cancel(metric_12.matrix[i12, i12]))
    R_over_g = sp.cancel(R_z / g_z)
    print(f"  ℛ_{name} = {R_z}")
    print(f"  ℛ_{name}/g_{name} = {R_over_g}")

# ===================================================================
# Part C: Express full 12d ℛ as form contributions + remainder
# ===================================================================
print("\n" + "="*70)
print("PART C: Full 12d equation — what sources ℛ_{MN}?")
print("="*70)

# D1 form in 12d (= G₄ = H₃ ∧ dz₂)
C_D1 = FormField(rank=3, dim=12)
C_D1[(0,1,3)] = 1/H  # C_{t,x1,z2} = H^{-1}
F_D1 = exterior_derivative(C_D1, coords_12)

FF_D1 = form_contraction(F_D1, metric_12)
S_D1 = form_norm_squared(F_D1, metric_12)

# From 10d: R_{mn} = (1/2)∂Φ² + (1/2)e^{-Φ}[FF_{H3}/2 - (1/4)|H₃|²g]
# The H₃ part uplifts to G₄. Let's see what the correct 12d formula is.

# Compute: for each block, find coefficients a, b such that
# ℛ_{MM} = a × FF_D1_{MM}/3! + b × S_D1 × g_{MM}
# (These are the only scalar densities available from the 4-form)

print("\nFitting ℛ_{MN} = (a/6) FF_D1_{MN} + b S_D1 g_{MN} for each block:")
print("-" * 60)

s_d1 = hf.substitute(sp.cancel(S_D1))

results = {}
for i, name in [(0,'t'), (1,'x1'), (2,'z1'), (3,'z2'), (4,'y0'), (5,'y1')]:
    R_val = hf.substitute(sp.cancel(Ric_12[i,i]))
    ff_val = hf.substitute(sp.cancel(FF_D1[i,i]))
    g_val = hf.substitute(sp.cancel(metric_12.matrix[i,i]))

    # ℛ = (a/6) FF + b S g
    # Two unknowns (a, b). One equation per block.
    # For blocks where FF = 0 (like z1, pure y-blocks): ℛ = b S g → b = ℛ/(S g)
    # For blocks where FF ≠ 0: need another equation.

    ff_coeff = sp.cancel(ff_val / 6)
    sg_coeff = sp.cancel(s_d1 * g_val)

    if ff_coeff == 0:
        b_val = sp.cancel(R_val / sg_coeff) if sg_coeff != 0 else 'N/A'
        print(f"  [{name}]: FF=0, b = ℛ/(S·g) = {b_val}")
        results[name] = {'ff': ff_coeff, 'sg': sg_coeff, 'R': R_val, 'b': b_val, 'a': 'free'}
    else:
        print(f"  [{name}]: FF/6 = {ff_coeff}, S·g = {sg_coeff}")
        print(f"         ℛ = {R_val}")
        # Can't determine a and b from one equation. Express as:
        # a = 6(ℛ - b·S·g) / FF
        a_of_b = sp.cancel(6 * (R_val - sp.Symbol('b') * sg_coeff) / ff_val)
        print(f"         a(b) = {a_of_b}")
        results[name] = {'ff': ff_coeff, 'sg': sg_coeff, 'R': R_val}

# Since z1 and y have FF=0, they pin down b. Let's check if they give the same b.
print("\n\nConsistency check — b from different blocks:")
b_z1 = results.get('z1', {}).get('b')
b_y0 = results.get('y0', {}).get('b')
b_y1 = results.get('y1', {}).get('b')
print(f"  b[z1] = {b_z1}")
print(f"  b[y0] = {b_y0}")
print(f"  b[y1] = {b_y1}")

if b_z1 is not None and b_y0 is not None:
    diff_b = sp.cancel(b_z1 - b_y0.subs(y0, 0))
    print(f"  b[z1] - b[y0](iso) = {diff_b}")

# ===================================================================
# Part D: What does the 10d equation become in 12d?
# ===================================================================
print("\n" + "="*70)
print("PART D: 10d equation rewritten in 12d")
print("="*70)

# From 10d:
# R_{mn}(10d) = (1/2)∂_mΦ∂_nΦ + (1/2)e^{-Φ}[FF_{H3,mn}/2! - (1/4)|H₃|²g_{mn}]
#
# The H₃ contraction/norm in 10d vs the D1 (= G₄) contraction/norm in 12d:
# G₄ = H₃ ∧ dz₂
#
# G₄²_{mn} = ∑_{PQR} G₄_{mPQR} G₄_n^{PQR}
#           = ∑_{k} G₄_{m,x1/t,z2,yk} G₄_n^{x1/t,z2,yk}  (schematic)
#           = g^{z2z2} × H₃²_{mn}  (for m,n both in 10d block)
#
# |G₄|² = g^{z2z2} × |H₃|²_{10d}

print("Relation between 10d H₃ and 12d G₄ contractions:")
print()

# Compute H₃ norm in 10d
s_h3_10 = hf.substitute(sp.cancel(S_H3))
print(f"  |H₃|²(10d) = {s_h3_10}")

# D1 norm in 12d
s_d1_12 = hf.substitute(sp.cancel(S_D1))
print(f"  |G₄|²(12d) = {s_d1_12}")

# Ratio
g_z2 = hf.substitute(sp.cancel(metric_12.matrix[3,3]))
g_z2_inv = sp.cancel(1/g_z2)
print(f"  g^{'{z2,z2}'} = {g_z2_inv}")
print(f"  |H₃|² × g^{'{z2,z2}'} = {sp.cancel(s_h3_10 * g_z2_inv)}")
print(f"  Match |G₄|²? {sp.cancel(s_h3_10 * g_z2_inv - s_d1_12) == 0}")

# Check FF contraction relation
print("\nFF contraction relation per 10d block:")
for i10, name in [(0,'t'), (2,'y0')]:
    i12 = idx_map_10_to_12[i10]

    ff_h3_10 = hf.substitute(sp.cancel(FF_H3[i10, i10]))
    ff_d1_12 = hf.substitute(sp.cancel(FF_D1[i12, i12]))

    # G₄²_{mm} should = g^{z2z2} × H₃²_{mm}
    expected = sp.cancel(ff_h3_10 * g_z2_inv)
    match = sp.cancel(expected - ff_d1_12) == 0

    print(f"  [{name}]: H₃²(10d)={ff_h3_10}, G₄²(12d)={ff_d1_12}")
    print(f"         g^z2z2 × H₃² = {expected}, match G₄²? {match}")

# Now express the 10d equation in 12d variables
print("\n\nExpress 10d T^{H₃} in terms of 12d G₄:")
print("  T^{H₃}_{mn} = (1/2)e^{-Φ}[FF_{H₃,mn}/2 - (1/4)|H₃|²g_{mn}]")
print("  Using: FF_{H₃} = g_{z2z2} × FF_{G₄}")
print("  Using: |H₃|² = g_{z2z2} × |G₄|²")
print("  Using: e^{-Φ} = H^{1/2} = g_{z1z1}  (since g_{z1}=H^{1/2})")
print()
print("  → T^{H₃}_{mn} = (1/2) g_{z1} [g_{z2} FF_{G₄,mn}/2 - (1/4) g_{z2} |G₄|² g_{mn}]")
print("                 = (g_{z1}g_{z2}/2) [FF_{G₄,mn}/2 - (1/4)|G₄|² g_{mn}]")
print()

# g_{z1}*g_{z2} = H^{1/2} * H^{-1/2} = 1
gz1_gz2 = sp.cancel(hf.substitute(metric_12.matrix[2,2]) * hf.substitute(metric_12.matrix[3,3]))
print(f"  g_z1 × g_z2 = {gz1_gz2}")

print("\n  → T^{H₃}_{mn} = (1/2)[FF_{G₄,mn}/2 - (1/4)|G₄|² g_{mn}]")
print("  This is the standard p=4 formula with λ = 1/4 (not 3/10!)")
print()
print("  But: the FULL 10d equation is:")
print("    R_{mn}(10d) = T^{H₃}_{mn} + (1/2)∂Φ²_{mn}")
print()
print("  And: ℛ_{mn}(12d) = R_{mn}(10d) + KK_{mn}")
print()
print("  So: ℛ_{mn} = T^{H₃}_{mn} + (1/2)∂Φ² + KK_{mn}")
print("             = (1/2)[FF_{G₄}/2 - (1/4)|G₄|²g] + (1/2)∂Φ² + KK_{mn}")

# ===================================================================
# Part E: Numerically verify the 10d → 12d equation
# ===================================================================
print("\n" + "="*70)
print("PART E: Verify ℛ₁₂ = T^{G₄}(λ=1/4) + (1/2)∂Φ² + KK")
print("="*70)

# Compute T^{G₄} with λ = 1/4 (the 10d-inherited value)
lam_10d = R(1,4)

for i12, name in [(0,'t'), (1,'x1'), (2,'z1'), (3,'z2'), (4,'y0'), (5,'y1')]:
    R12 = hf.substitute(sp.cancel(Ric_12[i12, i12]))

    ff = hf.substitute(sp.cancel(FF_D1[i12, i12]))
    s = hf.substitute(sp.cancel(S_D1))
    g = hf.substitute(sp.cancel(metric_12.matrix[i12, i12]))

    # T^{G₄} with λ = 1/4
    T_G4 = R(1,2) * (ff/6 - lam_10d * s * g)
    T_G4 = sp.cancel(T_G4)

    # Dilaton kinetic: (1/2)∂_MΦ∂_NΦ
    # Φ = -(1/2)ln H, only transverse derivatives non-zero
    if i12 >= 4:  # transverse
        coord_name = coords_12[i12]
        T_dil = R(1,2) * R(1,4) * hf.Hp**2 / hf.H**2 * coord_name**2 / hf.r**2
    else:
        T_dil = sp.Integer(0)

    # KK contribution (= ℛ₁₂ - R₁₀ for 10d indices, or ℛ₁₂ for z-indices)
    if i12 in [0, 1]:  # wv
        i10 = i12
        R10 = hf.substitute(sp.cancel(Ric_10[i10, i10]))
        KK = sp.cancel(R12 - R10)
    elif i12 >= 4:  # transverse
        i10 = i12 - 2
        R10 = hf.substitute(sp.cancel(Ric_10[i10, i10]))
        KK = sp.cancel(R12 - R10)
    else:
        # z-directions: no 10d counterpart
        KK = 'N/A (torus direction)'
        R10 = None

    # Residual: ℛ₁₂ - T^{G₄}(λ=1/4) - T_dil
    total_rhs = sp.cancel(T_G4 + T_dil) if T_dil != 0 else T_G4
    resid = sp.cancel(R12 - total_rhs)

    print(f"\n  [{name}] (12d idx={i12}):")
    print(f"    ℛ₁₂     = {R12}")
    print(f"    T^G₄    = {T_G4}")
    print(f"    T^Φ     = {T_dil}")
    if R10 is not None:
        print(f"    R₁₀     = {R10}")
        print(f"    KK      = {KK}")
    print(f"    ℛ-T^G₄-T^Φ = {resid}")
    if resid == 0:
        print(f"    ✓ EXACT MATCH (no KK needed)")
    elif i12 in [2, 3]:
        # For z-directions: ℛ_{zz} is entirely from KK. What sources it?
        # Try: ℛ_{z} = ? × S_D1 × g_z
        if s != 0:
            ratio = sp.cancel(R12 / (s * g))
            print(f"    ℛ_z/(S·g_z) = {ratio}")
    else:
        # Express residual / g in terms of H', H, r
        res_over_g = sp.cancel(resid / g)
        print(f"    residual/g = {res_over_g}")

# ===================================================================
# Part F: Try finding the correct 12d formula
# ===================================================================
print("\n" + "="*70)
print("PART F: Search for correct 12d formula")
print("="*70)

# The 10d equation uses:
# - FF contraction (from H₃)
# - norm squared (from H₃)
# - dilaton kinetic (from Φ)
#
# In 12d, the dilaton is in M_{ij}. The "dilaton kinetic energy" becomes
# part of the KK Ricci. But it also relates to derivatives of g_{z1}/g_{z2}.
#
# Key observation: g_{z1}/g_{z2} = H^{1/2}/H^{-1/2} = H
# So ln(g_{z1}/g_{z2}) = ln H = -2Φ
#
# Define: Ψ = (1/2) ln(g_{z1}/g_{z2}) = (1/2) ln H = -Φ
# ∂_{yk}Ψ = (1/2)(H'/H)(yk/r)
# (∂Ψ)² with 12d metric: g^{yk yk} (∂Ψ)² = (1/H^{1/4})(1/4)(H'/H)² × yk²/r²

# Compute Ψ-gradient squared in 12d for each block
print("Gradient of Ψ = (1/2)ln(g_{z1}/g_{z2}) = (1/2)ln H = -Φ:")
print("  (∂Ψ)_M only non-zero for transverse M=yk:")
print("  (∂Ψ)_{yk} = (1/2)(H'/H)(yk/r)")
print()

# Construct the tensor ∂Ψ_M ∂Ψ_N
# Only non-zero for M=N=yk: (1/4)(H'/H)²(yk²/r²)
# Full contraction |∂Ψ|² = g^{MN} ∂Ψ_M ∂Ψ_N = sum_k g^{ykyk} (1/4)(H'/H)²(yk²/r²)
#                         = (1/H^{1/4})(1/4)(H'/H)² = H'²/(4H^{9/4})

dPsi_sq = hf.Hp**2 / (4 * hf.H**R(9,4))
print(f"  |∂Ψ|²(12d) = {dPsi_sq}")

# Now check: is the residual ℛ - T^{G₄}(λ) expressible as:
#   α ∂Ψ_M ∂Ψ_N + β |∂Ψ|² g_{MN}
# for some α, β?

print("\nResidual analysis: ℛ - T^{G₄}(λ) = α(∂Ψ)²_{MN} + β|∂Ψ|²g_{MN}?")
print("Testing with standard λ = 3/10:")
print()

lam_std = R(3,10)
a_sym, b_sym = sp.symbols('alpha beta')

for i12, name in [(0,'t'), (2,'z1'), (3,'z2'), (4,'y0')]:
    R12 = hf.substitute(sp.cancel(Ric_12[i12, i12]))
    ff = hf.substitute(sp.cancel(FF_D1[i12, i12]))
    s = hf.substitute(sp.cancel(S_D1))
    g = hf.substitute(sp.cancel(metric_12.matrix[i12, i12]))

    T_std = R(1,2) * (ff/6 - lam_std * s * g)

    resid = sp.cancel(R12 - T_std)

    # ∂Ψ²_{MN}: zero for worldvolume/z, (1/4)(H'/H)²(yk²/r²) for transverse
    if i12 >= 4:
        dPsi_MN = R(1,4) * hf.Hp**2 / hf.H**2 * coords_12[i12]**2 / hf.r**2
    else:
        dPsi_MN = sp.Integer(0)

    # |∂Ψ|² g_{MN}
    dPsi_g = dPsi_sq * g

    # resid = α × dPsi_MN + β × dPsi_g
    # For wv/z (dPsi_MN=0): resid = β × dPsi_g → β = resid/(dPsi_g)
    # For transverse: resid = α × dPsi_MN + β × dPsi_g

    if dPsi_MN == 0:
        beta_val = sp.cancel(resid / dPsi_g) if dPsi_g != 0 else 'N/A'
        print(f"  [{name}] (non-transverse): β = {beta_val}")
    else:
        # Will solve once we know β from non-transverse
        print(f"  [{name}] (transverse): resid = {resid}")
        print(f"    ∂Ψ²_{name} = {dPsi_MN}")
        print(f"    |∂Ψ|²g_{name} = {sp.cancel(dPsi_g)}")

# Let's solve explicitly
# From t: resid_t = β × dPsi_sq × g_t
# → β = resid_t / (dPsi_sq × g_t)
R12_t = hf.substitute(sp.cancel(Ric_12[0,0]))
ff_t = hf.substitute(sp.cancel(FF_D1[0,0]))
g_t = hf.substitute(sp.cancel(metric_12.matrix[0,0]))
resid_t = sp.cancel(R12_t - R(1,2)*(ff_t/6 - lam_std*s_d1_12*g_t))
beta_from_t = sp.cancel(resid_t / (dPsi_sq * g_t))

print(f"\n  β (from t block) = {beta_from_t}")

# From z1: resid_z1 = β × dPsi_sq × g_z1
R12_z1 = hf.substitute(sp.cancel(Ric_12[2,2]))
ff_z1 = hf.substitute(sp.cancel(FF_D1[2,2]))
g_z1 = hf.substitute(sp.cancel(metric_12.matrix[2,2]))
resid_z1 = sp.cancel(R12_z1 - R(1,2)*(ff_z1/6 - lam_std*s_d1_12*g_z1))
beta_from_z1 = sp.cancel(resid_z1 / (dPsi_sq * g_z1))

print(f"  β (from z1 block) = {beta_from_z1}")

# From z2:
R12_z2 = hf.substitute(sp.cancel(Ric_12[3,3]))
ff_z2 = hf.substitute(sp.cancel(FF_D1[3,3]))
g_z2 = hf.substitute(sp.cancel(metric_12.matrix[3,3]))
resid_z2 = sp.cancel(R12_z2 - R(1,2)*(ff_z2/6 - lam_std*s_d1_12*g_z2))
beta_from_z2 = sp.cancel(resid_z2 / (dPsi_sq * g_z2))

print(f"  β (from z2 block) = {beta_from_z2}")

# From y0 (isotropic, y0=0):
R12_y0 = hf.substitute(sp.cancel(Ric_12[4,4])).subs(y0, 0)
ff_y0 = hf.substitute(sp.cancel(FF_D1[4,4])).subs(y0, 0)
g_y0 = hf.substitute(sp.cancel(metric_12.matrix[4,4]))
resid_y0 = sp.cancel(R12_y0 - R(1,2)*(ff_y0/6 - lam_std*s_d1_12*g_y0))
beta_from_y0 = sp.cancel(resid_y0 / (dPsi_sq * g_y0))

print(f"  β (from y0_iso) = {beta_from_y0}")

print(f"\n  β consistent? t vs z1: {sp.cancel(beta_from_t - beta_from_z1) == 0}")
print(f"  β consistent? t vs z2: {sp.cancel(beta_from_t - beta_from_z2) == 0}")
print(f"  β consistent? t vs y0: {sp.cancel(beta_from_t - beta_from_y0) == 0}")

# If β is not consistent, maybe the ansatz needs modification.
# Try: resid = α ∂Ψ²_{MN} + β |∂Ψ|² g_{MN} + γ T_{MN}^{correction}

# Actually, let me try a different decomposition:
# ℛ_{MN} = c₁ FF_{MN}/3! + c₂ S g_{MN} + c₃ ∂Ψ_M∂Ψ_N + c₄ |∂Ψ|² g_{MN}
# 4 unknowns, need 4 equations (t, z1, z2, y0_iso)
# (y0_aniso constrains c₃ independently)

print("\n\nFull 4-parameter fit:")
print("  ℛ_{MN} = c₁ FF_{MN}/6 + c₂ S g_{MN} + c₃ (∂Ψ)_M(∂Ψ)_N + c₄ |∂Ψ|² g_{MN}")
c1, c2, c3, c4 = sp.symbols('c1 c2 c3 c4')

eqs = []
for i12, name in [(0,'t'), (2,'z1'), (3,'z2')]:
    R12_val = hf.substitute(sp.cancel(Ric_12[i12, i12]))
    ff_val = hf.substitute(sp.cancel(FF_D1[i12, i12]))
    g_val = hf.substitute(sp.cancel(metric_12.matrix[i12, i12]))

    # ∂Ψ²_{MN} = 0 for non-transverse
    rhs = c1 * ff_val/6 + c2 * s_d1_12 * g_val + c4 * dPsi_sq * g_val
    eq = sp.cancel(rhs - R12_val)
    eqs.append(eq)

# y0 isotropic (y0=0): ∂Ψ²_{y0} = 0 at y0=0
R12_y0_iso = hf.substitute(sp.cancel(Ric_12[4,4])).subs(y0, 0)
ff_y0_iso = hf.substitute(sp.cancel(FF_D1[4,4])).subs(y0, 0)
g_y0_val = hf.substitute(sp.cancel(metric_12.matrix[4,4]))
rhs_y0 = c1 * ff_y0_iso/6 + c2 * s_d1_12 * g_y0_val + c4 * dPsi_sq * g_y0_val
eq_y0 = sp.cancel(rhs_y0 - R12_y0_iso)
eqs.append(eq_y0)

print("  Solving t, z1, z2, y0_iso for (c₁, c₂, c₃, c₄)...")
print("  (c₃ only appears in y0_aniso, not in these equations)")

sol = sp.solve(eqs, [c1, c2, c4], dict=True)
if sol:
    for s in sol:
        print(f"\n  Solution:")
        for var in [c1, c2, c4]:
            print(f"    {var} = {sp.cancel(s.get(var, 'free'))}")

        # Recover the formula
        c1v = s[c1]
        c2v = s[c2]
        c4v = s[c4]

        print(f"\n  12d formula (non-transverse blocks):")
        print(f"    ℛ_MN = ({c1v}) FF/6 + ({c2v}) S g + ({c4v}) |∂Ψ|² g")

        # Check: does λ_eff = -c₂/c₁ match anything?
        lam_eff = sp.cancel(-c2v / c1v)
        print(f"    Effective λ = -c₂/c₁ = {lam_eff}")
        print(f"    c₁/2 = {sp.cancel(c1v/2)} (compare with standard 1/2)")

        # Now determine c₃ from y0 anisotropic
        R12_y0_full = hf.substitute(sp.cancel(Ric_12[4,4]))
        ff_y0_full = hf.substitute(sp.cancel(FF_D1[4,4]))

        # y0² coefficient
        R_aniso = sp.cancel(sp.diff(R12_y0_full, y0, 2) / 2)
        ff_aniso = sp.cancel(sp.diff(ff_y0_full, y0, 2) / 2)

        # ∂Ψ²_{y0y0} has y0² coefficient: (1/4) H'²/(H² r²) (from ∂Ψ_{y0} = (1/2)(H'/H)(y0/r))
        dPsi_aniso = R(1,4) * hf.Hp**2 / (hf.H**2 * hf.r**2)

        # |∂Ψ|²g has no y0-dependence, so no aniso contribution

        # R_aniso = c₁ ff_aniso/6 + c₃ dPsi_aniso
        c3_val = sp.cancel((R_aniso - c1v * ff_aniso/6) / dPsi_aniso)
        print(f"\n    c₃ (from y0 aniso) = {c3_val}")

        # Final formula
        print(f"\n  *** COMPLETE 12d EINSTEIN EQUATION ***")
        print(f"  ℛ_MN = {sp.cancel(c1v)} FF_MN/6 + ({sp.cancel(c2v)}) S g_MN")
        print(f"       + {sp.cancel(c3_val)} ∂Ψ_M ∂Ψ_N + ({sp.cancel(c4v)}) |∂Ψ|² g_MN")
        print(f"  where Ψ = (1/2) ln(g_z1/g_z2) = (1/2) ln H = -Φ")

        # Verify: does this match ALL blocks?
        print(f"\n  Verification:")
        for i12, name in [(0,'t'), (1,'x1'), (2,'z1'), (3,'z2'), (4,'y0'), (5,'y1')]:
            R12_val = hf.substitute(sp.cancel(Ric_12[i12, i12]))
            ff_val = hf.substitute(sp.cancel(FF_D1[i12, i12]))
            g_val = hf.substitute(sp.cancel(metric_12.matrix[i12, i12]))

            if i12 >= 4:
                coord = coords_12[i12]
                dPsi_MN_val = R(1,4) * hf.Hp**2 / hf.H**2 * coord**2 / hf.r**2
            else:
                dPsi_MN_val = sp.Integer(0)

            rhs = (c1v * ff_val/6 + c2v * s_d1_12 * g_val
                   + c3_val * dPsi_MN_val + c4v * dPsi_sq * g_val)
            rhs = sp.cancel(rhs)

            diff = sp.cancel(R12_val - rhs)
            status = "✓" if diff == 0 else f"✗ diff={diff}"
            print(f"    [{name}] {status}")

else:
    print("  No solution found!")

print("\n" + "="*70)
print("EXPERIMENT 13 COMPLETE")
print("="*70)
