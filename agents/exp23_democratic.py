"""Experiment 23: Democratic formulation in 12d.

The "frame tension" from exp16-17: the Einstein-frame 12d metric is NOT
SL(2,Z) invariant, while the string-frame 12d metric IS, but the string-frame
equation has asymmetric dilaton couplings (e^{2Φ} for RR, none for NSNS).

Question: Can the democratic approach — using BOTH G₄ and *₁₂G₄ — yield a
manifestly SL(2,Z)-covariant 12d equation?

Key identity for p-forms in D dimensions:
  (*F)²_{MN}/(D-p-1)! = (-1)^{s+1}[F²_{MN}/(p-1)! - |F|²g_{MN}/p!]
where s = signature sign count.

For p=4, D=12 (Lorentzian, s=1):
  (*G₄)²_{MN}/7! = -[G₄²_{MN}/3! - |G₄|²g_{MN}/4!]

Democratic stress-energy (factor 1/2 for double-counting):
  T^{dem}_{MN} = (1/4)[G₄²_{MN}/3! + (*G₄)²_{MN}/7!]
               = (1/4)[G₄²_{MN}/3! - G₄²_{MN}/3! + |G₄|²g_{MN}/4!]
               = (1/4)|G₄|²g_{MN}/24

This is pure trace! The democratic approach for a non-self-dual form gives
T ∝ |F|²g — which cannot match a non-trivial Ricci tensor.

Plan:
A. Verify the (*F)²/(D-p-1)! identity numerically for F1's G₄
B. Compute the democratic stress-energy explicitly
C. Check if any modified democratic formula (with M-dependent factors) works
D. Conclusions about the frame tension
"""
import sys; sys.path.insert(0, 'src')
import sympy as sp
from sympy import Rational as R, ln, sqrt, cancel, simplify
from sugra import (HarmonicFunction, Metric, warped_product,
                   FormField, exterior_derivative, hodge_star,
                   form_contraction, form_norm_squared,
                   form_stress_energy)

print("=" * 70)
print("EXP 23: Democratic formulation in 12d")
print("=" * 70)

# ===================================================================
# Setup: F1-string in 12d (Einstein frame, from exp13)
# ===================================================================
wv_coords = list(sp.symbols('t x1', real=True))
z1s, z2s = sp.symbols('z1 z2', real=True)
trans_coords = list(sp.symbols('y0:8', real=True))
coords = wv_coords + [z1s, z2s] + trans_coords
D = 12

hf = HarmonicFunction(transverse_coords=trans_coords)
H = sp.Function('H')(hf.r_expr)

# 12d Einstein-frame metric
g_ef = sp.zeros(D, D)
g_ef[0, 0] = -H**R(-3, 4)   # t
g_ef[1, 1] = H**R(-3, 4)    # x1
g_ef[2, 2] = H**R(1, 2)     # z1
g_ef[3, 3] = H**R(-1, 2)    # z2
for k in range(8):
    g_ef[4 + k, 4 + k] = H**R(1, 4)  # transverse

metric_ef = Metric(g_ef, coords)

# G₄ = dC₃ where C₃_{t,x1,z2} = H^{-1} (the D1 potential in 12d)
C3 = FormField(rank=3, dim=D)
C3[(0, 1, 3)] = 1/H  # C_{t,x1,z2} = H^{-1}
G4 = exterior_derivative(C3, coords)

print("\nG₄ non-zero components:")
for idx, val in sorted(G4.nonzero_components.items()):
    names = [str(coords[i]) for i in idx]
    print(f"  G₄[{','.join(names)}] = {sp.cancel(val)}")

# ===================================================================
# Part A: Verify (*F)² identity
# ===================================================================
print("\n" + "=" * 70)
print("PART A: Hodge dual identity verification")
print("=" * 70)

print("\nComputing *₁₂G₄ (8-form)...")
star_G4 = hodge_star(G4, metric_ef, signature=-1)

print("Non-zero *G₄ components (first few):")
count = 0
for idx, val in sorted(star_G4.nonzero_components.items()):
    if count < 5:
        names = [str(coords[i]) for i in idx]
        print(f"  *G₄[{','.join(names)}] = {cancel(val)}")
    count += 1
print(f"  ... ({count} total non-zero components)")

# Compute G₄² and (*G₄)² contractions
print("\nComputing G₄² contractions...")
FF_G4 = form_contraction(G4, metric_ef)
S_G4 = form_norm_squared(G4, metric_ef)

print("\nComputing (*G₄)² contractions...")
FF_star = form_contraction(star_G4, metric_ef)
S_star = form_norm_squared(star_G4, metric_ef)

# The identity: (*F)²_{MN}/(D-p-1)! = -[F²_{MN}/(p-1)! - |F|²g_{MN}/p!]
# i.e., (*G₄)²_{MN}/7! = -[G₄²_{MN}/3! - |G₄|²g_{MN}/4!]
#
# Rearranging: G₄²_{MN}/3! + (*G₄)²_{MN}/7! = |G₄|²g_{MN}/4!

print("\n--- Verifying identity: G₄²/3! + (*G₄)²/7! = |G₄|² g/4! ---")

fac_G4 = 1  # 3! = 6, already divided in form_contraction? Need to check.
# form_contraction returns FF_{MN} = F_{MA...}F_N^{A...} (raw, no factorial)
# So FF/3! means dividing by (p-1)! = 3! = 6
# form_norm_squared returns |F|² = F_{A...}F^{A...} (raw, no factorial)
# So |F|²/4! means dividing by p! = 4! = 24

from sympy import factorial

p = 4
Dp = D - p  # = 8

blocks = {
    't':  (0, 0),
    'x1': (1, 1),
    'z1': (2, 2),
    'z2': (3, 3),
    'y0': (4, 4),
}

print(f"\n|G₄|² = {cancel(S_G4)}")
print(f"|*G₄|² = {cancel(S_star)}")
print(f"|G₄|²/p! = {cancel(S_G4 / factorial(p))}")

for name, (i, j) in blocks.items():
    # LHS: FF_{G4}/3! + FF_{*G4}/7!
    lhs = FF_G4[i, j] / factorial(p - 1) + FF_star[i, j] / factorial(Dp - 1)
    # RHS: |G₄|² g_{ij} / 4!
    rhs = S_G4 * g_ef[i, j] / factorial(p)
    diff = cancel(lhs - rhs)
    status = "✓" if diff == 0 else f"✗ (diff = {diff})"
    print(f"  [{name}] LHS = {cancel(lhs)}, RHS = {cancel(rhs)} → {status}")

# ===================================================================
# Part B: Democratic stress-energy
# ===================================================================
print("\n" + "=" * 70)
print("PART B: Democratic stress-energy")
print("=" * 70)

# If the identity holds, then:
# T^{dem}_{MN} = (1/4)[FF_{G4}/3! + FF_{*G4}/7!] = (1/4)|G₄|²g/4!
# This is pure trace: T^{dem} ∝ g_{MN}

print("\n12d Ricci tensor (needed for comparison)...")
Ric_12 = metric_ef.ricci_tensor(simplify_func=cancel)

# Check: is ℛ_{MN} ∝ g_{MN}? (It shouldn't be, so democratic T^{dem} fails.)
print("\nℛ_{MN}/g_{MN} by block:")
for name, (i, j) in blocks.items():
    ric_over_g = cancel(Ric_12[i, j] / g_ef[i, j])
    tdem_over_g = cancel(S_G4 / (4 * factorial(p)))
    print(f"  [{name}] ℛ/g = {ric_over_g},  T^dem/g = {tdem_over_g}")

print("\n★ CONCLUSION: Democratic T^{dem} for non-self-dual G₄ is pure trace,")
print("  which CANNOT match a non-trivial Ricci tensor (ℛ/g varies by block).")

# ===================================================================
# Part C: Modified democratic with M-dependent weights
# ===================================================================
print("\n" + "=" * 70)
print("PART C: Modified democratic with weighted G₄ and *G₄")
print("=" * 70)

# Could we use different weights: α·G₄² + β·(*G₄)²?
# T_{MN} = (1/2)[α·FF_{G4}/(p-1)! + β·FF_{*G4}/(D-p-1)!]
#
# For this to work, we need α, β such that T matches ℛ for all blocks.
# Since (*G₄)²/(D-p-1)! = -[G₄²/(p-1)! - |G₄|²g/p!]:
# T = (1/2)[(α-β)·G₄²/(p-1)! + β·|G₄|²g/p!]
# This is T = (α-β)/2 · form_part + β/2 · trace_part
# = standard formula with effective coefficient (α-β)/2 and λ_eff = β/(α-β)/p ...
#
# So adding the dual form just reparametrizes (coeff, λ). No new information.

# Let's verify: for the STANDARD formula T = (1/2)[FF/3! - λ|G₄|²g],
# what λ makes it work for EF?

print("\nStandard formula: T_{MN} = c·[FF_{MN}/3! - λ|G₄|²g_{MN}]")
print("Solving for (c, λ) per block:")

for name, (i, j) in blocks.items():
    ff_val = cancel(FF_G4[i, j] / factorial(p - 1))
    s_g = cancel(S_G4 * g_ef[i, j])
    ric = cancel(Ric_12[i, j])
    print(f"  [{name}] FF/3! = {ff_val}")
    print(f"  [{name}] |G₄|²g = {s_g}")
    print(f"  [{name}] ℛ     = {ric}")
    if ff_val != 0:
        # c·ff - c·λ·|G₄|²g = ℛ → two unknowns (c, c·λ)
        print(f"  [{name}] ℛ/(FF/3!) = {cancel(ric / ff_val)}")
    print()

# ===================================================================
# Part D: String-frame democratic formulation
# ===================================================================
print("=" * 70)
print("PART D: String-frame analysis")
print("=" * 70)

# Build 12d string-frame metric: g^S + M dz²
# g^S = H^{-1} ds²_{1,1} + ds²_8
g_sf = sp.zeros(D, D)
g_sf[0, 0] = -H**(-1)       # t (string frame)
g_sf[1, 1] = H**(-1)        # x1
g_sf[2, 2] = H**R(1, 2)     # z1 (torus, from M)
g_sf[3, 3] = H**R(-1, 2)    # z2
for k in range(8):
    g_sf[4 + k, 4 + k] = sp.Integer(1)  # transverse (flat in SF)

metric_sf = Metric(g_sf, coords)

print("\nComputing string-frame 12d Ricci...")
Ric_sf = metric_sf.ricci_tensor(simplify_func=cancel)

print("\nComputing G₄ contractions in string frame...")
FF_G4_sf = form_contraction(G4, metric_sf)
S_G4_sf = form_norm_squared(G4, metric_sf)

print(f"\n|G₄|²_SF = {cancel(S_G4_sf)}")

print("\nString-frame Ricci and G₄ terms by block:")
for name, (i, j) in blocks.items():
    ric = cancel(Ric_sf[i, j])
    ff = cancel(FF_G4_sf[i, j] / factorial(p - 1))
    sg = cancel(S_G4_sf * g_sf[i, j])
    ric_g = cancel(ric / g_sf[i, j]) if g_sf[i, j] != 0 else 'N/A'
    ff_g = cancel(ff / g_sf[i, j]) if g_sf[i, j] != 0 else 'N/A'
    print(f"  [{name:3s}] ℛ/g = {ric_g}")
    print(f"        FF/(3!g) = {ff_g}")
    print(f"        |G₄|²   = {cancel(S_G4_sf)}")

# Check: does any single (c, λ) work in string frame?
print("\n--- Can c·[FF/3! - λ|G₄|²g] = ℛ in string frame? ---")

# Extract the coefficient ratios
# For blocks where FF ≠ 0: c·FF/3! - c·λ·|G₄|²g = ℛ
# For blocks where FF = 0: -c·λ·|G₄|²g = ℛ → c·λ = -ℛ/(|G₄|²g)

for name, (i, j) in blocks.items():
    ff = cancel(FF_G4_sf[i, j] / factorial(p - 1))
    sg = cancel(S_G4_sf * g_sf[i, j])
    ric = cancel(Ric_sf[i, j])
    if ff == 0 and sg != 0:
        cl = cancel(-ric / sg)
        print(f"  [{name}] FF=0 → c·λ = {cl}")
    elif ff != 0:
        print(f"  [{name}] FF≠0 → c - c·λ·(|G|²g/FF_normed) = ℛ/FF_normed")
        print(f"           ℛ/FF = {cancel(ric / ff)}")

# ===================================================================
# Part E: The SL(2,Z) covariance of G₄ contraction in SF
# ===================================================================
print("\n" + "=" * 70)
print("PART E: G₄ contraction decomposition in string frame")
print("=" * 70)

# In SF: g_{12d} = g^S ⊕ M, with g^S SL(2,Z)-invariant.
# (G₄²)^{SF}_{mn}/3! = (M⁻¹)^{ab} (F₃^a · F₃^b)^{SF}_{mn}/2!
# For F1: only a=b=2 contributes: (M⁻¹)^{22} = e^{-Φ} = H^{1/2}
# So (G₄²)^{SF}/3! = H^{1/2} · (H₃²)^{SF}/2!
#
# The 10d SF equation: R^S + 2∇∂Φ = (1/4)(H₃²)^{SF}/2!
# → (H₃²)^{SF}/2! = 4(R^S + 2∇∂Φ)
# → (G₄²)^{SF}/3! = 4H^{1/2}(R^S + 2∇∂Φ) = 4e^{-Φ}(R^S + 2∇∂Φ)
#
# The e^{-Φ} = (M⁻¹)^{22} = (M⁻¹)^{ab}c_a c_b is SL(2,Z) invariant!
# But only for a specific charge vector. For general (p,q) it's different.

# Let's verify: (G₄²)^{SF}/3! = H^{1/2} · (H₃²)^{SF}/2!
# To get (H₃²)^{SF}, compute it in 10d SF metric

# 10d SF metric
coords_10 = wv_coords + trans_coords
g_sf_10 = sp.zeros(10, 10)
g_sf_10[0, 0] = -H**(-1)
g_sf_10[1, 1] = H**(-1)
for k in range(8):
    g_sf_10[2 + k, 2 + k] = sp.Integer(1)

metric_sf_10 = Metric(g_sf_10, coords_10)

# H₃ = dB₂ in 10d
B2_10 = FormField(rank=2, dim=10)
B2_10[(0, 1)] = 1/H  # B_{t,x1} = H^{-1}
H3_10 = exterior_derivative(B2_10, coords_10)

print("\nComputing H₃² in 10d string frame...")
FF_H3_sf = form_contraction(H3_10, metric_sf_10)
S_H3_sf = form_norm_squared(H3_10, metric_sf_10)
print(f"|H₃|²_SF = {cancel(S_H3_sf)}")

# Map 10d blocks to 12d: 10d index i → 12d index (i if i<2, i+2 if i>=2)
def map_10_to_12(i):
    return i if i < 2 else i + 2

print("\nVerifying: (G₄²)^{SF}_{mn}/3! = H^{1/2} · (H₃²)^{SF}_{mn}/2!")
blocks_10d = {
    't':  (0, 0),
    'x1': (1, 1),
    'y0': (2, 2),
}

for name, (i10, j10) in blocks_10d.items():
    i12 = map_10_to_12(i10)
    j12 = map_10_to_12(j10)
    g4_ff = cancel(FF_G4_sf[i12, j12] / factorial(3))  # /3!
    h3_ff = cancel(H**R(1, 2) * FF_H3_sf[i10, j10] / factorial(2))  # H^{1/2} × /2!
    diff = cancel(g4_ff - h3_ff)
    status = "✓" if diff == 0 else f"✗ (diff = {diff})"
    print(f"  [{name}]  G₄²/3! = {g4_ff}")
    print(f"        H^½·H₃²/2! = {h3_ff} → {status}")

# Also check norms
norm_check = cancel(S_G4_sf - H**R(1, 2) * S_H3_sf)
print(f"\n|G₄|²_SF - H^½|H₃|²_SF = {norm_check}  {'✓' if norm_check == 0 else '✗'}")

print("\n" + "=" * 70)
print("CONCLUSIONS")
print("=" * 70)
print("""
1. DEMOCRATIC FORMULA (pure form):
   The democratic stress-energy for a non-self-dual p-form reduces to
   T^{dem}_{MN} ∝ |F|²g_{MN} (pure trace), which cannot match a
   non-trivial Ricci tensor. This is USELESS for sourcing gravity.

2. WEIGHTED DEMOCRATIC:
   Using different weights α·G₄² + β·(*G₄)² just reparametrizes the
   standard (coefficient, λ) formula. No new information beyond exp10-12.

3. STRING-FRAME G₄ FACTORIZATION:
   (G₄²)^{SF}/3! = (c^T M⁻¹ c) · (F₃²)^{SF}/2!
   The factor c^T M⁻¹ c = e^{-Φ} (for F1) is SL(2,Z) invariant for
   a given charge vector, but the 10d SF equation needs (F₃²) without
   this factor. Extracting it requires knowing c and M separately.

4. FRAME TENSION IS FUNDAMENTAL:
   The democratic approach does NOT resolve the frame tension. The tension
   arises because:
   - Einstein frame: g^E changes under SL(2,Z) → 12d metric not invariant
   - String frame: form couplings are asymmetric (1 for NSNS, e^{2Φ} for RR)
   - Democratic formulation addresses field CONTENT (using F_p and *F_p),
     not the frame choice. It doesn't help.

   The 12d formulation is best understood as a KK decomposition:
   ℛ_{mn}(12d) = R^S_{mn}(10d) + (1/4)Tr(∂M⁻¹·∂M)
   where BOTH sides are independently SL(2,Z) invariant.
""")
