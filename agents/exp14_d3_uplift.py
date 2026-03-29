"""Experiment 14: D3-brane 12d uplift verification.

The D3-brane is unique: Φ=0, self-dual F₅, SL(2,Z) singlet.

12d metric:
  ds²₁₂ = H^{-1/2} ds²_{1,3} + dz₁² + dz₂² + H^{1/2} ds²_6

Since Φ=0: τ₂=1, torus is FLAT, Ψ=0.
The self-dual F₅ has no z-legs and |F₅|²=0 (self-duality).

Expected 12d equation:
  ℛ_{MN} = (1/2) F₅²_{MN}/4!
  (no trace term, no torus projection term)

This tests:
  1. Formula with self-dual 5-form (different from 4-form in exp13)
  2. Flat torus case (trivial KK, but important cross-check)
  3. That |F₅|²_{12d} = 0 (self-duality preserved in embedding)
"""
import sys; sys.path.insert(0, 'src')
import sympy as sp
from sympy import Rational as R
from sugra import (HarmonicFunction, warped_product,
                   FormField, exterior_derivative, hodge_star,
                   form_contraction, form_norm_squared)

# ===================================================================
# Part A: 10d D3-brane verification
# ===================================================================
print("=" * 70)
print("PART A: D3-brane in 10d")
print("=" * 70)

wv_coords_10 = list(sp.symbols('t x1 x2 x3', real=True))
harmonic_coords_10 = list(sp.symbols('y0:6', real=True))
coords_10 = wv_coords_10 + harmonic_coords_10

hf_10 = HarmonicFunction(transverse_coords=harmonic_coords_10)
H10 = sp.Function('H')(hf_10.r_expr)

metric_10 = warped_product(
    warp_factors=[H10**R(-1, 2), H10**R(1, 2)],
    block_dims=[4, 6],
    block_signatures=['lorentzian', 'euclidean'],
    coordinates=coords_10,
)

print("Computing 10d Ricci tensor...")
Ric_10 = metric_10.ricci_tensor(simplify_func=sp.cancel)
print("Done.")

# Self-dual F₅
C4 = FormField(rank=4, dim=10)
C4[(0, 1, 2, 3)] = 1 / H10
F_E = exterior_derivative(C4, coords_10)
F5_10 = (F_E + hodge_star(F_E, metric_10)) / sp.sqrt(2)

FF_10 = form_contraction(F5_10, metric_10)
S_10 = form_norm_squared(F5_10, metric_10)
S_10_val = hf_10.substitute(sp.cancel(S_10))
print(f"\n  |F₅|²₁₀ = {S_10_val}  (should be 0 for self-dual)")

# D3 has Φ=0, so T_{mn} = (1/2)[FF/4! - (3/8)|F₅|²g]
# But with |F₅|²=0, just T = (1/2)FF/4!
print("\n  Checking R_{mn} = (1/2) FF_{mn}/4!:")
all_ok_10 = True
for i, name in [(0, 't'), (1, 'x1'), (4, 'y0')]:
    R_val = hf_10.substitute(sp.cancel(Ric_10[i, i]))
    ff_val = hf_10.substitute(sp.cancel(FF_10[i, i]))
    T_val = sp.cancel(ff_val / 48)  # 1/(2·4!) = 1/48
    diff = sp.cancel(R_val - T_val)
    status = "✓" if diff == 0 else f"✗ diff={diff}"
    if diff != 0:
        all_ok_10 = False
    print(f"    [{name}] R={R_val}, T={T_val}  {status}")

print(f"\n  10d D3 verification: {'PASS ✓' if all_ok_10 else 'FAIL ✗'}")

# ===================================================================
# Part B: 12d D3-brane uplift
# ===================================================================
print("\n" + "=" * 70)
print("PART B: D3-brane 12d uplift")
print("=" * 70)

wv_coords = list(sp.symbols('t x1 x2 x3', real=True))
z1s, z2s = sp.symbols('z1 z2', real=True)
harmonic_coords = list(sp.symbols('y0:6', real=True))
coords = wv_coords + [z1s, z2s] + harmonic_coords  # 12d

hf = HarmonicFunction(transverse_coords=harmonic_coords)
H = sp.Function('H')(hf.r_expr)

# 12d metric: H^{-1/2} ds²_{1,3} + dz₁² + dz₂² + H^{1/2} ds²_6
# Flat torus (Φ=0 → τ₂=1 → warp = H^0 = 1)
metric_12 = warped_product(
    warp_factors=[H**R(-1, 2), H**R(0, 1), H**R(0, 1), H**R(1, 2)],
    block_dims=[4, 1, 1, 6],
    block_signatures=['lorentzian', 'euclidean', 'euclidean', 'euclidean'],
    coordinates=coords,
)

print("Computing 12d Ricci tensor...")
Ric_12 = metric_12.ricci_tensor(simplify_func=sp.cancel)
print("Done.")

# Check that torus Ricci components vanish
for i, name in [(4, 'z1'), (5, 'z2')]:
    R_val = hf.substitute(sp.cancel(Ric_12[i, i]))
    print(f"\n  ℛ[{name}] = {R_val}  (should be 0 for flat torus)")

# Embed F₅ in 12d (same form, no z-legs)
# Indices shift: in 12d, x₁=1,x₂=2,x₃=3 still, but y₀=6,...,y₅=11
# The worldvolume indices (0,1,2,3) are unchanged
# The transverse indices shift by 2 (z1=4, z2=5 inserted)
C4_12 = FormField(rank=4, dim=12)
C4_12[(0, 1, 2, 3)] = 1 / H
F_E_12 = exterior_derivative(C4_12, coords)

# Hodge star in 12d gives a 8-form, not 5-form!
# Instead, compute 10d Hodge star and embed
# *₁₀F₅^E has indices in transverse space (y₀...y₅)
# In 12d coords, these are indices 6..11
F_star_10 = hodge_star(F_E, metric_10)  # 5-form in 10d

# Embed *₁₀F₅^E into 12d by shifting transverse indices by 2
F_star_12 = FormField(rank=5, dim=12)
for idx, val in F_star_10.nonzero_components.items():
    # Map 10d indices to 12d: 0-3 stay, 4-9 → 6-11
    new_idx = tuple(i if i < 4 else i + 2 for i in idx)
    F_star_12[new_idx] = val

F5_12 = (F_E_12 + F_star_12) / sp.sqrt(2)

print("\nComputing 12d form contractions...")
FF_12 = form_contraction(F5_12, metric_12)
S_12 = form_norm_squared(F5_12, metric_12)
S_12_val = hf.substitute(sp.cancel(S_12))
print(f"  |F₅|²₁₂ = {S_12_val}  (should be 0)")

# ===================================================================
# Part C: Verify ℛ_{MN} = (1/2) FF_{MN}/4!
# ===================================================================
print("\n" + "=" * 70)
print("PART C: Verify ℛ_{MN} = (1/2) FF_{MN}/4!")
print("=" * 70)

all_match = True
for i in range(12):
    name = str(coords[i])
    R_val = hf.substitute(sp.cancel(Ric_12[i, i]))
    ff_val = hf.substitute(sp.cancel(FF_12[i, i]))
    g_val = hf.substitute(sp.cancel(metric_12.matrix[i, i]))

    # T = (1/2)[FF/4! - λ|F₅|²g]  with |F₅|²=0
    T_val = sp.cancel(ff_val / 48)
    diff = sp.cancel(R_val - T_val)

    if diff != 0:
        all_match = False

    if i <= 7 or diff != 0:
        print(f"\n  [{name}] (i={i}):")
        print(f"    ℛ      = {R_val}")
        print(f"    T      = {T_val}")
        print(f"    diff   = {diff}  {'✓' if diff == 0 else '✗'}")

if all_match:
    print(f"\n{'=' * 70}")
    print("★ ALL 12 COMPONENTS VERIFIED ★")
    print()
    print("D3-brane 12d uplift equation:")
    print("  ℛ_{MN} = (1/2) F₅²_{MN}/4!")
    print()
    print("Key features:")
    print("  • Torus is FLAT (Φ=0 → τ₂=1)")
    print("  • |F₅|² = 0 (self-duality) → trace term irrelevant")
    print("  • No torus projection term needed (Ψ=0)")
    print("  • F₅ is SL(2,Z) singlet — no doublet structure")
    print("  • Formula is trivially λ-independent (trace vanishes)")
    print(f"{'=' * 70}")
else:
    print("\n--- MISMATCH: D3 formula does NOT hold ---")

    # Diagnose: check effective λ per block
    print("\nDiagnostics — effective λ per block:")
    for i, name in [(0, 't'), (4, 'z1'), (5, 'z2'), (6, 'y0')]:
        R_val = hf.substitute(sp.cancel(Ric_12[i, i]))
        ff_val = hf.substitute(sp.cancel(FF_12[i, i]))
        g_val = hf.substitute(sp.cancel(metric_12.matrix[i, i]))
        # ℛ = (1/2)(ff/4! - λ S g) → λ = (ff/4! - 2ℛ)/(S g)
        if S_12_val != 0:
            lam_eff = sp.cancel((ff_val / 24 - 2 * R_val) / (S_12_val * g_val))
            print(f"  [{name}] λ_eff = {lam_eff}")
        else:
            residual = sp.cancel(R_val - ff_val / 48)
            print(f"  [{name}] |F₅|²=0, residual R-FF/48 = {residual}")
