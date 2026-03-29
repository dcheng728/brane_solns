"""Experiment 13b: Verify the 12d equation derived from the F1-string uplift.

Key finding from exp13:
  ℛ_{MN}(12d) = T^{G₄}_{MN}(λ=1/4) for ALL 10d-direction components
  ℛ_{zi}(12d) = T^{G₄}_{zi}(λ=1/4) + (1/2)|∂Ψ|² g_{zi} for torus components

where:
  T^{G₄}(λ) = (1/2)[FF/3! - λ|G₄|²g]
  Ψ = (1/2)ln H = -Φ  (the dilaton, geometrized in 12d)
  λ = 1/4 = (p-1)/(D-2) with p=3, D=10 (the 10D value, not 12D!)

Physical interpretation:
  The 12d metric is NOT a genuine 12d supergravity — it's a 10d IIB repackaging.
  The 10d block uses the 10d-inherited λ, and the torus block gets an extra
  scalar kinetic contribution from the modular field.

  Equivalently: ℛ_{MN} = (1/2)[FF/3! - (1/4)|G₄|²g] + (1/2)|∂Ψ|² P^{torus}_{MN}
  where P^{torus} is the projection onto torus directions.

This experiment verifies the formula for ALL 12 diagonal components.
"""
import sys; sys.path.insert(0, 'src')
import sympy as sp
from sympy import Rational as R
from sugra import (HarmonicFunction, warped_product,
                   FormField, exterior_derivative,
                   form_contraction, form_norm_squared)

# --- Setup ---
wv_coords = list(sp.symbols('t x1', real=True))
z1s, z2s = sp.symbols('z1 z2', real=True)
harmonic_coords = list(sp.symbols('y0:8', real=True))
coords = wv_coords + [z1s, z2s] + harmonic_coords

hf = HarmonicFunction(transverse_coords=harmonic_coords)
H = sp.Function('H')(hf.r_expr)

metric = warped_product(
    warp_factors=[H**R(-3,4), H**R(1,2), H**R(-1,2), H**R(1,4)],
    block_dims=[2, 1, 1, 8],
    block_signatures=['lorentzian', 'euclidean', 'euclidean', 'euclidean'],
    coordinates=coords,
)

print("Computing Ricci tensor...")
Ric = metric.ricci_tensor(simplify_func=sp.cancel)
print("Done.\n")

# G₄ = H₃ ∧ dz₂  (the D1 form)
C_D1 = FormField(rank=3, dim=12)
C_D1[(0,1,3)] = 1/H
F_D1 = exterior_derivative(C_D1, coords)

FF = form_contraction(F_D1, metric)
S = form_norm_squared(F_D1, metric)

# Ψ = (1/2) ln H
# ∂Ψ_M = (1/2)(H'/H)(y_k/r) δ_{M,y_k}
# |∂Ψ|²_{12d} = ∑_k g^{y_k y_k} (∂Ψ_{y_k})² = (1/H^{1/4})(1/4)(H'/H)² × ∑(y²/r²)
#             = (1/(4H^{1/4}))(H'/H)² = H'²/(4H^{9/4})

dPsi_sq = hf.Hp**2 / (4 * hf.H**R(9,4))
S_val = hf.substitute(sp.cancel(S))

# ===================================================================
# Verify: ℛ_{MN} = (1/2)[FF_{MN}/3! - (1/4)|G₄|²g_{MN}] + (1/2)|∂Ψ|²P_{MN}
# where P_{MN} = g_{MN} for M=z1,z2; 0 otherwise
# ===================================================================

print("="*70)
print("VERIFICATION: ℛ_{MN} = T^{G₄}(λ=1/4) + (1/2)|∂Ψ|² P^{torus}")
print("="*70)

lam = R(1,4)
all_match = True

for i in range(12):
    name = str(coords[i])

    R_val = hf.substitute(sp.cancel(Ric[i,i]))
    ff_val = hf.substitute(sp.cancel(FF[i,i]))
    g_val = hf.substitute(sp.cancel(metric.matrix[i,i]))

    # T^{G₄}(λ=1/4)
    T_G4 = R(1,2) * (ff_val/6 - lam * S_val * g_val)
    T_G4 = sp.cancel(T_G4)

    # Torus projection
    is_torus = (i == 2 or i == 3)  # z1, z2
    if is_torus:
        T_torus = R(1,2) * dPsi_sq * g_val
    else:
        T_torus = sp.Integer(0)

    T_total = sp.cancel(T_G4 + T_torus)
    diff = sp.cancel(R_val - T_total)

    status = "✓" if diff == 0 else f"✗ diff={diff}"
    if diff != 0:
        all_match = False

    # Only print details for representative components
    if i <= 5 or diff != 0:
        print(f"\n  [{name}] (i={i}):")
        print(f"    ℛ      = {R_val}")
        print(f"    T^G₄   = {T_G4}")
        if is_torus:
            print(f"    T^torus= {sp.cancel(T_torus)}")
        print(f"    Total  = {T_total}")
        print(f"    {status}")

if all_match:
    print(f"\n{'='*70}")
    print("★ ALL 12 COMPONENTS VERIFIED ★")
    print()
    print("The 12d Einstein equation for the F1-string uplift is:")
    print()
    print("  ℛ_{MN} = (1/2)[G₄²_{MN}/3! - (1/4)|G₄|² 𝒢_{MN}]")
    print("         + (1/2)|∂Ψ|² P^{torus}_{MN}")
    print()
    print("where:")
    print("  G₄ = H₃ ∧ dz₂ = dB₂ ∧ dz₂  (F1-string NSNS 3-form ∧ torus)")
    print("  Ψ = (1/2) ln(g_{z1}/g_{z2}) = -Φ  (dilaton from torus ratio)")
    print("  P^{torus}_{MN} = g_{MN} for M,N ∈ {z1,z2}; 0 otherwise")
    print("  λ = 1/4 = (p-1)/(D-2)|_{p=3,D=10}  (inherited from 10d!)")
    print()
    print("Physical interpretation:")
    print("  The equation decomposes by 10d origin:")
    print("  • 10d directions: ℛ = T^{G₄}(λ=1/4)  [10d Einstein, dilaton absorbed by KK]")
    print("  • Torus directions: ℛ = T^{G₄}(λ=1/4) + scalar kinetic  [dilaton EOM]")
    print()
    print("  This is NOT a standard 12d SUGRA. The theory is 10d IIB repackaged.")
    print("  The trace coefficient λ=1/4 is the 10D value, not the 12D value 3/10.")
    print("  The torus projection breaks 12d diffeomorphism invariance to 10d × T².")
    print(f"{'='*70}")
else:
    print("\n--- MISMATCH: formula does not hold for all components ---")

# ===================================================================
# Cross-check with standard λ=3/10 (should NOT match)
# ===================================================================
print("\n\nCross-check: does λ=3/10 (standard 12d) work?")
lam_std = R(3,10)
for i, name in [(0,'t'), (2,'z1'), (3,'z2'), (4,'y0')]:
    R_val = hf.substitute(sp.cancel(Ric[i,i]))
    ff_val = hf.substitute(sp.cancel(FF[i,i]))
    g_val = hf.substitute(sp.cancel(metric.matrix[i,i]))
    T_std = R(1,2) * (ff_val/6 - lam_std * S_val * g_val)
    diff = sp.cancel(R_val - T_std)
    print(f"  [{name}] R-T(λ=3/10) = {diff}  {'✓' if diff == 0 else '✗'}")

# ===================================================================
# Rewrite: Can the formula be expressed as a SINGLE covariant equation?
# ===================================================================
print("\n\nAlternative form: ℛ = (1/2)[FF/3! - λ|G₄|²g] with λ block-dependent?")
for i, name in [(0,'t'), (2,'z1'), (3,'z2'), (4,'y0')]:
    R_val = hf.substitute(sp.cancel(Ric[i,i]))
    ff_val = hf.substitute(sp.cancel(FF[i,i]))
    g_val = hf.substitute(sp.cancel(metric.matrix[i,i]))

    # ℛ = (1/2)(ff/6 - λ S g) → λ = (ff/6 - 2ℛ)/(S g)
    lam_eff = sp.cancel((ff_val/6 - 2*R_val) / (S_val * g_val))
    print(f"  [{name}] effective λ = {lam_eff}")
