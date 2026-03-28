"""
4-form flux F₄ for the M2-brane ansatz and its norm |F₄|².

  ds²₁₁ = H^{-2/3} ds²₃ + H^{1/3} ds²₈
  (C₃)_{μ₀μ₁μ₂} = ε_{μ₀μ₁μ₂} (H⁻¹ - 1)
  F₄ = dC₃  ⟹  F_{μ₀μ₁μ₂,m} = ε_{μ₀μ₁μ₂} ∂_m(H⁻¹)

Convention: |F|² = (1/p!) F_{M₁…Mₚ} F^{M₁…Mₚ}.
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import sympy as sp
from sugra import HarmonicFunction, warped_product, FormField, form_norm_squared

# ── Parameters (edit these) ───────────────────────────────────────────────────

D = 11                         # total dimension
d = 3                          # worldvolume dimension
a = sp.Rational(-2, 3)         # worldvolume warp: H^{a}
b = sp.Rational(1, 3)          # transverse warp:  H^{b}

D_perp = D - d

# ── Setup ─────────────────────────────────────────────────────────────────────

wv_coords = list(sp.symbols(
    ' '.join(['t'] + [f'x{i}' for i in range(1, d)]), real=True
))
hf = HarmonicFunction(transverse_dim=D_perp)
y = hf.transverse_coords
coords = wv_coords + y

H_func = sp.Function('H')(hf.r_expr)

metric = warped_product(
    warp_factors     = [H_func ** a, H_func ** b],
    block_dims       = [d, D_perp],
    block_signatures = ['lorentzian', 'euclidean'],
    coordinates      = coords,
    H                = H_func,
)

# ── 4-form ansatz: F_{μ₀μ₁μ₂, m} = ε_{μ₀μ₁μ₂} ∂_m(H⁻¹) ─────────────────

wv_indices = list(range(d))            # (0, 1, 2)
H_inv = 1 / H_func

F4 = FormField(rank=4, dim=D)
for k, ym in enumerate(y):
    dH_inv = sp.diff(H_inv, ym)
    if dH_inv != 0:
        F4[tuple(wv_indices) + (d + k,)] = dH_inv

# ── Compute & display ────────────────────────────────────────────────────────

print(f"ds² = H^({a}) ds²_{d}  +  H^({b}) ds²_{D_perp}")
print(f"F_{{μ₀μ₁μ₂, m}} = ε_{{μ₀μ₁μ₂}} ∂_m(H⁻¹)")
print()

# Nonzero components
for idx, val in F4.nonzero_components.items():
    names = ','.join(str(coords[i]) for i in idx)
    val_sub = hf.substitute(sp.cancel(val))
    print(f"  F[{names}] = {val_sub}")

print()

# |F₄|²
F_sq = form_norm_squared(F4, metric)
F_sq = hf.substitute(sp.cancel(F_sq))
print(f"  |F₄|² = {F_sq}")
