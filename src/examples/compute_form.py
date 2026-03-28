"""
4-form flux F₄ for the M2-brane ansatz and its norm |F₄|².

  ds²₁₁ = H^{-2/3} ds²₃ + H^{1/3} ds²₈

Specify one non-vanishing component of the gauge potential C₃,
then obtain all F₄ = dC₃ components via exterior derivative.

  (C₃)_{012} = −H⁻¹
  F₄ = dC₃   ⟹   F_{012,m} = ∂_m(H⁻¹)

Convention: |F|² = (1/p!) F_{M₁…Mₚ} F^{M₁…Mₚ}.
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import sympy as sp
from sugra import (HarmonicFunction, warped_product,
                   FormField, exterior_derivative, form_norm_squared)

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

# ── Gauge potential: one non-vanishing component ─────────────────────────────

C3 = FormField(rank=3, dim=D)
C3[(0, 1, 2)] = -1 / H_func          # (C₃)_{012} = −H⁻¹

# ── Field strength: all components from exterior derivative ──────────────────

F4 = exterior_derivative(C3, coords)

# ── Compute & display ────────────────────────────────────────────────────────

print(f"ds^2 = H^({a}) ds^2_{d}  +  H^({b}) ds^2_{D_perp}")
print(f"C3:  C_{{012}} = -H^{{-1}}")
print(f"F4 = dC3")
print()

# Nonzero components
for idx, val in F4.nonzero_components.items():
    names = ','.join(str(coords[i]) for i in idx)
    val_sub = hf.substitute(sp.cancel(val))
    print(f"  F[{names}] = {val_sub}")

print()

# |F4|^2
F_sq = form_norm_squared(F4, metric)
F_sq = hf.substitute(sp.cancel(F_sq))
print(f"  |F4|^2 = {F_sq}")
