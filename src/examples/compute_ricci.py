"""
Ricci tensor for a general (d, D-d) warped-product brane ansatz.

  ds²_D = H^{a}(r) ds²_d  +  H^{b}(r) ds²_{D-d}

H(r) harmonic in the (D-d)-dimensional transverse space.
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import sympy as sp
from sugra import HarmonicFunction, warped_product

# ── Parameters (edit these) ───────────────────────────────────────────────────

D = 11                         # total dimension
d = 3                          # worldvolume dimension
a = sp.Rational(-2,3)         # worldvolume warp: H^{a}
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

# ── Compute & display ────────────────────────────────────────────────────────

print(f"ds² = H^({a}) ds²_{d}  +  H^({b}) ds²_{D_perp}")
print(f"H'' + {D_perp - 1}/r H' = 0")
print()

R = metric.ricci_tensor(simplify_func=sp.cancel)

for i in range(D):
    expr = hf.substitute(sp.cancel(R[i, i]))
    label = "wv" if i < d else "tr"
    print(f"  R[{coords[i]},{coords[i]}]  ({label}) = {expr}")
