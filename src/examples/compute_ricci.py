"""
Ricci tensor for a general multi-block warped-product ansatz.

  ds^2_D = H^a ds^2_{wv} + H^c dz_1^2 + H^{-c} dz_2^2 + H^b ds^2_{harmonic}

H(r) harmonic in the transverse (harmonic) block.
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import sympy as sp
from sugra import HarmonicFunction, warped_product

# ── Parameters (edit these) ───────────────────────────────────────────────────

a = sp.Symbol('a')
b = sp.Symbol('b')
c = sp.Symbol('c')

d_wv = 2                          # worldvolume dimension
d_harm = 8                        # harmonic (transverse) dimension

# ── Setup Coordinates ─────────────────────────────────────────────────────────

wv_coords = list(sp.symbols(
    ' '.join(['t'] + [f'x{i}' for i in range(1, d_wv)]), real=True
))
z1, z2 = sp.symbols('z1 z2', real=True)
harmonic_coords = list(sp.symbols(f'y0:{d_harm}', real=True))

coords = wv_coords + [z1, z2] + harmonic_coords
D = len(coords)

# ── Build Harmonic Function ──────────────────────────────────────────────────

hf = HarmonicFunction(transverse_coords=harmonic_coords)
H_func = sp.Function('H')(hf.r_expr)

# ── Build Metric ─────────────────────────────────────────────────────────────
#
#   ds^2 = H^a ds^2_{wv} + H^c dz_1^2 + H^{-c} dz_2^2 + H^b ds^2_{harm}
#

metric = warped_product(
    warp_factors     = [H_func**a, H_func**c, H_func**(-c), H_func**b],
    block_dims       = [d_wv, 1, 1, d_harm],
    block_signatures = ['lorentzian', 'euclidean', 'euclidean', 'euclidean'],
    coordinates      = coords,
)

# ── Compute & display ────────────────────────────────────────────────────────

print(f"ds^2 = H^a ds^2_{d_wv} + H^c dz1^2 + H^(-c) dz2^2 + H^b ds^2_{d_harm}")
print(f"D = {D},  d_wv = {d_wv},  d_harm = {d_harm}")
print(f"H'' + {d_harm - 1}/r H' = 0")
print()

R = metric.ricci_tensor(simplify_func=sp.cancel)

for i in range(D):
    expr = hf.substitute(sp.cancel(R[i, i]))
    name = str(coords[i])
    if i < d_wv:
        label = "wv"
    elif i < d_wv + 2:
        label = "torus"
    else:
        label = "harm"
    print(f"  R[{name},{name}]  ({label}) = {expr}")
