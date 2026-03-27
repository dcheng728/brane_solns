"""
Verify the analytical Ricci tensor formulas from A_Ricci_tensor.tex
against brute-force symbolic computation.

Metric ansatz:
  ds²_D = H^a ds²_d + H^b ds²_{D-d-2} + H^c dz1² + H^{-c} dz2²

with H harmonic in (D-d-2)-dimensional y-space, and constraint da+(D-d-4)b=0.

Claims to verify (lines 155-162):
  R_μν = (a/2) η_μν H^{a-b} (H')²/H²
  R_z1z1 = +(c/2) H^{-b+c} (H')²/H²
  R_z2z2 = -(c/2) H^{-b-c} (H')²/H²
  R_mn = (H')²/H² [ (b/2) δ_mn + (da(b-a)/4 - c²/2) y_m y_n / r² ]
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import sympy as sp
from sugra import HarmonicFunction, warped_product

# ── Parameters ────────────────────────────────────────────────────────────────
# Pick concrete values satisfying constraint da + (D-d-4)b = 0

D_tot = 11
d_wv = 3
D_perp = D_tot - d_wv - 2   # = 6

a_val = sp.Rational(2, 3)
# constraint: d*a + (D-d-4)*b = 0 => b = -da/(D-d-4)
b_val = -d_wv * a_val / (D_tot - d_wv - 4)
c_val = sp.Rational(1, 3)

print(f"D={D_tot}, d={d_wv}, D_perp={D_perp}")
print(f"a={a_val}, b={b_val}, c={c_val}")
print(f"Constraint check: da+(D-d-4)b = {d_wv*a_val + (D_tot-d_wv-4)*b_val}")
print()

# ── Coordinates ───────────────────────────────────────────────────────────────

wv_names = ['t'] + [f'x{i}' for i in range(1, d_wv)]
wv_coords = list(sp.symbols(' '.join(wv_names), real=True))

hf = HarmonicFunction(transverse_dim=D_perp)
y = hf.transverse_coords

z1, z2 = sp.symbols('z1 z2', real=True)

coords = wv_coords + y + [z1, z2]

H_func = sp.Function('H')(hf.r_expr)

# ── Metric (4 blocks) ────────────────────────────────────────────────────────

metric = warped_product(
    warp_factors     = [a_val, b_val, c_val, -c_val],
    block_dims       = [d_wv, D_perp, 1, 1],
    block_signatures = ['lorentzian', 'euclidean', 'euclidean', 'euclidean'],
    coordinates      = coords,
    H                = H_func,
)

# ── Compute Ricci ────────────────────────────────────────────────────────────

print("Computing Ricci tensor (this may take a minute)...")
R = metric.ricci_tensor(simplify_func=sp.cancel)
print("Done.\n")

# ── Compare against claimed formulas ─────────────────────────────────────────

H = hf.H
Hp = hf.Hp
r = hf.r

def check(label, idx, claimed):
    computed = hf.substitute(sp.cancel(R[idx, idx]))
    diff = sp.cancel(computed - claimed)
    status = "OK" if diff == 0 else "MISMATCH"
    print(f"  {label}: {status}")
    if diff != 0:
        print(f"    computed = {computed}")
        print(f"    claimed  = {claimed}")
        print(f"    diff     = {diff}")

# R_tt (worldvolume time): claimed = (a/2) * (-1) * H^{a-b} * Hp²/H²
#   η_tt = -1
check("R_tt", 0,
      sp.Rational(1, 2) * a_val * (-1) * H**(a_val - b_val) * Hp**2 / H**2)

# R_x1x1 (worldvolume space): claimed = (a/2) * (+1) * H^{a-b} * Hp²/H²
check("R_x1x1", 1,
      sp.Rational(1, 2) * a_val * H**(a_val - b_val) * Hp**2 / H**2)

# R_y0y0 (transverse): claimed = Hp²/H² * [b/2 + (da(b-a)/4 - c²/2) y0²/r²]
y0 = y[0]
coeff_ym = d_wv * a_val * (b_val - a_val) / 4 - c_val**2 / 2
check("R_y0y0", d_wv,
      Hp**2 / H**2 * (b_val / 2 + coeff_ym * y0**2 / r**2))

# R_y1y1
y1 = y[1]
check("R_y1y1", d_wv + 1,
      Hp**2 / H**2 * (b_val / 2 + coeff_ym * y1**2 / r**2))

# R_z1z1: claimed = +(c/2) H^{-b+c} Hp²/H²
z1_idx = d_wv + D_perp
check("R_z1z1", z1_idx,
      sp.Rational(1, 2) * c_val * H**(-b_val + c_val) * Hp**2 / H**2)

# R_z2z2: claimed = -(c/2) H^{-b-c} Hp²/H²
z2_idx = d_wv + D_perp + 1
check("R_z2z2", z2_idx,
      -sp.Rational(1, 2) * c_val * H**(-b_val - c_val) * Hp**2 / H**2)

print()
print(f"Transverse y_my_n coefficient = da(b-a)/4 - c²/2 = {coeff_ym}")
