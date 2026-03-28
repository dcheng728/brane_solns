"""
Brane verification: R_{MN} = T_{MN}

Uncomment the desired brane parameterization below.
Everything below the separator is fixed verification code.
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import sympy as sp
from sugra import HarmonicFunction, warped_product
from sugra import FormField, exterior_derivative, hodge_star, form_stress_energy

# ── Parameters ───────────────────────────────────────────────────────────────
#
# General Dp-brane (electric):  a = -(7-p)/8,  b = (p+1)/8
#   dilaton_power = (3-p)/4,  alpha = (3-p)/2,  C_{0..p} = H^{-1}
#
# F1:   d_wv=2, d_harm=8, a=-3/4, b=1/4, dilaton_power=-1/2, alpha=-1
# D1:   d_wv=2, d_harm=8, a=-3/4, b=1/4, dilaton_power= 1/2, alpha= 1
# D2:   d_wv=3, d_harm=7, a=-5/8, b=3/8, dilaton_power= 1/4, alpha= 1/2
# D3:   d_wv=4, d_harm=6, a=-1/2, b=1/2, dilaton_power= 0,   alpha= 0
# D4:   d_wv=5, d_harm=5, a=-3/8, b=5/8, dilaton_power=-1/4, alpha=-1/2
# D5:   d_wv=6, d_harm=4, a=-1/4, b=3/4, dilaton_power=-1/2, alpha=-1
# NS5:  d_wv=6, d_harm=4, a=-1/4, b=3/4, dilaton_power= 1/2, alpha= 1

# ── F1 ──
# d_wv = 2;  d_harm = 8
# a = sp.Rational(-3, 4);  b = sp.Rational(1, 4)
# dilaton_power = sp.Rational(-1, 2);  alpha = sp.Rational(-1, 1)

# ── D1 ──
# d_wv = 2;  d_harm = 8
# a = sp.Rational(-3, 4);  b = sp.Rational(1, 4)
# dilaton_power = sp.Rational(1, 2);  alpha = sp.Rational(1, 1)

# ── D2 ──
# d_wv = 3;  d_harm = 7
# a = sp.Rational(-5, 8);  b = sp.Rational(3, 8)
# dilaton_power = sp.Rational(1, 4);  alpha = sp.Rational(1, 2)

# ── D3 ──
# d_wv = 4;  d_harm = 6
# a = sp.Rational(-1, 2);  b = sp.Rational(1, 2)
# dilaton_power = sp.Rational(0, 1);  alpha = sp.Rational(0, 1)

# ── D4 ──
# d_wv = 5;  d_harm = 5
# a = sp.Rational(-3, 8);  b = sp.Rational(5, 8)
# dilaton_power = sp.Rational(-1, 4);  alpha = sp.Rational(-1, 2)

# ── D5 ──
# d_wv = 6;  d_harm = 4
# a = sp.Rational(-1, 4);  b = sp.Rational(3, 4)
# dilaton_power = sp.Rational(-1, 2);  alpha = sp.Rational(-1, 1)

# ── NS5 ──
d_wv = 6;  d_harm = 4
a = sp.Rational(-1, 4);  b = sp.Rational(3, 4)
dilaton_power = sp.Rational(1, 2);  alpha = sp.Rational(1, 1)

D = d_wv + d_harm

# ── Setup Coordinates ─────────────────────────────────────────────────────────

wv_coords = list(sp.symbols(
    ' '.join(['t'] + [f'x{i}' for i in range(1, d_wv)]), real=True
))
harmonic_coords = list(sp.symbols(f'y0:{d_harm}', real=True))

coords = wv_coords + harmonic_coords

# ── Build Harmonic Function ──────────────────────────────────────────────────

hf = HarmonicFunction(transverse_coords=harmonic_coords)
H_func = sp.Function('H')(hf.r_expr)

# ── Build Metric ─────────────────────────────────────────────────────────────

metric = warped_product(
    warp_factors     = [H_func**a, H_func**b],
    block_dims       = [d_wv, d_harm],
    block_signatures = ['lorentzian', 'euclidean'],
    coordinates      = coords,
)

# ── Build F ────────────────────────────────────────────────────────
#
# Electric ansatz: C_{0..p} = H^{-1}, F = dC  (F1, D1, D2, D4, D5, NS5)
# Self-dual:       F = (F_E + *F_E) / sqrt(2)  (D3 only)

C_E = FormField(rank=d_wv, dim=D)
C_E[tuple(range(d_wv))] = 1 / H_func
F_E = exterior_derivative(C_E, coords)

# ── F1, D1, D2 (electric) ──
# F = F_E

# ── D3 (self-dual) ──
# F = (F_E + hodge_star(F_E, metric)) / sp.sqrt(2)

# ── D4, D5, NS5 (magnetic) ──
F = hodge_star(F_E, metric)

# ── Build dilaton ────────────────────────────────────────────────────────────

Phi = dilaton_power * sp.log(H_func)      # Phi = (1/4) ln H

# ═══════════════════════════════════════════════════════════════════════════
# Verify R_{MN} = T_{MN}
# Expects: metric, F, Phi, alpha, coords, hf
# ═══════════════════════════════════════════════════════════════════════════
print(f"Using F that is {F.rank}-form in {metric.dim} dimensions.")
R = metric.ricci_tensor(simplify_func=sp.cancel)
T = form_stress_energy(F, metric, dilaton=Phi, dilaton_coupling=alpha)

for i in range(metric.dim):
    R_expr = hf.substitute(sp.cancel(R[i, i]))
    T_expr = hf.substitute(sp.cancel(T[i, i]))
    diff = sp.cancel(R_expr - T_expr)
    name = str(coords[i])
    status = "OK" if diff == 0 else "MISMATCH"
    print(f"  [{status}] {name}")
    print(f"    R = {R_expr}")
    print(f"    T = {T_expr}")
    if diff != 0:
        print(f"    R - T = {diff}")
