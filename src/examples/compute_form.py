import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import sympy as sp
from sugra import HarmonicFunction, warped_product

# ── Parameters (edit these) ───────────────────────────────────────────────────

a = sp.Rational(-3, 4)
b = -sp.Rational(1, 3) * a #demanded for sensible solutions

d_wv = 2                          # worldvolume dimension
d_harm = 8                        # harmonic (transverse) dimension

# ── Setup Coordinates ─────────────────────────────────────────────────────────

wv_coords = list(sp.symbols(
    ' '.join(['t'] + [f'x{i}' for i in range(1, d_wv)]), real=True
))
harmonic_coords = list(sp.symbols(f'y0:{d_harm}', real=True))

coords = wv_coords + harmonic_coords
D = len(coords)

# ── Build Harmonic Function ──────────────────────────────────────────────────

hf = HarmonicFunction(transverse_coords=harmonic_coords)
H_func = sp.Function('H')(hf.r_expr)

# ── Build Metric ─────────────────────────────────────────────────────────────
#
#   ds^2 = H^a ds^2_{wv} + H^b ds^2_{harm}
#

metric = warped_product(
    warp_factors     = [H_func**a, H_func**b],
    block_dims       = [d_wv, d_harm],
    block_signatures = ['lorentzian', 'euclidean'],
    coordinates      = coords,
)

# -- Build 3-form field strength and dilaton for D1 ──────────────────────────────────────
#
#   (C_2)_{01} = H^{-1}
#   F_3 = dC_2  =>  (F_3)_{m,01} = -H^{-2} partial_m H * epsilon_{01}
#   e^\Phi = H^{1/2}
#

from sugra import FormField, exterior_derivative, form_norm_squared, form_contraction

C2 = FormField(rank=2, dim=D)
wv_idx = tuple(range(d_wv))              # (0, 1)
C2[wv_idx] = 1 / H_func                  # C_{01} = H^{-1}

F3 = exterior_derivative(C2, coords)

# -- Build dilaton ────────────────────────────────────────────────────────────
#
#   e^Phi = H^{1/2}  =>  Phi = (1/2) ln H
#

dilaton_power = sp.Rational(1, 2)         # exponent in e^Phi = H^p
Phi = dilaton_power * sp.log(H_func)      # Phi = p * ln H

# ── Display ─────────────────────────────────────────────────────────────────

print("F3 = dC2,  C_{01} = H^{-1}")
print()

for idx, val in F3.nonzero_components.items():
    names = ','.join(str(coords[i]) for i in idx)
    val_sub = hf.substitute(sp.cancel(val))
    print(f"  F3[{names}] = {val_sub}")

print()

F_sq = form_norm_squared(F3, metric)
F_sq = hf.substitute(sp.cancel(F_sq))
print(f"  |F3|^2 = {F_sq}")

# Dilaton gradient: partial_M Phi
print()
print(f"  Phi = {dilaton_power} * ln(H)")
print()

dPhi = [sp.diff(Phi, xi) for xi in coords]
for i, val in enumerate(dPhi):
    if val != 0:
        val_sub = hf.substitute(sp.cancel(val))
        print(f"  d_{{  {coords[i]}}} Phi = {val_sub}")

# ── Stress-energy tensor T_{MN} ──────────────────────────────────────────────
#
#   T_{MN} = 1/2 partial_M Phi partial_N Phi
#          + 1/2 e^{alpha*Phi} (|FF_{MN}| - (n-1)/(D-2) |F^2| g_{MN})
#
#   where |FF_{MN}| = 1/(n-1)! F_{M...} F_N^{...}
#         |F^2|     = 1/n!     F_{...}   F^{...}
#

alpha = 1#sp.Symbol('alpha')                # dilaton coupling
n = F3.rank                               # form rank
g = metric.matrix

# Dilaton prefactor
e_alpha_Phi = sp.exp(alpha * Phi)         # e^{alpha * Phi} = H^{alpha * p}

# Normalized form contraction: |FF_{MN}| = 1/(n-1)! F_{M P...} F_N^{P...}
FF_MN = form_contraction(F3, metric) / sp.factorial(n - 1)

# Normalized form norm squared: |F^2| = 1/n! F_{M...} F^{M...}
F_sq = form_norm_squared(F3, metric)

# Assemble T_{MN}
T = sp.zeros(D, D)
for M in range(D):
    for N in range(M, D):
        # Dilaton kinetic term
        dilaton_term = sp.Rational(1, 2) * dPhi[M] * dPhi[N]
        # Form field term
        form_term = sp.Rational(1, 2) * e_alpha_Phi * (
            FF_MN[M, N] - sp.Rational(n - 1, D - 2) * F_sq * g[M, N]
        )
        T[M, N] = dilaton_term + form_term
        if M != N:
            T[N, M] = T[M, N]

# ── Display T_{MN} ──────────────────────────────────────────────────────────

print()
print(f"T_{{MN}} = 1/2 d_M Phi d_N Phi + 1/(2(n-1)!) e^(alpha*Phi) "
      f"[FF_{{MN}} - (n-1)/(n(D-2)) F^2 g_{{MN}}]")
print(f"n = {n},  D = {D},  alpha = alpha,  Phi = {dilaton_power} ln(H)")
print()

for M in range(D):
    expr = hf.substitute(sp.cancel(T[M, M]))
    name = str(coords[M])
    if M < d_wv:
        label = "wv"
    else:
        label = "harm"
    print(f"  T[{name},{name}]  ({label}) = {expr}")
    print()