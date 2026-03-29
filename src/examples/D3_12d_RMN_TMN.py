"""
Brane verification: R_{MN} = T_{MN}

Uncomment the desired brane parameterization below.
Everything below the separator is fixed verification code.
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import sympy as sp
from sugra import HarmonicFunction, warped_product
from sugra import FormField, exterior_derivative, hodge_star, Verifier

# ── Parameters (edit these) ───────────────────────────────────────────────────

a = sp.Rational(-1, 2)
b = -sp.Rational(1, 2) * a #demanded for sensible solutions
c = 0

d_wv = 4                          # worldvolume dimension
d_harm = 6                        # harmonic (transverse) dimension

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

from sugra import FormField, exterior_derivative, form_stress_energy
d = sp.symbols('d', real=True)
C3 = FormField(rank=3, dim=D)
C3[(0, 1, 2)] = H_func ** (-1)               # C_{t, x1, z1} = H^{-1}

F4 = exterior_derivative(C3, coords)

C4 = FormField(rank=4, dim=D)
C4[(0, 1, 2, 3)] =  H_func ** (-1)
F5 = exterior_derivative(C4, coords)

F = F4

# ═══════════════════════════════════════════════════════════════════════════
# Verify R_{MN} = T_{MN}
# ═══════════════════════════════════════════════════════════════════════════

d,e = sp.symbols('d e', real=True)

soln = dict(
    metric = metric,
    forms  = [(F5, 0)],
    Phi    = sp.S(0),
    coords = coords,
    hf     = hf,
)

Verifier(soln).check()
