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

# ── Build brane solution ────────────────────────────────────────────────────

def build_brane(d_wv, d_harm, a, b, dilaton_power, alpha, form_type='electric'):
    D = d_wv + d_harm
    wv_coords = list(sp.symbols(
        ' '.join(['t'] + [f'x{i}' for i in range(1, d_wv)]), real=True
    ))
    harmonic_coords = list(sp.symbols(f'y0:{d_harm}', real=True))
    coords = wv_coords + harmonic_coords

    hf = HarmonicFunction(transverse_coords=harmonic_coords)
    H_func = sp.Function('H')(hf.r_expr)

    metric = warped_product(
        warp_factors     = [H_func**a, H_func**b],
        block_dims       = [d_wv, d_harm],
        block_signatures = ['lorentzian', 'euclidean'],
        coordinates      = coords,
    )

    C_E = FormField(rank=d_wv, dim=D)
    C_E[tuple(range(d_wv))] = 1 / H_func
    F_E = exterior_derivative(C_E, coords)

    if form_type == 'electric':
        F = F_E
    elif form_type == 'magnetic':
        F = hodge_star(F_E, metric)
    elif form_type == 'self-dual':
        F = (F_E + hodge_star(F_E, metric)) / sp.sqrt(2)

    Phi = dilaton_power * sp.log(H_func)

    return dict(metric=metric, forms=[[F, alpha]], Phi=Phi, coords=coords, hf=hf)

# ── Brane catalog ───────────────────────────────────────────────────────────
#
# General Dp-brane:  a = -(7-p)/8,  b = (p+1)/8
#   dilaton_power = (3-p)/4,  alpha = (3-p)/2,  C_{0..p} = H^{-1}

R = sp.Rational
branes = {
    'F1':  build_brane(2, 8, R(-3,4), R(1,4),  R(-1,2), R(-1,1), 'electric'),
    'D1':  build_brane(2, 8, R(-3,4), R(1,4),  R(1,2),  R(1,1),  'electric'),
    'D2':  build_brane(3, 7, R(-5,8), R(3,8),  R(1,4),  R(1,2),  'electric'),
    'D3':  build_brane(4, 6, R(-1,2), R(1,2),  R(0,1),  R(0,1),  'self-dual'),
    'D4':  build_brane(5, 5, R(-3,8), R(5,8),  R(-1,4), R(-1,2), 'electric'),
    'D5':  build_brane(6, 4, R(-1,4), R(3,4),  R(-1,2), R(-1,1), 'magnetic'),
    'NS5': build_brane(6, 4, R(-1,4), R(3,4),  R(1,2),  R(1,1),  'magnetic'),
}

soln = branes['D4']

# ═══════════════════════════════════════════════════════════════════════════
# Verify R_{MN} = T_{MN}
# ═══════════════════════════════════════════════════════════════════════════

Verifier(soln).check()
