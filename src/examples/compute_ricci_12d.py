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

a = sp.Rational(-3, 4)
b = -sp.Rational(1, 3) * a #demanded for sensible solutions
c = sp.Rational(-1,2) #sp.Symbol('c')

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

# ── Helper: organize expression by H'/H powers ──────────────────────────────

def format_ricci(expr, hf):
    """Rewrite expr as sum of terms: (rational in a,b,c,...) * H'^n / H^m / r^k.

    Separates each term so the prefactor contains no H or H',
    and the H/H'/r dependence is a single monomial.
    """
    Hp, H, r = hf.Hp, hf.H, hf.r
    expr = sp.cancel(expr)
    # Decompose into additive terms
    expr_expanded = sp.Add.make_args(sp.expand(expr))
    # Group by (H' power, H power, r power) monomial
    monomial_map = {}
    for term in expr_expanded:
        # Extract powers of Hp, H, r from this term
        hp_exp = 0
        h_exp = 0
        r_exp = 0
        coeff = term
        # Extract H' power
        powers = sp.Mul.make_args(coeff)
        new_powers = []
        for factor in powers:
            base, exp = factor.as_base_exp()
            if base == Hp:
                hp_exp += exp
            elif base == H:
                h_exp += exp
            elif base == r:
                r_exp += exp
            else:
                new_powers.append(factor)
        coeff = sp.Mul(*new_powers) if new_powers else sp.S(1)
        key = (hp_exp, h_exp, r_exp)
        monomial_map[key] = monomial_map.get(key, sp.S(0)) + coeff

    terms = []
    for key in sorted(monomial_map.keys(), key=lambda k: (-k[0], -k[1], -k[2])):
        coeff = sp.cancel(monomial_map[key])
        if coeff == 0:
            continue
        hp_exp, h_exp, r_exp = key
        # Build the H'^n * H^m * r^k monomial string
        parts = []
        if hp_exp != 0:
            parts.append(f"H'^({hp_exp})" if hp_exp != 1 else "H'")
        if h_exp != 0:
            parts.append(f"H^({h_exp})" if h_exp != 1 else "H")
        if r_exp != 0:
            parts.append(f"r^({r_exp})" if r_exp != 1 else "r")
        monomial_str = ' * '.join(parts) if parts else ''
        import re
        coeff_str = str(coeff)
        # Replace x**n with x^(n) to avoid ambiguity like y7**2/8
        coeff_str = re.sub(r'\*\*(\w+)', r'^(\1)', coeff_str)
        if monomial_str:
            terms.append(f"({coeff_str}) * {monomial_str}")
        else:
            terms.append(f"({coeff_str})")
    return '\n      + '.join(terms) if terms else '0'

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
    print(f"  R[{name},{name}]  ({label}) =")
    print(f"    {format_ricci(expr, hf)}")
    print()


# -- Build 4-form field strength ──────────────────────────────────────────
#
#   C_{t, x1, z1} = H^{-1}
#   F_4 = dC_3  =>  F_{t, x1, z1, m} = partial_m(H^{-1})
#

from sugra import FormField, exterior_derivative, form_stress_energy

C3 = FormField(rank=3, dim=D)
C3[(0, 1, 2)] = 1 / H_func               # C_{t, x1, z1} = H^{-1}

F4 = exterior_derivative(C3, coords)

# ── Compute & display T_{MN} ────────────────────────────────────────────────

T = form_stress_energy(F4, metric)

print("T_{MN} from F4 = dC3,  C_{t,x1,z1} = H^{-1}  (no dilaton)")
print()

for i in range(D):
    expr = hf.substitute(sp.cancel(T[i, i]))
    name = str(coords[i])
    if i < d_wv:
        label = "wv"
    elif i < d_wv + 2:
        label = "torus"
    else:
        label = "harm"
    print(f"  T[{name},{name}]  ({label}) =")
    print(f"    {format_ricci(expr, hf)}")
    print()