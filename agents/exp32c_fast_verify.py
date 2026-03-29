"""Experiment 32c: Fast numerical verification of 12d formula.

Use minimal transverse coordinates (2d) for speed.
Formula: ℛ_{MN} = (1/2)[FF_{MN}/3! - (1/4)|G₄|²g_{MN}] - K · M̃_{MN}

For magnetic branes (NS5/D5) with 2 transverse coords, the Hodge dual is simpler:
*₂(dH) = ε_{ij} ∂_jH dy_i  (a 1-form, not 3-form), and G₄ = *₂dH ∧ dz.
But wait — NS5 needs codim 4 for the *₄dH construction.

Alternative: use the HarmonicFunction abstraction with composite function
for the magnetic case, or reduce to codim-2 with 2-form G₄ = dH ∧ dz (wrong rank).

Better approach: test the formula using the KNOWN Ricci and building block values
from previous experiments (which used HarmonicFunction abstraction correctly).
We have analytic expressions per block — just verify the algebraic relation.
"""
import sys; sys.path.insert(0, 'src')
import sympy as sp
from sympy import Rational as R, cancel, sqrt, symbols, Function
from sugra import HarmonicFunction

# ===================================================================
# APPROACH: Use the known per-block values from the notes and code
# to check the formula algebraically.
#
# For each brane, we need: ℛ_{ii}/g_{ii}, FF_{ii}/g_{ii}, |G₄|², S_{ii}/g_{ii}, K, M̃_{ii}
# These are expressible as c × H'^2 / H^a for appropriate c, a.
# ===================================================================

print("="*60)
print("Algebraic verification of 12d formula per brane")
print("="*60)

# Common: H' = dH/dr, expressed as H'^2/H^a factors
Hp2 = sp.Symbol("Hp2", positive=True)  # stands for H'^2
H = sp.Symbol("H", positive=True)
y2_r2 = sp.Symbol("y2_r2")  # stands for y_k^2/r^2 (for anisotropic terms)

def check_formula(brane_name, blocks):
    """blocks: list of (name, Ric/g, FF/g, |G4|^2, K, M_tilde/g, g_power)
    where everything is in units of Hp2 and H powers.
    Ric/g means ℛ_{ii} / g_{ii}, etc.
    Actually, let's just use the raw components."""
    print(f"\n--- {brane_name} ---")
    all_pass = True
    for bname, ric, ff, norm, K_val, m_tilde in blocks:
        T_form = R(1,2) * (ff/6 - R(1,4) * norm * 1)  # 1 = g/g (already ℛ/g etc)
        # Wait, this doesn't work because we need actual ℛ, FF, g values.
        pass
    # This algebraic approach is getting complicated. Let me just do a focused
    # numerical test with small dimension.

# ===================================================================
# BETTER APPROACH: Use HarmonicFunction abstraction but with only 2 transverse
# coords for the F1 case. For NS5, use 4 transverse coords but with
# the HarmonicFunction abstraction (composite function) for Ricci,
# and construct G₄ more carefully.
# ===================================================================

# Actually, the cleanest approach: compute everything for NS5 using the
# HarmonicFunction abstraction (like exp19 did), then check the formula
# symbolically with the abstracted H'/H^a expressions.

print("\n" + "="*60)
print("NS5: Using HarmonicFunction abstraction")
print("="*60)

wv6 = list(sp.symbols('t x1 x2 x3 x4 x5', real=True))
z1s, z2s = sp.symbols('z1 z2', real=True)
hc4 = list(sp.symbols('y0:4', real=True))
coords_ns5 = wv6 + [z1s, z2s] + hc4

hf = HarmonicFunction(transverse_coords=hc4)
H_ns5 = sp.Function('H')(hf.r_expr)

from sugra import warped_product, FormField, exterior_derivative, form_contraction, form_norm_squared
from itertools import combinations
from sympy.combinatorics import Permutation

def levi_civita_4d(i, j, k, l):
    perm = [i, j, k, l]
    if len(set(perm)) < 4:
        return 0
    sign = 1
    for a in range(4):
        for b in range(a+1, 4):
            if perm[a] > perm[b]:
                sign *= -1
    return sign

# NS5 EF metric
m_ns5 = warped_product(
    warp_factors=[H_ns5**R(-1,4), H_ns5**R(3,4), H_ns5**R(-1,2), H_ns5**R(1,2)],
    block_dims=[6, 4, 1, 1],
    block_signatures=['lorentzian', 'euclidean', 'euclidean', 'euclidean'],
    coordinates=coords_ns5,
)

# Torus
M_ns5 = sp.Matrix([[H_ns5**R(-1,2), 0], [0, H_ns5**R(1,2)]])
Minv_ns5 = sp.Matrix([[H_ns5**R(1,2), 0], [0, H_ns5**R(-1,2)]])

# G₄ = *₄dH ∧ dz₂
z2_idx = 7
trans_start = 8
G4_ns5 = FormField(rank=4, dim=12)

for i, j, k in combinations(range(4), 3):
    l = [x for x in range(4) if x not in (i,j,k)][0]
    eps = levi_civita_4d(i, j, k, l)
    if eps == 0:
        continue
    dH = sp.diff(H_ns5, hc4[l])  # chain rule: H'·y_l/r
    original = [trans_start+i, trans_start+j, trans_start+k, z2_idx]
    sorted_idx = sorted(original)
    perm_map = [sorted_idx.index(x) for x in original]
    sort_sign = Permutation(perm_map).signature()
    G4_ns5[tuple(sorted_idx)] = sp.Rational(eps) * sort_sign * dH

# Print G₄ components
print("G₄ components:")
for idx in sorted(G4_ns5.nonzero_components.keys()):
    val = hf.substitute(cancel(G4_ns5[idx]))
    print(f"  G₄{idx} = {val}")

# Compute Ricci
print("\nComputing Ricci tensor...")
Ric = m_ns5.ricci_tensor(simplify_func=cancel)
print("Done.")

# Compute FF and |G₄|²
print("Computing form contraction...")
FF = form_contraction(G4_ns5, m_ns5)
norm = form_norm_squared(G4_ns5, m_ns5)
print("Done.")

# Compute K
print("Computing K...")
D = 12
K_val = sp.Integer(0)
for i in range(D):
    if coords_ns5[i] not in H_ns5.free_symbols:
        continue
    dM = M_ns5.diff(coords_ns5[i])
    if dM == sp.zeros(2,2):
        continue
    dMinv = -Minv_ns5 * dM * Minv_ns5
    trace = (dMinv * dM).trace()
    g_inv = hf.substitute(cancel(m_ns5.inv_matrix[i,i]))
    K_val += g_inv * R(1,4) * trace
K_val = hf.substitute(cancel(K_val))
print(f"K = {K_val}")

# Check formula per block
block_indices = [(0, 't', None), (1, 'x1', None),
                 (6, 'z1', 0), (7, 'z2', 1),
                 (8, 'y0', None), (9, 'y1', None)]

y0 = hc4[0]

print(f"\nFormula: ℛ = (1/2)[FF/3! - (1/4)|G₄|²g] - K·M̃")
for idx, name, torus_ab in block_indices:
    ric = hf.substitute(cancel(Ric[idx, idx]))
    ff = hf.substitute(cancel(FF[idx, idx]))
    g = hf.substitute(cancel(m_ns5.matrix[idx, idx]))
    n = hf.substitute(cancel(norm))

    T_form = cancel(R(1,2) * (ff/6 - R(1,4)*n*g))
    T_form = hf.substitute(T_form)

    if torus_ab is not None:
        m_tilde = hf.substitute(cancel(M_ns5[torus_ab, torus_ab]))
    else:
        m_tilde = sp.Integer(0)

    KM = hf.substitute(cancel(K_val * m_tilde))
    T_total = cancel(T_form - KM)
    residual = hf.substitute(cancel(ric - T_total))

    print(f"\n[{name}]:")
    print(f"  ℛ      = {ric}")
    print(f"  T_form = {T_form}")
    print(f"  K·M̃    = {KM}")
    print(f"  residual = {residual}")

# ===================================================================
# Also test with OPPOSITE sign: ℛ = T_form + K·M̃
# ===================================================================
print(f"\n\nAlternative: ℛ = (1/2)[FF/3! - (1/4)|G₄|²g] + K·M̃")
for idx, name, torus_ab in block_indices:
    ric = hf.substitute(cancel(Ric[idx, idx]))
    ff = hf.substitute(cancel(FF[idx, idx]))
    g = hf.substitute(cancel(m_ns5.matrix[idx, idx]))
    n = hf.substitute(cancel(norm))

    T_form = cancel(R(1,2) * (ff/6 - R(1,4)*n*g))
    T_form = hf.substitute(T_form)

    if torus_ab is not None:
        m_tilde = hf.substitute(cancel(M_ns5[torus_ab, torus_ab]))
    else:
        m_tilde = sp.Integer(0)

    KM = hf.substitute(cancel(K_val * m_tilde))
    T_total = cancel(T_form + KM)
    residual = hf.substitute(cancel(ric - T_total))

    if residual != 0:
        print(f"  [{name}]: residual = {residual}")
    else:
        print(f"  [{name}]: ✓")

# ===================================================================
# Also compute effective λ per block
# ===================================================================
print(f"\n\nEffective λ per block:")
lam = sp.Symbol('lambda')
for idx, name, torus_ab in block_indices:
    ric = hf.substitute(cancel(Ric[idx, idx]))
    ff = hf.substitute(cancel(FF[idx, idx]))
    g = hf.substitute(cancel(m_ns5.matrix[idx, idx]))
    n = hf.substitute(cancel(norm))

    eq = cancel(ric - R(1,2)*(ff/6 - lam*n*g))
    sol = sp.solve(eq, lam)
    if sol:
        print(f"  [{name}]: λ = {[cancel(s) for s in sol]}")
    else:
        print(f"  [{name}]: no solution for λ")

print("\n" + "="*60)
print("DONE")
print("="*60)
