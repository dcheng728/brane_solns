"""Experiment 32e: Verify the UNIVERSAL 12d formula.

Formula:
    ℛ_{MN} = (1/2)[FF_{MN}/3! - (1/4)|G₄|²g_{MN}] - (1/8)|G₄|²Π^{torus}_{MN}

where Π^{torus}_{MN} = g_{MN} for torus directions, 0 otherwise.

Equivalently: ℛ_{MN} = (1/2)[FF_{MN}/3! - λ_eff|G₄|²g_{MN}]
with λ_eff = 1/4 (10d) and 1/2 (torus).

Test on: F1, D1, NS5, D5, D3.
"""
import sys; sys.path.insert(0, 'src')
import sympy as sp
from sympy import Rational as R, cancel
from sugra import (HarmonicFunction, warped_product, FormField,
                   exterior_derivative, form_contraction, form_norm_squared)
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


def verify_universal(name, coords, hf, metric, G4, block_indices, is_torus):
    """Verify ℛ = (1/2)[FF/3! - λ_eff|G₄|²g] with λ=1/4(10d), 1/2(torus)."""
    print(f"\n{'='*60}")
    print(f"  {name}: Universal formula verification")
    print(f"{'='*60}")

    print("  Computing Ricci...")
    Ric = metric.ricci_tensor(simplify_func=cancel)

    if G4 is not None:
        print("  Computing form data...")
        FF = form_contraction(G4, metric)
        norm = form_norm_squared(G4, metric)
        norm_val = hf.substitute(cancel(norm))
        print(f"  |G₄|² = {norm_val}")
    else:
        FF = None
        norm_val = sp.Integer(0)

    all_pass = True
    for idx, label in block_indices:
        ric = hf.substitute(cancel(Ric[idx, idx]))
        g = hf.substitute(cancel(metric.matrix[idx, idx]))

        if FF is not None:
            ff = hf.substitute(cancel(FF[idx, idx]))
        else:
            ff = sp.Integer(0)

        lam = R(1,2) if is_torus(idx) else R(1,4)
        T = cancel(R(1,2) * (ff/6 - lam * norm_val * g))
        residual = cancel(ric - T)

        status = "✓" if residual == 0 else f"✗ res={residual}"
        print(f"  [{label:>4}] (λ={lam}): {status}")

        if residual != 0:
            all_pass = False

    result = 'ALL PASS ✓' if all_pass else 'SOME FAIL ✗'
    print(f"\n  {result}")
    return all_pass


# ===================================================================
# F1-BRANE (codim 8)
# ===================================================================
print("#"*60)
print("# F1-BRANE")
print("#"*60)

wv2 = list(sp.symbols('t x1', real=True))
z1s, z2s = sp.symbols('z1 z2', real=True)
hc8 = list(sp.symbols('y0:8', real=True))
coords_f1 = wv2 + [z1s, z2s] + hc8

hf_f1 = HarmonicFunction(transverse_coords=hc8)
H_f1 = sp.Function('H')(hf_f1.r_expr)

# F1: blocks [wv(2), z1(1), z2(1), trans(8)] matching coord order
m_f1 = warped_product(
    warp_factors=[H_f1**R(-3,4), H_f1**R(1,2), H_f1**R(-1,2), H_f1**R(1,4)],
    block_dims=[2, 1, 1, 8],
    block_signatures=['lorentzian', 'euclidean', 'euclidean', 'euclidean'],
    coordinates=coords_f1,
)

C3_f1 = FormField(rank=3, dim=12)
C3_f1[(0, 1, 3)] = 1/H_f1  # C_{t,x1,z2}
G4_f1 = exterior_derivative(C3_f1, coords_f1)

blocks_f1 = [(0,'t'), (1,'x1'), (2,'z1'), (3,'z2'), (4,'y0'), (5,'y1')]
is_torus_f1 = lambda idx: idx in (2, 3)

verify_universal("F1", coords_f1, hf_f1, m_f1, G4_f1, blocks_f1, is_torus_f1)


# ===================================================================
# D1-BRANE (codim 8)
# ===================================================================
print("\n" + "#"*60)
print("# D1-BRANE")
print("#"*60)

# D1: z₁↔z₂ swap vs F1
m_d1 = warped_product(
    warp_factors=[H_f1**R(-3,4), H_f1**R(-1,2), H_f1**R(1,2), H_f1**R(1,4)],
    block_dims=[2, 1, 1, 8],
    block_signatures=['lorentzian', 'euclidean', 'euclidean', 'euclidean'],
    coordinates=coords_f1,
)

C3_d1 = FormField(rank=3, dim=12)
C3_d1[(0, 1, 2)] = 1/H_f1  # C_{t,x1,z1}
G4_d1 = exterior_derivative(C3_d1, coords_f1)

verify_universal("D1", coords_f1, hf_f1, m_d1, G4_d1, blocks_f1, is_torus_f1)


# ===================================================================
# NS5-BRANE (codim 4)
# ===================================================================
print("\n" + "#"*60)
print("# NS5-BRANE")
print("#"*60)

wv6 = list(sp.symbols('t x1 x2 x3 x4 x5', real=True))
hc4 = list(sp.symbols('y0:4', real=True))
coords_ns5 = wv6 + [z1s, z2s] + hc4

hf_ns5 = HarmonicFunction(transverse_coords=hc4)
H_ns5 = sp.Function('H')(hf_ns5.r_expr)

# NS5: blocks [wv(6), z1(1), z2(1), trans(4)] matching coord order
m_ns5 = warped_product(
    warp_factors=[H_ns5**R(-1,4), H_ns5**R(-1,2), H_ns5**R(1,2), H_ns5**R(3,4)],
    block_dims=[6, 1, 1, 4],
    block_signatures=['lorentzian', 'euclidean', 'euclidean', 'euclidean'],
    coordinates=coords_ns5,
)

# Magnetic G₄ = *₄^{flat}dH ∧ dz₂
G4_ns5 = FormField(rank=4, dim=12)
z2_idx, trans_start = 7, 8
for i, j, k in combinations(range(4), 3):
    l = [x for x in range(4) if x not in (i,j,k)][0]
    eps = levi_civita_4d(i, j, k, l)
    if eps == 0:
        continue
    dH = sp.diff(H_ns5, hc4[l])
    original = [trans_start+i, trans_start+j, trans_start+k, z2_idx]
    sorted_idx = sorted(original)
    perm_map = [sorted_idx.index(x) for x in original]
    sort_sign = Permutation(perm_map).signature()
    G4_ns5[tuple(sorted_idx)] = sp.Rational(eps) * sort_sign * dH

blocks_5 = [(0,'t'), (1,'x1'), (6,'z1'), (7,'z2'), (8,'y0'), (9,'y1')]
is_torus_5 = lambda idx: idx in (6, 7)

verify_universal("NS5", coords_ns5, hf_ns5, m_ns5, G4_ns5, blocks_5, is_torus_5)


# ===================================================================
# D5-BRANE (codim 4, S-dual of NS5)
# ===================================================================
print("\n" + "#"*60)
print("# D5-BRANE")
print("#"*60)

m_d5 = warped_product(
    warp_factors=[H_ns5**R(-1,4), H_ns5**R(1,2), H_ns5**R(-1,2), H_ns5**R(3,4)],
    block_dims=[6, 1, 1, 4],
    block_signatures=['lorentzian', 'euclidean', 'euclidean', 'euclidean'],
    coordinates=coords_ns5,
)

# G₄ = *₄dH ∧ dz₁
G4_d5 = FormField(rank=4, dim=12)
z1_idx = 6
for i, j, k in combinations(range(4), 3):
    l = [x for x in range(4) if x not in (i,j,k)][0]
    eps = levi_civita_4d(i, j, k, l)
    if eps == 0:
        continue
    dH = sp.diff(H_ns5, hc4[l])
    original = [trans_start+i, trans_start+j, trans_start+k, z1_idx]
    sorted_idx = sorted(original)
    perm_map = [sorted_idx.index(x) for x in original]
    sort_sign = Permutation(perm_map).signature()
    G4_d5[tuple(sorted_idx)] = sp.Rational(eps) * sort_sign * dH

verify_universal("D5", coords_ns5, hf_ns5, m_d5, G4_d5, blocks_5, is_torus_5)


# ===================================================================
# D3-BRANE with F₅ (codim 6)
# ===================================================================
print("\n" + "#"*60)
print("# D3-BRANE (F₅ sector)")
print("#"*60)

wv4 = list(sp.symbols('t x1 x2 x3', real=True))
hc6 = list(sp.symbols('y0:6', real=True))
coords_d3 = wv4 + [z1s, z2s] + hc6

hf_d3 = HarmonicFunction(transverse_coords=hc6)
H_d3 = sp.Function('H')(hf_d3.r_expr)

m_d3 = warped_product(
    warp_factors=[H_d3**R(-1,2), sp.Integer(1), sp.Integer(1), H_d3**R(1,2)],
    block_dims=[4, 1, 1, 6],
    block_signatures=['lorentzian', 'euclidean', 'euclidean', 'euclidean'],
    coordinates=coords_d3,
)

# D3: F₅ = dC₄, self-dual → |F₅|² = 0
# Formula: ℛ = (1/2)FF_{F₅}/4! (no trace term since |F₅|²=0, no torus term since K=0)
C4 = FormField(rank=4, dim=12)
C4[(0, 1, 2, 3)] = 1/H_d3
F5 = exterior_derivative(C4, coords_d3)

print("\n  Computing D3 Ricci and F₅ contraction...")
Ric_d3 = m_d3.ricci_tensor(simplify_func=cancel)
FF_F5 = form_contraction(F5, m_d3)
norm_F5 = form_norm_squared(F5, m_d3)
norm_F5_val = hf_d3.substitute(cancel(norm_F5))
print(f"  |F₅|² = {norm_F5_val}")

# For D3, the formula is ℛ = (1/2)FF/4! since |F₅|²=0 and torus is flat
print(f"\n  D3 formula: ℛ = (1/2)FF_F₅/4!")
blocks_d3 = [(0,'t'), (1,'x1'), (4,'z1'), (5,'z2'), (6,'y0'), (7,'y1')]
all_pass = True
for idx, label in blocks_d3:
    ric = hf_d3.substitute(cancel(Ric_d3[idx, idx]))
    ff = hf_d3.substitute(cancel(FF_F5[idx, idx]))
    g = hf_d3.substitute(cancel(m_d3.matrix[idx, idx]))
    # For self-dual: |F₅|² = 0, but the computed norm is NOT zero
    # because we only have the electric half. Need F₅ = dC₄ + *dC₄.
    # The correct approach: T = (1/2)FF/4! with factor for self-duality.
    # From exp14: ℛ = (1/2)FF/4! with the electric F₅ only, DOUBLED.
    # Actually: for self-dual F₅ = F₅^e + *F₅^e, FF(F₅) = 2·FF(F₅^e)
    # and |F₅|² = 0.
    # So ℛ = (1/2)·2·FF^e/4! = FF^e/4!
    # But let's first check the naive formula.
    T_naive = cancel(R(1,2) * ff / 24)
    res = cancel(ric - T_naive)

    # Also check with factor 2 (self-dual doubling)
    T_double = cancel(ff / 24)
    res2 = cancel(ric - T_double)

    if res == 0:
        print(f"  [{label:>4}]: ✓ (naive 1/2)")
    elif res2 == 0:
        print(f"  [{label:>4}]: ✓ (doubled for self-duality)")
    else:
        print(f"  [{label:>4}]: ✗ naive_res={res}, doubled_res={res2}")
        all_pass = False

print(f"\n  {'ALL PASS ✓' if all_pass else 'SOME FAIL ✗'}")


print("\n" + "="*60)
print("COMPLETE — Universal 12d formula verification")
print("="*60)
