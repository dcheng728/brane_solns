"""Experiment 32d: Fixed verification of 12d formula on NS5, D5, D3.

BUG FIX: The warped_product block_dims must match coordinate ordering.
For coords = wv + [z1, z2] + trans, blocks must be [wv_dim, 1, 1, trans_dim].

Formula: ℛ_{MN} = (1/2)[FF_{MN}/3! - (1/4)|G₄|²g_{MN}] - K · M̃_{MN}
"""
import sys; sys.path.insert(0, 'src')
import sympy as sp
from sympy import Rational as R, cancel
from sugra import (HarmonicFunction, warped_product, Metric,
                   FormField, exterior_derivative,
                   form_contraction, form_norm_squared)
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


def compute_K(M, Minv, coords, metric, hf):
    D = len(coords)
    K = sp.Integer(0)
    free = set()
    for i in range(2):
        for j in range(2):
            free |= M[i,j].free_symbols
    for i in range(D):
        if coords[i] not in free:
            continue
        dM = M.diff(coords[i])
        if dM == sp.zeros(2,2):
            continue
        dMinv = -Minv * dM * Minv
        trace = (dMinv * dM).trace()
        g_inv = hf.substitute(cancel(metric.inv_matrix[i,i]))
        K += g_inv * R(1,4) * trace
    return hf.substitute(cancel(K))


def verify(name, coords, hf, metric, G4, M, Minv, block_indices):
    print(f"\n{'='*60}")
    print(f"  {name}")
    print(f"{'='*60}")

    H_func = sp.Function('H')(hf.r_expr)

    print("  Computing Ricci...")
    Ric = metric.ricci_tensor(simplify_func=cancel)

    if G4 is not None:
        print("  Computing form contraction...")
        FF = form_contraction(G4, metric)
        norm = form_norm_squared(G4, metric)
        norm_val = hf.substitute(cancel(norm))
    else:
        FF = None
        norm_val = sp.Integer(0)

    print("  Computing K...")
    K_val = compute_K(M, Minv, coords, metric, hf)
    print(f"  K = {K_val}")

    print(f"\n  Testing: ℛ = (1/2)[FF/3! - (1/4)|G₄|²g] - K·M̃")
    all_pass = True
    for idx, label, torus_ab in block_indices:
        ric = hf.substitute(cancel(Ric[idx, idx]))
        g = hf.substitute(cancel(metric.matrix[idx, idx]))

        if FF is not None:
            ff = hf.substitute(cancel(FF[idx, idx]))
            T_form = cancel(R(1,2) * (ff/6 - R(1,4) * norm_val * g))
        else:
            T_form = sp.Integer(0)

        if torus_ab is not None:
            m_tilde = hf.substitute(cancel(M[torus_ab, torus_ab]))
        else:
            m_tilde = sp.Integer(0)

        KM = cancel(K_val * m_tilde)
        T_total = cancel(T_form - KM)
        residual = cancel(ric - T_total)

        status = "✓" if residual == 0 else f"✗ residual={residual}"
        print(f"  [{label:>4}]: {status}")

        if residual != 0:
            all_pass = False
            print(f"         ℛ      = {ric}")
            print(f"         T_form = {T_form}")
            print(f"         K·M̃    = {KM}")

    result = 'ALL PASS ✓' if all_pass else 'SOME FAIL ✗'
    print(f"\n  {result}")
    return all_pass


# ===================================================================
# NS5-BRANE (codim 4)
# ===================================================================
print("#"*60)
print("# NS5-BRANE (FIXED block order)")
print("#"*60)

wv6 = list(sp.symbols('t x1 x2 x3 x4 x5', real=True))
z1s, z2s = sp.symbols('z1 z2', real=True)
hc4 = list(sp.symbols('y0:4', real=True))
coords = wv6 + [z1s, z2s] + hc4  # z before transverse

hf = HarmonicFunction(transverse_coords=hc4)
H = sp.Function('H')(hf.r_expr)

# FIXED: blocks match coordinate order [wv(6), z1(1), z2(1), trans(4)]
m_ns5 = warped_product(
    warp_factors=[H**R(-1,4), H**R(-1,2), H**R(1,2), H**R(3,4)],
    block_dims=[6, 1, 1, 4],
    block_signatures=['lorentzian', 'euclidean', 'euclidean', 'euclidean'],
    coordinates=coords,
)

# Verify metric
print("Metric diagonal:")
for i, name in [(0,'t'), (1,'x1'), (6,'z1'), (7,'z2'), (8,'y0')]:
    print(f"  g[{name}] = {m_ns5.matrix[i,i]}")

M_ns5 = sp.Matrix([[H**R(-1,2), 0], [0, H**R(1,2)]])
Minv_ns5 = sp.Matrix([[H**R(1,2), 0], [0, H**R(-1,2)]])

# G₄ = *₄^{flat}dH ∧ dz₂
z2_idx = 7
trans_start = 8

G4_ns5 = FormField(rank=4, dim=12)
for i, j, k in combinations(range(4), 3):
    l = [x for x in range(4) if x not in (i,j,k)][0]
    eps = levi_civita_4d(i, j, k, l)
    if eps == 0:
        continue
    dH = sp.diff(H, hc4[l])
    original = [trans_start+i, trans_start+j, trans_start+k, z2_idx]
    sorted_idx = sorted(original)
    perm_map = [sorted_idx.index(x) for x in original]
    sort_sign = Permutation(perm_map).signature()
    G4_ns5[tuple(sorted_idx)] = sp.Rational(eps) * sort_sign * dH

blocks = [
    (0, 't', None), (1, 'x1', None),
    (6, 'z1', 0), (7, 'z2', 1),
    (8, 'y0', None), (9, 'y1', None),
]

verify("NS5", coords, hf, m_ns5, G4_ns5, M_ns5, Minv_ns5, blocks)


# ===================================================================
# D5-BRANE (S-dual of NS5: z₁↔z₂ swap)
# ===================================================================
print("\n" + "#"*60)
print("# D5-BRANE")
print("#"*60)

# D5: H^{-1/4}ds²_{1,5} + H^{1/2}dz₁² + H^{-1/2}dz₂² + H^{3/4}ds²_4
m_d5 = warped_product(
    warp_factors=[H**R(-1,4), H**R(1,2), H**R(-1,2), H**R(3,4)],
    block_dims=[6, 1, 1, 4],
    block_signatures=['lorentzian', 'euclidean', 'euclidean', 'euclidean'],
    coordinates=coords,
)

M_d5 = sp.Matrix([[H**R(1,2), 0], [0, H**R(-1,2)]])
Minv_d5 = sp.Matrix([[H**R(-1,2), 0], [0, H**R(1,2)]])

# G₄ = *₄dH ∧ dz₁
z1_idx = 6
G4_d5 = FormField(rank=4, dim=12)
for i, j, k in combinations(range(4), 3):
    l = [x for x in range(4) if x not in (i,j,k)][0]
    eps = levi_civita_4d(i, j, k, l)
    if eps == 0:
        continue
    dH = sp.diff(H, hc4[l])
    original = [trans_start+i, trans_start+j, trans_start+k, z1_idx]
    sorted_idx = sorted(original)
    perm_map = [sorted_idx.index(x) for x in original]
    sort_sign = Permutation(perm_map).signature()
    G4_d5[tuple(sorted_idx)] = sp.Rational(eps) * sort_sign * dH

verify("D5", coords, hf, m_d5, G4_d5, M_d5, Minv_d5, blocks)


# ===================================================================
# D3-BRANE with F₅ (codim 6)
# ===================================================================
print("\n" + "#"*60)
print("# D3-BRANE (with F₅)")
print("#"*60)

wv4 = list(sp.symbols('t x1 x2 x3', real=True))
hc6 = list(sp.symbols('y0:6', real=True))
coords_d3 = wv4 + [z1s, z2s] + hc6

hf_d3 = HarmonicFunction(transverse_coords=hc6)
H_d3 = sp.Function('H')(hf_d3.r_expr)

# D3: H^{-1/2}ds²_{1,3} + dz₁² + dz₂² + H^{1/2}ds²_6
# blocks: [wv(4), z1(1), z2(1), trans(6)]
m_d3 = warped_product(
    warp_factors=[H_d3**R(-1,2), sp.Integer(1), sp.Integer(1), H_d3**R(1,2)],
    block_dims=[4, 1, 1, 6],
    block_signatures=['lorentzian', 'euclidean', 'euclidean', 'euclidean'],
    coordinates=coords_d3,
)

# D3: G₄ = 0, K = 0. Full equation is ℛ = (1/2)FF_{F₅}/4!
# Build F₅ = dC₄ with C₄ = H^{-1} vol_{1,3}
C4 = FormField(rank=4, dim=12)
C4[(0, 1, 2, 3)] = 1/H_d3
F5 = exterior_derivative(C4, coords_d3)

print("  Computing D3 Ricci...")
Ric_d3 = m_d3.ricci_tensor(simplify_func=cancel)

print("  Computing F₅ contraction...")
FF_F5 = form_contraction(F5, m_d3)
norm_F5 = form_norm_squared(F5, m_d3)
norm_F5_val = hf_d3.substitute(cancel(norm_F5))
print(f"  |F₅|² = {norm_F5_val}  (should be 0 for self-dual)")

print(f"\n  Testing: ℛ = (1/2)FF_F₅/4!")
blocks_d3 = [
    (0, 't'), (1, 'x1'),
    (4, 'z1'), (5, 'z2'),
    (6, 'y0'), (7, 'y1'),
]
all_pass = True
for idx, label in blocks_d3:
    ric = hf_d3.substitute(cancel(Ric_d3[idx, idx]))
    ff = hf_d3.substitute(cancel(FF_F5[idx, idx]))
    T = cancel(R(1,2) * ff / 24)  # 4! = 24
    residual = cancel(ric - T)
    status = "✓" if residual == 0 else f"✗ residual={residual}"
    print(f"  [{label:>4}]: {status}")
    if residual != 0:
        all_pass = False
        print(f"         ℛ = {ric}")
        print(f"         T = {T}")

print(f"\n  {'ALL PASS ✓' if all_pass else 'SOME FAIL ✗'}")


print("\n" + "="*60)
print("COMPLETE")
print("="*60)
