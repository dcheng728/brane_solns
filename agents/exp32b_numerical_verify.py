"""Experiment 32b: Numerical verification of the 12d formula for ALL branes.

Formula: ℛ_{MN} = (1/2)[FF_{MN}/3! - (1/4)|G₄|²g_{MN}] - K · M̃_{MN}

Use numerical H = 1 + a/r^n with specific coordinate values to verify.
"""
import sys; sys.path.insert(0, 'src')
import sympy as sp
from sympy import Rational as R, cancel, sqrt, diff, simplify
from sugra import (warped_product, Metric, FormField, exterior_derivative,
                   form_contraction, form_norm_squared)
from itertools import combinations
from sympy.combinatorics import Permutation


def levi_civita_4d(i, j, k, l):
    """Levi-Civita symbol in 4d."""
    perm = [i, j, k, l]
    if len(set(perm)) < 4:
        return 0
    sign = 1
    for a in range(4):
        for b in range(a+1, 4):
            if perm[a] > perm[b]:
                sign *= -1
    return sign


def compute_K_num(M, Minv, coords, metric):
    """Compute K = (1/4) g^{PQ} Tr(∂_P M⁻¹ · ∂_Q M), numerical."""
    D = len(coords)
    K = sp.Integer(0)
    free = set()
    for i in range(2):
        for j in range(2):
            free |= M[i,j].free_symbols
    for i in range(D):
        if coords[i] not in free:
            continue
        dM_i = M.diff(coords[i])
        if dM_i == sp.zeros(2, 2):
            continue
        dMinv_i = -Minv * dM_i * Minv
        trace = (dMinv_i * dM_i).trace()
        g_inv_ii = metric.inv_matrix[i, i]
        K += g_inv_ii * R(1, 4) * trace
    return K


def verify_brane_numerical(name, coords, metric, G4, M, Minv,
                           block_info, subs_dict):
    """Verify formula numerically with given substitutions."""
    print(f"\n{'='*60}")
    print(f"  {name}: Numerical verification")
    print(f"{'='*60}")

    # Ricci
    print("  Computing Ricci...")
    Ric = metric.ricci_tensor(simplify_func=cancel)

    # Form
    if G4 is not None:
        print("  Computing form contraction...")
        FF = form_contraction(G4, metric)
        norm = form_norm_squared(G4, metric)
    else:
        FF = None
        norm = sp.Integer(0)

    # K
    print("  Computing K...")
    K_val = compute_K_num(M, Minv, coords, metric)

    # Numerical evaluation
    print("  Substituting numerical values...")
    K_num = float(K_val.subs(subs_dict))
    print(f"  K = {K_num:.8f}")

    all_pass = True
    for idx, label, torus_ab in block_info:
        ric_num = float(Ric[idx, idx].subs(subs_dict))
        g_num = float(metric.matrix[idx, idx].subs(subs_dict))

        if FF is not None:
            ff_num = float(FF[idx, idx].subs(subs_dict))
            norm_num = float(norm.subs(subs_dict))
        else:
            ff_num = 0.0
            norm_num = 0.0

        T_form = 0.5 * (ff_num / 6 - 0.25 * norm_num * g_num)

        if torus_ab is not None:
            m_tilde = float(M[torus_ab, torus_ab].subs(subs_dict))
        else:
            m_tilde = 0.0

        T_total = T_form - K_num * m_tilde
        residual = ric_num - T_total

        rel_err = abs(residual) / max(abs(ric_num), 1e-15) if ric_num != 0 else abs(residual)
        status = "✓" if rel_err < 1e-8 else f"✗ rel_err={rel_err:.2e}"
        print(f"  [{label:>4}]: ℛ={ric_num:+.8f}  T={T_total:+.8f}  {status}")

        if rel_err >= 1e-8:
            all_pass = False
            print(f"         T_form={T_form:+.8f}  K·M̃={K_num*m_tilde:+.8f}")

    print(f"\n  {'ALL PASS ✓' if all_pass else 'SOME FAIL ✗'}")
    return all_pass


# ===================================================================
# Common setup: numerical H
# We use explicit H(y) = 1 + 1/r^n for appropriate n, at a specific point
# ===================================================================
# For string-like (codim 8): H(r) = 1 + a/r^6
# For 5-brane (codim 4): H(r) = 1 + a/r^2
# For D3 (codim 6): H(r) = 1 + a/r^4

z1s, z2s = sp.symbols('z1 z2', real=True)


# ===================================================================
# F1-BRANE (codim 8, sanity check)
# ===================================================================
print("#" * 60)
print("# F1-BRANE (sanity check)")
print("#" * 60)

wv2 = list(sp.symbols('t x1', real=True))
hc8 = list(sp.symbols('y0:8', real=True))
coords_f1 = wv2 + [z1s, z2s] + hc8
r_f1 = sqrt(sum(y**2 for y in hc8))
H_f1 = 1 + 1/r_f1**6

m_f1 = warped_product(
    warp_factors=[H_f1**R(-3,4), H_f1**R(1,2), H_f1**R(-1,2), H_f1**R(1,4)],
    block_dims=[2, 1, 1, 8],
    block_signatures=['lorentzian', 'euclidean', 'euclidean', 'euclidean'],
    coordinates=coords_f1,
)

M_f1 = sp.Matrix([[H_f1**R(1,2), 0], [0, H_f1**R(-1,2)]])
Minv_f1 = sp.Matrix([[H_f1**R(-1,2), 0], [0, H_f1**R(1,2)]])

C3_f1 = FormField(rank=3, dim=12)
C3_f1[(0, 1, 3)] = 1/H_f1
G4_f1 = exterior_derivative(C3_f1, coords_f1)

# Numerical point: y0=2, y1=1, y2=0.5, y3=0.3, rest=0.1
subs_f1 = {hc8[0]: 2.0, hc8[1]: 1.0, hc8[2]: 0.5, hc8[3]: 0.3,
           hc8[4]: 0.1, hc8[5]: 0.1, hc8[6]: 0.1, hc8[7]: 0.1}

blocks_f1 = [
    (0, 't', None), (1, 'x1', None),
    (2, 'z1', 0), (3, 'z2', 1),
    (4, 'y0', None), (5, 'y1', None),
]

verify_brane_numerical("F1", coords_f1, m_f1, G4_f1, M_f1, Minv_f1,
                       blocks_f1, subs_f1)


# ===================================================================
# D1-BRANE (codim 8)
# ===================================================================
print("\n" + "#" * 60)
print("# D1-BRANE")
print("#" * 60)

m_d1 = warped_product(
    warp_factors=[H_f1**R(-3,4), H_f1**R(-1,2), H_f1**R(1,2), H_f1**R(1,4)],
    block_dims=[2, 1, 1, 8],
    block_signatures=['lorentzian', 'euclidean', 'euclidean', 'euclidean'],
    coordinates=coords_f1,
)

M_d1 = sp.Matrix([[H_f1**R(-1,2), 0], [0, H_f1**R(1,2)]])
Minv_d1 = sp.Matrix([[H_f1**R(1,2), 0], [0, H_f1**R(-1,2)]])

C3_d1 = FormField(rank=3, dim=12)
C3_d1[(0, 1, 2)] = 1/H_f1  # C_{t,x1,z1}
G4_d1 = exterior_derivative(C3_d1, coords_f1)

verify_brane_numerical("D1", coords_f1, m_d1, G4_d1, M_d1, Minv_d1,
                       blocks_f1, subs_f1)


# ===================================================================
# NS5-BRANE (codim 4)
# ===================================================================
print("\n" + "#" * 60)
print("# NS5-BRANE")
print("#" * 60)

wv6 = list(sp.symbols('t x1 x2 x3 x4 x5', real=True))
hc4 = list(sp.symbols('y0:4', real=True))
coords_ns5 = wv6 + [z1s, z2s] + hc4
r_ns5 = sqrt(sum(y**2 for y in hc4))
H_ns5 = 1 + 1/r_ns5**2

m_ns5 = warped_product(
    warp_factors=[H_ns5**R(-1,4), H_ns5**R(3,4), H_ns5**R(-1,2), H_ns5**R(1,2)],
    block_dims=[6, 4, 1, 1],
    block_signatures=['lorentzian', 'euclidean', 'euclidean', 'euclidean'],
    coordinates=coords_ns5,
)

M_ns5 = sp.Matrix([[H_ns5**R(-1,2), 0], [0, H_ns5**R(1,2)]])
Minv_ns5 = sp.Matrix([[H_ns5**R(1,2), 0], [0, H_ns5**R(-1,2)]])

# G₄ = *₄dH ∧ dz₂ (magnetic)
z2_idx = 7  # z2 position in coords_ns5
trans_start = 8  # y0 starts here

G4_ns5 = FormField(rank=4, dim=12)
for i, j, k in combinations(range(4), 3):
    l = [x for x in range(4) if x not in (i,j,k)][0]
    eps = levi_civita_4d(i, j, k, l)
    if eps == 0:
        continue
    dH = diff(H_ns5, hc4[l])
    original = [trans_start+i, trans_start+j, trans_start+k, z2_idx]
    sorted_idx = sorted(original)
    perm_map = [sorted_idx.index(x) for x in original]
    sort_sign = Permutation(perm_map).signature()
    G4_ns5[tuple(sorted_idx)] = sp.Rational(eps) * sort_sign * dH

# Numerical point
subs_ns5 = {hc4[0]: 2.0, hc4[1]: 1.5, hc4[2]: 0.8, hc4[3]: 0.6}

blocks_ns5 = [
    (0, 't', None), (1, 'x1', None),
    (6, 'z1', 0), (7, 'z2', 1),
    (8, 'y0', None), (9, 'y1', None),
]

verify_brane_numerical("NS5", coords_ns5, m_ns5, G4_ns5, M_ns5, Minv_ns5,
                       blocks_ns5, subs_ns5)


# ===================================================================
# D5-BRANE (codim 4)
# ===================================================================
print("\n" + "#" * 60)
print("# D5-BRANE")
print("#" * 60)

m_d5 = warped_product(
    warp_factors=[H_ns5**R(-1,4), H_ns5**R(3,4), H_ns5**R(1,2), H_ns5**R(-1,2)],
    block_dims=[6, 4, 1, 1],
    block_signatures=['lorentzian', 'euclidean', 'euclidean', 'euclidean'],
    coordinates=coords_ns5,
)

M_d5 = sp.Matrix([[H_ns5**R(1,2), 0], [0, H_ns5**R(-1,2)]])
Minv_d5 = sp.Matrix([[H_ns5**R(-1,2), 0], [0, H_ns5**R(1,2)]])

# G₄ = *₄dH ∧ dz₁
z1_idx = 6
G4_d5 = FormField(rank=4, dim=12)
for i, j, k in combinations(range(4), 3):
    l = [x for x in range(4) if x not in (i,j,k)][0]
    eps = levi_civita_4d(i, j, k, l)
    if eps == 0:
        continue
    dH = diff(H_ns5, hc4[l])
    original = [trans_start+i, trans_start+j, trans_start+k, z1_idx]
    sorted_idx = sorted(original)
    perm_map = [sorted_idx.index(x) for x in original]
    sort_sign = Permutation(perm_map).signature()
    G4_d5[tuple(sorted_idx)] = sp.Rational(eps) * sort_sign * dH

verify_brane_numerical("D5", coords_ns5, m_d5, G4_d5, M_d5, Minv_d5,
                       blocks_ns5, subs_ns5)


# ===================================================================
# D3-BRANE (codim 6, with F₅)
# ===================================================================
print("\n" + "#" * 60)
print("# D3-BRANE")
print("#" * 60)

wv4 = list(sp.symbols('t x1 x2 x3', real=True))
hc6 = list(sp.symbols('y0:6', real=True))
coords_d3 = wv4 + [z1s, z2s] + hc6
r_d3 = sqrt(sum(y**2 for y in hc6))
H_d3 = 1 + 1/r_d3**4

m_d3 = warped_product(
    warp_factors=[H_d3**R(-1,2), sp.Integer(1), sp.Integer(1), H_d3**R(1,2)],
    block_dims=[4, 1, 1, 6],
    block_signatures=['lorentzian', 'euclidean', 'euclidean', 'euclidean'],
    coordinates=coords_d3,
)

M_d3 = sp.eye(2)
Minv_d3 = sp.eye(2)

# F₅ for D3: self-dual 5-form
# F₅ = dC₄ with C₄ = H^{-1} dt∧dx1∧dx2∧dx3
# The formula for D3 is: ℛ = (1/2)FF_{F₅}/4! (since |F₅|²=0, K=0)
# But our formula uses G₄ (which is 0 for D3) + F₅ separately.
# Let me just check ℛ = (1/2)FF_{F₅}/4! directly.

print("  D3 uses F₅ (self-dual). Testing ℛ = (1/2)FF_{F₅}/4!")
print("  Computing Ricci...")
Ric_d3 = m_d3.ricci_tensor(simplify_func=cancel)

# Build F₅: C₄ = H^{-1} vol_{1,3}
# dC₄ has components in (y_k, t, x1, x2, x3) directions
C4 = FormField(rank=4, dim=12)
C4[(0, 1, 2, 3)] = 1/H_d3  # t,x1,x2,x3 = indices 0,1,2,3
F5 = exterior_derivative(C4, coords_d3)

print("  Computing F₅ contraction...")
FF_F5 = form_contraction(F5, m_d3)

# D3 numerical point
subs_d3 = {hc6[0]: 2.0, hc6[1]: 1.0, hc6[2]: 0.7, hc6[3]: 0.5,
           hc6[4]: 0.3, hc6[5]: 0.2}

blocks_d3 = [
    (0, 't'), (1, 'x1'),
    (4, 'z1'), (5, 'z2'),
    (6, 'y0'), (7, 'y1'),
]

print(f"\n  Formula: ℛ = (1/2)FF_{'{F₅}'}/4!")
all_pass = True
for idx, label in blocks_d3:
    ric_num = float(Ric_d3[idx, idx].subs(subs_d3))
    ff_num = float(FF_F5[idx, idx].subs(subs_d3))
    T_num = 0.5 * ff_num / 24  # 4! = 24
    residual = ric_num - T_num
    rel_err = abs(residual) / max(abs(ric_num), 1e-15) if ric_num != 0 else abs(residual)
    status = "✓" if rel_err < 1e-8 else f"✗ rel_err={rel_err:.2e}"
    print(f"  [{label:>4}]: ℛ={ric_num:+.8f}  T={T_num:+.8f}  {status}")
    if rel_err >= 1e-8:
        all_pass = False

print(f"\n  {'ALL PASS ✓' if all_pass else 'SOME FAIL ✗'}")


print("\n" + "=" * 60)
print("COMPLETE")
print("=" * 60)
