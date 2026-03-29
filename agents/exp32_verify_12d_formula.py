"""Experiment 32: Verify the 12d formula on ALL fundamental branes.

The candidate formula from exp31 (verified on F1):
    ℛ_{MN} = (1/2)[FF_{MN}/3! - (1/4)|G₄|²g_{MN}] - K · M̃_{MN}

where:
    K = (1/4) g^{PQ} Tr(∂_P M⁻¹ · ∂_Q M)     [scalar kinetic scalar]
    M̃_{MN} = M_{ab} embedded in 12d             [torus projection]

Test on: NS5, D5, D3, D1, and (1,1)-string.
"""
import sys; sys.path.insert(0, 'src')
import sympy as sp
from sympy import Rational as R, cancel, simplify
from sugra import (HarmonicFunction, warped_product, Metric,
                   FormField, exterior_derivative,
                   form_contraction, form_norm_squared)


def compute_K(M, Minv, coords, metric, hf):
    """Compute K = (1/4) g^{PQ} Tr(∂_P M⁻¹ · ∂_Q M)."""
    D = len(coords)
    K = sp.Integer(0)
    for i in range(D):
        dM_i = M.diff(coords[i]) if coords[i] in M.free_symbols else None
        if dM_i is None or dM_i == sp.zeros(2, 2):
            continue
        dMinv_i = -Minv * dM_i * Minv
        dM_i_again = dM_i  # ∂_i M (same index, diagonal of S)
        trace = (dMinv_i * dM_i_again).trace()
        g_inv_ii = hf.substitute(cancel(metric.inv_matrix[i, i]))
        K += g_inv_ii * R(1, 4) * trace
    return hf.substitute(cancel(K))


def verify_formula(name, coords, hf, metric, G4, M, Minv, block_info, transverse_coords):
    """Verify ℛ = (1/2)[FF/3! - (1/4)|G₄|²g] - K·M̃ for all blocks.

    block_info: list of (idx, label, is_torus_idx)
    """
    print(f"\n{'='*60}")
    print(f"  {name}: Verifying 12d formula")
    print(f"{'='*60}")

    H = sp.Function('H')(hf.r_expr)

    # Ricci
    print("  Computing Ricci...")
    Ric = metric.ricci_tensor(simplify_func=cancel)

    # Form building blocks
    if G4 is not None:
        print("  Computing form contraction and norm...")
        FF = form_contraction(G4, metric)
        norm = form_norm_squared(G4, metric)
    else:
        FF = None
        norm = sp.Integer(0)

    # Scalar kinetic scalar K
    print("  Computing K...")
    K_val = compute_K(M, Minv, coords, metric, hf)
    print(f"  K = {K_val}")

    # Test formula per block
    print(f"\n  Formula: ℛ = (1/2)[FF/3! - (1/4)|G₄|²g] - K·M̃")
    print(f"  {'Block':<8} {'ℛ':>12} {'T_form':>12} {'K·M̃':>12} {'residual':>12}")
    print(f"  {'-'*56}")

    all_pass = True
    for idx, label, torus_ab in block_info:
        ric = hf.substitute(cancel(Ric[idx, idx]))
        g = hf.substitute(cancel(metric.matrix[idx, idx]))

        if FF is not None:
            ff = hf.substitute(cancel(FF[idx, idx]))
            n = hf.substitute(cancel(norm))
            T_form = R(1, 2) * (ff / 6 - R(1, 4) * n * g)
        else:
            T_form = sp.Integer(0)
        T_form = cancel(T_form)

        # Torus projection M̃
        if torus_ab is not None:
            m_tilde = hf.substitute(cancel(M[torus_ab, torus_ab]))
        else:
            m_tilde = sp.Integer(0)

        KM = cancel(K_val * m_tilde)
        T_total = cancel(T_form - KM)
        residual = cancel(ric - T_total)

        # For transverse coords, residual may have y_k terms that sum to r²
        # Substitute sum(y_k²) = r² explicitly
        if residual != 0 and transverse_coords:
            r = hf.r_expr
            sum_y2 = sum(y**2 for y in transverse_coords)
            # Try substituting r² = sum_y² (which is an identity)
            residual_sub = residual.subs(sum_y2, r**2) if hasattr(r, 'args') else residual
            residual_sub = cancel(residual_sub)
            # Also try: factor out and simplify
            if residual_sub != 0:
                # Check if residual is proportional to (r² - Σy²) which = 0
                test = cancel(residual * r**2)
                # replace r² with sum_y²
                test2 = test.subs(r**2, sum_y2)
                test2 = sp.expand(test2)
                if test2 == 0:
                    residual = sp.Integer(0)
                else:
                    residual = residual_sub

        status = "✓" if residual == 0 else f"✗ {residual}"
        print(f"  {label:<8} {status}")

        if residual != 0:
            all_pass = False
            # Print detailed breakdown
            print(f"         ℛ     = {ric}")
            print(f"         T_form= {T_form}")
            print(f"         K·M̃   = {KM}")

    print(f"\n  {'ALL PASS ✓' if all_pass else 'SOME FAIL ✗'}")
    return all_pass


# ===================================================================
# D3-BRANE
# ===================================================================
print("\n" + "#"*60)
print("# D3-BRANE")
print("#"*60)

wv4 = list(sp.symbols('t x1 x2 x3', real=True))
z1s, z2s = sp.symbols('z1 z2', real=True)
hc6 = list(sp.symbols('y0:6', real=True))
coords_d3 = wv4 + [z1s, z2s] + hc6
hf_d3 = HarmonicFunction(transverse_coords=hc6)
H_d3 = sp.Function('H')(hf_d3.r_expr)

# D3: ds² = H^{-1/2}ds²_{1,3} + dz₁² + dz₂² + H^{1/2}ds²_6
m_d3 = warped_product(
    warp_factors=[H_d3**R(-1,2), sp.Integer(1), sp.Integer(1), H_d3**R(1,2)],
    block_dims=[4, 1, 1, 6],
    block_signatures=['lorentzian', 'euclidean', 'euclidean', 'euclidean'],
    coordinates=coords_d3,
)

M_d3 = sp.eye(2)  # flat torus
Minv_d3 = sp.eye(2)

# D3 has F₅, not G₄. The formula uses G₄=0 for D3.
# F₅ self-dual → separate treatment. First check G₄=0 piece.
blocks_d3 = [
    (0, 't', None), (1, 'x1', None),
    (4, 'z1', 0), (5, 'z2', 1),
    (6, 'y0', None), (7, 'y1', None),
]

verify_formula("D3 (G₄=0 sector)", coords_d3, hf_d3, m_d3, None,
               M_d3, Minv_d3, blocks_d3, hc6)

# Now check D3 with F₅ stress-energy
print("\n  --- D3: Adding F₅ contribution ---")
print("  D3 uses self-dual F₅. Since |F₅|²=0, T^{F₅} = (1/2)FF_{F₅}/4!")
print("  (verified in exp14)")

# ===================================================================
# NS5-BRANE
# ===================================================================
print("\n" + "#"*60)
print("# NS5-BRANE")
print("#"*60)

wv6 = list(sp.symbols('t x1 x2 x3 x4 x5', real=True))
hc4 = list(sp.symbols('y0:4', real=True))
coords_ns5 = wv6 + [z1s, z2s] + hc4
hf_ns5 = HarmonicFunction(transverse_coords=hc4)
H_ns5 = sp.Function('H')(hf_ns5.r_expr)

# NS5: ds² = H^{-1/4}ds²_{1,5} + H^{3/4}ds²_4 + H^{-1/2}dz₁² + H^{1/2}dz₂²
m_ns5 = warped_product(
    warp_factors=[H_ns5**R(-1,4), H_ns5**R(3,4), H_ns5**R(-1,2), H_ns5**R(1,2)],
    block_dims=[6, 4, 1, 1],
    block_signatures=['lorentzian', 'euclidean', 'euclidean', 'euclidean'],
    coordinates=coords_ns5,
)

# NS5 torus: M = diag(H^{-1/2}, H^{1/2})  (Φ=+(1/2)lnH)
M_ns5 = sp.Matrix([[H_ns5**R(-1,2), 0], [0, H_ns5**R(1,2)]])
Minv_ns5 = sp.Matrix([[H_ns5**R(1,2), 0], [0, H_ns5**R(-1,2)]])

# G₄ for NS5 (magnetic): G₄ = *₄dH ∧ dz₂
# In practice, *₄dH has components (dy1∧dy2∧dy3, -dy0∧dy2∧dy3, dy0∧dy1∧dy3, -dy0∧dy1∧dy2)
# C₃ with legs in transverse + z₂: C_{yi,yj,z₂} = ...
# Actually, for the magnetic 3-form: H₃ = *₄dH = ε_{ijkl}∂_lH dy^i∧dy^j∧dy^k / 3!
# where *₄ is the flat Hodge dual in transverse R⁴.
#
# For H(r) in 4d: dH = H' y_k/r dy^k
# *₄dH has components proportional to ε_{ijk l} y_l/r × H'
#
# G₄ = H₃^{mag} ∧ dz₂:
# Components: G₄_{y_i, y_j, y_k, z₂} = ε_{ijkl} H' y_l / r

# Build G₄ for NS5 manually
from itertools import permutations

y0n, y1n, y2n, y3n = hc4
G4_ns5 = FormField(rank=4, dim=12)
trans_start = 6  # y0 is at index 6 in coords_ns5

# ε_{ijkl} (y_l / r) H'  for the 3-form part, then ∧ dz₂
# z₂ is at index 7 in coords_ns5
z2_idx = 7  # index of z2 in coords_ns5

# Actually z1 is at index 6, z2 at index 7, y0 at index 8, etc.
# Wait, coords_ns5 = wv6(0-5) + [z1(6), z2(7)] + hc4(8-11)
z1_idx_ns5 = 6
z2_idx_ns5 = 7
trans_start_ns5 = 8

# ε_{0123} = +1 (4d Levi-Civita in transverse coords y0,y1,y2,y3)
# G₄_{y_i,y_j,y_k,z₂} = ε_{ijkl} * (H'/r) * y_l
# For i<j<k, l is determined. The sign of ε_{ijkl} with l fixed:
Hp = sp.Function("H'")(hf_ns5.r_expr)  # stand-in

# Actually, let me just set the 4 independent components:
# (y0,y1,y2,z2): ε_{0123→3} → l=3, ε_{0123}=+1 → +H'y3/r
# Wait this is getting complicated. Let me use the exterior derivative approach.
# The magnetic potential is C₃ with legs in transverse space.
# For NS5: H₃ = *₄dH, then G₄ = H₃ ∧ dz₂.
# Let's construct it from a potential. In 4d: *dH = -d(*H) + ...
# Actually for a harmonic function in 4d, *₄dH is closed (since d*dH = ∇²H = 0).
# But constructing the potential is nontrivial.
#
# Alternative: just manually set the 4-form components.
# G₄_{y_i,y_j,y_k,z₂} = ε_{ijkl} H' y_l / r  (flat ε, normalized)

# In the 12d coord ordering, transverse y_a is at index trans_start_ns5 + a
def levi_civita_4d(i, j, k, l):
    """Levi-Civita symbol in 4d."""
    perm = [i, j, k, l]
    if len(set(perm)) < 4:
        return 0
    # Count transpositions
    n = 0
    p = list(perm)
    for ii in range(4):
        while p[ii] != ii:
            target = p[ii]
            p[ii], p[target] = p[target], p[ii]
            n += 1
    return (-1)**n

G4_ns5 = FormField(rank=4, dim=12)
r_ns5 = hf_ns5.r_expr
H_ns5_func = H_ns5

for i in range(4):
    for j in range(i+1, 4):
        for k in range(j+1, 4):
            # l is the remaining index
            l = [x for x in range(4) if x not in (i,j,k)][0]
            eps = levi_civita_4d(i, j, k, l)
            if eps == 0:
                continue
            # G₄_{y_i, y_j, y_k, z₂}
            idx_12d = tuple(sorted([trans_start_ns5+i, trans_start_ns5+j,
                                    trans_start_ns5+k, z2_idx_ns5]))
            # Need to account for the ordering — z2_idx_ns5=7 < trans_start_ns5=8
            # So z₂ comes before all y's in index ordering
            # idx_12d = (7, 8+i, 8+j, 8+k) already sorted since 7 < 8
            val = eps * hc4[l] / r_ns5  # H' * y_l / r, but H' comes from dH
            G4_ns5[idx_12d] = val

# The G₄ components should be: dH ∧ ... so they involve H'(r).
# Actually G₄ = H₃^{mag} ∧ dz₂ where H₃ = *₄dH.
# dH = H'(r) * Σ (y_k/r) dy^k
# *₄dH in flat R⁴: (*₄dH)_{ijk} = ε_{ijkl} (H'/r) y_l
# Then G₄ = H₃ ∧ dz₂: G₄_{ijk,z₂} = (H'/r) ε_{ijkl} y_l

# But we need to be careful: the FormField stores the component of the
# exterior derivative, not a potential. For magnetic forms, we can't easily
# use exterior_derivative. Let me just compute form_contraction manually.

# Actually, let me rethink. The form_contraction and form_norm_squared
# take a FormField. Let me set the components directly with H' included.

# Build G₄ using sp.diff(H, y_l) which applies chain rule for H(r(y))
from itertools import combinations
from sympy.combinatorics import Permutation

G4_ns5 = FormField(rank=4, dim=12)
for i, j, k in combinations(range(4), 3):
    l = [x for x in range(4) if x not in (i,j,k)][0]
    eps = levi_civita_4d(i, j, k, l)
    if eps == 0:
        continue
    dH = sp.diff(H_ns5, hc4[l])
    # 12d indices: y_a at trans_start_ns5+a, z₂ at z2_idx_ns5
    original = [trans_start_ns5+i, trans_start_ns5+j, trans_start_ns5+k, z2_idx_ns5]
    sorted_idx = sorted(original)
    perm_map = [sorted_idx.index(x) for x in original]
    sort_sign = Permutation(perm_map).signature()
    G4_ns5[tuple(sorted_idx)] = sp.Rational(eps) * sort_sign * dH

blocks_ns5 = [
    (0, 't', None), (1, 'x1', None),
    (6, 'z1', 0), (7, 'z2', 1),
    (8, 'y0', None), (9, 'y1', None),
]

verify_formula("NS5", coords_ns5, hf_ns5, m_ns5, G4_ns5,
               M_ns5, Minv_ns5, blocks_ns5, hc4)


# ===================================================================
# D5-BRANE
# ===================================================================
print("\n" + "#"*60)
print("# D5-BRANE")
print("#"*60)

# D5: ds² = H^{-1/4}ds²_{1,5} + H^{3/4}ds²_4 + H^{1/2}dz₁² + H^{-1/2}dz₂²
# (S-dual of NS5: z₁↔z₂ swapped)
m_d5 = warped_product(
    warp_factors=[H_ns5**R(-1,4), H_ns5**R(3,4), H_ns5**R(1,2), H_ns5**R(-1,2)],
    block_dims=[6, 4, 1, 1],
    block_signatures=['lorentzian', 'euclidean', 'euclidean', 'euclidean'],
    coordinates=coords_ns5,
)

# D5 torus: M = diag(H^{1/2}, H^{-1/2})  (Φ=-(1/2)lnH, S-dual of NS5)
M_d5 = sp.Matrix([[H_ns5**R(1,2), 0], [0, H_ns5**R(-1,2)]])
Minv_d5 = sp.Matrix([[H_ns5**R(-1,2), 0], [0, H_ns5**R(1,2)]])

# G₄ for D5: G₄ = *₄dH ∧ dz₁ (c=(1,0))
G4_d5 = FormField(rank=4, dim=12)
for i, j, k in combinations(range(4), 3):
    l = [x for x in range(4) if x not in (i,j,k)][0]
    eps = levi_civita_4d(i, j, k, l)
    if eps == 0:
        continue
    dH = sp.diff(H_ns5, hc4[l])
    original = [trans_start_ns5+i, trans_start_ns5+j, trans_start_ns5+k, z1_idx_ns5]
    sorted_idx = sorted(original)
    perm_map = [sorted_idx.index(x) for x in original]
    sort_sign = Permutation(perm_map).signature()
    G4_d5[tuple(sorted_idx)] = sp.Rational(eps) * sort_sign * dH

blocks_d5 = [
    (0, 't', None), (1, 'x1', None),
    (6, 'z1', 0), (7, 'z2', 1),
    (8, 'y0', None), (9, 'y1', None),
]

verify_formula("D5", coords_ns5, hf_ns5, m_d5, G4_d5,
               M_d5, Minv_d5, blocks_d5, hc4)


# ===================================================================
# D1-BRANE (quick check — S-dual of F1)
# ===================================================================
print("\n" + "#"*60)
print("# D1-BRANE")
print("#"*60)

wv2 = list(sp.symbols('t x1', real=True))
hc8 = list(sp.symbols('y0:8', real=True))
coords_d1 = wv2 + [z1s, z2s] + hc8
hf_d1 = HarmonicFunction(transverse_coords=hc8)
H_d1 = sp.Function('H')(hf_d1.r_expr)

# D1: ds² = H^{-3/4}ds²_{1,1} + H^{-1/2}dz₁² + H^{1/2}dz₂² + H^{1/4}ds²_8
# (z₁↔z₂ vs F1)
m_d1 = warped_product(
    warp_factors=[H_d1**R(-3,4), H_d1**R(-1,2), H_d1**R(1,2), H_d1**R(1,4)],
    block_dims=[2, 1, 1, 8],
    block_signatures=['lorentzian', 'euclidean', 'euclidean', 'euclidean'],
    coordinates=coords_d1,
)

# D1 torus: M = diag(H^{-1/2}, H^{1/2})  (Φ=+(1/2)lnH)
M_d1 = sp.Matrix([[H_d1**R(-1,2), 0], [0, H_d1**R(1,2)]])
Minv_d1 = sp.Matrix([[H_d1**R(1,2), 0], [0, H_d1**R(-1,2)]])

# G₄ for D1: G₄ = F₃ ∧ dz₁, C_{t,x1,z1} = H^{-1}
C3_d1 = FormField(rank=3, dim=12)
C3_d1[(0, 1, 2)] = 1/H_d1  # t=0, x1=1, z1=2
G4_d1 = exterior_derivative(C3_d1, coords_d1)

blocks_d1 = [
    (0, 't', None), (1, 'x1', None),
    (2, 'z1', 0), (3, 'z2', 1),
    (4, 'y0', None), (5, 'y1', None),
]

verify_formula("D1", coords_d1, hf_d1, m_d1, G4_d1,
               M_d1, Minv_d1, blocks_d1, hc8)


print("\n" + "="*60)
print("DONE — All brane tests complete")
print("="*60)
