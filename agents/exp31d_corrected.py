"""Experiment 31d: Verify 12d formula — corrected NS5 and D3.

FORMULA:
  ℛ_{MN} = (1/2)[(G₄²)_{MN}/3! - (1/4)|G₄|²g_{MN}] - K · M̃_{MN}
          + (for D3: F₅ contribution)
"""
import sys; sys.path.insert(0, 'src')
import sympy as sp
from sympy import Rational as R, cancel
from sympy.combinatorics import Permutation
from sugra import (HarmonicFunction, warped_product, Metric,
                   FormField, exterior_derivative,
                   form_contraction, form_norm_squared)
from itertools import combinations


def compute_K_scalar(M, Minv, coords, metric, hf):
    """K = (1/4) g^{PQ} Tr(∂_P M⁻¹ · ∂_Q M)."""
    K = sp.Integer(0)
    for mu in range(len(coords)):
        x = coords[mu]
        dM = M.diff(x) if x in M.free_symbols else None
        if dM is None or dM == sp.zeros(2, 2):
            continue
        dMinv = -Minv * dM * Minv
        trace_val = (dMinv * dM).trace()
        g_inv = metric.inv_matrix[mu, mu]
        K += R(1, 4) * g_inv * trace_val
    return hf.substitute(cancel(K))


def levi_civita_4(i, j, k, l):
    perm = [i, j, k, l]
    if len(set(perm)) < 4:
        return 0
    sign = 1
    for a in range(4):
        for b in range(a+1, 4):
            if perm[a] > perm[b]:
                sign *= -1
    return sign


def verify_formula(name, Ric, metric, FF, norm_val, K, M, hf,
                   torus_indices, labels):
    """Check ℛ = (1/2)[FF/3! - (1/4)|G₄|²g] - K·M̃ for all diagonal components."""
    print(f"\n{name}: ℛ = (1/2)[FF/3! - (1/4)|G₄|²g] - K·M̃")
    print("-"*60)
    all_pass = True
    for i12, label in labels:
        ric_val = hf.substitute(cancel(Ric[i12, i12]))
        ff_val = hf.substitute(cancel(FF[i12, i12]))
        g_val = hf.substitute(cancel(metric.matrix[i12, i12]))

        T_form = cancel(R(1,2) * (ff_val/6 - R(1,4)*norm_val*g_val))

        if i12 in torus_indices:
            a_idx = torus_indices.index(i12)
            m_tilde = hf.substitute(cancel(M[a_idx, a_idx]))
            T_torus = cancel(-K * m_tilde)
        else:
            T_torus = sp.Integer(0)

        T_total = cancel(T_form + T_torus)
        diff = cancel(ric_val - T_total)
        status = "✓" if diff == 0 else f"✗ residual={diff}"
        if diff != 0:
            all_pass = False
        print(f"  [{label}]: diff={diff}  {status}")

    # Off-diagonal torus
    i_z1, i_z2 = torus_indices
    ric_cross = hf.substitute(cancel(Ric[i_z1, i_z2]))
    ff_cross = hf.substitute(cancel(FF[i_z1, i_z2]))
    g_cross = hf.substitute(cancel(metric.matrix[i_z1, i_z2]))
    m_cross = hf.substitute(cancel(M[0, 1]))
    T_cross = cancel(R(1,2)*(ff_cross/6 - R(1,4)*norm_val*g_cross) - K*m_cross)
    diff_cross = cancel(ric_cross - T_cross)
    status = "✓" if diff_cross == 0 else f"✗ residual={diff_cross}"
    if diff_cross != 0:
        all_pass = False
    print(f"  [z1z2]: diff={diff_cross}  {status}")
    print(f"\n  {'ALL PASS ✓' if all_pass else 'SOME FAIL ✗'}")
    return all_pass


# ===================================================================
# NS5-BRANE (following exp19 pattern exactly)
# ===================================================================
print("="*70)
print("NS5-BRANE (magnetic)")
print("="*70)

wv_ns5 = list(sp.symbols('t x1 x2 x3 x4 x5', real=True))
z1s, z2s = sp.symbols('z1 z2', real=True)
trans_ns5 = list(sp.symbols('y0:4', real=True))
coords_ns5 = wv_ns5 + [z1s, z2s] + trans_ns5  # 12d

hf_ns5 = HarmonicFunction(transverse_coords=trans_ns5)
H = sp.Function('H')(hf_ns5.r_expr)

m_ns5 = warped_product(
    [H**R(-1,4), H**R(-1,2), H**R(1,2), H**R(3,4)],
    [6, 1, 1, 4], ['lorentzian','euclidean','euclidean','euclidean'],
    coords_ns5)

M_ns5 = sp.Matrix([[H**R(-1,2), 0], [0, H**R(1,2)]])
Minv_ns5 = sp.Matrix([[H**R(1,2), 0], [0, H**R(-1,2)]])

print("Computing Ricci...")
Ric_ns5 = m_ns5.ricci_tensor(simplify_func=cancel)
print("Done.")

K_ns5 = compute_K_scalar(M_ns5, Minv_ns5, coords_ns5, m_ns5, hf_ns5)
print(f"K(NS5) = {K_ns5}")

# G₄ = H₃^{flat} ∧ dz₂ with H₃_{ijk} = ε_{ijkl} ∂_{yl}H  (flat Hodge dual)
# Following exp19 exactly:
G4_ns5 = FormField(rank=4, dim=12)
trans_idx = list(range(4))  # internal 0,1,2,3
z2_pos = 7  # z2 is index 7 in 12d coords
y_base = 8  # y0 starts at index 8

for i, j, k in combinations(trans_idx, 3):
    l = [x for x in trans_idx if x not in (i, j, k)][0]
    eps = levi_civita_4(i, j, k, l)
    if eps == 0:
        continue
    y_l = trans_ns5[l]
    dH = sp.diff(H, y_l)  # = H'(r) · yl/r
    original = [y_base+i, y_base+j, y_base+k, z2_pos]
    sorted_idx = tuple(sorted(original))
    perm_map = [sorted_idx.index(x) for x in original]
    sort_sign = Permutation(perm_map).signature()
    G4_ns5[sorted_idx] = sp.Rational(eps) * sort_sign * dH

print("\nComputing NS5 form contraction and norm...")
FF_ns5 = form_contraction(G4_ns5, m_ns5)
norm_ns5 = form_norm_squared(G4_ns5, m_ns5)
norm_ns5_val = hf_ns5.substitute(cancel(norm_ns5))
print(f"|G₄|²(NS5) = {norm_ns5_val}")

labels_ns5 = [(0,'t'), (1,'x1'), (6,'z1'), (7,'z2'), (8,'y0'), (9,'y1')]
verify_formula("NS5", Ric_ns5, m_ns5, FF_ns5, norm_ns5_val,
               K_ns5, M_ns5, hf_ns5, [6,7], labels_ns5)


# ===================================================================
# D3-BRANE (self-dual F₅)
# ===================================================================
print("\n\n" + "="*70)
print("D3-BRANE (self-dual F₅)")
print("="*70)

wv_d3 = list(sp.symbols('t x1 x2 x3', real=True))
trans_d3 = list(sp.symbols('y0:6', real=True))
coords_d3 = wv_d3 + [z1s, z2s] + trans_d3

hf_d3 = HarmonicFunction(transverse_coords=trans_d3)
H_d3 = sp.Function('H')(hf_d3.r_expr)

m_d3 = warped_product(
    [H_d3**R(-1,2), sp.Integer(1), sp.Integer(1), H_d3**R(1,2)],
    [4, 1, 1, 6], ['lorentzian','euclidean','euclidean','euclidean'],
    coords_d3)

print("Computing D3 Ricci...")
Ric_d3 = m_d3.ricci_tensor(simplify_func=cancel)
print("Done.")

# For D3: flat torus, K=0, G₄=0.
# Only F₅. C₄ with C_{t,x1,x2,x3} = H^{-1}
C4 = FormField(rank=4, dim=12)
C4[(0, 1, 2, 3)] = 1/H_d3
F5_d3 = exterior_derivative(C4, coords_d3)

FF_F5 = form_contraction(F5_d3, m_d3)
norm_F5 = hf_d3.substitute(cancel(form_norm_squared(F5_d3, m_d3)))
print(f"|F₅|² (electric only) = {norm_F5}")

# For self-dual F₅, the trace |F₅|² = 0. The stress energy is:
# T_{MN} = (1/2) FF_{MN} / (p-1)! with p=5, so T = (1/2) FF/4!
# But for the self-dual case, the electric contribution is doubled:
# T^{SD}_{MN} = 2 × (1/2) FF^{elec}/4! = FF^{elec}/4!
# OR equivalently: T = (1/2) FF^{elec}/4! with NO trace subtraction
# (since the trace cancels between electric and magnetic).
#
# From exp14: ℛ = (1/2) FF^{F₅}/4! was verified.
# But exp14 used the FULL self-dual F₅ (both electric and magnetic halves).
# With only the electric half, we need to check the coefficient.

labels_d3 = [(0,'t'), (1,'x1'), (4,'z1'), (5,'z2'), (6,'y0'), (7,'y1')]

print(f"\nD3: Finding correct F₅ coefficient")
print("-"*60)

for i12, label in labels_d3:
    ric_val = hf_d3.substitute(cancel(Ric_d3[i12, i12]))
    ff_val = hf_d3.substitute(cancel(FF_F5[i12, i12]))
    g_val = hf_d3.substitute(cancel(m_d3.matrix[i12, i12]))

    # Try: ℛ = α·FF/4! + β·|F₅|²·g
    if ff_val != 0:
        # Try with no trace (self-dual):
        for coeff_name, coeff in [("1/4", R(1,4)), ("1/2", R(1,2)), ("1", sp.Integer(1))]:
            T = cancel(coeff * ff_val / 24)
            diff = cancel(ric_val - T)
            if diff == 0:
                print(f"  [{label}]: ℛ = ({coeff_name})·FF/4!  ✓")
                break
        else:
            # Solve for coefficient
            alpha = sp.Symbol('alpha')
            eq = cancel(ric_val - alpha * ff_val / 24)
            sol = sp.solve(eq, alpha)
            # Also try with trace term
            beta = sp.Symbol('beta')
            eq2 = cancel(ric_val - alpha * ff_val / 24 - beta * norm_F5 * g_val)
            sol2 = sp.solve(eq2, [alpha, beta])
            print(f"  [{label}]: ℛ={ric_val}, FF/24={cancel(ff_val/24)}")
            print(f"    α (no trace): {sol}")
            print(f"    (α,β): {sol2}")
    else:
        # FF = 0, need trace
        if norm_F5 != 0 and g_val != 0:
            beta = cancel(ric_val / (norm_F5 * g_val))
            print(f"  [{label}]: FF=0, ℛ/(|F₅|²·g) = {beta}")
        else:
            print(f"  [{label}]: ℛ={ric_val}, FF=0, |F₅|²·g=0")

# The D3 self-dual formula: with ONLY electric half of F₅,
# the correct D3 equation should be: ℛ = α·FF^{elec}/4! - β·|F₅^{elec}|²·g
# For the FULL self-dual F₅:
# T^{SD} = (1/2)[FF^{SD}/4!] (no trace, since |F₅^{SD}|²=0)
# FF^{SD} = FF^{elec} + FF^{mag}
# For the D3 metric, FF^{mag}_{tt} = FF^{elec}_{tt} (by self-duality),
# so FF^{SD}_{tt} = 2·FF^{elec}_{tt}.
# Thus: T^{SD}_{tt} = (1/2)·2·FF^{elec}_{tt}/4! = FF^{elec}_{tt}/4!

# But: for transverse directions, FF^{elec}_{yy} comes from F_{yPQRS} terms.
# The electric F₅ only has components with one transverse index (from dH ∧ vol₄).
# So FF^{elec}_{yy} = Σ F_{yPQRS}F_y^{PQRS} where P,Q,R,S are wv+z.
# And FF^{mag}_{yy} = Σ F_{yP'Q'R'S'}F_y^{...} where P'... are transverse+z.
# These are NOT equal in general.

# Let me compute the MAGNETIC F₅ = *₁₂F₅^{elec} and check.
print("\n\nComputing magnetic F₅ = *₁₂(F₅^{elec})...")
from sugra import hodge_star
F5_mag = hodge_star(F5_d3, m_d3, signature=-1)  # 12-5=7 form... wait

# *₁₂(5-form) = 7-form. For D=12: *F₅ is a 7-form, not a 5-form.
# Self-duality in 12d works for 6-forms (D/2 = 6).
# The D3 F₅ is actually self-dual in 10d (*₁₀F₅ = F₅), not in 12d.
# In 12d: F₅ lives in the 10d subspace and *₁₂F₅ = F₇ = F₅ ∧ dz₁ ∧ dz₂.

# So the correct approach for D3 in 12d is:
# ℛ = T^{F₅}(self-dual, 10d)  uplifted to 12d
# The F₅ stress-energy in 10d (self-dual): T_{mn} = (1/2)FF^{SD}/4!
# = (1/2)·2·FF^{elec}/4! = FF^{elec}/4! (for wv/trans directions)
# Torus: flat, no F₅ components along torus → T_{za} = 0.

# In 12d: the F₅ doesn't have z-legs (no F_{...z...}), so FF^{F₅}_{za} = 0.
# And |F₅|² is nonzero but trace |F₅^{SD}|² = 0.
# The "effective" T in 12d: ℛ_{mn} = FF^{elec}_{mn}/4! for 10d directions.

# But above, the fit shows ratio 1/4 for worldvolume. Let me check:
# FF^{elec}_{tt}/24. If ℛ = α·FF/24, then α=1/4.
# But the expected coefficient for the self-dual doubling is α=1.
# So either FF is off by a factor of 4, or the coefficient is different.

# Let me just print the raw values:
print("\nRaw values for D3:")
for i12, label in labels_d3:
    ric_val = hf_d3.substitute(cancel(Ric_d3[i12, i12]))
    ff_val = hf_d3.substitute(cancel(FF_F5[i12, i12]))
    print(f"  [{label}]: ℛ={ric_val}, FF={ff_val}")


# ===================================================================
# D7-BRANE check
# ===================================================================
print("\n\n" + "="*70)
print("D7-BRANE (vacuum)")
print("="*70)

wv_d7 = list(sp.symbols('t x1 x2 x3 x4 x5 x6 x7', real=True))
trans_d7 = list(sp.symbols('y0:2', real=True))
coords_d7 = wv_d7 + [z1s, z2s] + trans_d7

hf_d7 = HarmonicFunction(transverse_coords=trans_d7)
H_d7 = sp.Function('H')(hf_d7.r_expr)

# D7 (C₀=0): flat₈ + H·(dy0²+dy1²) + diag(1/H, H)·dz²
g_d7 = sp.zeros(12, 12)
g_d7[0, 0] = -1
for i in range(1, 8):
    g_d7[i, i] = 1
g_d7[8, 8] = 1/H_d7   # z1
g_d7[9, 9] = H_d7      # z2
g_d7[10, 10] = H_d7    # y0
g_d7[11, 11] = H_d7    # y1

m_d7 = Metric(g_d7, coords_d7)

print("Computing D7 Ricci...")
Ric_d7 = m_d7.ricci_tensor(simplify_func=cancel)
print("Done.")

K_d7 = compute_K_scalar(sp.Matrix([[1/H_d7,0],[0,H_d7]]),
                         sp.Matrix([[H_d7,0],[0,1/H_d7]]),
                         coords_d7, m_d7, hf_d7)
print(f"K(D7) = {K_d7}")

print("\nD7 Ricci components (should be 0 for vacuum):")
labels_d7 = [(0,'t'), (8,'z1'), (9,'z2'), (10,'y0'), (11,'y1')]
for i, label in labels_d7:
    ric_val = hf_d7.substitute(cancel(Ric_d7[i, i]))
    print(f"  [{label}]: ℛ = {ric_val}")

# Check off-diagonal
ric_01 = hf_d7.substitute(cancel(Ric_d7[10, 11]))
print(f"  [y0,y1]: ℛ = {ric_01}")

# D7 uses harmonic condition in 2d: H'' + H'/r = 0
# After substitution, ℛ should be 0 if the Cauchy-Riemann conditions are met.
# But for diagonal M = diag(1/H, H) with C₀=0, the CR conditions are trivially
# satisfied (τ₁=0, τ₂=H). The remaining condition is ∇²H = 0 (harmonic in 2d).
# With hf.substitute, H'' → -H'/r (2d harmonic condition).
# So ℛ should come out 0.

# If ℛ ≠ 0, it means the D7 metric setup is incorrect OR hf doesn't
# simplify correctly for 2d.

# The 2d harmonic condition: H'' + (D_perp-1)H'/r = 0 with D_perp = 2:
# H'' + H'/r = 0. The HarmonicFunction class should handle this.

print(f"\n2d harmonic params: D_perp={len(trans_d7)}, condition: H'' + {len(trans_d7)-1}H'/r = 0")
print(f"hf.Hpp expression after substitute = -{len(trans_d7)-1}*Hp/r")

print("\n" + "="*70)
print("SUMMARY")
print("="*70)
