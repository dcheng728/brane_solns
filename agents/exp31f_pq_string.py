"""Experiment 31f: Verify enhanced-metric formula on (1,1)-string (non-diagonal M).

The (1,1)-string is obtained from F1 via SL(2,Z) with Λ = [[1,0],[1,1]].
It has a NON-DIAGONAL torus metric M = Λ·M₀·Λ^T.

This is the strongest test: non-diagonal M, non-trivial G₄ doublet.
Uses 2 transverse coordinates for tractability with full (non-diagonal) Ricci.

FORMULA:
  ℛ_{MN} = (1/2) [(G₄²)_{MN}/3! - (1/4)|G₄|² ĝ_{MN}]
  ĝ_{MN} = g_{MN} + M̃_{MN}
"""
import sys; sys.path.insert(0, 'src')
import sympy as sp
from sympy import Rational as R, cancel
from sugra import (HarmonicFunction, Metric,
                   FormField, exterior_derivative,
                   form_contraction, form_norm_squared)

# ===================================================================
# (1,1)-string setup (following exp29 pattern)
# ===================================================================
print("="*70)
print("(1,1)-STRING (non-diagonal torus)")
print("="*70)

wv = list(sp.symbols('t x1', real=True))
z1s, z2s = sp.symbols('z1 z2', real=True)
trans = list(sp.symbols('y0:2', real=True))  # 2 transverse for speed
coords = wv + [z1s, z2s] + trans  # 6d total

hf = HarmonicFunction(transverse_coords=trans)
H = sp.Function('H')(hf.r_expr)

# F1 torus metric: M₀ = diag(H^{1/2}, H^{-1/2})
M0 = sp.Matrix([[H**R(1,2), 0], [0, H**R(-1,2)]])

# SL(2,Z) transformation: Λ = [[1,0],[1,1]]
Lambda = sp.Matrix([[1, 0], [1, 1]])
M = cancel(Lambda * M0 * Lambda.T)
# M = [[H^{1/2}, H^{1/2}], [H^{1/2}, H^{1/2}+H^{-1/2}]]
print(f"M = {M}")
print(f"det(M) = {cancel(M.det())}")

# String-frame metric (SL(2,Z) invariant): g^S = H^{-1}ds²_{wv} + ds²_trans
# Einstein-frame metric (NOT SL(2,Z) invariant) depends on τ₂.
# For (1,1)-string: τ₂ = det(M₀)^{1/2} / M_{11}^{...} ... complicated.
# Actually, the 12d metric is ds²₁₂ = g^{10d}_{mn}dx^mdx^n + M_{ab}dz^adz^b
# For the (1,1)-string in STRING FRAME:
#   g^S = H^{-1}ds²_{1,1} + ds²_trans  (SL(2,Z) invariant)

# Build 6d metric (string frame + torus)
g_6 = sp.zeros(6, 6)
g_6[0, 0] = -1/H   # t (string frame)
g_6[1, 1] = 1/H    # x1
g_6[2, 2] = M[0, 0]  # z1z1
g_6[2, 3] = M[0, 1]  # z1z2
g_6[3, 2] = M[1, 0]  # z2z1
g_6[3, 3] = M[1, 1]  # z2z2
g_6[4, 4] = 1        # y0
g_6[5, 5] = 1        # y1

metric = Metric(g_6, coords)

# The Metric class may not handle non-diagonal metrics correctly for Ricci
# (the bug from exp29). Force full computation.
# Monkey-patch: set is_diagonal to False
metric._is_diagonal = False

print("\nComputing 6d Ricci (full, non-diagonal)...")
Ric = metric.ricci_tensor(simplify_func=cancel)
print("Done.")

# G₄ for (1,1)-string:
# From exp29: G₄ = c_a F₃^a ∧ dz_a where c = (-1, 1) (charge vector for (1,1))
# G₄ = -H₃ ∧ dz₁ + H₃ ∧ dz₂
# where H₃ = dB₂ with B_{t,x1} = H^{-1}
# In 6d: C₃ with:
#   C_{t,x1,z1} = -H^{-1}  (from -H₃ ∧ dz₁)
#   C_{t,x1,z2} = +H^{-1}  (from +H₃ ∧ dz₂)

# But G₄ = dC₃, so we need to be careful. Actually:
# G₄ = -dB₂ ∧ dz₁ + dB₂ ∧ dz₂ = dB₂ ∧ (-dz₁ + dz₂)
# As a 4-form in 6d:
# G₄_{t,x1,z1,yk} = -∂_{yk}(H^{-1})  (from -H₃∧dz₁)
# G₄_{t,x1,z2,yk} = +∂_{yk}(H^{-1})  (from +H₃∧dz₂)

# Build via exterior derivative of C₃ = -B₂∧dz₁ + B₂∧dz₂
C3 = FormField(rank=3, dim=6)
C3[(0, 1, 2)] = -1/H  # C_{t,x1,z1} = -H^{-1}
C3[(0, 1, 3)] = 1/H   # C_{t,x1,z2} = +H^{-1}

G4 = exterior_derivative(C3, coords)

print("\nG₄ non-zero components:")
for idx, val in G4.nonzero_components.items():
    sval = hf.substitute(cancel(val))
    cnames = [str(coords[i]) for i in idx]
    print(f"  G₄[{','.join(cnames)}] = {sval}")

# For non-diagonal metrics, form_contraction assumes diagonal metric!
# Need to use a manual contraction.
print("\nComputing form contraction (manual, non-diagonal metric)...")

g_inv = metric.inv_matrix

def manual_form_contraction(form, g_inv, dim, hf_obj):
    """FF_{MN} = Σ_{P,Q,R} F_{MPQR} F_N^{PQR} for a 4-form."""
    FF = sp.zeros(dim, dim)
    rank = form.rank
    # Get all nonzero components
    nz = form.nonzero_components

    for M in range(dim):
        for N in range(M, dim):
            val = sp.Integer(0)
            # Sum over all (P,Q,R) triples
            for P in range(dim):
                for Q in range(P+1, dim):
                    for R in range(Q+1, dim):
                        # F_{MPQR}
                        f_lower = form[(M, P, Q, R)]
                        if f_lower == 0:
                            continue
                        # F_N^{PQR} = g^{PP'}g^{QQ'}g^{RR'} F_{NP'Q'R'}
                        # For the upper-index contraction, need full g^{-1}
                        f_upper = sp.Integer(0)
                        for Pp in range(dim):
                            if g_inv[P, Pp] == 0:
                                continue
                            for Qp in range(dim):
                                if g_inv[Q, Qp] == 0:
                                    continue
                                for Rp in range(dim):
                                    if g_inv[R, Rp] == 0:
                                        continue
                                    f_NpPQ = form[(N, Pp, Qp, Rp)]
                                    if f_NpPQ != 0:
                                        f_upper += g_inv[P,Pp]*g_inv[Q,Qp]*g_inv[R,Rp]*f_NpPQ

                        val += f_lower * f_upper
            val = hf_obj.substitute(cancel(val))
            FF[M, N] = val
            if M != N:
                FF[N, M] = val
    return FF


def manual_form_norm(form, g_inv, dim, hf_obj):
    """(1/p!) F_{M...} F^{M...}"""
    rank = form.rank
    val = sp.Integer(0)
    for combo in __import__('itertools').combinations(range(dim), rank):
        f_lower = form[combo]
        if f_lower == 0:
            continue
        # Raise all indices
        f_upper = sp.Integer(0)
        # This is expensive; use the relation |F|² = (1/p!) Tr(FF · g⁻¹) approach
        # Actually: |F|² = (1/p!) Σ_{i₁<...<iₚ} F_{i₁...iₚ} F^{i₁...iₚ}
        # F^{i₁...iₚ} = Σ g^{i₁j₁}...g^{iₚjₚ} F_{j₁...jₚ}

        i1, i2, i3, i4 = combo
        for j1 in range(dim):
            if g_inv[i1,j1] == 0: continue
            for j2 in range(dim):
                if g_inv[i2,j2] == 0: continue
                for j3 in range(dim):
                    if g_inv[i3,j3] == 0: continue
                    for j4 in range(dim):
                        if g_inv[i4,j4] == 0: continue
                        fj = form[(j1,j2,j3,j4)]
                        if fj != 0:
                            f_upper += g_inv[i1,j1]*g_inv[i2,j2]*g_inv[i3,j3]*g_inv[i4,j4]*fj
        val += f_lower * f_upper

    # Factor: no 1/p! because we sum over sorted indices only = already 1 representative per p! perms
    # Actually: |F|² = (1/p!) Σ_{ALL i} F·F = Σ_{sorted} F·F^{same sorted}
    # The sum over sorted combos with F·F^{sorted} gives |F|² directly.
    return hf_obj.substitute(cancel(val))


FF_G4 = manual_form_contraction(G4, g_inv, 6, hf)
norm_G4 = manual_form_norm(G4, g_inv, 6, hf)
print(f"|G₄|² = {norm_G4}")

# Enhanced metric: ĝ = g + M̃
# M̃ in 6d: indices 2,3 (z1,z2)
g_hat = sp.zeros(6, 6)
for i in range(6):
    for j in range(i, 6):
        g_val = hf.substitute(cancel(metric.matrix[i, j]))
        m_val = sp.Integer(0)
        if i in [2, 3] and j in [2, 3]:
            m_val = hf.substitute(cancel(M[i-2, j-2]))
        gh = cancel(g_val + m_val)
        g_hat[i, j] = gh
        if i != j:
            g_hat[j, i] = gh

# Verify formula
print(f"\n(1,1)-string: ℛ = (1/2)[FF/3! - (1/4)|G₄|²·ĝ]")
print("-"*60)

labels = [(0,'t'), (1,'x1'), (2,'z1'), (3,'z2'), (4,'y0'), (5,'y1')]
all_pass = True

for i, label in labels:
    ric_val = hf.substitute(cancel(Ric[i, i]))
    ff_val = FF_G4[i, i]
    gh_val = g_hat[i, i]

    T = cancel(R(1,2) * (ff_val/6 - R(1,4)*norm_G4*gh_val))
    diff = cancel(ric_val - T)
    status = "✓" if diff == 0 else f"✗ residual={diff}"
    if diff != 0:
        all_pass = False
    print(f"  [{label}]: ℛ={ric_val}, T={T}, diff={diff}  {status}")

# Off-diagonal: z1z2
ric_z12 = hf.substitute(cancel(Ric[2, 3]))
ff_z12 = FF_G4[2, 3]
gh_z12 = g_hat[2, 3]
T_z12 = cancel(R(1,2)*(ff_z12/6 - R(1,4)*norm_G4*gh_z12))
diff_z12 = cancel(ric_z12 - T_z12)
status = "✓" if diff_z12 == 0 else f"✗ residual={diff_z12}"
if diff_z12 != 0:
    all_pass = False
print(f"  [z1z2]: ℛ={ric_z12}, T={T_z12}, diff={diff_z12}  {status}")

# Off-diagonal: y0y1
ric_y01 = hf.substitute(cancel(Ric[4, 5]))
ff_y01 = FF_G4[4, 5]
gh_y01 = g_hat[4, 5]  # = 0
T_y01 = cancel(R(1,2)*(ff_y01/6 - R(1,4)*norm_G4*gh_y01))
diff_y01 = cancel(ric_y01 - T_y01)
status = "✓" if diff_y01 == 0 else f"✗ residual={diff_y01}"
if diff_y01 != 0:
    all_pass = False
print(f"  [y0y1]: ℛ={ric_y01}, T={T_y01}, diff={diff_y01}  {status}")

print(f"\n  {'ALL PASS ✓' if all_pass else 'SOME FAIL ✗'}")

if all_pass:
    print("""
★★★ FORMULA VERIFIED FOR NON-DIAGONAL TORUS ★★★

The enhanced-metric formula
  ℛ_{MN} = (1/2)[(G₄²)_{MN}/3! - (1/4)|G₄|²(g_{MN} + M̃_{MN})]
works for the (1,1)-string with non-diagonal M, confirming SL(2,Z) covariance.
""")
