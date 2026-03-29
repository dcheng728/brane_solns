"""Experiment 31: Find a self-contained 12d stress-energy formula.

Goal: Find T_{MN} expressed purely in 12d quantities such that ℛ_{MN} = T_{MN}
for ALL 12 components, without reducing to 10d.

Strategy:
1. Compute all 12d building blocks: ℛ, FF(G₄), |G₄|², scalar kinetic terms from M
2. Express ℛ as a linear combination of these building blocks
3. Find universal coefficients that work for ALL blocks

Building blocks available from G₄ and M:
  (a) FF_{MN} = (G₄²)_{MN} = form_contraction(G₄, g₁₂)
  (b) |G₄|² = form_norm_squared(G₄, g₁₂)
  (c) S_{MN} = (1/4)Tr(∂_M M⁻¹ · ∂_N M)    [scalar kinetic tensor]
  (d) K = g^{MN}S_{MN}                         [scalar kinetic scalar]
  (e) M̃_{MN} = M embedded in 12d              [torus metric]

Ansatz: ℛ_{MN} = a·FF_{MN}/3! + b·|G₄|²·g_{MN} + c·S_{MN} + d·K·g_{MN} + e·K·M̃_{MN}

Test on: F1-string (primary), then NS5, D3.
"""
import sys; sys.path.insert(0, 'src')
import sympy as sp
from sympy import Rational as R, cancel, Matrix
from sugra import (HarmonicFunction, warped_product, Metric,
                   FormField, exterior_derivative,
                   form_contraction, form_norm_squared)

# ===================================================================
# F1-BRANE SETUP (Einstein frame, 12d)
# ===================================================================
print("="*70)
print("F1-BRANE: Computing all 12d building blocks")
print("="*70)

wv_coords = list(sp.symbols('t x1', real=True))
z1s, z2s = sp.symbols('z1 z2', real=True)
harmonic_coords = list(sp.symbols('y0:8', real=True))
coords_12 = wv_coords + [z1s, z2s] + harmonic_coords
y0 = harmonic_coords[0]

hf = HarmonicFunction(transverse_coords=harmonic_coords)
H = sp.Function('H')(hf.r_expr)

# F1 12d metric (Einstein frame):
# ds² = H^{-3/4} ds²_{1,1} + H^{1/2} dz₁² + H^{-1/2} dz₂² + H^{1/4} ds²_8
metric_12 = warped_product(
    warp_factors=[H**R(-3,4), H**R(1,2), H**R(-1,2), H**R(1,4)],
    block_dims=[2, 1, 1, 8],
    block_signatures=['lorentzian', 'euclidean', 'euclidean', 'euclidean'],
    coordinates=coords_12,
)

# Torus metric M = diag(H^{1/2}, H^{-1/2})  [m₁=g_{z1z1}, m₂=g_{z2z2}]
M = sp.Matrix([[H**R(1,2), 0], [0, H**R(-1,2)]])
Minv = sp.Matrix([[H**R(-1,2), 0], [0, H**R(1,2)]])

print("Computing 12d Ricci tensor...")
Ric_12 = metric_12.ricci_tensor(simplify_func=cancel)
print("Done.")

# G₄ = H₃ ∧ dz₂  →  C₃ with C_{t,x1,z2} = H^{-1}
C3 = FormField(rank=3, dim=12)
C3[(0, 1, 3)] = 1/H  # t=0, x1=1, z2=3
G4 = exterior_derivative(C3, coords_12)

print("Computing G₄ form contraction and norm...")
FF_G4 = form_contraction(G4, metric_12)
norm_G4 = form_norm_squared(G4, metric_12)
print("Done.")

# ===================================================================
# BUILDING BLOCK (c): Scalar kinetic tensor S_{MN}
# S_{MN} = (1/4) Tr(∂_M M⁻¹ · ∂_N M)
# ===================================================================

def compute_S_MN(M, Minv, coords, g_inv_diag, hf_obj):
    """Compute S_{MN} = (1/4)Tr(∂_M M⁻¹ · ∂_N M) for all 12d indices.

    For diagonal M with det M = 1:
    ∂M⁻¹/∂x_M = -M⁻¹ (∂M/∂x_M) M⁻¹
    Tr(∂_MM⁻¹ · ∂_NM) = -Tr(M⁻¹∂_MM · M⁻¹∂_NM)
    """
    D = len(coords)
    S = sp.zeros(D, D)

    for M_idx in range(D):
        # ∂_{x_M} M  (2x2 matrix)
        dM_M = M.diff(coords[M_idx]) if coords[M_idx] in M.free_symbols else sp.zeros(2,2)
        if dM_M == sp.zeros(2,2):
            continue

        dMinv_M = -Minv * dM_M * Minv  # ∂_M(M⁻¹)

        for N_idx in range(M_idx, D):
            dM_N = M.diff(coords[N_idx]) if coords[N_idx] in M.free_symbols else sp.zeros(2,2)
            if dM_N == sp.zeros(2,2):
                continue

            # S_{MN} = (1/4) Tr(dMinv_M · dM_N)
            product = dMinv_M * dM_N
            trace = product.trace()
            s_val = R(1, 4) * trace
            s_val = hf_obj.substitute(cancel(s_val))
            S[M_idx, N_idx] = s_val
            if M_idx != N_idx:
                S[N_idx, M_idx] = s_val

    return S

print("\nComputing scalar kinetic tensor S_{MN}...")
S_MN = compute_S_MN(M, Minv, coords_12, None, hf)
print("Done.")

# Scalar K = g^{MN} S_{MN}
K = sp.Integer(0)
for i in range(12):
    if S_MN[i,i] != 0:
        g_inv_ii = hf.substitute(cancel(metric_12.inv_matrix[i,i]))
        K += g_inv_ii * S_MN[i,i]
K = cancel(K)
print(f"\nScalar K = g^{{MN}}S_{{MN}} = {K}")

# ===================================================================
# TABULATE ALL BUILDING BLOCKS PER BLOCK
# ===================================================================
print("\n" + "="*70)
print("BUILDING BLOCKS PER DIAGONAL COMPONENT (F1)")
print("="*70)

# Representative indices: t=0, x1=1, z1=2, z2=3, y0=4, y1=5
block_indices = [(0, 't'), (1, 'x1'), (2, 'z1'), (3, 'z2'), (4, 'y0'), (5, 'y1')]

# Store values for linear system
data = {}

for idx, name in block_indices:
    ric = hf.substitute(cancel(Ric_12[idx, idx]))
    ff = hf.substitute(cancel(FF_G4[idx, idx]))
    g = hf.substitute(cancel(metric_12.matrix[idx, idx]))
    s = S_MN[idx, idx]
    n_g4 = hf.substitute(cancel(norm_G4))

    # Torus metric embedding M̃
    if idx == 2:  # z1
        m_tilde = hf.substitute(cancel(M[0, 0]))
    elif idx == 3:  # z2
        m_tilde = hf.substitute(cancel(M[1, 1]))
    else:
        m_tilde = sp.Integer(0)

    data[name] = {
        'ric': ric,          # ℛ_{ii}
        'ff': ff,            # FF_{ii} = (G₄²)_{ii}
        'g': g,              # g_{ii}
        's': s,              # S_{ii}
        'norm': n_g4,        # |G₄|²
        'K': K,              # scalar K
        'm_tilde': m_tilde,  # M̃_{ii}
    }

    print(f"\n[{name}] (idx={idx}):")
    print(f"  ℛ         = {ric}")
    print(f"  FF        = {ff}")
    print(f"  |G₄|²     = {n_g4}")
    print(f"  g         = {g}")
    print(f"  S         = {s}")
    print(f"  K         = {K}")
    print(f"  M̃         = {m_tilde}")

# ===================================================================
# TEST CANDIDATE FORMULAS
# ===================================================================
print("\n" + "="*70)
print("CANDIDATE FORMULA TESTS")
print("="*70)

# --- Formula 1: Standard p-form T with λ=1/4 (10d value) ---
print("\n--- Formula 1: T = (1/2)[FF/3! - (1/4)|G₄|²g] ---")
for name in ['t', 'x1', 'z1', 'z2', 'y0', 'y1']:
    d = data[name]
    T = R(1,2) * (d['ff']/6 - R(1,4)*d['norm']*d['g'])
    T = cancel(T)
    diff = cancel(d['ric'] - T)
    status = "✓" if diff == 0 else f"✗ residual = {diff}"
    print(f"  [{name}]: {status}")

# --- Formula 2: T = (1/2)[FF/3! - (1/4)|G₄|²g] - K·M̃ ---
print("\n--- Formula 2: T = (1/2)[FF/3! - (1/4)|G₄|²g] - K·M̃ ---")
for name in ['t', 'x1', 'z1', 'z2', 'y0', 'y1']:
    d = data[name]
    T = R(1,2) * (d['ff']/6 - R(1,4)*d['norm']*d['g']) - d['K']*d['m_tilde']
    T = cancel(T)
    diff = cancel(d['ric'] - T)
    status = "✓" if diff == 0 else f"✗ residual = {diff}"
    print(f"  [{name}]: {status}")

# --- Formula 3: T = (1/2)[FF/3! - (3/10)|G₄|²g] + 2S (standard 12d action) ---
print("\n--- Formula 3: T = (1/2)[FF/3! - (3/10)|G₄|²g] + 2S (D=12 action) ---")
for name in ['t', 'x1', 'z1', 'z2', 'y0', 'y1']:
    d = data[name]
    T = R(1,2) * (d['ff']/6 - R(3,10)*d['norm']*d['g']) + 2*d['s']
    T = cancel(T)
    diff = cancel(d['ric'] - T)
    status = "✓" if diff == 0 else f"✗ residual = {diff}"
    print(f"  [{name}]: {status}")

# --- Formula 4: General fit ---
# ℛ = a·FF/3! + b·|G₄|²·g + c·S + d·K·g + e·K·M̃
# Solve for (a, b, c, d, e) using 6 blocks
print("\n--- Formula 4: General linear fit ---")
print("  ℛ = a·FF/6 + b·|G₄|²·g + c·S + d·K·g + e·K·M̃")
print()

# Build the linear system
a, b, c, d_coeff, e = sp.symbols('a b c d e')
equations = []
for name in ['t', 'x1', 'z1', 'z2', 'y0', 'y1']:
    dd = data[name]
    lhs = dd['ric']
    rhs = (a * dd['ff']/6 + b * dd['norm'] * dd['g'] + c * dd['s']
           + d_coeff * dd['K'] * dd['g'] + e * dd['K'] * dd['m_tilde'])
    rhs = sp.expand(rhs)
    eq = cancel(lhs - rhs)
    equations.append(eq)
    print(f"  [{name}]: ℛ - T = {eq}")

# Try to solve
print("\n  Solving system for (a, b, c, d, e)...")
try:
    sol = sp.solve(equations, [a, b, c, d_coeff, e], dict=True)
    if sol:
        for s in sol:
            print(f"  Solution: {s}")
    else:
        print("  No solution found.")
except Exception as ex:
    print(f"  Solve failed: {ex}")

# --- Formula 5: Reduced ansatz ---
# From analysis: S_{MN} = 0 for worldvolume (t, x1) and torus (z1, z2)
# S_{MN} ≠ 0 only for transverse (y_k).
# So c·S only affects y-blocks. Try: a·FF/6 + b·|G₄|²g + e·K·M̃
print("\n--- Formula 5: Reduced ansatz (no S, no K·g) ---")
print("  ℛ = a·FF/6 + b·|G₄|²·g + e·K·M̃")
eqs5 = []
for name in ['t', 'x1', 'z1', 'z2', 'y0']:
    dd = data[name]
    lhs = dd['ric']
    rhs = a * dd['ff']/6 + b * dd['norm'] * dd['g'] + e * dd['K'] * dd['m_tilde']
    eq = cancel(lhs - rhs)
    eqs5.append(eq)

try:
    sol5 = sp.solve(eqs5, [a, b, e], dict=True)
    if sol5:
        for s in sol5:
            print(f"  Solution: {s}")
            # Verify y1 (anisotropic block)
            dd = data['y1']
            rhs_y1 = s[a]*dd['ff']/6 + s[b]*dd['norm']*dd['g'] + s[e]*dd['K']*dd['m_tilde']
            diff_y1 = cancel(dd['ric'] - rhs_y1)
            print(f"  Verify y1: residual = {diff_y1}")
    else:
        print("  No solution.")
except Exception as ex:
    print(f"  Solve failed: {ex}")

# --- Formula 6: Per-block effective λ ---
# ℛ = (1/2)[FF/3! - λ|G₄|²g] → solve for λ per block
print("\n--- Formula 6: Effective λ per block ---")
print("  ℛ = (1/2)[FF/3! - λ|G₄|²·g]")
lam = sp.Symbol('lambda')
for name in ['t', 'x1', 'z1', 'z2', 'y0', 'y1']:
    dd = data[name]
    eq = cancel(dd['ric'] - R(1,2)*(dd['ff']/6 - lam*dd['norm']*dd['g']))
    sol_lam = sp.solve(eq, lam)
    if sol_lam:
        for sl in sol_lam:
            sl = cancel(sl)
            print(f"  [{name}]: λ = {sl}")
    else:
        print(f"  [{name}]: no solution for λ")

# ===================================================================
# INVERSE APPROACH: What exactly is the residual?
# ===================================================================
print("\n" + "="*70)
print("RESIDUAL ANALYSIS: ℛ - T(λ=1/4)")
print("="*70)

print("\nResidual D = ℛ - (1/2)[FF/3! - (1/4)|G₄|²g]:")
for name in ['t', 'x1', 'z1', 'z2', 'y0', 'y1']:
    dd = data[name]
    T_quarter = R(1,2) * (dd['ff']/6 - R(1,4)*dd['norm']*dd['g'])
    D_val = cancel(dd['ric'] - T_quarter)

    # Express residual / K / M̃ if possible
    if dd['m_tilde'] != 0 and dd['K'] != 0:
        ratio = cancel(D_val / (dd['K'] * dd['m_tilde']))
        print(f"  [{name}]: D = {D_val},  D/(K·M̃) = {ratio}")
    elif dd['g'] != 0 and dd['K'] != 0:
        ratio = cancel(D_val / (dd['K'] * dd['g']))
        print(f"  [{name}]: D = {D_val},  D/(K·g) = {ratio}")
    else:
        print(f"  [{name}]: D = {D_val}")

# Also check: is the residual proportional to S_{MN}?
print("\nResidual vs S_{MN}:")
for name in ['y0', 'y1']:
    dd = data[name]
    T_quarter = R(1,2) * (dd['ff']/6 - R(1,4)*dd['norm']*dd['g'])
    D_val = cancel(dd['ric'] - T_quarter)
    if dd['s'] != 0:
        ratio = cancel(D_val / dd['s'])
        print(f"  [{name}]: D/S = {ratio}")

print("\n" + "="*70)
print("DONE")
print("="*70)
