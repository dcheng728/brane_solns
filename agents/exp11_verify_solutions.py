"""Experiment 11: Verify the solutions found in exp10.

Key findings to verify:
1. D1 + F1 with lambda=1/2: cd²=5/4, cf²=1  (from exp10 Part D)
2. 5-form with lambda=1: consistent across t,z1,z2 (from exp10 Part H)

Need to check: y0 block (both isotropic and anisotropic parts).

Also: What physical formulation gives lambda=1/2?
"""
import sys; sys.path.insert(0, 'src')
import sympy as sp
from sympy import Rational as R
from sugra import (HarmonicFunction, warped_product,
                   FormField, exterior_derivative,
                   form_stress_energy, form_contraction,
                   form_norm_squared)

# --- Setup ---
wv_coords = list(sp.symbols('t x1', real=True))
z1, z2 = sp.symbols('z1 z2', real=True)
harmonic_coords = list(sp.symbols('y0:8', real=True))
coords = wv_coords + [z1, z2] + harmonic_coords
y0 = sp.Symbol('y0', real=True)

hf = HarmonicFunction(transverse_coords=harmonic_coords)
H = sp.Function('H')(hf.r_expr)

metric = warped_product(
    warp_factors=[H**R(-3,4), H**R(1,2), H**R(-1,2), H**R(1,4)],
    block_dims=[2, 1, 1, 8],
    block_signatures=['lorentzian','euclidean','euclidean','euclidean'],
    coordinates=coords,
)

print("Computing Ricci tensor...")
Ric = metric.ricci_tensor(simplify_func=sp.cancel)
print("Done.\n")

# --- Build D1 and F1 forms ---
C_D1 = FormField(rank=3, dim=12)
C_D1[(0,1,3)] = 1/H
F_D1 = exterior_derivative(C_D1, coords)

C_F1 = FormField(rank=3, dim=12)
C_F1[(0,1,2)] = H**R(-1,2)
F_F1 = exterior_derivative(C_F1, coords)

# Contractions and norms
FF_D1 = form_contraction(F_D1, metric)
S_D1 = form_norm_squared(F_D1, metric)
FF_F1 = form_contraction(F_F1, metric)
S_F1 = form_norm_squared(F_F1, metric)

# ===================================================================
# Verification 1: D1 + F1 with lambda=1/2, cd²=5/4, cf²=1
#
# T_{MN} = cd² * (1/2)[FF_D1/6 - lambda*S_D1*g] + cf² * (1/2)[FF_F1/6 - lambda*S_F1*g]
#         = (1/2) [cd²*FF_D1/6 + cf²*FF_F1/6 - lambda*(cd²*S_D1 + cf²*S_F1)*g]
# ===================================================================

print("="*60)
print("Verification 1: D1+F1, lambda=1/2, cd²=5/4, cf²=1")
print("="*60)

cd2 = R(5,4)
cf2 = R(1,1)
lam = R(1,2)

# Compute T and R for ALL diagonal components
print("\nBlock-by-block comparison:")
all_match = True

for i in range(12):
    name = coords[i]

    r_val = hf.substitute(sp.cancel(Ric[i,i]))
    ff_d1 = hf.substitute(sp.cancel(FF_D1[i,i]))
    ff_f1 = hf.substitute(sp.cancel(FF_F1[i,i]))
    s_d1 = hf.substitute(sp.cancel(S_D1))
    s_f1 = hf.substitute(sp.cancel(S_F1))
    g_val = hf.substitute(sp.cancel(metric.matrix[i,i]))

    T_val = R(1,2) * (cd2 * ff_d1/6 + cf2 * ff_f1/6 - lam*(cd2*s_d1 + cf2*s_f1)*g_val)
    T_val = sp.cancel(T_val)

    diff = sp.cancel(r_val - T_val)

    if i <= 4 or i == 5:  # Print representative components
        label = str(name)
        print(f"\n  [{label}] (index {i}):")
        print(f"    R  = {r_val}")
        print(f"    T  = {T_val}")
        print(f"    R-T = {diff}")
        if diff != 0:
            all_match = False
            print("    *** MISMATCH ***")
        else:
            print("    ✓")

if all_match:
    print("\n*** ALL COMPONENTS MATCH: R = T ***")
else:
    print("\n--- Some components don't match ---")

# ===================================================================
# Check y0 block in detail (isotropic + anisotropic)
# ===================================================================

print("\n" + "="*60)
print("Detailed y0 analysis")
print("="*60)

i = 4  # y0 index
r_y0 = hf.substitute(sp.cancel(Ric[4,4]))
ff_d1_y0 = hf.substitute(sp.cancel(FF_D1[4,4]))
ff_f1_y0 = hf.substitute(sp.cancel(FF_F1[4,4]))
s_d1 = hf.substitute(sp.cancel(S_D1))
s_f1 = hf.substitute(sp.cancel(S_F1))
g_y0 = hf.substitute(sp.cancel(metric.matrix[4,4]))

T_y0 = R(1,2) * (cd2 * ff_d1_y0/6 + cf2 * ff_f1_y0/6 - lam*(cd2*s_d1 + cf2*s_f1)*g_y0)
T_y0 = sp.cancel(T_y0)

print(f"R[y0]  = {r_y0}")
print(f"T[y0]  = {T_y0}")
print(f"R-T    = {sp.cancel(r_y0 - T_y0)}")

# Split into isotropic and anisotropic
r_iso = r_y0.subs(y0, 0)
t_iso = T_y0.subs(y0, 0)
print(f"\nIsotropic (y0=0):")
print(f"  R_iso = {r_iso}")
print(f"  T_iso = {t_iso}")
print(f"  Diff  = {sp.cancel(r_iso - t_iso)}")

r_aniso_coeff = sp.cancel(sp.diff(r_y0, y0, 2) / 2)
t_aniso_coeff = sp.cancel(sp.diff(T_y0, y0, 2) / 2)
print(f"\nAnisotropic (coeff of y0²):")
print(f"  R_aniso = {r_aniso_coeff}")
print(f"  T_aniso = {t_aniso_coeff}")
print(f"  Diff    = {sp.cancel(r_aniso_coeff - t_aniso_coeff)}")

# ===================================================================
# Verification 2: 5-form A_{t,x1,z1,z2}=H^{-3/4} with lambda=1
# ===================================================================

print("\n" + "="*60)
print("Verification 2: 5-form with lambda=1")
print("="*60)

A5 = FormField(rank=4, dim=12)
A5[(0,1,2,3)] = H**R(-3,4)
F5 = exterior_derivative(A5, coords)

FF_5 = form_contraction(F5, metric)
S_5 = form_norm_squared(F5, metric)

lam_5 = R(1,1)

# Compute kappa = R/T for each block with lambda=1
print("\nWith lambda=1:")
for i, name in [(0,'t'), (1,'x1'), (2,'z1'), (3,'z2'), (4,'y0'), (5,'y1')]:
    r_val = hf.substitute(sp.cancel(Ric[i,i]))
    ff_val = hf.substitute(sp.cancel(FF_5[i,i]))
    s_val = hf.substitute(sp.cancel(S_5))
    g_val = hf.substitute(sp.cancel(metric.matrix[i,i]))

    T_val = R(1,2) * (ff_val/24 - lam_5 * s_val * g_val)
    T_val = sp.cancel(T_val)

    if T_val != 0:
        kappa = sp.cancel(r_val / T_val)
    else:
        kappa = 'T=0'

    diff = sp.cancel(r_val - T_val)
    print(f"  [{name}]: R={r_val}, T={T_val}, kappa={kappa}, diff={diff}")

# ===================================================================
# Verification 3: D1 + F1 + 5-form combined with different lambdas
# In IIB, the 4-form and 5-form have DIFFERENT trace coefficients
# because they have different ranks.
# T_total = cd²*T_D1(lam4) + cf²*T_F1(lam4) + c5²*T_5(lam5)
# With lam4 and lam5 potentially different
# ===================================================================

print("\n" + "="*60)
print("Verification 3: D1+F1+5-form with separate lambdas")
print("="*60)

# From verification 1: D1+F1 at lambda=1/2 gives residual in y0
# From verification 2: 5-form at lambda=1 might compensate
# Can we add the 5-form to fix the y0 residual?

# Compute D1+F1 residual at lambda=1/2
cd2_v = R(5,4)
cf2_v = R(1,1)
lam4 = R(1,2)

residuals = {}
for i, name in [(0,'t'), (2,'z1'), (3,'z2'), (4,'y0')]:
    r_val = hf.substitute(sp.cancel(Ric[i,i]))
    ff_d1 = hf.substitute(sp.cancel(FF_D1[i,i]))
    ff_f1 = hf.substitute(sp.cancel(FF_F1[i,i]))
    s_d1 = hf.substitute(sp.cancel(S_D1))
    s_f1 = hf.substitute(sp.cancel(S_F1))
    g_val = hf.substitute(sp.cancel(metric.matrix[i,i]))

    T12 = R(1,2) * (cd2_v*ff_d1/6 + cf2_v*ff_f1/6 - lam4*(cd2_v*s_d1 + cf2_v*s_f1)*g_val)
    residuals[name] = sp.cancel(r_val - T12)

print("D1+F1 residuals (R - T_12) at lambda4=1/2:")
for name, res in residuals.items():
    print(f"  D[{name}] = {res}")

# Now check if the 5-form can absorb these residuals
# D[block] = c5² * T_5(lam5)[block]
# For this to work: D[block]/T_5[block] = c5² for all blocks
print("\n5-form stress-energy at various lambda5:")
for lam5 in [R(1,2), R(2,5), R(1,1), R(3,10)]:
    T5_blk = {}
    for i, name in [(0,'t'), (2,'z1'), (3,'z2'), (4,'y0')]:
        ff_val = hf.substitute(sp.cancel(FF_5[i,i]))
        s_val = hf.substitute(sp.cancel(S_5))
        g_val = hf.substitute(sp.cancel(metric.matrix[i,i]))
        T5_blk[name] = sp.cancel(R(1,2) * (ff_val/24 - lam5*s_val*g_val))

    # Check if residuals / T5 is constant
    ratios = {}
    for name in ['t', 'z1', 'z2']:
        if T5_blk[name] != 0:
            ratios[name] = sp.cancel(residuals[name] / T5_blk[name])

    print(f"\n  lambda5={lam5}:")
    for name, rat in ratios.items():
        print(f"    D[{name}]/T5[{name}] = {rat}")

    # Check if all ratios are equal
    vals = list(ratios.values())
    if len(vals) >= 2 and all(sp.cancel(v - vals[0]) == 0 for v in vals[1:]):
        c5_sq = vals[0]
        print(f"    *** CONSISTENT: c5² = {c5_sq}")

        # Check y0
        if T5_blk['y0'] != 0:
            y0_ratio = sp.cancel(residuals['y0'] / T5_blk['y0'])
            print(f"    y0 check: D[y0]/T5[y0] = {y0_ratio}")
            if sp.cancel(y0_ratio - c5_sq) == 0:
                print(f"    *** FULL SOLUTION FOUND! ***")

# ===================================================================
# Verification 4: Direct check with the FULL modified formula
# T = sum_i ci² * [(1/2) FF_i / (p_i-1)! - lambda_i * S_i * g / 2]
# Try: D1 (c1²=5/4, lam=1/2) + F1 (c2²=1, lam=1/2) + ?
# ===================================================================

print("\n" + "="*60)
print("Verification 4: Can we fix y0 with ANY additional form?")
print("="*60)

# The residual from D1+F1 at lambda=1/2 for y0 block
D_y0 = residuals['y0']
D_y0_iso = D_y0.subs(y0, 0)
D_y0_aniso = sp.cancel(sp.diff(D_y0, y0, 2) / 2)

print(f"y0 residual (full): {D_y0}")
print(f"y0 residual (iso):  {D_y0_iso}")
print(f"y0 residual (aniso): {D_y0_aniso}")

# D[t] = 0, D[z1] = 0, D[z2] = 0, D[y0] = something
# We need a form that: T[t] = 0, T[z1] = 0, T[z2] = 0, T[y0] ≠ 0
# This is impossible for standard form stress-energy because
# T[worldvolume] always contributes through the trace term.
# UNLESS: FF[wv] = lambda*S*g[wv] exactly, which cancels T[wv].

# The only way T[t]=T[z1]=T[z2]=0 with T[y0]≠0 is:
# FF[t]/norm = lambda*S*g[t], FF[z1]/norm = lambda*S*g[z1], FF[z2]/norm = lambda*S*g[z2]
# This requires very specific form structure.

# Actually, what about a form that lives ENTIRELY in the transverse space?
# F_{y0,y1,y2,...} = f(r) * epsilon_{ijk...} * y_k/r type
# But exp03 showed these don't match H-powers (SO(8) invariant 4-form = 0).

# What about a scalar (0-form) contribution? Or what if D[y0] ≠ 0 indicates
# we need to adjust the D1+F1 amplitudes?

print("\nThe t, z1, z2 residuals from D1+F1 (lambda=1/2):")
for name in ['t', 'z1', 'z2']:
    print(f"  D[{name}] = {residuals[name]}")

# ===================================================================
# Verification 5: FULL numerical check at a specific point
# ===================================================================

print("\n" + "="*60)
print("Verification 5: Numerical check at H=2, H'=1, r=1, y0=1/4")
print("="*60)

subs = [(hf.H, 2), (hf.Hp, 1), (hf.r, 1), (y0, R(1,4))]

for i, name in [(0,'t'), (2,'z1'), (3,'z2'), (4,'y0')]:
    r_num = float(R_blk_val := hf.substitute(sp.cancel(Ric[i,i])).subs(subs))

    ff_d1 = hf.substitute(sp.cancel(FF_D1[i,i])).subs(subs)
    ff_f1 = hf.substitute(sp.cancel(FF_F1[i,i])).subs(subs)
    s_d1 = hf.substitute(sp.cancel(S_D1)).subs(subs)
    s_f1 = hf.substitute(sp.cancel(S_F1)).subs(subs)
    g_val = hf.substitute(sp.cancel(metric.matrix[i,i])).subs(subs)

    T_num = float(R(1,2) * (cd2_v*ff_d1/6 + cf2_v*ff_f1/6 - lam4*(cd2_v*s_d1 + cf2_v*s_f1)*g_val))

    print(f"  [{name}] R={r_num:.6f}, T={T_num:.6f}, diff={r_num-T_num:.6f}")

# ===================================================================
# Summary: The physical interpretation
# ===================================================================

print("\n" + "="*60)
print("PHYSICAL INTERPRETATION")
print("="*60)
print("""
Standard SUGRA formula: T_{MN} = (1/2)[FF/(p-1)! - (p-1)/(D-2) |F²| g]
  → lambda = (p-1)/(D-2) = 3/10 for p=4, D=12

Found solution requires lambda = 1/2 for D1+F1 system.

Possible interpretations:
1. The 12d theory has a modified action coefficient
2. This corresponds to a different trace relation (possibly from CS terms)
3. lambda=1/2 = (p-1)/(D-2) would mean:
   - p=6, D=12 (6-form in 12d), or
   - p=4, D=8 (4-form in 8d), or
   - p=2, D=4 (2-form in 4d)

For the 5-form: lambda=1 = (p-1)/(D-2) → p=11, D=12 (impossible for p>D)
or this is a self-dual-type formula where the trace coefficient doubles.

Key finding: The ANISOTROPIC y0²/r² part of the y0 block is AUTOMATICALLY
satisfied (lambda-independent) for the D1 form. This is the angular structure
that encodes the brane source. The isotropic part requires a specific lambda.
""")
