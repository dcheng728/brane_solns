"""Experiment 09: Non-closed forms and modified field strengths.

Direction #2 and #6 from notes: What if F is not simply dC?

In IIB SUGRA, the self-dual 5-form is F_5 = dC_4 + B_2 ∧ F_3.
The Chern-Simons term means F is NOT the exterior derivative of a single potential.
Similarly, F_3 = dC_2 and H_3 = dB_2 are closed individually, but the
modified field strength F̃_3 = F_3 - C_0 * H_3 is not.

In the 12d uplift, the modified field strengths might introduce cross-terms
between different form components.

Strategy:
1. Specify F directly (not as dC) with free coefficient parameters
2. Compute T_{MN} and match to R_{MN}
3. Check if solutions exist that aren't exterior derivatives

Also: Test the specific ansatz F = dC_D1 + B ∧ dC_F1 (Chern-Simons-like).
"""
import sys; sys.path.insert(0, 'src')
import sympy as sp
from sympy import Rational as R
from sugra import (HarmonicFunction, warped_product,
                   FormField, exterior_derivative,
                   form_stress_energy, form_contraction,
                   form_norm_squared, Verifier)

# --- Setup ---
wv_coords = list(sp.symbols('t x1', real=True))
z1, z2 = sp.symbols('z1 z2', real=True)
harmonic_coords = list(sp.symbols('y0:8', real=True))
coords = wv_coords + [z1, z2] + harmonic_coords
D = 12

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

# Ricci blocks
R_blocks = {}
for i, name in [(0,'t'), (2,'z1'), (3,'z2'), (4,'y0')]:
    R_blocks[name] = hf.substitute(sp.cancel(Ric[i,i]))
    print(f"R[{name}] = {R_blocks[name]}")

# ===================================================================
# Part A: Inverse problem — what FF_{MN} and |F²| are required?
#
# T_{MN} = (1/2)[FF_{MN}/(p-1)! - (p-1)/(D-2)|F²|g_{MN}] = R_{MN}
# → FF_{MN}/(p-1)! = 2R_{MN} + (p-1)/(D-2)|F²|g_{MN}
#
# Trace: g^{MN}FF_{MN}/(p-1)! = |F²|/((p-1)!)  [by definition]
# → (1/p) |F²| = 2R + p(p-1)/(D-2)|F²|
# Wait, let me redo this properly.
#
# For a p-form: FF_{MN} = F_{M A1...A_{p-1}} F_N^{A1...A_{p-1}}
# |F²| = (1/p!) F_{A1...Ap} F^{A1...Ap} = g^{MN} FF_{MN} / (p!)... no.
# Actually: g^{MN} FF_{MN} = F^{M A1...} F_M^{A1...} = p! |F²|? No.
# g^{MN} FF_{MN} = g^{MN} F_{MA...} F_N^{A...} = p * |F|² * p!...
#
# Let me be precise. Define:
# S = |F|² = (1/p!) F_{M1...Mp} F^{M1...Mp}
# X_{MN} = FF_{MN} = F_{M A2...Ap} F_N^{A2...Ap}
# Then: g^{MN} X_{MN} = (p-1)! * p * S  ... actually:
# g^{MN} X_{MN} = g^{MN} F_{M A2...Ap} F_N^{A2...Ap}
#               = F^{N}_{A2...Ap} F_N^{A2...Ap}  (contracting M with g^{MN})
#               = p! * S  (by definition of S)
#
# So: Tr(X) = p! * S
#
# T_{MN} = (1/2) [X_{MN}/(p-1)! - (p-1)/(D-2) * S * g_{MN}] = R_{MN}
#
# Taking trace: Tr(T) = (1/2) [p!*S/(p-1)! - (p-1)*D/(D-2)*S]
#             = (1/2) S [p - p(p-1)/(D-2)]   ... wait
#             = (1/2) S [p - D(p-1)/(D-2)]
#
# R_trace = g^{MN} R_{MN} = Ricci scalar R
#
# (1/2) S [p - D(p-1)/(D-2)] = R
# S = 2R / [p - D(p-1)/(D-2)]
# S = 2R(D-2) / [p(D-2) - D(p-1)]
# S = 2R(D-2) / [pD - 2p - Dp + D]
# S = 2R(D-2) / [D - 2p]
#
# For p=4, D=12: S = 2R*10/(12-8) = 2R*10/4 = 5R
# For p=5, D=12: S = 2R*10/(12-10) = 2R*10/2 = 10R
# ===================================================================

print("\n" + "="*60)
print("Part A: Required X_{MN} from inverse problem")
print("="*60)

# Compute Ricci scalar
R_scalar = sp.S(0)
for i in range(12):
    r_val = hf.substitute(sp.cancel(Ric[i,i]))
    g_inv = sp.cancel(1/metric.matrix[i,i])
    g_inv_sub = hf.substitute(g_inv)
    R_scalar += sp.cancel(r_val * g_inv_sub)
R_scalar = sp.cancel(R_scalar)
print(f"Ricci scalar R = {R_scalar}")

for p in [4, 5]:
    S = sp.cancel(2*R_scalar*(D-2) / (D - 2*p))
    print(f"\nFor p={p}-form:")
    print(f"  Required |F|² = S = {S}")

    # Required X_{MN}/(p-1)! = 2R_{MN} + (p-1)/(D-2)*S*g_{MN}
    print(f"  Required X_{'{MN}'}/(p-1)! for each block:")
    for i, name in [(0,'t'), (2,'z1'), (3,'z2'), (4,'y0')]:
        g_val = hf.substitute(sp.cancel(metric.matrix[i,i]))
        x_val = sp.cancel(2*R_blocks[name] + R(p-1,D-2)*S*g_val)
        print(f"    X[{name}]/{sp.factorial(p-1)} = {x_val}")

# ===================================================================
# Part B: Direct 4-form construction
#
# Build F as a 4-form with free parameters for each independent component.
# The metric is block-diagonal, so FF_{MN} only gets contributions from
# components sharing the index M (or N).
#
# For the blocks to work, we need:
# - F with both t and z1 indices (for X[z1] ≠ 0) — key requirement from exp07
# - F with t and z2 indices (for D1-type contribution)
# - Transverse components for X[y0]
#
# The simplest ansatz: F has components
#   F_{t,x1,z2,yk} = f_d * H^{p_d} * H' * y_k/r   (D1-type, from dC)
#   F_{t,x1,z1,yk} = f_f * H^{p_f} * H' * y_k/r   (F1-type, from dC)
#   F_{z1,z2,yk,yl} = f_c * expr                     (cross/CS type)
#
# Actually, rather than guessing, let me work with the general structure.
# A 4-form F in 12d has C(12,4) = 495 independent components.
# But with the symmetry ansatz (all y_k equivalent under SO(8)), we reduce hugely.
#
# Key insight: The field strength components from dC have the form
#   F_{...,yk,...} = C' * H' * yk/r
# where the remaining indices are worldvolume/z-type.
# This is the ONLY radial dependence for closed forms dC with C = f(H).
#
# For a Chern-Simons term B ∧ F, the structure would be:
#   (B ∧ F)_{M1...M4} = B_{[M1 M2} F_{M3 M4]}
# where B has 2 indices and F has 2 indices (from dC_0 or similar).
# ===================================================================

print("\n" + "="*60)
print("Part B: Chern-Simons modification F = dC_D1 + lambda * B ∧ dC_0")
print("="*60)

# In IIB: F_3 = dC_2 - C_0 dB_2 (or + depending on convention)
# Uplifted to 12d: this might give F_4 = dC_3 + C_0 ∧ H_3 or similar
#
# But the Chern-Simons term is more naturally:
# F_5 = dC_4 + B_2 ∧ F_3
# Let me try building F_4 = dC_D1 + lambda * B_{z1,z2} ∧ F_2(transverse)

# More concretely: build a non-closed 4-form that has BOTH D1 and F1 type legs
# F_{t,x1,z2,yk} = f_d * H^{p_d}' * yk/r  (D1 type)
# F_{t,x1,z1,yk} = f_f * H^{p_f}' * yk/r  (F1 type)
# These are the components that dC_D1 and dC_F1 would give individually.

# Combined form: just ADD the two dC contributions
# This IS a valid (closed) 4-form since d(dC_D1 + dC_F1) = 0.
# But exp04 already showed no positive-weight combination works...

# The NEW thing is: what if F also has components with z1,z2 legs?
# F_{z1,z2,yk,yl} = f_c * expr  — these cannot come from dC of any 3-form
# with C ~ f(H), because:
# dC_{z1,z2,yl} = partial_{yk} C_{z1,z2,yl} * dyk + ...
# → F_{yk,z1,z2,yl} = dC term, which IS possible if C_{z1,z2,yl} = g(H) * yl/r

# Actually let me think about what Chern-Simons gives:
# If B_{z1,z2} = b * H^q (a 2-form along z1,z2)
# And F_2 = f * H^s * d(yk/r) (a 2-form along transverse)... this doesn't quite work.
#
# More simply: B_{z1,z2} ∧ F_{yk,yl} wouldn't be antisymmetric with 4 different indices.
# (B ∧ F)_{z1,z2,yk,yl} = B_{z1,z2} * F_{yk,yl} - ... = B_{z1,z2} * F_{yk,yl}
# (since z-indices don't overlap with y-indices, the antisymmetrization is trivial)

# Let me just BUILD the most general 4-form respecting SO(8) symmetry
# and the ansatz that everything depends only on r (through H).

# ===================================================================
# Part C: Most general SO(8)-symmetric 4-form ansatz
#
# With SO(8) symmetry, the independent 4-form structures are:
# Type 1: F_{t,x1,z2,yk} ~ yk    (1 function)  → D1-type
# Type 2: F_{t,x1,z1,yk} ~ yk    (1 function)  → F1-type
# Type 3: F_{t,z1,z2,yk} ~ yk    (1 function)  → 3-brane type
# Type 4: F_{x1,z1,z2,yk} ~ yk   (1 function)  → 3-brane type
# Type 5: F_{t,z2,yk,yl} ~ yk*yl (antisym→0) or ~ delta_{kl} (sym→not antisym)
#          Actually: F_{t,z2,yk,yl} with k≠l, ~ (yk*al - yl*ak) type?
#          Under SO(8): F_{t,z2,yk,yl} = f5 * epsilon_{kl...}? No, rank-2 antisym tensor
#          under SO(8) is the adjoint rep. For SO(8) invariance: no rank-2 antisym invariant.
#          So this type is ZERO under SO(8).
# Type 6: F_{z1,z2,yk,yl} ~ same argument, zero under SO(8)
# Type 7: F_{t,x1,yk,yl} ~ zero under SO(8) (same reason)
# Type 8: F_{z_a,yk,yl,ym} ~ zero under SO(8) (rank-3 antisym)
# Type 9: F_{yk,yl,ym,yn} ~ epsilon_{klmn...} (SO(8) 4-form invariant!)
#          Under SO(8), there IS an invariant 4-form if we consider the
#          Hodge dual structure, but for SO(8) on R^8, the unique invariant
#          4-form would need to be self-dual. Actually, there's no invariant
#          4-form under SO(8) — the 4-form rep is 70-dimensional, and the
#          SO(8) invariant subspace is 0 (the singlet doesn't appear in
#          the decomposition of ∧⁴(8)). So this is also zero.
#
# CONCLUSION: Under SO(8), only Types 1-4 survive! Each has the form
# F_{a,b,c,yk} = f(r) * yk/r  (from the SO(8) vector structure)
# with {a,b,c} being 3 of {t, x1, z1, z2}.
# ===================================================================

print("\n" + "="*60)
print("Part C: General SO(8)-symmetric 4-form (4 free functions)")
print("="*60)

# The 4 types with symbolic amplitudes
# F_{t,x1,z2,yk} = c1 * f1(H) * H' * yk/r  → D1 type
# F_{t,x1,z1,yk} = c2 * f2(H) * H' * yk/r  → F1 type
# F_{t,z1,z2,yk} = c3 * f3(H) * H' * yk/r  → mixed 1
# F_{x1,z1,z2,yk} = c4 * f4(H) * H' * yk/r → mixed 2
#
# These are the ONLY possibilities from dC with C ~ g(H).
# The H-power analysis from exp03 already found which f(H) powers match.
#
# From exp03 table:
# D1: C_{t,x1,z2}=H^{-1} → matches H-powers ✓
# F1: C_{t,x1,z1}=H^{-1/2} → matches H-powers ✓
# Type 3: C_{t,z1,z2}... checking from exp03 results
# Type 4: C_{x1,z1,z2}... checking from exp03 results

# Actually let me just check ALL four 3-form potentials and their H-power matching

print("\nChecking all 4 SO(8)-symmetric 3-form potentials:")
base_indices = [(0,1,3), (0,1,2), (0,2,3), (1,2,3)]
base_names = ['C_{t,x1,z2}', 'C_{t,x1,z1}', 'C_{t,z1,z2}', 'C_{x1,z1,z2}']

# For each, find the power p that gives matching H-powers
for base_idx, base_name in zip(base_indices, base_names):
    print(f"\n  {base_name}:")
    for p_num in range(-8, 4):
        p_val = R(p_num, 4)
        C = FormField(rank=3, dim=12)
        C[base_idx] = H**p_val
        F = exterior_derivative(C, coords)

        if len(F.nonzero_components) == 0:
            continue

        T = form_stress_energy(F, metric)

        # Check H-power matching: T[block]/R[block] should be H-independent
        all_match = True
        ratios = {}
        for i, name in [(0,'t'), (2,'z1'), (3,'z2')]:
            t_val = hf.substitute(sp.cancel(T[i,i]))
            r_val = R_blocks[name]
            if t_val == 0:
                all_match = False
                break
            ratio = sp.cancel(t_val / r_val)
            if ratio.has(hf.H):
                all_match = False
                break
            ratios[name] = ratio

        if all_match and ratios:
            print(f"    p={p_val}: H-powers match! T/R = {ratios}")

# ===================================================================
# Part D: With all 4 types identified, solve the 4-parameter system
# a₁²T₁ + a₂²T₂ + a₃²T₃ + a₄²T₄ = R
# ===================================================================

print("\n" + "="*60)
print("Part D: 4-parameter system with all SO(8)-symmetric forms")
print("="*60)

# Build forms for all matching cases found above
# From exp03: D1 (C_{t,x1,z2}=H^{-1}) and F1 (C_{t,x1,z1}=H^{-1/2}) match
# Let's find the others and build the full system

matching_forms = []

for base_idx, base_name in zip(base_indices, base_names):
    for p_num in range(-8, 4):
        p_val = R(p_num, 4)
        C = FormField(rank=3, dim=12)
        C[base_idx] = H**p_val
        F = exterior_derivative(C, coords)

        if len(F.nonzero_components) == 0:
            continue

        T = form_stress_energy(F, metric)

        all_match = True
        blocks = {}
        for i, name in [(0,'t'), (2,'z1'), (3,'z2'), (4,'y0')]:
            t_val = hf.substitute(sp.cancel(T[i,i]))
            r_val = R_blocks[name]
            if t_val == 0:
                all_match = False
                break
            ratio = sp.cancel(t_val / r_val)
            if ratio.has(hf.H):
                all_match = False
                break
            blocks[name] = ratio

        if all_match and blocks:
            matching_forms.append((base_name, p_val, blocks))

print(f"\nFound {len(matching_forms)} matching forms:")
for name, p, blocks in matching_forms:
    print(f"  {name}, p={p}: T/R = t:{blocks['t']}, z1:{blocks['z1']}, z2:{blocks['z2']}, y0:{blocks['y0']}")

# Now solve the linear system: sum_i wi * (T_i/R) = 1 for each block
# where wi = ai² ≥ 0
if len(matching_forms) >= 3:
    n_forms = len(matching_forms)
    w = [sp.Symbol(f'w{i}', positive=True) for i in range(n_forms)]

    # Build equations for t, z1, z2 blocks
    eqs = []
    block_names = ['t', 'z1', 'z2']
    for bname in block_names:
        eq = sum(w[i] * matching_forms[i][2][bname] for i in range(n_forms)) - 1
        eqs.append(eq)

    print(f"\nSolving {len(eqs)} equations in {n_forms} unknowns...")
    print("Equations:")
    for i, eq in enumerate(eqs):
        print(f"  {block_names[i]}: {eq} = 0")

    # For exactly 3 unknowns and 3 equations, solve directly
    if n_forms >= 3:
        for combo in __import__('itertools').combinations(range(n_forms), 3):
            w_sub = [w[i] for i in combo]
            eq_sub = []
            for bname in block_names:
                eq = sum(w_sub[j] * matching_forms[combo[j]][2][bname] for j in range(3)) - 1
                eq_sub.append(eq)
            sol = sp.solve(eq_sub, w_sub)
            if sol:
                names = [matching_forms[combo[j]][0] + f"(p={matching_forms[combo[j]][1]})" for j in range(3)]
                all_pos = all(sol[w_sub[j]] > 0 for j in range(3))
                print(f"\n  Combo {names}:")
                for j in range(3):
                    print(f"    w_{j} = {sol[w_sub[j]]}  {'✓' if sol[w_sub[j]] > 0 else '✗ NEGATIVE'}")

                if all_pos:
                    print("  *** ALL POSITIVE! Check y0 block:")
                    y0_check = sum(sol[w_sub[j]] * matching_forms[combo[j]][2]['y0'] for j in range(3))
                    print(f"    sum wi*T_i/R[y0] = {sp.cancel(y0_check)}")
                    if sp.cancel(y0_check - 1) == 0:
                        print("    *** PERFECT MATCH! ***")
                    else:
                        print(f"    Residual: {sp.cancel(y0_check - 1)}")

    # Also try ALL forms (overdetermined → least squares or check consistency)
    if n_forms >= 4:
        # Try 4-tuples
        for combo in __import__('itertools').combinations(range(n_forms), 4):
            w_sub = [w[i] for i in combo]
            eq_sub = []
            for bname in ['t', 'z1', 'z2', 'y0']:
                eq = sum(w_sub[j] * matching_forms[combo[j]][2][bname] for j in range(4)) - 1
                eq_sub.append(eq)
            sol = sp.solve(eq_sub, w_sub)
            if sol:
                names = [f"{matching_forms[combo[j]][0]}(p={matching_forms[combo[j]][1]})" for j in range(4)]
                all_pos = all(sol.get(w_sub[j], -1) > 0 for j in range(4))
                if all_pos:
                    print(f"\n  *** 4-form SOLUTION: {names}")
                    for j in range(4):
                        print(f"    w_{j} = {sol[w_sub[j]]}")

# ===================================================================
# Part E: Allow NEGATIVE coefficients (non-standard kinetic terms)
# In some SUGRA formulations, wrong-sign kinetic terms appear
# (ghost fields, or fields with imaginary coupling).
# Check: does dropping the wi > 0 constraint help?
# ===================================================================

print("\n" + "="*60)
print("Part E: Unrestricted weights (allow negative = ghost contributions)")
print("="*60)

if len(matching_forms) >= 3:
    # Use first 3 matching forms, solve without positivity constraint
    w_free = sp.symbols('w0 w1 w2')
    eq_free = []
    for bname in ['t', 'z1', 'z2']:
        eq = sum(w_free[j] * matching_forms[j][2][bname] for j in range(3)) - 1
        eq_free.append(eq)

    sol_free = sp.solve(eq_free, w_free)
    if sol_free:
        print(f"Solution with first 3 forms (unrestricted):")
        for j in range(3):
            name = f"{matching_forms[j][0]}(p={matching_forms[j][1]})"
            print(f"  w_{j} ({name}) = {sol_free[w_free[j]]}")

        # Check y0
        y0_val = sum(sol_free[w_free[j]] * matching_forms[j][2]['y0'] for j in range(3))
        print(f"  y0 check: {sp.cancel(y0_val)} (need 1)")

        # Check y0_aniso
        # Need to compute the anisotropic part separately
        print("\n  Checking anisotropic (y0²/r²) part separately...")
        y0s = sp.Symbol('y0', real=True)
        for j in range(3):
            base_idx = base_indices[[n for n, (bn, _, _) in enumerate(matching_forms) if n == j][0]] if j < len(base_indices) else None
            fname = matching_forms[j][0]
            pval = matching_forms[j][1]

            C = FormField(rank=3, dim=12)
            idx_map = {'C_{t,x1,z2}': (0,1,3), 'C_{t,x1,z1}': (0,1,2),
                       'C_{t,z1,z2}': (0,2,3), 'C_{x1,z1,z2}': (1,2,3)}
            if fname in idx_map:
                C[idx_map[fname]] = H**pval
                F = exterior_derivative(C, coords)
                T = form_stress_energy(F, metric)
                t_y0 = hf.substitute(sp.cancel(T[4,4]))
                r_y0 = R_blocks['y0']
                # Extract aniso coefficient (y0^2 part)
                t_aniso = sp.cancel(sp.diff(t_y0, y0s, 2) / 2)
                r_aniso = sp.cancel(sp.diff(r_y0, y0s, 2) / 2)
                if t_aniso != 0 and r_aniso != 0:
                    ratio_aniso = sp.cancel(t_aniso / r_aniso)
                    print(f"    {fname}(p={pval}): T_aniso/R_aniso = {ratio_aniso}")

# ===================================================================
# Part F: The 3-form potentials C_{t,z1,z2} and C_{x1,z1,z2}
# These were in the exp03 scan but might not have been fully explored.
# They have BOTH z1 AND z2 legs, which could be key.
# ===================================================================

print("\n" + "="*60)
print("Part F: Detailed analysis of C_{t,z1,z2} and C_{x1,z1,z2} forms")
print("="*60)

for base_idx, base_name in [((0,2,3), 'C_{t,z1,z2}'), ((1,2,3), 'C_{x1,z1,z2}')]:
    print(f"\n--- {base_name} ---")
    for p_num in range(-8, 4):
        p_val = R(p_num, 4)
        C = FormField(rank=3, dim=12)
        C[base_idx] = H**p_val
        F = exterior_derivative(C, coords)

        if len(F.nonzero_components) == 0:
            continue

        T = form_stress_energy(F, metric)

        # Quick check: is T nonzero?
        t_t = hf.substitute(sp.cancel(T[0,0]))
        if t_t == 0:
            continue

        # Check H-power matching
        all_match = True
        for i, name in [(0,'t'), (2,'z1'), (3,'z2')]:
            t_val = hf.substitute(sp.cancel(T[i,i]))
            r_val = R_blocks[name]
            if t_val == 0:
                all_match = False
                break
            ratio = sp.cancel(t_val / r_val)
            if ratio.has(hf.H):
                all_match = False
                break

        status = "H-match ✓" if all_match else "H-mismatch ✗"
        print(f"  p={p_val}: T[t]={t_t}, {status}")

        if all_match:
            for i, name in [(0,'t'), (2,'z1'), (3,'z2'), (4,'y0')]:
                t_val = hf.substitute(sp.cancel(T[i,i]))
                r_val = R_blocks[name]
                ratio = sp.cancel(t_val / r_val)
                print(f"    T/R[{name}] = {ratio}")

print("\n" + "="*60)
print("EXPERIMENT 09 COMPLETE")
print("="*60)
