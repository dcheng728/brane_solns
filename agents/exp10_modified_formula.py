"""Experiment 10: Modified stress-energy formula.

Direction #3 from notes: The standard T formula might not apply in 12d.

Key observation from exp05: D1 residual D = R - T has pattern (1,4,4,-1).
What if the trace subtraction coefficient differs from (p-1)/(D-2)?

General formula: T_{MN} = (1/2)[FF_{MN}/(p-1)! - lambda * |F²| * g_{MN}]
Standard: lambda = (p-1)/(D-2) = 3/10 for p=4, D=12

Question 1: Is there a lambda such that D1 ALONE satisfies R = T?
Question 2: Is D/g_{MN} (residual per metric component) constant across blocks?
Question 3: Two-parameter formula T = (1/2)[alpha*FF - beta*|F²|*g]?
Question 4: What about INCLUDING multiple forms with the modified formula?
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

# --- D1 form ---
C_D1 = FormField(rank=3, dim=12)
C_D1[(0,1,3)] = 1/H
F_D1 = exterior_derivative(C_D1, coords)

# Compute D1's FF contraction and norm squared
FF_D1 = form_contraction(F_D1, metric)
Fsq_D1 = form_norm_squared(F_D1, metric)

# Get simplified block values
R_blk = {}
FF_blk = {}
g_blk = {}
for i, name in [(0,'t'), (2,'z1'), (3,'z2'), (4,'y0')]:
    R_blk[name] = hf.substitute(sp.cancel(Ric[i,i]))
    FF_blk[name] = hf.substitute(sp.cancel(FF_D1[i,i]))
    g_blk[name] = hf.substitute(sp.cancel(metric.matrix[i,i]))

S = hf.substitute(sp.cancel(Fsq_D1))

print("D1 form components:")
print(f"  |F|² = S = {S}")
for name in ['t', 'z1', 'z2', 'y0']:
    print(f"  FF[{name}] = {FF_blk[name]}")
    print(f"  R[{name}]  = {R_blk[name]}")
    print(f"  g[{name}]  = {g_blk[name]}")
    print()

# ===================================================================
# Part A: Residual analysis D/g_{MN}
# T_{MN}(std) = (1/2)[FF_{MN}/3! - (3/10)*S*g_{MN}]
# D_{MN} = R_{MN} - T_{MN}(std)
# ===================================================================

print("="*60)
print("Part A: D1 residual D = R - T per metric component")
print("="*60)

lam_std = R(3,10)  # (p-1)/(D-2) for p=4, D=12

for name in ['t', 'z1', 'z2', 'y0']:
    T_std = R(1,2) * (FF_blk[name]/6 - lam_std * S * g_blk[name])
    D = sp.cancel(R_blk[name] - T_std)
    Dg = sp.cancel(D / g_blk[name])
    print(f"  D[{name}] = {D}")
    print(f"  D[{name}]/g[{name}] = {Dg}")
    print()

# ===================================================================
# Part B: Solve for lambda such that D1 + modified formula matches R
# R = (1/2)[FF/6 - lambda*S*g]
# → 2R = FF/6 - lambda*S*g
# → lambda = (FF/6 - 2R) / (S*g)
# ===================================================================

print("="*60)
print("Part B: Required lambda per block")
print("="*60)

for name in ['t', 'z1', 'z2']:
    lam_req = sp.cancel((FF_blk[name]/6 - 2*R_blk[name]) / (S * g_blk[name]))
    print(f"  lambda[{name}] = {lam_req}")

# y0: split into isotropic and anisotropic
r_y0_iso = R_blk['y0'].subs(y0, 0)
ff_y0_iso = FF_blk['y0'].subs(y0, 0)
lam_y0_iso = sp.cancel((ff_y0_iso/6 - 2*r_y0_iso) / (S * g_blk['y0']))
print(f"  lambda[y0_iso] = {lam_y0_iso}")

# Anisotropic: coefficient of y0^2
r_y0_aniso = sp.cancel(sp.diff(R_blk['y0'], y0, 2) / 2)
ff_y0_aniso = sp.cancel(sp.diff(FF_blk['y0'], y0, 2) / 2)
# S and g don't depend on y0, so:
# d²/dy0² [lambda*S*g] / 2 = 0 (since S*g is y0-independent)
# d²/dy0² [FF/6] / 2 = ff_y0_aniso/6
# d²/dy0² [2R] / 2 = r_y0_aniso
# So the aniso equation: ff_y0_aniso/6 = 2*r_y0_aniso (no lambda dependence!)
aniso_check = sp.cancel(ff_y0_aniso/6 - 2*r_y0_aniso)
print(f"\n  Anisotropic check: FF_aniso/6 - 2R_aniso = {aniso_check}")
if aniso_check == 0:
    print("  *** ANISOTROPIC PART MATCHES EXACTLY (lambda-independent)! ***")
else:
    print("  Anisotropic part does NOT match.")

# ===================================================================
# Part C: Two-parameter formula T = (1/2)[alpha*FF/6 - beta*S*g]
# Solve: alpha*FF/6 - beta*S*g = 2R for each block
# ===================================================================

print("\n" + "="*60)
print("Part C: Two-parameter formula alpha*FF/6 - beta*S*g = 2R")
print("="*60)

alpha, beta = sp.symbols('alpha beta')

# Build equations for t, z1, z2
eqs = []
for name in ['t', 'z1', 'z2']:
    eq = sp.Eq(alpha * FF_blk[name]/6 - beta * S * g_blk[name], 2*R_blk[name])
    eqs.append(eq)

# Solve t and z1 for alpha, beta
sol = sp.solve([eqs[0], eqs[1]], [alpha, beta])
if sol:
    a_val = sp.cancel(sol[alpha])
    b_val = sp.cancel(sol[beta])
    print(f"  From t,z1: alpha = {a_val}, beta = {b_val}")

    # Check z2
    z2_check = sp.cancel(a_val * FF_blk['z2']/6 - b_val * S * g_blk['z2'] - 2*R_blk['z2'])
    print(f"  z2 check: {z2_check}  {'✓' if z2_check == 0 else '✗'}")

    # Check y0 (isotropic)
    y0_check = sp.cancel(a_val * ff_y0_iso/6 - b_val * S * g_blk['y0'] - 2*r_y0_iso)
    print(f"  y0_iso check: {y0_check}  {'✓' if y0_check == 0 else '✗'}")

# Solve t and z2
sol2 = sp.solve([eqs[0], eqs[2]], [alpha, beta])
if sol2:
    a_val2 = sp.cancel(sol2[alpha])
    b_val2 = sp.cancel(sol2[beta])
    print(f"\n  From t,z2: alpha = {a_val2}, beta = {b_val2}")

    # Check z1
    z1_check = sp.cancel(a_val2 * FF_blk['z1']/6 - b_val2 * S * g_blk['z1'] - 2*R_blk['z1'])
    print(f"  z1 check: {z1_check}  {'✓' if z1_check == 0 else '✗'}")

# Solve z1 and z2
sol3 = sp.solve([eqs[1], eqs[2]], [alpha, beta])
if sol3:
    a_val3 = sp.cancel(sol3[alpha])
    b_val3 = sp.cancel(sol3[beta])
    print(f"\n  From z1,z2: alpha = {a_val3}, beta = {b_val3}")

    # Check t
    t_check = sp.cancel(a_val3 * FF_blk['t']/6 - b_val3 * S * g_blk['t'] - 2*R_blk['t'])
    print(f"  t check: {t_check}  {'✓' if t_check == 0 else '✗'}")

# ===================================================================
# Part D: Three-parameter: T = (1/2)[alpha*FF_D1/6 + gamma*FF_F1/6 - beta*S_total*g]
# This allows two forms with a SHARED trace term
# ===================================================================

print("\n" + "="*60)
print("Part D: D1+F1 with shared trace: alpha*FF_D1 + gamma*FF_F1 - beta*S*g = 2R")
print("="*60)

C_F1 = FormField(rank=3, dim=12)
C_F1[(0,1,2)] = H**R(-1,2)
F_F1 = exterior_derivative(C_F1, coords)
FF_F1 = form_contraction(F_F1, metric)
Fsq_F1 = form_norm_squared(F_F1, metric)

FF_F1_blk = {}
for i, name in [(0,'t'), (2,'z1'), (3,'z2'), (4,'y0')]:
    FF_F1_blk[name] = hf.substitute(sp.cancel(FF_F1[i,i]))

S_F1 = hf.substitute(sp.cancel(Fsq_F1))

print(f"F1 |F|² = {S_F1}")
for name in ['t', 'z1', 'z2']:
    print(f"  FF_F1[{name}] = {FF_F1_blk[name]}")

gamma = sp.Symbol('gamma')

# Total S depends on amplitudes: S_total = alpha * S_D1 + gamma * S_F1
# But actually, alpha and gamma multiply FF, so the trace relation is different.
# Let me use the more physical parametrization:
# T = (1/2) [cd² FF_D1/6 + cf² FF_F1/6 - lambda*(cd² S_D1 + cf² S_F1)*g]
# i.e., same lambda for both forms, different amplitudes

cd2, cf2 = sp.symbols('cd2 cf2')
lam = sp.Symbol('lambda')

eqs_3p = []
for name in ['t', 'z1', 'z2']:
    eq = cd2 * FF_blk[name]/6 + cf2 * FF_F1_blk[name]/6 - lam*(cd2*S + cf2*S_F1)*g_blk[name] - 2*R_blk[name]
    eqs_3p.append(sp.cancel(eq))

# 3 equations, 3 unknowns (cd2, cf2, lam)
print("\nSolving 3-parameter system...")
sol_3p = sp.solve(eqs_3p, [cd2, cf2, lam], dict=True)
if sol_3p:
    for s in sol_3p:
        print(f"\n  Solution:")
        for var in [cd2, cf2, lam]:
            val = sp.cancel(s[var]) if var in s else 'free'
            print(f"    {var} = {val}")
else:
    print("  No solution found. Trying parametric approach...")
    # Fix lam and solve for cd2, cf2
    for lam_val in [R(1,4), R(3,10), R(1,3), R(2,5), R(1,2), R(3,5)]:
        eqs_2p = []
        for name in ['t', 'z1', 'z2']:
            eq = cd2 * (FF_blk[name]/6 - lam_val*S*g_blk[name]) + cf2 * (FF_F1_blk[name]/6 - lam_val*S_F1*g_blk[name]) - 2*R_blk[name]
            eqs_2p.append(sp.cancel(eq))

        # Solve t and z1 for cd2, cf2
        try:
            sol_2p = sp.solve([eqs_2p[0], eqs_2p[1]], [cd2, cf2])
            if sol_2p and cd2 in sol_2p:
                cd2_v = sp.cancel(sol_2p[cd2])
                cf2_v = sp.cancel(sol_2p[cf2])
                z2_res = sp.cancel(eqs_2p[2].subs([(cd2, cd2_v), (cf2, cf2_v)]))
                match = "✓" if z2_res == 0 else f"✗ (residual={z2_res})"
                print(f"  lambda={lam_val}: cd²={cd2_v}, cf²={cf2_v}, z2 check: {match}")
                if z2_res == 0:
                    # Check y0
                    y0_eq = cd2_v * (ff_y0_iso/6 - lam_val*S*g_blk['y0']) + cf2_v * (FF_F1_blk['y0'].subs(y0,0)/6 - lam_val*S_F1*g_blk['y0']) - 2*r_y0_iso
                    y0_res = sp.cancel(y0_eq)
                    print(f"         y0_iso check: {y0_res}")
        except:
            pass

# ===================================================================
# Part E: What if D=12 theory has F^4 (not F^2) in the action?
# Some higher-dimensional theories have higher-order form terms.
# T_{MN} = FF_{MN}/6 * f(|F²|) - g(|F²|) * g_{MN}
# Too many free functions. Let me try the simplest modification:
# T_{MN} = (1/2n) [FF_{MN}/6 - ((p-1)/(D-2)) |F²| g_{MN}]
# where n is a normalization factor (like 2*kappa² or similar).
# This is just an overall rescaling — same as changing the form amplitude.
# So it can't help with the RATIO mismatch between blocks.
# ===================================================================

# ===================================================================
# Part F: What if the Einstein equation is NOT R_{MN} = T_{MN} but
# G_{MN} = T_{MN} (Einstein form)?
# G_{MN} = R_{MN} - (1/2) R g_{MN}
# Then: R_{MN} = T_{MN} + (1/(D-2)) T g_{MN}  (trace-reversed)
# where T = g^{MN} T_{MN}
#
# In the "Ricci form": R_{MN} = T_{MN} - (1/(D-2)) T_trace g_{MN}
# Wait, standard SUGRA uses: R_{MN} = T_{MN} (Ricci form)
# which is the trace-reversed version of G_{MN} = T_{MN}^{Einstein}
#
# The Ricci form IS what form_stress_energy computes.
# But what if the 12d theory uses the Einstein form?
# G_{MN} = T^E_{MN} = (1/2)[FF/(p-1)! - (1/2)|F²|g]
# Then R_{MN} = T^E_{MN} + (1/(D-2)) T^E g_{MN}
# = (1/2)[FF/6 - (1/2)Sg] + (1/(D-2)) * (1/2)[p*S - D*S/2] g
# = (1/2)FF/6 - S/4 g + S/(2(D-2))[p - D/2] g
# = (1/2)FF/6 + [-1/4 + (p - D/2)/(2(D-2))] S g
# = (1/2)FF/6 + [-1/4 + (4-6)/(20)] S g
# = (1/2)FF/6 + [-1/4 - 1/10] S g
# = (1/2)FF/6 - [7/20] S g
# = (1/2)[FF/6 - (7/10)Sg]
# i.e., lambda = 7/10 instead of 3/10
# ===================================================================

print("\n" + "="*60)
print("Part F: Einstein form G=T → R = T + trace term")
print("What if the equation is G_{MN} = T^E, not R_{MN} = T^Ricci?")
print("="*60)

# If G = T^E: with T^E = (1/2)[FF/6 - (1/2)|F²|g]
# Then R = (1/2)FF/6 + S*g * [(p-D/2)/(2(D-2)) - 1/4]
# For p=4, D=12: (4-6)/(20) - 1/4 = -1/10 - 1/4 = -7/20
# R = (1/2)FF/6 - (7/20) S g
# lambda_eff = 7/10  (after multiplying by 2)

# Wait, let me be more careful. The standard form kinetic term in the action is
# S = -1/(2p!) ∫ |F|² √g d^Dx
# The variation gives: T^E_{MN} = (1/2p!) [p F_{MA...}F_N^{A...} - (1/2)|F²|g_{MN}]
#                    = (1/2) [FF_{MN}/(p-1)! - (1/(2p))|F²| p g_{MN}]  ... hmm
# Actually: T^E_{MN} = (1/(2(p-1)!)) [FF_{MN} - (1/2) |F²| g_{MN} ... no

# Standard: from -1/(2p!) |F|² with |F|² = F_{M1..Mp}F^{M1..Mp}
# T^E_{MN} = (1/(p-1)!) F_{MA2..Ap}F_N^{A2..Ap} - (1/(2p!)) |F|² g_{MN}
# But our FF_{MN} includes the 1/(p-1)! normalization? No, form_contraction gives
# F_{M,A2,...Ap} F_N^{A2,...Ap} without the 1/(p-1)! factor.

# Let me just check: what lambda_eff values arise from different formulations?

print("\nChecking various effective lambda values:")
for lam_test_num, lam_test_den, label in [
    (3, 10, "standard Ricci form"),
    (7, 20, "Einstein form variant"),
    (1, 4, "1/4"),
    (1, 5, "1/5"),
    (1, 3, "1/3"),
    (2, 5, "2/5"),
    (7, 10, "7/10"),
    (1, 2, "1/2"),
]:
    lam_test = R(lam_test_num, lam_test_den)
    # For each lambda, compute what kappa (= R/T) would be per block
    kappas = {}
    for name in ['t', 'z1', 'z2', 'y0_iso']:
        if name == 'y0_iso':
            ff_val = ff_y0_iso
            r_val = r_y0_iso
            g_val = g_blk['y0']
        else:
            ff_val = FF_blk[name]
            r_val = R_blk[name]
            g_val = g_blk[name]
        T_val = R(1,2) * (ff_val/6 - lam_test * S * g_val)
        if T_val != 0:
            kappas[name] = sp.cancel(r_val / T_val)
        else:
            kappas[name] = 'inf'

    all_same = len(set(str(v) for v in kappas.values())) == 1
    print(f"\n  lambda={lam_test} ({label}):")
    for name, k in kappas.items():
        print(f"    kappa[{name}] = {k}")
    if all_same:
        print(f"    *** ALL EQUAL → SOLUTION with kappa = {list(kappas.values())[0]} ***")

# ===================================================================
# Part G: Solve for lambda that makes kappa[t] = kappa[z1] = kappa[z2]
# kappa[block] = R[block] / T[block](lambda)
# T[block](lambda) = (1/2)[FF[block]/6 - lambda*S*g[block]]
# Need: R_t / [(1/2)(FF_t/6 - lambda*S*g_t)] = R_z1 / [(1/2)(FF_z1/6 - lambda*S*g_z1)]
# ===================================================================

print("\n" + "="*60)
print("Part G: Solve for lambda analytically")
print("="*60)

lam_sym = sp.Symbol('lambda')

# kappa[t] = kappa[z1]
T_t = R(1,2) * (FF_blk['t']/6 - lam_sym * S * g_blk['t'])
T_z1 = R(1,2) * (FF_blk['z1']/6 - lam_sym * S * g_blk['z1'])
T_z2 = R(1,2) * (FF_blk['z2']/6 - lam_sym * S * g_blk['z2'])

# Cross-multiply: R_t * T_z1 = R_z1 * T_t
eq1 = sp.cancel(R_blk['t'] * T_z1 - R_blk['z1'] * T_t)
sol_lam_1 = sp.solve(eq1, lam_sym)
print(f"  kappa[t]=kappa[z1] → lambda = {sol_lam_1}")

# kappa[t] = kappa[z2]
eq2 = sp.cancel(R_blk['t'] * T_z2 - R_blk['z2'] * T_t)
sol_lam_2 = sp.solve(eq2, lam_sym)
print(f"  kappa[t]=kappa[z2] → lambda = {sol_lam_2}")

# kappa[z1] = kappa[z2]
eq3 = sp.cancel(R_blk['z1'] * T_z2 - R_blk['z2'] * T_z1)
sol_lam_3 = sp.solve(eq3, lam_sym)
print(f"  kappa[z1]=kappa[z2] → lambda = {sol_lam_3}")

# If all three give the same lambda, we have a solution!
if sol_lam_1 and sol_lam_2 and sol_lam_3:
    if sp.cancel(sol_lam_1[0] - sol_lam_2[0]) == 0:
        print(f"\n  *** ALL CONSISTENT: lambda = {sol_lam_1[0]} ***")
        lam_soln = sol_lam_1[0]
        T_t_soln = R(1,2) * (FF_blk['t']/6 - lam_soln * S * g_blk['t'])
        kappa_soln = sp.cancel(R_blk['t'] / T_t_soln)
        print(f"  kappa = {kappa_soln}")

        # Check y0
        T_y0_iso = R(1,2) * (ff_y0_iso/6 - lam_soln * S * g_blk['y0'])
        kappa_y0 = sp.cancel(r_y0_iso / T_y0_iso)
        print(f"  kappa[y0_iso] = {kappa_y0}")

        # Check aniso
        T_y0_aniso = R(1,2) * sp.cancel(sp.diff(FF_blk['y0'], y0, 2)/2) / 6
        kappa_aniso = sp.cancel(r_y0_aniso / T_y0_aniso) if T_y0_aniso != 0 else 'inf'
        print(f"  kappa[y0_aniso] = {kappa_aniso}")
    else:
        print(f"\n  INCONSISTENT: {sol_lam_1[0]} vs {sol_lam_2[0]} vs {sol_lam_3[0]}")
        print("  No single lambda works.")

        # What is the "best" lambda? Compute the one that minimizes residual
        # Or: what if the form is different from standard? Check the pattern.
        print("\n  Examining the structure more carefully:")
        print(f"  FF_D1[t]/g[t]   = {sp.cancel(FF_blk['t']/g_blk['t'])}")
        print(f"  FF_D1[z1]/g[z1] = {sp.cancel(FF_blk['z1']/g_blk['z1'])}")
        print(f"  FF_D1[z2]/g[z2] = {sp.cancel(FF_blk['z2']/g_blk['z2'])}")
        print(f"  R[t]/g[t]       = {sp.cancel(R_blk['t']/g_blk['t'])}")
        print(f"  R[z1]/g[z1]     = {sp.cancel(R_blk['z1']/g_blk['z1'])}")
        print(f"  R[z2]/g[z2]     = {sp.cancel(R_blk['z2']/g_blk['z2'])}")

# ===================================================================
# Part H: 5-form analysis with modified formula
# ===================================================================

print("\n" + "="*60)
print("Part H: Same analysis for 5-form A_{t,x1,z1,z2}=H^{-3/4}")
print("="*60)

A5 = FormField(rank=4, dim=12)
A5[(0,1,2,3)] = H**R(-3,4)
F5 = exterior_derivative(A5, coords)

FF_5 = form_contraction(F5, metric)
Fsq_5 = form_norm_squared(F5, metric)

FF_5_blk = {}
for i, name in [(0,'t'), (2,'z1'), (3,'z2'), (4,'y0')]:
    FF_5_blk[name] = hf.substitute(sp.cancel(FF_5[i,i]))
S_5 = hf.substitute(sp.cancel(Fsq_5))

print(f"  5-form |F|² = {S_5}")
for name in ['t', 'z1', 'z2']:
    print(f"  FF_5[{name}] = {FF_5_blk[name]}")

# For 5-form: standard lambda = (p-1)/(D-2) = 4/10 = 2/5
lam_5_std = R(2,5)

# Solve for lambda: kappa[t] = kappa[z1]
T5_t = R(1,2) * (FF_5_blk['t']/24 - lam_sym * S_5 * g_blk['t'])
T5_z1 = R(1,2) * (FF_5_blk['z1']/24 - lam_sym * S_5 * g_blk['z1'])
T5_z2 = R(1,2) * (FF_5_blk['z2']/24 - lam_sym * S_5 * g_blk['z2'])

eq5_1 = sp.cancel(R_blk['t'] * T5_z1 - R_blk['z1'] * T5_t)
sol5_1 = sp.solve(eq5_1, lam_sym)
eq5_2 = sp.cancel(R_blk['t'] * T5_z2 - R_blk['z2'] * T5_t)
sol5_2 = sp.solve(eq5_2, lam_sym)
eq5_3 = sp.cancel(R_blk['z1'] * T5_z2 - R_blk['z2'] * T5_z1)
sol5_3 = sp.solve(eq5_3, lam_sym)

print(f"\n  kappa[t]=kappa[z1] → lambda = {sol5_1}")
print(f"  kappa[t]=kappa[z2] → lambda = {sol5_2}")
print(f"  kappa[z1]=kappa[z2] → lambda = {sol5_3}")

if sol5_1 and sol5_2:
    if sp.cancel(sol5_1[0] - sol5_2[0]) == 0:
        print(f"  *** CONSISTENT: lambda = {sol5_1[0]} ***")
    else:
        print(f"  INCONSISTENT: {sol5_1[0]} vs {sol5_2[0]}")

# ===================================================================
# Part I: 3-form (B_{t,x1}=H^{-3/4}) with modified formula
# ===================================================================

print("\n" + "="*60)
print("Part I: 3-form B_{t,x1}=H^{-3/4}")
print("="*60)

B2 = FormField(rank=2, dim=12)
B2[(0,1)] = H**R(-3,4)
F3 = exterior_derivative(B2, coords)

FF_3 = form_contraction(F3, metric)
Fsq_3 = form_norm_squared(F3, metric)

FF_3_blk = {}
for i, name in [(0,'t'), (2,'z1'), (3,'z2'), (4,'y0')]:
    FF_3_blk[name] = hf.substitute(sp.cancel(FF_3[i,i]))
S_3 = hf.substitute(sp.cancel(Fsq_3))

print(f"  3-form |F|² = {S_3}")
for name in ['t', 'z1', 'z2']:
    print(f"  FF_3[{name}] = {FF_3_blk[name]}")

# For 3-form: (p-1)! = 2
T3_t = R(1,2) * (FF_3_blk['t']/2 - lam_sym * S_3 * g_blk['t'])
T3_z1 = R(1,2) * (FF_3_blk['z1']/2 - lam_sym * S_3 * g_blk['z1'])
T3_z2 = R(1,2) * (FF_3_blk['z2']/2 - lam_sym * S_3 * g_blk['z2'])

eq3_1 = sp.cancel(R_blk['t'] * T3_z1 - R_blk['z1'] * T3_t)
sol3_1 = sp.solve(eq3_1, lam_sym)
eq3_2 = sp.cancel(R_blk['t'] * T3_z2 - R_blk['z2'] * T3_t)
sol3_2 = sp.solve(eq3_2, lam_sym)
eq3_3 = sp.cancel(R_blk['z1'] * T3_z2 - R_blk['z2'] * T3_z1)
sol3_3 = sp.solve(eq3_3, lam_sym)

print(f"\n  kappa[t]=kappa[z1] → lambda = {sol3_1}")
print(f"  kappa[t]=kappa[z2] → lambda = {sol3_2}")
print(f"  kappa[z1]=kappa[z2] → lambda = {sol3_3}")

if sol3_1 and sol3_2:
    if sp.cancel(sol3_1[0] - sol3_2[0]) == 0:
        print(f"  *** CONSISTENT: lambda = {sol3_1[0]} ***")
        lam_3_soln = sol3_1[0]
        T3_t_soln = R(1,2) * (FF_3_blk['t']/2 - lam_3_soln * S_3 * g_blk['t'])
        kappa_3 = sp.cancel(R_blk['t'] / T3_t_soln)
        print(f"  kappa = {kappa_3}")
    else:
        print(f"  INCONSISTENT: {sol3_1[0]} vs {sol3_2[0]}")

# ===================================================================
# Part J: Combined D1 + B with modified lambda
# T = c_d² * T_D1(lambda) + c_b² * T_B(lambda)
# where lambda is a single shared parameter
# ===================================================================

print("\n" + "="*60)
print("Part J: D1 + B_{t,x1} combination with shared modified lambda")
print("="*60)

# For a given lambda, T_form[block] = (1/2)[FF/norm - lambda*S*g]
# Total: cd² * T_D1(lambda) + cb² * T_B(lambda) = R

# Parametric: for each lambda, solve for cd², cb²
lam_range = [R(i, 20) for i in range(-5, 15)]
for lam_v in lam_range:
    TD1 = {}
    TB = {}
    for name in ['t', 'z1', 'z2']:
        TD1[name] = R(1,2) * sp.cancel(FF_blk[name]/6 - lam_v * S * g_blk[name])
        TB[name] = R(1,2) * sp.cancel(FF_3_blk[name]/2 - lam_v * S_3 * g_blk[name])

    # Solve from t, z1
    cd2_sym, cb2_sym = sp.symbols('cd2 cb2')
    eq_t = cd2_sym * TD1['t'] + cb2_sym * TB['t'] - R_blk['t']
    eq_z1 = cd2_sym * TD1['z1'] + cb2_sym * TB['z1'] - R_blk['z1']
    try:
        s = sp.solve([eq_t, eq_z1], [cd2_sym, cb2_sym])
        if s and cd2_sym in s:
            cd2_v = sp.cancel(s[cd2_sym])
            cb2_v = sp.cancel(s[cb2_sym])
            # Check z2
            z2_res = sp.cancel(cd2_v * TD1['z2'] + cb2_v * TB['z2'] - R_blk['z2'])
            if z2_res == 0 and cd2_v > 0 and cb2_v > 0:
                print(f"  lambda={lam_v}: cd²={cd2_v}, cb²={cb2_v} — z2 ✓, ALL POSITIVE!")
                # Check y0
                TD1_y0 = R(1,2) * sp.cancel(ff_y0_iso/6 - lam_v * S * g_blk['y0'])
                TB_y0 = R(1,2) * sp.cancel(FF_3_blk['y0'].subs(y0,0)/2 - lam_v * S_3 * g_blk['y0'])
                y0_res = sp.cancel(cd2_v * TD1_y0 + cb2_v * TB_y0 - r_y0_iso)
                print(f"         y0_iso check: {y0_res}")
    except:
        pass

print("\n" + "="*60)
print("SUMMARY")
print("="*60)
