"""Experiment 12: Can a dilaton absorb the y0 residual from D1+F1 (lambda=1/2)?

From exp11: D1+F1 with lambda=1/2 gives R=T exactly for t,x1,z1,z2.
The y0 residual is: D[y_k] = H'²/(4H²) * (-1 + y_k²/r²)

A dilaton Phi = a*ln(H) contributes:
  T^dil_{y_k,y_k} = (1/2) a² H'²/H² * y_k²/r²  (aniso only)
  T^dil_{wv/z} = 0

This matches the aniso part if a² = 1/2. But the iso part -H'²/(4H²) remains.

Question: Can the dilaton coupling factor e^{alpha*Phi} shift the D1+F1
form part to absorb the isotropic residual too?

Full equation with dilaton:
  R = (1/2)(dPhi)² + e^{alpha*Phi} * T^form(lambda=1/2)
  = (1/2)(dPhi)² + H^{a*alpha} * T^form(lambda=1/2)

For t,x1,z1,z2: (dPhi)² = 0, so we need H^{a*alpha} * T^form = R.
Since T^form = R already (from exp11), we need H^{a*alpha} = 1, i.e., a*alpha = 0.

If a ≠ 0 (need dilaton for y0), then alpha = 0 (no coupling).
But with alpha = 0: T = (1/2)(dPhi)² + T^form(lambda=1/2).

For y0: R = T^form + T^dil = (R - residual) + (1/2)a² H'²/H² * y0²/r²
  → residual = (1/2) a² H'²/H² * y0²/r²
  → (-1 + y0²/r²) * H'²/(4H²) = (1/2) a² H'²/H² * y0²/r²

Matching y0² coefficient: 1/(4r²) = a²/(2r²) → a² = 1/2
Matching constant term: -1/4 ≠ 0  → IMPOSSIBLE with dilaton alone

BUT: what if we allow the dilaton to modify lambda?
With dilaton in the action, the trace structure changes.
The standard formula with dilaton is:
  T_{MN} = (1/2)(dPhi_M dPhi_N) + (1/2)e^{alpha Phi}[FF/6 - lambda*|F|²*g]

And the trace gives a DIFFERENT lambda than without dilaton because
the dilaton energy-momentum contributes to the trace:
  g^{MN}T_{MN} = (1/2)|dPhi|² + (D/2)*e^{a*alpha}*(-lambda*S) + (p/2)*e^{a*alpha}*S/p!... hmm

Actually, let's just do the computation directly.
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
z1s, z2s = sp.symbols('z1 z2', real=True)
harmonic_coords = list(sp.symbols('y0:8', real=True))
coords = wv_coords + [z1s, z2s] + harmonic_coords
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

# D1 and F1 forms
C_D1 = FormField(rank=3, dim=12)
C_D1[(0,1,3)] = 1/H
F_D1 = exterior_derivative(C_D1, coords)

C_F1 = FormField(rank=3, dim=12)
C_F1[(0,1,2)] = H**R(-1,2)
F_F1 = exterior_derivative(C_F1, coords)

# Contractions
FF_D1 = form_contraction(F_D1, metric)
S_D1 = form_norm_squared(F_D1, metric)
FF_F1 = form_contraction(F_F1, metric)
S_F1 = form_norm_squared(F_F1, metric)

# ===================================================================
# Part A: D1+F1 (lambda=1/2) + dilaton Phi = a*ln(H), no coupling (alpha=0)
# ===================================================================

print("="*60)
print("Part A: D1+F1 (lambda=1/2) + dilaton (alpha=0)")
print("="*60)

cd2 = R(5,4)
cf2 = R(1,1)
lam = R(1,2)

# Dilaton: Phi = a*ln(H), dPhi_M = a * H' * y_k / (H * r) for M = y_k
# T^dil_{y_k, y_l} = (1/2) a² H'^2 / H² * y_k*y_l / r²
# T^dil_{others} = 0

a_sq = sp.Symbol('a2', positive=True)  # a²

for i, name in [(0,'t'), (2,'z1'), (3,'z2'), (4,'y0')]:
    r_val = hf.substitute(sp.cancel(Ric[i,i]))

    # Form stress-energy (modified lambda)
    ff_d1 = hf.substitute(sp.cancel(FF_D1[i,i]))
    ff_f1 = hf.substitute(sp.cancel(FF_F1[i,i]))
    s_d1 = hf.substitute(sp.cancel(S_D1))
    s_f1 = hf.substitute(sp.cancel(S_F1))
    g_val = hf.substitute(sp.cancel(metric.matrix[i,i]))

    T_form = R(1,2) * (cd2*ff_d1/6 + cf2*ff_f1/6 - lam*(cd2*s_d1 + cf2*s_f1)*g_val)
    T_form = sp.cancel(T_form)

    # Dilaton kinetic
    if i >= 4:  # transverse
        T_dil = R(1,2) * a_sq * hf.Hp**2 / hf.H**2 * sp.Symbol(name, real=True)**2 / hf.r**2
    else:
        T_dil = 0

    T_total = sp.cancel(T_form + T_dil) if T_dil != 0 else T_form
    diff = sp.cancel(r_val - T_total)

    print(f"\n  [{name}]:")
    print(f"    R      = {r_val}")
    print(f"    T_form = {T_form}")
    print(f"    T_dil  = {T_dil}")
    print(f"    R-T    = {diff}")

# For y0: R - T = (-1 + y0²/r²)H'²/(4H²) - (1/2)a² H'²/H² * y0²/r²
# = -H'²/(4H²) + (1/(4r²) - a²/(2r²))*H'²*y0²/H²
# Set coeff of y0² to 0: 1/4 - a²/2 = 0 → a² = 1/2
# Then R-T = -H'²/(4H²) (pure isotropic)

print("\n--- Setting a² = 1/2 to kill anisotropic residual ---")
a2_val = R(1,2)

for i, name in [(0,'t'), (2,'z1'), (3,'z2'), (4,'y0')]:
    r_val = hf.substitute(sp.cancel(Ric[i,i]))

    ff_d1 = hf.substitute(sp.cancel(FF_D1[i,i]))
    ff_f1 = hf.substitute(sp.cancel(FF_F1[i,i]))
    s_d1 = hf.substitute(sp.cancel(S_D1))
    s_f1 = hf.substitute(sp.cancel(S_F1))
    g_val = hf.substitute(sp.cancel(metric.matrix[i,i]))

    T_form = R(1,2) * (cd2*ff_d1/6 + cf2*ff_f1/6 - lam*(cd2*s_d1 + cf2*s_f1)*g_val)
    T_form = sp.cancel(T_form)

    if i >= 4:
        T_dil = R(1,2) * a2_val * hf.Hp**2 / hf.H**2 * sp.Symbol(name, real=True)**2 / hf.r**2
    else:
        T_dil = 0

    T_total = sp.cancel(T_form + T_dil) if T_dil != 0 else T_form
    diff = sp.cancel(r_val - T_total)
    print(f"  [{name}] R-T = {diff}")

# ===================================================================
# Part B: The remaining isotropic residual -H'²/(4H²) in y0 block
# = -(1/4) * H'²/H² per transverse direction
#
# This has the structure of |dPhi|² contribution to the trace.
# |dPhi|² = g^{MN} dPhi_M dPhi_N = (1/H^{1/4}) * a² * H'^2/H^2 * yk^2/r^2
#          summed over k... actually:
# |dPhi|² = sum_k g^{yk yk} (dPhi_{yk})^2
#          = sum_k (1/H^{1/4}) * a^2 * H'^2/H^2 * yk^2/r^2
#          = (a^2/H^{1/4}) * H'^2/H^2 * (r^2/r^2)
#          = a^2 * H'^2 / H^{9/4}
#
# For a²=1/2: |dPhi|² = H'^2/(2H^{9/4})
#
# The isotropic residual is -H'²/(4H²) per y-direction.
# With g_{yk} = H^{1/4}: residual/g = -H'²/(4H^{9/4}) per direction.
# Summed: 8 × (-H'²/(4H^{9/4})) = -2H'²/H^{9/4}
#
# Meanwhile |dPhi|² = H'²/(2H^{9/4})
# And the dilaton trace contribution to T is (1/2)|dPhi|² = H'²/(4H^{9/4})
# per direction... no, the trace is the total, not per direction.
#
# The KEY QUESTION: is the remaining isotropic residual related to a
# trace term that we're not including?
# ===================================================================

print("\n" + "="*60)
print("Part B: Structure of the remaining isotropic residual")
print("="*60)

# The residual per transverse direction: -H'²/(4H²) (iso at y0=0)
# Note g[y0] = H^{1/4}, so residual/g = -H'²/(4H^{9/4})
# This is the SAME H-power combination as S (form norm squared)
# and |dPhi|² (dilaton kinetic norm).

# S_D1 = -H'^2/H^{9/4}
# S_F1 = -H'^2/(4H^{9/4})
# |dPhi|² = H'^2/(2H^{9/4})

print("Relevant scalar densities (all proportional to H'^2/H^{9/4}):")
print(f"  S_D1 = {hf.substitute(sp.cancel(S_D1))}")
print(f"  S_F1 = {hf.substitute(sp.cancel(S_F1))}")
print(f"  |dPhi|² (a²=1/2) = H'²/(2*H^(9/4))")

# The form part of T was: T^form = (1/2)[cd²*FF_D1/6 + cf²*FF_F1/6 - lambda*S_tot*g]
# where S_tot = cd²*S_D1 + cf²*S_F1
S_tot = sp.cancel(cd2*hf.substitute(sp.cancel(S_D1)) + cf2*hf.substitute(sp.cancel(S_F1)))
print(f"  S_total (cd²=5/4,cf²=1) = {S_tot}")

# So S_total = 5/4*(-H'^2/H^{9/4}) + 1*(-H'^2/(4H^{9/4})) = -5H'^2/(4H^{9/4}) - H'^2/(4H^{9/4}) = -6H'^2/(4H^{9/4}) = -3H'^2/(2H^{9/4})
print(f"  Simplified: {sp.cancel(S_tot)}")

# The isotropic residual per direction is -H'^2/(4H^2).
# With g = H^{1/4}: -H'^2/(4H^2) * (1/H^{1/4}) = ... no, residual is already
# the component R_{y0,y0} - T_{y0,y0} with indices down.

# Let me think of this differently. The residual is:
# D[y_k] = -H'^2/(4H^2) (isotropic part)
# What PRODUCES this? It would need to be a term that's the same for all y_k.
# Candidate: additional trace subtraction term c*S*g or c*|dPhi|²*g

# Try: add a correction term delta * g[y_k] to T[y_k]
# D[y_k] = -H'^2/(4H^2) = delta * H^{1/4}
# → delta = -H'^2/(4H^{9/4})

# And delta * g[t] should be 0 (since D[t]=0):
# delta * g[t] = -H'^2/(4H^{9/4}) * (-H^{-3/4}) = H'^2/(4H^3) ≠ 0

# So adding a cosmological-type term breaks the wv/z matching.

# The ONLY way to add something only in the transverse block is through
# a form that has FF_{wv/z} exactly cancelling its trace contribution there.

# ===================================================================
# Part C: Actually, let me reconsider the formula entirely.
# What if the 12d formula has a SEPARATE lambda for each form type?
# T = sum_i [(1/2) ci² * (FF_i/norm_i - lambda_i * S_i * g)]
#
# From exp10: kappa_i[block] = R/T_i depends on lambda_i.
# We need: sum_i ci² * T_i(lambda_i) = R for all blocks including y0.
#
# With 3 forms (D1, F1, maybe 5-form) and 3 lambdas + 3 amplitudes = 6 params,
# and 5 constraints (t, z1, z2, y0_iso, y0_aniso), the system is underdetermined.
# ===================================================================

print("\n" + "="*60)
print("Part C: Solve the FULL system with separate lambdas")
print("cd²_D1 T_D1(lam_D) + cf²_F1 T_F1(lam_F) = R for ALL blocks")
print("="*60)

# Parameters: cd2, cf2, lam_D, lam_F (4 unknowns)
# Equations: t, z1, z2, y0_iso, y0_aniso (5 equations)

cd2_s, cf2_s = sp.symbols('cd cf', positive=True)  # cd = cd², cf = cf²
ld, lf = sp.symbols('ld lf')  # lambda_D, lambda_F

# Build T blocks symbolically
eqs = []
for i, name in [(0,'t'), (2,'z1'), (3,'z2')]:
    ff_d = hf.substitute(sp.cancel(FF_D1[i,i]))
    ff_f = hf.substitute(sp.cancel(FF_F1[i,i]))
    s_d = hf.substitute(sp.cancel(S_D1))
    s_f = hf.substitute(sp.cancel(S_F1))
    g = hf.substitute(sp.cancel(metric.matrix[i,i]))
    r = hf.substitute(sp.cancel(Ric[i,i]))

    T = R(1,2) * (cd2_s * (ff_d/6 - ld*s_d*g) + cf2_s * (ff_f/6 - lf*s_f*g))
    eq = sp.cancel(T - r)
    eqs.append(eq)
    print(f"  eq[{name}] = {eq}")

# y0 isotropic (y0=0)
ff_d_y0 = hf.substitute(sp.cancel(FF_D1[4,4])).subs(y0, 0)
ff_f_y0 = hf.substitute(sp.cancel(FF_F1[4,4])).subs(y0, 0)
s_d = hf.substitute(sp.cancel(S_D1))
s_f = hf.substitute(sp.cancel(S_F1))
g_y0 = hf.substitute(sp.cancel(metric.matrix[4,4]))
r_y0_iso = hf.substitute(sp.cancel(Ric[4,4])).subs(y0, 0)

T_y0_iso = R(1,2) * (cd2_s * (ff_d_y0/6 - ld*s_d*g_y0) + cf2_s * (ff_f_y0/6 - lf*s_f*g_y0))
eq_y0_iso = sp.cancel(T_y0_iso - r_y0_iso)
eqs.append(eq_y0_iso)
print(f"  eq[y0_iso] = {eq_y0_iso}")

# y0 anisotropic (coeff of y0²)
ff_d_y0_full = hf.substitute(sp.cancel(FF_D1[4,4]))
ff_f_y0_full = hf.substitute(sp.cancel(FF_F1[4,4]))
r_y0_full = hf.substitute(sp.cancel(Ric[4,4]))

ff_d_aniso = sp.cancel(sp.diff(ff_d_y0_full, y0, 2) / 2)
ff_f_aniso = sp.cancel(sp.diff(ff_f_y0_full, y0, 2) / 2)
r_aniso = sp.cancel(sp.diff(r_y0_full, y0, 2) / 2)

# Lambda terms don't contribute to y0² coefficient (S and g are y0-independent)
T_aniso = R(1,2) * (cd2_s * ff_d_aniso/6 + cf2_s * ff_f_aniso/6)
eq_aniso = sp.cancel(T_aniso - r_aniso)
eqs.append(eq_aniso)
print(f"  eq[y0_aniso] = {eq_aniso}")

# Solve: 5 equations, 4 unknowns
# First: the aniso equation constrains cd2, cf2 only
print(f"\nAniso equation: {eq_aniso} = 0")
print("  This constrains cd and cf only (no lambda dependence)")

# From aniso: express cf2 in terms of cd2
cf2_from_aniso = sp.solve(eq_aniso, cf2_s)
print(f"  cf = {cf2_from_aniso}")

if cf2_from_aniso:
    cf2_expr = cf2_from_aniso[0]
    print(f"  cf = f(cd) = {cf2_expr}")

    # Substitute into t, z1, z2, y0_iso equations (4 eqs, 3 unknowns: cd, ld, lf)
    eqs_reduced = [eq.subs(cf2_s, cf2_expr) for eq in eqs[:4]]
    for j, name in enumerate(['t', 'z1', 'z2', 'y0_iso']):
        eqs_reduced[j] = sp.cancel(eqs_reduced[j])
        print(f"  reduced eq[{name}] = {eqs_reduced[j]}")

    # Solve t, z1, z2 for cd, ld, lf
    print("\n  Solving t, z1, z2 for (cd, ld, lf)...")
    sol = sp.solve(eqs_reduced[:3], [cd2_s, ld, lf], dict=True)
    if sol:
        for s in sol:
            print(f"\n  Solution:")
            for var in [cd2_s, ld, lf]:
                v = sp.cancel(s[var]) if var in s else 'free'
                print(f"    {var} = {v}")

            cd_val = s[cd2_s]
            cf_val = sp.cancel(cf2_expr.subs(cd2_s, cd_val))
            print(f"    cf = {cf_val}")

            # Check y0_iso
            y0_res = sp.cancel(eqs_reduced[3].subs(s))
            print(f"    y0_iso residual: {y0_res}")
            if y0_res == 0:
                print(f"    *** FULL SOLUTION! ***")
    else:
        print("  No solution from t,z1,z2.")

        # Try solving pairs
        for pair, pair_name in [([0,1],'t,z1'), ([0,2],'t,z2'), ([1,2],'z1,z2')]:
            print(f"\n  Solving {pair_name} for (cd, ld) with lf free...")
            sol_p = sp.solve([eqs_reduced[pair[0]], eqs_reduced[pair[1]]], [cd2_s, ld], dict=True)
            if sol_p:
                for s in sol_p:
                    cd_v = sp.cancel(s.get(cd2_s, cd2_s))
                    ld_v = sp.cancel(s.get(ld, ld))
                    print(f"    cd = {cd_v}, ld = {ld_v}")

# ===================================================================
# Part D: What if we include dilaton AND modified lambda?
# Full system: R = T^dil + T^form(lambda)
# T^dil only contributes to transverse block (y_k² terms)
# ===================================================================

print("\n" + "="*60)
print("Part D: D1+F1+dilaton with modified lambda")
print("R = (1/2)a²(dln H)² + T^form(cd², cf², lambda)")
print("="*60)

# For wv/z blocks: same as before (T^dil = 0)
# For y0 block:
#   R[y0] = (1/2)a²*H'²/H² * y0²/r² + T^form[y0]
#
# So: aniso part: (1/2)a²/r² * H'²/H² + (1/2)(cd²*ff_d_aniso/6 + cf²*ff_f_aniso/6)
#                = r_aniso
#     iso part:   0 + T^form_iso = r_iso
#     wv/z:       0 + T^form = R

# From wv/z: cd²=5/4, cf²=1, lambda=1/2 (exp10 solution)
# From aniso: (1/2)a²/r² + (1/2)(5/4*ff_d_aniso + 1*ff_f_aniso)/6 = r_aniso

# Let me compute the form's aniso contribution explicitly
form_aniso = sp.cancel(R(1,2) * (cd2*ff_d_aniso/6 + cf2*ff_f_aniso/6))
print(f"  Form aniso contribution: {form_aniso}")
print(f"  R aniso: {r_aniso}")
aniso_diff = sp.cancel(r_aniso - form_aniso)
print(f"  Needed from dilaton aniso: {aniso_diff}")
# (1/2) a² H'²/(H² r²) = aniso_diff
a_sq_needed = sp.cancel(aniso_diff * 2 * hf.H**2 * hf.r**2 / hf.Hp**2)
print(f"  → a² = {a_sq_needed}")

# From iso: T^form_iso = R_iso (no dilaton contribution to iso part)
form_iso = sp.cancel(R(1,2) * (cd2*ff_d_y0/6 + cf2*ff_f_y0/6 - lam*(cd2*s_d + cf2*s_f)*g_y0))
print(f"\n  Form iso contribution: {form_iso}")
print(f"  R iso: {r_y0_iso}")
iso_diff = sp.cancel(r_y0_iso - form_iso)
print(f"  Iso residual (must be 0): {iso_diff}")

if iso_diff != 0:
    print(f"\n  *** Iso residual is NOT zero: {iso_diff}")
    print("  The dilaton cannot contribute to the isotropic part.")
    print("  Need to adjust lambda to fix BOTH iso and wv/z blocks simultaneously.")

    # Try: allow lambda to vary and add dilaton
    # Equations:
    # t: cd²*(FF_D1[t]/6 - lam*S_D1*g[t]) + cf²*(FF_F1[t]/6 - lam*S_F1*g[t]) = 2R[t]
    # (same for z1, z2)
    # y0_iso: cd²*(ff_d_iso/6 - lam*S_D1*g_y0) + cf²*(ff_f_iso/6 - lam*S_F1*g_y0) = 2*r_y0_iso
    # y0_aniso: cd²*ff_d_aniso/6 + cf²*ff_f_aniso/6 + a²*H'^2/(H²r²) = 2*r_aniso
    #
    # 5 eqs (t, z1, z2, y0_iso, y0_aniso), 4 unknowns (cd², cf², lambda, a²)

    print("\n  Solving full 5-eq/4-unknown system with dilaton...")

    a2_s = sp.Symbol('asq', positive=True)

    full_eqs = []
    for i, name in [(0,'t'), (2,'z1'), (3,'z2')]:
        ff_d = hf.substitute(sp.cancel(FF_D1[i,i]))
        ff_f = hf.substitute(sp.cancel(FF_F1[i,i]))
        s_d = hf.substitute(sp.cancel(S_D1))
        s_f = hf.substitute(sp.cancel(S_F1))
        g = hf.substitute(sp.cancel(metric.matrix[i,i]))
        r = hf.substitute(sp.cancel(Ric[i,i]))

        eq = cd2_s*(ff_d/6 - ld*s_d*g) + cf2_s*(ff_f/6 - lf*s_f*g) - 2*r
        full_eqs.append(sp.cancel(eq))

    # y0 iso
    eq_iso = cd2_s*(ff_d_y0/6 - ld*s_d*g_y0) + cf2_s*(ff_f_y0/6 - lf*s_f*g_y0) - 2*r_y0_iso
    full_eqs.append(sp.cancel(eq_iso))

    # y0 aniso (with dilaton)
    eq_aniso_dil = cd2_s*ff_d_aniso/6 + cf2_s*ff_f_aniso/6 + a2_s*hf.Hp**2/(hf.H**2*hf.r**2) - 2*r_aniso
    full_eqs.append(sp.cancel(eq_aniso_dil))

    # Simplify: use lam_D = lam_F = lam (shared lambda, different from Part C)
    lam_s = sp.Symbol('lam')
    full_eqs_shared = [eq.subs([(ld, lam_s), (lf, lam_s)]) for eq in full_eqs]

    print("  With shared lambda (lam_D = lam_F = lam):")
    print("  5 equations, 4 unknowns: cd², cf², lam, a²")

    # Solve aniso for a² in terms of cd², cf²
    a2_from_aniso = sp.solve(full_eqs_shared[4], a2_s)
    if a2_from_aniso:
        a2_expr = a2_from_aniso[0]
        print(f"  From aniso: a² = {a2_expr}")

        # Now 4 equations (t, z1, z2, y0_iso) in 3 unknowns (cd², cf², lam)
        remaining = full_eqs_shared[:4]
        print("\n  Solving t, z1, z2, y0_iso for (cd², cf², lam)...")

        # Try solving t, z1, z2 first
        sol3 = sp.solve(remaining[:3], [cd2_s, cf2_s, lam_s], dict=True)
        if sol3:
            for s in sol3:
                print(f"\n  From t,z1,z2:")
                for var in [cd2_s, cf2_s, lam_s]:
                    print(f"    {var} = {sp.cancel(s[var])}")

                # Check y0_iso
                y0_check = sp.cancel(remaining[3].subs(s))
                print(f"    y0_iso check: {y0_check}")

                if y0_check == 0:
                    a2_v = sp.cancel(a2_expr.subs(s))
                    print(f"    a² = {a2_v}")
                    if a2_v > 0:
                        print(f"    *** FULL SOLUTION WITH DILATON! ***")
                        print(f"    cd² = {sp.cancel(s[cd2_s])}")
                        print(f"    cf² = {sp.cancel(s[cf2_s])}")
                        print(f"    lambda = {sp.cancel(s[lam_s])}")
                        print(f"    a² = {a2_v} → a = {sp.sqrt(a2_v)}")

        # Also try solving t, z1, y0_iso and checking z2
        sol_alt = sp.solve([remaining[0], remaining[1], remaining[3]], [cd2_s, cf2_s, lam_s], dict=True)
        if sol_alt:
            for s in sol_alt:
                z2_check = sp.cancel(remaining[2].subs(s))
                if z2_check == 0:
                    print(f"\n  From t,z1,y0_iso (checking z2):")
                    for var in [cd2_s, cf2_s, lam_s]:
                        print(f"    {var} = {sp.cancel(s[var])}")
                    print(f"    z2 check: ✓")
                    a2_v = sp.cancel(a2_expr.subs(s))
                    print(f"    a² = {a2_v}")

print("\n" + "="*60)
print("EXPERIMENT 12 COMPLETE")
print("="*60)
