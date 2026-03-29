"""Experiment 07: Try a D1 4-form combined with a symmetrized 5-form.

The 5-form is built from the sum of A_{x1,z1,z2,yk} = c*H^p for all k.
This should restore SO(8) symmetry in the transverse space.

Also try: A single 4-form with BOTH D1 and F1 components (varying amplitudes),
combined with various 5-forms.

Finally: try building T inversely — compute what T must be, then see if any
form produces it.
"""
import sys; sys.path.insert(0, 'src')
import sympy as sp
from sympy import Rational as R
from sugra import (HarmonicFunction, warped_product,
                   FormField, exterior_derivative,
                   form_stress_energy, Verifier)

# --- Setup ---
wv_coords = list(sp.symbols('t x1', real=True))
z1, z2 = sp.symbols('z1 z2', real=True)
harmonic_coords = list(sp.symbols('y0:8', real=True))
coords = wv_coords + [z1, z2] + harmonic_coords

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

# --- Test: D1 4-form (scaled) + symmetrized 5-form from A_{x1,z1,z2,yk} ---
print("="*60)
print("Building symmetrized 5-form: sum_k A_{x1,z1,z2,yk} = c*H^p")
print("="*60)

# Build the combined potential
for p_val in [R(-1,4), R(-1,2), R(-3,4), R(-1), R(1,4), R(1,2)]:
    A4 = FormField(rank=4, dim=12)
    for k in range(8):
        # A_{x1,z1,z2,yk} = H^p for each yk
        A4[(1, 2, 3, 4+k)] = H**p_val
    F5 = exterior_derivative(A4, coords)

    n_comp = len(F5.nonzero_components)
    if n_comp == 0:
        continue

    T5 = form_stress_energy(F5, metric)
    print(f"\np={p_val}, F5 has {n_comp} components")

    for i, name in [(0,'t'), (2,'z1'), (3,'z2'), (4,'y0')]:
        t_val = hf.substitute(sp.cancel(T5[i,i]))
        r_val = hf.substitute(sp.cancel(Ric[i,i]))
        print(f"  T5[{name}] = {t_val}")
        print(f"  R[{name}]  = {r_val}")
        print()

# --- Test: Try D1 + F1 combination at various amplitudes ---
print("\n" + "="*60)
print("Parametric scan: a²*T_D1 + b²*T_F1 = R")
print("="*60)

# Compute T for unit-amplitude D1 and F1
C_D1 = FormField(rank=3, dim=12)
C_D1[(0,1,3)] = 1/H
F_D1 = exterior_derivative(C_D1, coords)
T_D1 = form_stress_energy(F_D1, metric)

# F1 at various powers
for p_f1 in [R(-1), R(-1,2), R(-3,4), R(-1,4)]:
    C_F1 = FormField(rank=3, dim=12)
    C_F1[(0,1,2)] = H**p_f1
    F_F1 = exterior_derivative(C_F1, coords)
    T_F1 = form_stress_energy(F_F1, metric)

    # Check if F1 at this power has matching H-powers
    t_test = hf.substitute(sp.cancel(T_F1[0,0]))
    if t_test == 0:
        continue

    # For D1: T_D1 is computed
    # Try to solve: a²*T_D1[block] + b²*T_F1[block] = R[block]
    # for each pair of blocks, check consistency

    # Use symbolic a², b²
    a2, b2 = sp.symbols('a2 b2', positive=True)

    eqs = []
    for i, name in [(0,'t'), (2,'z1'), (3,'z2')]:
        r_val = hf.substitute(sp.cancel(Ric[i,i]))
        td = hf.substitute(sp.cancel(T_D1[i,i]))
        tf = hf.substitute(sp.cancel(T_F1[i,i]))
        eq = sp.Eq(a2*td + b2*tf, r_val)
        eqs.append(eq)

    # Solve using t and z1 equations
    try:
        sol = sp.solve([eqs[0], eqs[1]], [a2, b2])
        if sol:
            a2_val = sol[a2]
            b2_val = sol[b2]
            # Check z2
            td_z2 = hf.substitute(sp.cancel(T_D1[3,3]))
            tf_z2 = hf.substitute(sp.cancel(T_F1[3,3]))
            r_z2 = hf.substitute(sp.cancel(Ric[3,3]))
            check_z2 = sp.cancel(a2_val*td_z2 + b2_val*tf_z2 - r_z2)
            if check_z2 == 0:
                print(f"\n  F1 power={p_f1}: EXACT MATCH on t,z1,z2!")
                print(f"    a² = {a2_val}, b² = {b2_val}")
            else:
                print(f"\n  F1 power={p_f1}: a²={a2_val}, b²={b2_val}")
                print(f"    z2 residual: {check_z2}")
    except:
        pass

# --- Test: Direct inverse approach ---
# Compute what X_{MN} = FF_{MN}/(n-1)! must be for a 4-form
print("\n" + "="*60)
print("Inverse analysis: what must FF_{MN} be for a 4-form?")
print("="*60)

# From T_{MN} = (1/2)[X_{MN} - 3/10*S*g_{MN}] = R_{MN}
# Taking trace: 16S = 2R + 36S/5 → S = 5R/22

R_scalar = sp.S(0)
for i in range(12):
    r_val = hf.substitute(sp.cancel(Ric[i,i]))
    g_inv = sp.cancel(1/metric.matrix[i,i])
    R_scalar += sp.cancel(r_val * g_inv)
R_scalar = sp.cancel(R_scalar)
print(f"Ricci scalar R = {R_scalar}")

# For a 4-form: S = 5R/22
S_target = sp.cancel(5*R_scalar/22)
print(f"Required |F²| = S = 5R/22 = {S_target}")

# D1's |F²|
from sugra import form_norm_squared
F_sq_D1 = form_norm_squared(F_D1, metric)
F_sq_D1_sub = hf.substitute(sp.cancel(F_sq_D1))
print(f"D1's |F²| = {F_sq_D1_sub}")
print(f"Ratio D1/target = {sp.cancel(F_sq_D1_sub/S_target)}")

# Required X_{MN} = 2R_{MN} + (3/5)*S*g_{MN}
print("\nRequired X_{MN} (= FF_{MN}/(n-1)!):")
for i, name in [(0,'t'), (2,'z1'), (3,'z2'), (4,'y0')]:
    r_val = hf.substitute(sp.cancel(Ric[i,i]))
    g_val = metric.matrix[i,i]
    g_sub = hf.substitute(sp.cancel(g_val))
    x_val = sp.cancel(2*r_val + R(3,5)*S_target*g_sub)
    print(f"  X[{name}] = {x_val}")

# D1's X_{MN}
from sugra import form_contraction
FF_D1 = form_contraction(F_D1, metric)
for i, name in [(0,'t'), (2,'z1'), (3,'z2'), (4,'y0')]:
    ff_val = hf.substitute(sp.cancel(FF_D1[i,i] / sp.factorial(3)))
    print(f"  D1's X[{name}] = {ff_val}")
