"""Experiment 02: Test 5-form and two-form combinations.

Key findings from exp01:
- D1-type 4-form (z2 leg) matches all H-powers in R_{MN}
- Coefficient mismatch: T/R ratios are 14/15, 3/5, 7/5 for t, z1, z2
- Need a second form to fix coefficients

Test:
(a) 5-form with legs (t,x1,z1,z2,yi): A_{t,x1,z1,z2} = H^s
(b) Combined D1-type 4-form + 5-form
(c) Parameterize powers and coefficients to solve for exact match
"""
import sys; sys.path.insert(0, 'src')
import sympy as sp
from sugra import (HarmonicFunction, warped_product,
                   FormField, exterior_derivative,
                   form_stress_energy)

# --- Setup (same as exp01) ---
wv_coords = list(sp.symbols('t x1', real=True))
z1, z2 = sp.symbols('z1 z2', real=True)
harmonic_coords = list(sp.symbols('y0:8', real=True))
coords = wv_coords + [z1, z2] + harmonic_coords

hf = HarmonicFunction(transverse_coords=harmonic_coords)
H = sp.Function('H')(hf.r_expr)

metric = warped_product(
    warp_factors=[H**sp.Rational(-3,4), H**sp.Rational(1,2),
                  H**sp.Rational(-1,2), H**sp.Rational(1,4)],
    block_dims=[2, 1, 1, 8],
    block_signatures=['lorentzian','euclidean','euclidean','euclidean'],
    coordinates=coords,
)

print("Computing Ricci tensor...")
R = metric.ricci_tensor(simplify_func=sp.cancel)
print("Done.\n")

# Cache R values
R_vals = {}
for i, name in [(0,'t'), (2,'z1'), (3,'z2'), (4,'y0')]:
    R_vals[name] = hf.substitute(sp.cancel(R[i,i]))
    print(f"R[{name},{name}] = {R_vals[name]}")
print()

# === Test 2a: 5-form with all four special legs ===
print("=" * 60)
print("Test 2a: 5-form F5 from A_{t,x1,z1,z2} = H^{-1}")
print("=" * 60)

A4 = FormField(rank=4, dim=12)
A4[(0,1,2,3)] = 1/H  # A_{t,x1,z1,z2} = H^{-1}
F5 = exterior_derivative(A4, coords)

print(f"Nonzero F5 components: {len(F5.nonzero_components)}")
for idx, val in F5.nonzero_components.items():
    v = hf.substitute(sp.cancel(val))
    print(f"  F5{list(idx)} = {v}")

T5 = form_stress_energy(F5, metric)
print(f"\nStress-energy from 5-form:")
for i, name in [(0,'t'), (2,'z1'), (3,'z2'), (4,'y0')]:
    t_val = hf.substitute(sp.cancel(T5[i,i]))
    print(f"  T5[{name},{name}] = {t_val}")
    print(f"  R[{name},{name}]  = {R_vals[name]}")
    print(f"  R - T5       = {sp.cancel(R_vals[name] - t_val)}")
    print()

# === Test 2b: Combined D1-type 4-form + 5-form ===
print("=" * 60)
print("Test 2b: D1-type 4-form + 5-form (both H^{-1})")
print("=" * 60)

# D1-type 4-form
C3_D1 = FormField(rank=3, dim=12)
C3_D1[(0,1,3)] = 1/H
F4_D1 = exterior_derivative(C3_D1, coords)
T4 = form_stress_energy(F4_D1, metric)

print("T_total = T4(D1) + T5:")
for i, name in [(0,'t'), (2,'z1'), (3,'z2'), (4,'y0')]:
    t4 = hf.substitute(sp.cancel(T4[i,i]))
    t5 = hf.substitute(sp.cancel(T5[i,i]))
    total = sp.cancel(t4 + t5)
    diff = sp.cancel(R_vals[name] - total)
    print(f"  T4[{name}] = {t4}")
    print(f"  T5[{name}] = {t5}")
    print(f"  Total   = {total}")
    print(f"  R[{name}] = {R_vals[name]}")
    print(f"  R-T     = {diff}")
    print()

# === Test 2c: Parameterized combination ===
print("=" * 60)
print("Test 2c: a*F4(D1,H^p) + b*F5(H^q) — find a,b,p,q")
print("=" * 60)

# For a general potential H^p, the field strength is p*H^{p-1}*H'*yi/r
# T scales as (coefficient)^2 * p^2 relative to p=-1 case (where p^2=1)
# But H-power in T also changes. Let's compute T for general p.

# For the D1 4-form with C = H^p, the field is F ~ p*H^{p-1}*H'.
# T has terms like F*F ~ p^2*H^{2(p-1)}*H'^2.
# The metric contractions bring additional H-powers.
# Let me compute symbolically for a few powers.

print("\nD1-type 4-form T[block] as function of potential power p:")
for p_val in [-1, sp.Rational(-1,2), sp.Rational(-3,4), sp.Rational(-1,4), sp.Rational(-3,2)]:
    C3 = FormField(rank=3, dim=12)
    C3[(0,1,3)] = H**p_val
    F4 = exterior_derivative(C3, coords)
    T = form_stress_energy(F4, metric)
    t_t = hf.substitute(sp.cancel(T[0,0]))
    t_z1 = hf.substitute(sp.cancel(T[2,2]))
    t_z2 = hf.substitute(sp.cancel(T[3,3]))
    t_y0 = hf.substitute(sp.cancel(T[4,4]))
    print(f"  p={p_val}:")
    print(f"    T[t]  = {t_t}")
    print(f"    T[z1] = {t_z1}")
    print(f"    T[z2] = {t_z2}")
    print(f"    T[y0] = {t_y0}")
    print()

print("\n5-form T[block] as function of potential power q:")
for q_val in [-1, sp.Rational(-1,2), sp.Rational(-3,4), sp.Rational(-1,4), sp.Rational(-3,2)]:
    A4 = FormField(rank=4, dim=12)
    A4[(0,1,2,3)] = H**q_val
    F5 = exterior_derivative(A4, coords)
    T = form_stress_energy(F5, metric)
    t_t = hf.substitute(sp.cancel(T[0,0]))
    t_z1 = hf.substitute(sp.cancel(T[2,2]))
    t_z2 = hf.substitute(sp.cancel(T[3,3]))
    t_y0 = hf.substitute(sp.cancel(T[4,4]))
    print(f"  q={q_val}:")
    print(f"    T[t]  = {t_t}")
    print(f"    T[z1] = {t_z1}")
    print(f"    T[z2] = {t_z2}")
    print(f"    T[y0] = {t_y0}")
    print()
