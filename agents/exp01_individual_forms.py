"""Experiment 01: Test individual 4-forms from F1 and D1 potentials.

F1-type: C_{t,x1,z1} = H^{-a}  →  F4 with legs (t,x1,z1,yi)
D1-type: C_{t,x1,z2} = H^{-b}  →  F4 with legs (t,x1,z2,yi)

Goal: See what H-powers each produces in T_{MN} and compare to R_{MN}.
"""
import sys; sys.path.insert(0, 'src')
import sympy as sp
from sugra import (HarmonicFunction, warped_product,
                   FormField, exterior_derivative,
                   form_stress_energy)

# --- Coordinates ---
d_wv, d_harm = 2, 8
wv_coords = list(sp.symbols('t x1', real=True))
z1, z2 = sp.symbols('z1 z2', real=True)
harmonic_coords = list(sp.symbols('y0:8', real=True))
coords = wv_coords + [z1, z2] + harmonic_coords  # D=12

# --- Harmonic function ---
hf = HarmonicFunction(transverse_coords=harmonic_coords)
H = sp.Function('H')(hf.r_expr)

# --- 12d metric ---
metric = warped_product(
    warp_factors=[H**sp.Rational(-3,4), H**sp.Rational(1,2),
                  H**sp.Rational(-1,2), H**sp.Rational(1,4)],
    block_dims=[2, 1, 1, 8],
    block_signatures=['lorentzian','euclidean','euclidean','euclidean'],
    coordinates=coords,
)

# --- Compute Ricci tensor ---
print("Computing Ricci tensor...")
R = metric.ricci_tensor(simplify_func=sp.cancel)
print("Ricci tensor done.\n")

# Print R targets
print("=== Ricci tensor targets ===")
for i, name in [(0,'t'), (2,'z1'), (3,'z2'), (4,'y0')]:
    r = hf.substitute(sp.cancel(R[i,i]))
    print(f"  R[{name},{name}] = {r}")
print()

def test_form(label, potential_indices, potential_power):
    """Test a 3-form potential and its 4-form field strength."""
    print(f"=== {label}: C_{potential_indices} = H^({potential_power}) ===")
    C = FormField(rank=3, dim=12)
    C[potential_indices] = H**potential_power
    F = exterior_derivative(C, coords)

    print(f"  Nonzero F4 components: {len(F.nonzero_components)}")
    for idx, val in F.nonzero_components.items():
        v = hf.substitute(sp.cancel(val))
        print(f"    F{list(idx)} = {v}")

    T = form_stress_energy(F, metric)

    print(f"\n  Stress-energy T_{label}:")
    for i, name in [(0,'t'), (2,'z1'), (3,'z2'), (4,'y0')]:
        t_val = hf.substitute(sp.cancel(T[i,i]))
        r_val = hf.substitute(sp.cancel(R[i,i]))
        diff = sp.cancel(r_val - t_val)
        print(f"    T[{name},{name}] = {t_val}")
        print(f"    R[{name},{name}] = {r_val}")
        print(f"    R - T      = {diff}")
        print()
    return F, T

# --- Experiment 1a: F1-type with z1 leg ---
F1, T1 = test_form("F1", (0, 1, 2), -1)  # C_{t,x1,z1} = H^{-1}

# --- Experiment 1b: D1-type with z2 leg ---
F2, T2 = test_form("D1", (0, 1, 3), -1)  # C_{t,x1,z2} = H^{-1}

# --- Experiment 1c: Try different powers ---
print("\n" + "="*60)
print("Trying different H-powers for the potential...")
for p in [sp.Rational(-1,2), sp.Rational(-3,4), sp.Rational(-1,4)]:
    test_form(f"F1(H^{p})", (0, 1, 2), p)
