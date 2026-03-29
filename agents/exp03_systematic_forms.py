"""Experiment 03: Systematic scan of form index structures.

For each possible leg configuration, find the H-power p in the potential
that makes T have the correct H-powers as R, then record the T/R ratios.

We need two forms whose T/R ratios, when combined with weights a² and b²,
can solve a²*(T1/R) + b²*(T2/R) = 1 for all four blocks.
"""
import sys; sys.path.insert(0, 'src')
import sympy as sp
from sugra import (HarmonicFunction, warped_product,
                   FormField, exterior_derivative,
                   form_stress_energy)

# --- Setup ---
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

R_vals = {}
for i, name in [(0,'t'), (2,'z1'), (3,'z2'), (4,'y0')]:
    R_vals[name] = hf.substitute(sp.cancel(R[i,i]))
    print(f"R[{name}] = {R_vals[name]}")
print()

# --- Helper: test a potential and extract T/R ratios ---
def test_potential(label, rank, indices, power):
    """Test a potential form C_indices = H^power, compute F=dC, T(F)."""
    C = FormField(rank=rank, dim=12)
    C[indices] = H**power
    F = exterior_derivative(C, coords)
    if not F.nonzero_components:
        print(f"  {label}: F=0 (trivial)")
        return None
    T = form_stress_energy(F, metric)
    results = {}
    for i, name in [(0,'t'), (2,'z1'), (3,'z2'), (4,'y0')]:
        t_val = hf.substitute(sp.cancel(T[i,i]))
        results[name] = t_val
    return results

# --- Systematic scan ---
# Potential indices → (rank, indices, description)
# Convention: 0=t, 1=x1, 2=z1, 3=z2, 4+=yi

ansatze = [
    # 3-form potentials → 4-form F
    ("C_{t,x1,z2}", 3, (0,1,3)),      # D1-type
    ("C_{t,x1,z1}", 3, (0,1,2)),      # F1-type
    ("C_{t,z1,z2}", 3, (0,2,3)),      # mixed-1
    ("C_{x1,z1,z2}", 3, (1,2,3)),     # mixed-2
    ("C_{t,x1,y0}", 3, (0,1,4)),      # WV+transverse
    ("C_{t,z1,y0}", 3, (0,2,4)),      # WV+z1+transverse
    ("C_{t,z2,y0}", 3, (0,3,4)),      # WV+z2+transverse
    ("C_{z1,z2,y0}", 3, (2,3,4)),     # z1+z2+transverse

    # 4-form potentials → 5-form F
    ("A_{t,x1,z1,z2}", 4, (0,1,2,3)), # all four special legs
    ("A_{t,x1,z1,y0}", 4, (0,1,2,4)), # WV+z1+transverse
    ("A_{t,x1,z2,y0}", 4, (0,1,3,4)), # WV+z2+transverse
    ("A_{t,z1,z2,y0}", 4, (0,2,3,4)), # mixed
    ("A_{x1,z1,z2,y0}", 4, (1,2,3,4)), # mixed-2

    # 2-form potentials → 3-form F
    ("B_{t,x1}", 2, (0,1)),           # WV only
    ("B_{t,z1}", 2, (0,2)),
    ("B_{t,z2}", 2, (0,3)),
    ("B_{z1,z2}", 2, (2,3)),
    ("B_{t,y0}", 2, (0,4)),
    ("B_{z1,y0}", 2, (2,4)),
    ("B_{z2,y0}", 2, (3,4)),

    # 5-form potentials → 6-form F
    ("D_{t,x1,z1,z2,y0}", 5, (0,1,2,3,4)),
]

# For each ansatz, scan powers to find one that gives correct H-powers
print("=" * 70)
print("Scanning all ansatze with various H-powers in the potential...")
print("=" * 70)

powers_to_try = [sp.Integer(-1), sp.Rational(-1,2), sp.Rational(-3,4),
                 sp.Rational(-1,4), sp.Rational(-3,2), sp.Rational(1,4),
                 sp.Rational(1,2), sp.Rational(-5,4), sp.Rational(-7,4)]

# Target H-powers
target_hpowers = {'t': -3, 'z1': sp.Rational(-7,4), 'z2': sp.Rational(-11,4), 'y0': -2}

for label, rank, indices in ansatze:
    print(f"\n--- {label} (rank-{rank} potential → rank-{rank+1} field strength) ---")
    found_match = False
    for p in powers_to_try:
        result = test_potential(label, rank, indices, p)
        if result is None:
            break

        # Check if H-powers match targets
        # Extract H-power from expression like coeff*H'^2/H^n or coeff*H'^2*y0^2/(H^n*r^2)
        t_expr = result['t']
        if t_expr == 0:
            continue

        # For the t-block, R has H'^2/H^3. T should also have H'^2/H^n for some n.
        # If we substitute H'=1, H=H_sym, we can read off the H-power.
        H_sym = hf.H
        Hp_sym = hf.Hp
        r_sym = hf.r
        y0_sym = sp.Symbol('y0', real=True)

        # Substitute H'=1, r=1, y0=0 to isolate the H-power
        t_test = t_expr.subs([(Hp_sym, 1), (r_sym, 1), (y0_sym, 0)])
        # t_test should be ~ const * H^(-n)
        # Check if it's a single power of H
        if t_test == 0:
            continue

        # Check if all blocks have correct H-powers by comparing at two H values
        import random
        H_val1, H_val2 = sp.Rational(2), sp.Rational(3)

        def eval_at_H(expr, h_val):
            return expr.subs([(H_sym, h_val), (Hp_sym, 1), (r_sym, 1), (y0_sym, sp.Rational(1,3))])

        match = True
        ratios = {}
        for name in ['t', 'z1', 'z2', 'y0']:
            r_v1 = eval_at_H(R_vals[name], H_val1)
            r_v2 = eval_at_H(R_vals[name], H_val2)
            t_v1 = eval_at_H(result[name], H_val1)
            t_v2 = eval_at_H(result[name], H_val2)

            if t_v1 == 0 or r_v1 == 0:
                match = False
                break

            # If H-powers match, then T/R should be the same at both H values
            ratio1 = sp.cancel(t_v1 / r_v1)
            ratio2 = sp.cancel(t_v2 / r_v2)
            if sp.cancel(ratio1 - ratio2) != 0:
                match = False
                break
            ratios[name] = ratio1

        if match and ratios:
            found_match = True
            print(f"  ** p={p}: H-powers MATCH! T/R ratios = {ratios}")
            # Also print actual T values
            for name in ['t', 'z1', 'z2', 'y0']:
                print(f"     T[{name}] = {result[name]}")

    if not found_match:
        print(f"  No H-power match found in scanned range.")
