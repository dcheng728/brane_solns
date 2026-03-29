"""Experiment 04: Solve the linear system a²T₁ + b²T₂ + c²T₃ = R.

The 4-form can have BOTH D1 and F1 components (cross terms vanish).
The 5-form adds a third independent contribution.

We systematically check ALL combinations of H-power-matching forms.
Also try forms with exotic Y-leg structures that weren't in the initial scan.
"""
import sys; sys.path.insert(0, 'src')
import sympy as sp
from sympy import Rational as R
from sugra import (HarmonicFunction, warped_product,
                   FormField, exterior_derivative,
                   form_stress_energy)
from itertools import combinations

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

# Extract reference values
Hp, Hsym, rsym = hf.Hp, hf.H, hf.r
y0sym = sp.Symbol('y0', real=True)

R_exprs = {}
for i, name in [(0,'t'), (2,'z1'), (3,'z2'), (4,'y0')]:
    R_exprs[name] = hf.substitute(sp.cancel(Ric[i,i]))
print("R[t]  =", R_exprs['t'])
print("R[z1] =", R_exprs['z1'])
print("R[z2] =", R_exprs['z2'])
print("R[y0] =", R_exprs['y0'])

# For y0 block, separate into isotropic and anisotropic parts
# R[y0] = (1 - 4*y0²/r²) * H'^2/(8H²) → iso=1/8, aniso=-4/8=-1/2 (per H'^2/H²)

def extract_5_coeffs(T_exprs, R_target_Hpowers):
    """Extract the 5 scalar coefficients from T expression.

    For each block, T should have the form c * H'^2 / H^n (possibly with y0 terms).
    Returns (t_coeff, z1_coeff, z2_coeff, y0_iso_coeff, y0_aniso_coeff)
    or None if H-powers don't match.
    """
    coeffs = {}
    for name in ['t', 'z1', 'z2']:
        expr = T_exprs[name]
        if expr == 0:
            return None
        # Substitute H'=1, r=1, y0=0
        val = expr.subs([(Hp, 1), (rsym, 1), (y0sym, 0)])
        coeffs[name] = sp.cancel(val)

    # y0 block: separate iso and aniso
    y0_expr = T_exprs['y0']
    # iso: set y0=0
    iso = y0_expr.subs([(Hp, 1), (rsym, 1), (y0sym, 0)])
    # full: set y0=1/2, r=1 → iso + aniso*(1/4)
    full = y0_expr.subs([(Hp, 1), (rsym, 1), (y0sym, sp.Rational(1,2))])
    aniso = sp.cancel(4*(full - iso))  # aniso = coefficient of y0²/r²

    coeffs['y0_iso'] = sp.cancel(iso)
    coeffs['y0_aniso'] = sp.cancel(aniso)

    return coeffs

def check_hpower_match(T_exprs, R_exprs):
    """Check if T has the same H-power dependence as R in each block."""
    H_test = [R(2), R(3)]
    for name in ['t', 'z1', 'z2', 'y0']:
        r_vals = []
        t_vals = []
        for h in H_test:
            subs = [(Hsym, h), (Hp, 1), (rsym, 1), (y0sym, R(1,3))]
            rv = R_exprs[name].subs(subs)
            tv = T_exprs[name].subs(subs)
            r_vals.append(rv)
            t_vals.append(tv)

        if t_vals[0] == 0 or t_vals[1] == 0:
            return False

        ratio1 = sp.cancel(t_vals[0] / r_vals[0])
        ratio2 = sp.cancel(t_vals[1] / r_vals[1])
        if sp.cancel(ratio1 - ratio2) != 0:
            return False
    return True

# --- Broad scan of potentials ---
print("\n" + "="*70)
print("Extended scan: trying many potential forms and powers")
print("="*70)

# Powers to try (extended range)
powers = [R(p,4) for p in range(-10, 5)]  # -10/4 to 4/4

matching_forms = []

# 3-form potentials → 4-form F
potential_3forms = [
    ("C_{t,x1,z2}", (0,1,3)),
    ("C_{t,x1,z1}", (0,1,2)),
    ("C_{t,z1,z2}", (0,2,3)),
    ("C_{x1,z1,z2}", (1,2,3)),
]

# 4-form potentials → 5-form F
potential_4forms = [
    ("A_{t,x1,z1,z2}", (0,1,2,3)),
    ("A_{t,x1,z1,y0}", (0,1,2,4)),
    ("A_{t,x1,z2,y0}", (0,1,3,4)),
    ("A_{t,z1,z2,y0}", (0,2,3,4)),
    ("A_{x1,z1,z2,y0}", (1,2,3,4)),
    ("A_{t,x1,y0,y1}", (0,1,4,5)),
    ("A_{z1,z2,y0,y1}", (2,3,4,5)),
    ("A_{t,z1,y0,y1}", (0,2,4,5)),
    ("A_{t,z2,y0,y1}", (0,3,4,5)),
    ("A_{y0,y1,y2,y3}", (4,5,6,7)),
]

# 2-form potentials → 3-form F
potential_2forms = [
    ("B_{t,x1}", (0,1)),
    ("B_{t,z1}", (0,2)),
    ("B_{t,z2}", (0,3)),
    ("B_{z1,z2}", (2,3)),
    ("B_{t,y0}", (0,4)),
]

all_potentials = (
    [(name, 3, idx) for name, idx in potential_3forms] +
    [(name, 4, idx) for name, idx in potential_4forms] +
    [(name, 2, idx) for name, idx in potential_2forms]
)

for name, rank, indices in all_potentials:
    fs_rank = rank + 1  # field strength rank
    for p in powers:
        try:
            C = FormField(rank=rank, dim=12)
            C[indices] = H**p
            F = exterior_derivative(C, coords)
            if not F.nonzero_components:
                continue
            T = form_stress_energy(F, metric)

            T_exprs = {}
            for i, bname in [(0,'t'), (2,'z1'), (3,'z2'), (4,'y0')]:
                T_exprs[bname] = hf.substitute(sp.cancel(T[i,i]))

            if check_hpower_match(T_exprs, R_exprs):
                coeffs = extract_5_coeffs(T_exprs, None)
                if coeffs is not None:
                    matching_forms.append({
                        'name': f"{name}=H^({p})",
                        'rank': fs_rank,
                        'coeffs': coeffs,
                        'T_exprs': T_exprs
                    })
                    print(f"  MATCH: {name}=H^({p}) [F{fs_rank}]")
                    print(f"    coeffs = {coeffs}")
        except Exception as e:
            pass

print(f"\nTotal matching forms: {len(matching_forms)}")

# --- Now try to solve the linear system ---
print("\n" + "="*70)
print("Trying ALL pairs and triples to solve a²T₁ + b²T₂ [+ c²T₃] = R")
print("="*70)

# Target coefficients
R_coeffs = extract_5_coeffs(R_exprs, None)
print(f"R target: {R_coeffs}")

def try_solve(forms_list, R_target):
    """Try to solve sum_i w_i * T_i = R_target for w_i > 0."""
    n = len(forms_list)
    keys = ['t', 'z1', 'z2', 'y0_iso', 'y0_aniso']

    # Build matrix: each column is a form's coefficients
    A = sp.Matrix(len(keys), n, lambda i,j: forms_list[j]['coeffs'][keys[i]])
    b = sp.Matrix(len(keys), 1, lambda i,j: R_target[keys[i]])

    # Variables
    w = sp.symbols(f'w0:{n}', positive=True)

    # Solve the system Aw = b
    try:
        # For 2 unknowns, 5 equations: use least squares or check consistency
        if n <= len(keys):
            sol = A.solve(b)
            if sol is not None:
                # Check all positive
                all_pos = all(sp.cancel(s) > 0 for s in sol if s != 0)
                # Check residual
                residual = sp.cancel(A * sol - b)
                if residual.equals(sp.zeros(len(keys), 1)):
                    return sol, True
                return sol, False
    except Exception:
        pass

    # Manual solve for n=2: use first 2 equations, check rest
    if n == 2:
        A2 = A[:2, :]
        b2 = b[:2, :]
        try:
            sol2 = A2.solve(b2)
            residual = sp.cancel(A * sol2 - b)
            is_exact = residual.equals(sp.zeros(len(keys), 1))
            return sol2, is_exact
        except:
            return None, False

    if n == 3:
        A3 = A[:3, :]
        b3 = b[:3, :]
        try:
            sol3 = A3.solve(b3)
            residual = sp.cancel(A * sol3 - b)
            is_exact = residual.equals(sp.zeros(len(keys), 1))
            return sol3, is_exact
        except:
            return None, False

    return None, False

# Try all pairs
print("\n--- Pairs ---")
for i, j in combinations(range(len(matching_forms)), 2):
    fi, fj = matching_forms[i], matching_forms[j]
    sol, exact = try_solve([fi, fj], R_coeffs)
    if sol is not None and exact:
        w1, w2 = sol
        if sp.cancel(w1) > 0 and sp.cancel(w2) > 0:
            print(f"  SOLUTION! {fi['name']} (w²={w1}) + {fj['name']} (w²={w2})")

# Try all triples
print("\n--- Triples ---")
for i, j, k in combinations(range(len(matching_forms)), 3):
    fi, fj, fk = matching_forms[i], matching_forms[j], matching_forms[k]
    sol, exact = try_solve([fi, fj, fk], R_coeffs)
    if sol is not None and exact:
        w1, w2, w3 = sol
        if sp.cancel(w1) > 0 and sp.cancel(w2) > 0 and sp.cancel(w3) > 0:
            print(f"  SOLUTION! {fi['name']} + {fj['name']} + {fk['name']}")
            print(f"    weights: {w1}, {w2}, {w3}")
