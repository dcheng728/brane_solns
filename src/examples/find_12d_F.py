"""
Adaptive search for form field(s) F (no dilaton) balancing R_{MN} = T_{MN} for:

  ds^2_{12} = H^{-3/4} ds^2_{1,1} + H^{1/2} dz_1^2 + H^{-1/2} dz_2^2 + H^{1/4} ds^2_8

Strategy:
  1. Compute R_{MN} once.
  2. For each single-form ansatz with C = H^p (p = -1), compute T_{MN} and the
     ratio R/T to diagnose what's off.
  3. Try pairs of forms: if cross terms vanish, T = c1^2 T_1 + c2^2 T_2 = R
     is a linear system in (c1^2, c2^2).
  4. Verify any solution found with the full Verifier.
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import sympy as sp
from itertools import combinations
from sugra import (HarmonicFunction, warped_product,
                   FormField, exterior_derivative, hodge_star,
                   form_stress_energy, Verifier)

# ============================================================================
# Setup (fixed)
# ============================================================================

d_wv = 2
d_harm = 8

wv_coords = list(sp.symbols('t x1', real=True))
z1, z2 = sp.symbols('z1 z2', real=True)
harmonic_coords = list(sp.symbols(f'y0:{d_harm}', real=True))

coords = wv_coords + [z1, z2] + harmonic_coords
D = len(coords)   # 12

hf = HarmonicFunction(transverse_coords=harmonic_coords)
H_func = sp.Function('H')(hf.r_expr)

metric = warped_product(
    warp_factors     = [H_func**sp.Rational(-3, 4),
                        H_func**sp.Rational(1, 2),
                        H_func**sp.Rational(-1, 2),
                        H_func**sp.Rational(1, 4)],
    block_dims       = [d_wv, 1, 1, d_harm],
    block_signatures = ['lorentzian', 'euclidean', 'euclidean', 'euclidean'],
    coordinates      = coords,
)

# Representative index per block: t(wv), z1(torus1), z2(torus2), y0(transverse)
# x1 is redundant with t (same block), so we use 4 independent constraints.
REPS = [(0, 't'), (2, 'z1'), (3, 'z2'), (4, 'y0')]

# ============================================================================
# Phase 0: Compute R_{MN}
# ============================================================================

print("Phase 0: Computing Ricci tensor ...")
R_full = metric.ricci_tensor(simplify_func=sp.cancel)

R = {}
for idx, name in REPS:
    R[name] = hf.substitute(sp.cancel(R_full[idx, idx]))
    print(f"  R[{name},{name}] = {R[name]}")
print()

# ============================================================================
# Helpers
# ============================================================================

def make_F(legs, power, form_type='electric'):
    """Build F from potential C_{legs} = H^power."""
    C = FormField(rank=len(legs), dim=D)
    C[legs] = H_func**power
    F_E = exterior_derivative(C, coords)
    if form_type == 'electric':
        return F_E
    else:
        return hodge_star(F_E, metric)

def compute_T_reps(F):
    """Compute T_{ii} for representative indices, return as dict."""
    T_full = form_stress_energy(F, metric)
    T = {}
    for idx, name in REPS:
        T[name] = hf.substitute(sp.cancel(T_full[idx, idx]))
    return T

def diagnose(label, T):
    """Compare T to R for each representative, print ratios."""
    print(f"  {label}")
    ratios = {}
    for name in [n for _, n in REPS]:
        if T[name] == 0 and R[name] == 0:
            ratios[name] = 'both_zero'
            print(f"    {name}: T = 0, R = 0")
        elif T[name] == 0:
            ratios[name] = 'T_zero'
            print(f"    {name}: T = 0, R = {R[name]}  [T vanishes!]")
        else:
            ratio = sp.cancel(R[name] / T[name])
            ratios[name] = ratio
            diff = sp.cancel(R[name] - T[name])
            print(f"    {name}: R/T = {ratio}")
    return ratios

# ============================================================================
# Phase 1: Single forms with p = -1
# ============================================================================

# All subsets of {t, x1, z1, z2} that include t
def all_legs():
    others = [1, 2, 3]  # x1, z1, z2
    for r in range(len(others) + 1):
        for subset in combinations(others, r):
            yield (0,) + subset

print("Phase 1: Single forms, C = H^{-1}, diagnosing R/T ratios ...")
print()

single_T = {}  # (legs, type) -> T dict

for legs in all_legs():
    leg_names = ','.join(str(coords[i]) for i in legs)
    for ftype in ['electric', 'magnetic']:
        key = (legs, ftype)
        label_prefix = 'F' if ftype == 'electric' else '*F'
        label = f"{label_prefix} = dC({leg_names}), p=-1"

        F = make_F(legs, sp.S(-1), ftype)
        T = compute_T_reps(F)
        single_T[key] = T
        diagnose(label, T)
        print()

# ============================================================================
# Phase 2: Pairs of forms (no cross terms assumed)
# ============================================================================
#
# If F_1 and F_2 have non-overlapping index structures, cross terms in T vanish:
#   T_total = c1^2 * T_1 + c2^2 * T_2 = R
#
# This is a linear system in (c1^2, c2^2) with 4 equations (one per block).
# We need a consistent non-negative solution.

print("=" * 60)
print("Phase 2: Pairs of forms, solving c1^2 T_1 + c2^2 T_2 = R ...")
print("=" * 60)
print()

c1sq, c2sq = sp.symbols('c1sq c2sq', positive=True)
rep_names = [n for _, n in REPS]

keys = list(single_T.keys())
found_solutions = []

for i, key1 in enumerate(keys):
    for key2 in keys[i+1:]:
        T1 = single_T[key1]
        T2 = single_T[key2]

        # Skip if either has all zeros
        if all(T1[n] == 0 for n in rep_names):
            continue
        if all(T2[n] == 0 for n in rep_names):
            continue

        # Build the system: c1sq * T1[n] + c2sq * T2[n] = R[n]
        equations = []
        for name in rep_names:
            eq = c1sq * T1[name] + c2sq * T2[name] - R[name]
            eq = sp.cancel(eq)
            if eq != 0:
                equations.append(eq)

        if not equations:
            # Trivial match (shouldn't happen)
            continue

        try:
            sol = sp.solve(equations, [c1sq, c2sq], dict=True)
        except Exception:
            continue

        for s in sol:
            if c1sq in s and c2sq in s:
                v1 = s[c1sq]
                v2 = s[c2sq]
                # Check: both must be real, non-negative, and constant (no H, r, y dependence)
                free = v1.free_symbols | v2.free_symbols
                hf_syms = {hf.H, hf.Hp, hf.r} | set(hf._y)
                if free & hf_syms:
                    continue  # depends on H, r, etc. -- not valid
                if v1.is_negative or v2.is_negative:
                    continue

                legs1, type1 = key1
                legs2, type2 = key2
                l1 = ','.join(str(coords[j]) for j in legs1)
                l2 = ','.join(str(coords[j]) for j in legs2)
                prefix1 = '' if type1 == 'electric' else '*'
                prefix2 = '' if type2 == 'electric' else '*'

                print(f"  [SOLUTION CANDIDATE]")
                print(f"    F1: {prefix1}dC({l1}), c1^2 = {v1}")
                print(f"    F2: {prefix2}dC({l2}), c2^2 = {v2}")
                print()

                found_solutions.append((key1, v1, key2, v2))

# ============================================================================
# Phase 3: Verify solutions
# ============================================================================

if found_solutions:
    print("=" * 60)
    print("Phase 3: Full verification ...")
    print("=" * 60)
    print()

    for key1, v1, key2, v2 in found_solutions:
        legs1, type1 = key1
        legs2, type2 = key2

        c1 = sp.sqrt(v1)
        c2 = sp.sqrt(v2)

        F1 = make_F(legs1, sp.S(-1), type1) * c1
        F2 = make_F(legs2, sp.S(-1), type2) * c2
        F = F1 + F2

        l1 = ','.join(str(coords[j]) for j in legs1)
        l2 = ','.join(str(coords[j]) for j in legs2)
        prefix1 = '' if type1 == 'electric' else '*'
        prefix2 = '' if type2 == 'electric' else '*'
        print(f"  Verifying: {c1}*{prefix1}dC({l1}) + {c2}*{prefix2}dC({l2})")
        print()

        soln = dict(metric=metric, F=F, Phi=0, alpha=0, coords=coords, hf=hf)
        ok = Verifier(soln).check()
        if ok:
            print("  >>> SOLUTION FOUND <<<")
        print()
else:
    print()
    print("No pair solution found with p = -1.")
    print()

# ============================================================================
# Summary
# ============================================================================

print("=" * 60)
print("  SUMMARY")
print("=" * 60)
if found_solutions:
    print(f"  Found {len(found_solutions)} candidate(s).")
else:
    print("  No solution found from single forms or pairs with C = H^{-1}.")
    print("  Next steps to try:")
    print("    - Different powers p != -1")
    print("    - Three or more forms")
    print("    - Forms with transverse legs in the potential")
print()
