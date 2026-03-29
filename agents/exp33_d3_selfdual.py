"""Experiment 33: D3-brane with self-dual F₅ in 12d.

The D3 F₅ is self-dual: F₅ = *₁₀F₅. In the 12d formulation:
- F₇ = F₅ ∧ dz₁ ∧ dz₂, and F₅ = *₁₂F₇
- The self-dual condition means F₅^{SD} = F₅^e + *₁₀F₅^e

For stress-energy: (F₅^{SD})²_{MN}/4! = 2·(F₅^e)²_{MN}/4!
because the cross terms between F₅^e and *F₅^e vanish.

And |F₅^{SD}|² = 0 (self-dual in Lorentzian ⟹ null norm).

So: ℛ = (1/2)·2·(F₅^e)²/4! = (F₅^e)²/4! = FF^e/4!

Actually this factor of 2 needs careful verification. Let me also check
whether the enhanced-metric formula with F₅ gives the right result.
"""
import sys; sys.path.insert(0, 'src')
import sympy as sp
from sympy import Rational as R, cancel
from sugra import (HarmonicFunction, warped_product, FormField,
                   exterior_derivative, form_contraction, form_norm_squared)

# D3 setup
wv4 = list(sp.symbols('t x1 x2 x3', real=True))
z1s, z2s = sp.symbols('z1 z2', real=True)
hc6 = list(sp.symbols('y0:6', real=True))
coords = wv4 + [z1s, z2s] + hc6
D = 12

hf = HarmonicFunction(transverse_coords=hc6)
H = sp.Function('H')(hf.r_expr)

# D3 12d metric (EF): H^{-1/2}ds²_{1,3} + dz₁² + dz₂² + H^{1/2}ds²_6
m = warped_product(
    warp_factors=[H**R(-1,2), sp.Integer(1), sp.Integer(1), H**R(1,2)],
    block_dims=[4, 1, 1, 6],
    block_signatures=['lorentzian', 'euclidean', 'euclidean', 'euclidean'],
    coordinates=coords,
)

print("="*60)
print("D3-BRANE: Self-dual F₅ analysis")
print("="*60)

# Ricci
print("Computing Ricci...")
Ric = m.ricci_tensor(simplify_func=cancel)

# Electric F₅ = dC₄ where C₄ = H^{-1} vol_{1,3}
# This gives F₅ with components (y_k, t, x1, x2, x3)
C4 = FormField(rank=4, dim=D)
C4[(0, 1, 2, 3)] = 1/H  # t, x1, x2, x3
F5e = exterior_derivative(C4, coords)

print("Computing electric F₅ contraction...")
FFe = form_contraction(F5e, m)
norme = form_norm_squared(F5e, m)
norme_val = hf.substitute(cancel(norme))
print(f"|F₅^e|² = {norme_val}")

# Magnetic F₅ = *₁₂F₇ where F₇ = F₅^e ∧ dz₁ ∧ dz₂
# Components: F₅^m has legs entirely in transverse space
# For D3: F₅^m_{y_i,...,y_j,z₁,z₂} ... actually *₁₀F₅^e has 5 transverse legs

# Let me build *₁₀F₅^e directly.
# F₅^e has components F_{y_k,t,x1,x2,x3} = ∂H^{-1}/∂y_k
# In 10d (without torus): *₁₀F₅^e_{m₁...m₅} = (1/5!)ε^{10d}_{m₁...m₅n₁...n₅}F^{n₁...n₅}
# This is complicated. Let me instead use the 12d approach.

# In 12d: F₇ = F₅^e ∧ dz₁ ∧ dz₂ has 7 components: (y_k, t, x1, x2, x3, z₁, z₂)
# *₁₂F₇ should be a 5-form: *₁₂F₇ = F₅^m
# *₁₂F₇_{m₁...m₅} = (1/7!)√{-g}ε_{m₁...m₅n₁...n₇}F₇^{n₁...n₇}

# For the self-dual relation: F₅ = *₁₂F₇
# The self-dual F₅^{SD} has both electric and magnetic parts.

# Key property: for self-dual F₅, (F₅^{SD})² = 2·(F₅^e)²
# This is because F₅^e and *F₅^e have disjoint leg structure
# (electric has 4 wv + 1 trans, magnetic has 5 trans + 1 trans... wait)

# Actually in 12d: F₅^e has legs (wv, trans), F₅^m = *₁₂(F₅^e ∧ dz₁∧dz₂) has legs...
# Let me compute *₁₂F₇ explicitly.

# F₇ components: for each k ∈ {0,...,5}:
# F₇[t, x1, x2, x3, z1, z2, y_k] = F₅^e[t,x1,x2,x3,y_k] × 1 × 1 (wedge with dz1∧dz2)
# = ∂_{y_k}(H^{-1})

# Actually, F₇ = F₅^e ∧ dz₁ ∧ dz₂. So F₇[M₁...M₇] = F₅^e[M₁...M₅] when M₆=z₁, M₇=z₂
# (with appropriate antisymmetrization)

# *₁₂F₇ is a 5-form in 12d. Its components live in the complement of F₇'s legs.
# F₇ has legs (t,x1,x2,x3,z1,z2,y_k). The complement (the other 5 indices) is
# the remaining 5 transverse coords {y_0,...,y_5} \ {y_k}.

# So F₅^m has legs purely in transverse space! And F₅^e has legs in (wv + 1 trans).
# These are DISJOINT — they share at most the transverse index.

# More precisely: F₅^e_{y_k,t,x1,x2,x3} ∝ ∂_{y_k}H^{-1}
# F₅^m_{y_{i1},...,y_{i5}} ∝ ε_{i1...i5 k} ∂_{y_k}H^{-1} / (sqrt(g) factors)

# The contraction FF^{SD}_{MN} = (1/4!)g^{...}F^{SD}_{M...}F^{SD}_{N...}
# For M=t: F^{SD}_{t,...} = F^e_{t,...} (magnetic part has no t-leg)
# So FF^{SD}_{tt} = FF^e_{tt} + FF^m_{tt}
# But F^m has no t-leg → FF^m_{tt} = 0
# So FF^{SD}_{tt} = FF^e_{tt}

# Hmm, that means (F₅^{SD})²_{tt} = (F₅^e)²_{tt}, NOT 2×.
# The factor of 2 comes from contracting F₅^e with *F₅^e?

# Let me think more carefully. F₅^{SD} = F₅^e + F₅^m.
# (F₅^{SD})² = (F₅^e)² + (F₅^m)² + 2(F₅^e·F₅^m)

# Cross term: (F₅^e·F₅^m)_{MN} = (1/4!)g^{...}F^e_{M...}F^m_{N...}
# For M=t: F^e_{t,...} has remaining 4 indices in (x1,x2,x3,y_k)
# F^m has NO wv legs, so F^m_{N,...} with those 4 indices = 0
# → cross term = 0 for wv blocks ✓

# For M=y_0: F^e_{y_0,...} has remaining 4 indices in (t,x1,x2,x3)
# F^m_{N,...} with exactly (t,x1,x2,x3) = 0 (F^m has no wv legs)
# → cross term = 0 for transverse blocks ✓

# So cross terms vanish. Then (F₅^{SD})² = (F₅^e)² + (F₅^m)².

# For wv blocks: (F₅^m)²_{tt} = 0 → (F₅^{SD})²_{tt} = (F₅^e)²_{tt}
# For trans blocks: (F₅^e)²_{y0y0} ≠ 0 AND (F₅^m)²_{y0y0} ≠ 0
# For torus blocks: both have 0 (neither has torus legs for D3)

# Hmm wait, for trans: F₅^e_{y_0,t,x1,x2,x3} ≠ 0 → (F₅^e)²_{y0y0} ≠ 0
# And F₅^m_{y_0,y_1,...,y_4} ≠ 0 → (F₅^m)²_{y0y0} ≠ 0
# So the self-dual FF is NOT just 2× the electric FF for transverse.

# The issue is more subtle. In 10d: for self-dual F₅, the equation is
# R_{mn} = (1/4·4!)(F₅)²_{mn} where F₅ = *F₅.
# The normalization is chosen such that you use the FULL F₅ but with
# coefficient 1/4 instead of 1/2.

# In practice, for the D3-brane solution:
# R_{mn} = (1/4·4!)F_{mp₁p₂p₃p₄}F_n^{p₁p₂p₃p₄}
# where F₅ is the FULL self-dual field strength.

# Let me instead just build F₅^m and compute FF^m + FF^e = FF^{SD}.

# Build F₅^m = *₁₂F₇ where F₇ = F₅^e ∧ dz₁ ∧ dz₂
# F₇ components at sorted indices:
print("\nBuilding F₇ = F₅^e ∧ dz₁ ∧ dz₂...")
F7 = FormField(rank=7, dim=D)
for idx5, val5 in F5e.nonzero_components.items():
    # idx5 has 5 indices. Wedge with z1 (idx=4) and z2 (idx=5)
    z1_idx = 4  # index of z1 in coords
    z2_idx = 5  # index of z2 in coords
    if z1_idx in idx5 or z2_idx in idx5:
        continue  # would make it zero
    idx7 = tuple(sorted(list(idx5) + [z1_idx, z2_idx]))
    # Sign from rearranging to sorted order
    original = list(idx5) + [z1_idx, z2_idx]
    perm_map = [idx7.index(x) for x in original]
    from sympy.combinatorics import Permutation
    sign = Permutation(perm_map).signature()
    F7[idx7] = sign * val5

# Now compute *₁₂F₇ = F₅^m
# *F₇_{m1...m5} = (1/7!) √{-g} ε_{m1...m5 n1...n7} g^{n1a1}...g^{n7a7} F₇_{a1...a7}
# For diagonal metric: simple

print("Computing *₁₂F₇ (Hodge dual)...")

# Levi-Civita symbol in 12d
def levi_civita_12(*indices):
    """Levi-Civita symbol for 12 indices."""
    indices = list(indices)
    if len(set(indices)) < 12:
        return 0
    # Count inversions
    sign = 1
    for i in range(12):
        for j in range(i+1, 12):
            if indices[i] > indices[j]:
                sign *= -1
    return sign

from itertools import combinations

# Determinant of metric (diagonal)
det_g = sp.Integer(1)
for i in range(D):
    det_g *= m.matrix[i, i]
det_g = cancel(det_g)
sqrt_neg_g = cancel(sp.sqrt(-det_g))  # √{-g}
print(f"√(-g) = {hf.substitute(cancel(sqrt_neg_g))}")

F5m = FormField(rank=5, dim=D)
for combo5 in combinations(range(D), 5):
    val = sp.Integer(0)
    # Sum over all 7-tuples that complement combo5
    complement = [i for i in range(D) if i not in combo5]
    # There's only one 7-tuple (the complement), up to permutation
    combo7 = tuple(complement)
    # Get ε_{combo5 + combo7}
    full_12 = list(combo5) + list(combo7)
    eps = levi_civita_12(*full_12)
    if eps == 0:
        continue
    # Get F₇ at this 7-tuple
    sorted7 = tuple(sorted(combo7))
    if sorted7 not in F7.nonzero_components:
        continue
    f7_val = F7[sorted7]
    if f7_val == 0:
        continue
    # Sign from sorting combo7
    perm7 = [sorted7.index(x) for x in combo7]
    sign7 = Permutation(perm7).signature()
    # Raise indices with g^{nn}
    ginv_prod = sp.Integer(1)
    for n in combo7:
        ginv_prod *= m.inv_matrix[n, n]
    val = eps * sign7 * sqrt_neg_g * ginv_prod * f7_val
    val = cancel(val)
    if val != 0:
        F5m[combo5] = val

print(f"F₅^m nonzero components: {len(F5m.nonzero_components)}")
for idx in sorted(F5m.nonzero_components.keys())[:4]:
    val = hf.substitute(cancel(F5m[idx]))
    names = [str(coords[i]) for i in idx]
    print(f"  F₅^m[{','.join(names)}] = {val}")

# Self-dual F₅ = F₅^e + F₅^m
F5sd = FormField(rank=5, dim=D)
for idx, val in F5e.nonzero_components.items():
    old = F5sd.nonzero_components.get(idx, sp.Integer(0))
    F5sd[idx] = old + val
for idx, val in F5m.nonzero_components.items():
    old = F5sd.nonzero_components.get(idx, sp.Integer(0))
    F5sd[idx] = old + val

print(f"\nF₅^SD nonzero components: {len(F5sd.nonzero_components)}")

# Compute FF for self-dual F₅
print("Computing FF for self-dual F₅...")
FFsd = form_contraction(F5sd, m)
normsd = form_norm_squared(F5sd, m)
normsd_val = hf.substitute(cancel(normsd))
print(f"|F₅^SD|² = {normsd_val} (should be 0)")

# Also compute FF^m separately
print("Computing FF for magnetic F₅...")
FFm = form_contraction(F5m, m)
normm = form_norm_squared(F5m, m)
normm_val = hf.substitute(cancel(normm))
print(f"|F₅^m|² = {normm_val}")

# Test formulas
blocks = [(0,'t'), (1,'x1'), (4,'z1'), (5,'z2'), (6,'y0'), (7,'y1')]

print("\n--- Test 1: ℛ = (1/4·4!)(F₅^SD)² (standard self-dual formula) ---")
for idx, name in blocks:
    ric = hf.substitute(cancel(Ric[idx, idx]))
    ffsd = hf.substitute(cancel(FFsd[idx, idx]))
    T = cancel(R(1,4) * ffsd / 24)
    res = cancel(ric - T)
    status = "✓" if res == 0 else f"✗ res={res}"
    print(f"  [{name:>4}]: {status}")

print("\n--- Test 2: ℛ = (1/2·4!)(F₅^e)² (electric only, coeff 1/2) ---")
for idx, name in blocks:
    ric = hf.substitute(cancel(Ric[idx, idx]))
    ffe = hf.substitute(cancel(FFe[idx, idx]))
    T = cancel(R(1,2) * ffe / 24)
    res = cancel(ric - T)
    status = "✓" if res == 0 else f"✗ res={res}"
    print(f"  [{name:>4}]: {status}")

print("\n--- Test 3: ℛ = (1/4·4!)(F₅^e)² (electric only, coeff 1/4) ---")
for idx, name in blocks:
    ric = hf.substitute(cancel(Ric[idx, idx]))
    ffe = hf.substitute(cancel(FFe[idx, idx]))
    T = cancel(R(1,4) * ffe / 24)
    res = cancel(ric - T)
    status = "✓" if res == 0 else f"✗ res={res}"
    print(f"  [{name:>4}]: {status}")

# Compare FF^e vs FF^m per block
print("\n--- FF comparison: electric vs magnetic ---")
for idx, name in blocks:
    ffe = hf.substitute(cancel(FFe[idx, idx]))
    ffm = hf.substitute(cancel(FFm[idx, idx]))
    ffsd = hf.substitute(cancel(FFsd[idx, idx]))
    print(f"  [{name:>4}]: FF^e={ffe}, FF^m={ffm}, FF^SD={ffsd}")

print("\n" + "="*60)
print("DONE")
print("="*60)
