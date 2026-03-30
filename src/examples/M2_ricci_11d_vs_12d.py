"""
Compare Ricci tensor of the M2 brane in 11d vs 12d.

11d:  ds² = H^{-2/3} ds²_{1,2} + H^{1/3} ds²_8

12d:  ds² = H^{-2/3} ds²_{1,2} + H^{1/3} ds²_8 + H^c dz²
      with c as a free parameter (c=0 is flat, c≠0 is warped)

H(r) harmonic in the 8d transverse space: H'' + 7/r H' = 0.
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import sympy as sp
from sugra import Metric, HarmonicFunction

# ── Warp factor exponents ────────────────────────────────────────────────────

a = sp.Rational(-2, 3)   # worldvolume warp
b = sp.Rational(1, 3)    # transverse warp
c = sp.Symbol('c')       # extra direction warp (free parameter)

# ── Coordinates ──────────────────────────────────────────────────────────────

t, x1, x2 = sp.symbols('t x1 x2', real=True)
y = list(sp.symbols('y0:8', real=True))          # 8 transverse
z = sp.Symbol('z', real=True)                     # extra flat direction

wv   = [t, x1, x2]       # 3 worldvolume
d_wv = len(wv)
d_tr = len(y)             # 8 transverse

# ── Harmonic function in 8d ──────────────────────────────────────────────────

hf = HarmonicFunction(transverse_coords=y)
H  = sp.Function('H')(hf.r_expr)

# ══════════════════════════════════════════════════════════════════════════════
#  11d metric:  ds² = H^a ds²_{1,2} + H^b ds²_8
# ══════════════════════════════════════════════════════════════════════════════

coords_11 = wv + y
D_11 = len(coords_11)     # 11

g11 = sp.zeros(D_11, D_11)
g11[0, 0] = -H**a                                # -H^{-2/3} dt²
for i in range(1, d_wv):
    g11[i, i] = H**a                             #  H^{-2/3} dx_i²
for i in range(d_tr):
    g11[d_wv + i, d_wv + i] = H**b               #  H^{1/3} dy_i²

metric_11 = Metric(g11, coords_11)

# ══════════════════════════════════════════════════════════════════════════════
#  12d metric:  ds² = H^a ds²_{1,2} + H^b ds²_8 + H^c dz²
# ══════════════════════════════════════════════════════════════════════════════

coords_12 = wv + y + [z]
D_12 = len(coords_12)     # 12

g12 = sp.zeros(D_12, D_12)
g12[0, 0] = -H**a
for i in range(1, d_wv):
    g12[i, i] = H**a
for i in range(d_tr):
    g12[d_wv + i, d_wv + i] = H**b
g12[D_12 - 1, D_12 - 1] = H**c                    # H^c dz²

metric_12 = Metric(g12, coords_12)

# ══════════════════════════════════════════════════════════════════════════════
#  Compute Ricci tensors
# ══════════════════════════════════════════════════════════════════════════════

print("Computing 11d Ricci tensor ...")
R11 = metric_11.ricci_tensor(simplify_func=sp.cancel)

print("Computing 12d Ricci tensor ...")
R12 = metric_12.ricci_tensor(simplify_func=sp.cancel)

# ══════════════════════════════════════════════════════════════════════════════
#  Display and compare
# ══════════════════════════════════════════════════════════════════════════════

def show(expr):
    """Substitute harmonic condition and simplify."""
    return sp.cancel(hf.substitute(sp.cancel(expr)))

def label(i, has_z):
    if i < d_wv:
        return "wv"
    elif i < d_wv + d_tr:
        return "tr"
    else:
        return "z"

print()
print("=" * 70)
print("  11d M2 brane:  ds² = H^(-2/3) ds²_{1,2} + H^(1/3) ds²_8")
print("=" * 70)

for i in range(D_11):
    name = str(coords_11[i])
    expr = show(R11[i, i])
    print(f"  R[{name},{name}]  ({label(i, False)})  =  {expr}")

print()
print("=" * 70)
print("  12d:  ds² = H^(-2/3) ds²_{1,2} + H^(1/3) ds²_8 + H^c dz²")
print("=" * 70)

for i in range(D_12):
    name = str(coords_12[i])
    expr = show(R12[i, i])
    print(f"  R[{name},{name}]  ({label(i, True)})  =  {expr}")

# ── Direct comparison ────────────────────────────────────────────────────────

print()
print("=" * 70)
print("  Comparison:  R_12d[i,i] - R_11d[i,i]  for shared directions")
print("=" * 70)

for i in range(D_11):
    name = str(coords_11[i])
    diff = show(R12[i, i] - R11[i, i])
    status = "  (same)" if diff == 0 else f"  DIFFERS: {diff}"
    print(f"  delta R[{name},{name}]  ({label(i, False)})  =  {diff}{status}")

print(f"\n  R[z,z]  (z)   =  {show(R12[D_12-1, D_12-1])}")

# ── Off-diagonal components ──────────────────────────────────────────────────

print()
print("=" * 70)
print("  Off-diagonal components R[i,j], i < j")
print("=" * 70)

for dim_label, R, coords, D in [("11d", R11, coords_11, D_11),
                                  ("12d", R12, coords_12, D_12)]:
    nonzero = []
    for i in range(D):
        for j in range(i + 1, D):
            val = show(R[i, j])
            if val != 0:
                nonzero.append((str(coords[i]), str(coords[j]), val))
    if nonzero:
        print(f"\n  {dim_label}: {len(nonzero)} nonzero off-diagonal component(s):")
        for ci, cj, val in nonzero:
            print(f"    R[{ci},{cj}] = {val}")
    else:
        print(f"\n  {dim_label}: all off-diagonal components vanish  ✓")
