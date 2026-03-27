"""
D3 brane solution verification using the sugra package.

The D3 brane in Type IIB supergravity (Einstein frame) is:

  ds^2 = H^{-1/2}(-dt^2 + dx1^2 + dx2^2 + dx3^2)
       + H^{1/2}( dy0^2 + dy1^2 + dy2^2 + dy3^2 + dy4^2 + dy5^2 )

  phi = 0    (constant dilaton)

  F5 = dC4 + *dC4,   where C4 = H^{-1} dt ^ dx1 ^ dx2 ^ dx3
  F5 is self-dual:   F5 = *F5

H = H(r) is harmonic in the 6 transverse directions:
  H'' + (5/r) H' = 0   =>   H = 1 + L^4 / r^4

Coordinate index map (0-indexed):
  0 = t,  1 = x1,  2 = x2,  3 = x3,   (worldvolume)
  4 = y0, 5 = y1, 6 = y2, 7 = y3, 8 = y4, 9 = y5   (transverse)
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import sympy as sp
from sugra import (
    HarmonicFunction, warped_product,
    FormField, FluxField, ScalarField, Solution,
    hodge_star,
)

# ── Coordinates ────────────────────────────────────────────────────────────────

t, x1, x2, x3 = sp.symbols('t x1 x2 x3', real=True)

# HarmonicFunction manages the 6 transverse coords y0,...,y5
# and the radial symbol r = sqrt(y0^2 + ... + y5^2)
hf = HarmonicFunction(transverse_dim=6)
y  = hf.transverse_coords       # [y0, y1, y2, y3, y4, y5]

coords = [t, x1, x2, x3] + y   # full 10d coordinate list (length 10)

# H expressed as a SymPy function so SymPy can differentiate it w.r.t. y_a
H_func = sp.Function('H')(hf.r_expr)

# ── Metric ─────────────────────────────────────────────────────────────────────
#
# ds^2 = H^{-1/2} eta_{mu nu} dx^mu dx^nu  +  H^{1/2} delta_{ab} dy^a dy^b
#
# warped_product builds the diagonal metric block by block.
# Block 0: 4d worldvolume, Lorentzian, warp factor H^{-1/2}
# Block 1: 6d transverse,  Euclidean,  warp factor H^{+1/2}

metric = warped_product(
    warp_factors     = [sp.Rational(-1, 2), sp.Rational(1, 2)],
    block_dims       = [4, 6],
    block_signatures = ['lorentzian', 'euclidean'],
    coordinates      = coords,
    H                = H_func,
)

# ── Dilaton ────────────────────────────────────────────────────────────────────
#
# The D3 brane does not source the dilaton: phi = 0.
# We still register it so check_dilaton() runs (trivially Box(0) = 0).

dilaton = ScalarField(expr=sp.S(0), name='dilaton')

# ── RR 5-form F5 ───────────────────────────────────────────────────────────────
#
# The RR 4-form potential sourced by the D3 brane is:
#   C4 = H^{-1} dt ^ dx1 ^ dx2 ^ dx3
#
# Its field strength (worldvolume piece):
#   F5^{wv} = dC4 = d(H^{-1}) ^ dt ^ dx1 ^ dx2 ^ dx3
#
# Component by component:
#   F5^{wv}_{t,x1,x2,x3,ya} = d(H^{-1})/d(ya) = -H^{-2} H'(r) ya/r
#
# The full self-dual F5 is constructed by adding the Hodge dual piece:
#   F5 = F5^{wv} + *F5^{wv}
# which ensures F5 = *F5 by construction (using *(*F5^{wv}) = F5^{wv} in 10d Lorentzian).

H_inv = H_func ** (-1)

F5_wv = FormField(rank=5, dim=10)
for a, ya in enumerate(y):
    comp = sp.diff(H_inv, ya)   # = -H^{-2} H'(r) ya/r
    if comp != 0:
        # Sorted index: (t=0, x1=1, x2=2, x3=3, ya=4+a)
        F5_wv[(0, 1, 2, 3, 4 + a)] = comp

# Add the Hodge dual to get the full self-dual F5
# hodge_star uses signature=-1 (Lorentzian) by default
F5_dual = hodge_star(F5_wv, metric)
F5_full = F5_wv + F5_dual

flux = FluxField(form=F5_full, dilaton_coupling=0, self_dual=True, name='F5')

# ── Verification ───────────────────────────────────────────────────────────────

sol = Solution(
    metric      = metric,
    fields      = [dilaton, flux],
    coordinates = coords,
    harmonic    = hf,
)

print("=" * 60)
print("D3 brane — Type IIB supergravity (Einstein frame)")
print("Coordinates: (t, x1, x2, x3 | y0, y1, y2, y3, y4, y5)")
print("Equations of motion from action:")
print("  S = ∫ d^10x sqrt(-g) [R - 1/2(∂φ)² - 1/2 e^{αφ}|F5|²]")
print("=" * 60)
print()

# ── Inspect both sides of Einstein equation ───────────────────────────────────
# R_{MN}  =  (1/2)( ∂_M φ ∂_N φ  +  e^{αφ} S_{MN}|_{F5} )
#
# For a diagonal metric the only independent components are the diagonal ones.
# We apply hf.substitute() to trade H(y_a) derivatives for the symbols H, H', r.

# Compute R and T raw (no simplify_func), then substitute H(y_a) → H, H', r,
# and only then call cancel.  Calling cancel *before* substitution would ask
# SymPy to simplify expressions containing H(sqrt(y0²+...+y5²)), which it
# cannot do efficiently.
R = sol.ricci_tensor()
T = sol.stress_energy_tensor()

print("Left-hand side  R_MN  (after harmonic substitution):")
for i in range(10):
    expr = sp.cancel(hf.substitute(R[i, i]))
    print(f"  R[{i},{i}] = {expr}")

print()
print("Right-hand side  T_MN = (1/2)(∂φ∂φ + e^αφ S_MN|_F)  (after harmonic substitution):")
for i in range(10):
    expr = sp.cancel(hf.substitute(T[i, i]))
    print(f"  T[{i},{i}] = {expr}")

print()

# ── Run all checks ────────────────────────────────────────────────────────────
# Do not pass simplify_func here — the harmonic substitution inside each check
# already reduces residuals to H, H', r form before cancel is attempted.
report = sol.verify_all()
report.summary()
