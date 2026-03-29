"""Experiment 06: Sanity check with M2-brane in D=11.

The M2-brane solution:
  ds^2 = H^{-2/3} ds^2_{1,2} + H^{1/3} ds^2_8
  F4 = dC3, C_{012} = H^{-1}
  No dilaton

This should satisfy R_MN = T_MN exactly.
If it doesn't, there's a normalization issue in the code.
"""
import sys; sys.path.insert(0, 'src')
import sympy as sp
from sympy import Rational as R
from sugra import (HarmonicFunction, warped_product,
                   FormField, exterior_derivative, Verifier)

# M2-brane: D=11, p=2
# ds^2 = H^{-2/3}(-dt^2 + dx1^2 + dx2^2) + H^{1/3}(dy0^2+...+dy7^2)
# WV exponent: -(D-p-3)/(D-2) = -(11-2-3)/9 = -6/9 = -2/3 ✓
# Y exponent: (p+1)/(D-2) = 3/9 = 1/3 ✓

wv = list(sp.symbols('t x1 x2', real=True))
y = list(sp.symbols('y0:8', real=True))
coords = wv + y

hf = HarmonicFunction(transverse_coords=y)
H = sp.Function('H')(hf.r_expr)

metric = warped_product(
    warp_factors=[H**R(-2,3), H**R(1,3)],
    block_dims=[3, 8],
    block_signatures=['lorentzian', 'euclidean'],
    coordinates=coords,
)

C = FormField(rank=3, dim=11)
C[(0,1,2)] = 1/H
F = exterior_derivative(C, coords)

print("M2-brane in D=11:")
print(f"F4 components: {len(F.nonzero_components)}")

soln = dict(metric=metric, F=F, Phi=0, alpha=0, coords=coords, hf=hf)
result = Verifier(soln).check()
print(f"\nVerifier result: {'PASS' if result else 'FAIL'}")

# Now try M5-brane: D=11, p=5
# ds^2 = H^{-1/3}(-dt^2+dx1^2+...+dx5^2) + H^{2/3}(dy0^2+...+dy4^2)
# F4 = *dC6 → F4 is the magnetic dual, or equivalently dA3 from transverse
# Actually F7 = dA6 with A_{012345} = H^{-1}, and F4 = *F7

print("\n" + "="*60)
print("M5-brane in D=11 (using F7 = dA6):")
print("="*60)

wv5 = list(sp.symbols('t x1 x2 x3 x4 x5', real=True))
y5 = list(sp.symbols('y0:5', real=True))
coords5 = wv5 + y5

hf5 = HarmonicFunction(transverse_coords=y5)
H5 = sp.Function('H')(hf5.r_expr)

metric5 = warped_product(
    warp_factors=[H5**R(-1,3), H5**R(2,3)],
    block_dims=[6, 5],
    block_signatures=['lorentzian', 'euclidean'],
    coordinates=coords5,
)

A6 = FormField(rank=6, dim=11)
A6[(0,1,2,3,4,5)] = 1/H5
F7 = exterior_derivative(A6, coords5)

print(f"F7 components: {len(F7.nonzero_components)}")

soln5 = dict(metric=metric5, F=F7, Phi=0, alpha=0, coords=coords5, hf=hf5)
result5 = Verifier(soln5).check()
print(f"\nVerifier result: {'PASS' if result5 else 'FAIL'}")
