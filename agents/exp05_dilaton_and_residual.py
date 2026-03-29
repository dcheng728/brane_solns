"""Experiment 05: Check with dilaton and analyze residual structure.

Key idea: The 12d metric encodes the dilaton via g_{z1}/g_{z2} = H.
The F1-string has e^Phi = H^{-1/2}, alpha = -1 for NS-NS 3-form.
In 12d, maybe the form_stress_energy needs the dilaton coupling factor.

Also: try the Verifier directly, and check if scaling the ENTIRE
stress-energy formula by a constant fixes things.
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

# --- D1 4-form ---
C_D1 = FormField(rank=3, dim=12)
C_D1[(0,1,3)] = 1/H
F_D1 = exterior_derivative(C_D1, coords)

# --- Run Verifier with D1 form, no dilaton ---
print("=" * 60)
print("Test A: D1 4-form, no dilaton")
print("=" * 60)
soln = dict(metric=metric, F=F_D1, Phi=0, alpha=0, coords=coords, hf=hf)
Verifier(soln).check()

# --- Test B: D1 form WITH dilaton Phi = -(1/2)ln(H), alpha = -1 ---
# (This is the F1-string dilaton in type IIB)
print("\n" + "=" * 60)
print("Test B: D1 4-form + dilaton Phi=-(1/2)lnH, alpha=-1")
print("=" * 60)
Phi = -sp.Rational(1,2) * sp.log(H)
soln_dil = dict(metric=metric, F=F_D1, Phi=Phi, alpha=-1, coords=coords, hf=hf)
Verifier(soln_dil).check()

# --- Test C: F1 form WITH dilaton ---
print("\n" + "=" * 60)
print("Test C: F1 4-form (z1 leg) + dilaton Phi=-(1/2)lnH, alpha=-1")
print("=" * 60)
C_F1 = FormField(rank=3, dim=12)
C_F1[(0,1,2)] = 1/H
F_F1 = exterior_derivative(C_F1, coords)
soln_f1 = dict(metric=metric, F=F_F1, Phi=Phi, alpha=-1, coords=coords, hf=hf)
Verifier(soln_f1).check()

# --- Test D: Check what T/R ratio the D1 form gives ---
print("\n" + "=" * 60)
print("Test D: D1 form T/R ratios and residual analysis")
print("=" * 60)
Ric = metric.ricci_tensor(simplify_func=sp.cancel)
T_D1 = form_stress_energy(F_D1, metric)

for i, name in [(0,'t'), (2,'z1'), (3,'z2'), (4,'y0')]:
    r = hf.substitute(sp.cancel(Ric[i,i]))
    t = hf.substitute(sp.cancel(T_D1[i,i]))
    ratio = sp.cancel(t/r)
    # Evaluate at a specific point
    subs_vals = [(hf.H, 2), (hf.Hp, 1), (hf.r, 1), (sp.Symbol('y0',real=True), sp.Rational(1,3))]
    ratio_num = ratio.subs(subs_vals)
    print(f"  T/R [{name}] = {ratio}  (at test point: {ratio_num})")

# --- Test E: Try scaling T by a constant kappa ---
print("\n" + "=" * 60)
print("Test E: What constant kappa would make kappa*T_D1 = R?")
print("=" * 60)
for i, name in [(0,'t'), (2,'z1'), (3,'z2')]:
    r = hf.substitute(sp.cancel(Ric[i,i]))
    t = hf.substitute(sp.cancel(T_D1[i,i]))
    kappa = sp.cancel(r/t)
    print(f"  kappa[{name}] = {kappa}")

# y0 block: check aniso part
r_y0 = hf.substitute(sp.cancel(Ric[4,4]))
t_y0 = hf.substitute(sp.cancel(T_D1[4,4]))
# At y0=0
r0 = r_y0.subs(sp.Symbol('y0',real=True), 0)
t0 = t_y0.subs(sp.Symbol('y0',real=True), 0)
kappa_iso = sp.cancel(r0/t0)
print(f"  kappa[y0_iso] = {kappa_iso}")
# Aniso: coefficient of y0^2
y0s = sp.Symbol('y0', real=True)
r_aniso = sp.cancel(sp.diff(r_y0, y0s, 2) / 2)
t_aniso = sp.cancel(sp.diff(t_y0, y0s, 2) / 2)
kappa_aniso = sp.cancel(r_aniso/t_aniso)
print(f"  kappa[y0_aniso] = {kappa_aniso}")

# --- Test F: Try F1 + dilaton alpha=+1 ---
print("\n" + "=" * 60)
print("Test F: F1 4-form + dilaton Phi=-(1/2)lnH, alpha=+1")
print("=" * 60)
soln_f1p = dict(metric=metric, F=F_F1, Phi=Phi, alpha=1, coords=coords, hf=hf)
Verifier(soln_f1p).check()

# --- Test G: D1 + dilaton alpha=+1 ---
print("\n" + "=" * 60)
print("Test G: D1 4-form + dilaton Phi=-(1/2)lnH, alpha=+1")
print("=" * 60)
soln_d1p = dict(metric=metric, F=F_D1, Phi=Phi, alpha=1, coords=coords, hf=hf)
Verifier(soln_d1p).check()
