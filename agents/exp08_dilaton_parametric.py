"""Experiment 08: Parametric dilaton scan.

Exp05 only tested specific dilaton values (Phi = -1/2 ln H, alpha = ±1).
Here we solve the FULL parametric system with:
  Phi = a * ln(H)     (arbitrary power a)
  alpha = symbolic     (arbitrary coupling)
  D1 form amplitude = symbolic c_d
  F1 form amplitude = symbolic c_f  (optional)

The stress-energy is:
  T_{MN} = (1/2) dPhi_M dPhi_N + (1/2) e^{alpha*Phi} [FF_{MN} - (n-1)/(D-2) |F²| g_{MN}]

The dilaton kinetic term contributes ONLY to blocks that contain
transverse directions (since dH is transverse).
"""
import sys; sys.path.insert(0, 'src')
import sympy as sp
from sympy import Rational as R
from sugra import (HarmonicFunction, warped_product,
                   FormField, exterior_derivative,
                   form_stress_energy, form_contraction,
                   form_norm_squared, Verifier)

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
print("Done.\n")

# Collect Ricci components
R_blocks = {}
for i, name in [(0,'t'), (2,'z1'), (3,'z2'), (4,'y0')]:
    R_blocks[name] = hf.substitute(sp.cancel(Ric[i,i]))

print("Ricci tensor blocks:")
for name, val in R_blocks.items():
    print(f"  R[{name}] = {val}")

# ===================================================================
# Part A: D1 form with parametric dilaton Phi = a*ln(H), coupling alpha
# ===================================================================
print("\n" + "="*60)
print("Part A: D1 form + Phi = a*ln(H), arbitrary alpha")
print("="*60)

# D1 form
C_D1 = FormField(rank=3, dim=12)
C_D1[(0,1,3)] = 1/H
F_D1 = exterior_derivative(C_D1, coords)
T_D1_form = form_stress_energy(F_D1, metric)

# D1 form contraction and norm
FF_D1 = form_contraction(F_D1, metric)
Fsq_D1 = form_norm_squared(F_D1, metric)

# Compute T_D1 blocks (form part only, no dilaton)
T_D1_blocks = {}
for i, name in [(0,'t'), (2,'z1'), (3,'z2'), (4,'y0')]:
    T_D1_blocks[name] = hf.substitute(sp.cancel(T_D1_form[i,i]))

print("D1 stress-energy blocks (form part only, no dilaton, no exp coupling):")
for name, val in T_D1_blocks.items():
    print(f"  T_D1[{name}] = {val}")

# Dilaton Phi = a * ln(H)
# dPhi_M = a * (1/H) * dH_M = a * H'/H * (y_k/r) for transverse directions
# dPhi_M = 0 for worldvolume/z1/z2 directions
#
# So: (1/2) dPhi_M dPhi_N contributes ONLY to transverse block:
#   T^dil_{y_i, y_j} = (1/2) a^2 * (H'/H)^2 * (y_i * y_j / r^2)
#   T^dil_{others} = 0
#
# After substitution:
#   T^dil_{y0,y0} = (1/2) a^2 * H'^2 / H^2 * y0^2/r^2
#                  = (1/2) a^2 * y0^2/r^2 * H'^2/H^2
#
# But we need to account for the metric on the transverse block: g_{y0,y0} = H^{1/4}
# So dPhi_{y0} = a * (H'/(H)) * (y0/r) ... wait, this is the COORDINATE derivative.
# The stress-energy uses dPhi with indices down:
#   dPhi_M = partial_M Phi = a * partial_M ln(H)
# For M = y_k: partial_{yk} H = H' * y_k/r
# So dPhi_{yk} = a * H'/(H) * y_k/r
#
# T^dil_{yk, yl} = (1/2) * a^2 * H'^2/H^2 * y_k*y_l / r^2
# T^dil_{t,t} = T^dil_{z1} = T^dil_{z2} = 0

a, alpha = sp.symbols('a alpha', real=True)
y0 = sp.Symbol('y0', real=True)

# Dilaton kinetic contribution to y0,y0
T_dil_y0 = R(1,2) * a**2 * hf.Hp**2 / hf.H**2 * y0**2 / hf.r**2

print(f"\nDilaton kinetic T^dil[y0] = {T_dil_y0}")

# The full stress-energy with dilaton coupling:
# T_{MN} = T^dil_{MN} + (1/2) e^{alpha*Phi} [FF_{MN}/(p-1)! - (p-1)/(D-2)|F²|g_{MN}]
# For D1 4-form: p=4, D=12
# factor = (1/2) * e^{alpha*a*ln(H)} = (1/2) * H^{alpha*a}
#
# Actually form_stress_energy already computes the undilaton-coupled version:
# T^form_{MN} = (1/2) [FF_{MN}/(p-1)! - (p-1)/(D-2)|F²|g_{MN}]
# So the full T is: T^dil + H^{alpha*a} * T^form
#
# Wait, that's not right. Let me re-examine.
# The standard formula (hep-th/9701088):
# T_{MN} = (1/2)(dPhi)_M(dPhi)_N + (1/2) e^{alpha Phi} [1/(p!) F_{MA...} F_N^{A...} - (p-1)/((D-2)*p!) |F|^2 g_{MN}]
#
# Our form_stress_energy computes the non-dilaton part:
# T^form = (1/2)[FF/(p-1)! - (p-1)/(D-2) |F|^2 g]
# This is equal to the second term when alpha=0 (no dilaton coupling).
# With dilaton: the second term gets multiplied by exp(alpha * Phi) = H^{alpha*a}

# For each block, the equation R = T^dil + H^{a*alpha} * T^form becomes:
# R[block] = T^dil[block] + H^{a*alpha} * T^form[block]

# For t, z1, z2 blocks: T^dil = 0, so
# R[block] = H^{a*alpha} * T^form[block]

# This means: R[t]/T^form[t] = R[z1]/T^form[z1] = R[z2]/T^form[z2] = H^{a*alpha}
# We already know these ratios (kappa values from exp05):
# kappa[t] = 15/14, kappa[z1] = 5/3, kappa[z2] = 5/7
# These are DIFFERENT NUMBERS, not just H-power-different.
# For them to equal H^{a*alpha}, they'd need to be the same number.
# Unless the H-powers also differ...

# Actually wait, let me be more careful. The T^form values already have specific
# H-power dependence. The R values also have H-power dependence.
# The ratio R/T should be H-independent (a pure number) for block t, z1, z2
# because both R and T have the same H-power (that was verified in exp03).
# So H^{a*alpha} would just be a number multiplied by the ratio.

# Let's compute R/T^form for each block explicitly:
print("\n--- Checking R/T^form ratios ---")
for name in ['t', 'z1', 'z2']:
    ratio = sp.cancel(R_blocks[name] / T_D1_blocks[name])
    print(f"  R/T_D1[{name}] = {ratio}")

# For the y0 block: R = T^dil + H^{a*alpha} * T^form
# R[y0] = (1/2) a^2 H'^2/H^2 y0^2/r^2 + H^{a*alpha} * T_D1[y0]

# Since R_t/T_t ≠ R_z1/T_z1 ≠ R_z2/T_z2, we CANNOT satisfy
# R = H^{a*alpha} * T for all non-dilaton blocks simultaneously.
# The dilaton exponential is a SINGLE constant, can't fix different ratios.

print("\n*** CONCLUSION: Exponential dilaton coupling H^{a*alpha} is a single")
print("    constant multiplier on T^form. Since R/T^form ratios differ across")
print("    blocks (15/14, 5/3, 5/7), no single H^{a*alpha} can match all.")
print("    The dilaton coupling CANNOT resolve the mismatch.")

# ===================================================================
# Part B: But wait — what if a*alpha is chosen so H^{a*alpha} has a
# NONTRIVIAL H-dependence that shifts the H-powers?
# This requires that a*alpha contributes DIFFERENT effective H-powers
# to different blocks. But H^{a*alpha} is the same factor everywhere.
# The ONLY way it helps is if T^form has DIFFERENT H-powers per block
# and we want to shift them to match R.
# Since exp03 showed D1's H-powers ALREADY match R, H^{a*alpha} would
# BREAK the H-power matching unless a*alpha = 0.
# ===================================================================

print("\n" + "="*60)
print("Part B: Verify that a*alpha ≠ 0 breaks H-power matching")
print("="*60)

# R has H-powers: t→H^{-3}, z1→H^{-7/4}, z2→H^{-11/4}, y0→H^{-2}
# T_D1 has same H-powers (verified in exp03)
# H^{a*alpha} * T_D1 has H-powers shifted by a*alpha
# Match requires a*alpha = 0

print("T_D1 H-powers match R H-powers already.")
print("H^{a*alpha} shifts all by a*alpha → match requires a*alpha=0")
print("With a*alpha=0: no dilaton coupling effect at all.")
print("→ Single-form + dilaton CANNOT work for D1.")

# ===================================================================
# Part C: What about TWO forms with dilaton, where different couplings
# allow different H-power shifts?
# In IIB SUGRA: B_2 couples with e^{-Phi} and C_2 with e^{+Phi}
# So T = (1/2)(dPhi)^2 + (1/2)e^{-Phi}|H_3|^2 + (1/2)e^{+Phi}|F_3|^2
# In 12d: F_D1 = dC → couples with e^{+alpha1*Phi}
#          F_F1 = dB → couples with e^{-alpha2*Phi}
# If alpha1 ≠ alpha2, they have DIFFERENT H-power shifts!
# ===================================================================

print("\n" + "="*60)
print("Part C: Two forms with DIFFERENT dilaton couplings")
print("In IIB: B_2 couples as e^{-Phi}, C_2 couples as e^{+Phi}")
print("="*60)

# With Phi = a*ln(H):
# e^{+Phi} = H^a, e^{-Phi} = H^{-a}
# T = (1/2)a^2 (dln H)^2 + c_d^2 * H^a * T_D1 + c_f^2 * H^{-a} * T_F1
#
# For this to have correct H-powers:
# T_D1 has H-power h_D in each block
# H^a * T_D1 has H-power (h_D + a) → must match R's H-power h_R
# Similarly H^{-a} * T_F1 has H-power (h_F - a) → must also match h_R
#
# So: h_D + a = h_R AND h_F - a = h_R (for each block, or at least most)
# → a = h_R - h_D and a = h_F - h_R
# → h_D + h_F = 2*h_R

# Let me find which F1 power p_F gives h_F such that h_D + h_F = 2*h_R

# From exp03: D1 (C_{t,x1,z2}=H^{-1}) matches all H-powers
# So h_D = h_R already, meaning a = 0.
# If we want a ≠ 0, we need a D1-type form with SHIFTED H-powers.

# Alternative: Use F1 form (C_{t,x1,z1}) with various p values
# and find if h_D1 + a = h_R and h_F1 - a = h_R can be SIMULTANEOUSLY satisfied
# across all blocks.

# F1 at p=-1/2 matches H-powers (from exp03, form #2)
# Let's check: does h_F1 + h_D1 = 2*h_R for each block?

C_F1 = FormField(rank=3, dim=12)
C_F1[(0,1,2)] = H**R(-1,2)
F_F1 = exterior_derivative(C_F1, coords)
T_F1_form = form_stress_energy(F_F1, metric)

T_F1_blocks = {}
for i, name in [(0,'t'), (2,'z1'), (3,'z2'), (4,'y0')]:
    T_F1_blocks[name] = hf.substitute(sp.cancel(T_F1_form[i,i]))

print("\nF1(p=-1/2) stress-energy blocks:")
for name, val in T_F1_blocks.items():
    print(f"  T_F1[{name}] = {val}")

# Check H-powers in T_D1 and T_F1
# T_D1[t] ∝ H'²/H^3 → H-power = -3 (same as R[t])
# T_F1[t] at p=-1/2: let's extract
print("\nH-power analysis for two-form system:")
print("(Checking if h_D1 + h_F1 = 2*h_R for each block)")

# Extract H-powers by looking at the structure
# T ~ coeff * H'^2 * H^{power} * (y-dependent part)
# For the t-block at y0=0: T should be proportional to H'^2 * H^{power}

# D1 and F1 both match R H-powers (from exp03), so h_D1 = h_F1 = h_R
# → h_D1 + h_F1 = 2*h_R trivially.
# Then a = h_R - h_D1 = 0 again.

print("Both D1 and F1(p=-1/2) individually match R's H-powers.")
print("So h_D + h_F = 2*h_R is trivially satisfied with a=0.")
print("No gain from dual-coupling dilaton when BOTH forms already match H-powers.")

# ===================================================================
# Part D: What if we DON'T require H-power matching and let the dilaton
# coupling SHIFT the H-powers to match?
# Use D1 with some "wrong" H-power and let H^{a*alpha} fix it.
# But D1's H-power structure comes from its index structure and the
# metric warp factors — NOT from the potential power p.
# Actually the potential power p DOES affect the H-power.
# Let me check: what H-power does D1 at p=-3/4 give?
# ===================================================================

print("\n" + "="*60)
print("Part D: D1 at non-standard power + dilaton shift")
print("="*60)

for p_val in [R(-3,4), R(-1,2), R(-5,4), R(-3,2), R(-2)]:
    C_test = FormField(rank=3, dim=12)
    C_test[(0,1,3)] = H**p_val
    F_test = exterior_derivative(C_test, coords)
    T_test = form_stress_energy(F_test, metric)

    # Check H-power in t-block
    t_val = hf.substitute(sp.cancel(T_test[0,0]))
    r_val = R_blocks['t']

    if t_val == 0:
        print(f"  p={p_val}: T[t]=0, skip")
        continue

    ratio = sp.cancel(t_val / r_val)
    # This ratio should be H^{something} * constant if H-powers differ
    print(f"  D1(p={p_val}): T/R[t] = {ratio}")

    # Check all blocks
    for i, name in [(2,'z1'), (3,'z2')]:
        t_b = hf.substitute(sp.cancel(T_test[i,i]))
        r_b = R_blocks[name]
        if t_b != 0:
            ratio_b = sp.cancel(t_b / r_b)
            print(f"                T/R[{name}] = {ratio_b}")

# ===================================================================
# Part E: Systematic: For D1 + F1 + dilaton system,
# T = (1/2)a^2 (dlnH)^2 + c_d^2 H^{alpha_d*a} T^D1 + c_f^2 H^{alpha_f*a} T^F1
# with alpha_d = +1, alpha_f = -1 (IIB convention)
# Phi = a * ln(H)
#
# Write equations for each block and solve for (a, c_d, c_f).
# ===================================================================

print("\n" + "="*60)
print("Part E: Full IIB-style system with opposite dilaton couplings")
print("T = (1/2)a²(dlnH)² + c_d² H^a T_D1 + c_f² H^{-a} T_F1 = R")
print("="*60)

# Symbolic parameters
a_sym = sp.Symbol('a')
cd2 = sp.Symbol('cd2')  # c_d^2
cf2 = sp.Symbol('cf2')  # c_f^2

# For the t, z1, z2 blocks (no dilaton kinetic term):
# cd2 * H^a * T_D1[block] + cf2 * H^{-a} * T_F1[block] = R[block]

# Since T_D1, T_F1, and R all have the form coeff * H'^2 * H^{power},
# we need the H-power to match. T_D1 already matches R. So:
# H^a * H^{h_D1} = H^{h_R} → a = 0 (since h_D1 = h_R)
# H^{-a} * H^{h_F1} = H^{h_R} → a = h_F1 - h_R

# If h_F1 = h_R (as for F1 at p=-1/2), then a=0 from both.
# If h_F1 ≠ h_R, the first equation forces a=0, which means
# the F1 term has wrong H-powers. So only a=0 works.

# With a=0: back to exp04's result (no dilaton coupling).
# → The dilaton CANNOT help because D1 already has correct H-powers.

# BUT WAIT: What if neither form individually matches H-powers,
# but the sum of H^a*T_D1 + H^{-a}*T_F1 does?
# This requires cancellations between terms with different H-powers.
# E.g., T[t] = cd2 * A * H'^2/H^{3-a} + cf2 * B * H'^2/H^{3+a}
# For this to equal R[t] = (3/8) H'^2/H^3, we need either:
# (a) Both terms have H^{-3}: a=0
# (b) One term cancels, leaving only one: but both cd2, cf2 > 0

# Actually, if a≠0, the two terms have DIFFERENT H-powers.
# Their sum can only equal a single H-power term if one coefficient is zero.
# So either cd2=0 or cf2=0, reducing to single-form case.

print("Analysis: With D1 matching R's H-powers, the equation")
print("  cd2 * H^a * T_D1 + cf2 * H^{-a} * T_F1 = R")
print("splits by H-power into independent equations for each H-power.")
print("Since T_D1 has the same H-power as R:")
print("  - The H^{h_R} part: cd2 * H^a contributes only if a=0")
print("  - The H^{h_F1-a} part: cf2 term, matches h_R only if h_F1=h_R and a=0")
print("→ Two-form + dilaton reduces to a=0 (no coupling) case.")

# ===================================================================
# Part F: What about forms that DON'T individually match H-powers?
# Could two "wrong" forms with opposite dilaton couplings conspire?
# Try D1(p=p1) + F1(p=p2) + Phi=a*ln(H) with alpha_D=+1, alpha_F=-1
# ===================================================================

print("\n" + "="*60)
print("Part F: Search for D1(p1)+F1(p2)+dilaton with non-standard powers")
print("="*60)

# For D1(p=p1): T has H-power that depends on p1 and warp factors
# For F1(p=p2): T has H-power that depends on p2 and warp factors
# With dilaton shifts: effective powers are h_D1 + a and h_F1 - a

# Build T for several (p1, p2) combinations and check if there exists
# a value of 'a' such that BOTH h_D1+a = h_R AND h_F1-a = h_R.
# This means h_D1 + h_F1 = 2*h_R for EVERY block simultaneously.

# Let me compute h_D1(p) and h_F1(p) for each block by testing T ratios

def get_h_power(T_block, R_block):
    """Extract the H-power difference between T and R.
    If T = c * H'^2 * H^h and R = c' * H'^2 * H^{h_R},
    then T/R = (c/c') * H^{h-h_R}.
    Returns the exponent of H in the ratio, or None if zero."""
    ratio = sp.cancel(T_block / R_block)
    if ratio == 0:
        return None
    # Check if ratio is a pure number (no H dependence)
    if not ratio.has(hf.H):
        return 0  # Same H-power
    # Try to extract: ratio should be const * H^k
    # Substitute H=2 and H=3 to find k
    r2 = float(ratio.subs([(hf.H, 2), (hf.Hp, 1), (hf.r, 1), (sp.Symbol('y0',real=True), 0)]))
    r3 = float(ratio.subs([(hf.H, 3), (hf.Hp, 1), (hf.r, 1), (sp.Symbol('y0',real=True), 0)]))
    if abs(r2) < 1e-15 or abs(r3) < 1e-15:
        return None
    import math
    k = math.log(r3/r2) / math.log(3/2)
    return round(k*4)/4  # Round to nearest quarter

print("\nH-power shifts for D1(p) and F1(p) relative to R:")
print(f"{'p':>6} | {'D1 shift t':>10} {'z1':>6} {'z2':>6} | {'F1 shift t':>10} {'z1':>6} {'z2':>6}")
print("-"*65)

for p_num, p_den in [(-2,1),(-7,4),(-3,2),(-5,4),(-1,1),(-3,4),(-1,2),(-1,4)]:
    p_val = R(p_num, p_den)

    # D1
    C_d = FormField(rank=3, dim=12)
    C_d[(0,1,3)] = H**p_val
    F_d = exterior_derivative(C_d, coords)
    T_d = form_stress_energy(F_d, metric)

    # F1
    C_f = FormField(rank=3, dim=12)
    C_f[(0,1,2)] = H**p_val
    F_f = exterior_derivative(C_f, coords)
    T_f = form_stress_energy(F_f, metric)

    shifts_d = []
    shifts_f = []
    for i, name in [(0,'t'), (2,'z1'), (3,'z2')]:
        td = hf.substitute(sp.cancel(T_d[i,i]))
        tf = hf.substitute(sp.cancel(T_f[i,i]))
        rb = R_blocks[name]
        sd = get_h_power(td, rb) if td != 0 else 'X'
        sf = get_h_power(tf, rb) if tf != 0 else 'X'
        shifts_d.append(sd)
        shifts_f.append(sf)

    print(f"{p_val:>6} | {shifts_d[0]:>10} {shifts_d[1]:>6} {shifts_d[2]:>6} | {shifts_f[0]:>10} {shifts_f[1]:>6} {shifts_f[2]:>6}")

    # Check if h_D + h_F = 0 (which means h_D + a = 0 and h_F - a = 0 → a = h_F = -h_D)
    # i.e., shifts are equal and opposite
    if all(isinstance(s, (int, float)) for s in shifts_d + shifts_f):
        pairs = [(shifts_d[j] + shifts_f[j]) for j in range(3)]
        if all(abs(p) < 0.01 for p in pairs):
            print(f"    *** D1({p_val})+F1({p_val}) with a=-shift: H-powers match!")

# ===================================================================
# Part G: Try D1(p1) + F1(p2) with DIFFERENT powers
# ===================================================================
print("\n" + "="*60)
print("Part G: D1(p1) + F1(p2) with different powers")
print("="*60)

powers = [R(-2,1), R(-7,4), R(-3,2), R(-5,4), R(-1,1), R(-3,4), R(-1,2), R(-1,4)]

found = False
for p1 in powers:
    C_d = FormField(rank=3, dim=12)
    C_d[(0,1,3)] = H**p1
    F_d = exterior_derivative(C_d, coords)
    T_d = form_stress_energy(F_d, metric)

    shifts_d = []
    for i in [0, 2, 3]:
        td = hf.substitute(sp.cancel(T_d[i,i]))
        if td == 0:
            shifts_d = None
            break
        shifts_d.append(get_h_power(td, R_blocks[['t','','z1','z2'][i] if i < 4 else '']))

    if shifts_d is None:
        continue
    # Fix: properly index
    d_shifts = {}
    for idx, (i, name) in enumerate([(0,'t'), (2,'z1'), (3,'z2')]):
        td = hf.substitute(sp.cancel(T_d[i,i]))
        if td == 0:
            d_shifts = None
            break
        d_shifts[name] = get_h_power(td, R_blocks[name])
    if d_shifts is None:
        continue

    for p2 in powers:
        C_f = FormField(rank=3, dim=12)
        C_f[(0,1,2)] = H**p2
        F_f = exterior_derivative(C_f, coords)
        T_f = form_stress_energy(F_f, metric)

        f_shifts = {}
        for i, name in [(0,'t'), (2,'z1'), (3,'z2')]:
            tf = hf.substitute(sp.cancel(T_f[i,i]))
            if tf == 0:
                f_shifts = None
                break
            f_shifts[name] = get_h_power(tf, R_blocks[name])
        if f_shifts is None:
            continue

        # Check: d_shift + f_shift = same value for all blocks
        # (This would allow a single 'a' to fix both)
        sums = [d_shifts[n] + f_shifts[n] for n in ['t','z1','z2']
                if isinstance(d_shifts[n], (int,float)) and isinstance(f_shifts[n], (int,float))]
        if len(sums) == 3 and abs(sums[0] - sums[1]) < 0.01 and abs(sums[1] - sums[2]) < 0.01:
            a_needed = sums[0]/2 if sums[0] != 0 else 0
            print(f"  D1(p={p1}) + F1(p={p2}): shifts sum to {sums[0]:.2f} uniformly → a={a_needed:.4f}")
            found = True

            # Now check: with this 'a', do the coefficients work?
            # cd2 * coeff_D + cf2 * coeff_F = R_coeff for each block
            # Build the actual equations
            a_val = R(int(sums[0]*4), 8)  # approximate as rational
            print(f"    Approximate a = {a_val}")

if not found:
    print("  No (p1, p2) pair found with uniform H-power shift sum.")

print("\n" + "="*60)
print("SUMMARY")
print("="*60)
print("""
The dilaton coupling e^{alpha*Phi} multiplies the form stress-energy by
H^{a*alpha}, a SINGLE H-power shift applied uniformly to all blocks.

Since D1 already matches R's H-powers perfectly (shift=0), including
a dilaton coupling with alpha*a ≠ 0 would BREAK the matching.

For two forms with opposite couplings (IIB-style), the system reduces to
independent equations at each H-power, and cannot combine terms with
different H-powers to produce a single-power result (since squares are
positive).

CONCLUSION: Standard dilaton couplings cannot resolve the coefficient
mismatch found in exp04/exp05. The problem is in the COEFFICIENTS
(ratios 15/14, 5/3, 5/7) not in the H-power structure.
""")
