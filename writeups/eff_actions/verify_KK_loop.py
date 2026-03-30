"""
Verify the KK loop computation in 1__KK_loop.tex.

Computes the symmetrized rho integrals order by order, checks they are
proportional to sigma_2, sigma_3, sigma_2^2, then combines with Gamma
functions and verifies the final numerical coefficients.
"""
import sympy as sp
from itertools import permutations

s, t, u = sp.symbols('s t u')
r1, r2, r3 = sp.symbols('rho1 rho2 rho3')

# ── Schwinger integrand for a given ordering of (s, t, u) ───────────────────

def M_func(a, b, c):
    return a*r1*r2 + b*r2*r3 + c*r1*r3 + b*(r1 - r2)

def rho_integral(integrand):
    result = sp.integrate(integrand, (r1, 0, r2))
    result = sp.integrate(result, (r2, 0, r3))
    result = sp.integrate(result, (r3, 0, 1))
    return sp.expand(result)

# ── Symmetrized rho integrals ────────────────────────────────────────────────

perms = list(permutations([s, t, u]))
sigma2 = sp.expand((s*t + s*u + t*u).subs(u, -s - t))  # = -(s^2 + st + t^2)
sigma3 = sp.expand((s*t*u).subs(u, -s - t))              # = -(s^2*t + s*t^2)

print("=" * 60)
print("  Step 1: Symmetrized rho integrals")
print("=" * 60)

P_coeffs = {}  # k -> coefficient c_k such that I_k = c_k * sigma_basis

for k in range(5):
    total = sp.S(0)
    for perm in perms:
        integrand = sp.expand(((-M_func(*perm))**k) / sp.factorial(k))
        total += rho_integral(integrand)
    total = sp.cancel(total / 6)
    total_sub = sp.expand(total.subs(u, -s - t))

    if k == 0:
        assert total_sub == sp.Rational(1, 6), f"k=0 failed: {total_sub}"
        P_coeffs[0] = sp.Rational(1, 6)
        print(f"  k=0: I_0 = 1/6  ✓")
    elif k == 1:
        assert total_sub == 0, f"k=1 failed: {total_sub}"
        print(f"  k=1: I_1 = 0  (sigma_1 = 0)  ✓")
    elif k == 2:
        r = sp.cancel(total_sub / sigma2)
        assert r == sp.Rational(-1, 2160), f"k=2 ratio failed: {r}"
        P_coeffs[2] = sp.Rational(-1, 2160)
        print(f"  k=2: I_2 / sigma_2 = -1/2160  ✓")
    elif k == 3:
        r = sp.cancel(total_sub / sigma3)
        assert r == sp.Rational(1, 36288), f"k=3 ratio failed: {r}"
        P_coeffs[3] = sp.Rational(1, 36288)
        print(f"  k=3: I_3 / sigma_3 = 1/36288  ✓")
    elif k == 4:
        r = sp.cancel(total_sub / sigma2**2)
        assert r == sp.Rational(1, 1360800), f"k=4 ratio failed: {r}"
        P_coeffs[4] = sp.Rational(1, 1360800)
        print(f"  k=4: I_4 / sigma_2^2 = 1/1360800  ✓")

# ── Combine with Gamma functions ─────────────────────────────────────────────

print()
print("=" * 60)
print("  Step 2: Combine with Gamma(k - 1/2) and convert sigma -> power sums")
print("=" * 60)
print()
print("  sigma_2 = -(s^2+t^2+u^2)/2,  sigma_3 = (s^3+t^3+u^3)/3")
print("  sigma_2^2 = (s^2+t^2+u^2)^2/4")
print()

# k=2: c_2 * sigma_2 * Gamma(3/2) / M^3
#     = (-1/2160) * [-(s^2+t^2+u^2)/2] * (sqrt(pi)/2) / M^3
#     = sqrt(pi) * (s^2+t^2+u^2) / (2*2160*2) / M^3
#     = sqrt(pi) / 8640 * (s^2+t^2+u^2) / M^3
coeff_2 = sp.Rational(-1, 2160) * sp.Rational(-1, 2) * sp.gamma(sp.Rational(3, 2))
expected_2 = sp.sqrt(sp.pi) / 8640
assert sp.cancel(coeff_2 - expected_2) == 0, f"k=2 final coeff failed: {coeff_2}"
print(f"  k=2: coeff of (s^2+t^2+u^2)/M^3 = sqrt(pi)/8640  ✓")

# k=3: c_3 * sigma_3 * Gamma(5/2) / M^5
#     = (1/36288) * [(s^3+t^3+u^3)/3] * (3*sqrt(pi)/4) / M^5
#     = sqrt(pi) / (36288*4) * (s^3+t^3+u^3) / M^5
#     = sqrt(pi) / 145152 * (s^3+t^3+u^3) / M^5
coeff_3 = sp.Rational(1, 36288) * sp.Rational(1, 3) * sp.gamma(sp.Rational(5, 2))
expected_3 = sp.sqrt(sp.pi) / 145152
assert sp.cancel(coeff_3 - expected_3) == 0, f"k=3 final coeff failed: {coeff_3}"
print(f"  k=3: coeff of (s^3+t^3+u^3)/M^5 = sqrt(pi)/145152  ✓")

# k=4: c_4 * sigma_2^2 * Gamma(7/2) / M^7
#     = (1/1360800) * [(s^2+t^2+u^2)^2/4] * (15*sqrt(pi)/8) / M^7
#     = 15*sqrt(pi) / (1360800*4*8) * (s^2+t^2+u^2)^2 / M^7
#     = sqrt(pi) / 2903040 * (s^2+t^2+u^2)^2 / M^7
coeff_4 = sp.Rational(1, 1360800) * sp.Rational(1, 4) * sp.gamma(sp.Rational(7, 2))
expected_4 = sp.sqrt(sp.pi) / 2903040
assert sp.cancel(coeff_4 - expected_4) == 0, f"k=4 final coeff failed: {coeff_4}"
print(f"  k=4: coeff of (s^2+t^2+u^2)^2/M^7 = sqrt(pi)/2903040  ✓")

# ── Lattice sum ──────────────────────────────────────────────────────────────

print()
print("=" * 60)
print("  Step 3: Lattice sum")
print("=" * 60)
print()
print("  1/(M^2)^{p+1/2} = (R^2 tau_2^2)^{p+1/2} / |n+m*tau|^{2p+1}")
print("  sum = R^{2p+1} tau_2^{2p+1} sum 1/|n+m*tau|^{2p+1}")
print("       = R^{2p+1} tau_2^{p+1/2} E_{p+1/2}")
print()
print("  p=1: R^3 tau_2^{3/2} E_{3/2}   multiplies  (s^2+t^2+u^2)")
print("  p=2: R^5 tau_2^{5/2} E_{5/2}   multiplies  (s^3+t^3+u^3)")
print("  p=3: R^7 tau_2^{7/2} E_{7/2}   multiplies  (s^2+t^2+u^2)^2")

print()
print("All checks passed.")
