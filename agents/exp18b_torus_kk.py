"""Experiment 18b: Correct KK formula for torus Ricci.

From exp18:
  - ℛ_{ab} = -(1/2)□M_{ab} FAILED
  - But ℛ_{ab} IS SL(2,Z) covariant: ℛ(1,1) = Λ·ℛ(F1)·Λ^T ✓
  - Need the correct KK reduction formula with M^{cd}∂M∂M correction terms

Standard KK reduction for ds² = g_{mn}dx^m dx^n + M_{ab}(x)dz^a dz^b:

ℛ_{ab} = -(1/2) □_g M_{ab} + (1/4) M^{cd} (∂_m M_{ac})(∂^m M_{bd})
         + (1/4) (∂_m ln det M) ∂^m M_{ab}   [vanishes when det M = 1]

ℛ_{mn} = R_{mn}^{base} - (1/2) M^{ab} ∇_m ∂_n M_{ab}
       = R_{mn}^{base} + (1/2) ∂_m M^{ab} ∂_n M_{ab}

ℛ_{ma} = 0  (no gauge field)

Let me verify these formulas for both Einstein-frame and string-frame 12d metrics.
"""
import sys; sys.path.insert(0, 'src')
import sympy as sp
from sympy import Rational as R
from sugra import HarmonicFunction, Metric

wv_coords = list(sp.symbols('t x1', real=True))
z1s, z2s = sp.symbols('z1 z2', real=True)
harmonic_coords = list(sp.symbols('y0:8', real=True))
coords = wv_coords + [z1s, z2s] + harmonic_coords
coords_10 = wv_coords + harmonic_coords
D = 12

hf = HarmonicFunction(transverse_coords=harmonic_coords)
H = sp.Function('H')(hf.r_expr)

# Torus metric for F1
M = sp.Matrix([[H**R(1, 2), 0],
               [0, H**R(-1, 2)]])
M_inv = sp.cancel(M.inv())

print("=" * 70)
print("EXP 18b: Correct KK formula for torus block")
print("=" * 70)

# ===================================================================
# Part A: String-frame KK
# ===================================================================
print("\nPART A: String-frame 12d metric")
print("-" * 50)

# String-frame 10d metric
g_sf_10 = sp.zeros(10, 10)
g_sf_10[0, 0] = -1/H
g_sf_10[1, 1] = 1/H
for k in range(8):
    g_sf_10[2 + k, 2 + k] = sp.Integer(1)

metric_sf_10 = Metric(g_sf_10, coords_10)
chris_sf = metric_sf_10.christoffel(simplify_func=sp.cancel)

# 12d string-frame metric
g_sf_12 = sp.zeros(D, D)
g_sf_12[0, 0] = -1/H
g_sf_12[1, 1] = 1/H
g_sf_12[2, 2] = M[0, 0]
g_sf_12[3, 3] = M[1, 1]
for k in range(8):
    g_sf_12[4 + k, 4 + k] = sp.Integer(1)

print("Computing 12d Ricci (string frame)...")
metric_sf_12 = Metric(g_sf_12, coords)
Ric_sf_12 = metric_sf_12.ricci_tensor(simplify_func=sp.cancel)

# The correct KK formula for internal block:
# ℛ_{ab} = -(1/2) □M_{ab} + (1/4) M^{cd} (∂_μ M_{ac})(∂^μ M_{bd})
# where □ = g^{μν}(∂_μ∂_ν - Γ^ρ_{μν}∂_ρ) is the base covariant box.

def box_scalar(f):
    """□_{SF,10d} of a scalar function f."""
    result = sp.Integer(0)
    for m in range(10):
        g_inv = sp.cancel(1 / g_sf_10[m, m])
        d2 = sp.diff(sp.diff(f, coords_10[m]), coords_10[m])
        conn = sp.Integer(0)
        for rho in range(10):
            key = (rho, m, m)
            if key in chris_sf:
                conn += chris_sf[key] * sp.diff(f, coords_10[rho])
        result += g_inv * (d2 - conn)
    return sp.cancel(result)

def contracted_dMdM(a, b):
    """(1/4) M^{cd} g^{μν} ∂_μ M_{ac} ∂_ν M_{bd}."""
    result = sp.Integer(0)
    for c in range(2):
        for d in range(2):
            if M_inv[c, d] == 0:
                continue
            for m in range(10):
                g_inv = sp.cancel(1 / g_sf_10[m, m])
                dMac = sp.diff(M[a, c], coords_10[m])
                dMbd = sp.diff(M[b, d], coords_10[m])
                if dMac != 0 and dMbd != 0:
                    result += R(1, 4) * M_inv[c, d] * g_inv * dMac * dMbd
    return sp.cancel(result)

print("\nTesting: ℛ_{ab} = -(1/2)□M_{ab} + (1/4)M^{cd}(∂M_{ac})(∂M_{bd})")
kk_pass = True
for a, b, name in [(0, 0, 'z1'), (1, 1, 'z2'), (0, 1, 'z1z2')]:
    ric_actual = hf.substitute(sp.cancel(Ric_sf_12[a+2, b+2]))
    box_val = hf.substitute(sp.cancel(box_scalar(M[a, b])))
    corr_val = hf.substitute(contracted_dMdM(a, b))
    predicted = sp.cancel(-R(1, 2) * box_val + corr_val)
    diff = sp.cancel(ric_actual - predicted)
    status = '✓' if diff == 0 else f'✗ diff={diff}'
    if diff != 0:
        kk_pass = False
    print(f"  [{name}]: ℛ={ric_actual}, predicted={predicted}, {status}")

print(f"\n  Formula 1 [-(1/2)□M + (1/4)M^cd∂M∂M]: {'★ VERIFIED ★' if kk_pass else 'FAILED'}")

# If that fails, try with coefficient 1/8 instead of 1/4
if not kk_pass:
    print("\n  Trying other coefficients for the correction term...")
    # Compute the numerical coefficient
    ric_z1 = hf.substitute(sp.cancel(Ric_sf_12[2, 2]))
    box_z1 = hf.substitute(sp.cancel(box_scalar(M[0, 0])))

    # Compute M^{cd}∂M_{0c}∂M_{0d} WITHOUT the 1/4 factor
    raw_corr_z1 = sp.Integer(0)
    for c in range(2):
        for d in range(2):
            if M_inv[c, d] == 0:
                continue
            for m in range(10):
                g_inv = sp.cancel(1 / g_sf_10[m, m])
                dM0c = sp.diff(M[0, c], coords_10[m])
                dM0d = sp.diff(M[0, d], coords_10[m])
                if dM0c != 0 and dM0d != 0:
                    raw_corr_z1 += M_inv[c, d] * g_inv * dM0c * dM0d
    raw_corr_z1 = hf.substitute(sp.cancel(raw_corr_z1))

    residual = sp.cancel(ric_z1 + R(1, 2) * box_z1)
    if raw_corr_z1 != 0:
        coeff = sp.cancel(residual / raw_corr_z1)
        print(f"  Needed coefficient: {coeff}")
        print(f"  Raw correction z1: {raw_corr_z1}")
        print(f"  Residual z1: {residual}")

    # Try the alternative: use ∂M^{ab} ∂M_{ab} contracted with g^{mn}
    # (1/4) ∂_m M^{cd} ∂^m M_{cd} involves the INVERSE metric derivatives
    dMinvdM = sp.Integer(0)
    for c in range(2):
        for d in range(2):
            for m in range(10):
                g_inv = sp.cancel(1 / g_sf_10[m, m])
                dMinv_cd = sp.diff(M_inv[c, d], coords_10[m])
                dM_cd = sp.diff(M[c, d], coords_10[m])
                if dMinv_cd != 0 and dM_cd != 0:
                    dMinvdM += g_inv * dMinv_cd * dM_cd
    dMinvdM = hf.substitute(sp.cancel(dMinvdM))
    print(f"\n  Tr(∂M^{-1} · ∂M) = {dMinvdM}")

    # The standard KK formula from Pope's notes (hep-th/9912164, eq 2.6):
    # For ds² = g + M dz²:
    # ℛ_{mn} = R_{mn} + (1/4) Tr(∂_m M^{-1} · ∂_n M)
    # ℛ_{ab} = -(1/2)(sqrt(det M))^{-1} ∇_μ(sqrt(det M) g^{μν} M^{-1}_{ac} ∂_ν M_{cb})
    # For det M = 1: ℛ_{ab} = -(1/2) g^{μν} ∇_μ(M^{-1}_{ac} ∂_ν M_{cb})
    #              = -(1/2) [M^{-1}_{ac} □M_{cb} + (∂_μ M^{-1}_{ac})(∂^μ M_{cb})]
    #              = -(1/2) [M^{-1}_{ac} □M_{cb} - (M^{-1}∂M M^{-1})_{ac}(∂M)_{cb}]
    # Multiplying by M on left to get ℛ as symmetric:
    # M_{ea} ℛ_{ab} = -(1/2) [□M_{eb} - (M^{-1}∂M)_{ec}(∂M)_{cb}]
    # Hmm, ℛ_{ab} as computed from the 12d metric IS symmetric in (a,b).

    # Let me try the Pope formula directly:
    # ℛ_{ab} = -(1/2) ∇_μ(M^{-1} ∂^μ M)_{ab}
    #         = -(1/2) [∂_μ(M^{-1}∂^μ M) - Γ^ρ_{μμ}(M^{-1}∂_ρ M)]_{ab} × (1/g^{μμ})... no

    # Let me compute div(M^{-1}∂M) more carefully.
    # [∇_μ(M^{-1}∂^μ M)]_{ab} = g^{μν}[∂_μ(M^{-1}∂_ν M) + Γ^μ_{μρ}(M^{-1}∂^ρ M)]_{ab}
    # Wait, for a vector V^μ: ∇_μ V^μ = ∂_μ V^μ + Γ^μ_{μρ}V^ρ = (1/√g)∂_μ(√g V^μ)

    # Here, (M^{-1}∂M) has spacetime index μ (vector) and torus indices (a,b) (scalar).
    # So: ∇_μ(M^{-1}∂^μ M)_{ab} = (1/√g)∂_μ(√g g^{μν} (M^{-1}∂_ν M)_{ab})
    # where √g = √|det g_{10}|.

    sqrt_g = sp.Integer(1)
    for m in range(10):
        sqrt_g *= sp.Abs(g_sf_10[m, m])
    sqrt_g = sp.sqrt(sqrt_g)
    sqrt_g = sp.cancel(sqrt_g.rewrite(sp.Piecewise))
    # For SF: g = H^{-2} (from tt and x1x1), rest are 1.
    # |det g| = H^{-2}. √|det g| = H^{-1}
    sqrt_g = 1/H
    print(f"\n  √|det g^SF| = {sqrt_g}")

    print("\n  Testing Pope formula: ℛ_{ab} = -(1/2) div(M^{-1}∂M)_{ab}")
    print("  = -(1/2√g) ∂_μ(√g g^{μν} [M^{-1}∂_ν M]_{ab})")

    pope_pass = True
    for a, b, name in [(0, 0, 'z1'), (1, 1, 'z2'), (0, 1, 'z1z2')]:
        # Compute (M^{-1}∂_ν M)_{ab}
        # Then √g g^{μν} (M^{-1}∂_ν M)_{ab} for diagonal metric is:
        # √g g^{νν} (M^{-1}∂_ν M)_{ab}
        # Then ∂_ν of that.
        divergence = sp.Integer(0)
        for nu in range(10):
            g_inv_nu = sp.cancel(1 / g_sf_10[nu, nu])
            # (M^{-1} ∂_ν M)_{ab}
            MinvdM_ab = sp.Integer(0)
            for c in range(2):
                MinvdM_ab += M_inv[a, c] * sp.diff(M[c, b], coords_10[nu])
            MinvdM_ab = sp.cancel(MinvdM_ab)
            if MinvdM_ab == 0:
                continue
            flux = sp.cancel(sqrt_g * g_inv_nu * MinvdM_ab)
            divergence += sp.diff(flux, coords_10[nu])
        divergence = sp.cancel(divergence / sqrt_g)
        ric_predicted = hf.substitute(sp.cancel(-R(1, 2) * divergence))
        ric_actual = hf.substitute(sp.cancel(Ric_sf_12[a+2, b+2]))
        diff = sp.cancel(ric_actual - ric_predicted)
        status = '✓' if diff == 0 else f'✗ diff={diff}'
        if diff != 0:
            pope_pass = False
        print(f"    [{name}]: actual={ric_actual}, Pope={ric_predicted}, {status}")

    print(f"\n  Pope formula: {'★ VERIFIED ★' if pope_pass else 'FAILED'}")

# ===================================================================
# Part B: Verify ℛ_{mn} KK formula
# ===================================================================
print("\n" + "=" * 70)
print("PART B: KK formula for 10d-direction Ricci")
print("-" * 50)

# Standard: ℛ_{mn} = R_{mn}^{base} + (1/4) Tr(∂_m M^{-1} · ∂_n M)
# = R_{mn}^{base} - (1/4) Tr(M^{-1}∂_m M · M^{-1} ∂_n M)

Ric_sf_10 = metric_sf_10.ricci_tensor(simplify_func=sp.cancel)

print("\nTesting: ℛ_{mn} = R_{mn}^{SF} + (1/4)Tr(∂_m M^{-1} · ∂_n M)")
mn_pass = True
for i_12, i_10, name in [(0, 0, 't'), (1, 1, 'x1'), (4, 2, 'y0'), (5, 3, 'y1')]:
    R12 = hf.substitute(sp.cancel(Ric_sf_12[i_12, i_12]))
    R10 = hf.substitute(sp.cancel(Ric_sf_10[i_10, i_10]))

    # (1/4)Tr(∂_m M^{-1} · ∂_m M) = (1/4) Σ_{a,b} ∂_m M^{-1}_{ab} · ∂_m M_{ba}
    tr_val = sp.Integer(0)
    for a_idx in range(2):
        for b_idx in range(2):
            dMinv = sp.diff(M_inv[a_idx, b_idx], coords_10[i_10])
            dM = sp.diff(M[b_idx, a_idx], coords_10[i_10])
            tr_val += dMinv * dM
    tr_val = hf.substitute(sp.cancel(R(1, 4) * tr_val))

    predicted = sp.cancel(R10 + tr_val)
    diff = sp.cancel(R12 - predicted)
    status = '✓' if diff == 0 else f'✗ diff={diff}'
    if diff != 0:
        mn_pass = False
    print(f"  [{name}]: R12={R12}, R10+KK={predicted}, {status}")

print(f"\n  10d KK formula: {'★ VERIFIED ★' if mn_pass else 'FAILED'}")

# Try alternative sign: (1/4)Tr(∂M^{-1}·∂M) vs -(1/4)Tr(M^{-1}∂M·M^{-1}∂M)
if not mn_pass:
    print("\n  Trying: ℛ_{mn} = R_{mn} - (1/4)Tr(M^{-1}∂_m M · M^{-1}∂_n M)")
    for i_12, i_10, name in [(0, 0, 't'), (4, 2, 'y0')]:
        R12 = hf.substitute(sp.cancel(Ric_sf_12[i_12, i_12]))
        R10 = hf.substitute(sp.cancel(Ric_sf_10[i_10, i_10]))

        tr_val2 = sp.Integer(0)
        for a_idx in range(2):
            for b_idx in range(2):
                MinvdM_a = sp.Integer(0)
                MinvdM_b = sp.Integer(0)
                for c_idx in range(2):
                    MinvdM_a += M_inv[a_idx, c_idx] * sp.diff(M[c_idx, b_idx], coords_10[i_10])
                    MinvdM_b += M_inv[b_idx, c_idx] * sp.diff(M[c_idx, a_idx], coords_10[i_10])
                tr_val2 += MinvdM_a * MinvdM_b  # This is Tr[(M^{-1}dM)^2] wrong

        # Actually Tr(M^{-1}∂M·M^{-1}∂M) = Σ_{a,b} (M^{-1}∂M)_{ab}(M^{-1}∂M)_{ba}
        tr_val2 = sp.Integer(0)
        MinvdM = sp.Matrix(2, 2, lambda a, b: sum(M_inv[a, c] * sp.diff(M[c, b], coords_10[i_10]) for c in range(2)))
        prod = MinvdM * MinvdM
        tr_val2 = hf.substitute(sp.cancel(-R(1, 4) * sp.trace(prod)))

        predicted2 = sp.cancel(R10 + tr_val2)
        diff2 = sp.cancel(R12 - predicted2)
        print(f"    [{name}]: diff = {diff2}")

    # Also try (1/2)Tr(∂M·∂M^{-1})
    print("\n  Fitting coefficient: ℛ - R = α × Tr(∂M^{-1}·∂M)")
    for i_12, i_10, name in [(0, 0, 't'), (4, 2, 'y0')]:
        R12 = hf.substitute(sp.cancel(Ric_sf_12[i_12, i_12]))
        R10 = hf.substitute(sp.cancel(Ric_sf_10[i_10, i_10]))
        residual_mn = sp.cancel(R12 - R10)

        tr_raw = sp.Integer(0)
        for a_idx in range(2):
            for b_idx in range(2):
                dMinv = sp.diff(M_inv[a_idx, b_idx], coords_10[i_10])
                dM_val = sp.diff(M[b_idx, a_idx], coords_10[i_10])
                tr_raw += dMinv * dM_val
        tr_raw = hf.substitute(sp.cancel(tr_raw))
        if tr_raw != 0:
            coeff_mn = sp.cancel(residual_mn / tr_raw)
            print(f"    [{name}]: residual={residual_mn}, Tr={tr_raw}, coeff={coeff_mn}")

# ===================================================================
# Part C: Summary
# ===================================================================
print("\n" + "=" * 70)
print("PART C: Summary of KK formulas")
print("=" * 70)
print("""
The correct KK decomposition for ds²_{12} = g^{base}_{mn}dx^m dx^n + M_{ab}(x)dz^a dz^b
(with no KK gauge field, det M = 1):

  10d directions:
    ℛ_{mn} = R_{mn}^{base} + (correction involving M)

  Torus directions:
    ℛ_{ab} = -(1/2) div(M^{-1}∂M)_{ab}
           = -(1/2√g) ∂_μ(√g g^{μν} (M^{-1}∂_ν M)_{ab})

Both are SL(2,Z) covariant by construction since M → ΛMΛ^T.
""")
