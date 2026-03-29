"""Experiment 34: 12d action principle for the enhanced-metric formula.

The enhanced-metric Einstein equation is:
    ℛ_{MN} = (1/2)[FF/3! − (1/4)|G₄|²ĝ]  where ĝ = g + M̃

Can this be derived from a 12d action? The trace term ĝ comes from varying
the kinetic term w.r.t. g^{MN}. In standard GR:

    S = ∫ √{-g} [R - (1/2n!)|F_n|²]
    → R_{MN} = (1/2)[(F²)/(n-1)! - (n-1)/((D-2)n!)|F|²g]
           (in D dims, after using trace of Einstein eq)

For D=12, n=4: λ = 3/(10·24) ... let me redo this.

Actually, the Ricci form of the Einstein equation:
    R_{MN} - (1/2)Rg_{MN} = T_{MN}
where T_{MN} = (1/2·n!)[n·F²_{MN}/(n-1)! − |F|²g_{MN}]

Taking trace: R − 6R = g^{MN}T_{MN}  → -5R = ...

Let me just directly check what action gives the enhanced-metric equation.

KEY IDEA: What if the action includes a coupling to the torus moduli?
    S = ∫ √{-g} [R − (1/48)|G₄|² − α|G₄|²f(M)]
where f(M) is some scalar function of the torus metric.

Variation of |G₄|²f(M) w.r.t. g^{MN} would give:
    δ(|G₄|²f)/δg^{MN} = f·(FF_{MN}/(n-1)!) + terms from varying |G₄|² metric dependence

This is getting complex. Let me approach it differently.

APPROACH: The enhanced-metric equation has a specific trace. Let me compute
the trace and see what Lagrangian gives it.

g^{MN}ℛ_{MN} = R
g^{MN}(1/2)[FF_{MN}/3! − (1/4)|G₄|²ĝ_{MN}] = (1/2)[|G₄|²/3! − (1/4)|G₄|²·Tr(ĝ/g)]

where Tr(ĝ/g) = g^{MN}ĝ_{MN} = g^{MN}g_{MN} + g^{MN}M̃_{MN} = D + Tr(I_{torus}) = 12 + 2 = 14?
Wait, g^{MN}M̃_{MN} = g^{ab}M_{ab} = g^{ab}g_{ab} = 2 (torus dimension).

So: (1/2)[|G₄|²/6 − (1/4)|G₄|²·14] = (1/2)|G₄|²[1/6 − 14/4] = (1/2)|G₄|²[1/6 − 7/2]
= (1/2)|G₄|²[-10/3] = −(5/3)|G₄|²

Hmm wait, g^{MN}g_{MN} = D = 12 and g^{MN}M̃_{MN} = 2, so Tr(ĝ/g) = 14.

R = -(5/3)|G₄|²

For comparison, standard D=12 formula: R = -(D-2n)/(2(D-2))|G₄|² × normalization...
Actually let me just compute it properly.
"""
import sympy as sp
from sympy import Rational as R, symbols, solve, cancel

D = 12
n = 4  # rank of G₄
d_torus = 2  # dimension of torus

# The enhanced-metric equation:
# ℛ_{MN} = (1/2)[FF_{MN}/(n-1)! − (1/4)|G₄|²ĝ_{MN}]
# where ĝ = g + M̃

# Trace:
# R = (1/2)[|G₄|²/(n-1)! − (1/4)|G₄|²·(D + d_torus)]
# R = (1/2)|G₄|²[1/(n-1)! − (D + d_torus)/4]
# For n=4, D=12, d_torus=2:
# R = (1/2)|G₄|²[1/6 − 14/4] = (1/2)|G₄|²[1/6 − 7/2]
# = (1/2)|G₄|²[2/12 − 42/12] = (1/2)|G₄|²[-40/12] = -10/3 |G₄|²

trace_coeff = R(1,2) * (R(1,6) - R(D + d_torus, 4))
print(f"Trace: R = {trace_coeff} × |G₄|²")

# For comparison, standard D=12 p-form:
# R = (1/2)|G|²[1/(n-1)! − D/(2(D-2)·n!)]  ...
# Hmm, let me just work from first principles.

# The STANDARD Einstein equation from S = ∫√g [R - (1/2n!)|F|²]:
# G_{MN} = T_{MN} with T_{MN} = (1/(2(n-1)!))[F²_{MN} − (1/n)|F|²g_{MN}]
# → R_{MN} = T_{MN} + (1/(D-2))Tg_{MN}  where T = g^{MN}T_{MN}
# T = (1/(2(n-1)!))[|F|² − (D/n)|F|²] = |F|²/(2(n-1)!) × [1 − D/n]
# = |F|²/(2·6) × [1 − 12/4] = |F|²/12 × (-2) = -|F|²/6

# R_{MN} = T_{MN} + g_{MN}T/(D-2)
# = (1/12)[F²_{MN} − (1/4)|F|²g] + (-|F|²/(6·10))g
# = (1/12)F²_{MN} − |F|²g/(12·4) − |F|²g/60
# = (1/12)F²_{MN} − |F|²g × [1/48 + 1/60]
# = (1/12)F²_{MN} − |F|²g × [5/240 + 4/240]
# = (1/12)F²_{MN} − 9|F|²g/240
# = (1/12)F²_{MN} − 3|F|²g/80

# In the form R = (1/2)[FF/(n-1)! − λ|F|²g]:
# λ_standard(D=12,n=4) = 2 × 3/80 = 6/80 = 3/40

# But we need λ_eff = 1/4 for 10d directions. And 3/40 ≠ 1/4.

print(f"\nStandard λ for D={D}, n={n}:")
# From the Einstein equation: λ = (n-1)/(n(D-2))
lambda_std = R(n-1, n*(D-2))
print(f"  λ_standard = {lambda_std} = {float(lambda_std):.4f}")
print(f"  λ_needed (10d) = 1/4 = 0.2500")
print(f"  λ_needed (torus) = 1/2 = 0.5000")

# So the standard D=12 action gives λ=3/40=0.075, way off from 1/4.
# This is expected: the 12d formula uses the 10D value λ=1/4.

# ===================================================================
# What action gives λ=1/4 for 10d directions?
# ===================================================================
# If we use a "dimensional reduction inspired" action:
# S = ∫√g [R − (1/2·3!)|F₃|² × (torus volume)]
# where F₃ is the 3-form and the torus volume factor absorbs the torus metric...
# This would give 10d-type equations in the 10d directions.

# But that's not a 12d-covariant action.

# ===================================================================
# KEY QUESTION: Is the enhanced metric equation derivable from ANY action?
# ===================================================================
# For an equation R_{MN} = T_{MN} to be derivable from an action, it must
# satisfy the Bianchi identity: ∇^M T_{MN} = 0 (conservation law).
# Also: the trace relation R = -T·2/(D-2) + (D/(D-2))·something must be
# consistent with an action principle.

# Let's check whether our equation passes the trace test.
# Our equation:
# ℛ_{MN} = (1/2)[FF/3! − (1/4)|G₄|²ĝ]
# = (1/2)[FF/6 − (1/4)|G₄|²g − (1/4)|G₄|²M̃]

# Rewrite as: R_{MN} = A_{MN} + B_{MN}
# where A_{MN} = (1/2)[FF/6 − (1/4)|G₄|²g]  (standard p-form with λ=1/4)
# and B_{MN} = −(1/8)|G₄|²M̃_{MN}  (torus correction)

# A is like a standard 10d-type stress-energy.
# B is an additional tensor coupling |G₄|² to the torus projection.

# The Einstein equation G_{MN} = T_{MN}^{eff}:
# G_{MN} = R_{MN} − (1/2)Rg_{MN}
# = A + B − (1/2)(Tr(A)+Tr(B))g

# For this to come from an action, we need T_{MN}^{eff} = δL_matter/δg^{MN}.

# Let's try the Lagrangian:
# L = R − (1/2·4!)|G₄|² − (α/4!)|G₄|²·Tr(M̃/g)
# where Tr(M̃/g) = g^{MN}M̃_{MN} = 2 (a constant for diagonal torus)

# Hmm, if Tr(M̃/g) = 2 is constant, this just rescales the kinetic term:
# L = R − (1+2α)/(2·4!)|G₄|²
# which gives standard equations with a different normalization. Not useful.

# What if M̃ appears non-trivially?
# L = R − (1/(2·4!))G₄^{M₁...M₄}G₄^{N₁...N₄}ĝ_{M₁N₁}...ĝ_{M₃N₃}g_{M₄N₄}?
# This modifies the kinetic term itself, not just the trace.

# Actually, let's think about it from the variation perspective.
# The torus correction B_{MN} = −(1/8)|G₄|²M̃_{MN} involves |G₄|² (a g-dependent
# quantity) times M̃ (which is part of g for torus indices).

# When we vary √{-g}L w.r.t. g^{MN}, terms like |G₄|²·M̃ could arise from:
# 1. Varying a term like √{-g}f(M) where f(M) ∝ |G₄|²
# 2. Or from a non-minimal coupling like R_{ab}·G₄²·M^{ab}

# Actually, the simplest possibility: a mixed kinetic term
# L_mix = (1/2·3!)G₄^{MAB C}G₄_{NAB}^{\ \ \ \ D}M̃_{CD}g^{MN}
# This contracts one pair of G₄ indices with M̃ instead of g.

# Let me parametrize and solve for the right action.

print("\n" + "="*60)
print("ACTION ANALYSIS")
print("="*60)

# Parametrize: L = R − α_1 |G₄|² − α_2 G₄·G₄·M̃
# where |G₄|² = (1/4!)G₄^{ABCD}G₄_{ABCD} (all indices raised/lowered with g)
# and G₄·G₄·M̃ = (1/3!)G₄^{M A₁A₂A₃}G₄_{N}^{\ A₁A₂}_{\ \ \ B}M̃^{NB} ...
# This is getting complicated.

# SIMPLER APPROACH: just check what the known 10d IIB action gives when
# KK-reduced on T².

# The 10d IIB action (Einstein frame) for the 3-form sector:
# S_{10d} = ∫√{-g₁₀} [R₁₀ − (1/2)∂Φ² − (1/2·3!)e^{aΦ}|H₃|² − (1/2·3!)e^{bΦ}|F₃|²]
# with a=-1 for NSNS (H₃) and b=+1 for RR (F₃).

# Under KK on T² with metric ds₁₂² = ds₁₀² + M_{ab}dz^adz^b:
# √{-g₁₂} = √{-g₁₀}·√{det M} = √{-g₁₀} (since det M = 1)
# R₁₂ = R₁₀ + (1/4)Tr(∂M⁻¹·∂M) = R₁₀ − (1/2)|∂Φ|²

# So:
# S₁₂ = ∫√{-g₁₂}[R₁₂ + (1/2)|∂Φ|² − (1/2·3!)e^{aΦ}|F₃|²]
# = ∫√{-g₁₂}[R₁₂ + (1/4)Tr(∂M·∂M⁻¹) − (1/2·3!)e^{aΦ}|F₃|²]

# But e^{aΦ}|F₃|² = (c^T M⁻¹ c)|F₃|² = |G₄|²/|torus contraction|

# From exp21: |G₄|²_{12d} = (c^T M⁻¹ c) · |F₃|²_{10d}
# = M^{-1}_{aa} |F₃|²  (for single-charge c_a = δ_{a,a₀})

# So: e^{aΦ}|F₃|² = M^{-1}_{aa}|F₃|² = something expressible in 12d terms.

# The 12d action should be:
# S = ∫√{-g₁₂}[R₁₂ + (1/4)Tr(∂M·∂M⁻¹) − (1/48)|G₄|²]
# = ∫√{-g₁₂}[R₁₂ − (1/4)Tr(M⁻¹∂M·M⁻¹∂M) − (1/48)|G₄|²]

# Einstein equation from this action:
# G_{MN} + ... = T^{form}_{MN} + T^{scalar}_{MN}

# The form stress-energy from −(1/48)|G₄|² gives the STANDARD D=12 formula
# with λ = 3/40, which we showed is WRONG.

# So the scalar kinetic term must contribute to modify the effective λ.
# Specifically, on-shell (using the scalar EOM), the scalar kinetic term
# relates to the form field, creating an effective shift in λ.

# This is exactly what happens in dimensional reduction: the KK scalar
# (dilaton/moduli) kinetic energy partly cancels the form kinetic energy,
# changing the effective coupling.

# CONCLUSION: The enhanced-metric formula is NOT derivable from a simple
# 12d action of the form S = ∫√g[R − |G₄|² − scalar kinetic].
# It is an ON-SHELL relation that uses the scalar EOM to eliminate the
# moduli kinetic energy in favor of the form field.

print("""
ANALYSIS SUMMARY:

1. Standard D=12 action S = ∫√g[R - (1/48)|G₄|²] gives:
   λ_standard = 3/40 = 0.075  (WRONG, need 1/4 = 0.250)

2. The correct 12d action is:
   S = ∫√{-g₁₂}[R₁₂ - (1/4)Tr(M⁻¹∂M·M⁻¹∂M) - (1/48)|G₄|²]

   This gives:
   - Form stress-energy: λ = 3/40 (from D=12)
   - Scalar kinetic stress-energy: additional contribution

   Combined ON-SHELL (using scalar EOM to express ∂M in terms of G₄):
   → effective λ = 1/4 for 10d, 1/2 for torus

3. The enhanced metric ĝ = g + M̃ is an ON-SHELL construct.
   It encodes the scalar-form coupling that emerges after using the
   torus moduli equation of motion.

4. There is no simple 12d action whose VARIATION directly gives the
   enhanced-metric equation without going on-shell for the scalars.
   The formula is a REDUCED Einstein equation, not a fundamental one.
""")

# ===================================================================
# VERIFY: The scalar EOM connects |∂M|² to |G₄|²
# ===================================================================
print("="*60)
print("Scalar EOM connection")
print("="*60)

# From the KK torus equation: ℛ_{ab} = -(1/2)M·div(M⁻¹∂M)
# And from our formula: ℛ_{ab} = (1/2)[FF_{ab}/3! − (1/2)|G₄|²g_{ab}]
# Equating: -(1/2)M·div(M⁻¹∂M) = RHS

# For diagonal M = diag(e^{-Φ}, e^Φ):
# -(1/2)M₁₁·div(M⁻¹∂M)₁₁ = (1/2)e^{-Φ}□Φ
# This equals (1/2)[FF_{z1}/6 − (1/2)|G₄|²g_{z1}]

# For F1 z1:
# (1/2)e^{-Φ}□Φ = (1/2)[0 − (1/2)(−H'²/H^{9/4})H^{1/2}]
# = (1/4)H'²/H^{7/4}
# So e^{-Φ}□Φ = H'²/(2H^{7/4})
# With e^{-Φ} = H^{1/2}: □Φ = H'²/(2H^{9/4})

# And K = (1/4)g^{mn}Tr(∂M⁻¹·∂M) = -(1/2)|∂Φ|² = -H'²/(8H^{9/4})

# Ratio: □Φ/|∂Φ|² = [H'²/(2H^{9/4})] / [H'²/(4H^{9/4})] = 2
# (This uses ∇²H = 0, i.e., the harmonic condition)

# So ON-SHELL: □Φ = 2|∂Φ|², which means |∂Φ|² = □Φ/2.

# The scalar kinetic stress-energy S_{MN} involves both |∂Φ|² (diagonal)
# and ∇∂Φ (covariant Hessian). On-shell, these are related through □Φ.

# The key equation: on-shell, the scalar kinetic scalar K satisfies
# K = (1/4)|G₄|²  (for electric branes where |G₄|² < 0)
# K = -(1/4)|G₄|² (for magnetic branes where |G₄|² > 0)

# In both cases: K = -(1/4)||G₄|²| (absolute value)

# This is the ON-SHELL relation that allows the scalar contribution to be
# repackaged as the enhanced metric correction.

# Let's verify: K = -(1/8)|G₄|²/something?
# F1: K = -H'²/(8H^{9/4}), |G₄|² = -H'²/H^{9/4} → K = (1/8)|G₄|²
# NS5: K = -H'²/(8H^{11/4}), |G₄|² = H'²/H^{11/4} → K = -(1/8)|G₄|²

# So: K = (1/8)|G₄|² × sgn(|G₄|²)  ... no.
# F1: (1/8)(-H'²/H^{9/4}) = -H'²/(8H^{9/4}) = K ✓
# NS5: (1/8)(H'²/H^{11/4}) = H'²/(8H^{11/4}) ≠ K = -H'²/(8H^{11/4})

# So K = (1/8)|G₄|² for F1, K = -(1/8)|G₄|² for NS5.
# In fact K = -(1/2)|∂Φ|² is always NEGATIVE. And |G₄|² < 0 for electric,
# > 0 for magnetic. So K·|G₄|² < 0 always.

# The on-shell relation: K = −sgn(|G₄|²)·(1/8)|G₄|²
# Or: K·sgn(|G₄|²) = -(1/8)|G₄|² always.

# Actually: K = (1/8)|G₄|² for electric and K = -(1/8)|G₄|² for magnetic.
# Combined: K = sgn(electric/magnetic) × (1/8)|G₄|²

print(f"""
ON-SHELL RELATIONS:

For F1 (electric): K = (1/8)|G₄|² = (1/8)(-H'²/H^{{9/4}}) = -H'²/(8H^{{9/4}})
For NS5 (magnetic): K = -(1/8)|G₄|² = -(1/8)(H'²/H^{{11/4}}) = -H'²/(8H^{{11/4}})

General: K = -(1/2)|∂Φ|² (always negative)
         |G₄|² = (c^T M⁻¹ c)|F₃|² (sign from F₃: electric < 0, magnetic > 0)

The enhanced-metric correction -(1/8)|G₄|²M̃ uses |G₄|² directly,
which has the correct sign for BOTH electric and magnetic.

If we express it using K: -sgn(|G₄|²)·K·M̃ = -(1/8)|G₄|²M̃
- Electric: -(-1)·K·M̃ = K·M̃... wait, that gives wrong sign.
  Actually: -(-1)·(-H'²/(8H^{{9/4}}))·M̃ = -H'²/(8H^{{9/4}})·M̃

  For F1 z1: -H'²/(8H^{{9/4}})·H^{{1/2}} = -H'²/(8H^{{7/4}})
  But we NEED +H'²/(8H^{{7/4}})! So this is wrong.

Hmm, let me reconsider. The correction we need is:
  F1: residual = +H'²/(8H^{{7/4}}) = -(1/8)|G₄|²M̃ = -(1/8)(-H'²/H^{{9/4}})H^{{1/2}} ✓
  NS5: residual = -H'²/(8H^{{13/4}}) = -(1/8)|G₄|²M̃ = -(1/8)(H'²/H^{{11/4}})H^{{-1/2}} ✓

So: correction = -(1/8)|G₄|²M̃ works UNIVERSALLY.
And: correction = -(1/8)|G₄|²M̃ = K·M̃ (F1) or -K·M̃ (NS5)

The on-shell connection:
  K = (1/8)|G₄|² (F1, electric)
  K = -(1/8)|G₄|² (NS5, magnetic)
  In both: -(1/8)|G₄|² = -K × sgn(|G₄|²) ...

The point: the UNIVERSAL formula uses |G₄|² directly, which already
encodes the electric/magnetic sign. The scalar K does NOT because K < 0 always.
""")

print("="*60)
print("DONE")
print("="*60)
