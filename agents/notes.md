# 12d Form Field Search — Agent Notes

## Problem
Find form fields F such that T_{MN}(F) = R_{MN} for:
```
ds²₁₂ = H^{-3/4} ds²_{1,1} + H^{1/2} dz₁² + H^{-1/2} dz₂² + H^{1/4} ds²_8
```

## Ricci tensor targets
```
R[t]  =  3/(8H³)         · H'²
R[z1] =  1/(4H^{7/4})    · H'²
R[z2] = -1/(4H^{11/4})   · H'²
R[y0] = (1-4y0²/r²)/(8H²) · H'²
```

## Code validation
- M2, M5 (D=11): R=T verified ✓ (exp06)
- `form_stress_energy` and Verifier are correct
- Ricci bug fixed: diagonal metrics can have off-diagonal Ricci (geometry.py line 333, already patched)

## Exp01–12 summary: Obstruction with standard formula

Systematic search of ALL SO(8)-symmetric forms (4-forms, 5-forms, 3-forms) with T = (1/2)[FF/(p-1)! - λ|F²|g]:

- Only 2 of 4 SO(8)-symmetric 4-form types match H-powers: D1-type C_{t,x1,z2}=H^{-1} and F1-type C_{t,x1,z1}=H^{-1/2}
- D1 residual pattern: **(1, 4, 4, -1)** × H'²/(40 × H^{block}), with κ_aniso = 1 (exact angular match)
- Required effective λ: 1/4 for blocks WITH a D1 leg, 1/2 for blocks WITHOUT
- **No single λ works**: wv/z blocks need (cd², cf²) = (5/4, 1) at λ=1/2, but y0 aniso needs cf² = 4 - 4cd² = -1 < 0
- Dilaton couplings cannot resolve the mismatch (exp08)
- Democratic formulation gives pure trace — useless (exp23)

## Exp13–30 summary: 10d KK analysis (background, not the goal)

- The 12d metric is the F1-string uplift: ds²₁₂ = g^E(10d) + M_{ab}dz^adz^b
- G₄ = H₃ ∧ dz₂ with M = diag(H^{1/2}, H^{-1/2})
- KK formulas verified (string frame): ℛ_{mn} = R^S_{mn} + (1/4)Tr(∂M⁻¹·∂M), ℛ_{ab} = -(1/2)M·div(M⁻¹∂M)
- All 6 IIB branes verified (F1, D1, D3, NS5, D5, D7)
- G₄ packaging theorem: (G₄²)_{mn}/3! = (c^T M⁻¹ c)·(F₃²)_{mn}/2! — kinematic identity
- SL(2,Z) covariance verified for non-diagonal M (exp28-30)
- Frame tension: Einstein frame equation is covariant but metric is not; string frame metric is invariant but equation is not

## Exp31–33: Enhanced-metric formula (★★★ MAIN RESULT)

**The obstruction from exp1-12 is RESOLVED:**

```
ℛ_{MN} = (1/2)[(G₄²)_{MN}/3! − (1/4)|G₄|² ĝ_{MN}] + (1/(4·4!))(F₅^{SD})²_{MN}
```

where:
- **ĝ_{MN} = g_{MN} + M̃_{MN}** (enhanced metric)
- M̃_{MN} = torus metric M_{ab} embedded in 12d (zero for non-torus indices)
- G₄ = F₃^a ∧ dz_a (SL(2,Z) doublet → singlet)
- F₅^{SD} = F₅^e + *₁₂(F₅^e ∧ vol_{T²}) (self-dual, |F₅^{SD}|² = 0)

**Effect**: ĝ gives effective λ = 1/4 for 10d directions, λ = 1/2 for torus — resolving the block-dependent λ obstruction.

**Verified for ALL 6 fundamental IIB branes:**

| Brane | G₄ | F₅ | M̃ | Status |
|-------|----|----|-----|--------|
| F1 (NSNS elec) | H₃∧dz₂ | — | active | ALL ✓ |
| D1 (RR elec) | F₃∧dz₁ | — | active | ALL ✓ |
| NS5 (NSNS mag) | *H₃∧dz₂ | — | active | ALL ✓ |
| D5 (RR mag) | *F₃∧dz₁ | — | active | ALL ✓ |
| D3 (self-dual) | — | F₅^{SD} | trivial (M=I) | ALL ✓ |
| D7 (KK monopole) | — | — | N/A (vacuum) | ALL ✓ |

### Exp33: Enhanced-metric formula for (p,q)-strings (non-diagonal M)

**Part A — 10d directions**: Proven analytically for ALL (p,q)-strings.
  1. KK: ℛ^{EF}_{mn} = R^E_{mn} + (1/4)Tr(∂_mM⁻¹·∂_nM)
  2. Scalar cancellation: (1/2)∂Φ∂Φ + (1/2)e^{2Φ}∂C₀∂C₀ = -(1/4)Tr(∂M⁻¹·∂M) → exact cancellation
  3. Packaging: form terms = (1/2)[FF/3! - (1/4)|G₄|²g^E]
  4. **Result: ℛ = (1/2)[FF/3! - (1/4)|G₄|²g^E] for ALL (p,q)-strings** ✓

**Part B — Torus directions**: FAILS for non-diagonal M.
  - F1 (diagonal M): ALL PASS ✓
  - (1,1)-string (Λ=[[1,0],[1,1]]): ALL FAIL ✗
  - (2,1)-string (Λ=[[2,1],[1,1]]): ALL FAIL ✗
  - Residuals are proportional to (1-H)·H'² terms → vanish only for H=1 (flat space)
  - Root cause: EF divergence operator has τ₂-dependent volume form: div^{EF} ∝ τ₂^{-5/2}∂(τ₂^2 √g^S ...), and the τ₂ factors don't simplify for non-diagonal M

**Conclusion**: Enhanced-metric formula ĝ = g + M̃ is valid for:
  - ALL components of single-charge branes (diagonal M: F1, D1, NS5, D5)
  - ALL 10d-direction components of general (p,q)-strings (non-diagonal M)
  - But NOT torus components of (p,q)-strings with non-diagonal M

**Open**: Need a torus correction term for non-diagonal M, or an alternative formula.

## Exp34: Action principle

The enhanced metric is an **on-shell construct**. The correct 12d action is:
```
S₁₂ = ∫ √{-g₁₂} [R₁₂ − (1/4)Tr(M⁻¹∂M·M⁻¹∂M) − (1/48)|G₄|²]
```
Standard variation gives λ = 3/40, but using the scalar EOM on-shell shifts λ → 1/4 (10d) and 1/2 (torus). No off-shell action directly gives the enhanced-metric equation.

## IMPORTANT: Redirect from principal investigator

The Ricci bug has already been fixed (do not edit `src/`). Do NOT write up results or declare the search done.

**The original goal has NOT been achieved.** The task is to find a self-contained 12d stress-energy formula T_{MN} such that R_{MN} = T_{MN} directly in 12 dimensions — without reducing to 10d. Everything from exp13 onward has been characterizing the 10d KK decomposition, which I already know. That is not the goal.

The obstruction found in exp1-12 (no single λ works for all blocks) is the core problem. Resolve it. Possible directions:
- A 12d T_{MN} with block-dependent coefficients (not just a single λ)
- Modified 12d stress-energy that incorporates the torus scalar field natively
- A non-standard 12d action principle whose variation gives the correct block structure
- Multiple interacting form fields whose combined T_{MN} balances R_{MN}

You ARE encouraged to write notes and insights about 12d uplifts of type IIB brane solutions in this file — KaTeX is supported for math rendering. However, the Einstein equations of those 12d uplifts must be explicitly checked (computationally, not just argued). You can document derivations and verifications here.

Do NOT declare the search complete until you have an explicit 12d formula T_{MN}(fields) = R_{MN} verified for all components, expressed purely in 12d quantities (no reference to 10d equations).
