# Simplifying $t8t8RRDPbDP$
One realizes that
$$
\begin{equation}
\begin{aligned}
t_8 t_8 R^2 |DS^2|
&=\frac{1}{4}{}R^{m n r s} R_{m n r s} S^{t u v w} \bar{S}{}_{t u v w}
+R^{m n r s} R_{m n}{}^{t u} S_{r s}{}^{v w} \bar{S}{}_{t u v w}\\
&-4{}R^{m n r s} R_{m n r}{}^{t} S_{s}{}^{u v w} \bar{S}{}_{t u v w}
-2{}R^{m n r s} R_{m n}{}^{t u} S_{r t}{}^{v w} \bar{S}{}_{s u v w}\\
&+\frac{1}{2}{}R^{m n r s} R^{t u v w} S_{m n r s} \bar{S}{}_{t u v w}
+\frac{1}{2}{}R^{m n r s} R^{t u v w} S_{m n t u} \bar{S}{}_{r s v w}\\
&-4{}R^{m n r s} R_{m}{}^{t u v} S_{n}{}^{w}{}_{u v} \bar{S}{}_{r s t w}
-4{}R^{m n r s} R_{m}{}^{t u v} S_{n}{}^{w}{}_{r s} \bar{S}{}_{t w u v}\\
&-4{}R^{m n r s} R^{t u v w} S_{m n r t} \bar{S}{}_{s u v w}
+8{}R^{m n r s} R_{m}{}^{t}{}_{r}{}^{u} S_{n}{}^{v}{}_{s}{}^{w} \bar{S}{}_{t v u w}\\
&+8{}R^{m n r s} R_{m}{}^{t}{}_{r}{}^{u} S_{n}{}^{v}{}_{u}{}^{w} \bar{S}{}_{s w t v}
+16{}R^{m n r s} R_{m}{}^{t u v} S_{n}{}^{w}{}_{r u} \bar{S}{}_{s v t w}\\
&+4{}R^{m n r s} R^{t u v w} S_{m t r v} \bar{S}{}_{n u s w}
\end{aligned}
\end{equation}
$$

Letting $S_{mnrs} = DP_{mnrs}$, Using Cadabra, we arrive at, here we are using $K_{mn}$ to denote $D_m P_n$, which is symmetric in $m,n$, same for $\bar{K}_{mn}$.

$$
\begin{equation}
\begin{aligned}
t_8 t_8 R^2 |DP^2|
&=\frac{1}{4}{}R^{m n r s} R_{m n r s} K^{t u} \bar{K}{}_{t u}
+\frac{1}{16}{}R^{m n r s} R_{m n r s} K^{t}{}_{t} \bar{K}{}^{u}{}_{u}\\
&-2{}R^{m n r s} R_{m n r}{}^{t} K_{s}{}^{u} \bar{K}{}_{t u}
-R^{m n r s} R_{m n}{}^{t u} K_{r t} \bar{K}{}_{s u}\\
&-\frac{1}{4}{}R^{m n r s} R_{m n r}{}^{t} K^{u}{}_{u} \bar{K}{}_{s t}
-\frac{1}{4}{}R^{m n r s} R_{m n r}{}^{t} K_{s t} \bar{K}{}^{u}{}_{u}\\
&+R^{m n}{}_{m}{}^{r} R^{s t}{}_{s}{}^{u} K_{n r} \bar{K}{}_{t u}
+6{}R^{m n r s} R_{m}{}^{t}{}_{r}{}^{u} K_{n t} \bar{K}{}_{s u}\\
&-2{}R^{m n}{}_{m}{}^{r} R_{n}{}^{s t u} K_{r t} \bar{K}{}_{s u}
+2{}R^{m n}{}_{m}{}^{r} R_{n}{}^{s t u} K_{s t} \bar{K}{}_{r u}\\
&+2{}R^{m n r s} R_{m}{}^{t}{}_{r}{}^{u} K_{n s} \bar{K}{}_{t u}
-R^{m n}{}_{m}{}^{r} R^{s t}{}_{s}{}^{u} K_{n t} \bar{K}{}_{r u}\\
&+\frac{1}{2}{}R^{m n}{}_{m}{}^{r} R_{n}{}^{s}{}_{r s} K^{t u} \bar{K}{}_{t u}
-2{}R^{m n}{}_{m}{}^{r} R_{n}{}^{s}{}_{r}{}^{t} K_{s}{}^{u} \bar{K}{}_{t u}\\
&+\frac{1}{2}{}R^{m n}{}_{m}{}^{r} R_{n}{}^{s}{}_{r}{}^{t} K^{u}{}_{u} \bar{K}{}_{s t}
+\frac{1}{2}{}R^{m n}{}_{m}{}^{r} R_{n}{}^{s}{}_{r}{}^{t} K_{s t} \bar{K}{}^{u}{}_{u}\\
&+\frac{1}{2}{}R^{m n}{}_{m n} R^{r s t u} K_{r t} \bar{K}{}_{s u}
\end{aligned}
\end{equation}
$$

Sorting

$$
\begin{equation}
\begin{aligned}
t_8 t_8 R^2 |DP^2|
%--- 4 free indices on K \bar{K} ---
&=-{}R^{m n r s} R_{m n}{}^{t u} K_{r t} \bar{K}{}_{s u}
+R^{m n}{}_{m}{}^{r} R^{s t}{}_{s}{}^{u} K_{n r} \bar{K}{}_{t u}\\
&+6{}R^{m n r s} R_{m}{}^{t}{}_{r}{}^{u} K_{n t} \bar{K}{}_{s u}
-2{}R^{m n}{}_{m}{}^{r} R_{n}{}^{s t u} K_{r t} \bar{K}{}_{s u}\\
&+2{}R^{m n}{}_{m}{}^{r} R_{n}{}^{s t u} K_{s t} \bar{K}{}_{r u}
+\cancel{2{}R^{m n r s} R_{m}{}^{t}{}_{r}{}^{u} K_{n s} \bar{K}{}_{t u}}\\
&-R^{m n}{}_{m}{}^{r} R^{s t}{}_{s}{}^{u} K_{n t} \bar{K}{}_{r u}
+\frac{1}{2}{}R^{m n}{}_{m n} R^{r s t u} K_{r t} \bar{K}{}_{s u}\\
%--- 2 free indices on K \bar{K} ---
&-2{}R^{m n r s} R_{m n r}{}^{t} K_{s}{}^{u} \bar{K}{}_{t u}
-\frac{1}{4}{}R^{m n r s} R_{m n r}{}^{t} K^{u}{}_{u} \bar{K}{}_{s t}\\
&-\frac{1}{4}{}R^{m n r s} R_{m n r}{}^{t} K_{s t} \bar{K}{}^{u}{}_{u}
-2{}R^{m n}{}_{m}{}^{r} R_{n}{}^{s}{}_{r}{}^{t} K_{s}{}^{u} \bar{K}{}_{t u}\\
&+\frac{1}{2}{}R^{m n}{}_{m}{}^{r} R_{n}{}^{s}{}_{r}{}^{t} K^{u}{}_{u} \bar{K}{}_{s t}
+\frac{1}{2}{}R^{m n}{}_{m}{}^{r} R_{n}{}^{s}{}_{r}{}^{t} K_{s t} \bar{K}{}^{u}{}_{u}\\
%--- 0 free indices on K \bar{K} ---
&+\frac{1}{4}{}R^{m n r s} R_{m n r s} K^{t u} \bar{K}{}_{t u}
+\frac{1}{16}{}R^{m n r s} R_{m n r s} K^{t}{}_{t} \bar{K}{}^{u}{}_{u}\\
&+\frac{1}{2}{}R^{m n}{}_{m}{}^{r} R_{n}{}^{s}{}_{r s} K^{t u} \bar{K}{}_{t u}
\end{aligned}
\end{equation}
$$

One can show that define $K^{mnrs} \equiv D^m\bar{P}^n D^r P^s$.
For any contraction with a real tensor $T$, one has the effective symmetry

$$
\begin{equation}
    K^{mnrs} \cong K^{rsmn}.
\end{equation}
$$

Proof is done by conjugation.
It follows that

$$
\begin{equation}
    K^{mnrs}R_{m}{}^{t}{}_{n}{}^{u}{}R_{rtsu} = 0.
\end{equation}
$$
Apply the first Bianchi identity $R_{r[tsu]}=0$ to the second Riemann factor.

Furthermore,

$$
\begin{equation}
    K^{mnrs}{}R_{mr}{}^{tu}{}R_{nstu}
    = 2{}K^{mnrs}{}R_{m}{}^{t}{}_{r}{}^{u}{}R_{nstu}.
\end{equation}
$$

Proof by First Bianchi identity.

From the properties we have derived, we know that $2{}R^{m n r s} R_{m}{}^{t}{}_{r}{}^{u} K_{n s} \bar{K}{}_{t u} = 0$.

## Agent Simplification — Step 1
*Generated: 2026-04-12 19:29:01*

### Proposal

Looking at the current expression term by term, I'll label them in order of appearance (skipping the already-cancelled T6):

- T10: $-\frac{1}{4}R^{mnrs}R_{mnr}{}^{t}\,K^{u}{}_{u}\,\bar{K}_{st}$
- T11: $-\frac{1}{4}R^{mnrs}R_{mnr}{}^{t}\,K_{st}\,\bar{K}^{u}{}_{u}$

These two have identical Riemann structure but swapped $K \leftrightarrow \bar{K}$ roles.

```
SIMPLIFICATION_FOUND

TARGET: T10 and T11

CLAIM: T10 = T11 by the K-swap symmetry, so T10 + T11 = 2 T11

DERIVATION:

1. Write T10 in terms of the composite object:
   T10 = -1/4 · R^{mnrs} R_{mnr}^t · bar{K}_{st} K^u_u
       = -1/4 · K^{abcd} T_{abcd}
   where K^{abcd} = bar{K}^{ab} K^{cd}, and
   T_{abcd} = R^{mnr}{}_a R_{mnrb} delta_{cd}
   (here a=s, b=t carry the Riemann contraction; c=d are traced).

2. T_{abcd} is real (built entirely from the real metric's Riemann tensor
   and the metric delta). This is the condition needed to apply the
   K-swap symmetry.

3. Apply the effective symmetry K^{abcd} T_{abcd} = K^{cdab} T_{abcd}:
   K^{cdab} T_{abcd} = bar{K}^{cd} K^{ab} · R^{mnr}{}_a R_{mnrb} · delta_{cd}
                      = K^{ab} R^{mnr}{}_a R_{mnrb} · bar{K}^c_c
                      = R^{mnrs} R_{mnr}^t K_{st} bar{K}^u_u

4. Therefore:
   T10 = -1/4 · R^{mnrs} R_{mnr}^t · K_{st} bar{K}^u_u = T11

5. T10 + T11 = 2 T11 = -1/2 · R^{mnrs} R_{mnr}^t K_{st} bar{K}^u_u

SIMPLIFIED_EXPRESSION:
t_8 t_8 R^2 |DP^2|
  = - R^{mnrs} R_{mn}^{tu} K_{rt} bar{K}_{su}                         [T1]
  + R^{mn}_{m}^{r} R^{st}_{s}^{u} K_{nr} bar{K}_{tu}                  [T2]
  + 6 R^{mnrs} R_{m}^{t}_{r}^{u} K_{nt} bar{K}_{su}                   [T3]
  - 2 R^{mn}_{m}^{r} R_{n}^{stu} K_{rt} bar{K}_{su}                   [T4]
  + 2 R^{mn}_{m}^{r} R_{n}^{stu} K_{st} bar{K}_{ru}                   [T5]
  - R^{mn}_{m}^{r} R^{st}_{s}^{u} K_{nt} bar{K}_{ru}                  [T7]
  + 1/2 R^{mn}_{mn} R^{rstu} K_{rt} bar{K}_{su}                       [T8]
  - 2 R^{mnrs} R_{mnr}^{t} K_{s}^{u} bar{K}_{tu}                      [T9]
  - 1/2 R^{mnrs} R_{mnr}^{t} K_{st} bar{K}^{u}_{u}                    [T10+T11]
  - 2 R^{mn}_{m}^{r} R_{n}^{s}_{r}^{t} K_{s}^{u} bar{K}_{tu}         [T12]
  + 1/2 R^{mn}_{m}^{r} R_{n}^{s}_{r}^{t} K^{u}_{u} bar{K}_{st}       [T13]
  + 1/2 R^{mn}_{m}^{r} R_{n}^{s}_{r}^{t} K_{st} bar{K}^{u}_{u}       [T14]
  + 1/4 R^{mnrs} R_{mnrs} K^{tu} bar{K}_{tu}                          [T15]
  + 1/16 R^{mnrs} R_{mnrs} K^{t}_{t} bar{K}^{u}_{u}                   [T16]
  + 1/2 R^{mn}_{m}^{r} R_{n}^{s}_{rs} K^{tu} bar{K}_{tu}              [T17]
```

### Verification

```
VERDICT: ACCEPT

CHECK:

Step 1 (Writing T10 in terms of the composite object):
  T10 = -1/4 · R^{mnrs} R_{mnr}^t · bar{K}_{st} · K^u_u

  Define S_{ab} = R^{mnr}_{a} R_{mnrb} (a real rank-2 tensor built from Riemann).
  Define T_{abcd} = S_{ab} g_{cd}.

  Then K^{abcd} T_{abcd} = bar{K}^{ab} K^{cd} S_{ab} g_{cd}
                          = (bar{K}^{ab} S_{ab})(K^c_c)
                          = S^{st} bar{K}_{st} · K^u_u.

  So T10 = -1/4 · K^{abcd} T_{abcd}. CORRECT.

Step 2 (T_{abcd} is real):
  T_{abcd} = R^{mnr}_{a} R_{mnrb} g_{cd} is built entirely from the
  Riemann tensor and the metric, both real. CORRECT.

Step 3 (Applying the K-swap symmetry):
  The effective symmetry gives K^{abcd} T_{abcd} = K^{cdab} T_{abcd}.

  K^{cdab} T_{abcd} = bar{K}^{cd} K^{ab} · S_{ab} g_{cd}
                     = (K^{ab} S_{ab})(bar{K}^{cd} g_{cd})
                     = S^{st} K_{st} · bar{K}^u_u.

  CORRECT.

Step 4 (Concluding T10 = T11):
  From step 3: -1/4 S^{st} bar{K}_{st} K^u_u = -1/4 S^{st} K_{st} bar{K}^u_u.
  The RHS is exactly T11 = -1/4 R^{mnrs} R_{mnr}^t K_{st} bar{K}^u_u.
  CORRECT.

Step 5 (Combined coefficient):
  T10 + T11 = 2 · T11 = -1/2 R^{mnrs} R_{mnr}^t K_{st} bar{K}^u_u.
  Arithmetic: 2 × (-1/4) = -1/2. CORRECT.

Final expression check:
  All 15 remaining terms (T1–T5, T7–T9, T10+T11, T12–T17) are present
  with correct coefficients and index structures, matching the sorted
  expression with T6 cancelled and T10+T11 combined. Term count: 15 lines
  in the output, consistent with 17 original terms minus T6 (cancelled)
  minus one line from merging T10+T11. CORRECT.
```

---

## Agent Simplification — Step 2
*Generated: 2026-04-12 19:30:44*

### Proposal

Looking at the current expression, I'll check for the same type of K-swap symmetry that combined T10+T11, applied to another pair.

```
SIMPLIFICATION_FOUND

TARGET: T13 and T14

CLAIM: T13 = T14 by the K-swap symmetry, so T13 + T14 = 2 T14

DERIVATION:

1. Write T13 explicitly:
   T13 = 1/2 · R^{mn}_{m}^{r} R_{n}^{s}_{r}^{t} · K^{u}_{u} · bar{K}_{st}

2. Define the real rank-2 tensor S_{ab} = R^{mn}_{m a} R_{n}^{s}_{a}^{t}
   evaluated as S_{st} = R^{mn}_{m}^{r} R_{n s r t} (contracted with
   the metric to raise/lower as needed). More precisely, define
   T_{abcd} = S_{ab} g_{cd}, where S_{ab} = R^{mn}_{ma} R_{n}^{e}_{ae} delta^...

   Let me be more careful. Define:
   T_{abcd} = R^{mn}_{m}^{r} R_{narb} · g_{cd}

   This is real (built from Riemann and the metric).

3. Compute K^{abcd} T_{abcd}:
   K^{abcd} T_{abcd} = bar{K}^{ab} K^{cd} · R^{mn}_{m}^{r} R_{narb} · g_{cd}
                      = bar{K}^{ab} · R^{mn}_{m}^{r} R_{narb} · K^{c}_{c}

   Identifying a = s, b = t: this gives
   R^{mn}_{m}^{r} R_{n}^{s}_{r}^{t} · bar{K}_{st} · K^{u}_{u}
   which (up to the 1/2 prefactor) is T13.

4. Apply the effective symmetry K^{abcd} T_{abcd} = K^{cdab} T_{abcd}:
   K^{cdab} T_{abcd} = bar{K}^{cd} K^{ab} · R^{mn}_{m}^{r} R_{narb} · g_{cd}
                      = K^{ab} · R^{mn}_{m}^{r} R_{narb} · bar{K}^{c}_{c}

   This gives R^{mn}_{m}^{r} R_{n}^{s}_{r}^{t} · K_{st} · bar{K}^{u}_{u},
   which (up to the 1/2 prefactor) is T14.

5. Therefore T13 = T14, and:
   T13 + T14 = 2 · T14 = R^{mn}_{m}^{r} R_{n}^{s}_{r}^{t} K_{st} bar{K}^{u}_{u}

SIMPLIFIED_EXPRESSION:
t_8 t_8 R^2 |DP^2|
  = - R^{mnrs} R_{mn}^{tu} K_{rt} bar{K}_{su}                         [T1]
  + R^{mn}_{m}^{r} R^{st}_{s}^{u} K_{nr} bar{K}_{tu}                  [T2]
  + 6 R^{mnrs} R_{m}^{t}_{r}^{u} K_{nt} bar{K}_{su}                   [T3]
  - 2 R^{mn}_{m}^{r} R_{n}^{stu} K_{rt} bar{K}_{su}                   [T4]
  + 2 R^{mn}_{m}^{r} R_{n}^{stu} K_{st} bar{K}_{ru}                   [T5]
  - R^{mn}_{m}^{r} R^{st}_{s}^{u} K_{nt} bar{K}_{ru}                  [T7]
  + 1/2 R^{mn}_{mn} R^{rstu} K_{rt} bar{K}_{su}                       [T8]
  - 2 R^{mnrs} R_{mnr}^{t} K_{s}^{u} bar{K}_{tu}                      [T9]
  - 1/2 R^{mnrs} R_{mnr}^{t} K_{st} bar{K}^{u}_{u}                    [T10+T11]
  - 2 R^{mn}_{m}^{r} R_{n}^{s}_{r}^{t} K_{s}^{u} bar{K}_{tu}         [T12]
  + R^{mn}_{m}^{r} R_{n}^{s}_{r}^{t} K_{st} bar{K}^{u}_{u}           [T13+T14]
  + 1/4 R^{mnrs} R_{mnrs} K^{tu} bar{K}_{tu}                          [T15]
  + 1/16 R^{mnrs} R_{mnrs} K^{t}_{t} bar{K}^{u}_{u}                   [T16]
  + 1/2 R^{mn}_{m}^{r} R_{n}^{s}_{rs} K^{tu} bar{K}_{tu}              [T17]
```

### Verification

```
VERDICT: ACCEPT

CHECK:

Step 1 (Writing T13 in terms of the composite object):
  T13 = 1/2 · R^{mn}_{m}^{r} R_{n}^{s}_{r}^{t} · K^{u}_{u} · bar{K}_{st}

  Rewrite using index lowering: R_{n}^{s}_{r}^{t} bar{K}_{st} = R_{narb} bar{K}^{ab}
  (since R_{n}^{s}_{r}^{t} = g^{sa} g^{tb} R_{narb}, contracting with bar{K}_{st}
  gives R_{narb} bar{K}^{ab}).

  Define S_{ab} = R^{mn}_{m}^{r} R_{narb}. Then:
  T13 = 1/2 · S_{ab} bar{K}^{ab} · K^{c}_{c}
      = 1/2 · K^{abcd} T_{abcd}
  where K^{abcd} = bar{K}^{ab} K^{cd} and T_{abcd} = S_{ab} g_{cd}.

  Verification: K^{abcd} T_{abcd} = bar{K}^{ab} K^{cd} S_{ab} g_{cd}
                                   = (bar{K}^{ab} S_{ab})(K^{c}_{c})
                                   = S_{ab} bar{K}^{ab} · K^{u}_{u}.   CORRECT.

Step 2 (T_{abcd} is real):
  S_{ab} = R^{mn}_{m}^{r} R_{narb} is built entirely from the Riemann tensor
  and metric. T_{abcd} = S_{ab} g_{cd} likewise. Both are real.   CORRECT.

Step 3 (Applying the K-swap symmetry):
  K^{abcd} T_{abcd} = K^{cdab} T_{abcd}.

  K^{cdab} T_{abcd} = bar{K}^{cd} K^{ab} · S_{ab} · g_{cd}
                     = (K^{ab} S_{ab})(bar{K}^{c}_{c}).

  Now K^{ab} S_{ab} = K^{ab} R^{mn}_{m}^{r} R_{narb}.
  Since R_{narb} K^{ab} = R_{narb} g_{aa'} g_{bb'} K^{a'b'}... more directly:
  contracting K^{ab} with R_{narb} in positions (2,4) gives
  R_{n}^{s}_{r}^{t} K_{st} (relabeling dummy indices).

  So K^{cdab} T_{abcd} = R^{mn}_{m}^{r} R_{n}^{s}_{r}^{t} K_{st} · bar{K}^{u}_{u}.
  This is exactly T14 (without the 1/2 prefactor).   CORRECT.

Step 4 (Concluding T13 = T14):
  From steps 1 and 3: both T13 and T14 equal 1/2 · K^{abcd} T_{abcd}
  (one via the original pairing, one via the swapped pairing, which are
  equal by the K-swap symmetry). Therefore T13 = T14.   CORRECT.

Step 5 (Combined coefficient):
  T13 + T14 = 2 · (1/2) · R^{mn}_{m}^{r} R_{n}^{s}_{r}^{t} K_{st} bar{K}^{u}_{u}
            = 1 · R^{mn}_{m}^{r} R_{n}^{s}_{r}^{t} K_{st} bar{K}^{u}_{u}.
  Arithmetic: 2 × (1/2) = 1.   CORRECT.

Final expression check:
  Starting from the 15-term expression (post T10+T11 merger), merging T13+T14
  reduces to 14 terms. The proposed output lists: T1, T2, T3, T4, T5, T7, T8,
  T9, T10+T11, T12, T13+T14, T15, T16, T17 — that is 14 lines. All coefficients
  and index structures match the input expression, with T13+T14 correctly shown
  with coefficient 1.   CORRECT.
```

---

## Agent Simplification — Step 1
*Generated: 2026-04-12 19:37:29*

### Proposal

Looking at the current 14-term expression, I'll check whether T4 and T5 are related by the $\mathcal{K}$-swap symmetry.

```
SIMPLIFICATION_FOUND

TARGET: T4 and T5

CLAIM: T4 = T5 by K-swap symmetry plus dummy index relabeling, so T4 + T5 = 2 T4

DERIVATION:

1. Write out the two terms:
   T4 = -2 R^{mn}_{m}^{r} R_{n}^{stu} K_{rt} bar{K}_{su}
   T5 = +2 R^{mn}_{m}^{r} R_{n}^{stu} K_{st} bar{K}_{ru}

2. Define C = R^{mn}_{m}^{r} R_{n}^{stu} K_{rt} bar{K}_{su}, so T4 = -2C.
   Write C in composite notation: C = K^{abcd} T_{abcd} where
   K^{abcd} = bar{K}^{ab} K^{cd} with a=s, b=u, c=r, d=t,
   and T_{surt} = R^{mn}_{mr} R_{nsut} (all-lower form of the
   Riemann bilinear). T_{surt} is real, being built entirely
   from the Riemann tensor.

3. Apply the K-swap symmetry K^{abcd} T_{abcd} = K^{cdab} T_{abcd}:
   C = R^{mn}_{m}^{r} R_{n}^{stu} bar{K}_{rt} K_{su}

4. Relabel dummy indices t <-> u (both are fully contracted):
   C = R^{mn}_{m}^{r} R_{n}^{sut} bar{K}_{ru} K_{st}

5. Apply Riemann antisymmetry R_{abcd} = -R_{abdc} in the last
   two indices: R_{n}^{sut} = -R_{n}^{stu}.
   Therefore:
   C = -R^{mn}_{m}^{r} R_{n}^{stu} K_{st} bar{K}_{ru}

6. This means R^{mn}_{m}^{r} R_{n}^{stu} K_{st} bar{K}_{ru} = -C.
   Therefore T5 = +2(-C) = -2C = T4.

7. Combined: T4 + T5 = -2C + (-2C) = -4C
   = -4 R^{mn}_{m}^{r} R_{n}^{stu} K_{rt} bar{K}_{su}
   Arithmetic: 2 × (-2) = -4. CORRECT.

SIMPLIFIED_EXPRESSION:
t_8 t_8 R^2 |DP^2|
  = - R^{mnrs} R_{mn}^{tu} K_{rt} bar{K}_{su}                         [T1]
  + R^{mn}_{m}^{r} R^{st}_{s}^{u} K_{nr} bar{K}_{tu}                  [T2]
  + 6 R^{mnrs} R_{m}^{t}_{r}^{u} K_{nt} bar{K}_{su}                   [T3]
  - 4 R^{mn}_{m}^{r} R_{n}^{stu} K_{rt} bar{K}_{su}                   [T4+T5]
  - R^{mn}_{m}^{r} R^{st}_{s}^{u} K_{nt} bar{K}_{ru}                  [T7]
  + 1/2 R^{mn}_{mn} R^{rstu} K_{rt} bar{K}_{su}                       [T8]
  - 2 R^{mnrs} R_{mnr}^{t} K_{s}^{u} bar{K}_{tu}                      [T9]
  - 1/2 R^{mnrs} R_{mnr}^{t} K_{st} bar{K}^{u}_{u}                    [T10+T11]
  - 2 R^{mn}_{m}^{r} R_{n}^{s}_{r}^{t} K_{s}^{u} bar{K}_{tu}         [T12]
  + R^{mn}_{m}^{r} R_{n}^{s}_{r}^{t} K_{st} bar{K}^{u}_{u}           [T13+T14]
  + 1/4 R^{mnrs} R_{mnrs} K^{tu} bar{K}_{tu}                          [T15]
  + 1/16 R^{mnrs} R_{mnrs} K^{t}_{t} bar{K}^{u}_{u}                   [T16]
  + 1/2 R^{mn}_{m}^{r} R_{n}^{s}_{rs} K^{tu} bar{K}_{tu}              [T17]
```

### Verification

```
VERDICT: ACCEPT

CHECK:

Step 1 (Writing out T4 and T5): Matches the 14-term expression from the previous step. CORRECT.

Step 2 (Defining C and composite notation):
  C = R^{mn}_{m}^{r} R_{n}^{stu} K_{rt} bar{K}_{su}, so T4 = -2C.
  In composite notation: K^{abcd} = bar{K}^{ab} K^{cd} with a=s, b=u, c=r, d=t.
  T_{surt} is built from the Riemann bilinear R^{mn}_{mr} R_{nsut}, which is real.
  CORRECT.

Step 3 (K-swap symmetry):
  The symmetry K^{surt} T_{surt} = K^{rtsu} T_{surt} swaps which pair
  carries bar{K} vs K:
    Original: bar{K}_{su} K_{rt} · (Riemann structure)
    Swapped:  bar{K}_{rt} K_{su} · (same Riemann structure)
  So C = R^{mn}_{m}^{r} R_{n}^{stu} bar{K}_{rt} K_{su}. CORRECT.

Step 4 (Dummy relabeling t <-> u):
  Both t and u are fully contracted, so the relabeling is a valid bijection.
  Free indices: none beyond m,n,r,s,t,u which are all summed.
  After t <-> u:
    bar{K}_{rt} -> bar{K}_{ru}
    K_{su} -> K_{st}
    R_{n}^{stu} -> R_{n}^{sut}  (t in slot 3 becomes u, u in slot 4 becomes t)
  Result: C = R^{mn}_{m}^{r} R_{n}^{sut} K_{st} bar{K}_{ru}. CORRECT.

Step 5 (Riemann antisymmetry in last two indices):
  R_{abcd} = -R_{abdc}, so R_{n}^{sut} = -R_{n}^{stu}
  (slots 3 and 4 are swapped, picking up a minus sign).
  Therefore: C = -R^{mn}_{m}^{r} R_{n}^{stu} K_{st} bar{K}_{ru}. CORRECT.

Step 6 (Concluding T5 = T4):
  From step 5: R^{mn}_{m}^{r} R_{n}^{stu} K_{st} bar{K}_{ru} = -C.
  T5 = +2 · R^{mn}_{m}^{r} R_{n}^{stu} K_{st} bar{K}_{ru} = +2(-C) = -2C.
  T4 = -2C.
  So T5 = T4. CORRECT.

Step 7 (Combined coefficient):
  T4 + T5 = -2C + (-2C) = -4C
  = -4 R^{mn}_{m}^{r} R_{n}^{stu} K_{rt} bar{K}_{su}.
  Arithmetic: (-2) + (-2) = -4. CORRECT.

Final expression check:
  Starting from the 14-term expression, merging T4+T5 reduces to 13 terms.
  The output lists: T1, T2, T3, T4+T5, T7, T8, T9, T10+T11, T12, T13+T14,
  T15, T16, T17 — 13 lines. The T4+T5 coefficient is -4, matching -2C-2C.
  All other terms carry over unchanged with correct coefficients and index
  structures. CORRECT.
```

---
