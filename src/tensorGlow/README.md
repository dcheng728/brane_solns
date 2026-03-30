# tensorGlow

Abstract tensor algebra package for differential geometry and theoretical physics. Computes with abstract index notation — no coordinates, pure symbolic tensor algebra.

## Goal

A self-contained Python/SymPy package that can:

1. Manipulate tensor expressions with abstract indices (contractions, symmetries, canonicalization)
2. Compute covariant derivatives with Leibniz rule and commutators producing Riemann tensors
3. Compute Christoffel symbols and Riemann tensors from the defining formulas given a metric ansatz
4. Perform dimensional reduction: given a higher-dimensional metric as blocks, derive all lower-dimensional curvature quantities

The primary use case is string theory / supergravity: KK reductions (11d → IIA, F-theory 12d → 10d), conformal transformations, and higher-derivative effective actions like the t8 t8 R^4 computation.

## What's done

### Core algebra (complete)
- **Index system**: `IndexType`, `Index` with `-a` flipping variance, `indices('a b c', L)` factory
- **Tensor declarations**: `TensorHead` with symmetry — `R(a, b, -c, -d)` creates expressions
- **Expression tree**: `TensorAtom`, `TensorProduct`, `TensorSum`, `ScalarExpr` with `+`, `*`, scalar multiplication
- **Symmetry**: `TensorSymmetry` wrapping SymPy's BSGS format — `.riemann()`, `.fully_symmetric()`, `.fully_antisymmetric()`
- **Canonicalization**: Butler-Portugal algorithm via `sympy.combinatorics.tensor_can` — full index remapping, sign extraction

### Metric operations (complete)
- **MetricTensor**: raising/lowering indices, `contract_metrics()` absorbs metrics into other tensors
- Kronecker delta trace: δ_a^a → dim

### Derivatives and curvature (complete)
- **CovDerivative**: Leibniz rule, metric compatibility (D_a g_{bc} = 0), commutator [D_a, D_b] → Riemann terms
- **PartialDerivative**: bare partial derivative without connection — needed for computing Christoffel from the metric
- **RiemannGeometry**: declares Christoffel, Riemann, Ricci, Weyl, Einstein, scalar curvature with correct symmetries
- Weyl decomposition formula, Ricci from Riemann contraction, scalar from Ricci

### Block metric and dimensional reduction (in progress)
- **IndexSplit**: parent IndexType → list of child types (e.g. 11d → 10d + S1)
- **ScalarField**: abstract scalar fields (e.g. dilaton) with chain rule differentiation
- **BlockMetric**: metric expressed as blocks in a KK ansatz — computes Christoffel symbols from the defining formula
  - Tested on 11d = 10d + S1 with dilaton-dependent circle metric
  - Handles scalar field chain rule, truncation (∂_z = 0), Leibniz rule
  - No assumption of block-diagonality — cross-blocks (KK gauge fields) supported

### LaTeX output (complete)
- `to_latex()` with grouped index rendering: T^{ab}_{cd}

## What's next

### Near-term
1. **Riemann from Christoffel** in `BlockMetric` — implement the formula R^A_{BCD} = ∂_C Γ^A_{DB} - ∂_D Γ^A_{CB} + Γ^A_{CE} Γ^E_{DB} - Γ^A_{DE} Γ^E_{CB}, computed block-by-block
2. **Ricci tensor and scalar** from contracting Riemann blocks
3. **Dummy index management** — clean up index names in Christoffel output so formulas are readable
4. **Non-diagonal KK ansatz** — test with KK gauge field C_m in the metric (off-diagonal blocks)

### Medium-term
5. **Conformal transformations** — compute (not substitute) the transformed Riemann under g → e^{2ω} g
6. **Full 11d → type IIA reduction** — reproduce the standard dimensional reduction of the 11d Ricci scalar
7. **F-theory 12d → 10d+T^2** — reproduce equation B.10 of arXiv:1506.06756 (t8 t8 R^4 expansion)

### Longer-term
8. **Hodge duality** with epsilon tensor in abstract notation
9. **Substitution rules and pattern matching** for simplifying with Bianchi identities
10. **Automatic inverse metric** via Schur complement (currently user-provided)

## Usage

```python
from tensorGlow import *
import sympy as sp

# Declare index type and indices
L = IndexType('Lorentz', dim=4, dummy_prefix='L')
a, b, c, d = indices('a b c d', L)

# Metric and curvature
g = MetricTensor(L)
geom = RiemannGeometry(g, L)
R = geom.Riemann
D = geom.D

# Riemann symmetries
R(-a,-b,-d,-c).canon_bp()          # → -R(-a,-b,-c,-d)
(R(-a,-b,-c,-d) + R(-a,-b,-d,-c)).canon_bp()  # → 0

# Covariant derivative
V = TensorHead('V', [L])
D(-a, V(b))                         # → DV(-a, b)
D(-a, g(-b, -c))                    # → 0  (metric compatibility)
D(-c, V(a) * V(b))                  # → DV(-c,a)*V(b) + V(a)*DV(-c,b)

# Commutator
D.commutator(a, b, V(c))            # → R(c,-L_0,-a,-b)*V(L_0)

# Metric contraction
g.contract_metrics(g(-a,-b) * V(b)) # → V(a)

# LaTeX
to_latex(R(a,b,-c,-d) * R(-a,-b,c,d))  # → R^{a b}_{c d}\,R_{a b}^{c d}

# Block metric (KK reduction)
D11 = IndexType('11d', 11, 'M')
D10 = IndexType('10d', 10, 'm')
S1  = IndexType('S1', 1, 'z')
split = IndexSplit(D11, [D10, S1])

tau2 = ScalarField('tau2', depends_on={D10})
bm = BlockMetric(split, blocks={...}, inv_blocks={...},
                 scalar_fields=[tau2], truncation={S1})
christoffel = bm.christoffel()  # all Christoffel blocks from defining formula
```

## Architecture

```
tensorGlow/
├── index.py                 # IndexType, Index, indices()
├── symmetry.py              # TensorSymmetry (BSGS wrapper)
├── tensor_head.py           # TensorHead declarations
├── expr.py                  # TensorAtom, TensorProduct, TensorSum, ScalarExpr
├── canonicalize.py          # Butler-Portugal bridge to SymPy
├── metric.py                # MetricTensor, raising/lowering
├── derivative.py            # CovDerivative, PartialDerivative
├── curvature.py             # RiemannGeometry
├── dimensional_reduction.py # IndexSplit
├── block_metric.py          # BlockMetric, ScalarField
├── latex.py                 # to_latex()
├── utils.py                 # dummy name generation, permutation helpers
└── __init__.py              # public API
```

## Design decisions

- **Own expression tree, SymPy for canonicalization**: SymPy's `tensor.tensor` module is hard to extend. We own the API and expression representation; SymPy's `tensor_can.canonicalize` (Butler-Portugal) does the heavy lifting underneath.
- **No coordinates**: everything is abstract index algebra. For concrete coordinate computations, use the separate `sugra/` package.
- **Computation from defining formulas**: conformal transformations and KK reductions are computed, not substitution rules. The BlockMetric computes Christoffel from Γ = ½ g^{-1}(∂g + ∂g - ∂g).
- **No block-diagonality assumption**: cross-blocks (KK gauge fields, off-diagonal metric components) are fully supported.

## Dependencies

- Python 3.10+
- SymPy 1.14+
