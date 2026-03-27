# sugra — Symbolic Supergravity Solution Verifier

A Python package for defining supergravity solutions symbolically and checking whether they satisfy their equations of motion.

## What it does

Given a brane-like ansatz (a metric, a dilaton, and flux fields), `sugra` will:

1. Compute the Ricci tensor and scalar curvature from the metric
2. Compute the stress-energy tensor from the flux fields
3. Check that the Einstein equations, Bianchi identities, dilaton equation, and self-duality constraints are satisfied
4. Report which checks passed and which failed, using both symbolic simplification and numerical spot-checks as fallbacks

## Modules

### `geometry.py` — Metrics and curvature

**`Metric`** — Represents a pseudo-Riemannian metric. Given a list of diagonal components (as SymPy expressions) and a list of coordinate symbols, it computes:
- Christoffel symbols Γ^ρ_{μν}
- Ricci tensor R_{μν}
- Ricci scalar R

Diagonal metrics are handled with an optimized O(D²) algorithm. Expensive computations are cached.

**`HarmonicFunction`** — Represents the harmonic function H(r) that appears in brane solutions. It enforces the harmonic condition:

```
H'' + (D_perp - 1)/r · H' = 0
```

where `D_perp` is the number of transverse dimensions. Used to simplify expressions by substituting H'' → -(D_perp-1)/r · H'.

**`warped_product(coords, blocks)`** — Constructs a warped product metric. Each block is `(metric_components, warp_factor)`, producing:

```
ds² = H^a₁ ds₁² + H^a₂ ds₂² + ...
```

---

### `forms.py` — Differential form fields

**`FormField(p, D, components)`** — A sparse antisymmetric p-form on a D-dimensional spacetime. Only independent components are stored; antisymmetry is handled automatically when you access a permuted index.

Key operations (all take a `Metric` and return a new `FormField` or SymPy expression):

| Function | Output |
|---|---|
| `exterior_derivative(F, coords)` | dF (a (p+1)-form) |
| `hodge_star(F, metric)` | ⋆F (the Hodge dual) |
| `form_norm_squared(F, metric)` | \|F\|² = (1/p!) F_{M₁…Mₚ} F^{M₁…Mₚ} |
| `form_contraction(F, metric)` | F_{M P₁…} F_N^{P₁…} (used in stress-energy) |
| `form_stress_energy(F, metric, p)` | T_{MN} contribution from a p-form |

---

### `verify.py` — Equations of motion checks

**`ScalarField(expr, symbol)`** — Wraps a scalar field (e.g. the dilaton φ).

**`FluxField(form, dilaton_coupling, self_dual)`** — Wraps a p-form flux, with an optional dilaton coupling e^{αφ} and a flag for self-duality constraints.

**`Solution(metric, scalars, fluxes, harmonic)`** — The central object. Takes a `Metric`, a list of `ScalarField`s, a list of `FluxField`s, and optionally a `HarmonicFunction`. Provides:

| Method | Checks |
|---|---|
| `check_einstein()` | R_{MN} = T_{MN} |
| `check_bianchi()` | dF = 0 for each flux |
| `check_self_duality()` | F = ⋆F for self-dual fluxes |
| `check_dilaton()` | Dilaton equation of motion |
| `verify_all()` | All of the above, returns a `VerificationReport` |

**`VerificationReport`** — A collection of `CheckResult` objects. Prints a table showing which checks passed, which failed, and by what method (symbolic or numerical).

#### How verification works

For each expression that should vanish, `check_expression()` tries:
1. Fast symbolic simplification (`cancel`, `factor`, `radsimp`)
2. Full `simplify()` (slower)
3. Numerical evaluation at random points consistent with the harmonic condition

If any method confirms vanishing, the check passes.

---

## Typical workflow

```python
import sympy as sp
from sugra import Metric, HarmonicFunction, warped_product, FormField, FluxField, ScalarField, Solution

# 1. Define coordinates and the harmonic function
r, H_sym = sp.symbols('r H', positive=True)
coords = [t, x1, x2, x3, x4, x5, x6, x7, x8, x9]  # 10d coordinates

harmonic = HarmonicFunction(H_sym, r, D_perp=3)

# 2. Build the metric as a warped product
H = harmonic.H
g_components = [...]   # list of diagonal metric entries
metric = Metric(g_components, coords)

# 3. Define flux fields
F_components = {...}   # {(μ,ν,...): expr} sparse components
F = FormField(p=5, D=10, components=F_components)
flux = FluxField(form=F, dilaton_coupling=0, self_dual=True)

# 4. Define the dilaton
phi = ScalarField(expr=sp.log(H)/4, symbol=sp.Symbol('phi'))

# 5. Build and verify the solution
sol = Solution(metric=metric, scalars=[phi], fluxes=[flux], harmonic=harmonic)
report = sol.verify_all()
print(report)
```

---

## Dependencies

- [SymPy](https://www.sympy.org/) — all symbolic computation
- Python standard library (`itertools`, `dataclasses`, `functools`, `random`)
