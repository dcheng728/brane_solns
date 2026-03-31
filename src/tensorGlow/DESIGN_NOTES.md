# tensorGlow Design Notes

Lessons from v1 and principles for v2.

## Core principle

**Everything is an operator.** Scalars, tensors, derivatives, coordinates, multiplication — all operators. An expression is a nested composition of operators.

## Key lessons from v1

1. **`*` means "apply to."** Left acts on right. For tensors this is multiplication. For derivatives this is differentiation. For composed actions like `U(b) * D(-b) * V(a)`, evaluation is right-to-left.

2. **D is defined from partial + Gamma, not hardcoded.** Metric compatibility is derived, not assumed. `[D_a, D_b]` produces Riemann by algebraic expansion.

3. **Free indices only for connection terms.** When D acts on an expression, Gamma terms attach to FREE indices only, not internal dummies.

4. **Lazy vs eager expansion.** D should produce an unevaluated node by default. `.expand()` gives `partial + Gamma` terms. This prevents the outer D from seeing internal structure of the inner D.

5. **Coordinates belong to geometry.** Not a separate class.

6. **Every expression should have `.eval()`.** The abstract formula is defined once; `.eval()` evaluates it with specific coordinates. No duplicate implementations.

7. **Two stages:** (a) derive identities abstractly, (b) substitute and evaluate concretely. The same expression tree serves both.

8. **Console IO:** Plain `(a, -b)` notation or LaTeX. Don't try Unicode superscripts.

9. **Butler-Portugal canonicalization** from SymPy works. Group identical tensor heads with `n>1`, pass `msym` as single int, generators have size `rank+2`.

## What needs redesign

- Expression tree: unify everything under "operator applied to operator"
- Remove distinction between abstract geometry and concrete geometry — one class, `.eval()` for concrete
- `derivative.py` is legacy — `operator.py` replaces it but both coexist
- `block_metric.py` and `dimensional_reduction.py` need rethinking in operator framework

## v2 architecture

Start from the operator principle:

```
Operator (base)
  - every operator can be applied to another operator via *
  - every operator has .eval(coords, metric) for concrete evaluation
  - every operator has .free_indices for its index structure

Scalar(value)          — a number, acts by multiplication
Coord(name)            — a coordinate, acts as a scalar field
Partial(coord)         — d/d(coord), acts by differentiation
Tensor(head, indices)  — a named tensor with indices
Metric(coords, components) — a specific metric, defines the geometry
D(metric)              — covariant derivative, defined from Partial + Christoffel(metric)

Expression = nested application of operators
  A * B means "A applied to B"
  (A + B) * C = A*C + B*C (linearity)
```

The Christoffel formula, Riemann tensor, geodesic equation, etc. are all expressions built from these primitives. Abstract when you don't specify a metric. Concrete when you do.
