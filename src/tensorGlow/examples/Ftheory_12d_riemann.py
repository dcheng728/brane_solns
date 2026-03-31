"""
Riemann tensor of the 12d F-theory metric, computed from first principles.

The 12d metric is block-diagonal:

    G_MN = g_{mn}  oplus  M_{ab}

where g_{mn} is the 10d SL(2,R)-invariant Einstein frame metric
and M_{ab} = (1/tau_2) * [[1, tau_1], [tau_1, |tau|^2]] is the
axio-dilaton matrix on T^2.

We split 12d -> 10d + T^2, with truncation d_a = 0 on all fields.

BlockMetric computes everything from the defining formulas:
  1. Christoffel symbol blocks  Gamma^A_{BC}
  2. Riemann tensor blocks      R_{ABCD} (all lower)
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from tensorGlow import *
from tensorGlow.core.expr import replace_index


# ═══════════════════════════════════════════════════════════════
# Index types
# ═══════════════════════════════════════════════════════════════

L = IndexType('10d', dim=10, dummy_prefix='L')
T = IndexType('T2', dim=2, dummy_prefix='T')


# ═══════════════════════════════════════════════════════════════
# Abstract metrics
# ═══════════════════════════════════════════════════════════════

g = MetricTensor(L, name='g')      # 10d Einstein-frame metric
M = MetricTensor(T, name='M')      # T^2 axio-dilaton matrix


# ═══════════════════════════════════════════════════════════════
# 12d block metric with truncation d_a = 0
# ═══════════════════════════════════════════════════════════════

parent = IndexType('12d', dim=12, dummy_prefix='P')
split = IndexSplit(parent, [L, T])

bm = BlockMetric(
    split,
    blocks     = {(L, L): g,     (T, T): M},
    inv_blocks = {(L, L): g.inv, (T, T): M.inv},
    truncation = {T},
)


# ═══════════════════════════════════════════════════════════════
# Relabeling: canonical indices → readable names
# ═══════════════════════════════════════════════════════════════

def _name_gen(start_names):
    """Yield names from start_names, then prefix_0, prefix_1, ..."""
    for n in start_names:
        yield n
    i = 0
    while True:
        yield f'{start_names[0]}{i}'
        i += 1


def relabel(expr):
    """Replace all canonical index names (_ca, _ra, _g0, ...) with
    readable physics indices (m, n, p, q, ... for 10d; a, b, c, d, ... for T²).

    Free indices are assigned first (in order of appearance),
    then dummies get subsequent names.
    """
    pools = {
        id(L): _name_gen('mnpqrstuvw'),
        id(T): _name_gen('abcdef'),
    }
    assigned = {}   # (name, type_id) → readable name

    def assign(name, itype):
        key = (name, id(itype))
        if key not in assigned:
            assigned[key] = next(pools[id(itype)])

    # 1. Assign free indices first (ensures they get primary names)
    for idx in expr.free_indices:
        assign(idx.name, idx.index_type)

    # 2. Walk the expression tree to assign remaining (dummy) names
    def visit(e):
        if isinstance(e, TensorAtom):
            for idx in e.indices:
                assign(idx.name, idx.index_type)
        elif isinstance(e, TensorProduct):
            for f in e.factors:
                visit(f)
        elif isinstance(e, TensorSum):
            for t in e.terms:
                visit(t)

    visit(expr)

    # 3. Apply all renames
    result = expr
    for (name, type_id), new_name in assigned.items():
        itype = L if type_id == id(L) else T
        result = replace_index(result,
                               Index(name, itype, True),
                               Index(new_name, itype, True))
        result = replace_index(result,
                               Index(name, itype, False),
                               Index(new_name, itype, False))
    return result


def _expr_is_zero(expr):
    """Check if a block expression is zero."""
    if isinstance(expr, TensorProduct) and expr.coeff == 0:
        return True
    if isinstance(expr, TensorSum) and len(expr.terms) == 0:
        return True
    return False


def _block_label(types, up_slots=(0,)):
    """Build an index label like R^{m}_{a n p} from child types."""
    L_gen = _name_gen('mnpq')
    T_gen = _name_gen('abcd')
    names = []
    for t in types:
        names.append(next(L_gen) if t is L else next(T_gen))
    ups = [names[i] for i in range(len(names)) if i in up_slots]
    downs = [names[i] for i in range(len(names)) if i not in up_slots]
    label = ''
    if ups:
        label += '^{' + ' '.join(ups) + '}'
    if downs:
        label += '_{' + ' '.join(downs) + '}'
    return label


# ═══════════════════════════════════════════════════════════════
# 1. Christoffel symbol blocks
# ═══════════════════════════════════════════════════════════════

print("=" * 64)
print("CHRISTOFFEL BLOCKS  Gamma^A_{BC}")
print("  from:  Gamma^A_{BC} = 1/2 g^{AD}(d_B g_{DC} + d_C g_{DB} - d_D g_{BC})")
print("=" * 64)

christoffel = bm.christoffel()

for (tA, tB, tC), expr in christoffel.items():
    if not _expr_is_zero(expr):
        label = 'Gamma' + _block_label([tA, tB, tC], up_slots=(0,))
        nice = relabel(expr)
        print(f"\n  {label}  =  {nice}")

print(f"\n  All other blocks  =  0")


# ═══════════════════════════════════════════════════════════════
# 2. Riemann tensor (all lower)
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 64)
print("RIEMANN BLOCKS  R_{ABCD}  (all lower, from defining formula)")
print("  R^A_{BCD} = d_C Gamma^A_{DB} - d_D Gamma^A_{CB}")
print("            + Gamma^A_{CE} Gamma^E_{DB} - Gamma^A_{DE} Gamma^E_{CB}")
print("  R_{ABCD} = g_{AE} R^E_{BCD}")
print("=" * 64)

riemann = bm.riemann(christoffel)
riemann_low = bm.riemann_lower(riemann)

for (tA, tB, tC, tD), expr in riemann_low.items():
    if not _expr_is_zero(expr):
        label = 'R' + _block_label([tA, tB, tC, tD], up_slots=())
        nice = relabel(expr)
        print(f"\n  {label}  =  {nice}")

print(f"\n  Odd T^2 index count  =>  0  (block-diagonal + truncation)")
print()
