"""
Riemann tensor of the 12d F-theory metric, computed from first principles.

The 12d metric is block-diagonal:

    G_MN = g_{mn}  oplus  M_{ab}

where g_{mn} is the 10d SL(2,R)-invariant Einstein frame metric
and M_{ab} = (1/tau_2) * [[1, tau_1], [tau_1, |tau|^2]] is the
axio-dilaton matrix on T^2.

We split 12d -> 10d + T^2, with truncation d_a = 0 on all fields.

RiemannGeometry computes everything from the defining formulas,
then decomposes into blocks via index splitting:
  1. Christoffel symbol blocks  Gamma^A_{BC}
  2. Riemann tensor blocks      R^A_{BCD}
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from tensorGlow import *
from tensorGlow.core.expr import Apply, Tensor, Prod, Sum, Scalar


# ═══════════════════════════════════════════════════════════════
# Index types
# ═══════════════════════════════════════════════════════════════

L = IndexType('10d', dim=10, dummy_prefix='L')
T = IndexType('T2', dim=2, dummy_prefix='T')

parent = IndexType('12d', dim=12, dummy_prefix='P')
split = IndexSplit(parent, [L, T])


# ═══════════════════════════════════════════════════════════════
# 12d metric with block-diagonal structure and truncation
# ═══════════════════════════════════════════════════════════════

G = MetricTensor(parent, name='G')

geom = RiemannGeometry(G, parent)
geom.split_indices(split)
# Block-diagonal: cross-blocks vanish
geom.set_metric_zero(L, T)
geom.set_metric_zero(T, L)
# Truncation: fields are T²-independent
geom.set_partial_zero(T)


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


def _rename_idx(idx, rename_map):
    """Rename a single index if it's in the map."""
    key = (idx.name, idx.index_type)
    if key in rename_map:
        return Index(rename_map[key], idx.index_type, idx.is_up)
    return idx


def _bulk_rename(expr, rename_map):
    """Rename all indices in one pass using rename_map: (name, IndexType) → new_name."""
    if isinstance(expr, Scalar):
        return expr

    if isinstance(expr, Tensor):
        new_indices = tuple(_rename_idx(idx, rename_map) for idx in expr.indices)
        return Tensor(expr.head, new_indices)

    if isinstance(expr, Apply):
        new_di = _rename_idx(expr.deriv_index, rename_map)
        new_operand = _bulk_rename(expr.operand, rename_map)
        return Apply(expr.op, new_di, new_operand)

    if isinstance(expr, Prod):
        new_factors = tuple(_bulk_rename(f, rename_map) for f in expr.factors)
        return Prod(expr.coeff, new_factors)

    if isinstance(expr, Sum):
        new_terms = tuple(_bulk_rename(t, rename_map) for t in expr.terms)
        return Sum(new_terms)

    return expr


def relabel(expr):
    """Replace all canonical index names with readable physics indices.

    10d: m, n, p, q, ...
    T²: a, b, c, d, ...
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
        elif isinstance(e, Apply):
            assign(e.deriv_index.name, e.deriv_index.index_type)
            visit(e.operand)
        elif isinstance(e, TensorProduct):
            for f in e.factors:
                visit(f)
        elif isinstance(e, TensorSum):
            for t in e.terms:
                visit(t)

    visit(expr)

    # 3. Build rename map: (old_name, IndexType) → new_name
    rename_map = {}
    for (name, type_id), new_name in assigned.items():
        itype = L if type_id == id(L) else T
        rename_map[(name, itype)] = new_name

    # 4. Apply all renames in a single pass (avoids sequential collision)
    return _bulk_rename(expr, rename_map)


def _expr_is_zero(expr):
    """Check if a block expression is zero."""
    if isinstance(expr, TensorProduct) and expr.coeff == 0:
        return True
    if isinstance(expr, TensorSum) and len(expr.terms) == 0:
        return True
    return False


# ── Symmetry filters: pick one canonical rep per equivalence class ──

_type_rank = {id(L): 0, id(T): 1}

def _tkey(types):
    """Sortable tuple for a block key."""
    return tuple(_type_rank[id(t)] for t in types)


def _christoffel_is_canonical(key):
    """Gamma^A_{BC} is symmetric in (B, C).  Keep B <= C."""
    _, tB, tC = key
    return _type_rank[id(tB)] <= _type_rank[id(tC)]


def _riemann_is_canonical(key):
    """R_{ABCD} has antisym in (A,B), antisym in (C,D), pair sym (AB)↔(CD).

    Two blocks are "the same component" if related by any of the 8
    Riemann symmetry permutations (regardless of sign).  Pick the
    lex-smallest block-key in the orbit as the canonical representative.
    """
    t1, t2, t3, t4 = key
    variants = [
        (t1, t2, t3, t4),      # identity
        (t2, t1, t3, t4),      # swap first pair
        (t1, t2, t4, t3),      # swap second pair
        (t2, t1, t4, t3),      # swap both pairs
        (t3, t4, t1, t2),      # pair exchange
        (t4, t3, t1, t2),      # pair exchange + swap first
        (t3, t4, t2, t1),      # pair exchange + swap second
        (t4, t3, t2, t1),      # pair exchange + swap both
    ]
    canon = min(variants, key=_tkey)
    return key == canon


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
# Free indices for the formulas (parent type)
# ═══════════════════════════════════════════════════════════════

A, B, C, D = (Index(n, parent, is_up=False) for n in ['A', 'B', 'C', 'D'])
A_up = Index('A', parent, is_up=True)


# ═══════════════════════════════════════════════════════════════
# 1. Christoffel symbol blocks
# ═══════════════════════════════════════════════════════════════

print("=" * 64)
print("CHRISTOFFEL BLOCKS  Gamma^A_{BC}")
print("  from:  Gamma^A_{BC} = 1/2 G^{AD}(d_B G_{DC} + d_C G_{DB} - d_D G_{BC})")
print("=" * 64)

christoffel_expr = geom.christoffel_formula(A_up, -B, -C)
christoffel_blocks = geom.decompose(christoffel_expr)

for block_key, expr in sorted(christoffel_blocks.items(), key=lambda x: str(x[0])):
    if not _expr_is_zero(expr) and _christoffel_is_canonical(block_key):
        label = 'Gamma' + _block_label(list(block_key), up_slots=(0,))
        nice = relabel(expr)
        print(f"\n  {label}  =  {nice}")

print(f"\n  All other blocks  =  0  (or related by Gamma^A_{{BC}} = Gamma^A_{{CB}})")


# ═══════════════════════════════════════════════════════════════
# 2. Riemann tensor blocks
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 64)
print("RIEMANN BLOCKS  R_{ABCD}")
print("  R_{ABCD} = G_{AE} R^E_{BCD}")
print("           = G_{AE} (d_C Gamma^E_{DB} - d_D Gamma^E_{CB}")
print("                    + Gamma^E_{CF} Gamma^F_{DB} - Gamma^E_{DF} Gamma^F_{CB})")
print("=" * 64)

riemann_expr = geom.riemann_formula(-A, -B, -C, -D)
riemann_blocks = geom.decompose(riemann_expr)

for block_key, expr in sorted(riemann_blocks.items(), key=lambda x: str(x[0])):
    if not _expr_is_zero(expr) and _riemann_is_canonical(block_key):
        label = 'R' + _block_label(list(block_key), up_slots=())
        nice = relabel(expr)
        print(f"\n  {label}  =  {nice}")

print(f"\n  Others  =  0, or related by R_{{ABCD}} = -R_{{BACD}} = -R_{{ABDC}} = R_{{CDAB}}")
print()
