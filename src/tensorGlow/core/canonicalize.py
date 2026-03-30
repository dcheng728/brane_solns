"""Canonicalization engine: bridge to SymPy's Butler-Portugal algorithm.

Slot convention
---------------
tensor_can.canonicalize works with permutations of size ``n_slots + 2``
where ``n_slots`` is the total number of index slots across all tensor
factors.  The last two entries encode the sign:

    perm[n_slots] == n_slots      → positive sign
    perm[n_slots] == n_slots + 1  → negative sign

The permutation ``g`` maps slot position → canonical index position.
Free indices occupy the first positions in canonical order, followed
by dummy pairs listed as [d0_up, d0_down, d1_up, d1_down, ...].

Each tensor type provides a BSGS (base, strong generating set) whose
generators act on ``rank_i + 2`` elements (its own slots plus 2 sign bits).
``canonicalize`` is told ``(base, gens, n_tensors_of_this_type, sym)``
for each tensor type, where ``n_tensors_of_this_type`` groups identical
tensors.  For our use, each atom is its own type with n=1.
"""

from sympy.combinatorics import Permutation
from sympy.combinatorics.tensor_can import canonicalize as _canon_bp

from .index import Index
from .expr import TensorAtom, TensorProduct, TensorSum, ScalarExpr


def canonicalize_expr(expr):
    """Canonicalize any TensorExpr."""
    if isinstance(expr, TensorProduct):
        return canonicalize_product(expr)
    if isinstance(expr, TensorSum):
        return canonicalize_sum(expr)
    if isinstance(expr, ScalarExpr):
        return expr
    raise TypeError(f"Cannot canonicalize {type(expr)}")


def canonicalize_sum(tensor_sum):
    """Canonicalize each term, then collect like terms."""
    new_terms = []
    for term in tensor_sum.terms:
        canon = canonicalize_product(term)
        if isinstance(canon, TensorProduct) and canon.coeff != 0:
            new_terms.append(canon)
        elif isinstance(canon, TensorSum):
            new_terms.extend(canon.terms)
    result = TensorSum(tuple(new_terms))
    return result.collect()


def canonicalize_product(product):
    """Canonicalize a TensorProduct using Butler-Portugal.

    Returns a TensorProduct in canonical form (or with coeff=0 if zero).
    """
    if not product.atoms or product.coeff == 0:
        return product

    # ── Collect all slot indices ─────────────────────────────────────
    all_indices = []  # flat list of Index, one per slot
    for atom in product.atoms:
        all_indices.extend(atom.indices)

    n_slots = len(all_indices)
    if n_slots == 0:
        return product
    size = n_slots + 2

    # ── Classify indices as free or dummy ────────────────────────────
    # Build a map: (name, index_type_id) → list of (slot_pos, Index)
    idx_groups = {}
    for slot, idx in enumerate(all_indices):
        key = (idx.name, id(idx.index_type))
        idx_groups.setdefault(key, []).append((slot, idx))

    free_slots = []     # slot positions of free indices
    dummy_pairs = []    # list of (slot_up, slot_down)
    dummy_msym = []     # metric symmetry for each dummy pair

    for key, entries in idx_groups.items():
        if len(entries) == 1:
            free_slots.append(entries[0][0])
        elif len(entries) == 2:
            (s0, i0), (s1, i1) = entries
            if i0.is_up and not i1.is_up:
                dummy_pairs.append((s0, s1))
            elif i1.is_up and not i0.is_up:
                dummy_pairs.append((s1, s0))
            else:
                # Same variance — treat as free (shouldn't happen for valid exprs)
                free_slots.extend([s0, s1])
                continue
            dummy_msym.append(i0.index_type.metric_symmetry)
        else:
            raise ValueError(
                f"Index {entries[0][1].name} appears {len(entries)} times "
                f"(expected 1 or 2)"
            )

    # Sort free slots by their current index name for canonical ordering
    free_slots.sort(key=lambda s: all_indices[s].name)

    # ── Build the permutation g ──────────────────────────────────────
    # g maps: slot position → canonical position
    # Canonical positions: free indices first (sorted), then dummy pairs
    g = [0] * size

    # Free index positions
    for canon_pos, slot in enumerate(free_slots):
        g[slot] = canon_pos

    # Dummy pair positions (after free indices)
    n_free = len(free_slots)
    dummies_flat = []
    for pair_idx, (slot_up, slot_down) in enumerate(dummy_pairs):
        canon_up = n_free + 2 * pair_idx
        canon_down = n_free + 2 * pair_idx + 1
        g[slot_up] = canon_up
        g[slot_down] = canon_down
        dummies_flat.extend([canon_up, canon_down])

    # Sign bits (positive)
    g[n_slots] = n_slots  # not a real slot, but needed for size
    # Actually we build g as size n_slots + 2
    g_arr = [0] * size
    for slot in range(n_slots):
        g_arr[slot] = g[slot]
    g_arr[n_slots] = n_slots
    g_arr[n_slots + 1] = n_slots + 1
    g_perm = Permutation(g_arr)

    # ── Build v-list (one entry per atom) ────────────────────────────
    # Each atom is its own tensor type with n=1.
    # The BSGS generators are on (rank + 2) elements.
    # canonicalize handles embedding into the full size.
    v_args = []
    for atom in product.atoms:
        sym = atom.head.symmetry
        v_args.append((sym.base, sym.generators, 1, 0))

    # ── Call tensor_can ──────────────────────────────────────────────
    try:
        canon = _canon_bp(g_perm, dummies_flat, dummy_msym, *v_args)
    except Exception:
        return product  # fallback: return unchanged

    if canon == 0:
        return TensorProduct(0, ())

    # ── Interpret the result ─────────────────────────────────────────
    canon_list = list(canon)

    # Extract sign
    sign = 1 if canon_list[n_slots] == n_slots else -1

    # Build the inverse map: canonical position → slot position
    # canon_list[slot] = canonical_position, so we need the inverse
    inv_canon = [0] * n_slots
    for slot in range(n_slots):
        inv_canon[canon_list[slot]] = slot

    # Reconstruct: for each canonical position, which index goes there?
    # The canonical ordering has free indices first, then dummies.
    # We need to rebuild the index list for each atom.

    # Build the canonical index at each canonical position:
    # Free positions: the sorted free index names
    canon_indices = [None] * n_slots
    for canon_pos, slot in enumerate(free_slots):
        canon_indices[canon_pos] = all_indices[slot]

    for pair_idx, (slot_up, slot_down) in enumerate(dummy_pairs):
        canon_up = n_free + 2 * pair_idx
        canon_down = n_free + 2 * pair_idx + 1
        canon_indices[canon_up] = all_indices[slot_up]
        canon_indices[canon_down] = all_indices[slot_down]

    # Now: the canonical permutation tells us that slot s should hold
    # the index that was at canonical position canon_list[s].
    # So new_index_at_slot[s] = canon_indices[canon_list[s]]
    new_all_indices = [None] * n_slots
    for s in range(n_slots):
        new_all_indices[s] = canon_indices[canon_list[s]]

    # Rebuild atoms
    new_atoms = []
    offset = 0
    for atom in product.atoms:
        rank = atom.head.rank
        new_idx = tuple(new_all_indices[offset:offset + rank])
        new_atoms.append(TensorAtom(atom.head, new_idx))
        offset += rank

    return TensorProduct(product.coeff * sign, tuple(new_atoms))


def rename_dummies(product, existing_names=None):
    """Rename dummy indices to canonical form (prefix_0, prefix_1, ...).

    Parameters
    ----------
    product : TensorProduct
    existing_names : set, optional
        Names to avoid.
    """
    if not isinstance(product, TensorProduct) or not product.atoms:
        return product

    if existing_names is None:
        existing_names = set()

    pairs = product.dummy_pairs
    if not pairs:
        return product

    rename = {}
    counter = {}

    for (ai, si), (aj, sj) in pairs:
        idx = product.atoms[ai].indices[si]
        prefix = idx.index_type.dummy_prefix
        if prefix not in counter:
            counter[prefix] = 0
        new_name = f"{prefix}_{counter[prefix]}"
        while new_name in existing_names:
            counter[prefix] += 1
            new_name = f"{prefix}_{counter[prefix]}"
        counter[prefix] += 1
        rename[idx.name] = new_name

    if not rename:
        return product

    new_atoms = []
    for atom in product.atoms:
        new_indices = tuple(
            Index(rename.get(idx.name, idx.name), idx.index_type, idx.is_up)
            for idx in atom.indices
        )
        new_atoms.append(TensorAtom(atom.head, new_indices))

    return TensorProduct(product.coeff, tuple(new_atoms))
