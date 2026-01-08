from itertools import product, permutations
import sympy as sp
from sympy.combinatorics.permutations import Permutation

def epsilon(n, idx):
    """
    Returns the Levi-Civita symbol in the given number of dimensions with the specified indices.
    The indices should be a tuple or list of integers.

    Examples:
    >>> epsilon(3, [0, 1, 2]) # returns 1
    >>> epsilon(3, [1, 0, 2]) # returns -1
    >>> epsilon(3, [0, 2, 2]) # returns 0
    """
    if len(idx) != n:
        raise ValueError("Number of indices must match the number of dimensions.")
    if any(i < 0 or i >= n for i in idx):
        raise ValueError("Indices must be within the valid range.")
    if len(set(idx)) < n:
        return 0  # repeated index → 0

    # parity via inversion count
    inv = 0
    a = tuple(idx)
    for i in range(n - 1):
        ai = a[i]
        inv += sum(ai > a[j] for j in range(i + 1, n))
    return -1 if (inv & 1) else 1


def perm_with_sign(t):
    """
    Generate all permutations of the input tuple `t` along with their sign.

    Usage:
    t = (1, 2, 3, 11)
    for tup, sgn in perm_with_sign(t):
        print(tup, sgn)
    """
    n = len(t)
    for p in permutations(range(n)):                # p is a perm of indices
        yield tuple(t[i] for i in p), Permutation(p).signature()

def find_nonzero_indices(_tensor):
    """
    Given a sympy tensor F, find all indices where F is nonzero.
    Returns a list of index tuples.
    """
    nonzeros = []
    for idx in product(*(range(s) for s in _tensor.shape)):
        val = _tensor[idx]
        if val != 0:                    # fast, but purely syntactic
            nonzeros.append(idx)
    return nonzeros

def antisymmetrize_tensor(__tensor):
    """
    Given a sympy tensor F, antisymmetrize it by filling in all permutations
    of its nonzero components with the appropriate sign.
    """
    nonzeros = find_nonzero_indices(__tensor)

    # iterate over the nonzero indices, apply antisymmetry
    for nonzero_index in nonzeros:
        for permuted_indices, sgn in perm_with_sign(nonzero_index):
            __tensor[permuted_indices] = sgn * __tensor[nonzero_index]
    return __tensor

def get_independent_form_indices(num_dimensions, rank):
    """
    Get a list of independent index combinations for an antisymmetric form
    of given rank in the specified number of dimensions.

    Example:
    >>> get_independent_form_indices(4, 2)
    [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    """
    from itertools import combinations
    return list(combinations(range(num_dimensions), rank))

def hodge_dual(tensor, metric, signature=-1):
    """
    we assume the metric is diagonal and that tensor is a form field
    """
    abs_det_g = signature * sp.det(metric)
    sqrt_abs_det_g = sp.sqrt(abs_det_g)
    inv_metric = sp.simplify(metric.inv())

    rank = len(tensor.shape)
    num_dims = metric.shape[0]
    dual_rank = num_dims - rank
    dual_tensor = sp.MutableDenseNDimArray.zeros(*([num_dims] * dual_rank))

    # find independent components of the dual tensor, upto antisymmetry
    independent_dual_indices = get_independent_form_indices(num_dims, dual_rank)

    for dual_idx in independent_dual_indices:
        remaining_indices = [i for i in range(num_dims) if i not in dual_idx]
        # print(dual_idx,remaining_indices)

        # sum over all permutations of the remaining indices
        running_sum = sp.S(0)
        remaining_indices_perms = permutations(remaining_indices)
        for remaining_perm in remaining_indices_perms:
            # print(dual_idx + remaining_perm)

            contribution = sqrt_abs_det_g * epsilon(num_dims, dual_idx + remaining_perm) * tensor[remaining_perm]
            for summed_index in remaining_perm:
                # raise index using inverse metric
                contribution *= inv_metric[summed_index, summed_index]
            running_sum += contribution

        dual_tensor[dual_idx] = running_sum / sp.factorial(rank)

    dual_tensor = antisymmetrize_tensor(dual_tensor)

    return dual_tensor

def compute_form_squared(form_field, inv_metric):
    """
    Compute the squared norm of a form field using the inverse metric.
    Assumes the form field is fully antisymmetric, and the metric is diagonal
    """
    num_dims = inv_metric.shape[0]
    rank = len(form_field.shape)
    
    print(f"Computing squared norm of a rank {rank} form in {num_dims} dimensions.")

    squared_norm = sp.S(0)

    # # Get independent components of the tensor
    independent_indices = get_independent_form_indices(num_dims, rank)

    for idx in independent_indices:
        component = form_field[idx]
        if component != 0:
            contribution = component**2
            for i in idx:
                contribution *= inv_metric[i, i]
            squared_norm += contribution * sp.factorial(rank)

    return squared_norm

def compute_FF_MN(form_field, inv_metric):
    """
    Compute the tensor F_{MPQRS} F_N^{PQRS}
    Assumes the form field is fully antisymmetric, and the metric is diagonal
    """
    num_dims = inv_metric.shape[0]
    rank = len(form_field.shape)
    
    print(f"Computing F_M PQRS F_N^PQRS for a rank {rank} form in {num_dims} dimensions.")

    toret = sp.MutableDenseNDimArray.zeros(num_dims, num_dims)
    
    # we are essentially summing over a worth of rank 4 antismmetric indices
    independent_sub_indices = get_independent_form_indices(num_dims, rank-1)

    for M in range(num_dims):
        for N in range(num_dims):
            sum_term = sp.S(0)
            for idx in independent_sub_indices:
                if M in idx or N in idx: # vanishes due to antisymmetry
                    continue
                component_M = form_field[(M,) + idx]
                component_N = form_field[(N,) + idx]
                contribution = component_M * component_N
                for i in idx:
                    contribution *= inv_metric[i, i]
                sum_term += contribution
            toret[M, N] = sp.simplify(sum_term)

    return toret * sp.factorial(rank - 1)

def compute_S_MN(metric, rank, ff_MN, ff):
    """
    Given FF_{MN} and F^2, compute S_{MN} as defined in eq 2.2b of hep-th/9701088
    """
    n_dims = metric.shape[0]
    toret = sp.MutableDenseNDimArray.zeros(n_dims, n_dims)

    for M in range(n_dims):
        for N in range(n_dims):
            toret[M, N] = sp.Rational(1,2 * sp.factorial(rank - 1) ) * (ff_MN[M, N] - sp.Rational(rank-1,rank * (n_dims-2)) * metric[M, N] * ff)
    return toret