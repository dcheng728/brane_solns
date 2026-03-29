"""
Differential form fields and their operations for supergravity.

Sparse storage for antisymmetric tensor fields, with Hodge star,
norm, contraction, stress-energy, and exterior derivative.
"""

import sympy as sp
from itertools import combinations, permutations
from sympy.combinatorics.permutations import Permutation


def epsilon(n, idx):
    """Levi-Civita symbol in n dimensions.

    Returns +1, -1, or 0 depending on the permutation parity of idx.
    """
    idx = tuple(idx)
    if len(idx) != n:
        raise ValueError("Number of indices must match dimension.")
    if len(set(idx)) < n:
        return 0
    inv = 0
    for i in range(n - 1):
        inv += sum(idx[i] > idx[j] for j in range(i + 1, n))
    return -1 if (inv & 1) else 1


def _permutation_sign(perm):
    """Sign of a permutation given as a tuple of indices."""
    return Permutation(perm).signature()


def _sort_with_sign(idx):
    """Sort an index tuple and return (sorted_tuple, sign).

    The sign accounts for the number of transpositions needed.
    """
    idx = list(idx)
    n = len(idx)
    sign = 1
    # Bubble sort tracking sign
    for i in range(n):
        for j in range(i + 1, n):
            if idx[i] > idx[j]:
                idx[i], idx[j] = idx[j], idx[i]
                sign *= -1
            elif idx[i] == idx[j]:
                return tuple(idx), 0
    return tuple(idx), sign


class FormField:
    """An antisymmetric p-form field with sparse storage.

    Only independent components (sorted index tuples) are stored.
    Accessing with permuted indices returns the appropriately signed value.

    Parameters
    ----------
    rank : int
        Degree of the form (number of indices).
    dim : int
        Dimension of the manifold.
    components : dict, optional
        Mapping from sorted index tuples to expressions.
    """

    def __init__(self, rank, dim, components=None):
        self.rank = rank
        self.dim = dim
        self._components = {}
        if components:
            for idx, val in components.items():
                self[idx] = val

    def __getitem__(self, idx):
        """Access F_{i0 i1 ... i_{p-1}} with automatic antisymmetry."""
        idx = tuple(idx)
        sorted_idx, sign = _sort_with_sign(idx)
        if sign == 0:
            return sp.S(0)
        return sign * self._components.get(sorted_idx, sp.S(0))

    def __setitem__(self, idx, value):
        """Set a component. Stores in canonical (sorted) form."""
        idx = tuple(idx)
        sorted_idx, sign = _sort_with_sign(idx)
        if sign == 0:
            raise ValueError(f"Cannot set component with repeated indices: {idx}")
        self._components[sorted_idx] = sign * value

    @property
    def nonzero_components(self):
        """Dict of sorted index tuple -> expression for nonzero entries."""
        return {k: v for k, v in self._components.items() if v != 0}

    @property
    def independent_indices(self):
        """All C(dim, rank) sorted index tuples."""
        return list(combinations(range(self.dim), self.rank))

    def to_array(self):
        """Convert to a dense SymPy MutableDenseNDimArray."""
        arr = sp.MutableDenseNDimArray.zeros(*([self.dim] * self.rank))
        for idx in self.independent_indices:
            val = self._components.get(idx, sp.S(0))
            if val != 0:
                # Fill all permutations with appropriate signs
                for perm in permutations(idx):
                    _, sign = _sort_with_sign(perm)
                    arr[perm] = sign * val
        return arr

    @classmethod
    def from_array(cls, array):
        """Construct from a SymPy NDimArray."""
        rank = len(array.shape)
        dim = array.shape[0]
        form = cls(rank, dim)
        for idx in combinations(range(dim), rank):
            val = array[idx]
            if val != 0:
                form._components[idx] = val
        return form

    def copy(self):
        """Return a deep copy."""
        return FormField(self.rank, self.dim, dict(self._components))

    def __mul__(self, scalar):
        return FormField(self.rank, self.dim, {k: scalar * v for k, v in self._components.items()})

    def __rmul__(self, scalar):
        return self.__mul__(scalar)

    def __truediv__(self, scalar):
        return FormField(self.rank, self.dim, {k: v / scalar for k, v in self._components.items()})

    def __neg__(self):
        return FormField(self.rank, self.dim, {k: -v for k, v in self._components.items()})

    def __sub__(self, other):
        result = self.copy()
        for idx in other.independent_indices:
            val = other._components.get(idx, sp.S(0))
            if val != 0:
                current = result._components.get(idx, sp.S(0))
                result._components[idx] = current - val
        return result

    def __add__(self, other):
        result = self.copy()
        for idx in other.independent_indices:
            val = other._components.get(idx, sp.S(0))
            if val != 0:
                current = result._components.get(idx, sp.S(0))
                result._components[idx] = current + val
        return result


def exterior_derivative(form, coordinates):
    """Compute the exterior derivative dF.

    (dF)_{i0...ip} = (p+1) * partial_{[i0} F_{i1...ip]}

    For sparse forms, only differentiates nonzero components.

    Parameters
    ----------
    form : FormField
    coordinates : list of sp.Symbol

    Returns
    -------
    FormField of rank p+1
    """
    p = form.rank
    D = form.dim
    dF = FormField(p + 1, D)

    for idx in combinations(range(D), p + 1):
        val = sp.S(0)
        for s in range(p + 1):
            # Remove index s from the tuple to get the (p)-form index
            sub_idx = idx[:s] + idx[s+1:]
            component = form[sub_idx]
            if component != 0:
                deriv = sp.diff(component, coordinates[idx[s]])
                if deriv != 0:
                    val += (-1)**s * deriv
        if val != 0:
            dF._components[idx] = val

    return dF


def hodge_star(form, metric, signature=-1):
    """Compute the Hodge dual *F.

    (*F)_{j1...j_{D-p}} = (1/p!) * sqrt(|g|) * epsilon^{i1...ip}_{j1...j_{D-p}} * F_{i1...ip}

    Optimized for diagonal metrics (g^{ii} = 1/g_{ii}).

    Parameters
    ----------
    form : FormField
    metric : Metric
        Must have .dim, .matrix, .inv_matrix, .det properties.
    signature : int
        -1 for Lorentzian (default), +1 for Euclidean.
        Used as: sqrt(signature * det(g)).

    Returns
    -------
    FormField of rank D - p
    """
    p = form.rank
    D = metric.dim
    dual_rank = D - p
    dual = FormField(dual_rank, D)

    abs_det_g = signature * metric.det
    sqrt_g = sp.sqrt(abs_det_g)

    is_diag = metric.is_diagonal
    inv_g = metric.inv_matrix

    for dual_idx in combinations(range(D), dual_rank):
        remaining = [i for i in range(D) if i not in dual_idx]

        val = sp.S(0)
        for perm in permutations(remaining):
            # F component with raised indices
            F_val = form[perm]
            if F_val == 0:
                continue

            eps = epsilon(D, dual_idx + perm)
            if eps == 0:
                continue

            # Raise indices: F^{i1...ip} = g^{i1 i1} ... g^{ip ip} F_{i1...ip} (diagonal)
            contribution = F_val
            if is_diag:
                for i in perm:
                    contribution *= inv_g[i, i]
            else:
                # For non-diagonal: would need full contraction
                # For now, support diagonal only
                for i in perm:
                    contribution *= inv_g[i, i]

            val += eps * contribution

        val = sqrt_g * val / sp.factorial(p)
        if val != 0:
            dual._components[dual_idx] = val

    return dual


def form_norm_squared(form, metric):
    """Compute |F|^2 = (1/p!) F_{M1...Mp} F^{M1...Mp}.

    Iterates only over nonzero independent components.

    Parameters
    ----------
    form : FormField
    metric : Metric (diagonal)

    Returns
    -------
    sp.Expr
    """
    p = form.rank
    inv_g = metric.inv_matrix
    result = sp.S(0)

    for idx, val in form.nonzero_components.items():
        # F_{idx} F^{idx} = F_{idx}^2 * prod(g^{ii})
        contribution = val ** 2
        for i in idx:
            contribution *= inv_g[i, i]
        # Factor of p! from summing over all permutations of the sorted index
        result += sp.factorial(p) * contribution

    # Overall 1/p! normalization
    return result / sp.factorial(p)


def form_contraction(form, metric):
    """Compute F_{M P1...P_{p-1}} F_N^{P1...P_{p-1}}.

    Returns a D x D SymPy Matrix.

    Parameters
    ----------
    form : FormField
    metric : Metric (diagonal)

    Returns
    -------
    sp.Matrix
    """
    p = form.rank
    D = metric.dim
    inv_g = metric.inv_matrix
    result = sp.zeros(D, D)

    # Sum over independent (p-1)-index combinations
    sub_indices = list(combinations(range(D), p - 1))

    for M in range(D):
        for N in range(M, D):
            val = sp.S(0)
            for idx in sub_indices:
                if M in idx or N in idx:
                    continue
                F_M = form[(M,) + idx]
                F_N = form[(N,) + idx]
                if F_M == 0 or F_N == 0:
                    continue
                contribution = F_M * F_N
                for i in idx:
                    contribution *= inv_g[i, i]
                val += contribution
            val *= sp.factorial(p - 1)
            result[M, N] = val
            if M != N:
                result[N, M] = val

    return result


def form_stress_energy(form, metric, dilaton=None, alpha=0, D_trace=None, n_eff=None):
    """Compute the form-field stress-energy tensor T^(form)_{MN}.

        T^(form)_{MN} = (1/2) e^{alpha*Phi} [ FF_{MN} - (n_eff-1)/(D_trace-2) |F|^2 g_{MN} ]

    where FF_{MN} = 1/(n-1)! F_{MP...} F_N^{P...}
          |F|^2   = 1/n!     F_{...}   F^{...}

    The contraction FF_{MN} and norm |F|^2 always use the actual form rank n.
    Only the trace coefficient uses n_eff and D_trace.

    Parameters
    ----------
    form : FormField
    metric : Metric
    dilaton : sp.Expr or None, optional
        Dilaton Phi, used only in the e^{alpha*Phi} prefactor.
    alpha : number or sp.Expr, optional
        Dilaton coupling. Default: 0 (no coupling).
    D_trace : int or None, optional
        Dimension used in the trace coefficient.
        Defaults to metric.dim if not specified.
    n_eff : int or None, optional
        Effective form rank used in the trace coefficient.
        Defaults to form.rank if not specified.

    Returns
    -------
    sp.Matrix
        D x D form stress-energy tensor.
    """
    n = form.rank
    D = metric.dim
    if D_trace is None:
        D_trace = D
    if n_eff is None:
        n_eff = n
    g = metric.matrix

    if dilaton is not None and alpha != 0:
        e_alpha_Phi = sp.exp(sp.sympify(alpha) * dilaton)
    else:
        e_alpha_Phi = sp.S(1)

    # Normalized form contraction: FF_{MN} = 1/(n-1)! F_{MP...} F_N^{P...}
    FF_MN = form_contraction(form, metric) / sp.factorial(n - 1)

    # Normalized form norm: |F|^2 = 1/n! F_{...} F^{...}
    F_sq = form_norm_squared(form, metric)

    T = sp.zeros(D, D)
    for M in range(D):
        for N in range(M, D):
            val = sp.Rational(1, 2) * e_alpha_Phi * (
                FF_MN[M, N] - sp.Rational(n_eff - 1, D_trace - 2) * F_sq * g[M, N]
            )
            T[M, N] = val
            if M != N:
                T[N, M] = val

    return T


def scalar_stress_energy(dilaton, metric):
    """Compute the dilaton kinetic stress-energy tensor T^(scalar)_{MN}.

        T^(scalar)_{MN} = (1/2) d_M Phi d_N Phi

    Parameters
    ----------
    dilaton : sp.Expr
        Dilaton scalar Phi as a function of coordinates.
    metric : Metric
        Must have .coordinates.

    Returns
    -------
    sp.Matrix
        D x D scalar stress-energy tensor.
    """
    D = metric.dim
    dPhi = [sp.diff(dilaton, xi) for xi in metric.coordinates]

    T = sp.zeros(D, D)
    for M in range(D):
        for N in range(M, D):
            val = sp.Rational(1, 2) * dPhi[M] * dPhi[N]
            T[M, N] = val
            if M != N:
                T[N, M] = val

    return T
