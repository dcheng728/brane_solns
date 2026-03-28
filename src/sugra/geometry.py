"""
Metric and curvature computations for supergravity solutions.

Provides diagonal-optimized Christoffel symbol and Ricci tensor computation,
harmonic function utilities, and warped product metric construction.
"""

import sympy as sp
from functools import cached_property


class Metric:
    """A pseudo-Riemannian metric with diagonal-optimized curvature computation.

    For diagonal metrics, exploits closed-form Christoffel formulas:
      Gamma^i_ii = (d_i g_ii) / (2 g_ii)
      Gamma^i_ij = (d_j g_ii) / (2 g_ii)       (i != j)
      Gamma^i_jj = -(d_i g_jj) / (2 g_ii)      (i != j)
    giving O(D^2) instead of O(D^4) work.

    For non-diagonal metrics, falls back to the full computation.

    Parameters
    ----------
    matrix : sp.Matrix
        D x D symmetric metric matrix.
    coordinates : list of sp.Symbol
        D coordinate symbols.
    """

    def __init__(self, matrix, coordinates):
        if isinstance(matrix, list):
            # Interpret as diagonal entries
            matrix = sp.diag(*matrix)
        self._matrix = sp.Matrix(matrix)
        self._coordinates = list(coordinates)
        self._dim = self._matrix.shape[0]
        assert self._matrix.shape == (self._dim, self._dim)
        assert len(self._coordinates) == self._dim
        # Detect diagonal structure
        self._is_diagonal = all(
            self._matrix[i, j] == 0
            for i in range(self._dim)
            for j in range(self._dim)
            if i != j
        )

    @property
    def dim(self):
        return self._dim

    @property
    def is_diagonal(self):
        return self._is_diagonal

    @property
    def matrix(self):
        return self._matrix

    @property
    def coordinates(self):
        return self._coordinates

    @cached_property
    def inv_matrix(self):
        if self._is_diagonal:
            return sp.diag(*[1 / self._matrix[i, i] for i in range(self._dim)])
        return self._matrix.inv()

    @cached_property
    def det(self):
        if self._is_diagonal:
            result = sp.S(1)
            for i in range(self._dim):
                result *= self._matrix[i, i]
            return result
        return self._matrix.det()

    @cached_property
    def sqrt_abs_det(self):
        return sp.sqrt(sp.Abs(self.det))

    def christoffel(self, simplify_func=None):
        """Compute nonzero Christoffel symbols Gamma^i_{jk}.

        Returns dict mapping (i, j, k) -> expression.
        Symmetric in the lower two indices: only (i, j, k) with j <= k is stored.

        Parameters
        ----------
        simplify_func : callable, optional
            Applied to each component (e.g., sp.cancel, sp.factor). Default: None.
        """
        if self._is_diagonal:
            return self._christoffel_diagonal(simplify_func)
        return self._christoffel_general(simplify_func)

    def _christoffel_diagonal(self, simplify_func):
        D = self._dim
        g = self._matrix
        x = self._coordinates
        Gamma = {}

        for i in range(D):
            g_ii = g[i, i]
            inv_g_ii = 1 / g_ii

            # Gamma^i_ii = (d_i g_ii) / (2 g_ii)
            dg = sp.diff(g_ii, x[i])
            if dg != 0:
                val = sp.Rational(1, 2) * inv_g_ii * dg
                if simplify_func:
                    val = simplify_func(val)
                Gamma[(i, i, i)] = val

            for j in range(D):
                if j == i:
                    continue

                # Gamma^i_ij = (d_j g_ii) / (2 g_ii)
                dg_ij = sp.diff(g_ii, x[j])
                if dg_ij != 0:
                    val = sp.Rational(1, 2) * inv_g_ii * dg_ij
                    if simplify_func:
                        val = simplify_func(val)
                    key = (i, i, j) if i <= j else (i, j, i)
                    Gamma[key] = val

                # Gamma^i_jj = -(d_i g_jj) / (2 g_ii)
                g_jj = g[j, j]
                dg_ji = sp.diff(g_jj, x[i])
                if dg_ji != 0:
                    val = -sp.Rational(1, 2) * inv_g_ii * dg_ji
                    if simplify_func:
                        val = simplify_func(val)
                    Gamma[(i, j, j)] = val

        return Gamma

    def _christoffel_general(self, simplify_func):
        D = self._dim
        g = self._matrix
        ginv = self.inv_matrix
        x = self._coordinates
        Gamma = {}

        for i in range(D):
            for j in range(D):
                for k in range(j, D):  # symmetry: j <= k
                    val = sp.S(0)
                    for l in range(D):
                        val += sp.Rational(1, 2) * ginv[i, l] * (
                            sp.diff(g[l, j], x[k])
                            + sp.diff(g[l, k], x[j])
                            - sp.diff(g[j, k], x[l])
                        )
                    if simplify_func:
                        val = simplify_func(val)
                    if val != 0:
                        Gamma[(i, j, k)] = val

        return Gamma

    def _get_christoffel(self, Gamma, i, j, k):
        """Look up Gamma^i_{jk} from the stored dict, using symmetry in (j,k)."""
        key = (i, j, k) if j <= k else (i, k, j)
        return Gamma.get(key, sp.S(0))

    def _detect_composite_functions(self):
        """Detect composite function applications in metric components.

        If metric components involve f(ξ(x1,...,xn)), returns helpers that
        replace f(ξ) with plain symbols for fast polynomial Ricci arithmetic
        and a chain-rule--aware derivative.  Returns None otherwise.
        """
        # Collect applied undefined functions (e.g. H(sqrt(Σyi²)))
        func_apps = set()
        for i in range(self._dim):
            func_apps.update(self._matrix[i, i].atoms(sp.core.function.AppliedUndef))
        if not func_apps:
            return None

        # Require a single shared composite argument that depends on coordinates
        coord_set = set(self._coordinates)
        args = {app.args[0] for app in func_apps
                if len(app.args) == 1 and app.args[0].free_symbols & coord_set}
        if len(args) != 1:
            return None

        xi = args.pop()                                # e.g. sqrt(Σyi²)
        xi_sq = xi.base if isinstance(xi, sp.Pow) and xi.exp == sp.Rational(1, 2) else xi**2
        r = sp.Dummy('r', positive=True)

        # Abstract symbols for each function: f, f', f''
        func_info = {}
        for app in func_apps:
            fc = app.func
            if fc not in func_info:
                nm = fc.__name__
                func_info[fc] = (sp.Dummy(nm, positive=True),
                                 sp.Dummy(nm + 'p'),
                                 sp.Dummy(nm + 'pp'))

        abstract_syms = {r}
        for f, fp, fpp in func_info.values():
            abstract_syms |= {f, fp, fpp}

        # Pre-compute ∂ξ/∂xk (in abstract form, e.g. yk/r)
        dxi = {}
        for xk in self._coordinates:
            d = sp.diff(xi, xk)
            if d != 0:
                dxi[xk] = d.subs(xi_sq, r**2).subs(xi, r)

        # ── Substitution helpers ──────────────────────────────────────

        def to_abstract(expr):
            """f(ξ) → plain symbols, resolve Subs objects."""
            out = expr.subs(xi_sq, r**2).subs(xi, r).doit()
            for fc, (f, fp, fpp) in func_info.items():
                fr = fc(r)
                out = out.subs(sp.Derivative(fr, (r, 2)), fpp)
                out = out.subs(sp.Derivative(fr, r), fp)
                out = out.subs(fr, f)
            return out

        def to_original(expr):
            """Plain symbols → f(ξ) form."""
            out = expr
            for fc, (f, fp, fpp) in func_info.items():
                fr = fc(r)
                out = out.subs(fpp, sp.Derivative(fr, (r, 2)))
                out = out.subs(fp, sp.Derivative(fr, r))
                out = out.subs(f, fr)
            return out.subs(r**2, xi_sq).subs(r, xi)

        # ── Chain-rule derivative ─────────────────────────────────────

        # Per-coordinate rules: ∂s/∂xk for each abstract symbol s
        chain_rules = {}
        for xk, dxi_dxk in dxi.items():
            rules = {r: dxi_dxk}
            for fc, (f, fp, fpp) in func_info.items():
                rules[f]  = fp  * dxi_dxk
                rules[fp] = fpp * dxi_dxk
            chain_rules[xk] = rules

        def chain_diff(expr, xk):
            """Differentiate expr (in abstract symbols) w.r.t. coordinate xk."""
            if expr == 0:
                return sp.S(0)
            present = expr.free_symbols & abstract_syms
            if not present and xk not in expr.free_symbols:
                return sp.S(0)

            # Direct partial holding abstract symbols constant
            if present:
                subs = {s: sp.Dummy() for s in present}
                inv = {v: k for k, v in subs.items()}
                direct = sp.diff(expr.subs(subs), xk).subs(inv)
            else:
                direct = sp.diff(expr, xk)

            # Chain-rule contributions
            if xk in chain_rules:
                for s, ds_dxk in chain_rules[xk].items():
                    if s in present:
                        df_ds = sp.diff(expr, s)
                        if df_ds != 0:
                            direct += df_ds * ds_dxk
            return direct

        return {
            'to_abstract': to_abstract,
            'to_original': to_original,
            'chain_diff': chain_diff,
        }

    def ricci_tensor(self, simplify_func=None):
        """Compute the Ricci tensor R_{MN}.

        For diagonal metrics, R_{MN} is diagonal — only D components are computed.
        When the metric involves composite functions like H(sqrt(Σyi²)), an
        automatic abstraction layer replaces them with plain symbols so that
        all intermediate algebra is polynomial-fast.

        Parameters
        ----------
        simplify_func : callable, optional
            Applied to each R_{MN} component (e.g., sp.cancel). Default: None.

        Returns
        -------
        sp.Matrix
            D x D matrix of Ricci tensor components.
        """
        Gamma = self.christoffel(simplify_func=simplify_func)

        # If the metric contains composite functions like H(sqrt(Σyi²)),
        # lift Christoffel symbols into an abstract polynomial representation,
        # compute Ricci there, and restore at the end.
        abstract = self._detect_composite_functions()
        if abstract is not None:
            to_abs, to_orig = abstract['to_abstract'], abstract['to_original']
            diff_func = abstract['chain_diff']
            Gamma = {k: to_abs(v) for k, v in Gamma.items()}
            if simplify_func:
                Gamma = {k: simplify_func(v) for k, v in Gamma.items() if v != 0}
        else:
            to_orig = None
            diff_func = sp.diff

        return self._ricci_loop(Gamma, diff_func, simplify_func, to_orig)

    def _ricci_loop(self, Gamma, diff_func, simplify_func, post_func):
        """Core Ricci tensor loop.

        Parameters
        ----------
        Gamma : dict
            Christoffel symbols (possibly in abstract form).
        diff_func : callable(expr, x) -> expr
            Differentiation function (sp.diff or chain-rule aware variant).
        simplify_func : callable or None
            Applied after each P-iteration and to each final component.
        post_func : callable or None
            Applied to each final component to restore original form.
        """
        D = self._dim
        x = self._coordinates
        R = sp.zeros(D, D)

        index_pairs = [(i, i) for i in range(D)] if self._is_diagonal else \
            [(i, j) for i in range(D) for j in range(i, D)]

        for M, N in index_pairs:
            val = sp.S(0)
            for P in range(D):
                # Derivative terms
                G_PMN = self._get_christoffel(Gamma, P, M, N)
                if G_PMN != 0:
                    dG = diff_func(G_PMN, x[P])
                    if dG != 0:
                        val += dG
                G_PMP = self._get_christoffel(Gamma, P, M, P)
                if G_PMP != 0:
                    dG = diff_func(G_PMP, x[N])
                    if dG != 0:
                        val -= dG
                # Quadratic terms
                for Q in range(D):
                    G_PPQ = self._get_christoffel(Gamma, P, P, Q)
                    G_QMN = self._get_christoffel(Gamma, Q, M, N)
                    if G_PPQ != 0 and G_QMN != 0:
                        val += G_PPQ * G_QMN
                    G_PNQ = self._get_christoffel(Gamma, P, N, Q)
                    G_QMP = self._get_christoffel(Gamma, Q, M, P)
                    if G_PNQ != 0 and G_QMP != 0:
                        val -= G_PNQ * G_QMP
                if simplify_func and val != 0:
                    val = simplify_func(val)

            if post_func is not None:
                val = post_func(val)
                if simplify_func:
                    val = simplify_func(val)
            R[M, N] = val
            if M != N:
                R[N, M] = val
        return R

    def ricci_scalar(self, simplify_func=None):
        """Compute the Ricci scalar R = g^{MN} R_{MN}."""
        R_tensor = self.ricci_tensor(simplify_func)
        ginv = self.inv_matrix
        R = sp.S(0)
        for M in range(self._dim):
            for N in range(self._dim):
                if ginv[M, N] != 0 and R_tensor[M, N] != 0:
                    R += ginv[M, N] * R_tensor[M, N]
        if simplify_func:
            R = simplify_func(R)
        return R

class HarmonicFunction:
    """Manages a harmonic function H(r) and substitution rules for brane ansatze.

    In a D_perp-dimensional transverse space, H satisfies:
        H'' + (D_perp - 1)/r * H' = 0

    This class provides symbolic symbols for H, H', H'', r and methods to
    substitute these into expressions, eliminating the boilerplate that is
    otherwise copy-pasted across every notebook.

    Parameters
    ----------
    transverse_dim : int
        Dimension of the transverse space.
    transverse_coords : list of sp.Symbol, optional
        Explicit transverse coordinate symbols. If None, creates y0, y1, ...
    """

    def __init__(self, transverse_dim, transverse_coords=None):
        self._d_perp = transverse_dim

        if transverse_coords is None:
            self._y = list(sp.symbols(f'y0:{transverse_dim}', real=True))
        else:
            self._y = list(transverse_coords)
            assert len(self._y) == transverse_dim

        # Radial coordinate
        self._r_expr = sp.sqrt(sum(yi**2 for yi in self._y))
        self.r = sp.Symbol('r', positive=True)

        # H and its derivatives as abstract symbols
        self.H = sp.Symbol('H', positive=True)
        self.Hp = sp.Symbol("H'", real=True)
        self.Hpp = sp.Symbol("H''", real=True)

        # H as a function of r (for differentiation)
        self._H_func = sp.Function('H')(self._r_expr)

    @property
    def transverse_dim(self):
        return self._d_perp

    @property
    def transverse_coords(self):
        return self._y

    @property
    def r_expr(self):
        """The radial coordinate as sqrt(sum yi^2)."""
        return self._r_expr

    def substitute(self, expr):
        """
        CAUTION: This method should only be invoked at the last stage of cleaning up algebraic expressions. It is destructive to further differentiation.

        Apply the full substitution chain to an expression.

        1. Replace sqrt(sum yi^2) -> r
        2. Replace H''(r) -> Hpp, H'(r) -> Hp, H(r) -> H
        3. Impose harmonic condition: Hpp -> -(D_perp - 1)/r * Hp

        Parameters
        ----------
        expr : sp.Expr
            Expression involving H(sqrt(sum yi^2)) and its derivatives.

        Returns
        -------
        sp.Expr
        """
        result = expr

        # Step 1: replace the radial expression and its square with r symbols.
        # Must substitute sum(yi^2) = r^2 *before* sqrt(sum yi^2) = r, otherwise
        # the sum-of-squares pattern lingers in denominators.
        y_sq = sum(yi**2 for yi in self._y)
        result = result.subs(y_sq, self.r**2)
        result = result.subs(self._r_expr, self.r)
        # After the r substitutions SymPy may leave Subs(Derivative(H(...), ...))
        # objects unevaluated (e.g. Subs(H'(_xi), _xi, r)).  .doit() resolves them
        # into ordinary Derivative(H(r), r) forms that the next step can match.
        result = result.doit()

        # Step 2: replace function evaluations with algebraic symbols
        H_func_r = sp.Function('H')(self.r)
        H_func_r_d1 = sp.Derivative(sp.Function('H')(self.r), self.r)
        H_func_r_d2 = sp.Derivative(sp.Function('H')(self.r), (self.r, 2))

        result = result.subs(H_func_r_d2, self.Hpp)
        result = result.subs(H_func_r_d1, self.Hp)
        result = result.subs(H_func_r, self.H)

        # Step 3: eliminate any remaining yi² using r² = sum(yi²).
        # After .doit() the sum is fragmented into individual yi² terms with
        # distributed coefficients, so a direct subs(y_sq, r²) misses them.
        # Strategy: try eliminating each yi in turn via yi² = r² - Σ(others)²
        # and pick the result with the fewest remaining transverse symbols.
        # This preserves the "natural" variable for non-isotropic expressions
        # (e.g. R[ya,ya] keeps ya) while fully simplifying isotropic ones.
        if any(yi in result.free_symbols for yi in self._y):
            best = result
            best_count = sum(1 for yi in self._y if yi in result.free_symbols)
            for yk in self._y:
                others_sq = sum(yi**2 for yi in self._y if yi is not yk)
                trial = result.subs(yk**2, self.r**2 - others_sq)
                trial = sp.cancel(trial.expand())
                count = sum(1 for yi in self._y if yi in trial.free_symbols)
                if count < best_count:
                    best = trial
                    best_count = count
                    if count == 0:
                        break
            result = best

        # Step 4: impose harmonic condition
        result = result.subs(self.Hpp, -(self._d_perp - 1) / self.r * self.Hp)

        return result

    def random_values(self, rng=None):
        """Generate random numerical values consistent with the harmonic condition.

        Returns a dict mapping all symbols (H, H', H'', r, y0, y1, ...) to
        float values, with H'' = -(D_perp - 1)/r * H' enforced.

        Parameters
        ----------
        rng : random.Random, optional
            Random number generator. Uses a new one if None.
        """
        import random
        if rng is None:
            rng = random.Random()

        r_val = rng.uniform(0.5, 3.0)
        H_val = rng.uniform(0.5, 3.0)
        Hp_val = rng.uniform(-1.0, 1.0)
        Hpp_val = -(self._d_perp - 1) / r_val * Hp_val

        # Generate yi values on the sphere of radius r_val
        raw = [rng.gauss(0, 1) for _ in range(self._d_perp)]
        norm = (sum(x**2 for x in raw)) ** 0.5
        y_vals = [r_val * x / norm for x in raw]

        values = {
            self.r: r_val,
            self.H: H_val,
            self.Hp: Hp_val,
            self.Hpp: Hpp_val,
        }
        for yi, vi in zip(self._y, y_vals):
            values[yi] = vi

        return values


def warped_product(warp_factors, block_dims, block_signatures, coordinates, H=None):
    """Build a warped product metric.

    ds^2 = w_1 ds_{d_1}^2 + w_2 ds_{d_2}^2 + ...

    Parameters
    ----------
    warp_factors : list of sp.Expr
        Warp factor for each block (e.g., H**a, H**b).
    block_dims : list of int
        Dimension of each block.
    block_signatures : list of str
        'lorentzian' or 'euclidean' for each block.
    coordinates : list of sp.Symbol
        All D coordinates (must have sum(block_dims) entries).
    H : sp.Expr, optional
        Unused; kept for backward compatibility.

    Returns
    -------
    Metric
    """
    assert len(warp_factors) == len(block_dims) == len(block_signatures)
    D = sum(block_dims)
    assert len(coordinates) == D

    diagonal = []
    for warp, d, sig in zip(warp_factors, block_dims, block_signatures):
        if sig == 'lorentzian':
            diagonal.append(-warp)  # time component
            for _ in range(d - 1):
                diagonal.append(warp)
        elif sig == 'euclidean':
            for _ in range(d):
                diagonal.append(warp)
        else:
            raise ValueError(f"Unknown signature: {sig}. Use 'lorentzian' or 'euclidean'.")

    return Metric(diagonal, coordinates)
