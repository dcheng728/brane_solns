"""Concrete geometry: coordinates + metric -> Christoffel, Riemann, etc.

A Geometry is defined by coordinates and a metric. Everything else
is computed from the defining formulas.

    geom = Geometry(
        coordinates=('r', 'theta'),
        metric={('r', 'r'): 1, ('theta', 'theta'): r**2}
    )
    geom.christoffel()
    geom.riemann()
    geom.ricci()
    geom.ricci_scalar()
    geom.geodesic_equations()
"""

import sympy as sp


class Geometry:
    """A concrete Riemannian geometry: coordinates + metric.

    All curvature quantities computed from the defining formulas.

    Parameters
    ----------
    coordinates : tuple of str or tuple of sp.Symbol
        Coordinate names.
    metric : dict
        Maps (name_i, name_j) -> SymPy expression.
        Symmetric: only need upper triangle.
    """

    def __init__(self, coordinates, metric):
        # Parse coordinates
        if all(isinstance(c, str) for c in coordinates):
            self.names = list(coordinates)
            self.x = [sp.Symbol(n) for n in self.names]
        else:
            self.x = list(coordinates)
            self.names = [str(s) for s in self.x]

        self.dim = len(self.names)

        # Build metric matrix
        self.g = sp.zeros(self.dim, self.dim)
        for (ni, nj), val in metric.items():
            i = self._idx(ni)
            j = self._idx(nj)
            self.g[i, j] = val
            self.g[j, i] = val

        self.g_inv = self.g.inv()

        # Caches
        self._gamma = None
        self._riem = None
        self._ric = None

    def _idx(self, name):
        """Convert coordinate name or symbol to integer index."""
        if isinstance(name, str):
            return self.names.index(name)
        return self.x.index(name) if name in self.x else self.names.index(str(name))

    # ── Christoffel ─────────────────────────────────────────────────

    def christoffel(self):
        """Gamma^c_{ab} = 1/2 g^{cd}(d_a g_{db} + d_b g_{ad} - d_d g_{ab}).

        Returns
        -------
        dict : (c, a, b) -> SymPy expression (nonzero only, names as strings)
        """
        gamma = self._christoffel_array()
        result = {}
        for c in range(self.dim):
            for a in range(self.dim):
                for b in range(a, self.dim):
                    val = gamma[c][a][b]
                    if val != 0:
                        key = (self.names[c], self.names[a], self.names[b])
                        result[key] = val
                        if a != b:
                            result[(self.names[c], self.names[b], self.names[a])] = val
        return result

    # ── Riemann ─────────────────────────────────────────────────────

    def riemann(self):
        """R^a_{bcd} from Christoffel symbols.

        Returns
        -------
        dict : (a, b, c, d) -> SymPy expression (nonzero only)
        """
        gamma = self._christoffel_array()
        result = {}
        for a in range(self.dim):
            for b in range(self.dim):
                for c in range(self.dim):
                    for d in range(self.dim):
                        val = (sp.diff(gamma[a][d][b], self.x[c])
                               - sp.diff(gamma[a][c][b], self.x[d]))
                        for e in range(self.dim):
                            val += (gamma[a][c][e] * gamma[e][d][b]
                                    - gamma[a][d][e] * gamma[e][c][b])
                        val = sp.simplify(val)
                        if val != 0:
                            key = tuple(self.names[i] for i in [a, b, c, d])
                            result[key] = val
        return result

    # ── Ricci ───────────────────────────────────────────────────────

    def ricci(self):
        """Ric_{ab} = R^c_{acb}.

        Returns
        -------
        dict : (a, b) -> SymPy expression (nonzero only)
        """
        riem = self._riemann_array()
        result = {}
        for a in range(self.dim):
            for b in range(a, self.dim):
                val = sp.S(0)
                for c in range(self.dim):
                    val += riem[c][a][c][b]
                val = sp.simplify(val)
                if val != 0:
                    key = (self.names[a], self.names[b])
                    result[key] = val
                    if a != b:
                        result[(self.names[b], self.names[a])] = val
        return result

    # ── Ricci scalar ────────────────────────────────────────────────

    def ricci_scalar(self):
        """R = g^{ab} Ric_{ab}.

        Returns
        -------
        SymPy expression
        """
        ric = self._ricci_array()
        val = sp.S(0)
        for a in range(self.dim):
            for b in range(self.dim):
                val += self.g_inv[a, b] * ric[a][b]
        return sp.simplify(val)

    # ── Geodesic equations ──────────────────────────────────────────

    def geodesic_equations(self):
        """ddot{x}^a + Gamma^a_{bc} dot{x}^b dot{x}^c = 0.

        Returns
        -------
        list of (coord_name, SymPy expression) — each expr = 0
        """
        gamma = self._christoffel_array()
        x_dot = [sp.Symbol(f'dot_{n}') for n in self.names]
        x_ddot = [sp.Symbol(f'ddot_{n}') for n in self.names]

        eqs = []
        for a in range(self.dim):
            eq = x_ddot[a]
            for b in range(self.dim):
                for c in range(self.dim):
                    eq += gamma[a][b][c] * x_dot[b] * x_dot[c]
            eqs.append((self.names[a], sp.simplify(eq)))
        return eqs

    # ── Internal ────────────────────────────────────────────────────

    def _christoffel_array(self):
        if self._gamma is not None:
            return self._gamma
        n = self.dim
        gamma = [[[sp.S(0)] * n for _ in range(n)] for _ in range(n)]
        for c in range(n):
            for a in range(n):
                for b in range(n):
                    val = sp.S(0)
                    for d in range(n):
                        val += sp.Rational(1, 2) * self.g_inv[c, d] * (
                            sp.diff(self.g[d, b], self.x[a])
                            + sp.diff(self.g[a, d], self.x[b])
                            - sp.diff(self.g[a, b], self.x[d])
                        )
                    gamma[c][a][b] = sp.simplify(val)
        self._gamma = gamma
        return gamma

    def _riemann_array(self):
        if self._riem is not None:
            return self._riem
        gamma = self._christoffel_array()
        n = self.dim
        riem = [[[[sp.S(0)] * n for _ in range(n)] for _ in range(n)] for _ in range(n)]
        for a in range(n):
            for b in range(n):
                for c in range(n):
                    for d in range(n):
                        val = (sp.diff(gamma[a][d][b], self.x[c])
                               - sp.diff(gamma[a][c][b], self.x[d]))
                        for e in range(n):
                            val += (gamma[a][c][e] * gamma[e][d][b]
                                    - gamma[a][d][e] * gamma[e][c][b])
                        riem[a][b][c][d] = sp.simplify(val)
        self._riem = riem
        return riem

    def _ricci_array(self):
        if self._ric is not None:
            return self._ric
        riem = self._riemann_array()
        n = self.dim
        ric = [[sp.S(0)] * n for _ in range(n)]
        for a in range(n):
            for b in range(n):
                val = sp.S(0)
                for c in range(n):
                    val += riem[c][a][c][b]
                ric[a][b] = sp.simplify(val)
        self._ric = ric
        return ric

    def __repr__(self):
        return f"Geometry({', '.join(self.names)})"
