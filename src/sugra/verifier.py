"""
Einstein equation verifier: R_{MN} = T_{MN}.
"""

import sympy as sp
from .forms import form_stress_energy


class Verifier:
    """Verify R_{MN} = T_{MN} for a brane solution.

    Parameters
    ----------
    soln : dict
        Keys: metric, F, Phi, alpha, coords, hf
    """

    def __init__(self, soln):
        self.g      = soln['metric']
        self.F      = soln['F']
        self.Phi    = soln['Phi']
        self.alpha  = soln['alpha']
        self.coords = soln['coords']
        self.hf     = soln['hf']

        self._R = None
        self._T = None

    # ── Formatting helpers ──────────────────────────────────────────────────

    @staticmethod
    def _sub(n):
        t = str.maketrans('0123456789',
                          '\u2080\u2081\u2082\u2083\u2084\u2085\u2086\u2087\u2088\u2089')
        return str(n).translate(t)

    @staticmethod
    def _sup(n):
        t = str.maketrans('0123456789',
                          '\u2070\u00b9\u00b2\u00b3\u2074\u2075\u2076\u2077\u2078\u2079')
        return str(n).translate(t)

    @staticmethod
    def _header(width=56):
        return '\u2550' * width

    @staticmethod
    def _rule(width=56):
        return '\u2500' * width

    def _coord_label(self, i):
        return str(self.coords[i])

    def _coord_width(self):
        return max(len(self._coord_label(i)) for i in range(self.g.dim))

    def _display(self, expr):
        """Substitute H(sqrt(..)) -> H and sqrt(..) -> r for readable display.

        Unlike hf.substitute(), this does NOT impose the harmonic condition
        or eliminate transverse coordinates, so it's safe for presentation only.
        """
        hf = self.hf
        r = hf.r
        H = hf.H
        Hp = hf.Hp

        # Replace sqrt(sum yi^2) -> r, sum yi^2 -> r^2
        y_sq = sum(yi**2 for yi in hf._y)
        result = expr.subs(y_sq, r**2).subs(hf.r_expr, r)
        result = result.doit()

        # Replace H(r) -> H, H'(r) -> H'
        H_func_r = sp.Function('H')(r)
        result = result.subs(sp.Derivative(H_func_r, r), Hp)
        result = result.subs(H_func_r, H)

        return result

    # ── Display: action ─────────────────────────────────────────────────────

    def action_str(self):
        D = self.g.dim
        n = self.F.rank
        a = self.alpha
        if a == 0:
            coupling = ''
        elif a == 1:
            coupling = 'e^\u03a6 '
        elif a == -1:
            coupling = 'e^(-\u03a6) '
        else:
            coupling = f'e^({a}\u03a6) '
        return (f"S = \u222b d{self._sup(D)}x \u221a(-g) \u00b7 "
                f"[R - \u00bd(\u2202\u03a6)\u00b2 - \u00bd{coupling}|F{self._sub(n)}|\u00b2]")

    def print_action(self):
        print(f"  {self.action_str()}")

    # ── Display: metric ─────────────────────────────────────────────────────

    def print_metric(self):
        g = self.g
        x = self.coords
        if g.is_diagonal:
            print("  Metric (diagonal):")
            for i in range(g.dim):
                print(f"    g[{x[i]},{x[i]}] = {self._display(g.matrix[i, i])}")
        else:
            print("  Metric:")
            for i in range(g.dim):
                for j in range(i, g.dim):
                    val = g.matrix[i, j]
                    if val != 0:
                        print(f"    g[{x[i]},{x[j]}] = {self._display(val)}")

    # ── Display: form field ─────────────────────────────────────────────────

    def _distinct_components(self):
        """Return one representative per distinct component value (up to sign)."""
        comps = self.F.nonzero_components
        seen = {}
        for idx, val in comps.items():
            display_val = self._display(val)
            canonical = sp.cancel(sp.Abs(display_val))
            if canonical not in seen:
                seen[canonical] = (idx, display_val)
        return seen

    def print_form(self):
        comps = self.F.nonzero_components
        distinct = self._distinct_components()
        x = self.coords
        print(f"  F ({self.F.rank}-form), "
              f"{len(comps)} independent component(s), "
              f"{len(distinct)} distinct:")
        for idx, val in distinct.values():
            labels = ','.join(str(x[i]) for i in idx)
            print(f"    F[{labels}] = {self._display(val)}")

    # ── Display: dilaton ────────────────────────────────────────────────────

    def print_dilaton(self):
        if self.Phi == 0:
            print(f"  \u03a6 = 0 (no dilaton)")
        else:
            print(f"  \u03a6 = {self._display(self.Phi)}")
        print(f"  \u03b1 = {self.alpha}")

    # ── Display: full ansatz ────────────────────────────────────────────────

    def print_ansatz(self):
        self.print_metric()
        print()
        self.print_form()
        print()
        self.print_dilaton()

    # ── Compute ─────────────────────────────────────────────────────────────

    def compute(self):
        self._R = self.g.ricci_tensor(simplify_func=sp.cancel)
        self._T = form_stress_energy(
            self.F, self.g, dilaton=self.Phi, dilaton_coupling=self.alpha
        )

    # ── Display: component-wise results ─────────────────────────────────────

    def print_results(self):
        if self._R is None:
            self.compute()

        D = self.g.dim
        w = self._coord_width()
        pad = ' ' * (w + 4)

        all_ok = True
        for i in range(D):
            R_expr = self.hf.substitute(sp.cancel(self._R[i, i]))
            T_expr = self.hf.substitute(sp.cancel(self._T[i, i]))
            diff = sp.cancel(R_expr - T_expr)
            name = self._coord_label(i).ljust(w)
            ok = (diff == 0)
            if not ok:
                all_ok = False
            tag = "\u2713" if ok else "\u2717"
            print(f"  {tag} {name}   R = {R_expr}")
            print(f"    {pad} T = {T_expr}")
            if not ok:
                print(f"    {pad} \u0394 = {diff}")

        return all_ok

    # ── Main entry point ────────────────────────────────────────────────────

    def check(self):
        print()
        print(f"  {self._header()}")
        print(f"  Verifying  R_MN = T_MN")
        print(f"  {self._header()}")
        print()
        self.print_action()
        print()
        self.print_ansatz()
        print()
        print(f"  {self._rule()}")
        print(f"  Computing R_MN and T_MN ...")
        print()

        all_ok = self.print_results()

        print()
        print(f"  {self._rule()}")
        verdict = ("\u2713 All components match." if all_ok
                   else "\u2717 MISMATCH detected.")
        print(f"  {verdict}")
        print()

        return all_ok
