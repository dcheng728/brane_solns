"""
Einstein equation verifier: R_{MN} = T_{MN}.
"""

import sympy as sp
from .forms import form_stress_energy, scalar_stress_energy


class Verifier:
    """Verify R_{MN} = T_{MN} for a brane solution.

    Parameters
    ----------
    soln : dict
        Keys: metric, forms, Phi, coords, hf
        ``forms`` is a list of [FormField, alpha] pairs, e.g. [[F4, 0], [F5, 0]].
    """

    def __init__(self, soln):
        self.g      = soln['metric']
        self.forms  = [tuple(pair) for pair in soln['forms']]
        self.Phi    = soln['Phi']
        self.coords = soln['coords']
        self.hf     = soln['hf']
        self.D_trace = soln.get('D_trace', None)

        self._R = None
        self._T = None

    # -- Formatting helpers --------------------------------------------------

    @staticmethod
    def _header(width=56):
        return '=' * width

    @staticmethod
    def _rule(width=56):
        return '-' * width

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

    # -- Display: action -----------------------------------------------------

    @staticmethod
    def _coupling_str(alpha):
        if alpha == 0:
            return ''
        elif alpha == 1:
            return 'e^Phi '
        elif alpha == -1:
            return 'e^(-Phi) '
        else:
            return f'e^({alpha} Phi) '

    def action_str(self):
        D = self.g.dim
        form_terms = ' '.join(
            f"- 1/2 {self._coupling_str(entry[1])}|F_{entry[0].rank}|^2"
            for entry in self.forms
        )
        dilaton_term = '- 1/2 (dPhi)^2 ' if self.Phi != 0 else ''
        return f"S = int d^{D}x sqrt(-g) [R {dilaton_term}{form_terms}]"

    def print_action(self):
        print(f"  {self.action_str()}")

    # -- Display: metric -----------------------------------------------------

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

    # -- Display: form fields ------------------------------------------------

    def _distinct_components(self, F):
        """Return one representative per distinct component value (up to sign)."""
        seen = {}
        for idx, val in F.nonzero_components.items():
            display_val = self._display(val)
            canonical = sp.cancel(sp.Abs(display_val))
            if canonical not in seen:
                seen[canonical] = (idx, display_val)
        return seen

    def print_forms(self):
        x = self.coords
        for k, entry in enumerate(self.forms):
            F, alpha = entry[0], entry[1]
            n_eff = entry[2] if len(entry) > 2 else F.rank
            comps = F.nonzero_components
            distinct = self._distinct_components(F)
            coupling = f", alpha={alpha}" if alpha != 0 else ""
            neff_str = f", n_eff={n_eff}" if n_eff != F.rank else ""
            print(f"  F{k+1} ({F.rank}-form{coupling}{neff_str}), "
                  f"{len(comps)} independent component(s), "
                  f"{len(distinct)} distinct:")
            for idx, val in distinct.values():
                labels = ','.join(str(x[i]) for i in idx)
                print(f"    F[{labels}] = {val}")

    # -- Display: dilaton ----------------------------------------------------

    def print_dilaton(self):
        if self.Phi == 0:
            print("  Phi = 0 (no dilaton)")
        else:
            print(f"  Phi = {self._display(self.Phi)}")

    # -- Display: full ansatz ------------------------------------------------

    def print_ansatz(self):
        self.print_metric()
        print()
        self.print_forms()
        print()
        self.print_dilaton()

    # -- Compute -------------------------------------------------------------

    def compute(self):
        self._R = self.g.ricci_tensor(simplify_func=sp.cancel)
        D = self.g.dim
        self._T = sp.zeros(D, D)
        for entry in self.forms:
            F, alpha = entry[0], entry[1]
            n_eff = entry[2] if len(entry) > 2 else None
            self._T += form_stress_energy(F, self.g, dilaton=self.Phi, alpha=alpha,
                                          D_trace=self.D_trace, n_eff=n_eff)
        if self.Phi != 0:
            self._T += scalar_stress_energy(self.Phi, self.g)

    # -- Display: component-wise results -------------------------------------

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
            tag = "OK" if ok else "XX"
            print(f"  [{tag}] {name}   R = {R_expr}")
            print(f"        {pad} T = {T_expr}")
            if not ok:
                print(f"        {pad} D = {diff}")

        return all_ok

    # -- Main entry point ----------------------------------------------------

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
        verdict = "All components match." if all_ok else "MISMATCH detected."
        print(f"  {verdict}")
        print()

        return all_ok
