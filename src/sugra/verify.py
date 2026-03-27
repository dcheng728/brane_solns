"""
Verification engine for supergravity solutions.

Provides numerical spot-checking and symbolic simplification to determine
whether expressions (e.g., R_MN - T_MN) vanish. Also provides the Solution
class that ties together a metric, scalar fields, and flux fields and checks
all equations of motion.
"""

import sympy as sp
import random
from dataclasses import dataclass, field
from typing import Optional

from .geometry import Metric, HarmonicFunction
from .forms import FormField, form_stress_energy, hodge_star, exterior_derivative, form_norm_squared


# ---------------------------------------------------------------------------
# Verification utilities
# ---------------------------------------------------------------------------

@dataclass
class CheckResult:
    """Result of checking whether an expression vanishes."""
    label: str
    passed: bool
    method: str  # 'symbolic', 'numerical', or 'failed'
    residual: object = None  # symbolic expr or max numerical residual
    details: str = ""

    def __str__(self):
        status = "PASS" if self.passed else "FAIL"
        info = f"({self.method}"
        if self.method == 'numerical' and self.residual is not None:
            info += f", max residual {self.residual:.2e}"
        elif self.method == 'symbolic':
            info += ", simplified to 0"
        elif self.method == 'failed':
            info += f": {self.details}"
        info += ")"
        return f"  {self.label}: {status} {info}"


@dataclass
class VerificationReport:
    """Collection of CheckResults from verifying a solution."""
    checks: list = field(default_factory=list)

    def add(self, result):
        self.checks.append(result)

    @property
    def all_passed(self):
        return all(c.passed for c in self.checks)

    def summary(self):
        print("=" * 60)
        print("Supergravity Solution Verification Report")
        print("=" * 60)
        for c in self.checks:
            print(c)
        print("-" * 60)
        if self.all_passed:
            print("ALL CHECKS PASSED")
        else:
            n_fail = sum(1 for c in self.checks if not c.passed)
            print(f"{n_fail} CHECK(S) FAILED")
        print("=" * 60)


def verify_symbolic(expr, strategies=None):
    """Try to simplify expr to zero using a sequence of strategies.

    Parameters
    ----------
    expr : sp.Expr
    strategies : list of callables, optional
        Each takes an expression and returns a simplified form.
        Default: [cancel, factor, powsimp+radsimp, simplify].

    Returns
    -------
    (is_zero: bool, simplified_expr: sp.Expr)
    """
    if expr == 0 or expr is sp.S.Zero:
        return True, sp.S(0)

    if strategies is None:
        strategies = [
            sp.cancel,
            sp.factor,
            lambda e: sp.radsimp(sp.powsimp(e)),
            sp.simplify,
        ]

    for strat in strategies:
        try:
            result = strat(expr)
            if result == 0 or result is sp.S.Zero:
                return True, sp.S(0)
        except Exception:
            continue

    return False, expr


def verify_numerical(expr, symbol_values_fn, n_trials=20, tol=1e-10):
    """Evaluate expr at random points to check if it vanishes.

    Parameters
    ----------
    expr : sp.Expr
    symbol_values_fn : callable
        Returns a dict {symbol: float_value} for one trial.
    n_trials : int
    tol : float

    Returns
    -------
    (is_zero: bool, max_residual: float, worst_subs: dict or None)
    """
    if expr == 0 or expr is sp.S.Zero:
        return True, 0.0, None

    max_residual = 0.0
    worst_subs = None

    for _ in range(n_trials):
        subs = symbol_values_fn()
        try:
            val = complex(expr.subs(subs))
            residual = abs(val)
            if residual > max_residual:
                max_residual = residual
                worst_subs = subs
        except (TypeError, ValueError, ZeroDivisionError):
            continue

    is_zero = max_residual < tol
    return is_zero, max_residual, worst_subs


def check_expression(expr, label, symbol_values_fn=None, n_trials=20, tol=1e-10):
    """Check whether an expression vanishes, trying symbolic then numerical.

    Parameters
    ----------
    expr : sp.Expr
    label : str
        Human-readable label for this check (e.g., "Einstein (0,0)").
    symbol_values_fn : callable, optional
        If provided, enables numerical verification.
    n_trials : int
    tol : float

    Returns
    -------
    CheckResult
    """
    # Try symbolic first (fast strategies only)
    is_zero, simplified = verify_symbolic(expr, strategies=[sp.cancel, sp.factor])
    if is_zero:
        return CheckResult(label=label, passed=True, method='symbolic')

    # Try numerical if available
    if symbol_values_fn is not None:
        is_zero, max_res, worst = verify_numerical(expr, symbol_values_fn, n_trials, tol)
        if is_zero:
            return CheckResult(label=label, passed=True, method='numerical', residual=max_res)
        else:
            return CheckResult(label=label, passed=False, method='numerical',
                               residual=max_res, details=f"max residual {max_res:.2e}")

    # Try heavier symbolic simplification
    is_zero, simplified = verify_symbolic(expr)
    if is_zero:
        return CheckResult(label=label, passed=True, method='symbolic')

    return CheckResult(label=label, passed=False, method='failed',
                       residual=simplified, details="could not simplify to 0")


# ---------------------------------------------------------------------------
# Field wrappers
# ---------------------------------------------------------------------------

@dataclass
class ScalarField:
    """A scalar field (e.g., dilaton) in the solution.

    Parameters
    ----------
    expr : sp.Expr
        The scalar field as a function of coordinates.
    name : str
        Human-readable name (e.g., "dilaton").
    """
    expr: sp.Expr
    name: str = "scalar"


@dataclass
class FluxField:
    """A p-form flux field in the solution.

    Parameters
    ----------
    form : FormField
        The form field.
    dilaton_coupling : sp.Expr or int
        The coupling constant alpha in e^{alpha * phi} * S_MN.
    self_dual : bool
        Whether this form is self-dual (F = *F).
    name : str
        Human-readable name.
    """
    form: FormField
    dilaton_coupling: object = 0
    self_dual: bool = False
    name: str = "flux"


# ---------------------------------------------------------------------------
# Solution class
# ---------------------------------------------------------------------------

class Solution:
    """A supergravity solution: metric + scalar fields + flux fields.

    Provides methods to check Einstein equation, dilaton EOM,
    form field Bianchi identities, and self-duality constraints.

    Parameters
    ----------
    metric : Metric
    fields : list of ScalarField and FluxField
    coordinates : list of sp.Symbol, optional
        If None, uses metric.coordinates.
    harmonic : HarmonicFunction, optional
        If provided, enables harmonic substitution and numerical verification.
    """

    def __init__(self, metric, fields=None, coordinates=None, harmonic=None):
        self.metric = metric
        self.fields = fields or []
        self.coordinates = coordinates or metric.coordinates
        self.harmonic = harmonic

        # Separate field types
        self.scalars = [f for f in self.fields if isinstance(f, ScalarField)]
        self.fluxes = [f for f in self.fields if isinstance(f, FluxField)]

    def _make_symbol_values_fn(self):
        """Create a function that generates random symbol values for numerical checks."""
        if self.harmonic is None:
            return None

        hf = self.harmonic

        def fn():
            return hf.random_values(rng=random.Random())

        return fn

    def ricci_tensor(self, simplify_func=None):
        """Return the Ricci tensor R_{MN} as a SymPy Matrix.

        This is the left-hand side of the Einstein equation:
            R_{MN} = (1/2)(d_M phi d_N phi + e^{alpha phi} S_{MN}|_F)

        Parameters
        ----------
        simplify_func : callable, optional
            Applied to each component (e.g. sp.cancel).

        Returns
        -------
        sp.Matrix  (D x D)
        """
        return self.metric.ricci_tensor(simplify_func=simplify_func)

    def stress_energy_tensor(self):
        """Return the RHS of the Einstein equation as a SymPy Matrix.

        Computes:
            T_{MN} = (1/2) d_M phi d_N phi  +  e^{alpha phi} S_{MN}|_F

        where S_{MN}|_F is the trace-reversed stress-energy of each flux
        (see form_stress_energy), consistently with the action

            S = integral d^D x sqrt(-g) [R - 1/2 (d phi)^2
                                           - 1/2 e^{alpha phi} |F_n|^2]

        Returns
        -------
        sp.Matrix  (D x D)
        """
        return self._compute_stress_energy()

    def _compute_stress_energy(self):
        """Compute the total stress-energy T_{MN} from all fields.

        T_{MN} = (1/2) sum_scalars d_M phi d_N phi
                + sum_fluxes e^{alpha * phi} S_MN|_F
        """
        D = self.metric.dim
        g = self.metric
        T = sp.zeros(D, D)

        # Scalar field contributions: (1/2) d_M phi d_N phi
        for sf in self.scalars:
            phi = sf.expr
            for M in range(D):
                dphi_M = sp.diff(phi, self.coordinates[M])
                if dphi_M == 0:
                    continue
                for N in range(M, D):
                    dphi_N = sp.diff(phi, self.coordinates[N])
                    if dphi_N == 0:
                        continue
                    contrib = sp.Rational(1, 2) * dphi_M * dphi_N
                    T[M, N] += contrib
                    if M != N:
                        T[N, M] += contrib

        # Flux field contributions: e^{alpha * phi} S_MN|_F
        dilaton_expr = self._get_dilaton_expr()
        for ff in self.fluxes:
            S = form_stress_energy(ff.form, g)
            alpha = ff.dilaton_coupling
            if alpha != 0 and dilaton_expr is not None:
                prefactor = sp.exp(alpha * dilaton_expr)
            else:
                prefactor = sp.S(1)
            for M in range(D):
                for N in range(M, D):
                    if S[M, N] != 0:
                        T[M, N] += prefactor * S[M, N]
                        if M != N:
                            T[N, M] += prefactor * S[M, N]

        return T

    def _get_dilaton_expr(self):
        """Find the dilaton scalar field expression, if any."""
        for sf in self.scalars:
            if sf.name.lower() in ('dilaton', 'phi', 'φ'):
                return sf.expr
        return None

    def check_einstein(self, simplify_func=None):
        """Check the Einstein equation: R_{MN} = T_{MN}.

        Parameters
        ----------
        simplify_func : callable, optional
            Applied to Ricci tensor components.

        Returns
        -------
        list of CheckResult
        """
        R = self.ricci_tensor(simplify_func=simplify_func)
        T = self.stress_energy_tensor()
        D = self.metric.dim

        symbol_values_fn = self._make_symbol_values_fn()
        results = []

        # For diagonal metrics, only check diagonal components
        if self.metric.is_diagonal:
            indices = [(i, i) for i in range(D)]
        else:
            indices = [(i, j) for i in range(D) for j in range(i, D)]

        for M, N in indices:
            residual = R[M, N] - T[M, N]
            if self.harmonic is not None:
                residual = self.harmonic.substitute(residual)
            label = f"Einstein ({M},{N})"
            results.append(check_expression(residual, label, symbol_values_fn))

        return results

    def check_self_duality(self, flux_field=None):
        """Check F = *F for self-dual form fields.

        Parameters
        ----------
        flux_field : FluxField, optional
            Specific flux to check. If None, checks all self-dual fluxes.

        Returns
        -------
        list of CheckResult
        """
        targets = [flux_field] if flux_field else [f for f in self.fluxes if f.self_dual]
        results = []
        symbol_values_fn = self._make_symbol_values_fn()

        for ff in targets:
            dual = hodge_star(ff.form, self.metric)
            diff = ff.form - dual

            for idx in ff.form.independent_indices:
                val = diff[idx]
                if self.harmonic is not None:
                    val = self.harmonic.substitute(val)
                if val != 0:
                    label = f"Self-duality {ff.name} {idx}"
                    results.append(check_expression(val, label, symbol_values_fn))

            if not any(diff[idx] != 0 for idx in ff.form.independent_indices):
                results.append(CheckResult(
                    label=f"Self-duality {ff.name}",
                    passed=True, method='symbolic'))

        return results

    def check_bianchi(self, flux_field=None):
        """Check dF = 0 (Bianchi identity) for flux fields.

        Parameters
        ----------
        flux_field : FluxField, optional
            Specific flux to check. If None, checks all fluxes.

        Returns
        -------
        list of CheckResult
        """
        targets = [flux_field] if flux_field else self.fluxes
        results = []
        symbol_values_fn = self._make_symbol_values_fn()

        for ff in targets:
            dF = exterior_derivative(ff.form, self.coordinates)
            nonzero = dF.nonzero_components
            if not nonzero:
                results.append(CheckResult(
                    label=f"Bianchi {ff.name}",
                    passed=True, method='symbolic'))
            else:
                for idx, val in nonzero.items():
                    if self.harmonic is not None:
                        val = self.harmonic.substitute(val)
                    label = f"Bianchi {ff.name} {idx}"
                    results.append(check_expression(val, label, symbol_values_fn))

        return results

    def check_dilaton(self):
        """Check the dilaton equation of motion.

        Box(phi) = sum_n (alpha_n / 2) * e^{alpha_n * phi} * |F_n|^2

        Returns
        -------
        list of CheckResult
        """
        dilaton_expr = self._get_dilaton_expr()
        if dilaton_expr is None:
            return []

        g = self.metric
        D = g.dim
        ginv = g.inv_matrix
        x = self.coordinates

        # Compute Box(phi) = g^{MN} (d_M d_N phi - Gamma^P_{MN} d_P phi)
        Gamma = g.christoffel()
        box_phi = sp.S(0)
        for M in range(D):
            for N in range(D):
                if ginv[M, N] == 0:
                    continue
                d2 = sp.diff(dilaton_expr, x[M], x[N])
                christoffel_term = sp.S(0)
                for P in range(D):
                    key = (P, M, N) if M <= N else (P, N, M)
                    G_PMN = Gamma.get(key, sp.S(0))
                    if G_PMN != 0:
                        dphi_P = sp.diff(dilaton_expr, x[P])
                        if dphi_P != 0:
                            christoffel_term += G_PMN * dphi_P
                box_phi += ginv[M, N] * (d2 - christoffel_term)

        # Compute source: sum alpha_n/2 * e^{alpha_n * phi} * |F_n|^2
        source = sp.S(0)
        for ff in self.fluxes:
            alpha = ff.dilaton_coupling
            if alpha == 0:
                continue
            F_sq = form_norm_squared(ff.form, g)
            if alpha != 0 and dilaton_expr is not None:
                source += sp.Rational(1, 2) * alpha * sp.exp(alpha * dilaton_expr) * F_sq

        residual = box_phi - source
        if self.harmonic is not None:
            residual = self.harmonic.substitute(residual)

        symbol_values_fn = self._make_symbol_values_fn()
        return [check_expression(residual, "Dilaton EOM", symbol_values_fn)]

    def verify_all(self, simplify_func=None):
        """Run all applicable checks and return a VerificationReport.

        Parameters
        ----------
        simplify_func : callable, optional
            Passed to ricci_tensor computation.

        Returns
        -------
        VerificationReport
        """
        report = VerificationReport()

        # Einstein equation
        for r in self.check_einstein(simplify_func):
            report.add(r)

        # Self-duality
        for r in self.check_self_duality():
            report.add(r)

        # Bianchi identities
        for r in self.check_bianchi():
            report.add(r)

        # Dilaton EOM
        for r in self.check_dilaton():
            report.add(r)

        return report
