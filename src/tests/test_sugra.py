"""Tests for the sugra library."""

import unittest
import sympy as sp
from sympy import symbols, Rational, Function, sqrt, diff, diag, simplify, cancel

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from sugra.geometry import Metric, HarmonicFunction, warped_product
from sugra.forms import (
    FormField, epsilon, hodge_star, form_norm_squared,
    form_contraction, form_stress_energy, exterior_derivative,
)
from sugra.verify import (
    verify_symbolic, verify_numerical, check_expression,
    Solution, ScalarField, FluxField,
)


# =========================================================================
# Geometry tests
# =========================================================================

class TestMetricBasics(unittest.TestCase):
    """Test Metric construction and properties."""

    def test_diagonal_detection(self):
        g = Metric([1, 1, 1], symbols('x y z'))
        self.assertTrue(g.is_diagonal)
        self.assertEqual(g.dim, 3)

    def test_general_metric(self):
        x, y = symbols('x y')
        m = sp.Matrix([[1, x], [x, 1 + x**2]])
        g = Metric(m, [x, y])
        self.assertFalse(g.is_diagonal)

    def test_inv_metric_diagonal(self):
        a, b = symbols('a b', positive=True)
        g = Metric([a, b], symbols('x y'))
        self.assertEqual(g.inv_matrix[0, 0], 1/a)
        self.assertEqual(g.inv_matrix[1, 1], 1/b)

    def test_det_diagonal(self):
        a, b, c = symbols('a b c', positive=True)
        g = Metric([a, b, c], symbols('x y z'))
        self.assertEqual(g.det, a * b * c)


class TestFlatSpaceRicci(unittest.TestCase):
    """Flat Minkowski/Euclidean metrics should have R_MN = 0."""

    def test_minkowski_4d(self):
        coords = symbols('t x y z')
        g = Metric([-1, 1, 1, 1], coords)
        R = g.ricci_tensor()
        for i in range(4):
            for j in range(4):
                self.assertEqual(R[i, j], 0)

    def test_euclidean_3d(self):
        coords = symbols('x y z')
        g = Metric([1, 1, 1], coords)
        R = g.ricci_tensor()
        for i in range(3):
            self.assertEqual(R[i, i], 0)

    def test_minkowski_10d(self):
        coords = symbols('x0:10')
        g = Metric([-1] + [1]*9, coords)
        R = g.ricci_tensor()
        for i in range(10):
            self.assertEqual(R[i, i], 0)


class TestSchwarzschildRicci(unittest.TestCase):
    """Schwarzschild metric should have R_MN = 0 (vacuum solution)."""

    def test_schwarzschild_4d(self):
        """Full 4D Schwarzschild: ds^2 = -f dt^2 + f^{-1} dr^2 + r^2 dOmega^2.

        R_MN = 0 for all components.
        """
        t, r, theta, phi = symbols('t r theta phi')
        M = symbols('M', positive=True)
        f = 1 - 2*M/r

        g = Metric([-f, 1/f, r**2, r**2 * sp.sin(theta)**2],
                   [t, r, theta, phi])
        R = g.ricci_tensor(simplify_func=cancel)
        for i in range(4):
            self.assertEqual(cancel(R[i, i]), 0, f"R[{i},{i}] != 0: {R[i,i]}")


class TestHarmonicFunction(unittest.TestCase):
    """Test HarmonicFunction utilities."""

    def test_basic_creation(self):
        hf = HarmonicFunction(transverse_dim=6)
        self.assertEqual(hf.transverse_dim, 6)
        self.assertEqual(len(hf.transverse_coords), 6)

    def test_random_values_harmonic_condition(self):
        """Check that random values satisfy H'' = -(n-1)/r * H'."""
        hf = HarmonicFunction(transverse_dim=6)
        for _ in range(10):
            vals = hf.random_values()
            hpp = vals[hf.Hpp]
            hp = vals[hf.Hp]
            r = vals[hf.r]
            expected = -5.0 / r * hp  # (d_perp - 1) = 5
            self.assertAlmostEqual(hpp, expected, places=12)

    def test_substitute_basic(self):
        """Test that substitution works on a simple expression."""
        hf = HarmonicFunction(transverse_dim=3)
        # After substitution, Hpp should become -(3-1)/r * Hp = -2/r * Hp
        expr = hf.Hpp
        result = hf.substitute(expr)
        self.assertEqual(result, -2 / hf.r * hf.Hp)


class TestWarpedProduct(unittest.TestCase):
    """Test warped product metric construction."""

    def test_d3_brane_metric_shape(self):
        coords = symbols('t x1 x2 x3 y1 y2 y3 y4 y5 y6')
        H = symbols('H', positive=True)
        g = warped_product(
            [Rational(-1, 2), Rational(1, 2)],
            [4, 6],
            ['lorentzian', 'euclidean'],
            coords, H
        )
        self.assertEqual(g.dim, 10)
        self.assertTrue(g.is_diagonal)
        # Time component should be -H^{-1/2}
        self.assertEqual(g.matrix[0, 0], -H**Rational(-1, 2))
        # Spatial worldvolume
        self.assertEqual(g.matrix[1, 1], H**Rational(-1, 2))
        # Transverse
        self.assertEqual(g.matrix[4, 4], H**Rational(1, 2))


# =========================================================================
# Forms tests
# =========================================================================

class TestEpsilon(unittest.TestCase):
    """Test Levi-Civita symbol."""

    def test_identity_permutation(self):
        self.assertEqual(epsilon(3, (0, 1, 2)), 1)

    def test_odd_permutation(self):
        self.assertEqual(epsilon(3, (0, 2, 1)), -1)

    def test_repeated_index(self):
        self.assertEqual(epsilon(3, (0, 0, 2)), 0)

    def test_4d(self):
        self.assertEqual(epsilon(4, (0, 1, 2, 3)), 1)
        self.assertEqual(epsilon(4, (1, 0, 2, 3)), -1)


class TestFormField(unittest.TestCase):
    """Test FormField sparse storage and antisymmetry."""

    def test_antisymmetry(self):
        F = FormField(2, 4)
        F[0, 1] = sp.Symbol('a')
        self.assertEqual(F[1, 0], -sp.Symbol('a'))
        self.assertEqual(F[0, 0], 0)

    def test_3form(self):
        F = FormField(3, 5)
        F[0, 1, 2] = sp.Symbol('b')
        self.assertEqual(F[0, 2, 1], -sp.Symbol('b'))
        self.assertEqual(F[2, 0, 1], sp.Symbol('b'))
        self.assertEqual(F[2, 1, 0], -sp.Symbol('b'))

    def test_nonzero_components(self):
        F = FormField(2, 3)
        F[0, 1] = 5
        F[1, 2] = 3
        nz = F.nonzero_components
        self.assertEqual(len(nz), 2)
        self.assertIn((0, 1), nz)
        self.assertIn((1, 2), nz)

    def test_to_array_roundtrip(self):
        F = FormField(2, 3)
        F[0, 1] = sp.Symbol('a')
        F[0, 2] = sp.Symbol('b')
        arr = F.to_array()
        F2 = FormField.from_array(arr)
        for idx in F.independent_indices:
            self.assertEqual(F[idx], F2[idx])


class TestHodgeStar(unittest.TestCase):
    """Test Hodge dual computation."""

    def test_hodge_2form_in_3d_euclidean(self):
        """In 3D Euclidean, *dx^1 ∧ dx^2 = dx^3 (up to sign)."""
        coords = symbols('x y z')
        g = Metric([1, 1, 1], coords)
        F = FormField(2, 3)
        F[0, 1] = 1
        star_F = hodge_star(F, g, signature=1)
        # *F should be a 1-form
        self.assertEqual(star_F.rank, 1)
        self.assertEqual(star_F[2,], 1)

    def test_double_hodge_euclidean(self):
        """In D-dim Euclidean, **F = (-1)^{p(D-p)} F."""
        coords = symbols('x y z w')
        g = Metric([1, 1, 1, 1], coords)
        F = FormField(2, 4)
        F[0, 1] = sp.Symbol('a')
        F[2, 3] = sp.Symbol('b')

        star_F = hodge_star(F, g, signature=1)
        star_star_F = hodge_star(star_F, g, signature=1)

        # p=2, D=4: (-1)^{2*2} = 1, so **F = F
        for idx in F.independent_indices:
            self.assertEqual(simplify(star_star_F[idx] - F[idx]), 0)


class TestExteriorDerivative(unittest.TestCase):
    """Test exterior derivative computation."""

    def test_d_of_0form(self):
        """d(f) for a scalar should give the gradient as a 1-form."""
        x, y = symbols('x y')
        f = FormField(0, 2, {(): x**2 + y})
        df = exterior_derivative(f, [x, y])
        self.assertEqual(df.rank, 1)
        self.assertEqual(df[0,], 2*x)
        self.assertEqual(df[1,], 1)

    def test_d_squared_is_zero(self):
        """d(d(omega)) = 0 for any 1-form."""
        x, y, z = symbols('x y z')
        omega = FormField(1, 3)
        omega[0,] = x*y
        omega[1,] = y*z
        omega[2,] = x*z

        d_omega = exterior_derivative(omega, [x, y, z])
        dd_omega = exterior_derivative(d_omega, [x, y, z])
        # dd_omega should be zero (rank 3 in 3D has one component)
        for idx in dd_omega.independent_indices:
            self.assertEqual(dd_omega[idx], 0)


class TestFormStressEnergy(unittest.TestCase):
    """Test form norm and stress-energy computation."""

    def test_norm_of_1form_euclidean(self):
        """A 1-form (1, 0, 0) in flat 3D Euclidean should have |F|^2 = 1."""
        g = Metric([1, 1, 1], symbols('x y z'))
        F = FormField(1, 3)
        F[0,] = 1
        self.assertEqual(form_norm_squared(F, g), 1)

    def test_norm_of_2form(self):
        """A 2-form F_{01} = a in flat 4D Minkowski."""
        g = Metric([-1, 1, 1, 1], symbols('t x y z'))
        F = FormField(2, 4)
        a = sp.Symbol('a')
        F[0, 1] = a
        # |F|^2 = (1/2!) * 2 * a^2 * g^{00} g^{11} = a^2 * (-1)(1) = -a^2
        self.assertEqual(form_norm_squared(F, g), -a**2)


# =========================================================================
# Verify tests
# =========================================================================

class TestVerifySymbolic(unittest.TestCase):
    """Test symbolic verification."""

    def test_zero(self):
        is_zero, _ = verify_symbolic(sp.S(0))
        self.assertTrue(is_zero)

    def test_cancellation(self):
        x = sp.Symbol('x')
        expr = (x**2 - 1) / (x - 1) - (x + 1)
        is_zero, _ = verify_symbolic(expr)
        self.assertTrue(is_zero)

    def test_nonzero(self):
        x = sp.Symbol('x')
        is_zero, _ = verify_symbolic(x + 1)
        self.assertFalse(is_zero)


class TestVerifyNumerical(unittest.TestCase):
    """Test numerical verification."""

    def test_zero_expression(self):
        x = sp.Symbol('x')
        expr = (x + 1)**2 - x**2 - 2*x - 1
        is_zero, max_res, _ = verify_numerical(
            expr,
            lambda: {x: __import__('random').uniform(0.5, 3.0)},
        )
        self.assertTrue(is_zero)
        self.assertLess(max_res, 1e-10)

    def test_nonzero_expression(self):
        x = sp.Symbol('x')
        expr = x + 1
        is_zero, _, _ = verify_numerical(
            expr,
            lambda: {x: __import__('random').uniform(0.5, 3.0)},
        )
        self.assertFalse(is_zero)


class TestSolutionFlat(unittest.TestCase):
    """Test Solution verification on flat space (trivial case)."""

    def test_flat_space_no_fields(self):
        coords = symbols('t x y z')
        g = Metric([-1, 1, 1, 1], coords)
        sol = Solution(metric=g)
        report = sol.verify_all()
        self.assertTrue(report.all_passed)


# =========================================================================
# Integration test: Schwarzschild as a vacuum solution
# =========================================================================

class TestSchwarzschildSolution(unittest.TestCase):
    """Schwarzschild should satisfy vacuum Einstein equation (R_MN = 0)."""

    def test_4d(self):
        t, r, theta, phi = symbols('t r theta phi')
        M = symbols('M', positive=True)
        f = 1 - 2*M/r
        g = Metric([-f, 1/f, r**2, r**2 * sp.sin(theta)**2],
                   [t, r, theta, phi])
        sol = Solution(metric=g)
        results = sol.check_einstein(simplify_func=cancel)
        for res in results:
            self.assertTrue(res.passed, f"{res.label} failed: {res.residual}")


if __name__ == '__main__':
    unittest.main()
