"""Tests for tensorGlow core algebra (Phase 1)."""

import unittest
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import sympy as sp
from tensorGlow import (
    IndexType, Index, indices, TensorHead, TensorSymmetry,
    TensorProduct, TensorSum, rename_dummies, to_latex,
)


class TestIndex(unittest.TestCase):
    def setUp(self):
        self.L = IndexType('Lorentz', dim=4, dummy_prefix='L')

    def test_create_index(self):
        a = Index('a', self.L, is_up=True)
        self.assertEqual(a.name, 'a')
        self.assertTrue(a.is_up)
        self.assertTrue(a.is_contravariant)

    def test_negate_flips_variance(self):
        a = Index('a', self.L, is_up=True)
        a_down = -a
        self.assertEqual(a_down.name, 'a')
        self.assertFalse(a_down.is_up)
        self.assertTrue(a_down.is_covariant)

    def test_matches(self):
        a_up = Index('a', self.L, is_up=True)
        a_down = Index('a', self.L, is_up=False)
        b_up = Index('b', self.L, is_up=True)
        self.assertTrue(a_up.matches(a_down))
        self.assertTrue(a_down.matches(a_up))
        self.assertFalse(a_up.matches(b_up))
        self.assertFalse(a_up.matches(a_up))  # same variance

    def test_indices_factory(self):
        a, b, c = indices('a b c', self.L)
        self.assertEqual(a.name, 'a')
        self.assertEqual(b.name, 'b')
        self.assertEqual(c.name, 'c')
        self.assertIs(a.index_type, self.L)

    def test_single_index_factory(self):
        a = indices('a', self.L)
        self.assertIsInstance(a, Index)
        self.assertEqual(a.name, 'a')


class TestTensorHead(unittest.TestCase):
    def setUp(self):
        self.L = IndexType('Lorentz', dim=4)
        self.a, self.b, self.c, self.d = indices('a b c d', self.L)

    def test_create_vector(self):
        V = TensorHead('V', [self.L])
        self.assertEqual(V.rank, 1)

    def test_call_creates_product(self):
        V = TensorHead('V', [self.L])
        expr = V(self.a)
        self.assertIsInstance(expr, TensorProduct)
        self.assertEqual(len(expr.atoms), 1)
        self.assertEqual(expr.atoms[0].head, V)

    def test_wrong_number_of_indices(self):
        V = TensorHead('V', [self.L])
        with self.assertRaises(ValueError):
            V(self.a, self.b)

    def test_wrong_index_type(self):
        L2 = IndexType('Other', dim=3)
        V = TensorHead('V', [self.L])
        i = Index('i', L2)
        with self.assertRaises(TypeError):
            V(i)


class TestExprAlgebra(unittest.TestCase):
    def setUp(self):
        self.L = IndexType('Lorentz', dim=4)
        self.a, self.b, self.c, self.d = indices('a b c d', self.L)
        self.V = TensorHead('V', [self.L])
        self.A = TensorHead('A', [self.L, self.L])

    def test_product_free_indices(self):
        expr = self.V(self.a)
        self.assertEqual(len(expr.free_indices), 1)
        self.assertEqual(expr.free_indices[0].name, 'a')

    def test_contraction_dummy_pairs(self):
        # V^a * V_a (trace)
        expr = self.V(self.a) * self.V(-self.a)
        self.assertEqual(len(expr.dummy_pairs), 1)
        self.assertEqual(len(expr.free_indices), 0)

    def test_mixed_contraction(self):
        # A^{ab} * V_b — one dummy, one free
        expr = self.A(self.a, self.b) * self.V(-self.b)
        self.assertEqual(len(expr.dummy_pairs), 1)
        self.assertEqual(len(expr.free_indices), 1)
        self.assertEqual(expr.free_indices[0].name, 'a')

    def test_scalar_multiplication(self):
        expr = self.V(self.a)
        doubled = 2 * expr
        self.assertIsInstance(doubled, TensorProduct)
        self.assertEqual(doubled.coeff, 2)

    def test_scalar_multiplication_right(self):
        expr = self.V(self.a) * sp.Rational(1, 2)
        self.assertEqual(expr.coeff, sp.Rational(1, 2))

    def test_addition(self):
        expr = self.V(self.a) + self.V(self.a)
        self.assertIsInstance(expr, TensorSum)
        self.assertEqual(len(expr.terms), 2)

    def test_subtraction(self):
        expr = self.V(self.a) - self.V(self.a)
        self.assertIsInstance(expr, TensorSum)

    def test_negation(self):
        expr = -self.V(self.a)
        self.assertIsInstance(expr, TensorProduct)
        self.assertEqual(expr.coeff, -1)

    def test_repr(self):
        expr = self.V(self.a)
        self.assertIn('V', repr(expr))


class TestRenameDummies(unittest.TestCase):
    def setUp(self):
        self.L = IndexType('Lorentz', dim=4, dummy_prefix='L')
        self.a, self.b = indices('a b', self.L)
        self.V = TensorHead('V', [self.L])

    def test_rename(self):
        expr = self.V(self.a) * self.V(-self.a)
        renamed = rename_dummies(expr)
        # The dummy should be renamed to L_0
        all_names = [idx.name for idx in renamed.all_indices]
        self.assertIn('L_0', all_names)

    def test_no_dummies_unchanged(self):
        expr = self.V(self.a)
        renamed = rename_dummies(expr)
        self.assertEqual(renamed.all_indices[0].name, 'a')


class TestLaTeX(unittest.TestCase):
    def setUp(self):
        self.L = IndexType('Lorentz', dim=4)
        self.a, self.b, self.c, self.d = indices('a b c d', self.L)

    def test_vector_latex(self):
        V = TensorHead('V', [self.L])
        expr = V(self.a)
        latex = to_latex(expr)
        self.assertIn('V', latex)
        self.assertIn('a', latex)

    def test_mixed_tensor_latex(self):
        T = TensorHead('T', [self.L, self.L])
        expr = T(self.a, -self.b)
        latex = to_latex(expr)
        self.assertIn('^{a}', latex)
        self.assertIn('_{b}', latex)

    def test_sum_latex(self):
        V = TensorHead('V', [self.L])
        expr = V(self.a) + 2 * V(self.a)
        latex = to_latex(expr)
        self.assertIn('+', latex)

    def test_zero_sum_latex(self):
        latex = to_latex(TensorSum(()))
        self.assertEqual(latex, '0')


class TestSymmetry(unittest.TestCase):
    def test_no_symmetry(self):
        sym = TensorSymmetry.no_symmetry(2)
        self.assertEqual(sym.rank, 2)

    def test_fully_symmetric(self):
        sym = TensorSymmetry.fully_symmetric(3)
        self.assertEqual(sym.rank, 3)
        self.assertTrue(len(sym.generators) > 0)

    def test_fully_antisymmetric(self):
        sym = TensorSymmetry.fully_antisymmetric(3)
        self.assertEqual(sym.rank, 3)

    def test_riemann(self):
        sym = TensorSymmetry.riemann()
        self.assertEqual(sym.rank, 4)


class TestCanonicalization(unittest.TestCase):
    """Test the Butler-Portugal canonicalization bridge."""

    def setUp(self):
        self.L = IndexType('Lorentz', dim=4, dummy_prefix='L')
        self.a, self.b, self.c, self.d = indices('a b c d', self.L)
        self.e, self.f = indices('e f', self.L)

    def test_antisymmetric_swap(self):
        """A_{ba} = -A_{ab}."""
        A = TensorHead('A', [self.L, self.L], TensorSymmetry.fully_antisymmetric(2))
        expr = A(-self.b, -self.a)
        canon = expr.canon_bp()
        # Should be -A_{ab}
        self.assertEqual(canon.coeff, -1)
        idx = canon.atoms[0].indices
        self.assertEqual(idx[0].name, 'a')
        self.assertEqual(idx[1].name, 'b')

    def test_symmetric_swap(self):
        """g_{ba} = g_{ab} (symmetric tensor is already canonical)."""
        g = TensorHead('g', [self.L, self.L], TensorSymmetry.fully_symmetric(2))
        expr = g(-self.b, -self.a)
        canon = expr.canon_bp()
        self.assertEqual(canon.coeff, 1)
        idx = canon.atoms[0].indices
        self.assertEqual(idx[0].name, 'a')
        self.assertEqual(idx[1].name, 'b')

    def test_riemann_antisymmetry_last_pair(self):
        """R_{abdc} = -R_{abcd}."""
        R = TensorHead('R', [self.L]*4, TensorSymmetry.riemann())
        expr = R(-self.a, -self.b, -self.d, -self.c)
        canon = expr.canon_bp()
        self.assertEqual(canon.coeff, -1)
        names = [idx.name for idx in canon.atoms[0].indices]
        self.assertEqual(names, ['a', 'b', 'c', 'd'])

    def test_riemann_pair_symmetry(self):
        """R_{cdab} = R_{abcd} (pair exchange)."""
        R = TensorHead('R', [self.L]*4, TensorSymmetry.riemann())
        expr = R(-self.c, -self.d, -self.a, -self.b)
        canon = expr.canon_bp()
        self.assertEqual(canon.coeff, 1)
        names = [idx.name for idx in canon.atoms[0].indices]
        self.assertEqual(names, ['a', 'b', 'c', 'd'])

    def test_antisym_sum_vanishes(self):
        """A_{ab} + A_{ba} = 0."""
        A = TensorHead('A', [self.L, self.L], TensorSymmetry.fully_antisymmetric(2))
        expr = A(-self.a, -self.b) + A(-self.b, -self.a)
        canon = expr.canon_bp()
        # Should simplify to 0 after collecting
        total = sum(t.coeff for t in canon.terms) if isinstance(canon, TensorSum) else canon.coeff
        self.assertEqual(total, 0)

    def test_riemann_antisym_sum_vanishes(self):
        """R_{abcd} + R_{abdc} = 0."""
        R = TensorHead('R', [self.L]*4, TensorSymmetry.riemann())
        expr = R(-self.a, -self.b, -self.c, -self.d) + R(-self.a, -self.b, -self.d, -self.c)
        canon = expr.canon_bp()
        total = sum(t.coeff for t in canon.terms) if isinstance(canon, TensorSum) else canon.coeff
        self.assertEqual(total, 0)


if __name__ == '__main__':
    unittest.main()
