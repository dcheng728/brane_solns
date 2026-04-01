"""Tests for tensorGlow.core.metric."""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..'))

import unittest
from tensorGlow import (
    IndexType, Index, MetricTensor, TensorHead, TensorSymmetry,
    PartialDerivative,
)
from tensorGlow.core.expr import Prod, Sum, Scalar


L = IndexType('V', dim=4, dummy_prefix='x')
g = MetricTensor(L, name='g')
T = TensorHead('T', [L, L], TensorSymmetry.no_symmetry(2))
R = TensorHead('R', [L, L, L, L], TensorSymmetry.no_symmetry(4))
V = TensorHead('V', [L], TensorSymmetry.no_symmetry(1))

a = Index('a', L, is_up=True)
b = Index('b', L, is_up=True)
c = Index('c', L, is_up=True)
d = Index('d', L, is_up=True)
e = Index('e', L, is_up=True)



class TestMetricContraction(unittest.TestCase):
    """contract_metrics should raise/lower indices correctly."""

    # ── Basic raising and lowering ───────────────────────────────
    def test_identity_contraction1(self):
        """
        $$g_{ab}g^{bc} = \delta_a^c$$
        """
        result = g.contract_metrics(g(-a, -b) * g(b, c))
        fi = result.free_indices
        self.assertEqual(len(fi), 2)
        by_name = {idx.name: idx.is_up for idx in fi}
        self.assertFalse(by_name['a'])
        self.assertTrue(by_name['c'])


    def test_lower_vector(self):
        """g_{ab} V^a = V_b: single free index should be b_down."""
        result = g.contract_metrics(g(-a, -b) * V(a))
        fi = result.free_indices
        self.assertEqual(len(fi), 1)
        self.assertEqual(fi[0].name, 'b')
        self.assertFalse(fi[0].is_up)

    def test_lower_one_index_of_rank2(self):
        """g_{ab} T^{ac} = T_b^c: b should be down, c stays up."""
        result = g.contract_metrics(g(-a, -b) * T(a, c))
        fi = result.free_indices
        self.assertEqual(len(fi), 2)
        by_name = {idx.name: idx.is_up for idx in fi}
        self.assertFalse(by_name['b'])
        self.assertTrue(by_name['c'])

    def test_raise_vector(self):
        """g^{ab} V_a = V^b: single free index should be b_up."""
        result = g.contract_metrics(g.inv(a, b) * V(-a))
        fi = result.free_indices
        self.assertEqual(len(fi), 1)
        self.assertEqual(fi[0].name, 'b')
        self.assertTrue(fi[0].is_up)

    def test_raise_one_index_of_rank2(self):
        """g^{ab} T_{ac} = T^b_c: b should be up, c stays down."""
        result = g.contract_metrics(g.inv(a, b) * T(-a, -c))
        fi = result.free_indices
        self.assertEqual(len(fi), 2)
        by_name = {idx.name: idx.is_up for idx in fi}
        self.assertTrue(by_name['b'])
        self.assertFalse(by_name['c'])

    # ── Metric inverse identity ──────────────────────────────────

    def test_metric_inverse(self):
        """g_{ab} g^{bc}: a should be down, c should be up."""
        result = g.contract_metrics(g(-a, -b) * g.inv(b, c))
        fi = result.free_indices
        by_name = {idx.name: idx.is_up for idx in fi}
        self.assertFalse(by_name['a'])
        self.assertTrue(by_name['c'])

    def test_metric_inverse_other_order(self):
        """g^{ab} g_{bc}: a should be up, c should be down."""
        result = g.contract_metrics(g.inv(a, b) * g(-b, -c))
        fi = result.free_indices
        by_name = {idx.name: idx.is_up for idx in fi}
        self.assertTrue(by_name['a'])
        self.assertFalse(by_name['c'])

    # ── Trace ────────────────────────────────────────────────────

    def test_delta_trace(self):
        """delta^a_a = dim = 4."""
        result = g.delta(a, -a)
        self.assertIsInstance(result, Scalar)
        self.assertEqual(result.expr, 4)

    def test_metric_full_trace(self):
        """g_{ab} g^{ab} = dim = 4 (scalar, no free indices)."""
        result = g.contract_metrics(g(-a, -b) * g.inv(a, b))
        fi = result.free_indices
        self.assertEqual(len(fi), 0)

    # ── No contraction possible ──────────────────────────────────

    def test_no_contraction_unchanged(self):
        """g_{ab} T^{cd}: no shared index, all four indices survive."""
        result = g.contract_metrics(g(-a, -b) * T(c, d))
        self.assertEqual(len(result.free_indices), 4)

    # ── Scalar coefficient preserved ─────────────────────────────

    def test_scalar_coefficient(self):
        """3 * g_{ab} V^a should preserve coefficient 3."""
        result = g.contract_metrics(3 * g(-a, -b) * V(a))
        self.assertIsInstance(result, Prod)
        self.assertEqual(result.coeff, 3)

    # ── Sum distribution ─────────────────────────────────────────

    def test_sum_contracts_each_term(self):
        """g_{ab} (T^{ac} + T^{ac}): both terms should contract."""
        expr = g(-a, -b) * (T(a, c) + T(a, c))
        result = g.contract_metrics(expr)
        self.assertIsInstance(result, Sum)
        for term in result.terms:
            by_name = {idx.name: idx.is_up for idx in term.free_indices}
            self.assertFalse(by_name['b'])
            self.assertTrue(by_name['c'])

    # ── Non-Tensor factors (Apply) must survive ──────────────────

    def test_apply_factor_preserved(self):
        """contract_metrics must not drop Apply (partial derivative) factors.

        # g^{ab} * \partial_c(g_{de})
        : no contraction possible,
        so both the metric and the Apply should survive (5 free indices).
        """
        partial = PartialDerivative(L, 'partial')
        deriv = partial(-c, g(-d, -e))
        result = g.contract_metrics(g.inv(a, b) * deriv)
        self.assertEqual(len(result.free_indices), 5)

    # ── raise_index / lower_index roundtrip ──────────────────────

    def test_raise_then_contract(self):
        """raise_index introduces g^{ab}; contract_metrics should absorb it."""
        raised = g.raise_index(V(-a), -a)
        result = g.contract_metrics(raised)
        fi = result.free_indices
        self.assertEqual(len(fi), 1)
        self.assertTrue(fi[0].is_up)

    def test_lower_then_contract(self):
        """lower_index introduces g_{ab}; contract_metrics should absorb it."""
        lowered = g.lower_index(V(a), a)
        result = g.contract_metrics(lowered)
        fi = result.free_indices
        self.assertEqual(len(fi), 1)
        self.assertFalse(fi[0].is_up)


if __name__ == '__main__':
    unittest.main()
