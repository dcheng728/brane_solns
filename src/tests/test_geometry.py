"""
Tests for sugra.geometry: Metric, HarmonicFunction, warped_product.

The brane Ricci tests compare computed Ricci tensors against the analytical
formulas from hep-th/9803116, equations (59)-(63):

  ds^2_D = H^a ds^2_d + H^b ds^2_{D-d},   H harmonic in (D-d) dims.

  R_uv = -eta_uv * a*(s-2)/4 * (H')^2 * H^{a-b-2}

  R_mn = (H')^2/(4H^2) * [-b(s-2) delta_mn
         - y_m y_n/r^2 (ad(a-b) - s(b+2))]
       + s/(2rH) * H' * [-delta_mn + y_m y_n/r^2 (2+d_tilde)]

where s = ad + b*d_tilde, d_tilde = D - d - 2.
"""

import unittest
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import sympy as sp
from sympy import Rational, cancel, symbols

from sugra.geometry import Metric, HarmonicFunction, warped_product


# =========================================================================
# Metric basics
# =========================================================================

class TestMetric(unittest.TestCase):

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
        self.assertEqual(g.inv_matrix[0, 0], 1 / a)
        self.assertEqual(g.inv_matrix[1, 1], 1 / b)

    def test_det_diagonal(self):
        a, b, c = symbols('a b c', positive=True)
        g = Metric([a, b, c], symbols('x y z'))
        self.assertEqual(g.det, a * b * c)


# =========================================================================
# Flat-space and vacuum Ricci
# =========================================================================

class TestRicciFlat(unittest.TestCase):

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
        g = Metric([-1] + [1] * 9, coords)
        R = g.ricci_tensor()
        for i in range(10):
            self.assertEqual(R[i, i], 0)

    def test_schwarzschild_4d(self):
        t, r, theta, phi = symbols('t r theta phi')
        M = symbols('M', positive=True)
        f = 1 - 2 * M / r
        g = Metric(
            [-f, 1 / f, r**2, r**2 * sp.sin(theta)**2],
            [t, r, theta, phi],
        )
        R = g.ricci_tensor(simplify_func=cancel)
        for i in range(4):
            self.assertEqual(cancel(R[i, i]), 0,
                             f"R[{i},{i}] != 0: {R[i, i]}")


# =========================================================================
# HarmonicFunction
# =========================================================================

class TestHarmonicFunction(unittest.TestCase):

    def test_basic_creation(self):
        y = list(sp.symbols('y0:6', real=True))
        hf = HarmonicFunction(transverse_coords=y)
        self.assertEqual(hf.transverse_dim, 6)
        self.assertEqual(len(hf.transverse_coords), 6)

    def test_random_values_harmonic_condition(self):
        y = list(sp.symbols('y0:6', real=True))
        hf = HarmonicFunction(transverse_coords=y)
        for _ in range(10):
            vals = hf.random_values()
            hpp = vals[hf.Hpp]
            hp = vals[hf.Hp]
            r = vals[hf.r]
            expected = -5.0 / r * hp
            self.assertAlmostEqual(hpp, expected, places=12)

    def test_substitute_basic(self):
        y = list(sp.symbols('y0:3', real=True))
        hf = HarmonicFunction(transverse_coords=y)
        expr = hf.Hpp
        result = hf.substitute(expr)
        self.assertEqual(result, -2 / hf.r * hf.Hp)


# =========================================================================
# Warped product construction
# =========================================================================

class TestWarpedProduct(unittest.TestCase):

    def test_d3_brane_metric_shape(self):
        coords = symbols('t x1 x2 x3 y1 y2 y3 y4 y5 y6')
        H = symbols('H', positive=True)
        g = warped_product(
            [H**Rational(-1, 2), H**Rational(1, 2)],
            [4, 6],
            ['lorentzian', 'euclidean'],
            coords,
        )
        self.assertEqual(g.dim, 10)
        self.assertTrue(g.is_diagonal)
        self.assertEqual(g.matrix[0, 0], -H**Rational(-1, 2))
        self.assertEqual(g.matrix[1, 1], H**Rational(-1, 2))
        self.assertEqual(g.matrix[4, 4], H**Rational(1, 2))


# =========================================================================
# Brane Ricci tensor vs analytical formulas (hep-th/9803116)
# =========================================================================

def _expected_ricci_wv(a, b, d, dtilde, eta_ii, H, Hp):
    """Worldvolume R_{uu} from eq (61), after harmonic condition.

    R_uv = -eta_uv * a*(s-2)/4 * (H')^2 * H^{a-b-2}
    """
    s = a * d + b * dtilde
    return -eta_ii * a * (s - 2) / 4 * Hp**2 * H**(a - b - 2)


def _expected_ricci_tr(a, b, d, dtilde, yk, H, Hp, r):
    """Transverse R_{yy} from eq (61), after harmonic condition.

    R_mn = (H')^2/(4H^2) * [-b(s-2) delta_mn
           - y_m y_n/r^2 (ad(a-b) - s(b+2))]
         + s/(2rH) * H' * [-delta_mn + y_m y_n/r^2 (2+d_tilde)]
    """
    s = a * d + b * dtilde
    result = Hp**2 / (4 * H**2) * (
        -b * (s - 2) - yk**2 / r**2 * (a * d * (a - b) - s * (b + 2))
    )
    result += s / 2 * Hp / (r * H) * (
        -1 + yk**2 / r**2 * (2 + dtilde)
    )
    return result


# Brane configurations: (D_total, d_worldvolume, a, b)
BRANE_CONFIGS = {
    # 11d M-theory
    'M2': (11, 3, Rational(-2, 3), Rational(1, 3)),
    'M5': (11, 6, Rational(-1, 3), Rational(2, 3)),
    # 10d type II (Einstein frame)
    'D1': (10, 2, Rational(-3, 4), Rational(1, 4)),
    'D2': (10, 3, Rational(-5, 8), Rational(3, 8)),
    'D3': (10, 4, Rational(-1, 2), Rational(1, 2)),
    'D4': (10, 5, Rational(-3, 8), Rational(5, 8)),
    'D5': (10, 6, Rational(-1, 4), Rational(3, 4)),
}


class TestBraneRicci(unittest.TestCase):
    """Compare computed Ricci tensor against analytical brane formulas."""

    def _check_brane(self, name, D_tot, d, a, b):
        D_perp = D_tot - d
        dtilde = D_tot - d - 2

        # Verify constraint b = -ad/d_tilde
        self.assertEqual(b, -a * d / dtilde,
                         f"{name}: b = -ad/d_tilde not satisfied")

        # Coordinates
        wv_names = ['t'] + [f'x{i}' for i in range(1, d)]
        wv_coords = list(sp.symbols(' '.join(wv_names), real=True))

        y = list(sp.symbols(f'y0:{D_perp}', real=True))
        coords = wv_coords + y
        hf = HarmonicFunction(transverse_coords=y)

        H_func = sp.Function('H')(hf.r_expr)

        metric = warped_product(
            warp_factors=[H_func**a, H_func**b],
            block_dims=[d, D_perp],
            block_signatures=['lorentzian', 'euclidean'],
            coordinates=coords,
        )

        R = metric.ricci_tensor(simplify_func=cancel)

        H = hf.H
        Hp = hf.Hp
        r = hf.r

        # -- Worldvolume time: eta_tt = -1 --
        computed = hf.substitute(cancel(R[0, 0]))
        expected = _expected_ricci_wv(a, b, d, dtilde, -1, H, Hp)
        diff = cancel(computed - expected)
        self.assertEqual(diff, 0,
                         f"{name} R_tt: computed={computed}, expected={expected}")

        # -- Worldvolume space: eta_{x1 x1} = +1 --
        computed = hf.substitute(cancel(R[1, 1]))
        expected = _expected_ricci_wv(a, b, d, dtilde, +1, H, Hp)
        diff = cancel(computed - expected)
        self.assertEqual(diff, 0,
                         f"{name} R_x1x1: computed={computed}, expected={expected}")

        # -- Transverse y0 --
        y0 = y[0]
        idx0 = d
        computed = hf.substitute(cancel(R[idx0, idx0]))
        expected = _expected_ricci_tr(a, b, d, dtilde, y0, H, Hp, r)
        diff = cancel(computed - expected)
        self.assertEqual(diff, 0,
                         f"{name} R_y0y0: computed={computed}, expected={expected}")

        # -- Transverse y1 --
        y1 = y[1]
        idx1 = d + 1
        computed = hf.substitute(cancel(R[idx1, idx1]))
        expected = _expected_ricci_tr(a, b, d, dtilde, y1, H, Hp, r)
        diff = cancel(computed - expected)
        self.assertEqual(diff, 0,
                         f"{name} R_y1y1: computed={computed}, expected={expected}")

    # -- 11d branes --

    def test_M2(self):
        self._check_brane('M2', *BRANE_CONFIGS['M2'])

    def test_M5(self):
        self._check_brane('M5', *BRANE_CONFIGS['M5'])

    # -- 10d branes --

    def test_D1(self):
        self._check_brane('D1', *BRANE_CONFIGS['D1'])

    def test_D2(self):
        self._check_brane('D2', *BRANE_CONFIGS['D2'])

    def test_D3(self):
        self._check_brane('D3', *BRANE_CONFIGS['D3'])

    def test_D4(self):
        self._check_brane('D4', *BRANE_CONFIGS['D4'])

    def test_D5(self):
        self._check_brane('D5', *BRANE_CONFIGS['D5'])


if __name__ == '__main__':
    unittest.main()
