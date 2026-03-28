"""
Tests for sugra.verifier: Verify R_{MN} = T_{MN} for standard brane solutions.

Each test case explicitly constructs the metric, form field, and dilaton
for a specific brane, matching the conventions from hep-th/9701088.
"""

import unittest
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import sympy as sp
from sympy import Rational as R

from sugra.geometry import HarmonicFunction, warped_product
from sugra.forms import FormField, exterior_derivative, hodge_star
from sugra.verifier import Verifier


def make_coords(d_wv, d_harm):
    """Set up coordinates and harmonic function for a (d_wv + d_harm)-dim brane."""
    wv_coords = list(sp.symbols(
        ' '.join(['t'] + [f'x{i}' for i in range(1, d_wv)]), real=True
    ))
    harmonic_coords = list(sp.symbols(f'y0:{d_harm}', real=True))
    coords = wv_coords + harmonic_coords
    hf = HarmonicFunction(transverse_coords=harmonic_coords)
    H = sp.Function('H')(hf.r_expr)
    return coords, hf, H


class TestVerifierF1(unittest.TestCase):
    """F1 (fundamental string): H_3 = dB_2, e^Phi = H^{-1/2}, alpha = -1."""

    def test_einstein(self):
        coords, hf, H = make_coords(d_wv=2, d_harm=8)
        D = 10

        # Metric: ds^2 = H^{-3/4} ds^2_2 + H^{1/4} ds^2_8
        metric = warped_product(
            warp_factors=[H**R(-3, 4), H**R(1, 4)],
            block_dims=[2, 8],
            block_signatures=['lorentzian', 'euclidean'],
            coordinates=coords,
        )

        # B_{01} = H^{-1},  H_3 = dB_2
        B2 = FormField(rank=2, dim=D)
        B2[(0, 1)] = 1 / H
        F = exterior_derivative(B2, coords)

        Phi = R(-1, 2) * sp.log(H)
        alpha = R(-1, 1)

        soln = dict(metric=metric, F=F, Phi=Phi, alpha=alpha, coords=coords, hf=hf)
        self.assertTrue(Verifier(soln).check())


class TestVerifierD1(unittest.TestCase):
    """D1-brane: F_3 = dC_2, e^Phi = H^{1/2}, alpha = 1."""

    def test_einstein(self):
        coords, hf, H = make_coords(d_wv=2, d_harm=8)
        D = 10

        metric = warped_product(
            warp_factors=[H**R(-3, 4), H**R(1, 4)],
            block_dims=[2, 8],
            block_signatures=['lorentzian', 'euclidean'],
            coordinates=coords,
        )

        C2 = FormField(rank=2, dim=D)
        C2[(0, 1)] = 1 / H
        F = exterior_derivative(C2, coords)

        Phi = R(1, 2) * sp.log(H)
        alpha = R(1, 1)

        soln = dict(metric=metric, F=F, Phi=Phi, alpha=alpha, coords=coords, hf=hf)
        self.assertTrue(Verifier(soln).check())


class TestVerifierD2(unittest.TestCase):
    """D2-brane: F_4 = dC_3, e^Phi = H^{1/4}, alpha = 1/2."""

    def test_einstein(self):
        coords, hf, H = make_coords(d_wv=3, d_harm=7)
        D = 10

        metric = warped_product(
            warp_factors=[H**R(-5, 8), H**R(3, 8)],
            block_dims=[3, 7],
            block_signatures=['lorentzian', 'euclidean'],
            coordinates=coords,
        )

        C3 = FormField(rank=3, dim=D)
        C3[(0, 1, 2)] = 1 / H
        F = exterior_derivative(C3, coords)

        Phi = R(1, 4) * sp.log(H)
        alpha = R(1, 2)

        soln = dict(metric=metric, F=F, Phi=Phi, alpha=alpha, coords=coords, hf=hf)
        self.assertTrue(Verifier(soln).check())


class TestVerifierD3(unittest.TestCase):
    """D3-brane: self-dual F_5, Phi = 0, alpha = 0."""

    def test_einstein(self):
        coords, hf, H = make_coords(d_wv=4, d_harm=6)
        D = 10

        metric = warped_product(
            warp_factors=[H**R(-1, 2), H**R(1, 2)],
            block_dims=[4, 6],
            block_signatures=['lorentzian', 'euclidean'],
            coordinates=coords,
        )

        C4 = FormField(rank=4, dim=D)
        C4[(0, 1, 2, 3)] = 1 / H
        F_E = exterior_derivative(C4, coords)
        F = (F_E + hodge_star(F_E, metric)) / sp.sqrt(2)

        Phi = sp.S(0)
        alpha = R(0, 1)

        soln = dict(metric=metric, F=F, Phi=Phi, alpha=alpha, coords=coords, hf=hf)
        self.assertTrue(Verifier(soln).check())


class TestVerifierD4(unittest.TestCase):
    """D4-brane: F_6 = dC_5, e^Phi = H^{-1/4}, alpha = -1/2."""

    def test_einstein(self):
        coords, hf, H = make_coords(d_wv=5, d_harm=5)
        D = 10

        metric = warped_product(
            warp_factors=[H**R(-3, 8), H**R(5, 8)],
            block_dims=[5, 5],
            block_signatures=['lorentzian', 'euclidean'],
            coordinates=coords,
        )

        C5 = FormField(rank=5, dim=D)
        C5[(0, 1, 2, 3, 4)] = 1 / H
        F = exterior_derivative(C5, coords)

        Phi = R(-1, 4) * sp.log(H)
        alpha = R(-1, 2)

        soln = dict(metric=metric, F=F, Phi=Phi, alpha=alpha, coords=coords, hf=hf)
        self.assertTrue(Verifier(soln).check())


class TestVerifierD5(unittest.TestCase):
    """D5-brane: magnetic F_3 = *dC_6, e^Phi = H^{-1/2}, alpha = -1."""

    def test_einstein(self):
        coords, hf, H = make_coords(d_wv=6, d_harm=4)
        D = 10

        metric = warped_product(
            warp_factors=[H**R(-1, 4), H**R(3, 4)],
            block_dims=[6, 4],
            block_signatures=['lorentzian', 'euclidean'],
            coordinates=coords,
        )

        # Electric C_6, then Hodge dual to get magnetic F_3
        C6 = FormField(rank=6, dim=D)
        C6[(0, 1, 2, 3, 4, 5)] = 1 / H
        F_E = exterior_derivative(C6, coords)
        F = hodge_star(F_E, metric)

        Phi = R(-1, 2) * sp.log(H)
        alpha = R(-1, 1)

        soln = dict(metric=metric, F=F, Phi=Phi, alpha=alpha, coords=coords, hf=hf)
        self.assertTrue(Verifier(soln).check())


class TestVerifierNS5(unittest.TestCase):
    """NS5-brane: magnetic H_3 = *dB_6, e^Phi = H^{1/2}, alpha = 1."""

    def test_einstein(self):
        coords, hf, H = make_coords(d_wv=6, d_harm=4)
        D = 10

        metric = warped_product(
            warp_factors=[H**R(-1, 4), H**R(3, 4)],
            block_dims=[6, 4],
            block_signatures=['lorentzian', 'euclidean'],
            coordinates=coords,
        )

        # Electric B_6, then Hodge dual to get magnetic H_3
        B6 = FormField(rank=6, dim=D)
        B6[(0, 1, 2, 3, 4, 5)] = 1 / H
        F_E = exterior_derivative(B6, coords)
        F = hodge_star(F_E, metric)

        Phi = R(1, 2) * sp.log(H)
        alpha = R(1, 1)

        soln = dict(metric=metric, F=F, Phi=Phi, alpha=alpha, coords=coords, hf=hf)
        self.assertTrue(Verifier(soln).check())

class TestVerifierM2(unittest.TestCase):
    """M2-brane: F_4 = dC_3, no dilaton, alpha = 0."""

    def test_einstein(self):
        coords, hf, H = make_coords(d_wv=3, d_harm=8)
        D = 11

        # Metric: ds^2 = H^{-3/4} ds^2_2 + H^{1/4} ds^2_8
        metric = warped_product(
            warp_factors=[H**R(-2, 3), H**R(1, 3)],
            block_dims=[3, 8],
            block_signatures=['lorentzian', 'euclidean'],
            coordinates=coords,
        )

        # B_{01} = H^{-1},  H_3 = dB_2
        C3 = FormField(rank=3, dim=D)
        C3[(0, 1, 2)] = 1 / H
        F = exterior_derivative(C3, coords)

        Phi = 0
        alpha = 0

        soln = dict(metric=metric, F=F, Phi=Phi, alpha=alpha, coords=coords, hf=hf)
        self.assertTrue(Verifier(soln).check())

class TestVerifierM5(unittest.TestCase):
    """M5-brane: F_4 = *dC_6, no dilaton, alpha = 0."""

    def test_einstein(self):
        coords, hf, H = make_coords(d_wv=6, d_harm=5)
        D = 11

        # Metric: ds^2 = H^{-3/4} ds^2_2 + H^{1/4} ds^2_8
        metric = warped_product(
            warp_factors=[H**R(-1, 3), H**R(2, 3)],
            block_dims=[6, 5],
            block_signatures=['lorentzian', 'euclidean'],
            coordinates=coords,
        )

        # B_{01} = H^{-1},  H_3 = dB_2
        C6 = FormField(rank=6, dim=D)
        C6[(0, 1, 2, 3, 4, 5)] = 1 / H
        F_E = exterior_derivative(C6, coords)
        F = hodge_star(F_E, metric)

        Phi = 0
        alpha = 0

        soln = dict(metric=metric, F=F, Phi=Phi, alpha=alpha, coords=coords, hf=hf)
        self.assertTrue(Verifier(soln).check())

if __name__ == '__main__':
    unittest.main()
