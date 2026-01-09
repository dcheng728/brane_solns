import unittest
import time
from helpers.form_helpers import epsilon, perm_with_sign

class TestPermWithSign(unittest.TestCase):
    def test_perm_with_sign_3_elements(self):
        elements = (0, 1, 2)
        expected_outputs = {
            (0, 1, 2): 1,
            (0, 2, 1): -1,
            (1, 0, 2): -1,
            (1, 2, 0): 1,
            (2, 0, 1): 1,
            (2, 1, 0): -1,
        }
        for tup, sgn in perm_with_sign(elements):
            self.assertEqual(expected_outputs[tup], sgn)

    def test_perm_with_sign_4_elements(self):
        elements = (0, 1, 2, 3)
        expected_outputs = {
            (0, 1, 2, 3): 1,
            (0, 1, 3, 2): -1,
            (0, 2, 1, 3): -1,
            (0, 2, 3, 1): 1,
            (0, 3, 1, 2): 1,
            (0, 3, 2, 1): -1,
            (1, 0, 2, 3): -1,
            (1, 0, 3, 2): 1,
            (1, 2, 0, 3): 1,
            (1, 2, 3, 0): -1,
            (1, 3, 0, 2): -1,
            (1, 3, 2, 0): 1,
            (2, 0, 1, 3): 1,
            (2, 0, 3, 1): -1,
            (2, 1, 0, 3): -1,
            (2, 1, 3, 0): 1,
            (2, 3, 0, 1): 1,
            (2, 3, 1, 0): -1,
            (3, 0, 1, 2): -1,
            (3, 0, 2, 1): 1,
            (3, 1, 0, 2): 1,
            (3, 1, 2, 0): -1,
            (3, 2, 0, 1): -1,
            (3, 2, 1, 0): 1,
        }
        for tup, sgn in perm_with_sign(elements):
            self.assertEqual(expected_outputs[tup], sgn)

class TestEpsilon(unittest.TestCase):
    def test_epsilon_3d(self):
        # Test the epsilon function in 3 dimensions
        num_dims = 3
        expected = {
            (0, 1, 2): 1,
            (0, 2, 1): -1,
            (1, 0, 2): -1,
            (1, 2, 0): 1,
            (2, 0, 1): 1,
            (2, 1, 0): -1,
        }
        start_time = time.time()
        for indices, value in expected.items():
            self.assertEqual(epsilon(num_dims, indices), value)
        end_time = time.time()
        print(f"3D epsilon test completed in {(end_time - start_time)/len(expected):.6f} seconds/call")

    def test_epsilon_4d(self):
        num_dims = 4
        expected = {
            (0, 1, 2, 3): 1,
            (0, 1, 3, 2): -1,
            (0, 2, 1, 3): -1,
            (0, 2, 3, 1): 1,
            (0, 3, 1, 2): 1,
            (0, 3, 2, 1): -1,
            (1, 0, 2, 3): -1,
            (1, 0, 3, 2): 1,
            (1, 2, 0, 3): 1,
            (1, 2, 3, 0): -1,
            (1, 3, 0, 2): -1,
            (1, 3, 2, 0): 1,
            (2, 0, 1, 3): 1,
            (2, 0, 3, 1): -1,
            (2, 1, 0, 3): -1,
            (2, 1, 3, 0): 1,
            (2, 3, 0, 1): 1,
            (2, 3, 1, 0): -1,
            (3, 0, 1, 2): -1,
            (3, 0, 2, 1): 1,
            (3, 1, 0, 2): 1,
            (3, 1, 2, 0): -1,
            (3, 2, 0, 1): -1,
            (3, 2, 1, 0): 1,
        }
        start_time = time.time()
        for indices, value in expected.items():
            self.assertEqual(epsilon(num_dims, indices), value)
        end_time = time.time()
        print(f"4D epsilon test completed in {(end_time - start_time)/len(expected):.6f} seconds/call")


if __name__ == '__main__':
    unittest.main()