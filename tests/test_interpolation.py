import unittest
import numpy as np
from wagl.interpolation import bilinear, subdivide, interpolate_block, interpolate_grid


class TestBilinearFnc(unittest.TestCase):
    def test_bilinear(self):
        """
        Simple test case for the bilinear interpolation function
        """
        expected = np.arange(0, 100, dtype=np.float64).reshape(10, 10)

        in_arr = np.zeros(100).reshape(10, 10)

        result = bilinear(in_arr.shape, 0, 9, 99, 90)
        self.assertTrue(np.allclose(result, expected))


class TestSubdivideFnc(unittest.TestCase):
    def test_subdivide(self):
        """
        Simple test case for the subdivide function
        """
        expected = {
            "UL": [(0, 0), (0, 4), (4, 0), (4, 4)],
            "LL": [(4, 0), (4, 4), (7, 0), (7, 4)],
            "UR": [(0, 4), (0, 7), (4, 4), (4, 7)],
            "LR": [(4, 4), (4, 7), (7, 4), (7, 7)],
        }
        in_arr = np.zeros(64).reshape(8, 8)
        result = subdivide((0, 0), in_arr.shape)
        self.assertEqual(result, expected)


class TestInterpolateBlock(unittest.TestCase):
    def test_interpolate_block(self):
        """
        Simple test case for the interpolate block function
        """
        test_fnc = lambda y, x: 8 * y + x

        expected = np.arange(0, 16).reshape(4, 4) * 2

        in_arr = np.zeros(16).reshape(4, 4)

        in_arr[0, 3] = 3
        in_arr[3, 0] = 12
        in_arr[3, 3] = 15

        result = interpolate_block(
            (0, 0), shape=in_arr.shape, eval_func=test_fnc, grid=in_arr
        )


class TestInterpolateGrid(unittest.TestCase):
    def test_interpolate_grid(self):
        """
        Simple test case for interpolate grid
        """
        result = np.arange(64).reshape(8, 8)
        in_arr = result.copy()  # Interpolation is performed in place
        eval_func = lambda y, x: 8 * y + x
        depth = 3
        interpolate_grid(result, eval_func, depth)
        self.assertTrue(np.allclose(result, in_arr))

    def test_small_grid(self):
        """
        Test grid too small to calculate bilinear interpolation
        """
        with self.assertRaises(ValueError):
            in_arr = np.zeros(1).reshape(1, 1)
            eval_func = lambda y, x: y * 1 + x
            depth = 7
            interpolate_grid(in_arr, eval_func, depth)

    def test_interpolate_default_max_depth(self):
        """
        Test that the wrapper defaults to max depth
        """
        result = np.arange(64).reshape(8, 8)
        in_arr = result.copy()  # Interpolation is performed in place
        eval_func = lambda y, x: 8 * y + x
        depth = 10
        interpolate_grid(result, eval_func, depth)
        self.assertTrue(np.allclose(result, in_arr))
