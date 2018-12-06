#!/usr/bin/env python

"""
Unit test for the bilinear recursive bisection function found in
wagl.interpolation.bilinear
"""

from __future__ import absolute_import
import unittest
import numpy
import math
from wagl.interpolation import bilinear, subdivide, indices, interpolate_block


class BLRBTest(unittest.TestCase):

    def setUp(self):
        self.origin = (0, 0)
        self.shape = (16, 32)  # (nrows, ncols)

    def test_indices(self):
        t = indices(self.origin, self.shape)
        self.assertEqual(t, (0, 15, 0, 31))
        t = indices((2, 3), (3, 4))
        self.assertEqual(t, (2, 4, 3, 6))

    def test_subdivide(self):
        d = subdivide(self.origin, self.shape)
        self.assertEqual(sorted(d.keys()), ['LL', 'LR', 'UL', 'UR'])
        self.assertEqual(d['UL'], [(0, 0), (0, 16), (8, 0), (8, 16)])
        self.assertEqual(d['UR'], [(0, 16), (0, 31), (8, 16), (8, 31)])
        self.assertEqual(d['LL'], [(8, 0), (8, 16), (15, 0), (15, 16)])
        self.assertEqual(d['LR'], [(8, 16), (8, 31), (15, 16), (15, 31)])

    def test_bilinear_0(self):
        x = math.pi
        a = bilinear((5, 5), x, x, x, x)
        self.assertEqual(a[0, 0], x)
        self.assertEqual(a[0, 4], x)
        self.assertEqual(a[4, 4], x)
        self.assertEqual(a[2, 2], x)
        self.assertEqual(a[4, 0], x)

    def test_bilinear_1(self):
        a = bilinear((5, 5), 0.0, 1.0, 1.0, 0.0)
        # print '\n', a
        self.assertEqual(a[0, 0], 0.0)
        self.assertEqual(a[0, 4], 1.0)
        self.assertEqual(a[4, 0], 0.0)
        self.assertEqual(a[2, 2], 0.5)
        self.assertEqual(a[4, 4], 1.0)

    def test_bilinear_2(self):
        a = bilinear((5, 5), 0.0, 1.0, 2.0, 1.0)
        # print '\n', a
        self.assertEqual(a[0, 0], 0.0)
        self.assertEqual(a[0, 4], 1.0)
        self.assertEqual(a[4, 0], 1.0)
        self.assertEqual(a[2, 2], 1.0)
        self.assertEqual(a[4, 4], 2.0)

    def test_interpolate_block_0(self):

        def f(i, j):
            return float(i * j)

        b = interpolate_block((0, 0), (5, 5), f, grid=None)
        # print '\n', b

        self.assertEqual(b[0, 0], 0.0)
        self.assertEqual(b[0, 4], 0.0)
        self.assertEqual(b[4, 0], 0.0)
        self.assertEqual(b[2, 2], 4.0)
        self.assertEqual(b[4, 4], 16.0)

    def test_interpolate_block_1(self):

        def f(i, j):
            return float(i * j)

        b = interpolate_block((0, 0), (5, 11), f, grid=None)
        # print '\n', b

        self.assertEqual(b[0, 0], 0.0)
        self.assertEqual(b[0, -1], 0.0)
        self.assertEqual(b[-1, 0], 0.0)
        self.assertEqual(b[-1, -1], 40.0)
        self.assertEqual(b[-2, -1], 30.0)
        self.assertEqual(b[-1, -2], 36.0)

    def test_interpolate_block_2(self):

        def f(i, j):
            return float(i + j)

        origin = (0, 0)
        shape = (3, 5)  # nrows, ncols

        fUL = f(0, 0)   # 0.0
        fUR = f(0, 4)   # 4.0
        fLL = f(2, 0)   # 2.0
        fLR = f(2, 4)   # 6.0

        a = bilinear(shape, fUL, fUR, fLR, fLL)
        # print '\n', a

        self.assertEqual(a[0, 0], fUL)
        self.assertEqual(a[0, -1], fUR)
        self.assertEqual(a[-1, -1], fLR)
        self.assertEqual(a[-1, 0], fLL)

        b = interpolate_block(origin, shape, f, grid=None)
        # print '\n', b

        self.assertTrue(numpy.max(b - a) < 0.000001)

    def test_interpolate_block_3(self):

        def f(i, j):
            return float(i)**2 * math.sqrt(float(j))

        origin = (0, 0)
        shape = (3, 5)  # nrows, ncols

        fUL = f(0, 0)
        fUR = f(0, 4)
        fLL = f(2, 0)
        fLR = f(2, 4)

        a = bilinear(shape, fUL, fUR, fLR, fLL)
        # print '\n', a

        self.assertEqual(a[0, 0], fUL)
        self.assertEqual(a[0, -1], fUR)
        self.assertEqual(a[-1, -1], fLR)
        self.assertEqual(a[-1, 0], fLL)

        b = interpolate_block(origin, shape, f, grid=None)
        # print '\n', b

        self.assertTrue(numpy.max(b - a) < 0.000001)


def the_suite():
    """Returns a test suite of all the tests in this module."""

    test_classes = [
        BLRBTest
    ]

    suite_list = map(unittest.defaultTestLoader.loadTestsFromTestCase,
                     test_classes)

    suite = unittest.TestSuite(suite_list)

    return suite

#
# Run unit tests if in __main__
#

if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(the_suite())
