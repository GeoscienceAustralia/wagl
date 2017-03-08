#!/usr/bin/env python

import unittest

from gaip.tests import unittesting_tools as ut

class TestRandomPixelLocations(unittest.TestCase):

    def test_non_2D(self):
        """
        Test that specifying a non 2D tuple raises a TypeEror.
        """

        dims = (3,100,100)
        self.assertRaises(TypeError, ut.randomPixelLocations, dims)

    def test_nPixels(self):
        """
        Test that the correct number of random pixels are returned.
        The default return value is 100 pixels.
        """

        dims = (100,100)
        idx = ut.randomPixelLocations(dims)
        n = idx[0].shape[0]
        self.assertTrue(n == 100)

if __name__ == '__main__':
    unittest.main()
