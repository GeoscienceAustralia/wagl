#!/usr/bin/env python

import unittest

import numpy
import numpy.testing as npt

from gaip import unmiximage

class TestUnmix(unittest.TestCase):
    """
    Test the unmiximge Fortran module.
    """

    def test_unimx(self):
        """Test that the unmix module returns the expected result."""

        m = 6 # number of bands.
        n = 3 # number of endmembers.
        
        # number of rows and columns in image.
        numRows = 3
        numCols = 2
        
        # create an image with 3 rows 2 columns and 6 layers.
        image = numpy.zeros(m*numRows*numCols)
        image.shape = (m, numRows, numCols)
        
        # initialize the pixels in the image - left to right, top to bottom.
        image[:,0,0] = numpy.array([0.08, 0.1655, 0.1415, 0.225, 0.29, 0.22])
        image[:,0,1] = numpy.array([0.06, 0.14, 0.15, 0.23, 0.29, 0.20])
        image[:,1,0] = numpy.array([0.07, 0.16, 0.14, 0.11, 0.31, 0.25])
        image[:,1,1] = numpy.array([0.10, 0.20, 0.21, 0.19, 0.29, 0.21])
        image[:,2,0] = numpy.array([0.07, 0.19, 0.11, 0.26, 0.22, 0.29])
        image[:,2,1] = numpy.array([0.0, 0.0 , 0.0, 0.0, 0.0, 0.0])
        
        # A is a mxn array of endmembers
        A = numpy.array([[0.1, 0.075, 0.05],
                         [0.1, 0.1,  0.2],
                         [0.2, 0.15, 0.15],
                         [0.25, 0.20, 0.40],
                         [0.40, 0.25, 0.1],
                         [0.35, 0.25, 0.1]], dtype='float64')

        # Compute the fraction components
        fractions = unmiximage.unmiximage(image, A, 0.0, -10.0)

        # The reference image
        ref = numpy.array([0.44815772, 0.25296948, 0.21086879, 0.06933006,
                           0.58656233, 0, 0.2432272 , 0.05688818, 0.66786787,
                           0.06846847, 0, 0.11460207, 0, 0.99988303, 0.107227,
                           0.12129509, 0, 0.8780015, 0.23776422, 0.09579347,
                           -10, -10, -10, -10], dtype='float64')

        ref = ref.reshape(numRows, numCols, 4)
        ref = ref.transpose(2,0,1)

        self.assertIsNone(npt.assert_almost_equal(fractions, ref, decimal=6))

if __name__ == '__main__':

    unittest.main()
