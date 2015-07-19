#!/usr/bin/env python

import argparse
import glob
import os
from os.path import join as pjoin, abspath
import unittest

import numpy.testing as npt
import rasterio

from gaip.tests.unittesting_tools import ParameterisedTestCaseFiles


class TestFCProducts(ParameterisedTestCaseFiles):

    """
    Test and compare the FC product outputs.
    """

    def test_photosynthetic_veg_band(self):
        """
        Test and compare the differences between the
        reference and test photosynthetic vegetation
        output bands.
        """
        with rasterio.open(self.reference_fname, 'r') as ref_ds:
            ref_img = ref_ds.read_band(2, masked=False)

        with rasterio.open(self.test_fname, 'r') as test_ds:
            test_img = test_ds.read_band(2, masked=False)

        # Precision
        ip = self.integer_precision

        self.assertIsNone(npt.assert_all_close(test_img, ref_img, atol=ip))


    def test_non_photosynthetic_veg_band(self):
        """
        Test and compare the differences between the
        reference and test non-photosynthetic vegetation
        output bands.
        """
        with rasterio.open(self.reference_fname, 'r') as ref_ds:
            ref_img = ref_ds.read_band(3, masked=False)

        with rasterio.open(self.test_fname, 'r') as test_ds:
            test_img = test_ds.read_band(3, masked=False)

        # Precision
        ip = self.integer_precision

        self.assertIsNone(npt.assert_all_close(test_img, ref_img, atol=ip))


    def test_bare_soil_band(self):
        """
        Test and compare the differences between the
        reference and test bare soil output bands.
        """
        with rasterio.open(self.reference_fname, 'r') as ref_ds:
            ref_img = ref_ds.read_band(1, masked=False)

        with rasterio.open(self.test_fname, 'r') as test_ds:
            test_img = test_ds.read_band(1, masked=False)

        # Precision
        ip = self.integer_precision

        self.assertIsNone(npt.assert_all_close(test_img, ref_img, atol=ip))


    def test_unmixing_error_band(self):
        """
        Test and compare the differences between the
        reference and test unmixing error output bands.
        """
        with rasterio.open(self.reference_fname, 'r') as ref_ds:
            ref_img = ref_ds.read_band(4, masked=False)

        with rasterio.open(self.test_fname, 'r') as test_ds:
            test_img = test_ds.read_band(4, masked=False)

        # Precision
        ip = self.integer_precision

        self.assertIsNone(npt.assert_all_close(test_img, ref_img, atol=ip))


if __name__ == '__main__':

    desc = 'Test and compare the fractional cover outputs.'
    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument('--reference_filename', required=True,
                        help='The filename for the referenc dataset.')
    parser.add_argument('--test_filename', required=True,
                        help=('The filename for the test dataset.'))
    parser.add_argument('--int_precision', default=1, type=int,
                        help='The integer precision used for array comparison')

    parsed_args = parser.parse_args()

    reference_fname = parsed_args.reference_filename
    test_fname = parsed_args.test_filename
    int_precision = parsed_args.int_precision


    print "Testing the Fractional cover products."
    suite = unittest.TestSuite()
    suite.addTest(ParameterisedTestCase.parameterise(TestFCProducts,
                  reference_fname=reference_fname, test_fname=test_fname,
                  integer_precision=int_precision))
    unittest.TextTestRunner(verbosity=2).run(suite)
