#!/usr/bin/env python

import argparse
import os
from os.path import join
from os.path import dirname
import unittest

import numpy as np
import numpy.testing as npt

from gaip import acquisitions
from gaip import filter_dsm
from gaip import find_file
from gaip import gridded_geo_box
from gaip import read_img
from gaip import write_img
from gaip.tests.unittesting_tools import ParameterisedTestCase


def calculate_smoothed_dsm(geobox, ref_dir, outdir):
    """
    Creates a smoothed DSM.
    """

    # Check and load the required files from disk
    fname_dsm = 'region_dsm_image.img'

    dsm = (read_img(find_file(ref_dir, fname_dsm))).astype('float32')

    smoothed_dsm = filter_dsm(dsm)

    # Write out the smoothed dsm file
    out_fname = os.path.join(outdir, 'region_dsm_image_smoothed.img')
    write_img(smoothed_dsm, out_fname, geobox=geobox)


class TestFilterFileNames(ParameterisedTestCase):
    """
    Unittests will occur for the following files:
    region_dsm_image_smoothed.img
    """

    # File of interest
    ParameterisedTestCase.fname_smoothed_dsm = 'region_dsm_image_smoothed.img'

    def test_smoothed_dsm_ref(self):
        """
        Check that the smoothed dsm reference file exists.
        """

        fname = os.path.join(self.reference_dir, self.fname_smoothed_dsm)
        self.assertIs(os.path.exists(fname), True,
                      'Reference file does not exist: %s'%fname)

    def test_smoothed_dsm_tst(self):
        """
        Check that the smoothed dsm test file exists.
        """

        fname = os.path.join(self.test_dir, self.fname_smoothed_dsm)
        self.assertIs(os.path.exists(fname), True,
                      'Test file does not exist: %s'%fname)


class TestFilterOutputs(ParameterisedTestCase):
    """
    Unittests will occur for the following files:
    region_dsm_image_smoothed.img
    """

    # File of interest
    ParameterisedTestCase.fname_smoothed_dsm = 'region_dsm_image_smoothed.img'

    def test_smoothed_dsm(self):
        """
        Test the smoothed dsm image against the reference image.
        """

        # Get the filenames for both the reference and test files
        ref_fname  = find_file(self.reference_dir, self.fname_smoothed_dsm)
        test_fname = find_file(self.test_dir, self.fname_smoothed_dsm)

        # Get the image data
        ref_img  = read_img(ref_fname)
        test_img = read_img(test_fname)

        # Precision
        dp = self.decimal_precision

        self.assertIsNone(npt.assert_almost_equal(test_img, ref_img,
                                                  decimal=dp))


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(description='Performs unittests against the smoothed DSM and optionally calculates the smoothed DSM.')

    parser.add_argument('--L1T_dir', required=True, help='A directory path of a L1T scene.')
    parser.add_argument('--nbar_work_dir', required=True, help='A directory path to the associated NBAR working directory.')
    parser.add_argument('--outdir', required=True, help='A directory path that will contain the output files.')
    parser.add_argument('--dec_precision', default=3, type=int, help='The decimal precision used for array comparison')
    parser.add_argument('--int_precision', default=1, type=int, help='The integer precision used for array comparison')
    parser.add_argument('--compute', action='store_true', help='If set then a smoothed dsm image will be created.')

    parsed_args = parser.parse_args()

    L1T_dir       = parsed_args.L1T_dir
    nbar_work_dir = parsed_args.nbar_work_dir
    outdir        = parsed_args.outdir
    dec_precision = parsed_args.dec_precision
    int_precision = parsed_args.int_precision
    compute       = parsed_args.compute

    if compute:
        # Check the output directory
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        # Get the current directory
        cwd = os.getcwd()

        # Change to the output directory that will contain the results
        os.chdir(outdir)

        # Open the L1T dataset
        acqs = acquisitions(L1T_dir)

        # Get a geobox of the 1st acquisition
        geobox = gridded_geo_box(acqs[0])

        # Compute the smoothed dsm
        calculate_smoothed_dsm(geobox, nbar_work_dir, outdir)

        # Close the L1T dataset
        acqs = None

        # Change back to the original directory
        os.chdir(cwd)

    print "Checking that we have all the reference and test data files neccessary."
    suite = unittest.TestSuite()
    suite.addTest(ParameterisedTestCase.parameterise(TestFilterFileNames,
                  reference_dir=nbar_work_dir, test_dir=outdir,
                  decimal_precision=dec_precision,
                  integer_precision=int_precision))
    unittest.TextTestRunner(verbosity=2).run(suite)

    print "Comparing the reference and test smoothed dsm output files."
    suite = unittest.TestSuite()
    suite.addTest(ParameterisedTestCase.parameterise(TestFilterOutputs,
                  reference_dir=nbar_work_dir, test_dir=outdir,
                  decimal_precision=dec_precision,
                  integer_precision=int_precision))
    unittest.TextTestRunner(verbosity=2).run(suite)

