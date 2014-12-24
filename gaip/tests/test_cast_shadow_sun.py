#!/usr/bin/env python

import argparse
import os
from os.path import join as pjoin
from os.path import exists as pexists
import unittest

import numpy.testing as npt

from gaip import acquisitions
from gaip import find_file
from gaip import read_img
from gaip.tests.unittesting_tools import ParameterisedTestCase

#TODO Filename to be determined from the nbar.cfg file

def calculate_cast_shadow_sun(acquisition, ref_dir, outdir, buffer=250,
        block_height=500, block_width=500):
    """
    Calculates the cast shadow from the direction of the sun.
    """
    # TC_Intermediates directory
    tc_dir = pjoin(ref_dir, 'TC_Intermediates')
    tc_outdir = pjoin(outdir, 'TC_Intermediates')

    # Check and load the required files from disk
    fname_solar_zenith  = find_file(ref_dir, 'SOLAR_ZENITH.bin')
    fname_solar_azimuth = find_file(ref_dir, 'SOLAR_AZIMUTH.bin')
    fname_smoothed_dsm  = find_file(tc_dir, 'dsm_subset_smoothed.bin')

    # Define the output filename
    outfname = pjoin(tc_outdir, 'cast_shadow_sun.bin')

    # Compute the satellite view cast shadow mask
    calculate_cast_shadow(acquisition, fname_smoothed_dsm, buffer,
        block_height, block_width, fname_solar_zenith,
        fname_solar_azimuth, outfname)


class TestCasthadowSunFileNames(ParameterisedTestCase):
    """
    Unittests will occur for the following files:
    cast_shadow_sun.bin
    """

    # File of interest
    ParameterisedTestCase.fname_cast_shadow_sun = 'cast_shadow_sun.bin'

    def test_cast_shadow_sun_ref(self):
        """
        Check that the cast shadow sun reference file exists.
        """
        # TC_Intermediates directory
        tc_dir = pjoin(self.reference_dir, 'TC_Intermediates')

        fname = pjoin(tc_dir, self.fname_cast_shadow_sun)
        msg = 'Reference file does not exist: {fname}'.format(fname=fname)
        self.assertIs(pexists(fname), True,
                      'Reference file does not exist: %s'%fname)

    def test_cast_shadow_sun_tst(self):
        """
        Check that the cast shadow test file exists.
        """
        # TC_Intermediates directory
        tc_dir = pjoin(self.test_dir, 'TC_Intermediates')

        fname = pjoin(tc_dir, self.fname_cast_shadow_sun)
        msg = 'Test file does not exist: {fname}'.format(fname=fname)
        self.assertIs(pexists(fname), True,
                      'Test file does not exist: %s'%fname)


class TestCastShadowSunOutputs(ParameterisedTestCase):
    """
    Unittests will occur for the following files:
    cast_shadow_sun.bin
    """

    # File of interest
    ParameterisedTestCase.fname_cast_shadow_sun = 'cast_shadow_sun.bin'

    def test_cast_shadow_sun(self):
        """
        Test the cast shadow sun mask against the reference mask.
        """
        # TC_Intermediates directory (reference and test)
        tc_ref_dir = pjoin(self.reference_dir, 'TC_Intermediates')
        tc_tst_dir = pjoin(self.test_dir, 'TC_Intermediates')

        # Get the filenames for both the reference and test files
        ref_fname  = find_file(tc_ref_dir, self.fname_cast_shadow_sun)
        test_fname = find_file(tc_tst_dir, self.fname_cast_shadow_sun)

        # Get the image data
        ref_img  = read_img(ref_fname)
        test_img = read_img(test_fname)

        # Precision
        dp = self.decimal_precision

        self.assertIsNone(npt.assert_almost_equal(test_img, ref_img,
                                                  decimal=dp))


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(description='Perform unittesting for the cast shadow sun mask and optionaly calculates the cast shadow sun mask.')

    parser.add_argument('--L1T_dir', required=True, help='A directory path of a L1T scene.')
    parser.add_argument('--nbar_work_dir', required=True, help='A directory path to the associated NBAR working directory.')
    parser.add_argument('--outdir', required=True, help='A directory path that will contain the output files.')
    parser.add_argument('--dec_precision', default=4, type=int, help='The decimal precision used for array comparison')
    parser.add_argument('--int_precision', default=1, type=int, help='The integer precision used for array comparison')
    parser.add_argument('--compute', action='store_true', help='If set then the self shadow array will be computed before running the unittests.')
    parser.add_argument('--buffer', default=250, help='The buffer in pixels to be used in calculating the cast shadow sun mask.')
    parser.add_argument('--block_x', default=500, help='The x block size in pixels (Twice the buffer).')
    parser.add_argument('--block_y', default=500, help='The y block size in pixels (Twice the buffer)..')

    parsed_args = parser.parse_args()

    L1T_dir       = parsed_args.L1T_dir
    nbar_work_dir = parsed_args.nbar_work_dir
    outdir        = parsed_args.outdir
    dec_precision = parsed_args.dec_precision
    int_precision = parsed_args.int_precision
    compute       = parsed_args.compute
    buffer        = parsed_args.buffer
    block_x       = parsed_args.block_x
    block_y       = parsed_args.block_y

    if compute:
        # Check the output directory
        if not pexists(outdir):
            os.makedirs(outdir)

        # Get the current directory
        cwd = os.getcwd()

        # Change to the output directory that will contain the results
        os.chdir(outdir)

        # Open the L1T dataset
        acqs = acquisitions(L1T_dir)

        # Compute the angles
        calculate_self_shadow(acqs[0], nbar_work_dir, outdir, buffer,
                              block_y, block_x)

        # Close the L1T dataset
        acqs = None

        # Change back to the original directory
        os.chdir(cwd)

    print "Checking that we have all the reference and test data files neccessary."
    suite = unittest.TestSuite()
    suite.addTest(ParameterisedTestCase.parameterise(
                  TestCasthadowSunFileNames,
                  reference_dir=nbar_work_dir, test_dir=outdir,
                  decimal_precision=dec_precision,
                  integer_precision=int_precision))
    unittest.TextTestRunner(verbosity=2).run(suite)

    print "Comparing the reference and test cast shadow sun masks."
    suite = unittest.TestSuite()
    suite.addTest(ParameterisedTestCase.parameterise(TestCastShadowSunOutputs,
                  reference_dir=nbar_work_dir, test_dir=outdir,
                  decimal_precision=dec_precision,
                  integer_precision=int_precision))
    unittest.TextTestRunner(verbosity=2).run(suite)

