#!/usr/bin/env python

import argparse
import os
from os.path import join as pjoin
from os.path import exists as pexists
import unittests

import numpy.testing as npt

from gaip import acquisitions
from gaip import slope_aspect_arrays
from gaip import find_file
from gaip import read_img
from gaip.tests.unittesting_tools import ParameterisedTestCase

#TODO Filename to be determined from the nbar.cfg file


def compute_slope_aspect(acquisition, ref_dir, out_dir, margins):
    """
    A small wrapper for executing the slope and aspect function prior to
    running unittests.
    """
    # TC_Intermediates directory
    tc_dir = pjoin(ref_dir, 'TC_Intermediates')
    tc_outdir = pjoin(out_dir, 'TC_Intermediates')

    # Get the smoothed dsm
    dsm_fname = find_file(ref_dir, 'dsm_subset_smoothed.bin')

    # Define the output file names
    slope_out_fname = pjoin(tc_outdir, 'slope.bin')
    aspect_out_fname = pjoin(tc_outdir, 'aspect.bin')

    slope_aspect_arrays(acquisition, dsm_fname, margins, slope_out_fname,
                        aspect_out_fname)


class TestSlopeAspectFileNames(ParameterisedTestCase):

    """
    Unittests will occur for the following files:
        * slope.bin
        * aspect.bin
    """

    ParameterisedTestCase.fname_slope = 'slope.bin'
    ParameterisedTestCase.fname_aspect = 'aspect.bin'


    def test_slope_ref(self):
        """
        Check that the slope reference file exists.
        """
        # TC_Intermediates directory
        tc_dir = pjoin(self.reference_dir, 'TC_Intermediates')

        fname = pjoin(tc_dir, self.fname_slope)
        msg = 'Reference file does not exist: {fname}'.format(fname=fname)
        self.assertIs(pexists(fname), True, msg)


    def test_slope_tst(self):
        """
        Check that the slope test file exists.
        """
        # TC_Intermediates directory
        tc_dir = pjoin(self.reference_dir, 'TC_Intermediates')

        fname = pjoin(tc_dir, self.fname_slope)
        msg = 'Reference file does not exist: {fname}'.format(fname=fname)
        self.assertIs(pexists(fname), True, msg)


    def test_aspect_ref(self):
        """
        Check that the aspect reference file exists.
        """
        # TC_Intermediates directory
        tc_dir = pjoin(self.reference_dir, 'TC_Intermediates')

        fname = pjoin(tc_dir, self.fname_aspect)
        msg = 'Reference file does not exist: {fname}'.format(fname=fname)
        self.assertIs(pexists(fname), True, msg)


    def test_aspect_tst(self):
        """
        Check that the aspect test file exists.
        """
        # TC_Intermediates directory
        tc_dir = pjoin(self.reference_dir, 'TC_Intermediates')

        fname = pjoin(tc_dir, self.fname_aspect)
        msg = 'Reference file does not exist: {fname}'.format(fname=fname)
        self.assertIs(pexists(fname), True, msg)


class TestSlopeAspectOutputs(ParameterisedTestCase):

    """
    Unittests will occur for the following files:
        * slope.bin
        * aspect.bin
    """

    ParameterisedTestCase.fname_slope = 'slope.bin'
    ParameterisedTestCase.fname_aspect = 'aspect.bin'


    def test_slope(self):
        """
        Test the slope image against the reference image.
        """
        # TC_Intermediates directory (reference and test)
        tc_ref_dir = pjoin(self.reference_dir, 'TC_Intermediates')
        tc_tst_dir = pjoin(self.test_dir, 'TC_Intermediates')

        # Get the filenames for both the reference and test files
        ref_fname  = find_file(tc_ref_dir, self.fname_slope)
        test_fname = find_file(tc_tst_dir, self.fname_slope)

        # Get the image data
        ref_img  = read_img(ref_fname)
        test_img = read_img(test_fname)

        # Precision
        dp = self.decimal_precision

        self.assertIsNone(npt.assert_almost_equal(test_img, ref_img,
                                                  decimal=dp))


    def test_aspect(self):
        """
        Test the aspect image against the reference image.
        """
        # TC_Intermediates directory (reference and test)
        tc_ref_dir = pjoin(self.reference_dir, 'TC_Intermediates')
        tc_tst_dir = pjoin(self.test_dir, 'TC_Intermediates')

        # Get the filenames for both the reference and test files
        ref_fname  = find_file(tc_ref_dir, self.fname_aspect)
        test_fname = find_file(tc_tst_dir, self.fname_aspect)

        # Get the image data
        ref_img  = read_img(ref_fname)
        test_img = read_img(test_fname)

        # Precision
        dp = self.decimal_precision

        self.assertIsNone(npt.assert_almost_equal(test_img, ref_img,
                                                  decimal=dp))


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(description=('Perform unittesting for '
                                                  'the slope and aspect '
                                                  'angle arrays.'))

    parser.add_argument('--L1T_dir', required=True,
                        help='A directory path of a L1T scene.')
    parser.add_argument('--nbar_work_dir', required=True,
                        help=('A directory path to the associated NBAR '
                              'working directory.'))
    parser.add_argument('--outdir', required=True,
                        help=('A directory path that will contain the output '
                              'files.'))
    parser.add_argument('--dec_precision', default=4, type=int,
                        help='The decimal precision used for array comparison')
    parser.add_argument('--int_precision', default=1, type=int,
                        help='The integer precision used for array comparison')
    parser.add_argument('--compute', action='store_true',
                        help=('If set then the self shadow array will be '
                              'computed before running the unittests.'))
    parser.add_argument('--margins', default=250,
                        help=('The buffer in pixels to be used in calculating '
                              'the cast shadow sun mask.'))

    parsed_args = parser.parse_args()

    L1T_dir = parsed_args.L1T_dir
    nbar_work_dir = parsed_args.nbar_work_dir
    outdir = parsed_args.outdir
    dec_precision = parsed_args.dec_precision
    int_precision = parsed_args.int_precision
    compute = parsed_args.compute
    margins = parsed_args.margins

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
        compute_slope_aspect(acqs[0], nbar_work_dir, outdir, margins)

        # Close the L1T dataset
        acqs = None

        # Change back to the original directory
        os.chdir(cwd)

    print "Checking that we have the reference and test data files neccessary."
    suite = unittest.TestSuite()
    suite.addTest(ParameterisedTestCase.parameterise(
                  compute_slope_aspect,
                  reference_dir=nbar_work_dir, test_dir=outdir,
                  decimal_precision=dec_precision,
                  integer_precision=int_precision))
    unittest.TextTestRunner(verbosity=2).run(suite)

    print "Comparing the reference and test slope and aspect outputs."
    suite = unittest.TestSuite()
    suite.addTest(ParameterisedTestCase.parameterise(
                  TestSlopeAspectOutputs,
                  reference_dir=nbar_work_dir, test_dir=outdir,
                  decimal_precision=dec_precision,
                  integer_precision=int_precision))
    unittest.TextTestRunner(verbosity=2).run(suite)
