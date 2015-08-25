#!/usr/bin/env python

import argparse
import os
from os.path import join as pjoin
from os.path import exists as pexists
import unittest

import numpy.testing as npt

from gaip import acquisitions
from gaip import self_shadow
from gaip import find_file
from gaip import read_img
from gaip.tests.unittesting_tools import ParameterisedTestCase

#TODO Filename to be determined from the nbar.cfg file

def compute_self_shadow(ref_dir, out_dir):
    """
    A small wrapper for running the slope and self shadow function,
    specifically for the unittests.
    """
    # TC_Intermediates directory
    tc_dir = pjoin(ref_dir, 'TC_Intermediates')
    tc_outdir = pjoin(out_dir, 'TC_Intermediates')
    if not pexists(tc_outdir):
        os.makedirs(tc_outdir)

    # Check and load the required files from disk
    fname_incident_angle = find_file(tc_dir, 'incident_angle.bin')
    fname_exiting_angle = find_file(tc_dir, 'exiting_angle.bin')

    # Define the output filenames
    out_fname = pjoin(tc_outdir, 'self_shadow_mask.bin')

    self_shadow(fname_incident_angle, fname_exiting_angle, out_fname)


class TestSelfShadowMaskFileNames(ParameterisedTestCase):

    """
    Unittests will occur for the following files:
        * self_shadow_mask.bin
    """

    ParameterisedTestCase.fname_self_shadow_mask = 'self_shadow_mask.bin'


    def test_self_shadow_mask_ref(self):
        """
        Check that the self shadow mask reference file exists.
        """
        # TC_Intermediates directory
        tc_dir = pjoin(self.reference_dir, 'TC_Intermediates')

        fname = pjoin(tc_dir, self.fname_self_shadow_mask)
        msg = 'Reference file does not exist: {fname}'.format(fname=fname)
        self.assertIs(pexists(fname), True, msg)


    def test_self_shadow_mask_tst(self):
        """
        Check that the self shadow mask test file exists.
        """
        # TC_Intermediates directory
        tc_dir = pjoin(self.test_dir, 'TC_Intermediates')

        fname = pjoin(tc_dir, self.fname_self_shadow_mask)
        msg = 'Reference file does not exist: {fname}'.format(fname=fname)
        self.assertIs(pexists(fname), True, msg)


class TestSelfShadowMaskOutput(ParameterisedTestCase):

    """
    Unittests will occur for the following files:
        * self_shadow_mask.bin
    """

    ParameterisedTestCase.fname_self_shadow_mask = 'self_shadow_mask.bin'

    def test_self_shadow_mask(self):
        """
        Test the self shadow mask against the reference mask.
        """
        # TC_Intermediates directory (reference and test)
        tc_ref_dir = pjoin(self.reference_dir, 'TC_Intermediates')
        tc_tst_dir = pjoin(self.test_dir, 'TC_Intermediates')

        # Get the filenames for both the reference and test files
        ref_fname  = find_file(tc_ref_dir, self.fname_self_shadow_mask)
        test_fname = find_file(tc_tst_dir, self.fname_self_shadow_mask)

        # Get the image data
        ref_img  = read_img(ref_fname)
        test_img = read_img(test_fname)

        # Precision
        int_prec = self.integer_precision

        self.assertIsNone(npt.assert_allclose(test_img, ref_img,
                                              atol=int_prec))


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(description=('Perform unittesting for '
                                                  'the self shadow masks.'))
    parser.add_argument('--L1T_dir', required=True,
                        help='A directory path of a L1T scene.')
    parser.add_argument('--nbar_work_dir', required=True,
                        help=('A directory path to the associated NBAR '
                              'working directory.'))
    parser.add_argument('--outdir', required=True,
                        help=('A directory path that will contain the '
                              'output files.'))
    parser.add_argument('--dec_precision', default=4, type=int,
                        help='The decimal precision used for array comparison')
    parser.add_argument('--int_precision', default=1, type=int,
                        help='The integer precision used for array comparison')
    parser.add_argument('--compute', action='store_true',
                        help=('If set then the self shadow array will be '
                              'computed before running the unittests.'))

    parsed_args = parser.parse_args()

    L1T_dir = parsed_args.L1T_dir
    nbar_work_dir = parsed_args.nbar_work_dir
    outdir = parsed_args.outdir
    dec_precision = parsed_args.dec_precision
    int_precision = parsed_args.int_precision
    compute = parsed_args.compute

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
        compute_self_shadow(nbar_work_dir, outdir)

        # Close the L1T dataset
        acqs = None

        # Change back to the original directory
        os.chdir(cwd)

    print "Checking that we have the reference and test self shadow masks."
    suite = unittest.TestSuite()
    suite.addTest(ParameterisedTestCase.parameterise(
                  TestSelfShadowMaskFileNames,
                  reference_dir=nbar_work_dir, test_dir=outdir,
                  decimal_precision=dec_precision,
                  integer_precision=int_precision))
    unittest.TextTestRunner(verbosity=2).run(suite)

    print "Comparing the reference and test self shadow masks."
    suite = unittest.TestSuite()
    suite.addTest(ParameterisedTestCase.parameterise(
                  TestSelfShadowMaskOutput,
                  reference_dir=nbar_work_dir, test_dir=outdir,
                  decimal_precision=dec_precision,
                  integer_precision=int_precision))
    unittest.TextTestRunner(verbosity=2).run(suite)
