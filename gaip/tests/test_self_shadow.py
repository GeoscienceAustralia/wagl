#!/usr/bin/env python

import argparse
import os
from os.path import join as pjoin
from os.path import exists as pexists
import unittest

import numpy.testing as npt

from gaip import acquisitions
from gaip import calculate_self_shadow
from gaip import find_file
from gaip import read_img
from gaip.tests.unittesting_tools import ParameterisedTestCase

#TODO Filename to be determined from the nbar.cfg file

def calculate_self_shadow(acquisition, ref_dir, outdir, buffer=250):
    """
    A small wrapper for running the slope and self shadow function,
    specifically for the unittests.
    """
    # TC_Intermediates directory
    tc_dir = pjoin(ref_dir, 'TC_Intermediates')
    tc_outdir = pjoin(outdir, 'TC_Intermediates')

    # Check and load the required files from disk
    fname_solar_zenith  = find_file(ref_dir, 'SOLAR_ZENITH.bin')
    fname_solar_azimuth = find_file(ref_dir, 'SOLAR_AZIMUTH.bin')
    fname_satellite_view   = find_file(ref_dir, 'SATELLITE_VIEW.bin')
    fname_satellte_azimuth = find_file(ref_dir, 'SATELLITE_AZIMUTH.bin')
    fname_smoothed_dsm  = find_file(tc_dir, 'dsm_subset_smoothed.bin')

    # Define the output filenames
    fnames = ['self_shadow_mask.bin', 'slope.bin', 'aspect.bin',
              'incident_angle.bin', 'exiting_angle.bin',
              'azimuth_incident_angle.bin', 'azimuth_exiting_angle.bin',
              'relative_slope.bin']

    outfnames = [pjoin(tc_outdir, name) for name in fnames]

    calculate_self_shadow(acquisition, fname_smoothed_dsm, buffer,
         fname_solar_zenith, fname_solar_azimuth, fname_satellite_view,
         fname_satellte_azimuth, outfnames)


class TestSelfShadowSlopeFileNames(ParameterisedTestCase):
    """
    Unittests will occur for the following files:
        * self_shadow_mask.bin
        * slope.bin
        * aspect.bin
        * incident_angle.bin
        * exiting_angle.bin
        * azimuth_incident_angle.bin
        * azimuth_exiting_angle.bin
        * relative_slope.bin
    """

    ParameterisedTestCase.fname_self_shadow_mask = 'self_shadow_mask.bin'
    ParameterisedTestCase.fname_slope = 'slope.bin'
    ParameterisedTestCase.fname_aspect = 'aspect.bin'
    ParameterisedTestCase.fname_incident_angle = 'incident_angle.bin'
    ParameterisedTestCase.fname_exiting_angle = 'exiting_angle.bin'
    ParameterisedTestCase.fname_azi_inci_angle = 'azimuth_incident_angle.bin'
    ParameterisedTestCase.fname_azi_exit_angle = 'azimuth_exiting_angle.bin'
    ParameterisedTestCase.fname_relative_slope = 'relative_slope.bin'

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
        tc_dir = pjoin(self.reference_dir, 'TC_Intermediates')

        fname = pjoin(tc_dir, self.fname_self_shadow_mask)
        msg = 'Reference file does not exist: {fname}'.format(fname=fname)
        self.assertIs(pexists(fname), True, msg)

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

    def test_incident_angle_ref(self):
        """
        Check that the incident angle reference file exists.
        """
        # TC_Intermediates directory
        tc_dir = pjoin(self.reference_dir, 'TC_Intermediates')

        fname = pjoin(tc_dir, self.fname_incident_angle)
        msg = 'Reference file does not exist: {fname}'.format(fname=fname)
        self.assertIs(pexists(fname), True, msg)

    def test_incident_angle_tst(self):
        """
        Check that the incident angle test file exists.
        """
        # TC_Intermediates directory
        tc_dir = pjoin(self.reference_dir, 'TC_Intermediates')

        fname = pjoin(tc_dir, self.fname_incident_angle)
        msg = 'Reference file does not exist: {fname}'.format(fname=fname)
        self.assertIs(pexists(fname), True, msg)

    def test_exiting_angle_ref(self):
        """
        Check that the exiting angle reference file exists.
        """
        # TC_Intermediates directory
        tc_dir = pjoin(self.reference_dir, 'TC_Intermediates')

        fname = pjoin(tc_dir, self.fname_exiting_angle)
        msg = 'Reference file does not exist: {fname}'.format(fname=fname)
        self.assertIs(pexists(fname), True, msg)

    def test_exiting_angle_tst(self):
        """
        Check that the exiting angle test file exists.
        """
        # TC_Intermediates directory
        tc_dir = pjoin(self.reference_dir, 'TC_Intermediates')

        fname = pjoin(tc_dir, self.fname_exiting_angle)
        msg = 'Reference file does not exist: {fname}'.format(fname=fname)
        self.assertIs(pexists(fname), True, msg)

    def test_azimuth_incident_angle_ref(self):
        """
        Check that the azimuth incident angle reference file exists.
        """
        # TC_Intermediates directory
        tc_dir = pjoin(self.reference_dir, 'TC_Intermediates')

        fname = pjoin(tc_dir, self.fname_azi_inci_angle)
        msg = 'Reference file does not exist: {fname}'.format(fname=fname)
        self.assertIs(pexists(fname), True, msg)

    def test_azimuth_incident_angle_tst(self):
        """
        Check that the azimuth incident angle test file exists.
        """
        # TC_Intermediates directory
        tc_dir = pjoin(self.reference_dir, 'TC_Intermediates')

        fname = pjoin(tc_dir, self.fname_azi_inci_angle)
        msg = 'Reference file does not exist: {fname}'.format(fname=fname)
        self.assertIs(pexists(fname), True, msg)

    def test_azimuth_exiting_angle_ref(self):
        """
        Check that the azimuth exiting angle reference file exists.
        """
        # TC_Intermediates directory
        tc_dir = pjoin(self.reference_dir, 'TC_Intermediates')

        fname = pjoin(tc_dir, self.fname_azi_exit_angle)
        msg = 'Reference file does not exist: {fname}'.format(fname=fname)
        self.assertIs(pexists(fname), True, msg)

    def test_azimuth_exiting_angle_tst(self):
        """
        Check that the azimuth exiting angle test file exists.
        """
        # TC_Intermediates directory
        tc_dir = pjoin(self.reference_dir, 'TC_Intermediates')

        fname = pjoin(tc_dir, self.fname_azi_exit_angle)
        msg = 'Reference file does not exist: {fname}'.format(fname=fname)
        self.assertIs(pexists(fname), True, msg)

    def test_relative_slope_ref(self):
        """
        Check that the relative slope reference file exists.
        """
        # TC_Intermediates directory
        tc_dir = pjoin(self.reference_dir, 'TC_Intermediates')

        fname = pjoin(tc_dir, self.fname_relative_slope)
        msg = 'Reference file does not exist: {fname}'.format(fname=fname)
        self.assertIs(pexists(fname), True, msg)

    def test_relative_slope_tst(self):
        """
        Check that the relative slope test file exists.
        """
        # TC_Intermediates directory
        tc_dir = pjoin(self.reference_dir, 'TC_Intermediates')

        fname = pjoin(tc_dir, self.fname_relative_slope)
        msg = 'Reference file does not exist: {fname}'.format(fname=fname)
        self.assertIs(pexists(fname), True, msg)


class TestSelfShadowSlopeOutputs(ParameterisedTestCase):
    """
    Unittests will occur for the following files:
        * self_shadow_mask.bin
        * slope.bin
        * aspect.bin
        * incident_angle.bin
        * exiting_angle.bin
        * azimuth_incident_angle.bin
        * azimuth_exiting_angle.bin
        * relative_slope.bin
    """

    ParameterisedTestCase.fname_self_shadow_mask = 'self_shadow_mask.bin'
    ParameterisedTestCase.fname_slope = 'slope.bin'
    ParameterisedTestCase.fname_aspect = 'aspect.bin'
    ParameterisedTestCase.fname_incident_angle = 'incident_angle.bin'
    ParameterisedTestCase.fname_exiting_angle = 'exiting_angle.bin'
    ParameterisedTestCase.fname_azi_inci_angle = 'azimuth_incident_angle.bin'
    ParameterisedTestCase.fname_azi_exit_angle = 'azimuth_exiting_angle.bin'
    ParameterisedTestCase.fname_relative_slope = 'relative_slope.bin'

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

    def test_incident_angle(self):
        """
        Test the incident angle image against the reference image.
        """
        # TC_Intermediates directory (reference and test)
        tc_ref_dir = pjoin(self.reference_dir, 'TC_Intermediates')
        tc_tst_dir = pjoin(self.test_dir, 'TC_Intermediates')

        # Get the filenames for both the reference and test files
        ref_fname  = find_file(tc_ref_dir, self.fname_incident_angle)
        test_fname = find_file(tc_tst_dir, self.fname_incident_angle)

        # Get the image data
        ref_img  = read_img(ref_fname)
        test_img = read_img(test_fname)

        # Precision
        dp = self.decimal_precision

        self.assertIsNone(npt.assert_almost_equal(test_img, ref_img,
                                                  decimal=dp))

    def test_exiting_angle(self):
        """
        Test the exiting angle image against the reference image.
        """
        # TC_Intermediates directory (reference and test)
        tc_ref_dir = pjoin(self.reference_dir, 'TC_Intermediates')
        tc_tst_dir = pjoin(self.test_dir, 'TC_Intermediates')

        # Get the filenames for both the reference and test files
        ref_fname  = find_file(tc_ref_dir, self.fname_exiting_angle)
        test_fname = find_file(tc_tst_dir, self.fname_)

        # Get the image data
        ref_img  = read_img(ref_fname)
        test_img = read_img(test_fname)

        # Precision
        dp = self.decimal_precision

        self.assertIsNone(npt.assert_almost_equal(test_img, ref_img,
                                                  decimal=dp))

    def test_azimuth_incident_angle(self):
        """
        Test the azimuth incident angle image against the reference
        image.
        """
        # TC_Intermediates directory (reference and test)
        tc_ref_dir = pjoin(self.reference_dir, 'TC_Intermediates')
        tc_tst_dir = pjoin(self.test_dir, 'TC_Intermediates')

        # Get the filenames for both the reference and test files
        ref_fname  = find_file(tc_ref_dir, self.fname_azi_inci_angle)
        test_fname = find_file(tc_tst_dir, self.fname_azi_inci_angle)

        # Get the image data
        ref_img  = read_img(ref_fname)
        test_img = read_img(test_fname)

        # Precision
        dp = self.decimal_precision

        self.assertIsNone(npt.assert_almost_equal(test_img, ref_img,
                                                  decimal=dp))

    def test_azimuth_exiting_angle(self):
        """
        Test the azimuth exiting angle image against the reference
        image.
        """
        # TC_Intermediates directory (reference and test)
        tc_ref_dir = pjoin(self.reference_dir, 'TC_Intermediates')
        tc_tst_dir = pjoin(self.test_dir, 'TC_Intermediates')

        # Get the filenames for both the reference and test files
        ref_fname  = find_file(tc_ref_dir, self.fname_azi_exit_angle)
        test_fname = find_file(tc_tst_dir, self.fname_azi_exit_angle)

        # Get the image data
        ref_img  = read_img(ref_fname)
        test_img = read_img(test_fname)

        # Precision
        dp = self.decimal_precision

        self.assertIsNone(npt.assert_almost_equal(test_img, ref_img,
                                                  decimal=dp))

    def test_relative_slope(self):
        """
        Test the relative_slope image against the reference image.
        """
        # TC_Intermediates directory (reference and test)
        tc_ref_dir = pjoin(self.reference_dir, 'TC_Intermediates')
        tc_tst_dir = pjoin(self.test_dir, 'TC_Intermediates')

        # Get the filenames for both the reference and test files
        ref_fname  = find_file(tc_ref_dir, self.fname_relative_slope)
        test_fname = find_file(tc_tst_dir, self.fname_relative_slope)

        # Get the image data
        ref_img  = read_img(ref_fname)
        test_img = read_img(test_fname)

        # Precision
        dp = self.decimal_precision

        self.assertIsNone(npt.assert_almost_equal(test_img, ref_img,
                                                  decimal=dp))

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(description='Perform unittesting for the self shadow and various slope angles.')

    parser.add_argument('--L1T_dir', required=True, help='A directory path of a L1T scene.')
    parser.add_argument('--nbar_work_dir', required=True, help='A directory path to the associated NBAR working directory.')
    parser.add_argument('--outdir', required=True, help='A directory path that will contain the output files.')
    parser.add_argument('--dec_precision', default=4, type=int, help='The decimal precision used for array comparison')
    parser.add_argument('--int_precision', default=1, type=int, help='The integer precision used for array comparison')
    parser.add_argument('--compute', action='store_true', help='If set then the self shadow array will be computed before running the unittests.')
    parser.add_argument('--buffer', default=250, help='The buffer in pixels to be used in calculating the cast shadow sun mask.')

    parsed_args = parser.parse_args()

    L1T_dir       = parsed_args.L1T_dir
    nbar_work_dir = parsed_args.nbar_work_dir
    outdir        = parsed_args.outdir
    dec_precision = parsed_args.dec_precision
    int_precision = parsed_args.int_precision
    compute       = parsed_args.compute
    buffer        = parsed_args.buffer

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
        calculate_self_shadow(acqs[0], nbar_work_dir, outdir, buffer)

        # Close the L1T dataset
        acqs = None

        # Change back to the original directory
        os.chdir(cwd)

    print "Checking that we have all the reference and test data files neccessary."
    suite = unittest.TestSuite()
    suite.addTest(ParameterisedTestCase.parameterise(
                  TestSelfShadowSlopeFileNames,
                  reference_dir=nbar_work_dir, test_dir=outdir,
                  decimal_precision=dec_precision,
                  integer_precision=int_precision))
    unittest.TextTestRunner(verbosity=2).run(suite)

    print "Comparing the reference and test self shadow masks and angle outputs."
    suite = unittest.TestSuite()
    suite.addTest(ParameterisedTestCase.parameterise(
                  TestSelfShadowSlopeOutputs,
                  reference_dir=nbar_work_dir, test_dir=outdir,
                  decimal_precision=dec_precision,
                  integer_precision=int_precision))
    unittest.TextTestRunner(verbosity=2).run(suite)
