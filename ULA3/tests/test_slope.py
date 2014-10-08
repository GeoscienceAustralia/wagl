#!/usr/bin/env python

import argparse
import os
import unittest

import numpy.testing as npt

from EOtools.DatasetDrivers import SceneDataset
from ULA3._gdal_tools import Buffers
from ULA3.tc import run_slope
from unittesting_tools import find_file
from unittesting_tools import ParameterisedTestCase
from unittesting_tools import read_img


def calculate_slope_angles(scene_dataset, ref_dir, outdir, pixel_buffer=250):
    """
    Calculates all the slope angles.
    Output files produced from the calculations are:
    mask_self.img
    slope.img
    aspect.img
    incident.img
    exiting.img
    azi_incident.img
    azi_exiting.img
    rela_slope.img
    """

    # Define the pixel buffering object
    pixel_buf = Buffers(pixel_buffer)

    # Check and load the required files from disk
    fname_sat_v  = 'SAT_V.bin'
    fname_sat_az = 'SAT_AZ.bin'
    fname_sol_z  = 'SOL_Z.bin'
    fname_sol_az = 'SOL_AZ.bin'
    fname_dsm    = 'region_dsm_image_smoothed.img'

    # Read the reference images
    view_angle  = read_img(find_file(ref_dir, fname_sat_v))
    azi_angle   = read_img(find_file(ref_dir, fname_sat_az))
    solar_angle = read_img(find_file(ref_dir, fname_sol_z))
    sazi_angle  = read_img(find_file(ref_dir, fname_sol_az))
    dsm_data    = read_img(find_file(ref_dir, fname_dsm))

    # Are we looking at a UTM product?
    is_utm =  not scene_dataset.IsGeographic()

    slope_results = run_slope(scene_dataset, dsm_data, solar_angle, view_angle,
                              sazi_angle, azi_angle, pixel_buf, is_utm)

    slope_results.dump_arrays(outdir, scene_dataset, "ENVI", ".img")


class TestSlopeFilenames(ParameterisedTestCase):
    """
    Unittests will occur for the following files:
    mask_self.img
    slope.img
    aspect.img
    incident.img
    exiting.img
    azi_incident.img
    azi_exiting.img
    rela_slope.img
    """

    # Files of interest
    ParameterisedTestCase.fname_mask_self     = 'mask_self.img'
    ParameterisedTestCase.fname_slope         = 'slope.img'
    ParameterisedTestCase.fname_aspect        = 'aspect.img'
    ParameterisedTestCase.fname_incident      = 'incident.img'
    ParameterisedTestCase.fname_exiting       = 'exiting.img'
    ParameterisedTestCase.fname_azi_indcident = 'azi_incident.img'
    ParameterisedTestCase.fname_azi_exiting   = 'azi_exiting.img'
    ParameterisedTestCase.fname_rela_slope    = 'rela_slope.img'

    def test_mask_self_ref(self):
        """
        Check that the mask self reference file exists.
        """

        fname = os.path.join(self.reference_dir, self.fname_mask_self)
        self.assertIs(os.path.exists(fname), True,
                      'Reference file does not exist: %s'%fname)

    def test_mask_self_tst(self):
        """
        Check that the mask self test file exists.
        """

        fname = os.path.join(self.test_dir, self.fname_mask_self)
        self.assertIs(os.path.exists(fname), True,
                      'Test file does not exist: %s'%fname)

    def test_slope_ref(self):
        """
        Check that the slope reference file exists.
        """

        fname = os.path.join(self.reference_dir, self.fname_slope)
        self.assertIs(os.path.exists(fname), True,
                      'Reference file does not exist: %s'%fname)

    def test_slope_tst(self):
        """
        Check that the slope test file exists.
        """

        fname = os.path.join(self.test_dir, self.fname_slope)
        self.assertIs(os.path.exists(fname), True,
                      'Test file does not exist: %s'%fname)

    def test_aspect_ref(self):
        """
        Check that the aspect reference file exists.
        """

        fname = os.path.join(self.reference_dir, self.fname_aspect)
        self.assertIs(os.path.exists(fname), True,
                      'Reference file does not exist: %s'%fname)

    def test_aspect_tst(self):
        """
        Check that the aspect test file exists.
        """

        fname = os.path.join(self.test_dir, self.fname_aspect)
        self.assertIs(os.path.exists(fname), True,
                      'Test file does not exist: %s'%fname)

    def test_incident_ref(self):
        """
        Check that the incident reference file exists.
        """

        fname = os.path.join(self.reference_dir, self.fname_incident)
        self.assertIs(os.path.exists(fname), True,
                      'Reference file does not exist: %s'%fname)

    def test_incident_tst(self):
        """
        Check that the mask self test file exists.
        """

        fname = os.path.join(self.test_dir, self.fname_incident)
        self.assertIs(os.path.exists(fname), True,
                      'Test file does not exist: %s'%fname)

    def test_exiting_ref(self):
        """
        Check that the exiting reference file exists.
        """

        fname = os.path.join(self.reference_dir, self.fname_exiting)
        self.assertIs(os.path.exists(fname), True,
                      'Reference file does not exist: %s'%fname)

    def test_exiting_tst(self):
        """
        Check that the exiting self test file exists.
        """

        fname = os.path.join(self.test_dir, self.fname_exiting)
        self.assertIs(os.path.exists(fname), True,
                      'Test file does not exist: %s'%fname)

    def test_azi_incident_ref(self):
        """
        Check that the azimuth incident reference file exists.
        """

        fname = os.path.join(self.reference_dir, self.fname_azi_incident)
        self.assertIs(os.path.exists(fname), True,
                      'Reference file does not exist: %s'%fname)

    def test_azi_incident_tst(self):
        """
        Check that the azimuth incident test file exists.
        """

        fname = os.path.join(self.test_dir, self.fname_azi_incident)
        self.assertIs(os.path.exists(fname), True,
                      'Test file does not exist: %s'%fname)

    def test_azi_exiting_ref(self):
        """
        Check that the azimuth exiting reference file exists.
        """

        fname = os.path.join(self.reference_dir, self.fname_azi_exiting)
        self.assertIs(os.path.exists(fname), True,
                      'Reference file does not exist: %s'%fname)

    def test_azi_exiting_tst(self):
        """
        Check that the azimuth exiting test file exists.
        """

        fname = os.path.join(self.test_dir, self.fname_azi_exiting)
        self.assertIs(os.path.exists(fname), True,
                      'Test file does not exist: %s'%fname)

    def test_rela_slope_ref(self):
        """
        Check that the relative slope reference file exists.
        """

        fname = os.path.join(self.reference_dir, self.fname_rela_slope)
        self.assertIs(os.path.exists(fname), True,
                      'Reference file does not exist: %s'%fname)

    def test_rela_slope_tst(self):
        """
        Check that the relative slope test file exists.
        """

        fname = os.path.join(self.test_dir, self.fname_rela_slope)
        self.assertIs(os.path.exists(fname), True,
                      'Test file does not exist: %s'%fname)


class TestSlopeAngles(ParameterisedTestCase):
    """
    Unittests will occur for the following files:
    mask_self.img
    slope.img
    aspect.img
    incident.img
    exiting.img
    azi_incident.img
    azi_exiting.img
    rela_slope.img
    """

    ParameterisedTestCase.fname_mask_self     = 'mask_self.img'
    ParameterisedTestCase.fname_slope         = 'slope.img'
    ParameterisedTestCase.fname_aspect        = 'aspect.img'
    ParameterisedTestCase.fname_incident      = 'incident.img'
    ParameterisedTestCase.fname_exiting       = 'exiting.img'
    ParameterisedTestCase.fname_azi_indcident = 'azi_incident.img'
    ParameterisedTestCase.fname_azi_exiting   = 'azi_exiting.img'
    ParameterisedTestCase.fname_rela_slope    = 'rela_slope.img'

    def test_mask_self(self):
        """
        Test the self mask array.
        """

        # Get the filenames for both the reference and test files
        ref_fname  = find_file(self.reference_dir, self.fname_mask_self)
        test_fname = find_file(self.test_dir, self.fname_mask_self)

        # Get the image data
        ref_img  = read_img(ref_fname)
        test_img = read_img(test_fname)

        # Precision
        dp = self.decimal_precision

        self.assertIsNone(npt.assert_almost_equal(test_img, ref_img,
                                                  decimal=dp))

    def test_slope(self):
        """
        Test the slope array.
        """

        # Get the filenames for both the reference and test files
        ref_fname  = find_file(self.reference_dir, self.fname_slope)
        test_fname = find_file(self.test_dir, self.fname_slope)

        # Get the image data
        ref_img  = read_img(ref_fname)
        test_img = read_img(test_fname)

        # Precision
        dp = self.decimal_precision

        self.assertIsNone(npt.assert_almost_equal(test_img, ref_img,
                                                  decimal=dp))

    def test_aspect(self):
        """
        Test the aspect array.
        """

        # Get the filenames for both the reference and test files
        ref_fname  = find_file(self.reference_dir, self.fname_aspect)
        test_fname = find_file(self.test_dir, self.fname_aspect)

        # Get the image data
        ref_img  = read_img(ref_fname)
        test_img = read_img(test_fname)

        # Precision
        dp = self.decimal_precision

        self.assertIsNone(npt.assert_almost_equal(test_img, ref_img,
                                                  decimal=dp))

    def test_incident(self):
        """
        Test the incident array.
        """

        # Get the filenames for both the reference and test files
        ref_fname  = find_file(self.reference_dir, self.fname_incident)
        test_fname = find_file(self.test_dir, self.fname_incident)

        # Get the image data
        ref_img  = read_img(ref_fname)
        test_img = read_img(test_fname)

        # Precision
        dp = self.decimal_precision

        self.assertIsNone(npt.assert_almost_equal(test_img, ref_img,
                                                  decimal=dp))

    def test_exiting(self):
        """
        Test the exiting array.
        """

        # Get the filenames for both the reference and test files
        ref_fname  = find_file(self.reference_dir, self.fname_exiting)
        test_fname = find_file(self.test_dir, self.fname_exiting)

        # Get the image data
        ref_img  = read_img(ref_fname)
        test_img = read_img(test_fname)

        # Precision
        dp = self.decimal_precision

        self.assertIsNone(npt.assert_almost_equal(test_img, ref_img,
                                                  decimal=dp))

    def test_azimuth_indcident(self):
        """
        Test the azimuth incident array.
        """

        # Get the filenames for both the reference and test files
        ref_fname = find_file(self.reference_dir, self.fname_azi_indcident)
        test_fname = find_file(self.test_dir, self.fname_azi_indcident)

        # Get the image data
        ref_img = read_img(ref_fname)
        test_img = read_img(test_fname)

        # Precision
        dp = self.decimal_precision

        self.assertIsNone(npt.assert_almost_equal(test_img, ref_img,
                                                  decimal=dp))

    def test_azimuth_exiting(self):
        """
        Test the azimuth exiting array.
        """

        # Get the filenames for both the reference and test files
        ref_fname  = find_file(self.reference_dir, self.fname_azi_exiting)
        test_fname = find_file(self.test_dir, self.fname_azi_exiting)

        # Get the image data
        ref_img  = read_img(ref_fname)
        test_img = read_img(test_fname)

        # Precision
        dp = self.decimal_precision

        self.assertIsNone(npt.assert_almost_equal(test_img, ref_img,
                                                  decimal=dp))

    def test_relative_slope(self):
        """
        Test the relative slope array.
        """

        # Get the filenames for both the reference and test files
        ref_fname  = find_file(self.reference_dir, self.fname_rela_slope)
        test_fname = find_file(self.test_dir, self.fname_rela_slope)

        # Get the image data
        ref_img  = read_img(ref_fname)
        test_img = read_img(test_fname)

        # Precision
        dp = self.decimal_precision

        self.assertIsNone(npt.assert_almost_equal(test_img, ref_img,
                                                  decimal=dp))


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(description='Calculates slope, aspect, incident, exiting & relative slope.')

    parser.add_argument('--L1T_dir', required=True, help='A directory path of a L1T scene.')
    parser.add_argument('--nbar_work_dir', required=True, help='A directory path to the associated NBAR working directory.')
    parser.add_argument('--outdir', required=True, help='A directory path that will contain the output files.')
    parser.add_argument('--dec_precision', default=4, help='The decimal precision used for array comparison')
    parser.add_argument('--int_precision', default=1, help='The integer precision used for array comparison')
    parser.add_argument('--compute', action='store_true', help='If set then the slope, aspect, incident, exiting & relative angles will be computed before running the unittests.')
    parser.add_argument('--buffer', default=250, help='The buffer in pixels to be used in calculating the slope.')

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
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        # Get the current directory
        cwd = os.getcwd()

        # Change to the output directory that will contain the results
        os.chdir(outdir)

        # Open the L1T dataset
        ds = SceneDataset(L1T_dir)

        # Compute the angles
        calculate_slope_angles(ds, nbar_work_dir, outdir, buffer)

        # Close the L1T dataset
        ds = None

        # Change back to the original directory
        os.chdir(cwd)

    print "Checking that we have all the reference and test data files neccessary."
    suite = unittest.TestSuite()
    suite.addTest(ParameterisedTestCase.parameterise(TestSlopeFilenames,
                  reference_dir=nbar_work_dir, test_dir=outdir,
                  decimal_precision=dec_precision,
                  integer_precision=int_precision))
    unittest.TextTestRunner(verbosity=2).run(suite)

    print "Comparing the reference and test slope output files."
    suite = unittest.TestSuite()
    suite.addTest(ParameterisedTestCase.parameterise(TestSlopeAngles,
                  reference_dir=nbar_work_dir, test_dir=outdir,
                  decimal_precision=dec_precision,
                  integer_precision=int_precision))
    unittest.TextTestRunner(verbosity=2).run(suite)

