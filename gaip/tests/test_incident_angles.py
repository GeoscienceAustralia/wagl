#!/usr/bin/env python

import argparse
import os
from os.path import join as pjoin
from os.path import exists as pexists
import unittests

import numpy.testing as npt

from gaip import acquisitions
from gaip import incident_angles
from gaip import find_file
from gaip import read_img
from gaip.tests.unittesting_tools import ParameterisedTestCase

#TODO Filename to be determined from the nbar.cfg file

def compute_incident_angles(ref_dir, out_dir):
    """
    
    """
    # TC_Intermediates directory
    tc_dir = pjoin(ref_dir, 'TC_Intermediates')
    tc_outdir = pjoin(out_dir, 'TC_Intermediates')

    # Check and load the required files from disk
    fname_solar_zenith  = find_file(ref_dir, 'SOLAR_ZENITH.bin')
    fname_solar_azimuth = find_file(ref_dir, 'SOLAR_AZIMUTH.bin')
    fname_slope = find_file(tc_dir, 'slope.bin')
    fname_aspect = find_file(tc_dir, 'aspect.bin')

    # Output filenames
    incident_out_fname = pjoin(tc_outdir, 'incident_angle.bin')
    azimuth_incident_out_fname = pjoin(tc_outdir, 'azimuth_incident_angle.bin')
    
    incident_angles(solar_zenith_fname, solar_azimuth_fname, slope_fname,
                    aspect_fname, incident_out_fname,
                    azimuth_incident_out_fname)


class TestIncidentAngleFileNames(ParameterisedTestCase):

    """
    Unittests will occur for the following files:
        * incident_angle.bin
        * azimuth_incident_angle.bin
    """

    ParameterisedTestCase.fname_incident_angle = 'incident_angle.bin'
    ParameterisedTestCase.fname_azi_inci_angle = 'azimuth_incident_angle.bin'


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


class TestIncidentAngleOutputs(ParameterisedTestCase):

    """
    Unittests will occur for the following files:
        * incident_angle.bin
        * azimuth_incident_angle.bin
    """

    ParameterisedTestCase.fname_incident_angle = 'incident_angle.bin'
    ParameterisedTestCase.fname_azi_inci_angle = 'azimuth_incident_angle.bin'


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


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(description=('Perform unittesting for '
                                                  'the incident angles.'))

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

    parsed_args = parser.parse_args()

    L1T_dir       = parsed_args.L1T_dir
    nbar_work_dir = parsed_args.nbar_work_dir
    outdir        = parsed_args.outdir
    dec_precision = parsed_args.dec_precision
    int_precision = parsed_args.int_precision
    compute       = parsed_args.compute

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
        compute_incident_angles(nbar_work_dir, outdir)

        # Close the L1T dataset
        acqs = None

        # Change back to the original directory
        os.chdir(cwd)

    print "Checking that we have reference and test data files neccessary."
    suite = unittest.TestSuite()
    suite.addTest(ParameterisedTestCase.parameterise(
                  TestIncidentAngleFileNames,
                  reference_dir=nbar_work_dir, test_dir=outdir,
                  decimal_precision=dec_precision,
                  integer_precision=int_precision))
    unittest.TextTestRunner(verbosity=2).run(suite)

    print "Comparing the reference and test incident angle outputs."
    suite = unittest.TestSuite()
    suite.addTest(ParameterisedTestCase.parameterise(
                  TestIncidentAngleOutputs,
                  reference_dir=nbar_work_dir, test_dir=outdir,
                  decimal_precision=dec_precision,
                  integer_precision=int_precision))
    unittest.TextTestRunner(verbosity=2).run(suite)
