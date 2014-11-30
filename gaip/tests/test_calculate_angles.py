#!/usr/bin/env python

import argparse
import os
import unittest

import numpy
import numpy.testing as npt

from EOtools.DatasetDrivers import SceneDataset

from gaip import acquisitions
from gaip import calculate_angles as ca
from gaip import find_file
from gaip import read_img
from gaip import write_img
from unittesting_tools import ParameterisedTestCase


def compute_angles(L1T_path, lon_fname, lat_fname, npoints=12):
    """
    Creates the satellite and solar angle arrays as well as the time
    array.
    """
    # Retrieve an acquisitions object
    acqs = acquisitions(L1T_path)

    # Get the datetime of acquisition
    acqs[0].scene_center_datetime

    # create the geo_box
    geobox = gridded_geo_box(acqs[0])

    # Get the array dimensions
    cols = acqs[0].samples
    rows = acqs[0].lines

    # Initialise the satellite maximum view angle
    view_max = 9.0

    # Define the output file names
    sat_view_zenith_fname  = 'SAT_V.bin'
    sat_azimuth_fname      = 'SAT_AZ.bin'
    solar_zenith_fname     = 'SOL_Z.bin'
    solar_azimuth_fname    = 'SOL_AZ.bin'
    relative_azimuth_fname = 'REL_AZ.bin'
    time_fname             = 'TIME.bin'

    out_fnames = [sat_view_zenith_fname, sat_azimuth_fname, solar_zenith_fname,
                  solar_azimuth_fname, relative_azimuth_fname, time_fname]

    # Get the angles, time, & satellite track coordinates
    (satellite_zenith, satellite_azimuth, solar_zenith, 
     solar_azimuth, relative_azimuth, time,
     Y_cent, X_cent, N_cent) = ca.calculate_angles(Datetime, geobox, lon_fname,
                                   lat_fname, npoints=12, to_disk=out_fnames)

    print "Writing out the centreline file"
    # Write the centreline to disk
    outf = open('CENTRELINE', 'w')

    outf.write('%f\n'%view_max)
    outf.write('%i, %i\n'%(rows, cols))
    for r in range(rows):
        outf.write('%i, %i, %f\n'%(Y_cent[r], X_cent[r], N_cent[r]))

    outf.close()


class TestAngleFilenames(ParameterisedTestCase):
    """
    Tests to ensure we have all the correct files on disk before
    proceeding with array comparisons.
    Unittests will occur for the following files:
    SAT_V.bin
    SAT_AZ.bin
    SOL_Z.bin
    SOL_AZ.bin
    REL_AZ.bin
    TIME.bin
    CENTRELINE
    """

    # Files of interest
    ParameterisedTestCase.fname_sat_v      = 'SAT_V.bin'
    ParameterisedTestCase.fname_sat_az     = 'SAT_AZ.bin'
    ParameterisedTestCase.fname_sol_z      = 'SOL_Z.bin'
    ParameterisedTestCase.fname_sol_az     = 'SOL_AZ.bin'
    ParameterisedTestCase.fname_rel_az     = 'REL_AZ.bin'
    ParameterisedTestCase.fname_time       = 'TIME.bin'
    ParameterisedTestCase.fname_centreline = 'CENTRELINE'

    def test_centreline_ref(self):
        """
        Check that reference CENTRELINE text file exists.
        """

        fname = os.path.join(self.reference_dir, self.fname_centreline)
        self.assertIs(os.path.exists(fname), True,
                      'Reference file does not exist: %s'%fname)

    def test_centreline_tst(self):
        """
        Check that test CENTRELINE text file exists.
        """

        fname = os.path.join(self.test_dir, self.fname_centreline)
        self.assertIs(os.path.exists(fname), True, 
                      'Test file does not exist: %s'%fname)

    def test_sat_zenith_ref(self):
        """
        Check that the reference SAT_V.bin image file exists.
        """

        fname = os.path.join(self.reference_dir, self.fname_sat_v)
        self.assertIs(os.path.exists(fname), True,
                      'Reference file does not exist: %s'%fname)

    def test_sat_zenith_tst(self):
        """
        Check that the test SAT_V.bin image file exists.
        """

        fname = os.path.join(self.test_dir, self.fname_sat_v)
        self.assertIs(os.path.exists(fname), True,
                      'Test file does not exist: %s'%fname)

    def test_sat_azimuth_ref(self):
        """
        Check that the reference SAT_AZ.bin image file exists.
        """

        fname = os.path.join(self.reference_dir, self.fname_sat_az)
        self.assertIs(os.path.exists(fname), True,
                      'Reference file does not exist: %s'%fname)

    def test_sat_azimuth_tst(self):
        """
        Check that the test SAT_AZ.bin image file exists.
        """

        fname = os.path.join(self.test_dir, self.fname_sat_az)
        self.assertIs(os.path.exists(fname), True,
                      'Test file does not exist: %s'%fname)

    def test_sol_zenith_ref(self):
        """
        Check that the reference SOL_Z.bin image file exists.
        """

        fname = os.path.join(self.reference_dir, self.fname_sol_z)
        self.assertIs(os.path.exists(fname), True,
                      'Reference file does not exist: %s'%fname)

    def test_sol_zenith_tst(self):
        """
        Check that the test SOL_Z.bin image file exists.
        """

        fname = os.path.join(self.test_dir, self.fname_sol_z)
        self.assertIs(os.path.exists(fname), True,
                      'Test file does not exist: %s'%fname)

    def test_sol_azimuth_ref(self):
        """
        Check that the reference SOL_AZ.bin image file exists.
        """

        fname = os.path.join(self.reference_dir, self.fname_sol_az)
        self.assertIs(os.path.exists(fname), True,
                      'Reference file does not exist: %s'%fname)

    def test_sol_azimuth_tst(self):
        """
        Check that the test SOL_AZ.bin image file exists.
        """

        fname = os.path.join(self.test_dir, self.fname_sol_az)
        self.assertIs(os.path.exists(fname), True,
                      'Test file does not exist: %s'%fname)

    def test_rel_azimuth_ref(self):
        """
        Check that the reference REL_AZ.bin image file exists.
        """

        fname = os.path.join(self.reference_dir, self.fname_rel_az)
        self.assertIs(os.path.exists(fname), True,
                      'Reference file does not exist: %s'%fname)

    def test_rel_azimuth_tst(self):
        """
        Check that the test REL_AZ.bin image file exists.
        """

        fname = os.path.join(self.test_dir, self.fname_rel_az)
        self.assertIs(os.path.exists(fname), True,
                      'Test file does not exist: %s'%fname)

    def test_time_ref(self):
        """
        Check that the reference TIME.bin image file exists.
        """

        fname = os.path.join(self.reference_dir, self.fname_time)
        self.assertIs(os.path.exists(fname), True,
                      'Reference file does not exist: %s'%fname)

    def test_time_tst(self):
        """
        Check that the test TIME.bin image file exists.
        """

        fname = os.path.join(self.test_dir, self.fname_time)
        self.assertIs(os.path.exists(fname), True,
                      'Test file does not exist: %s'%fname)


#class AnglesOutputsTester(unittest.TestCase):
class TestSatSolAngles(ParameterisedTestCase):
    """
    Unittesting for the satellite and solar angle calculations.
    Unittests will occur for the following files:
    SAT_V.bin
    SAT_AZ.bin
    SOL_Z.bin
    SOL_AZ.bin
    REL_AZ.bin
    TIME.bin
    CENTRELINE
    """

    # Files of interest
    ParameterisedTestCase.fname_sat_v      = 'SAT_V.bin'
    ParameterisedTestCase.fname_sat_az     = 'SAT_AZ.bin'
    ParameterisedTestCase.fname_sol_z      = 'SOL_Z.bin'
    ParameterisedTestCase.fname_sol_az     = 'SOL_AZ.bin'
    ParameterisedTestCase.fname_rel_az     = 'REL_AZ.bin'
    ParameterisedTestCase.fname_time       = 'TIME.bin'
    ParameterisedTestCase.fname_centreline = 'CENTRELINE'

    # Read and store the centreline data in memory
    ParameterisedTestCase.centreline_ref = None
    ParameterisedTestCase.centreline_test = None

    def _read_centreline_files(self):
        """
        Read the centreline reference and test text files.
        """

        # Get the filenames for both the reference and test files
        ref_fname  = find_file(self.reference_dir, self.fname_centreline)
        test_fname = find_file(self.test_dir, self.fname_centreline)

        # Open and read the reference data
        f        = open(ref_fname)
        ref_data = f.readlines()
        f.close()

        # Open and read the test data
        f         = open(test_fname)
        test_data = f.readlines()
        f.close()

        self.centreline_ref  = ref_data
        self.centreline_test = test_data

    def test_satellite_view(self):
        """
        Test the satellite view angle array.
        """

        # Get the filenames for both the reference and test files
        ref_fname  = find_file(self.reference_dir, self.fname_sat_v)
        test_fname = find_file(self.test_dir, self.fname_sat_v)

        # Get the image data
        ref_img  = read_img(ref_fname)
        test_img = read_img(test_fname)

        # Precision
        dp = self.decimal_precision

        self.assertIsNone(npt.assert_almost_equal(test_img, ref_img,
                                                  decimal=dp))

    def test_satellite_azimuth(self):
        """
        Test the satellite azimuth angle array.
        """

        # Get the filenames for both the reference and test files
        ref_fname  = find_file(self.reference_dir, self.fname_sat_az)
        test_fname = find_file(self.test_dir, self.fname_sat_az)

        # Get the image data
        ref_img  = read_img(ref_fname)
        test_img = read_img(test_fname)

        # Precision
        dp = self.decimal_precision

        self.assertIsNone(npt.assert_almost_equal(test_img, ref_img,
                                                  decimal=dp))

    def test_solar_zenith(self):
        """
        Test the solar zenith angle array.
        """

        # Get the filenames for both the reference and test files
        ref_fname  = find_file(self.reference_dir, self.fname_sol_z)
        test_fname = find_file(self.test_dir, self.fname_sol_z)

        # Get the image data
        ref_img  = read_img(ref_fname)
        test_img = read_img(test_fname)

        # Precision
        dp = self.decimal_precision

        self.assertIsNone(npt.assert_almost_equal(test_img, ref_img,
                                                  decimal=dp))

    def test_solar_azimuth(self):
        """
        Test the solar azimuth angle array.
        """

        # Get the filenames for both the reference and test files
        ref_fname  = find_file(self.reference_dir, self.fname_sol_az)
        test_fname = find_file(self.test_dir, self.fname_sol_az)

        # Get the image data
        ref_img  = read_img(ref_fname)
        test_img = read_img(test_fname)

        # Precision
        dp = self.decimal_precision

        self.assertIsNone(npt.assert_almost_equal(test_img, ref_img,
                                                  decimal=dp))

    def test_relative_azimuth(self):
        """
        Test the relative azimuth angle array.
        """

        # Get the filenames for both the reference and test files
        ref_fname  = find_file(self.reference_dir, self.fname_rel_az)
        test_fname = find_file(self.test_dir, self.fname_rel_az)

        # Get the image data
        ref_img  = read_img(ref_fname)
        test_img = read_img(test_fname)

        # Precision
        dp = self.decimal_precision

        self.assertIsNone(npt.assert_almost_equal(test_img, ref_img,
                                                  decimal=dp))

    def test_time_array(self):
        """
        Test the time array.
        """

        # Get the filenames for both the reference and test files
        ref_fname  = find_file(self.reference_dir, self.fname_time)
        test_fname = find_file(self.test_dir, self.fname_time)

        # Get the image data
        ref_img  = read_img(ref_fname)
        test_img = read_img(test_fname)

        # Precision
        dp = self.decimal_precision

        self.assertIsNone(npt.assert_almost_equal(test_img, ref_img,
                                                  decimal=dp))

    def test_centreline_max_view_angle(self):
        """
        Test the maximum satellite view angle.
        This is the first line of the centreline file.
        """

        # Read the CENTRELINE file if needed
        if self.centreline_ref is None:
            self._read_centreline_files()

        ref_data  = float(self.centreline_ref[0].split()[0])
        test_data = float(self.centreline_test[0].split()[0])

        self.assertEqual(ref_data, test_data)

    def test_centreline_rows(self):
        """
        Test the number of rows.
        This is the first element of the second line in the
        centreline file.
        """

        # Read the CENTRELINE file if needed
        if self.centreline_ref is None:
            self._read_centreline_files()

        ref_rows  = int(self.centreline_ref[1].split()[0])
        test_rows = int(self.centreline_test[1].split()[0])

        self.assertEqual(ref_rows, test_rows)

    def test_centreline_columns(self):
        """
        Test the number of columns.
        This is the second element of the second line in the
        centreline file.
        """

        # Read the CENTRELINE file if needed
        if self.centreline_ref is None:
            self._read_centreline_files()

        ref_cols  = int(self.centreline_ref[1].split()[1])
        test_cols = int(self.centreline_test[1].split()[1])

        self.assertEqual(ref_cols, test_cols)

    def test_number_centreline_points(self):
        """
        Test that the number of centreline points are same.
        These start at the third line of the centreline file and
        contain 3 elements.
        """

        # Read the CENTRELINE file if needed
        if self.centreline_ref is None:
            self._read_centreline_files()

        ref_data  = self.centreline_ref[2:]
        test_data = self.centreline_test[2:]

        self.assertEqual(len(ref_data), len(test_data))

    def test_centreline_points(self):
        """
        Test that the centreline points are roughly the same.
        These start at the third line of the centreline file and
        contain 3 elements. We'll allow for variation by 1 integer.

        eg ['1742','4624','1.00000']
        """

        # Read the CENTRELINE file if needed
        if self.centreline_ref is None:
            self._read_centreline_files()

        ref_data  = self.centreline_ref[2:]
        test_data = self.centreline_test[2:]

        int_prec = self.integer_precision

        ref_points  = numpy.zeros((len(ref_data), 3), dtype='int')
        test_points = numpy.zeros((len(ref_data), 3))

        for i in range(len(ref_data)):
            rx, ry, rz = ref_data[i].split()
            tx, ty, tz = test_data[i].split()
            ref_points[i,0], ref_points[i,1], ref_points[i,2] = rx, ry, float(rz)
            test_points[i,0], test_points[i,1], test_points[i,2] = tx, ty, float(tz)

        self.assertIsNone(npt.assert_allclose(test_points, ref_points, atol=int_prec))

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(description='Calculates satellite and solar angles.')

    parser.add_argument('--L1T_dir', required=True, help='A directory path of a L1T scene.')
    parser.add_argument('--nbar_work_dir', required=True, help='A directory path to the associated NBAR working directory.')
    parser.add_argument('--outdir', required=True, help='A directory path that will contain the output files.')
    parser.add_argument('--dec_precision', default=4, help='The decimal precision used for array comparison')
    parser.add_argument('--int_precision', default=1, help='The integer precision used for array comparison')
    parser.add_argument('--compute', action='store_true', help='If set then the solar and sateliite angles will be computed as will the CENTRELINE text file before running the unittests.')

    parsed_args = parser.parse_args()

    #f = '/g/data1/v10/NBAR_validation_reference/Nov2013/L1T_Input/LS7_90-84_2013-10-03/EQR/LS7_ETM_OTH_P51_GALPGS04-002_090_084_20131003'
    #lon_f = '/short/v10/jps547/nbar/test_skew/work/LON.tif'
    #lat_f = '/short/v10/jps547/nbar/test_skew/work/LAT.tif'

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

        # Find and open the longitude and lattitude files
        lon_fname = find_file(nbar_work_dir, 'LON.tif')
        lat_fname = find_file(nbar_work_dir, 'LAT.tif')

        print "Computing satellite & solar angle grids."
        compute_angles(L1T_dir, lon_fname, lat_fname)

        # Change back to the original directory
        os.chdir(cwd)

    print "Checking that we have all the reference and test data files neccessary."
    suite = unittest.TestSuite()
    suite.addTest(ParameterisedTestCase.parameterise(TestAngleFilenames,
                  reference_dir=nbar_work_dir, test_dir=outdir,
                  decimal_precision=dec_precision,
                  integer_precision=int_precision))
    unittest.TextTestRunner(verbosity=2).run(suite)

    print "Comparing the reference and test angle grids and CENTRELINE file."
    suite = unittest.TestSuite()
    suite.addTest(ParameterisedTestCase.parameterise(TestSatSolAngles,
                  reference_dir=nbar_work_dir, test_dir=outdir,
                  decimal_precision=dec_precision,
                  integer_precision=int_precision))
    unittest.TextTestRunner(verbosity=2).run(suite)

