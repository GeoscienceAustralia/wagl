#!/usr/bin/env python

import argparse
import os
import unittest

import numpy
import numpy.testing as npt
import h5py

from gaip import read_table
from gaip.tests.unittesting_tools import ParameterisedTestCase


class TestCalculateAngles(ParameterisedTestCase):
    """
    Unittesting for the satellite and solar angles computation
    found in `gaip.calculate_angles`.

    Unittests will occur for the following datasets:

        * satellite-view
        * satellite-azimuth
        * solar-zenith
        * solar-azimuth
        * relative-azimuth
        * acquisition-time
        * centreline
        * boxline
        * coordinator
    """

    reference_fid = h5py.File(self.reference_fname, 'r')
    test_fid = h5py.File(self.test_fname, 'r')

    def test_satellite_view(self):
        """
        Test the satellite view angle array.
        """
        ref_dset = self.reference_fid['satellite-view']
        test_dset = self.test_fid['satellite-view']

        npt.assert_almost_equal(test_dset, ref_dset,
                                decimal=self.decimal_precision)

    def test_satellite_azimuth(self):
        """
        Test the satellite azimuth angle array.
        """
        ref_dset = self.reference_fid['satellite-azimuth']
        test_dset = self.test_fid['satellite-azimuth']

        npt.assert_almost_equal(test_dset, ref_dset,
                                decimal=self.decimal_precision)

    def test_solar_zenith(self):
        """
        Test the solar zenith angle array.
        """
        ref_dset = self.reference_fid['solar-zenith']
        test_dset = self.test_fid['solar-zenith']

        npt.assert_almost_equal(test_dset, ref_dset,
                                decimal=self.decimal_precision)

    def test_solar_azimuth(self):
        """
        Test the solar azimuth angle array.
        """
        ref_dset = self.reference_fid['solar-azimuth']
        test_dset = self.test_fid['solar-azimuth']

        npt.assert_almost_equal(test_dset, ref_dset,
                                decimal=self.decimal_precision)

    def test_relative_azimuth(self):
        """
        Test the relative azimuth angle array.
        """
        ref_dset = self.reference_fid['relative-azimuth']
        test_dset = self.test_fid['relative-azimuth']

        npt.assert_almost_equal(test_dset, ref_dset,
                                decimal=self.decimal_precision)

    def test_time_array(self):
        """
        Test the time array.
        """
        ref_dset = self.reference_fid['acquisition-time']
        test_dset = self.test_fid['acquisition-time']

        npt.assert_almost_equal(test_dset, ref_dset,
                                decimal=self.decimal_precision)

    def test_centreline_points(self):
        """
        Test that the centreline points are roughly the same.
        These start at the third line of the centreline file and
        contain 5 elements. We'll allow for variation by 1 integer.

        eg ['1742','4624','1.00000', '-33.608469649266', '150.080204768921']
        We'll only test the 1st three elements.
        """
        # TODO:
        # ref_df = read_table(self.reference_fid, 'centreline')
        # test_df = read_table(self.test_fid, 'centreline')
        # test_df.equals(ref_df) ???

        # self.assertIsNone(npt.assert_allclose(test_points, ref_points,
        #     atol=int_prec))


    def test_centreline_lonlat_points(self):
        """
        Test that the centreline longitude and latitude points are
        roughly the same.
        These start at the third line of the centreline file and
        contain 5 elements. We'll allow for variation at the 6th
        decimal place.

        eg ['1742','4624','1.00000', '-33.608469649266', '150.080204768921']
        We'll only test the last two elements.
        """
        # self.assertIsNone(npt.assert_almost_equal(test_points, ref_points,
        #     decimal=6))

    # TODO
    # def test_boxline(self):
    # def test_coordinator(self):

if __name__ == '__main__':
    description = ("Unittests for `gaip.calculate_angles` function.\n"
                   "Comparisons tests will occur for the following "
                   "datasets: \n"
                   "\t* satellite-view-angles\n"
                   "\t* satellite-azimuth-angles\n"
                   "\t* solar-zenith-angles\n"
                   "\t* solar-azimuth-angles\n"
                   "\t* relative-azimuth-angles\n"
                   "\t* acquisition-times\n"
                   "\t* centreline\n"
                   "\t* boxline\n"
                   "\t* coordinator\n")
                   
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--reference_fname', requried=True,
                        help=('The filename containing the reference datasets '
                              'to be used as a baseline.'))
    parser.add_argument('--test_fname', require=True,
                        help=('The filename containing the test datasets '
                              'to be used in comparing against the '
                              'base/reference datasets.'))
    parser.add_argument('--decimal_precision', default=4, type=int,
                        help=('The decimal precision used for the comparison '
                              'of images.'))
    parser.add_argument('--integer_precision', default=1, type=int,
                        help=('The integer precision used for the comparison '
                              'of images.'))

    parsed_args = parser.parse_args()

    reference_fname = parsed_args.reference_fname
    test_fname = parsed_args.test_fname
    decimal_precision = parsed_args.decimal_precision
    integer_precision = parsed_args.integer_precision

    suite = unittest.TestSuite()
    suite.addTest(ParameterisedTestCase.parameterise(TestCalculateAngles,
                  reference_fname=reference_fname, test_fname=test_fname,
                  decimal_precision=decimal_precision,
                  integer_precision=integer_precision))
    unittest.TextTestRunner(verbosity=2).run(suite)
