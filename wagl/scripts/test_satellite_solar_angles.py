#!/usr/bin/env python
"""
Unittesting framework for the `wagl.satellite_solar_angles` function.
"""

from __future__ import absolute_import, print_function, unicode_literals
import unittest
import argparse
from argparse import RawTextHelpFormatter

import numpy.testing as npt
import h5py

from wagl.hdf5 import read_h5_table
from wagl.unittesting_tools import ParameterisedTestCase


class TestCalculateAngles(ParameterisedTestCase):
    """
    Unittesting for the satellite and solar angles computation
    found in `wagl.satellite_solar_angles`.

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

    We're not explicitly testing the function, but implicitly
    testing the function by comparing against an existing
    dataset.
    """

    def test_satellite_view(self):
        """
        Test the satellite view angle array.
        """
        with h5py.File(self.reference_fname, "r") as reference_fid, h5py.File(
            self.test_fname, "r"
        ) as test_fid:

            ref_dset = reference_fid["satellite-view"]
            test_dset = test_fid["satellite-view"]

            npt.assert_almost_equal(test_dset, ref_dset, decimal=self.decimal_precision)

    def test_satellite_azimuth(self):
        """
        Test the satellite azimuth angle array.
        """
        with h5py.File(self.reference_fname, "r") as reference_fid, h5py.File(
            self.test_fname, "r"
        ) as test_fid:

            ref_dset = reference_fid["satellite-azimuth"]
            test_dset = test_fid["satellite-azimuth"]

            npt.assert_almost_equal(test_dset, ref_dset, decimal=self.decimal_precision)

    def test_solar_zenith(self):
        """
        Test the solar zenith angle array.
        """
        with h5py.File(self.reference_fname, "r") as reference_fid, h5py.File(
            self.test_fname, "r"
        ) as test_fid:

            ref_dset = reference_fid["solar-zenith"]
            test_dset = test_fid["solar-zenith"]

            npt.assert_almost_equal(test_dset, ref_dset, decimal=self.decimal_precision)

    def test_solar_azimuth(self):
        """
        Test the solar azimuth angle array.
        """
        with h5py.File(self.reference_fname, "r") as reference_fid, h5py.File(
            self.test_fname, "r"
        ) as test_fid:

            ref_dset = reference_fid["solar-azimuth"]
            test_dset = test_fid["solar-azimuth"]

            npt.assert_almost_equal(test_dset, ref_dset, decimal=self.decimal_precision)

    def test_relative_azimuth(self):
        """
        Test the relative azimuth angle array.
        """
        with h5py.File(self.reference_fname, "r") as reference_fid, h5py.File(
            self.test_fname, "r"
        ) as test_fid:

            ref_dset = reference_fid["relative-azimuth"]
            test_dset = test_fid["relative-azimuth"]

            npt.assert_almost_equal(test_dset, ref_dset, decimal=self.decimal_precision)

    def test_time_array(self):
        """
        Test the time array.
        """
        with h5py.File(self.reference_fname, "r") as reference_fid, h5py.File(
            self.test_fname, "r"
        ) as test_fid:

            ref_dset = reference_fid["acquisition-time"]
            test_dset = test_fid["acquisition-time"]

            npt.assert_almost_equal(test_dset, ref_dset, decimal=self.decimal_precision)

    def test_centreline_dataset(self):
        """
        Test the centreline dataset exactly using the
        `pandas.DataFrame` equality test.
        """
        with h5py.File(self.reference_fname, "r") as reference_fid, h5py.File(
            self.test_fname, "r"
        ) as test_fid:

            ref_data = read_h5_table(reference_fid, "centreline")
            test_data = read_h5_table(test_fid, "centreline")

            self.assertTrue(test_data.equals(ref_data))

    def test_centreline_row_index(self):
        """
        Test the `row_index` column of the centreline dataset.
        """
        with h5py.File(self.reference_fname, "r") as reference_fid, h5py.File(
            self.test_fname, "r"
        ) as test_fid:

            ref_data = reference_fid["centreline"]["row_index"]
            test_data = test_fid["centreline"]["row_index"]

            npt.assert_allclose(test_data, ref_data, atol=self.integer_precision)

    def test_centreline_col_index(self):
        """
        Test the `col_index` column of the centreline dataset.
        """
        with h5py.File(self.reference_fname, "r") as reference_fid, h5py.File(
            self.test_fname, "r"
        ) as test_fid:

            ref_data = reference_fid["centreline"]["col_index"]
            test_data = test_fid["centreline"]["col_index"]

            npt.assert_allclose(test_data, ref_data, atol=self.integer_precision)

    def test_centreline_n_pixels(self):
        """
        Test the `n_pixels` column of the centreline dataset.
        """
        with h5py.File(self.reference_fname, "r") as reference_fid, h5py.File(
            self.test_fname, "r"
        ) as test_fid:

            ref_data = reference_fid["centreline"]["n_pixels"]
            test_data = test_fid["centreline"]["n_pixels"]

            npt.assert_allclose(test_data, ref_data, atol=self.integer_precision)

    def test_centreline_latitude(self):
        """
        Test the `latitude` column of the centreline dataset.
        """
        with h5py.File(self.reference_fname, "r") as reference_fid, h5py.File(
            self.test_fname, "r"
        ) as test_fid:

            ref_data = reference_fid["centreline"]["latitude"]
            test_data = test_fid["centreline"]["latitude"]

            npt.assert_almost_equal(test_data, ref_data, decimal=self.decimal_precision)

    def test_centreline_longitude(self):
        """
        Test the `longitude` column of the centreline dataset.
        """
        with h5py.File(self.reference_fname, "r") as reference_fid, h5py.File(
            self.test_fname, "r"
        ) as test_fid:

            ref_data = reference_fid["centreline"]["longitude"]
            test_data = test_fid["centreline"]["longitude"]

            npt.assert_almost_equal(test_data, ref_data, decimal=self.decimal_precision)

    def test_boxline_dataset(self):
        """
        Test the boxline dataset exactly using the
        `pandas.DataFrame` equality test.
        """
        with h5py.File(self.reference_fname, "r") as reference_fid, h5py.File(
            self.test_fname, "r"
        ) as test_fid:

            ref_data = read_h5_table(reference_fid, "boxline")
            test_data = read_h5_table(test_fid, "boxline")

            self.assertTrue(test_data.equals(ref_data))

    def test_boxline_row_index(self):
        """
        Test the `row_index` column of the boxline dataset.
        """
        with h5py.File(self.reference_fname, "r") as reference_fid, h5py.File(
            self.test_fname, "r"
        ) as test_fid:

            ref_data = reference_fid["boxline"]["row_index"]
            test_data = test_fid["boxline"]["row_index"]

            npt.assert_allclose(test_data, ref_data, atol=self.integer_precision)

    def test_boxline_bisection_index(self):
        """
        Test the `bisection_index` column of the boxline dataset.
        """
        with h5py.File(self.reference_fname, "r") as reference_fid, h5py.File(
            self.test_fname, "r"
        ) as test_fid:

            ref_data = reference_fid["boxline"]["bisection_index"]
            test_data = test_fid["boxline"]["bisection_index"]

            npt.assert_allclose(test_data, ref_data, atol=self.integer_precision)

    def test_boxline_start_index(self):
        """
        Test the `start_index` column of the boxline dataset.
        """
        with h5py.File(self.reference_fname, "r") as reference_fid, h5py.File(
            self.test_fname, "r"
        ) as test_fid:

            ref_data = reference_fid["boxline"]["start_index"]
            test_data = test_fid["boxline"]["start_index"]

            npt.assert_allclose(test_data, ref_data, atol=self.integer_precision)

    def test_boxline_end_index(self):
        """
        Test the `end_index` column of the boxline dataset.
        """
        with h5py.File(self.reference_fname, "r") as reference_fid, h5py.File(
            self.test_fname, "r"
        ) as test_fid:

            ref_data = reference_fid["boxline"]["end_index"]
            test_data = test_fid["boxline"]["end_index"]

            npt.assert_allclose(test_data, ref_data, atol=self.integer_precision)

    def test_coordinator_dataset(self):
        """
        Test the boxline dataset exactly using the
        `pandas.DataFrame` equality test.
        """
        with h5py.File(self.reference_fname, "r") as reference_fid, h5py.File(
            self.test_fname, "r"
        ) as test_fid:

            ref_data = read_h5_table(reference_fid, "coordinator")
            test_data = read_h5_table(test_fid, "coordinator")

            self.assertTrue(test_data.equals(ref_data))

    def test_coordinator_row_index(self):
        """
        Test the `row_index` column of the boxline dataset.
        """
        with h5py.File(self.reference_fname, "r") as reference_fid, h5py.File(
            self.test_fname, "r"
        ) as test_fid:

            ref_data = reference_fid["coordinator"]["row_index"]
            test_data = test_fid["coordinator"]["row_index"]

            npt.assert_allclose(test_data, ref_data, atol=self.integer_precision)

    def test_coordinator_col_index(self):
        """
        Test the `col_index` column of the boxline dataset.
        """
        with h5py.File(self.reference_fname, "r") as reference_fid, h5py.File(
            self.test_fname, "r"
        ) as test_fid:

            ref_data = reference_fid["coordinator"]["col_index"]
            test_data = test_fid["coordinator"]["col_index"]

            npt.assert_allclose(test_data, ref_data, atol=self.integer_precision)


def _parser():
    """Argument parser."""
    description = (
        "Unittests for `wagl.satellite_solar_angles` function.\n"
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
        "\t* coordinator\n"
    )

    parser = argparse.ArgumentParser(
        description=description, formatter_class=RawTextHelpFormatter
    )
    parser.add_argument(
        "--reference_fname",
        required=True,
        help=(
            "The filename containing the reference datasets " "to be used as a baseline."
        ),
    )
    parser.add_argument(
        "--test_fname",
        required=True,
        help=(
            "The filename containing the test datasets "
            "to be used in comparing against the "
            "base/reference datasets."
        ),
    )
    parser.add_argument(
        "--decimal_precision",
        default=4,
        type=int,
        help=("The decimal precision used for the comparison " "of images."),
    )
    parser.add_argument(
        "--integer_precision",
        default=1,
        type=int,
        help=("The integer precision used for the comparison " "of images."),
    )

    return parser


def main():
    """Main execution."""
    parser = _parser()
    args = parser.parse_args()

    reference_fname = args.reference_fname
    test_fname = args.test_fname
    decimal_precision = args.decimal_precision
    integer_precision = args.integer_precision

    suite = unittest.TestSuite()
    test_case = ParameterisedTestCase()
    suite.addTest(
        test_case.parameterise(
            TestCalculateAngles,
            reference_fname=reference_fname,
            test_fname=test_fname,
            decimal_precision=decimal_precision,
            integer_precision=integer_precision,
        )
    )
    unittest.TextTestRunner(verbosity=2).run(suite)
