#!/usr/bin/env python

import argparse
import os
from os.path import join
from os.path import dirname
import unittest

import numpy as np
import numpy.testing as npt
from osgeo import gdal
from osgeo import gdalconst

from EOtools.DatasetDrivers import SceneDataset
from ULA3.filtering import filter_float as filter
from ULA3.filtering import read_array
from unittesting_tools import find_file
from unittesting_tools import ParameterisedTestCase
from unittesting_tools import read_img
from unittesting_tools import write_img


"""
Since we are lacking a portable test harness, these tests are configured to run on the NCI.

Hopefully, in the near future, we will have a test harness that can be more portable, and testing will be
possible elsewhere.
"""

class FilterTestCase(unittest.TestCase):
    """
    Test case for running Fuquin's filter.
    """
    def setUp(self):
        self.eps = 1e-4
        pass

    def tearDown(self):
        pass

    def runTest(self):
        #input_dataset = gdal.Open(
        #    join(dirname(__file__), "filter_data", "test_filter_input.img"),
        #    gdalconst.GA_ReadOnly)
        #filtered_data = filter(input_dataset.ReadAsArray().astype(np.float32))
        #print '(' + str(input_dataset.RasterYSize), str(input_dataset.RasterXSize) + ')'
        #input_dataset = None

        nrow = 153
        ncol = 131
        input_data = read_array(join(dirname(__file__), "filter_data", "test_filter_input.img"), nrow, ncol)
        filtered_data = filter(input_data)

        comparison_dataset = gdal.Open(
            join(dirname(__file__), "filter_data", "test_filter_target_output.img"),
            gdalconst.GA_ReadOnly)
        comparison_data = comparison_dataset.ReadAsArray().astype(np.float32)
        print '(' + str(comparison_dataset.RasterYSize), str(comparison_dataset.RasterXSize) + ')'
        comparison_dataset = None

        differences = filtered_data - comparison_data

        #for fd, cd in zip(filtered_data, comparison_data):
        #    print ', '.join([':'.join([str(x) for x in x]) for x in zip(fd, cd)])

        #for d in differences:
        #    print ', '.join([str(x) for x in d])

        self.assertTrue((differences < self.eps).all())





def suite_old():
    return unittest.TestSuite((FilterTestCase(),))


def calculate_smoothed_dsm(scene_dataset, ref_dir, outdir):
    """
    Creates a smoothed DSM.
    """

    # Check and load the required files from disk
    fname_dsm = 'region_dsm_image.img'

    dsm = (read_img(find_file(ref_dir, fname_dsm))).astype('float32')

    smoothed_dsm = filter(dsm)

    # Image projection, geotransform
    prj = scene_dataset.GetProjection()
    geoT = scene_dataset.GetGeoTransform()

    # Write out the smoothed dsm file
    out_fname = os.path.join(outdir, 'region_dsm_image_smoothed.img')
    write_img(smoothed_dsm, out_fname, projection=prj, geotransform=geoT)


class TestFilterFileNames(ParameterisedTestCase):
    """
    Unittests will occur for the following files:
    region_dsm_image_smoothed.img
    """

    def __init__(self):
        """
        """

        self.fname_smoothed_dsm = 'region_dsm_image_smoothed.img'

    def test_smoothed_dsm_ref(self):
        """
        Check that the smoothed dsm reference file exists.
        """
        fname = os.path.join(self.reference_dir, self.fname_smoothed_dsm)
        self.assertIs(os.path.exists(fname), True,
                      '_reference file does not exist: %s'%fname)

    def test_smoothed_dsm_tst(self):
        """
        Check that the smoothed dsm test file exists.
        """
        fname = os.path.join(self.test_dir, self.fname_smoothed_dsm)
        self.assertIs(os.path.exists(fname), True,
                      '_reference file does not exist: %s'%fname)


class TestFilterOutputs(ParameterisedTestCase):
    """
    Unittests will occur for the following files:
    region_dsm_image_smoothed.img
    """

    def __init__(self):
        """
        """

        self.fname_smoothed_dsm = 'region_dsm_image_smoothed.img'

    def test_smoothed_dsm(self):
        """
        Test the smoothed dsm image against the reference image.
        """

        # Get the filenames for both the reference and test files
        ref_fname = find_file(self.reference_dir, self.fname_smoothed_dsm)
        test_fname = find_file(self.test_dir, self.fname_smoothed_dsm)

        # Get the image data
        ref_img = read_img(ref_fname)
        test_img = read_img(test_fname)

        self.assertIsNone(npt.assert_almost_equal(test_img, ref_img,
                                                  decimal=self.dec_precision))


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(description='Performs unittests against the smoothed DSM and optionally calculates the smoothed DSM.')

    parser.add_argument('--L1T_dir', required=True, help='A directory path of a L1T scene.')
    parser.add_argument('--nbar_work_dir', required=True, help='A directory path to the associated NBAR working directory.')
    parser.add_argument('--outdir', required=True, help='A directory path that will contain the output files.')
    parser.add_argument('--dec_precision', default=4, help='The decimal precision used for array comparison')
    parser.add_argument('--int_precision', default=1, help='The integer precision used for array comparison')
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
        ds = SceneDataset(L1T_dir)

        # Compute the smoothed dsm
        calculate_smoothed_dsm(ds, nbar_work_dir, outdir)

        # Close the L1T dataset
        ds = None

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

