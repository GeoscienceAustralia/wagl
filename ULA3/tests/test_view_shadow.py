#!/usr/bin/env python

import argparse
import os
import unittest

import numpy.testing as npt

from EOtools.DatasetDrivers import SceneDataset
from ULA3._gdal_tools import Buffers
from ULA3.tc import run_castshadow
from unittesting_tools import find_file
from unittesting_tools import ParameterisedTestCase
from unittesting_tools import read_img
from unittesting_tools import write_img

def calculate_view_shadow(scene_dataset, ref_dir, outdir, pixel_buffer=250,
                          block_height=500, block_width=500):
    """
    Calculates the view shadow array.
    """

    # Image projection, geotransform, UTM(True/False)
    prj    = scene_dataset.GetProjection()
    geoT   = scene_dataset.GetGeoTransform()
    is_utm =  not scene_dataset.IsGeographic()

    # Define the pixel buffering object
    pixel_buf = Buffers(pixel_buffer)

    # Check and load the required files from disk
    fname_satellite_zenith  = 'SAT_V.bin'
    fname_satellte_azimuth  = 'SAT_AZ.bin'
    fname_smoothed_dsm      = 'region_dsm_image_smoothed.img'

    zen_angle = read_img(find_file(ref_dir, fname_satellite_zenith))
    azi_angle = read_img(find_file(ref_dir, fname_satellte_azimuth))

    dsm = (read_img(find_file(ref_dir, fname_smoothed_dsm))).astype('float32')

    # Compute the self shadow
    shadow_self = run_castshadow(scene_dataset, dsm, zen_angle, azi_angle,
                                 pixel_buf, block_height, block_width, is_utm)

    # Write the self shadow result to disk
    outfname = os.path.join(outdir, 'shadow_v.img')
    write_img(shadow_self, outfname, projection=prj, geotransform=geoT)


class TestViewShadowFileNames(ParameterisedTestCase):
    """
    Unittests will occur for the following files:
    shadow_v.img
    """

    ParameterisedTestCase.fname_view_shadow = 'shadow_v.img'

    def test_view_shadow_ref(self):
        """
        Check that the view shadow reference file exists.
        """

        fname = os.path.join(self.reference_dir, self.fname_view_shadow)
        self.assertIs(os.path.exists(fname), True,
                      'Reference file does not exist: %s'%fname)

    def test_view_shadow_tst(self):
        """
        Check that the view shadow test file exists.
        """

        fname = os.path.join(self.test_dir, self.fname_view_shadow)
        self.assertIs(os.path.exists(fname), True,
                      'Test file does not exist: %s'%fname)


class TestViewShadowOutputs(ParameterisedTestCase):
    """
    Unittests will occur for the following files:
    shadow_v.img
    """

    ParameterisedTestCase.fname_view_shadow = 'shadow_v.img'

    def test_view_shadow(self):
        """
        Test the view shadow image against the reference image.
        """

        # Get the filenames for both the reference and test files
        ref_fname  = find_file(self.reference_dir, self.fname_view_shadow)
        test_fname = find_file(self.test_dir, self.fname_view_shadow)

        # Get the image data
        ref_img  = read_img(ref_fname)
        test_img = read_img(test_fname)

        # Precision
        dp = self.dec_precision

        self.assertIsNone(npt.assert_almost_equal(test_img, ref_img,
                                                  decimal=dp))


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(description='Perform unittesting for the view shadow array and optionaly calculates the view shadow array.')

    parser.add_argument('--L1T_dir', required=True, help='A directory path of a L1T scene.')
    parser.add_argument('--nbar_work_dir', required=True, help='A directory path to the associated NBAR working directory.')
    parser.add_argument('--outdir', required=True, help='A directory path that will contain the output files.')
    parser.add_argument('--dec_precision', default=4, help='The decimal precision used for array comparison')
    parser.add_argument('--int_precision', default=1, help='The integer precision used for array comparison')
    parser.add_argument('--compute', action='store_true', help='If set then the view shadow array will be computed before running the unittests.')
    parser.add_argument('--buffer', default=250, help='The buffer in pixels to be used in calculating the view shadow.')
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
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        # Get the current directory
        cwd = os.getcwd()

        # Change to the output directory that will contain the results
        os.chdir(outdir)

        # Open the L1T dataset
        ds = SceneDataset(L1T_dir)

        # Compute the angles
        calculate_view_shadow(ds, nbar_work_dir, outdir, buffer,
                              block_y, block_x)

        # Close the L1T dataset
        ds = None

        # Change back to the original directory
        os.chdir(cwd)

    print "Checking that we have all the reference and test data files neccessary."
    suite = unittest.TestSuite()
    suite.addTest(ParameterisedTestCase.parameterise(TestViewShadowFileNames,
                  reference_dir=nbar_work_dir, test_dir=outdir,
                  decimal_precision=dec_precision,
                  integer_precision=int_precision))
    unittest.TextTestRunner(verbosity=2).run(suite)

    print "Comparing the reference and test view shadow output files."
    suite = unittest.TestSuite()
    suite.addTest(ParameterisedTestCase.parameterise(TestViewShadowOutputs,
                  reference_dir=nbar_work_dir, test_dir=outdir,
                  decimal_precision=dec_precision,
                  integer_precision=int_precision))
    unittest.TextTestRunner(verbosity=2).run(suite)

