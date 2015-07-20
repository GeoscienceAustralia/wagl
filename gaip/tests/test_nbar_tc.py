#!/usr/bin/env python

import argparse
import glob
import os
from os.path import join as pjoin, abspath
import unittest

import numpy.testing as npt

from gaip import read_img
from gaip.tests.unittesting_tools import ParameterisedTestCase


class TestProductFileNames(ParameterisedTestCase):

    """
    Tests that the reference and the current datasets
    have matching files.
    """

    def test_lambert_files(self):
        """
        Check that the same number of reference and test files
        exist for the lambertian output.
        """
        reflectance_ref_dir = pjoin(self.reference_dir, 'Reflectance_Outputs')
        reflectance_test_dir = pjoin(self.test_dir, 'Reflectance_Outputs')

        cwd = os.getcwd()

        # Get the reference files
        os.chdir(reflectance_ref_dir)
        files = glob.glob('*.bin')
        ref_files = [abspath(f) for f in files if 'lambertian' in f]

        # Get the test files
        os.chdir(reflectance_test_dir)
        files = glob.glob('*.bin')
        test_files = [abspath(f) for f in files if 'lambertian' in f]

        # Change back to the original directory
        os.chdir(cwd)

        self.assertEqual(len(ref_files), len(test_files))


    def test_brdf_files(self):
        """
        Check that the same number of reference and test files
        exist for the brdf output.
        """
        reflectance_ref_dir = pjoin(self.reference_dir, 'Reflectance_Outputs')
        reflectance_test_dir = pjoin(self.test_dir, 'Reflectance_Outputs')

        cwd = os.getcwd()

        # Get the reference files
        os.chdir(reflectance_ref_dir)
        files = glob.glob('*.bin')
        ref_files = [abspath(f) for f in files if 'brdf' in f]

        # Get the test files
        os.chdir(reflectance_test_dir)
        files = glob.glob('*.bin')
        test_files = [abspath(f) for f in files if 'brdf' in f]

        # Change back to the original directory
        os.chdir(cwd)

        self.assertEqual(len(ref_files), len(test_files))


    def test_tc_files(self):
        """
        Check that the same number of reference and test files
        exist for the tc output.
        """
        reflectance_ref_dir = pjoin(self.reference_dir, 'Reflectance_Outputs')
        reflectance_test_dir = pjoin(self.test_dir, 'Reflectance_Outputs')

        cwd = os.getcwd()

        # Get the reference files
        os.chdir(reflectance_ref_dir)
        files = glob.glob('*.bin')
        ref_files = [abspath(f) for f in files if 'terrain' in f]

        # Get the test files
        os.chdir(reflectance_test_dir)
        files = glob.glob('*.bin')
        test_files = [abspath(f) for f in files if 'terrain' in f]

        # Change back to the original directory
        os.chdir(cwd)

        self.assertEqual(len(ref_files), len(test_files))


class TestProductDifference(ParameterisedTestCase):

    """
    Test the numerical difference of the product outputs.
    """

    def test_compare_lmbrt_files(self):
        """
        Check that the lambertian outputs are roughly equal.
        """
        reflectance_ref_dir = pjoin(self.reference_dir, 'Reflectance_Outputs')
        reflectance_test_dir = pjoin(self.test_dir, 'Reflectance_Outputs')

        cwd = os.getcwd()

        # Get the reference files
        os.chdir(reflectance_ref_dir)
        files = glob.glob('*.bin')
        ref_files = [abspath(f) for f in files if 'lambertian' in f]

        # Get the test files
        os.chdir(reflectance_test_dir)
        files = glob.glob('*.bin')
        test_files = [abspath(f) for f in files if 'lambertian' in f]

        # Change back to the original directory
        os.chdir(cwd)

        ref_files.sort()
        test_files.sort()

        # Precision
        ip = self.integer_precision

        for i in range(len(ref_files)):
            ref_fname = ref_files[i]
            test_fname = test_files[i]
            print "Testing:\n {}\n{}".format(ref_fname, test_fname)

            # Get the reference data
            ref_img = read_img(ref_fname)

            # Get the test data
            test_img = read_img(test_fname)

            self.assertIsNone(npt.assert_allclose(test_img, ref_img, atol=ip))


    def test_compare_brdf_files(self):
        """
        Check that the brdf outputs are roughly equal.
        """
        reflectance_ref_dir = pjoin(self.reference_dir, 'Reflectance_Outputs')
        reflectance_test_dir = pjoin(self.test_dir, 'Reflectance_Outputs')

        cwd = os.getcwd()

        # Get the reference files
        os.chdir(reflectance_ref_dir)
        files = glob.glob('*.bin')
        ref_files = [abspath(f) for f in files if 'brdf' in f]

        # Get the test files
        os.chdir(reflectance_test_dir)
        files = glob.glob('*.bin')
        test_files = [abspath(f) for f in files if 'brdf' in f]

        # Change back to the original directory
        os.chdir(cwd)

        ref_files.sort()
        test_files.sort()

        # Precision
        ip = self.integer_precision

        for i in range(len(ref_files)):
            ref_fname = ref_files[i]
            test_fname = test_files[i]
            print "Testing:\n {}\n{}".format(ref_fname, test_fname)

            # Get the reference data
            ref_img = read_img(ref_fname)

            # Get the test data
            test_img = read_img(test_fname)

            self.assertIsNone(npt.assert_allclose(test_img, ref_img, atol=ip))


    def test_compare_tc_files(self):
        """
        Check that the tc outputs are roughly equal.
        """
        reflectance_ref_dir = pjoin(self.reference_dir, 'Reflectance_Outputs')
        reflectance_test_dir = pjoin(self.test_dir, 'Reflectance_Outputs')

        cwd = os.getcwd()

        # Get the reference files
        os.chdir(reflectance_ref_dir)
        files = glob.glob('*.bin')
        ref_files = [abspath(f) for f in files if 'terrain' in f]

        # Get the test files
        os.chdir(reflectance_test_dir)
        files = glob.glob('*.bin')
        test_files = [abspath(f) for f in files if 'terrain' in f]

        # Change back to the original directory
        os.chdir(cwd)

        ref_files.sort()
        test_files.sort()

        # Precision
        ip = self.integer_precision

        for i in range(len(ref_files)):
            ref_fname = ref_files[i]
            test_fname = test_files[i]
            print "Testing:\n {}\n{}".format(ref_fname, test_fname)

            # Get the reference data
            ref_img = read_img(ref_fname)

            # Get the test data
            test_img = read_img(test_fname)

            self.assertIsNone(npt.assert_allclose(test_img, ref_img, atol=ip))


if __name__ == '__main__':

    desc = 'Test and compare the Lambertian, BRDF and TC outputs.'
    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument('--reference_dir', required=True,
                        help='A directory path of a nbar reference output.')
    parser.add_argument('--test_dir', required=True,
                        help=('A directory path to the test output.'))
    parser.add_argument('--int_precision', default=1, type=int,
                        help='The integer precision used for array comparison')

    parsed_args = parser.parse_args()

    reference_dir = parsed_args.reference_dir
    test_dir = parsed_args.test_dir
    int_precision = parsed_args.int_precision


    print "\nChecking that we have all the reference and test data files."
    suite = unittest.TestSuite()
    suite.addTest(ParameterisedTestCase.parameterise(TestProductFileNames,
                  reference_dir=reference_dir, test_dir=test_dir,
                  integer_precision=int_precision))
    unittest.TextTestRunner(verbosity=2).run(suite)

    print "\nTesting the numerical precision on each product output."
    suite = unittest.TestSuite()
    suite.addTest(ParameterisedTestCase.parameterise(TestProductDifference,
                  reference_dir=reference_dir, test_dir=test_dir,
                  integer_precision=int_precision))
    unittest.TextTestRunner(verbosity=2).run(suite)
