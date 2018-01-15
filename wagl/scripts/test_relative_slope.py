#!/usr/bin/env python
"""
Unittesting framework for the `wagl.relative_azimuth_slope` function.
"""

from __future__ import absolute_import, print_function, unicode_literals
import unittest
import argparse
from argparse import RawTextHelpFormatter

import numpy.testing as npt
import h5py

from wagl.unittesting_tools import ParameterisedTestCase


class TestRelativeSlope(ParameterisedTestCase):
    """
    Unittesting for the relative azimuth slope computation
    found in `wagl.relative_azimuth_slope`.

    Unittests will occur for the following datasets:

        * relative-slope

    We're not explicitly testing the function, but implicitly
    testing the function by comparing against an existing
    dataset.
    """

    def test_relative_slope(self):
        """
        Test the relative slope array.
        """
        with h5py.File(self.reference_fname, 'r') as reference_fid,\
            h5py.File(self.test_fname, 'r') as test_fid:

            ref_dset = reference_fid['relative-slope']
            test_dset = test_fid['relative-slope']

            npt.assert_almost_equal(test_dset, ref_dset,
                                    decimal=self.decimal_precision)


def _parser():
    """ Argument parser. """
    description = ("Unittests for `wagl.relative_azimuth_slope` function.\n"
                   "Comparisons tests will occur for the following "
                   "datasets: \n"
                   "\t* relative-slope\n")
                   
    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=RawTextHelpFormatter)
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

    return parser


def main():
    """ Main execution. """
    parser = _parser()
    args = parser.parse_args()

    reference_fname = args.reference_fname
    test_fname = args.test_fname
    decimal_precision = args.decimal_precision
    integer_precision = args.integer_precision

    suite = unittest.TestSuite()
    test_case = ParameterisedTestCase()
    suite.addTest(test_case.parameterise(TestRelativeSlope,
                                         reference_fname=reference_fname,
                                         test_fname=test_fname,
                                         decimal_precision=decimal_precision,
                                         integer_precision=integer_precision))
    unittest.TextTestRunner(verbosity=2).run(suite)
