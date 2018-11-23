#!/usr/bin/env python
"""
Unittesting framework for the `wagl.terrain_shadow_masks` module.
"""

from __future__ import absolute_import, print_function, unicode_literals
import unittest
import argparse
from argparse import RawTextHelpFormatter

import numpy
import h5py

from wagl.unittesting_tools import ParameterisedTestCase


class TestShadowMasks(ParameterisedTestCase):
    """
    Unittesting for the self shadow, cast shadow (from both
    the sun and satellite directions), and the combined shadow
    masks, found in `wagl.terrain_shadow_masks`.

    Unittests will occur for the following datasets:

        * self-shadow
        * cast-shadow-sun
        * cast-shadow-satellite
        * combined-shadow

    We're not explicitly testing the function, but implicitly
    testing the function by comparing against an existing
    dataset.
    """

    def test_self_shadow(self):
        """
        Test the self shadow mask.
        """
        with h5py.File(self.reference_fname, 'r') as reference_fid,\
            h5py.File(self.test_fname, 'r') as test_fid:

            ref_dset = reference_fid['self-shadow']
            test_dset = test_fid['self-shadow']

            self.assertTrue(numpy.array_equal(ref_dset, test_dset))

    def test_cast_shadow_sun(self):
        """
        Test the cast shadow mask (sun direction).
        """
        with h5py.File(self.reference_fname, 'r') as reference_fid,\
            h5py.File(self.test_fname, 'r') as test_fid:

            ref_dset = reference_fid['cast-shadow-sun']
            test_dset = test_fid['cast-shadow-sun']

            self.assertTrue(numpy.array_equal(ref_dset, test_dset))

    def test_cast_shadow_satellite(self):
        """
        Test the cast shadow mask (satellite direction).
        """
        with h5py.File(self.reference_fname, 'r') as reference_fid,\
            h5py.File(self.test_fname, 'r') as test_fid:

            ref_dset = reference_fid['cast-shadow-satellite']
            test_dset = test_fid['cast-shadow-satellite']

            self.assertTrue(numpy.array_equal(ref_dset, test_dset))

    def test_combined_shadow(self):
        """
        Test the combined shadow mask.
        """
        with h5py.File(self.reference_fname, 'r') as reference_fid,\
            h5py.File(self.test_fname, 'r') as test_fid:

            ref_dset = reference_fid['combined-shadow']
            test_dset = test_fid['combined-shadow']

            self.assertTrue(numpy.array_equal(ref_dset, test_dset))


def _parser():
    """ Argument parser. """
    description = ("Unittests for `wagl.terrain_shadow_masks` module.\n"
                   "Comparisons tests will occur for the following "
                   "datasets: \n"
                   "\t* self-shadow\n"
                   "\t* cast-shadow-sun\n"
                   "\t* cast-shadow-satellite\n"
                   "\t* combined-shadow\n")
                   
    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=RawTextHelpFormatter)
    parser.add_argument('--reference_fname', required=True,
                        help=('The filename containing the reference datasets '
                              'to be used as a baseline.'))
    parser.add_argument('--test_fname', required=True,
                        help=('The filename containing the test datasets '
                              'to be used in comparing against the '
                              'base/reference datasets.'))

    return parser


def main():
    """ Main execution. """
    parser = _parser()
    args = parser.parse_args()

    reference_fname = args.reference_fname
    test_fname = args.test_fname

    suite = unittest.TestSuite()
    test_case = ParameterisedTestCase()
    suite.addTest(test_case.parameterise(TestShadowMasks,
                                         reference_fname=reference_fname,
                                         test_fname=test_fname))
    unittest.TextTestRunner(verbosity=2).run(suite)
