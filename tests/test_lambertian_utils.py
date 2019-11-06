#!/usr/bin/env python

"""
Test the various utilites contained in the wagl.reflectance module.
"""

import unittest
import numpy

from wagl.convolution import _sequential_valid_rows, _fill_nulls


class SequentialRowsTest(unittest.TestCase):

    def setUp(self):
        # define a 'mask' where True are bad rows
        self.sequential = numpy.zeros((10, 10), dtype='bool')
        self.sequential[0] = True   # set first row
        self.sequential[-1] = True  # set last row

        self.non_sequential = numpy.zeros((10, 10), dtype='bool')
        self.non_sequential[0] = True   # set first row
        self.non_sequential[-1] = True  # set last row
        self.non_sequential[5] = True   # set middle row

    def test_sequential_start(self):
        start_row, end_row = _sequential_valid_rows(self.sequential)

        self.assertTrue(start_row, 1)

    def test_sequential_end(self):
        start_row, end_row = _sequential_valid_rows(self.sequential)

        self.assertTrue(end_row, 8)

    def test_non_sequential(self):
        with self.assertRaises(Exception):
            _sequential_valid_rows(self.non_sequential)


class FillNullsTest(unittest.TestCase):

    def setUp(self):
        # data array with some holes
        self.data = numpy.random.ranf((10, 10))
        self.data[0, 0] = numpy.nan
        self.data[3, 3] = numpy.nan
        self.data[5, 5] = numpy.nan
        self.data[7, 7] = numpy.nan

        # data mask
        self.data_mask = ~numpy.isfinite(self.data)  # noqa # pylint: disable

    def test_row_zero(self):
        # take a copy as the array is modified in-place
        copy = self.data.copy()
        xbar = numpy.nanmean(copy[0])
        _fill_nulls(copy, self.data_mask, None)
        self.assertTrue(xbar, copy[0, 0])

    def test_row_four(self):
        # take a copy as the array is modified in-place
        copy = self.data.copy()
        xbar = numpy.nanmean(copy[3])
        _fill_nulls(copy, self.data_mask, None)
        self.assertTrue(xbar, copy[3, 3])

    def test_row_six(self):
        # take a copy as the array is modified in-place
        copy = self.data.copy()
        xbar = numpy.nanmean(copy[5])
        _fill_nulls(copy, self.data_mask, None)
        self.assertTrue(xbar, copy[5, 5])

    def test_row_eight(self):
        # take a copy as the array is modified in-place
        copy = self.data.copy()
        xbar = numpy.nanmean(copy[7])
        _fill_nulls(copy, self.data_mask, None)
        self.assertTrue(xbar, copy[7, 7])


if __name__ == '__main__':
    unittest.main()
