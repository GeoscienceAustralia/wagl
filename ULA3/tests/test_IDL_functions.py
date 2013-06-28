#! /usr/bin/env python

import unittest
import numpy
import _test_IDL_Histogram

from IDL_functions import IDL_Histogram


class IDL_Hist_Tester(unittest.TestCase):

    def setUp(self):
        self.array1 = numpy.arange(10)
        self.control1 = self.array1 > 5
        self.array2 = numpy.arange(256)
        self.array3 = numpy.arange(10,20,0.5)
        self.array4 = numpy.random.ranf(1000)
        self.array5 = numpy.random.randint(0,11,(100,100))
        self.array6 = numpy.random.randint(0,256,(256))

    def test_true_false_a(self):
        '''
        Test that TRUE = 1 and FALSE = 0.
        '''
        # Using an array 0->9, test how many are > 5
        bool_ = numpy.zeros((10), dtype='int8')
        _test_IDL_Histogram.test_bool(self.array1, 10, bool_)
        self.assertEqual(self.control1.sum(), 4)

    def test_true_false_b(self):
        '''
        Test that the boolean array returned by Fortran gives the same
        result as that given by numpy.
        '''
        # Using an array 0->9, test how many are > 5
        bool_ = numpy.zeros((10), dtype='int8')
        _test_IDL_Histogram.test_bool(self.array1, 10, bool_)
        eq = self.control1 == bool_
        self.assertEqual(eq.sum(), 10)

    def test_hist(self):
        '''
        Test that the histogram works. Default binsize is 1, so there
        should be 256 bins.
        '''
        h = IDL_Histogram(self.array2)
        # Should be 256 elements, and the value 1 contained within each.
        self.assertEqual(h['histogram'].shape[0], 256)
        self.assertEqual((h['histogram'] == 1).sum(), 256)

    def test_hist_max(self):
        '''
        Test that the max keyword works.
        '''
        # Using an array 0->255, check that 255 gets omitted
        h = IDL_Histogram(self.array2, max=254)
        self.assertEqual(h['histogram'].shape[0], 255)
        self.assertEqual((h['histogram'] == 1).sum(), 255)

    def test_hist_min(self):
        '''
        Test that the min keyword works.
        '''
        # Using an array 0->255, check that 0 gets omitted
        h = IDL_Histogram(self.array2, min=1)
        self.assertEqual(h['histogram'].shape[0], 255)
        self.assertEqual((h['histogram'] == 1).sum(), 255)

    def test_omin(self):
        '''
        Test that the omin keyword works.
        '''
        # The output should be the same. Using an array 0->255
        h = IDL_Histogram(self.array2, omin='omin')
        self.assertEqual(h['omin'], 0)

    def test_omax(self):
        '''
        Test that the omin keyword works.
        '''
        # Using an array 0->255
        # The returned value should be the same as the derived max, unless
        # the nbins keyword is set, in which case the max gets rescaled by
        # nbins*binsize+min in order to maintain equal bin widths.
        h = IDL_Histogram(self.array2, omax='omax')
        self.assertEqual(h['omax'], 255)

    def test_nan(self):
        '''
        Test that the NaN keyword works.
        '''
        a = self.array2.astype('float64')
        a[0] = numpy.NaN
        h = IDL_Histogram(a, NaN=True)
        # The histogram will fail if array contains NaN's and NaN isn't set.
        # One element is excluded (the NaN), so test the length.
        self.assertEqual(h['histogram'].shape[0], 255)

    def test_binsize(self):
        '''
        Test that the binsize keyword works.
        '''
        h = IDL_Histogram(self.array3, binsize=0.5)
        # should be 20 bins to contain the values 10 -> 19.5
        self.assertEqual(h['histogram'].shape[0], 20)

    def test_default_binsize(self):
        '''
        Test that the default binsize is 1 and works accordingly.
        '''
        # Using an array of values in range 0->1
        h = IDL_Histogram(self.array4)
        self.assertEqual(h['histogram'].shape[0], 1)
        # All values should be in the first bin.
        self.assertEqual(self.array4.shape[0], h['histogram'][0])

    def test_nbins(self):
        '''
        Test that the nbins keyword works.
        '''
        h = IDL_Histogram(self.array4, nbins=256)
        # There should be 256 bins
        self.assertEqual(h['histogram'].shape[0], 256)

    def test_two_dimensional(self):
        '''
        Test that inputing a 2D array will raise an error.
        '''
        self.assertRaises(Exception, IDL_Histogram, self.array5)

    def test_reverse_indices1(self):
        '''
        Test that the reverse indices keyword works.
        '''
        # Make a copy then shuffle the array. Elements are in a random order.
        a = self.array2.copy()
        numpy.random.shuffle(a)
        h = IDL_Histogram(a, reverse_indices='ri')
        # Let's see if we can access the correct element. As we are dealing with
        # int's (and the binsize is one), pick a random element and the value
        # of the element represents the bin.
        # If reverse indices works, then the reeturned value should equal data.
        element = numpy.random.randint(0,256,(1))[0]
        bin = a[element]
        ri = h['ri']
        data = a[ri[ri[bin]:ri[bin+1]]]
        self.assertEqual(bin, data)

    def test_reverse_indices2(self):
        '''
        Test whether mulitple values in a single bin are correctly returned
        by the reverse indices.
        '''
        # Make a copy then shuffle the array. Elements are in a random order.
        a = self.array2.copy()
        numpy.random.shuffle(a)
        h = IDL_Histogram(a, reverse_indices='ri', binsize=5)
        # Using an array in the range 0->255, find data >=100<105
        # This should be bin 21 (20th if start from the 0th bin)
        ri = h['ri']
        # We know that each ri has adjacent groups (no empty bin), so no need
        # to check that ri[21] > ri[21]
        data = a[ri[ri[20]:ri[21]]]
        # The order should be the same as well. If not then numpy has changed.
        control = a[(a >= 100) & (a < 105)]
        self.assertEqual((control - data).sum(), 0)

    def test_reverse_indices3(self):
        '''
        Test that the reverse indices keyword works across multiple bins
        and values.
        '''
        # A random floating array in range 0-20
        a = (self.array4)*20
        # Specifying min=0 should give bin start points 0, 2.5, 5, 7.5 etc
        h = IDL_Histogram(a, reverse_indices='ri', min=0, binsize=2.5)
        # Find values >= 7.5 < 17.5
        control = numpy.sort(a[(a >= 7.5) & (a < 17.5)])
        ri = h['ri']
        # If the locations keyword was set then the starting locations of each
        # bin would be:
        # [  0. ,   2.5,   5. ,   7.5,  10. ,  12.5,  15. ,  17.5]
        # So we want bins 3, 4, 5, 6. (Bins start at 0)
        # Sort the arrays; so we can do an element by element difference
        data = numpy.sort(a[ri[ri[3]:ri[7]]])
        self.assertEqual((control - data).sum(), 0)

    def test_input1(self):
        '''
        Test that the input keyword works and has the same length as the
        number of expected bins.
        '''
        # Output should be of the same length as the number of bins.
        # We are using the default binsize, for values in range [0,255]
        a = self.array2
        b = self.array6
        h = IDL_Histogram(a, input=b)
        self.assertEqual(h['histogram'].shape[0], 256)

    def test_input2(self):
        '''
        Test that the input keyword works and that the input is correctly
        added to the histogram.
        The expected length of the histogram is the same size as the input
        array.
        '''
        # We are using the default binsize, for values in range [0,255]
        a = self.array2
        b = self.array6
        h = IDL_Histogram(a, input=b)['histogram']
        hcontrol = IDL_Histogram(a)['histogram'] + b
        self.assertEqual((h - hcontrol).sum(), 0)

    def test_input3(self):
        '''
        Test that the input keyword works and that the expected length of the
        histogram should take the length of the input array.
        '''
        # We are using the default binsize, for values in range [0,10)
        # Without using the input keyword, the histogram size should be 10.
        # However, the size of the input array is 256, so the output histogram
        # should be 256.
        a = self.array1
        b = self.array6
        h = IDL_Histogram(a, input=b)
        self.assertEqual(h['histogram'].shape[0], b.shape[0])

    def test_input4(self):
        '''
        Test that the input keyword works and that the input is correctly
        added to the histogram.
        The length of the histogram without setting the input keyword, is 10.
        However, an input array of length 256 will be used, thereby the length
        of the histogram will be 256.
        '''
        # We are using the default binsize, for values in range [0,10)
        a = self.array1
        b = self.array6
        h = IDL_Histogram(a, input=b)['histogram']
        # Elements h[0:9] should equal 1, and elements h[10:-1] should equal 0
        diff = h - b
        self.assertEqual(diff.sum(),10)

if __name__ == '__main__':
    unittest.main()
