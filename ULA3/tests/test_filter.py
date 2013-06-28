import numpy as np
import unittest as ut
from pprint import pprint
from os.path import join, dirname
from osgeo import gdal, gdalconst
from ULA3.filter import filter_float as filter
from ULA3.filter import read_array

"""
Since we are lacking a portable test harness, these tests are configured to run on the NCI.

Hopefully, in the near future, we will have a test harness that can be more portable, and testing will be
possible elsewhere.
"""

class FilterTestCase(ut.TestCase):
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





def suite():
    return ut.TestSuite((FilterTestCase(),))
