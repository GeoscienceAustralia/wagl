#!/usr/bin/env python

import os
import sys
import unittests

import numpy
from osgeo import gdal

class ParameterisedTestCase(unittest.TestCase):
    """ TestCase classes that want to be parameterised should
        inherit from this class.

        Source for this code taken from:
        http://eli.thegreenplace.net/2011/08/02/python-unit-testing-parametrized-test-cases/

        Modified to suit our given parameters.
    """
    def __init__(self, methodName='runTest', reference_dir=None, test_dir=None,
                 decimal_precision=4, integer_precision=1):
        super(ParameterisedTestCase, self).__init__(methodName)

        # Reference and Testing directories
        self.reference_dir = reference_dir
        self.test_dir      = test_dir

        # Allowable numerical precision
        self.decimal_precision = decimal_precision
        self.integer_precision = integer_precision


    @staticmethod
    def parameterise(testcase_klass, reference_dir=None, test_dir=None,
                    decimal_precision=4, integer_precision=1):
        """ Create a suite containing all tests taken from the given
            subclass, passing them the parameters 'reference_dir,
            test_dir, decimal_precision, integer_precision'.
        """
        testloader = unittest.TestLoader()
        testnames = testloader.getTestCaseNames(testcase_klass)
        suite = unittest.TestSuite()
        for name in testnames:
            suite.addTest(testcase_klass(name, reference_dir=reference_dir,
                          test_dir=test_dir,
                          decimal_precision=decimal_precision,
                          integer_precision=integer_precision))
        return suite


def write_img(array, name='', format='ENVI', projection=None, geotransform=None):
    """
    A small and simple function to write single or multiband 2D arrays.
    Practical, lightweight and independent for the unittesting code.
    """

    dtypes = {
              'uint8'     : 1,
              'uint16'    : 2,
              'int16'     : 3,
              'uint32'    : 4,
              'int32'     : 5,
              'float32'   : 6,
              'float64'   : 7,
              'complex64' : 8,
              'complex64' : 9,
              'complex64' : 10,
              'complex128': 11,
              'bool'      : 1
             }

    dims = array.shape

    if (len(dims) == 2):
        samples = dims[1]
        lines   = dims[0]
        bands   = 1
    elif (len(dims) == 3):
        samples = dims[2]
        lines   = dims[1]
        bands   = dims[0]
    else:
        print 'Input array is not of 2 or 3 dimensions!!!'
        print 'Array dimensions: ', len(dims)
        return

    dtype = dtypes.get(array.dtype.name, 7)

    drv = gdal.GetDriverByName(format)
    outds = drv.Create(name, samples, lines, bands, dtype)

    if prj:
        outds.SetProjection(prj)

    if geot:
        outds.SetGeoTransform(geot)

    if bands > 1:
        for i in range(bands):
            band = outds.GetRasterBand(i+1)
            band.WriteArray(array)
            band.FlushCache()
    else:
        band = outds.GetRasterBand(1)
        band.WriteArray(array)
        band.FlushCache()

    outds = None


def read_img(fname):
    """
    A small and simple routine to read a GDAL compliant image.
    This is only intended for reading the raw file into a NumPy memory
    variable. It suits the purposes of unittesting.
    """

    ds = gdal.Open(fname)

    img = ds.ReadAsArray()

    ds = None

    return img


def find_file(dir, file):
    """
    A simple routine for checking existance of files on disk.
    No error catching, it'll bail out of the main level program
    as it is designed for the unittests.
    """
    fname = os.path.join(dir, file)
    if os.path.isfile(fname):
        return fname
    else:
        print "Error! file not found: %s" %fname
        sys.exit(2)

