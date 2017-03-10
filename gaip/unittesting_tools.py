#!/usr/bin/env python

from __future__ import absolute_import
import unittest

import numpy
import osr

from gaip.geobox import GriddedGeoBox

# GDA94/ MGA Zone 55
CRS = "EPSG:28355"

class ParameterisedTestCase(unittest.TestCase):

    """ 
    TestCase classes that want to be parameterised should
    inherit from this class.

    Source for this code taken from:
    http://eli.thegreenplace.net/2011/08/02/python-unit-testing-parametrized-test-cases/

    Modified to suit our given parameters.
    """

    def __init__(self, method_name='runTest', reference_fname=None,
                 test_fname=None, tolerance=3, decimal_precision=4,
                 integer_precision=1):
        super(ParameterisedTestCase, self).__init__(method_name)

        # Reference and Testing directories
        self.reference_fname = reference_fname
        self.test_fname = test_fname

        # Allowable tollerance
        self.tolerance = tolerance

        # Allowable numerical precision
        self.decimal_precision = decimal_precision
        self.integer_precision = integer_precision

    @staticmethod
    def parameterise(testcase_klass, reference_fname=None, test_fname=None,
                     tolerance=3, decimal_precision=4, integer_precision=1):
        """
        Create a suite containing all tests taken from the given
        subclass, passing them the parameters 'reference_fname,
        test_fname, decimal_precision, integer_precision'.

        :param testcase_klass:
            A unittest.TestCase Class

        :param reference_fname:
            A full file pathname oo the file containing the
            reference data.

        :param test_fname:
            A full file pathname to the file containing the
            test data.

        :param tolerance:
            Default is 3, i.e. tolerance is 97% must be the same in
            order to Pass.

        :param decimal_precision:
            The decimal precision to be used during array comparison.
            Default is 4, i.e. values must be correct up to 4 d.p. in
            order to Pass.

        :param integer_precision:
            The intger precision to be used during array comparison.
            Default is 1, i.e. values must be correct within 1 integer
            in order to Pass.
        """
        testloader = unittest.TestLoader()
        testnames = testloader.getTestCaseNames(testcase_klass)
        suite = unittest.TestSuite()
        for name in testnames:
            suite.addTest(testcase_klass(name, reference_fname=reference_fname,
                          test_fname=test_fname,
                          tolerance=tolerance,
                          decimal_precision=decimal_precision,
                          integer_precision=integer_precision))
        return suite


def create_test_image(dimensions=(1000,1000), geotransform=None,
                      projection=None, resolution=(25.0,25.0), dtype='uint8'):
    """
    Creates an image with geo-location information.

    :param dimensions:
        A tuple containing the (y, x) dimensions of the 2D image to
        be generated.

    :param geotransform:
        A list or tuple containing the upper left co-ordinate of the image.
        This info can be retrieved from gdal. Otherwise create your own using
        the following as a guide. Must have 6 elements.
        geoT = (635000.0, 25.0, 0.0, 6277000.0, 0.0, 25.0)
        geoT[0] is top left x co-ordinate.
        geoT[1] is west to east pixel size.
        geoT[2] is image rotation (0 if image is north up).
        geoT[3] is to left y co-ordinate.
        geoT[4] is image rotation (0 if image is north up).
        geoT[5] is north to south pixel size.

        If either the geotransform or projection keywords are None,
        then geotransform will be set to:
        (635000.0, 25.0, 0.0, 6277000.0, 0.0, 25.0)
        and the projection will be set to EPSG:28355.

    :param projection:
        An osr compliant projection input such as WKT.
        If either the projection or geotransform keywords are None,
        then the projection will be set to EPSG:28355 and the
        geotransform will be set to:
        (635000.0, 25.0, 0.0, 6277000.0, 0.0, 25.0)

    :param resolution:
        A tuple containing the (x, y) pixel resolution/size. Default
        is (25.0, 25.0).

    :return:
        A tuple of two elements.
        The 1st element contains a random 8bit Unsigned Integer of
        dimensions (y, x), containing values in the range [0,256).
        The 2nd element contains an instance of a GriddedGeoBox.
    """
    img = numpy.random.randint(0, 256, dimensions).astype(dtype)

    if (geotransform is None) or (projection is None):
        geotransform = (635000.0, 25.0, 0.0, 6277000.0, 0.0, 25.0)

        sr = osr.SpatialReference()
        # GDA94/ MGA Zone 55
        sr.SetFromUserInput(CRS)
        projection = sr.ExportToWkt()
        resolution = (geotransform[1], geotransform[5])

    UL = (geotransform[0], geotransform[3])

    geobox = GriddedGeoBox(shape=dimensions, origin=UL, pixelsize=resolution,
        crs=projection)

    return (img, geobox)

def random_pixel_locations(dimensions, npixels=100):
    """
    Given a tuple of (y, x) dimensions, generate a random index of
    pixel locations.
    Returns a standard NumPy tuple index i.e. (y, x) where y & x are
    NumPy arrays of length determined by the nPixels keyword.

    :param dimensions:
        A tuple containing the (y, x) dimensions of the 2D image of
        interest. The dimensions are used to confine the pixel
        indices to the 2D space.

    :param npixels:
        An integer representing the desired number of random pixels.
        Default is 100 pixels.

    :return:
        A tuple of (y, x) where both y & x are 1D NumPy arrays.
    """
    if len(dimensions) != 2:
        raise TypeError('2 Dimensions must be specified.')

    cols = dimensions[1]
    rows = dimensions[0]
    x = numpy.random.randint(0, cols, npixels)
    y = numpy.random.randint(0, rows, npixels)

    index = (y,x)
    return index

