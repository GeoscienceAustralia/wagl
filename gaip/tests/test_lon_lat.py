#!/usr/bin/env python

import argparse
import os
import unittest

import numpy
import numpy.testing as npt
import osr
import rasterio

from EOtools import blrb

from gaip import GriddedGeoBox
from gaip import create_lon_lat_grids

# GDA94/ MGA Zone 55
CRS1 = "EPSG:28355"

# WGS84
CRS2 = "EPSG:4326"

def createTestImage(dimensions=(1000,1000), geotransform=None,
        projection=None, resolution=(25.0,25.0)):
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

    img = numpy.random.randint(0, 256, dimensions).astype('uint8')

    if (geotransform is None) or (projection is None):
        geotransform = (635000.0, 25.0, 0.0, 6277000.0, 0.0, 25.0)

        sr = osr.SpatialReference()
        # GDA94/ MGA Zone 55
        sr.SetFromUserInput(CRS1)
        projection = sr.ExportToWkt()
        resolution = (geotransform[1], geotransform[5])

    UL = (geotransform[0], geotransform[3])

    geobox = GriddedGeoBox(shape=dimensions, origin=UL, pixelsize=resolution,
        crs=projection)

    return (img, geobox)

def randomPixelLocations(dimensions, nPixels=10):
    """
    Given a tuple of (y, x) dimensions, generate a random index of
    pixel locations.
    Returns a standard NumPy tuple index i.e. (y, x) where y & x are
    NumPy arrays of length determined by the nPixels keyword.

    :param dimensions:
        A tuple containing the (y, x) dimensions of the 2D image of
        interest. The dimensions are used to confine the pixel
        indices to the 2D space.

    :param nPixels:
        An integer representing the desired number of random pixels.

    :return:
        A tuple of (y, x) where both y & x are 1D NumPy arrays.
    """
    cols = dimensions[1]
    rows = dimensions[0]
    x = numpy.random.randint(0, cols, nPixels)
    y = numpy.random.randint(0, rows, nPixels)

    index = (y,x)
    return index

class TestLonLatArrays(unittest.TestCase):

    def test_lon_array(self):
        """
        Test that the interpolated longitude array has sensible
        values.

        :notes:
        The default image that is used is a 1000 by 1000. This might
        be too small for a recursive depth of 7 (default within the
        NBAR framework). As such a depth of 3 is used and produced
        correct results. If a larger image is used, then the depth
        level should probably be increased.
        """
        # Initialise the test data
        img, geobox = createTestImage()
        lon, _ = create_lon_lat_grids(geobox, depth=3, to_disk=False)
        ids = randomPixelLocations(img.shape)

        # We'll transform (reproject) to WGS84
        sr = osr.SpatialReference()
        sr.SetFromUserInput(CRS2)

        # Get a list of x reprojected co-ordinates
        reprj = []
        for i in range(ids[0].shape[0]):
            # Get pixel (x,y)
            xy = (ids[1][i], ids[0][i])
            # Convert pixel to map
            mapXY = geobox.convert_coordinates(xy, centre=True)
            # Transform map to another crs
            x, _ = geobox.transform_coordinates(mapXY, to_crs=sr)
            reprj.append(x)

        reprj = numpy.array(reprj)

        lon_values = lon[ids]

        # A decimal degree co-ordinate should be correct up to the 6th dp
        self.assertIsNone(npt.assert_almost_equal(lon_values, reprj,
            decimal=6))

    def test_lat_array(self):
        """
        Test that the interpolated latitude array has sensible
        values.

        :notes:
        The default image that is used is a 1000 by 1000. This might
        be too small for a recursive depth of 7 (default within the
        NBAR framework). As such a depth of 3 is used and produced
        correct results. If a larger image is used, then the depth
        level should probably be increased.
        """
        # Initialise the test data
        img, geobox = createTestImage()
        _, lat = create_lon_lat_grids(geobox, depth=3, to_disk=False)
        ids = randomPixelLocations(img.shape)

        # We'll transform (reproject) to WGS84
        sr = osr.SpatialReference()
        sr.SetFromUserInput(CRS2)

        # Get a list of the y reprojected co-ordinates
        reprj = []
        for i in range(ids[0].shape[0]):
            # Get pixel (x,y)
            xy = (ids[1][i], ids[0][i])
            # Convert pixel to map
            mapXY = geobox.convert_coordinates(xy, centre=True)
            # Transform map to another crs
            _, y = geobox.transform_coordinates(mapXY, to_crs=sr)
            reprj.append(y)

        reprj = numpy.array(reprj)

        lat_values = lat[ids]

        # A decimal degree co-ordinate should be correct up to the 6th dp
        self.assertIsNone(npt.assert_almost_equal(lat_values, reprj,
            decimal=6))

if __name__ == '__main__':
    fname = '/g/data1/v10/NBAR_validation_reference/Nov2013/NBAR/LS5_90-84_1996-08-25/UTM/output/LS5_TM_NBAR_P54_GANBAR01-002_090_084_19960825/scene01/LS5_TM_NBAR_P54_GANBAR01-002_090_084_19960825_B40.tif'
    ds = rasterio.open(fname)
    geobox = GriddedGeoBox.from_rio_dataset(ds)
    outdir = '/g/data1/v10/testing_ground/jps547/test_lon_lat'
    create_lon_lat_grids(geobox, work_dir=outdir, lon_fname='LON_new2.tif', lat_fname='LAT_new2.tif')
    unittest.main()

