#!/usr/bin/env python

from __future__ import absolute_import
import unittest

from posixpath import join as ppjoin
import numpy
import numpy.testing as npt
import osr

from gaip.constants import DatasetName
from gaip.longitude_latitude_arrays import create_lon_lat_grids
from gaip import unittesting_tools as ut

# WGS84
CRS = "EPSG:4326"

@unittest.skip('Requires refactoring')
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
        img, geobox = ut.create_test_image()
        fid = create_lon_lat_grids(geobox, depth=3)
        dataset_name = ppjoin(GroupName.lon_lat_group.value, DatasetName.lon)
        lon = fid[dataset_name][:]
        ids = ut.random_pixel_locations(img.shape)

        # We'll transform (reproject) to WGS84
        sr = osr.SpatialReference()
        sr.SetFromUserInput(CRS)

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
        img, geobox = ut.create_test_image()
        fid = create_lon_lat_grids(geobox, depth=3)
        dataset_name = ppjoin(GroupName.lon_lat_group.value, DatasetName.lat)
        lat = fid[dataset_name][:]
        ids = ut.random_pixel_locations(img.shape)

        # We'll transform (reproject) to WGS84
        sr = osr.SpatialReference()
        sr.SetFromUserInput(CRS)

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
    unittest.main()
