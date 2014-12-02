#!/bin/env python

import unittest
from gaip import GriddedGeoBox
import affine
import rasterio as rio
from osgeo import osr
from osgeo import gdal

affine.EPSILON=1e-9
affine.EPSILON2=1e-18

def getFlindersIsletGGB():
    flindersOrigin = (150.927659, -34.453309)
    flindersCorner = (150.931697, -34.457915)
       
    return GriddedGeoBox.from_corners(flindersOrigin, flindersCorner)

class TestGriddedGeoBox(unittest.TestCase):

    def test_create(self):
        scale = 0.00025
        shape = (3,2)
        origin = (150.0, -34.0)
        corner = (shape[1]*scale+origin[0], origin[1]-shape[0]*scale)
        ggb = GriddedGeoBox(shape, origin)
        assert ggb is not None
        self.assertEqual(shape, ggb.shape)
        self.assertEqual(origin, ggb.origin)
        self.assertEqual(corner, ggb.corner)
 

    def test_create_small_offset_using_corners(self):
        # create small GGB centred on (150.00025,-34.00025)
        expectedShape = (2,2)
        scale = 0.00025
        delta = 0.00012
        origin = (150.0+delta, -34.0-delta)
        corner = (150.0+2*scale-delta, -34.0-2*scale+delta)
        # the grid should "pixel align"
        expectedOrigin = (150.0, -34.0)
        expectedCorner = (150.0+expectedShape[1]*scale, -34.0-expectedShape[0]*scale)
        ggb = GriddedGeoBox.from_corners(origin, corner)
        self.assertEqual(expectedShape, ggb.shape)
        self.assertEqual(expectedOrigin, ggb.origin)
        self.assertEqual(expectedCorner, ggb.corner)

    def test_create_unit_GGB_using_corners(self):
        # create small GGB centred on (150.00025,-34.00025)
        expectedShape = (1,1)
        scale = 0.00025
        origin = (150.0, -34.0)
        corner = (150.0+scale, -34.0-scale)
        ggb = GriddedGeoBox.from_corners(origin, corner)
        self.assertEqual(expectedShape, ggb.shape)
        self.assertEqual(corner, ggb.corner)

    def test_grid_alignment(self):
        scale = 0.00025
        shape = (20,30)
        origin = (150.00012, -34.00014)
        originShouldBe = (150.0, -34.0)
       
        corner = (shape[1]*scale+150, -34-shape[0]*scale)
        ggb = GriddedGeoBox(shape, origin)
        assert ggb is not None
        self.assertEqual(shape, ggb.shape)
        self.assertEqual(originShouldBe, ggb.origin)
        self.assertEqual(corner, ggb.corner)

    def test_real_world(self):
        # Flinders Islet, NSW
        flindersOrigin = (150.927659, -34.453309)
        flindersCorner = (150.931697, -34.457915)
        originShouldBe = (150.9275, -34.45325)
        cornerShouldBe = (150.93175, -34.45800)
        shapeShouldBe = (19, 17)
       
        ggb = GriddedGeoBox.from_corners(flindersOrigin, flindersCorner)
        assert ggb is not None
        self.assertEqual(shapeShouldBe, ggb.shape)
        self.assertAlmostEqual(originShouldBe[0], ggb.origin[0])
        self.assertAlmostEqual(originShouldBe[1], ggb.origin[1])
        self.assertAlmostEqual(cornerShouldBe[0], ggb.corner[0])
        self.assertAlmostEqual(cornerShouldBe[1], ggb.corner[1])

    def _land_sea_GGB_asserts(self, aGGB):

        assert aGGB is not None
        # print "aGGB.origin=",aGGB.origin
        # print "aGGB.corner=",aGGB.corner
        # print "aGGB.pixelsize=",aGGB.pixelsize
        # print "aGGB.shape=",aGGB.shape
        self.assertTrue(aGGB.origin == (-221573.33728165628, 9850031.808687024))
        self.assertTrue(aGGB.corner == (1088726.6627183438, 936806.8086870238))
        # self.assertTrue(aGGB.pixelsize == (25.0, 25.0))
        self.assertTrue(isinstance(aGGB.crs, osr.SpatialReference))

    def test_ggb_from_rio_dataset(self):
        # get Land/Sea data file for this bounding box
        utmZone = 56
        utmDataPath = '/g/data/v10/eoancillarydata/Land_Sea_Rasters/WORLDzone%d.tif' % (utmZone,)

        # read the data for the Flinders islet region
        with rio.open(utmDataPath) as ds:

            # get the gridded box for the full data extent
            datasetGGB = GriddedGeoBox.from_dataset(ds)
            # print
            # print "For Land/Sea read via rasterio:"
            self._land_sea_GGB_asserts(datasetGGB)

    def test_ggb_from_gdal_dataset(self):
        # get Land/Sea data file for this bounding box
        utmZone = 56
        utmDataPath = '/g/data/v10/eoancillarydata/Land_Sea_Rasters/WORLDzone%d.tif' % (utmZone,)

        # read the data for the Flinders islet region
        ds = gdal.Open(utmDataPath)

        # get the gridded box for the full data extent
        datasetGGB = GriddedGeoBox.from_dataset(ds)
        # print
        # print "For Land/Sea read via GDAL:"
        self._land_sea_GGB_asserts(datasetGGB)

    def no_test_ggb_subsets(self):
        # get Land/Sea data file for this bounding box
        utmZone = 56
        utmDataPath = '/g/data/v10/eoancillarydata/Land_Sea_Rasters/WORLDzone%d.tif' % (utmZone,)

        # read the data for the Flinders islet region
        with rio.open(utmDataPath) as ds:

            # get the gridded box for the full data extent
            datasetGGB = GriddedGeoBox.from_dataset(ds)

            flindersGGB = getFlindersIsletGGB()
            scale = 0.00025
            shapeYX = (3,3)
            origin = (150.0, -34.0)
            corner = (shapeYX[1]*scale+origin[0], origin[1]-shapeYX[0]*scale)
            ggb = GriddedGeoBox(shapeYX, origin)
            assert ggb is not None
            self.assertEqual(shapeYX, ggb.shape)
            self.assertEqual(origin, ggb.origin)
            self.assertEqual(corner, ggb.corner)


    def test_ggb_copy_and_rereference(self):
        # get Land/Sea data file for this bounding box
        utmZone = 56
        utmDataPath = '/g/data/v10/eoancillarydata/Land_Sea_Rasters/WORLDzone%d.tif' % (utmZone,)

        # read the data for the Flinders islet region
        with rio.open(utmDataPath) as ds:

            # get the gridded box for the full data extent
            datasetGGB = GriddedGeoBox.from_dataset(ds)

            # dataset will be UTM

            # print datasetGGB

            # copy to WSG84

            newGGB = datasetGGB.copy(crs='EPSG:4326')
            # print newGGB


    #TODO: define this test -- currently there are projection problems
    def no_test_ggb_subsets(self):
        # get Land/Sea data file for this bounding box
        utmZone = 56
        utmDataPath = '/g/data/v10/eoancillarydata/Land_Sea_Rasters/WORLDzone%d.tif' % (utmZone,)

        # read the data for the Flinders islet region
        with rio.open(utmDataPath) as ds:

            # get the gridded box for the full data extent
            datasetGGB = GriddedGeoBox.from_dataset(ds)

            flindersGGB = getFlindersIsletGGB()

            window = datasetGGB.window(flindersGGB)
            print 'window=', window
            
            windowShape=(window[1][1]-window[1][0], window[0][1]-window[0][0])
            print windowShape

            # transform the origin of the window

            windowOriginXY = (window[1][0], window[0][0])
            print 'windowOriginXY=', windowOriginXY

            windowOriginUTM = datasetGGB.affine * windowOriginXY
            print 'windowOriginUTM=', windowOriginUTM

            utm2wgs84 = osr.CoordinateTransformation( datasetGGB.crs, flindersGGB.crs)

            windowOrigin = utm2wgs84.TransformPoint(windowOriginUTM[0], windowOriginUTM[1])
            print 'windowOrigin=', windowOrigin

            flindersOrigin = (150.927659, -34.453309)
            print 'flindersOrigin=', flindersOrigin

            diff = (flindersOrigin[0]-windowOrigin[0], flindersOrigin[1]-windowOrigin[1])
            print 'diff=', diff

            wgs842utm = osr.CoordinateTransformation(flindersGGB.crs, datasetGGB.crs)

            (xUtm, yUtm, zzz) = wgs842utm.TransformPoint(flindersOrigin[0], flindersOrigin[1])

            (xWgs84, yWgs84, zzz) = utm2wgs84.TransformPoint(xUtm, yUtm)

            print '(xUtm, yUtm)=    ', (xUtm, yUtm)
            print '(xWgs84, yWgs84)=', (xWgs84, yWgs84)

  

            
            if not windowShape == flindersGGB.shape:
                self.fail("%s == %s" % (windowShape, flindersGGB.shape))


    def test_convert_coordinate_to_map(self):
        """
        Test that an input image/array co-ordinate is correctly
        converted to a map co-cordinate.
        Simple case: The first pixel.
        """
        # get Land/Sea data file for this bounding box
        utmZone = 56
        utmDataPath = '/g/data/v10/eoancillarydata/Land_Sea_Rasters/WORLDzone%d.tif' % (utmZone,)

        # read the data for the Flinders islet region
        with rio.open(utmDataPath) as ds:

            # get the gridded box for the full data extent
            datasetGGB = GriddedGeoBox.from_dataset(ds)

            xmap, ymap = datasetGGB.convert_coordinates((0,0))

            self.assertTrue(datasetGGB.origin == (xmap, ymap))


    def test_convert_coordinate_to_image(self):
        """
        Test that an input image/array co-ordinate is correctly
        converted to a map co-cordinate.
        Simple case: The first pixel.
        """
        # get Land/Sea data file for this bounding box
        utmZone = 56
        utmDataPath = '/g/data/v10/eoancillarydata/Land_Sea_Rasters/WORLDzone%d.tif' % (utmZone,)

        # read the data for the Flinders islet region
        with rio.open(utmDataPath) as ds:

            # get the gridded box for the full data extent
            datasetGGB = GriddedGeoBox.from_dataset(ds)

            ximg, yimg = datasetGGB.convert_coordinates(datasetGGB.origin, to_map=False)

            self.assertTrue((0, 0) == (ximg, yimg))


    def test_convert_coordinate_to_map_offset(self):
        """
        Test that an input image/array co-ordinate is correctly
        converted to a map co-cordinate using a pixel centre offset.
        Simple case: The first pixel.
        """
        # get Land/Sea data file for this bounding box
        utmZone = 56
        utmDataPath = '/g/data/v10/eoancillarydata/Land_Sea_Rasters/WORLDzone%d.tif' % (utmZone,)

        # read the data for the Flinders islet region
        with rio.open(utmDataPath) as ds:

            # get the gridded box for the full data extent
            datasetGGB = GriddedGeoBox.from_dataset(ds)

            xmap, ymap = datasetGGB.convert_coordinates((0,0), centre=True)

            # Get the actual centre co-ordinate of the first pixel
            xcentre, ycentre = datasetGGB.convert_coordinates((0.5, 0.5))

            self.assertTrue((xcentre, ycentre) == (xmap, ymap))


if __name__ == '__main__':
    unittest.main()
