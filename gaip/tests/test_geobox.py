#!/bin/env python

from __future__ import absolute_import
from __future__ import print_function
import unittest
import h5py
import gdal
import affine
import rasterio as rio
from osgeo import osr
from osgeo import gdal
from gaip.geobox import GriddedGeoBox
from gaip import unittesting_tools as ut

affine.EPSILON = 1e-9
affine.EPSILON2 = 1e-18

def getFlindersIsletGGB():
    flindersOrigin = (150.927659, -34.453309)
    flindersCorner = (150.931697, -34.457915)
       
    return GriddedGeoBox.from_corners(flindersOrigin, flindersCorner)

class TestGriddedGeoBox(unittest.TestCase):

    def test_create(self):
        scale = 0.00025
        shape = (3, 2)
        origin = (150.0, -34.0)
        corner = (shape[1] * scale + origin[0], origin[1] - shape[0] * scale)
        ggb = GriddedGeoBox(shape, origin)
        self.assertEqual(shape, ggb.shape)
        self.assertEqual(origin, ggb.origin)
        self.assertEqual(corner, ggb.corner)
 

    def test_create_unit_GGB_using_corners(self):
        # create small GGB centred on (150.00025,-34.00025)
        expectedShape = (1, 1)
        scale = 0.00025
        origin = (150.0, -34.0)
        corner = (150.0 + scale, -34.0 - scale)
        ggb = GriddedGeoBox.from_corners(origin, corner)
        self.assertEqual(expectedShape, ggb.shape)
        self.assertEqual(corner, ggb.corner)


    def test_agdc_copy(self):
        scale = 0.00025
        shape = (4000, 4000)
        origin = (150.0, -34.0)
       
        corner = (shape[1] * scale + 150, -34 - shape[0] * scale)
        ggb = GriddedGeoBox(shape, origin, pixelsize=(scale, scale))
        self.assertEqual(shape, ggb.shape)
        self.assertEqual(corner, ggb.corner)

        # now get UTM equilavent

        utm_ggb = ggb.copy(crs="EPSG:32756")

        self.assertAlmostEqual(utm_ggb.origin[0], 222908.70452156663)
        self.assertAlmostEqual(utm_ggb.origin[1], 6233785.283900621)
        self.assertAlmostEqual(utm_ggb.corner[0], 317483.90638409054)
        self.assertAlmostEqual(utm_ggb.corner[1], 6125129.365269075)
        self.assertAlmostEqual(utm_ggb.pixelsize[0], 23.643800465630978)
        self.assertAlmostEqual(utm_ggb.pixelsize[1], 27.16397965788655)


    def test_real_world(self):
        # Flinders Islet, NSW
        flindersOrigin = (150.927659, -34.453309)
        flindersCorner = (150.931697, -34.457915)
        originShouldBe = flindersOrigin
        shapeShouldBe = (19, 17)
        cornerShouldBe = (flindersOrigin[0] + shapeShouldBe[1] * 0.00025, \
            flindersOrigin[1] - shapeShouldBe[0] * 0.00025)
       
        ggb = GriddedGeoBox.from_corners(flindersOrigin, flindersCorner)
        self.assertEqual(shapeShouldBe, ggb.shape)
        self.assertAlmostEqual(originShouldBe[0], ggb.origin[0])
        self.assertAlmostEqual(originShouldBe[1], ggb.origin[1])
        self.assertAlmostEqual(cornerShouldBe[0], ggb.corner[0])
        self.assertAlmostEqual(cornerShouldBe[1], ggb.corner[1])


    def test_ggb_from_rio_dataset(self):
        img, geobox = ut.create_test_image()
        kwargs = {'driver': 'Memory',
                  'width': img.shape[1],
                  'height': img.shape[0],
                  'count': 1,
                  'transform': geobox.affine,
                  'crs': geobox.crs.ExportToWkt(),
                  'dtype': img.dtype.name}

        with rio.open('tmp.tif', 'w', **kwargs) as ds:
            new_geobox = GriddedGeoBox.from_rio_dataset(ds)

            self.assrtTrue(new_geobox.affine == geobox.affine)
            self.assrtTrue(new_geobox.crs.ExportToWkt() ==
                           geobox.crs.ExportToWkt())
            self.assrtTrue(new_geobox.shape == img.shape)


    def test_ggb_from_gdal_dataset(self):
        img, geobox = ut.create_test_image()
        drv = gdal.GetDriverByName('MEM')
        ds = drv.Create('tmp.tif', img.shape[1], img.shape[0], 1, 1)
        ds.SetGeoTransform(geobox.affine.to_gdal())
        ds.SetProjection(geobox.crs.ExportToWkt())

        new_geobox = GriddedGeoBox.from_gdal_dataset(ds)
        self.assrtTrue(new_geobox.affine == geobox.affine)
        self.assrtTrue(new_geobox.crs.ExportToWkt() ==
                       geobox.crs.ExportToWkt())
        self.assrtTrue(new_geobox.shape == img.shape)
        drv = None
        ds = None


    def test_ggb_from_h5_dataset(self):
        img, geobox = ut.create_test_image()
        with h5py.File('tmp.h5', driver='core', backing_store=False) as fid:
            ds = fid.create_dataset('test', data=img)
            ds.attrs['geotransform'] = geobox.affine.to_gdal()
            ds.attrs['crs_wkt'] = geobox.crs.ExportToWkt()

            new_geobox = GriddedGeoBox.from_h5_dataset(ds)
            self.assrtTrue(new_geobox.affine == geobox.affine)
            self.assrtTrue(new_geobox.crs.ExportToWkt() ==
                           geobox.crs.ExportToWkt())
            self.assrtTrue(new_geobox.shape == img.shape)


    def test_convert_coordinate_to_map(self):
        """
        Test that an input image/array co-ordinate is correctly
        converted to a map co-cordinate.
        Simple case: The first pixel.
        """
        _, geobox = ut.create_test_image()
        xmap, ymap = geobox.convert_coordinates((0, 0))
        self.assertTrue(geobox.origin == (xmap, ymap))


    def test_convert_coordinate_to_image(self):
        """
        Test that an input image/array co-ordinate is correctly
        converted to a map co-cordinate.
        Simple case: The first pixel.
        """
        _, geobox = ut.create_test_image()
        ximg, yimg = geobox.convert_coordinates(geobox.origin, to_map=False)
        self.assertTrue((0, 0) == (ximg, yimg))


    def test_convert_coordinate_to_map_offset(self):
        """
        Test that an input image/array co-ordinate is correctly
        converted to a map co-cordinate using a pixel centre offset.
        Simple case: The first pixel.
        """
        _, geobox = ut.create_test_image()
        xmap, ymap = geobox.convert_coordinates((0, 0), centre=True)

        # Get the actual centre co-ordinate of the first pixel
        xcentre, ycentre = geobox.convert_coordinates((0.5, 0.5))
        self.assertTrue((xcentre, ycentre) == (xmap, ymap))


    def test_pixelscale_metres(self):
        scale = 0.00025
        shape = (4000, 4000)
        origin = (150.0, -34.0)
        ggb = GriddedGeoBox(shape, origin, pixelsize=(scale, scale))
        (size_x, size_y) = ggb.get_pixelsize_metres(xy=(0, 0))
        self.assertAlmostEqual(size_x, 23.0962, places=4)
        self.assertAlmostEqual(size_y, 27.7306, places=4)


    def test_all_pixelscale_metres(self):
        scale = 0.00025
        shape = (4000, 4000)
        origin = (150.0, -34.0)
        ggb = GriddedGeoBox(shape, origin, pixelsize=(scale, scale))
        size_array = ggb.get_all_pixelsize_metres()
        
        self.assertEqual(len(size_array), 4000)
        (size_x, size_y) = size_array[0]
        self.assertAlmostEqual(size_x, 23.0962, places=4)
        self.assertAlmostEqual(size_y, 27.7306, places=4)
        (size_x, size_y) = size_array[3999]
        self.assertAlmostEqual(size_x, 22.8221, places=4)
        self.assertAlmostEqual(size_y, 27.7351, places=4)


if __name__ == '__main__':
    unittest.main()
