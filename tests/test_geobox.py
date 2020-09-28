#!/bin/env python

from __future__ import absolute_import
from __future__ import print_function
import unittest
import numpy
import numpy.testing as npt
import h5py
import gdal
import affine
import rasterio as rio
from osgeo import osr
from osgeo import gdal

from wagl.acquisition import acquisitions
from wagl.geobox import GriddedGeoBox
from wagl import unittesting_tools as ut

from .data import LS8_SCENE1

affine.EPSILON = 1e-9
affine.EPSILON2 = 1e-18


def getFlindersIsletGGB():
    flindersOrigin = (150.927659, -34.453309)
    flindersCorner = (150.931697, -34.457915)

    return GriddedGeoBox.from_corners(flindersOrigin, flindersCorner)


class TestGriddedGeoBox(unittest.TestCase):
    def test_create_shape(self):
        scale = 0.00025
        shape = (3, 2)
        origin = (150.0, -34.0)
        corner = (shape[1] * scale + origin[0], origin[1] - shape[0] * scale)
        ggb = GriddedGeoBox(shape, origin)
        self.assertEqual(shape, ggb.shape)

    def test_get_shape_xy(self):
        scale = 0.00025
        shape = (3, 2)
        shape_xy = (2, 3)
        origin = (150.0, -34.0)
        corner = (shape[1] * scale + origin[0], origin[1] - shape[0] * scale)
        ggb = GriddedGeoBox(shape, origin)
        self.assertEqual(shape_xy, ggb.get_shape_xy())

    def test_get_shape_yx(self):
        scale = 0.00025
        shape = (3, 2)
        origin = (150.0, -34.0)
        corner = (shape[1] * scale + origin[0], origin[1] - shape[0] * scale)
        ggb = GriddedGeoBox(shape, origin)
        self.assertEqual(shape, ggb.get_shape_yx())

    def test_x_size(self):
        scale = 0.00025
        shape = (3, 2)
        origin = (150.0, -34.0)
        corner = (shape[1] * scale + origin[0], origin[1] - shape[0] * scale)
        ggb = GriddedGeoBox(shape, origin)
        self.assertEqual(shape[1], ggb.x_size())

    def test_y_size(self):
        scale = 0.00025
        shape = (3, 2)
        origin = (150.0, -34.0)
        corner = (shape[1] * scale + origin[0], origin[1] - shape[0] * scale)
        ggb = GriddedGeoBox(shape, origin)
        self.assertEqual(shape[0], ggb.y_size())

    def test_create_origin(self):
        scale = 0.00025
        shape = (3, 2)
        origin = (150.0, -34.0)
        corner = (shape[1] * scale + origin[0], origin[1] - shape[0] * scale)
        ggb = GriddedGeoBox(shape, origin)
        self.assertEqual(origin, ggb.origin)

    def test_create_corner(self):
        scale = 0.00025
        shape = (3, 2)
        origin = (150.0, -34.0)
        corner = (shape[1] * scale + origin[0], origin[1] - shape[0] * scale)
        ggb = GriddedGeoBox(shape, origin)
        self.assertEqual(corner, ggb.corner)

    def test_shape_create_unit_GGB_using_corners(self):
        # create small GGB centred on (150.00025,-34.00025)
        expectedShape = (1, 1)
        scale = 0.00025
        origin = (150.0, -34.0)
        corner = (150.0 + scale, -34.0 - scale)
        ggb = GriddedGeoBox.from_corners(origin, corner)
        self.assertEqual(expectedShape, ggb.shape)

    def test_corner_create_unit_GGB_using_corners(self):
        # create small GGB centred on (150.00025,-34.00025)
        expectedShape = (1, 1)
        scale = 0.00025
        origin = (150.0, -34.0)
        corner = (150.0 + scale, -34.0 - scale)
        ggb = GriddedGeoBox.from_corners(origin, corner)
        self.assertEqual(corner, ggb.corner)

    def test_real_world_shape(self):
        # Flinders Islet, NSW
        flindersOrigin = (150.927659, -34.453309)
        flindersCorner = (150.931697, -34.457915)
        originShouldBe = flindersOrigin
        shapeShouldBe = (19, 17)
        cornerShouldBe = (
            flindersOrigin[0] + shapeShouldBe[1] * 0.00025,
            flindersOrigin[1] - shapeShouldBe[0] * 0.00025,
        )

        ggb = GriddedGeoBox.from_corners(flindersOrigin, flindersCorner)
        self.assertEqual(shapeShouldBe, ggb.shape)

    def test_real_world_origin_lon(self):
        # Flinders Islet, NSW
        flindersOrigin = (150.927659, -34.453309)
        flindersCorner = (150.931697, -34.457915)
        originShouldBe = flindersOrigin
        shapeShouldBe = (19, 17)
        cornerShouldBe = (
            flindersOrigin[0] + shapeShouldBe[1] * 0.00025,
            flindersOrigin[1] - shapeShouldBe[0] * 0.00025,
        )

        ggb = GriddedGeoBox.from_corners(flindersOrigin, flindersCorner)
        self.assertEqual(shapeShouldBe, ggb.shape)
        self.assertAlmostEqual(originShouldBe[0], ggb.origin[0])

    def test_real_world_origin_lat(self):
        # Flinders Islet, NSW
        flindersOrigin = (150.927659, -34.453309)
        flindersCorner = (150.931697, -34.457915)
        originShouldBe = flindersOrigin
        shapeShouldBe = (19, 17)
        cornerShouldBe = (
            flindersOrigin[0] + shapeShouldBe[1] * 0.00025,
            flindersOrigin[1] - shapeShouldBe[0] * 0.00025,
        )

        ggb = GriddedGeoBox.from_corners(flindersOrigin, flindersCorner)
        self.assertAlmostEqual(originShouldBe[1], ggb.origin[1])

    def test_real_world_corner_lon(self):
        # Flinders Islet, NSW
        flindersOrigin = (150.927659, -34.453309)
        flindersCorner = (150.931697, -34.457915)
        originShouldBe = flindersOrigin
        shapeShouldBe = (19, 17)
        cornerShouldBe = (
            flindersOrigin[0] + shapeShouldBe[1] * 0.00025,
            flindersOrigin[1] - shapeShouldBe[0] * 0.00025,
        )

        ggb = GriddedGeoBox.from_corners(flindersOrigin, flindersCorner)
        self.assertAlmostEqual(cornerShouldBe[0], ggb.corner[0])

    def test_real_world_corner_lat(self):
        # Flinders Islet, NSW
        flindersOrigin = (150.927659, -34.453309)
        flindersCorner = (150.931697, -34.457915)
        originShouldBe = flindersOrigin
        shapeShouldBe = (19, 17)
        cornerShouldBe = (
            flindersOrigin[0] + shapeShouldBe[1] * 0.00025,
            flindersOrigin[1] - shapeShouldBe[0] * 0.00025,
        )

        ggb = GriddedGeoBox.from_corners(flindersOrigin, flindersCorner)
        self.assertAlmostEqual(cornerShouldBe[1], ggb.corner[1])

    def test_ggb_transform_from_rio_dataset(self):
        img, geobox = ut.create_test_image()
        kwargs = {
            "driver": "MEM",
            "width": img.shape[1],
            "height": img.shape[0],
            "count": 1,
            "transform": geobox.transform,
            "crs": geobox.crs.ExportToWkt(),
            "dtype": img.dtype.name,
        }

        with rio.open("tmp.tif", "w", **kwargs) as ds:
            new_geobox = GriddedGeoBox.from_rio_dataset(ds)

            self.assertTrue(new_geobox.transform == geobox.transform)

    def test_ggb_crs_from_rio_dataset(self):
        img, geobox = ut.create_test_image()
        kwargs = {
            "driver": "MEM",
            "width": img.shape[1],
            "height": img.shape[0],
            "count": 1,
            "transform": geobox.transform,
            "crs": geobox.crs.ExportToWkt(),
            "dtype": img.dtype.name,
        }

        with rio.open("tmp.tif", "w", **kwargs) as ds:
            new_geobox = GriddedGeoBox.from_rio_dataset(ds)

            self.assertTrue(new_geobox.crs.ExportToWkt() == geobox.crs.ExportToWkt())

    def test_ggb_shape_from_rio_dataset(self):
        img, geobox = ut.create_test_image()
        kwargs = {
            "driver": "MEM",
            "width": img.shape[1],
            "height": img.shape[0],
            "count": 1,
            "transform": geobox.transform,
            "crs": geobox.crs.ExportToWkt(),
            "dtype": img.dtype.name,
        }

        with rio.open("tmp.tif", "w", **kwargs) as ds:
            new_geobox = GriddedGeoBox.from_rio_dataset(ds)

            self.assertTrue(new_geobox.shape == img.shape)

    def test_ggb_transform_from_gdal_dataset(self):
        img, geobox = ut.create_test_image()
        drv = gdal.GetDriverByName("MEM")
        ds = drv.Create("tmp.tif", img.shape[1], img.shape[0], 1, 1)
        ds.SetGeoTransform(geobox.transform.to_gdal())
        ds.SetProjection(geobox.crs.ExportToWkt())

        new_geobox = GriddedGeoBox.from_gdal_dataset(ds)
        self.assertTrue(new_geobox.transform == geobox.transform)
        drv = None
        ds = None

    def test_ggb_crs_from_gdal_dataset(self):
        img, geobox = ut.create_test_image()
        drv = gdal.GetDriverByName("MEM")
        ds = drv.Create("tmp.tif", img.shape[1], img.shape[0], 1, 1)
        ds.SetGeoTransform(geobox.transform.to_gdal())
        ds.SetProjection(geobox.crs.ExportToWkt())

        new_geobox = GriddedGeoBox.from_gdal_dataset(ds)
        self.assertTrue(new_geobox.crs.ExportToWkt() == geobox.crs.ExportToWkt())
        drv = None
        ds = None

    def test_ggb_shape_from_gdal_dataset(self):
        img, geobox = ut.create_test_image()
        drv = gdal.GetDriverByName("MEM")
        ds = drv.Create("tmp.tif", img.shape[1], img.shape[0], 1, 1)
        ds.SetGeoTransform(geobox.transform.to_gdal())
        ds.SetProjection(geobox.crs.ExportToWkt())

        new_geobox = GriddedGeoBox.from_gdal_dataset(ds)
        self.assertTrue(new_geobox.shape == img.shape)
        drv = None
        ds = None

    def test_ggb_transform_from_h5_dataset(self):
        img, geobox = ut.create_test_image()
        with h5py.File("tmp.h5", "w", driver="core", backing_store=False) as fid:
            ds = fid.create_dataset("test", data=img)
            ds.attrs["geotransform"] = geobox.transform.to_gdal()
            ds.attrs["crs_wkt"] = geobox.crs.ExportToWkt()

            new_geobox = GriddedGeoBox.from_h5_dataset(ds)
            self.assertTrue(new_geobox.transform == geobox.transform)

    def test_ggb_crs_from_h5_dataset(self):
        img, geobox = ut.create_test_image()
        with h5py.File("tmp.h5", "w", driver="core", backing_store=False) as fid:
            ds = fid.create_dataset("test", data=img)
            ds.attrs["geotransform"] = geobox.transform.to_gdal()
            ds.attrs["crs_wkt"] = geobox.crs.ExportToWkt()

            new_geobox = GriddedGeoBox.from_h5_dataset(ds)
            self.assertTrue(new_geobox.crs.ExportToWkt() == geobox.crs.ExportToWkt())

    def test_ggb_shape_from_h5_dataset(self):
        img, geobox = ut.create_test_image()
        with h5py.File("tmp.h5", "w", driver="core", backing_store=False) as fid:
            ds = fid.create_dataset("test", data=img)
            ds.attrs["geotransform"] = geobox.transform.to_gdal()
            ds.attrs["crs_wkt"] = geobox.crs.ExportToWkt()

            new_geobox = GriddedGeoBox.from_h5_dataset(ds)
            self.assertTrue(new_geobox.shape == img.shape)

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

    def test_project_extents(self):
        """
        Test that the geobox extents are correctly converted.
        """
        # values that have been pre-computed
        values = [
            148.53163652556793,
            -35.75238475020376,
            151.2210802166511,
            -33.501139639003675,
        ]
        truth = numpy.array(values)

        # set up CRS; WGS84 Geographics
        crs = osr.SpatialReference()
        crs.SetFromUserInput("EPSG:4326")

        # load the acquisition, get geobox, compute extents
        acq_cont = acquisitions(LS8_SCENE1)
        acq = acq_cont.get_acquisitions("RES-GROUP-1")[0]
        gb = acq.gridded_geo_box()
        exts = gb.project_extents(crs)

        self.assertIsNone(npt.assert_almost_equal(truth, exts, decimal=7))

    def test_lonlat_to_utm(self):
        """
        Test that coordinates in one CRS are correctly transformed
        from lonlat into utm
        This came about with GDAL3 and Proj6, where the native axis
        mapping of a CRS is respected, meaning that an x,y input and
        output is dependent on the CRS axis.
        We'll instead enforce the x,y axis mapping strategy.
        """
        shape = (3, 2)
        origin = (150.0, -34.0)
        ggb = GriddedGeoBox(shape, origin)

        lon = 148.862561
        lat = -35.123064

        to_crs = osr.SpatialReference()
        to_crs.ImportFromEPSG(32755)

        easting, northing = ggb.transform_coordinates((lon, lat), to_crs)

        self.assertAlmostEqual(easting, 669717.105361586)
        self.assertAlmostEqual(northing, 6111722.038508673)

    def test_utm_to_lonlat(self):
        """
        Test that coordinates in one CRS are correctly transformed
        from utm into lonlat
        This came about with GDAL3 and Proj6, where the native axis
        mapping of a CRS is respected, meaning that an x,y input and
        output is dependent on the CRS axis.
        We'll instead enforce the x,y axis mapping strategy.
        """
        shape = (3, 2)
        origin = (669700, 6111700)
        ggb = GriddedGeoBox(shape, origin, crs="EPSG:32755")

        lon = 148.862561
        lat = -35.123064
        easting = 669717.105361586
        northing = 6111722.038508673

        to_crs = osr.SpatialReference()
        to_crs.ImportFromEPSG(4326)

        lon, lat = ggb.transform_coordinates((easting, northing), to_crs)

        self.assertAlmostEqual(lon, 148.862561)
        self.assertAlmostEqual(lat, -35.123064)

    # def test_pixelscale_metres(self):
    #     scale = 0.00025
    #     shape = (4000, 4000)
    #     origin = (150.0, -34.0)
    #     ggb = GriddedGeoBox(shape, origin, pixelsize=(scale, scale))
    #     (size_x, size_y) = ggb.get_pixelsize_metres(xy=(0, 0))
    #     self.assertAlmostEqual(size_x, 23.0962, places=4)
    #     self.assertAlmostEqual(size_y, 27.7306, places=4)

    # def test_all_pixelscale_metres(self):
    #    scale = 0.00025
    #    shape = (4000, 4000)
    #    origin = (150.0, -34.0)
    #    ggb = GriddedGeoBox(shape, origin, pixelsize=(scale, scale))
    #    size_array = ggb.get_all_pixelsize_metres()
    #
    #    self.assertEqual(len(size_array), 4000)
    #    (size_x, size_y) = size_array[0]
    #    self.assertAlmostEqual(size_x, 23.0962, places=4)
    #    self.assertAlmostEqual(size_y, 27.7306, places=4)
    #    (size_x, size_y) = size_array[3999]
    #    self.assertAlmostEqual(size_x, 22.8221, places=4)
    #    self.assertAlmostEqual(size_y, 27.7351, places=4)


if __name__ == "__main__":
    unittest.main()
