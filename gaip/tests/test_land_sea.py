#!/bin/env python

import unittest
import gaip
import os
import rasterio


class TestLandSea(unittest.TestCase):

    def test_flinders_islet(self):
        flindersOrigin = (150.927659, -34.453309)
        flindersCorner = (150.931697, -34.457915)

        # geobox with crs = WGS84
        flindersGGB = gaip.GriddedGeoBox.from_corners(flindersOrigin, flindersCorner)
        # print flindersGGB.affine

        # get the land sea mask for this small area

        mask = gaip.get_land_sea_mask(flindersGGB)
       
        total_pixels = mask.shape[0] * mask.shape[1] 
        land_pixels = sum(sum(mask.astype("uint32")))
        sea_pixels = total_pixels - land_pixels

        # print mask.shape
        self.assertEquals(total_pixels, 17*19)
        self.assertEquals(land_pixels, 182)
        self.assertEquals(sea_pixels, 141)
        # print mask

    def test_whole_AGDC_cell_no_sea(self):
        origin = (150.0, -33.0)
        corner = (151.0, -34.0)

        cellGGB = gaip.GriddedGeoBox.from_corners(origin, corner)
        mask = gaip.get_land_sea_mask(cellGGB)
        self.assertTrue(mask.shape == (4000,4000))
       
        total_pixels = mask.shape[0] * mask.shape[1] 
        land_pixels = sum(sum(mask.astype("uint32")))
        sea_pixels = total_pixels - land_pixels

        # print mask.shape
        self.assertEquals(total_pixels, 4000 * 4000)
        self.assertEquals(land_pixels, 4000*4000)
        self.assertEquals(sea_pixels, 0)
        # print mask

    def test_whole_AGDC_cell_land_and_sea(self):
        origin = (151.0, -33.0)
        corner = (152.0, -34.0)

        cellGGB = gaip.GriddedGeoBox.from_corners(origin, corner)
        mask = gaip.get_land_sea_mask(cellGGB)
        self.assertTrue(mask.shape == (4000,4000))
       
        total_pixels = mask.shape[0] * mask.shape[1] 
        land_pixels = sum(sum(mask.astype("uint32")))
        sea_pixels = total_pixels - land_pixels

        # print mask.shape
        self.assertEquals(total_pixels, 4000 * 4000)
        self.assertEquals(land_pixels, 6920586)
        self.assertEquals(sea_pixels, 9079414)
        # print mask

    def test_whole_AGDC_cell_all_sea(self):
        origin = (145.0, -39.0)
        corner = (146.0, -40.0)

        cellGGB = gaip.GriddedGeoBox.from_corners(origin, corner)
        mask = gaip.get_land_sea_mask(cellGGB)
        self.assertTrue(mask.shape == (4000,4000))
       
        total_pixels = mask.shape[0] * mask.shape[1] 
        land_pixels = sum(sum(mask.astype("uint32")))
        sea_pixels = total_pixels - land_pixels

        # print mask.shape
        self.assertEquals(total_pixels, 4000 * 4000)
        self.assertEquals(land_pixels, 0)
        self.assertEquals(sea_pixels, 4000 * 4000)
        # print mask
 
    def test_l1t_scene(self):
        p = 'data/NBAR/2013-06/LS8_OLI_TIRS_NBAR_P54_GANBAR01-032_103_070_20130617/scene01/LS8_OLI_TIRS_NBAR_P54_GANBAR01-032_103_070_20130617_B1.tif'
        with rasterio.open(p) as ds:
            geo_box = gaip.GriddedGeoBox.from_dataset(ds)
            mask = gaip.get_land_sea_mask(geo_box)
            self.assertTrue(mask.shape == (8801, 9041))
       
            total_pixels = mask.shape[0] * mask.shape[1] 
            land_pixels = sum(sum(mask.astype("uint32")))
            sea_pixels = total_pixels - land_pixels

            self.assertEquals(total_pixels, 8801 * 9041)
            self.assertEquals(land_pixels, 76137796)
            self.assertEquals(sea_pixels, 3432045)
            # print mask
 


if __name__ == '__main__':
    # need this hack until set_epsilon patch v1.0.5
    # affine.set_epsilon(1e-9)
    import affine
    affine.EPSILON=1e-9
    affine.EPSILON2=1e-18

    unittest.main()
