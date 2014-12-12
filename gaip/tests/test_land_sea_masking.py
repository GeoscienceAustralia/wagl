#!/bin/env python

import unittest
import rasterio as rio
import gaip
import os


class TestLandSea(unittest.TestCase):

    def test_L8(self):
        path = 'data/L1T/LS8_90_84_2013-10-11/UTM/LS8_OLITIRS_OTH_P51_GALPGS01-002_090_084_20131011'
        acqs = gaip.acquisitions(path)
        stack, geo_box = acqs[0].data_and_box()
        mask = gaip.calc_land_sea_mask(geo_box).astype('uint32')

        total_pixels = geo_box.shape[1]*geo_box.shape[0]
        land_pixels = sum(sum(mask))
        sea_pixels = total_pixels - land_pixels
        sea_pct = 100.0 * sea_pixels / total_pixels
        land_pct = 100.0 * land_pixels / total_pixels

        self.assertEqual(land_pixels, 75976737)
        self.assertEqual(sea_pixels, 9393744)
        self.assertEqual(total_pixels, 85370481)
        print "land=%f%%, sea=%f%%" % (land_pct, sea_pct)

#        gaip.write_img(mask.astype('uint8'), "./test.tif", format="GTiff", geobox=geo_box)

if __name__ == '__main__':
    # need this hack until set_epsilon patch v1.0.5
    # affine.set_epsilon(1e-9)
    import affine
    affine.EPSILON=1e-9
    affine.EPSILON2=1e-18

    unittest.main()
