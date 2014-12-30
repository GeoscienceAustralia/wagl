#!/bin/env python

import unittest
import rasterio as rio
import gaip
import os


class TestLandSea(unittest.TestCase):

    def no_test_L8(self):
        path = 'data/L1T/LS8_90_84_2013-10-11/UTM/LS8_OLITIRS_OTH_P51_GALPGS01-002_090_084_20131011'
        acqs = gaip.acquisitions(path)
        stack, geo_box = acqs[0].data_and_box()
        print geo_box
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
 
    def no_test_l1t_scene(self):
        p = '/g/data1/v10/projects/Luigi_work_flow_test/L1T/LS8_OLITIRS_OTH_P51_GALPGS01-002_115_075_20141014/scene01/LC81150752014287ASA00_B1.TIF'
        with rio.open(p) as ds:
            geo_box = gaip.GriddedGeoBox.from_dataset(ds)
            print geo_box
            mask = gaip.calc_land_sea_mask(geo_box).astype('uint32')

            total_pixels = geo_box.shape[1]*geo_box.shape[0]
            land_pixels = sum(sum(mask))
            sea_pixels = total_pixels - land_pixels
            sea_pct = 100.0 * sea_pixels / total_pixels
            land_pct = 100.0 * land_pixels / total_pixels

#            self.assertEqual(land_pixels, 75976737)
#            self.assertEqual(sea_pixels, 9393744)
#            self.assertEqual(total_pixels, 85370481)
            print "land=%f%%, sea=%f%%" % (land_pct, sea_pct)

    def no_test_agdc_cell(self):
        scale = 0.00025
        shape = (4000,4000)
        origin = (150.0, -34.0)

        corner = (shape[1]*scale+150, -34-shape[0]*scale)
        ggb = gaip.GriddedGeoBox(shape, origin, pixelsize=(scale, scale))

        # now get UTM equilavent

        geo_box = ggb.copy(crs="EPSG:32756")
        print geo_box


        # and get the mask

        mask = gaip.calc_land_sea_mask(geo_box).astype('uint32')

        total_pixels = geo_box.shape[1]*geo_box.shape[0]
        land_pixels = sum(sum(mask))
        sea_pixels = total_pixels - land_pixels
        sea_pct = 100.0 * sea_pixels / total_pixels
        land_pct = 100.0 * land_pixels / total_pixels

        print "land_pixels=%d" % land_pixels
        print "sea_pixels=%d" % sea_pixels
        self.assertEqual(land_pixels, 14554858)
        self.assertEqual(sea_pixels, 1445142)
        self.assertEqual(total_pixels, 16000000)

        print "land=%f%%, sea=%f%%" % (land_pct, sea_pct)

    def test_agdc_cell_darwin(self):
        scale = 0.00025
        shape = (4000,4000)
        origin = (130.0, -18.0)

        corner = (shape[1]*scale+origin[0], origin[1]-shape[0]*scale)
        ggb = gaip.GriddedGeoBox(shape, origin, pixelsize=(scale, scale))
        print ggb

        # now get UTM equilavent

        geo_box = ggb.copy(crs="EPSG:32752")
        print geo_box


        # and get the mask

        mask = gaip.calc_land_sea_mask(geo_box).astype('uint32')

        total_pixels = geo_box.shape[1]*geo_box.shape[0]
        land_pixels = sum(sum(mask))
        sea_pixels = total_pixels - land_pixels
        sea_pct = 100.0 * sea_pixels / total_pixels
        land_pct = 100.0 * land_pixels / total_pixels

        print "geobox_shape=%s" % str(geo_box.shape)
        print "mask_shape=%s" % str(mask.shape)
        print "total_pixels=%d" % total_pixels
        print "land_pixels=%d" % land_pixels
        print "sea_pixels=%d" % sea_pixels
        # self.assertEqual(land_pixels, 14554858)
        # self.assertEqual(sea_pixels, 1445142)
        # self.assertEqual(total_pixels, 16000000)

        print "land=%f%%, sea=%f%%" % (land_pct, sea_pct)



if __name__ == '__main__':
    # need this hack until set_epsilon patch v1.0.5
    # affine.set_epsilon(1e-9)
    import affine
    affine.EPSILON=1e-9
    affine.EPSILON2=1e-18

    unittest.main()
