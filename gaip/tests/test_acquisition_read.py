from __future__ import absolute_import
import unittest
import gaip
import datetime
import rasterio
import os
from osgeo import osr

DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')

L5_MTL = os.path.join(DATA_DIR, 'L5090081_08120090407_MTL.txt')
L5_DIR = os.path.join(DATA_DIR, 'L1T', 'LS5_90-84_1996-08-25', 'UTM',
                      'LS5_TM_OTH_P51_GALPGS01-002_090_084_19960825')

L7_MTL = os.path.join(DATA_DIR, 'L71090081_08120090415_MTL.txt')
L7_DIR = os.path.join(DATA_DIR, 'L1T', 'LS7_90-84_2000-09-13', 'UTM',
                      'LS7_ETM_OTH_P51_GALPGS01-002_090_084_20000913')

L8_MTL = os.path.join(DATA_DIR, 'LO80900842013284ASA00_MTL.txt')
L8_DIR = os.path.join(DATA_DIR, 'L1T', 'LS8_90_84_2013-10-11', 'UTM',
                      'LS8_OLITIRS_OTH_P51_GALPGS01-002_090_084_20131011')


class Landsat5AcquisitionTest(unittest.TestCase):

    def setUp(self):
        self.acqs = gaip.acquisitions(L5_DIR)

    def test_multi_band_read(self):
        acqs_subset, bands, geo_box = gaip.stack_data(self.acqs, \
            filter=(lambda acq: acq.band_type != gaip.PAN \
            and acq.band_type != gaip.THM))
        self.assertEqual(bands.shape, (6, 8801, 9721))

    def test_multi_band_read(self):
        acqs_subset, bands, geo_box = gaip.stack_data(self.acqs, \
            filter=(lambda acq: acq.band_type != gaip.PAN))
        self.assertEqual(bands.shape, (7, 8801, 9721))

    def test_read_with_no_acqs_selected(self):
        acqs_subset, bands, geo_box = gaip.stack_data(self.acqs, \
            filter=(lambda acq: False))
        self.assertEqual(len(acqs_subset), 0)
        self.assertTrue(bands == None)

    def test_read_with_no_acqs_input(self):
        acqs_subset, bands, geo_box = gaip.stack_data([], \
            filter=(lambda acq: False))
        self.assertEqual(len(acqs_subset), 0)
        self.assertTrue(bands == None)



class Landsat7AcquisitionTest(unittest.TestCase):

    def setUp(self):
        self.acqs = gaip.acquisitions(L7_DIR)

    def test_multi_band_read(self):
        acqs_subset, bands, geo_box = gaip.stack_data(self.acqs, \
            filter=(lambda acq: acq.band_type != gaip.PAN))
        self.assertEqual(bands.shape, (8, 8761, 9761))

    def test_multi_band_read(self):
        acqs_subset, bands, geo_box = gaip.stack_data(self.acqs, \
            filter=(lambda acq: acq.band_type != gaip.PAN))
        self.assertEqual(bands.shape, (8, 8761, 9761))

    def test_read_with_no_acqs_selected(self):
        acqs_subset, bands, geo_box = gaip.stack_data(self.acqs, \
            filter=(lambda acq: False))
        self.assertEqual(len(acqs_subset), 0)
        self.assertTrue(bands == None)

    def test_read_with_no_acqs_input(self):
        acqs_subset, bands, geo_box = gaip.stack_data([], \
            filter=(lambda acq: False))
        self.assertEqual(len(acqs_subset), 0)
        self.assertTrue(bands == None)

    def test_read_with_all_selected(self):
        try:
            acqs_subset, bands, geo_box = gaip.stack_data(self.acqs)
            self.fail("Should have got ValueError exception "
                "because Panchromatic band is wrong size" )
        except ValueError:
            pass


class Landsat8AcquisitionTest(unittest.TestCase):

    def setUp(self):
        self.acqs = gaip.acquisitions(L8_DIR)

    def test_gridded_geo_box(self):
        box = self.acqs[0].gridded_geo_box()
        self.assertEqual(type(box), gaip.GriddedGeoBox)
        self.assertEqual(box.origin, (644000.0, 6283000.0))
        self.assertEqual(box.corner, (879025.0, 6055975.0))
        self.assertEqual(box.shape, (9081, 9401))
        self.assertEqual(box.pixelsize, (25.0, 25.0))

    def test_gridded_geo_box_crs(self):
        box = self.acqs[0].gridded_geo_box()
        crs = box.crs
        self.assertEqual(type(crs), osr.SpatialReference)
        self.assertEqual(crs.IsProjected(), 1)
        self.assertEqual(crs.IsGeographic(), 0)
        self.assertEqual(crs.IsGeocentric(), 0)
        self.assertEqual(crs.GetUTMZone(), -55)
        self.assertEqual(crs.GetLinearUnits(), 1.0)
        self.assertEqual(crs.GetLinearUnitsName(), 'metre')
        self.assertAlmostEqual(crs.GetAngularUnits(), 0.0174532925199433)
        self.assertEqual(crs.ExportToProj4(), "+proj=utm +zone=55 +south +ellps=GRS80 " \
             "+towgs84=0,0,0,0,0,0,0 +units=m +no_defs " )

    def test_single_band_read(self):
        acq = self.acqs[0]
        band_data = acq.data()
        self.assertEqual(band_data.shape, (9081, 9401))

    def test_single_band_read_with_gridded_geo_box(self):
        acq = self.acqs[0]
        band_data, box = acq.data_and_box()
        self.assertEqual(band_data.shape, (9081, 9401))
        self.assertEqual(type(box), gaip.GriddedGeoBox)
        self.assertEqual(box.origin, (644000.0, 6283000.0))
        self.assertEqual(box.corner, (879025.0, 6055975.0))
        self.assertEqual(box.shape, (9081, 9401))
        self.assertEqual(box.pixelsize, (25.0, 25.0))

    def test_multi_band_read(self):
        acqs_subset, bands, geo_box = gaip.stack_data(self.acqs, \
            filter=(lambda acq: acq.band_type != gaip.PAN))
        self.assertEqual(bands.shape, (9, 9081, 9401))

    def test_multi_band_read(self):
        acqs_subset, bands, geo_box = gaip.stack_data(self.acqs, \
            filter=(lambda acq: acq.band_type != gaip.PAN))
        self.assertEqual(bands.shape, (9, 9081, 9401))

    def test_read_with_no_acqs_selected(self):
        acqs_subset, bands, geo_box = gaip.stack_data(self.acqs, \
            filter=(lambda acq: False))
        self.assertEqual(len(acqs_subset), 0)
        self.assertTrue(bands == None)
        self.assertTrue(geo_box == None)

    def test_read_with_no_acqs_input(self):
        acqs_subset, bands, geo_box = gaip.stack_data([], \
            filter=(lambda acq: False))
        self.assertEqual(len(acqs_subset), 0)
        self.assertTrue(bands == None)

    def test_read_with_all_selected(self):
        try:
            acqs_subset, bands, geo_box = gaip.stack_data(self.acqs)
            self.fail("Should have got ValueError exception "
                "because Panchromatic band is wrong size" )
        except ValueError:
            pass





class L1TDataTest(unittest.TestCase):

    def test_acquisition(self):
        self.assertTrue(os.path.exists(DATA_DIR))


if __name__ == '__main__':
    unittest.main()
