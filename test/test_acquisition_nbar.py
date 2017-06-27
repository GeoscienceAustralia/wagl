from __future__ import absolute_import
import unittest
import gaip
import datetime
import rasterio
import os
from osgeo import osr

DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')

NBAR_DIR_LS5 = os.path.join(DATA_DIR, 'NBAR', '2009-01',
    'LS5_TM_NBAR_P54_GANBAR01-002_092_086_20090115')
NBAR_DIR_LS7 = os.path.join(DATA_DIR, 'NBAR', '2004-03',
    'LS7_ETM_NBAR_P54_GANBAR01-002_112_083_20040311')
NBAR_DIR_LS8 = os.path.join(DATA_DIR, 'NBAR', '2013-06',
    'LS8_OLI_TIRS_NBAR_P54_GANBAR01-032_090_082_20130605')


class AcquisitionTest(unittest.TestCase):

    def test_load_acquisitions(self):
        acqs = gaip.acquisitions(NBAR_DIR_LS5)
        self.assertEqual(len(acqs), 6)

    def test_acquisition(self):
        acq = gaip.acquisitions(NBAR_DIR_LS5)[0]

        self.assertTrue(isinstance(acq, gaip.Landsat5Acquisition))

        self.assertEqual(acq.band_type, gaip.REF)
        self.assertEqual(acq.band_name, 'band1')
        self.assertEqual(acq.band_num, 1)

        self.assertFalse(hasattr(acq, 'band1_sl_gain_change'))
        self.assertEqual(acq.sensor_id, 'TM')
        self.assertEquals(acq.no_data, -999, 'no_data value should be -999')      

    def test_acquisition_LS7(self):
        acqs = gaip.acquisitions(NBAR_DIR_LS7)
        self.assertEqual(len(acqs), 6)

        acq = acqs[0]
        self.assertTrue(isinstance(acq, gaip.Landsat7Acquisition))
        self.assertEqual(acq.band_type, gaip.REF)
        self.assertEqual(acq.band_name, 'band1')
        self.assertEqual(acq.band_num, 1)

        acq = acqs[5]
        self.assertEqual(acq.band_name, 'band7')
        self.assertEqual(acq.band_num, 7)
        self.assertEqual(acq.sensor_id, 'ETM+')
        self.assertEquals(acq.no_data, -999, 'no_data value should be -999')      

    def test_acquisition_LS8(self):
        acqs = gaip.acquisitions(NBAR_DIR_LS8)

        acq = acqs[0]
        self.assertTrue(isinstance(acq, gaip.Landsat8Acquisition))
        self.assertEqual(acq.band_type, gaip.REF)
        self.assertEqual(acq.band_name, 'band1')
        self.assertEqual(acq.band_num, 1)

        acq = acqs[6]
        self.assertEqual(acq.band_name, 'band7')
        self.assertEqual(acq.band_num, 7)
        self.assertEqual(acq.sensor_id, 'OLI_TIRS')
        self.assertEquals(acq.no_data, -999, 'no_data value should be -999')      

    def test_acquisition_LS5(self):
        acqs = gaip.acquisitions(NBAR_DIR_LS5)

        acq = acqs[0]
        self.assertTrue(isinstance(acq, gaip.Landsat5Acquisition))
        self.assertEqual(acq.band_type, gaip.REF)
        self.assertEqual(acq.band_name, 'band1')
        self.assertEqual(acq.band_num, 1)

        acq = acqs[1]
        self.assertEqual(acq.band_type, gaip.REF)
        self.assertEqual(acq.band_name, 'band2')
        self.assertEqual(acq.band_num, 2)
        self.assertEqual(acq.sensor_id, 'TM')
        self.assertEquals(acq.no_data, -999, 'no_data value should be -999')      

    def test_acquisition_LS8_no_data(self):
        acqs = gaip.acquisitions(NBAR_DIR_LS8)
        acq = acqs[0]
        self.assertTrue(isinstance(acq, gaip.Landsat8Acquisition))
        self.assertEqual(acq.band_type, gaip.REF)
        self.assertEquals(acq.no_data, -999, 'no_data value should be -999')      

    def test_single_band_read(self):
        acqs = gaip.acquisitions(NBAR_DIR_LS5)
        acq = acqs[0]
        band_data = acq.data()
        self.assertEqual(band_data.shape, (8561, 9641))

    def no_test_single_band_read_with_gridded_geo_box(self):
        acq = self.acqs[0]
        band_data, box = acq.data_and_box()
        self.assertEqual(band_data.shape, (9081, 9401))
        self.assertEqual(type(box), gaip.GriddedGeoBox)
        self.assertEqual(box.origin, (644000.0, 6283000.0))
        self.assertEqual(box.corner, (879025.0, 6055975.0))
        self.assertEqual(box.shape, (9081, 9401))
        self.assertEqual(box.pixelsize, (25.0, 25.0))

    def test_LS5_stack_read(self):
        acqs = gaip.acquisitions(NBAR_DIR_LS5)
        (acqs_read, stack, geobox) = gaip.stack_data(acqs, filter=(lambda acq: acq.band_type == gaip.REF))
        self.assertEqual(stack.shape, (6, 8561, 9641))



if __name__ == '__main__':
    unittest.main()
