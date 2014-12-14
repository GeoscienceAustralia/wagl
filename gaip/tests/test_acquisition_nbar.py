import unittest
import gaip
import datetime
import rasterio
import os
from osgeo import osr

DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')

NBAR_DIR_LS5 = os.path.join(DATA_DIR, 'NBAR', '2009-01',
    'LS5_TM_NBAR_P54_GANBAR01-002_092_086_20090115')
NBAR_DIR_LS8 = os.path.join(DATA_DIR, 'NBAR', '2013-06',
    'LS8_OLI_TIRS_NBAR_P54_GANBAR01-032_090_082_20130605')


class AcquisitionTest(unittest.TestCase):

    def test_load_acquisitions(self):
        acqs = gaip.acquisitions(NBAR_DIR_LS5)
        self.assertEqual(len(acqs), 6)

    def test_acquisition(self):
        acq = gaip.acquisitions(NBAR_DIR_LS5)[0]

        self.assertEqual(acq.band_name, 'band1')
        self.assertEqual(acq.band_num, 1)

        self.assertFalse(hasattr(acq, 'band1_sl_gain_change'))
        self.assertEqual(acq.sensor_id, 'TM')

    def test_acquisition_LS8(self):
        acqs = gaip.acquisitions(NBAR_DIR_LS8)

        acq = acqs[0]
        self.assertEqual(acq.band_name, 'band1')
        self.assertEqual(acq.band_num, 1)

        acq = acqs[6]
        self.assertEqual(acq.band_name, 'band7')
        self.assertEqual(acq.band_num, 7)
        self.assertEqual(acq.sensor_id, 'OLI_TIRS')

if __name__ == '__main__':
    unittest.main()
