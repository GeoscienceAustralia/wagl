import unittest
import gaip
import os

DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')

L5_DIR = os.path.join(DATA_DIR, 'L1T', 'LS5_90-84_1996-08-25', 'UTM',
                      'LS5_TM_OTH_P51_GALPGS01-002_090_084_19960825')

L7_DIR = os.path.join(DATA_DIR, 'L1T', 'LS7_90-84_2000-09-13', 'UTM',
                      'LS7_ETM_OTH_P51_GALPGS01-002_090_084_20000913')

L8_DIR = os.path.join(DATA_DIR, 'L1T', 'LS8_90_84_2013-10-11', 'UTM',
                      'LS8_OLITIRS_OTH_P51_GALPGS01-002_090_084_20131011')


class Landsat5DataTest(unittest.TestCase):

    def setUp(self):
        self.acqs = gaip.acquisitions(L5_DIR)

    def test_load_band(self):
        acqs = self.acqs
        for acq in acqs:
            band_data = gaip.data(acq)
            w, h = band_data.shape
            self.assertEqual(h, acq.samples)
            self.assertEqual(w, acq.lines)

    def test_load_band_alternative(self):
        acqs = self.acqs
        for acq in acqs:
            band_data = acq.data()
            w, h = band_data.shape
            self.assertEqual(h, acq.samples)
            self.assertEqual(w, acq.lines)


class Landsat7DataTest(unittest.TestCase):

    def setUp(self):
        self.acqs = gaip.acquisitions(L7_DIR)

    def test_load_band(self):
        acqs = self.acqs
        for acq in acqs:
            band_data = gaip.data(acq)
            w, h = band_data.shape
            self.assertEqual(h, acq.samples)
            self.assertEqual(w, acq.lines)

    def test_load_band_alternative(self):
        acqs = self.acqs
        for acq in acqs:
            band_data = acq.data()
            w, h = band_data.shape
            self.assertEqual(h, acq.samples)
            self.assertEqual(w, acq.lines)


class Landsat8DataTest(unittest.TestCase):

    def setUp(self):
        self.acqs = gaip.acquisitions(L8_DIR)

    def test_load_band(self):
        acqs = self.acqs
        for acq in acqs:
            band_data = gaip.data(acq)
            w, h = band_data.shape
            self.assertEqual(h, acq.samples)
            self.assertEqual(w, acq.lines)

    def test_load_band_alternative(self):
        acqs = self.acqs
        for acq in acqs:
            band_data = acq.data()
            w, h = band_data.shape
            self.assertEqual(h, acq.samples)
            self.assertEqual(w, acq.lines)


if __name__ == '__main__':
    unittest.main()
