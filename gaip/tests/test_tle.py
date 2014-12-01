import unittest
import gaip
import ephem
import os

DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')

L5_DIR = os.path.join(DATA_DIR, 'L1T', 'LS5_90-84_1996-08-25', 'UTM',
                      'LS5_TM_OTH_P51_GALPGS01-002_090_084_19960825')
L7_DIR = os.path.join(DATA_DIR, 'L1T', 'LS7_90-84_2000-09-13', 'UTM',
                      'LS7_ETM_OTH_P51_GALPGS01-002_090_084_20000913')
L8_DIR = os.path.join(DATA_DIR, 'L1T', 'LS8_90_84_2013-10-11', 'UTM',
                      'LS8_OLITIRS_OTH_P51_GALPGS01-002_090_084_20131011')

TLE_DIR = '/g/data/v10/eoancillarydata/ephemeris_kill'

class TLELoadingTest(unittest.TestCase):

    def test_load_tle_landsat5(self):
        acq = gaip.acquisitions(L5_DIR)[0]
        data = gaip.load_tle(acq, TLE_DIR)
        self.assertIsInstance(data, ephem.EarthSatellite)

    def test_load_tle_landsat7(self):
        acq = gaip.acquisitions(L7_DIR)[0]
        data = gaip.load_tle(acq, TLE_DIR)
        self.assertIsInstance(data, ephem.EarthSatellite)

    def test_load_tle_landsat8(self):
        acq = gaip.acquisitions(L8_DIR)[0]
        data = gaip.load_tle(acq, TLE_DIR)
        self.assertIsInstance(data, ephem.EarthSatellite)

if __name__ == '__main__':
    unittest.main()
