from __future__ import absolute_import
from os.path import join as pjoin, abspath, dirname
import unittest
import ephem
from gaip.acquisition import acquisitions
from gaip.tle import load_tle

DATA_DIR = TLE_DIR = pjoin(dirname(abspath(__file__)), 'data')
LS5_SCENE1 = pjoin(DATA_DIR, 'LANDSAT5', 'LS5_TM_OTH_P51_GALPGS01-002_090_081_20090407')
LS7_SCENE1 = pjoin(DATA_DIR, 'LANDSAT7', 'LS7_ETM_OTH_P51_GALPGS01-002_090_081_20090415')
LS8_SCENE1 = pjoin(DATA_DIR, 'LANDSAT8', 'LS8_OLITIRS_OTH_P51_GALPGS01-032_090_084_20131011')


class TLELoadingTest(unittest.TestCase):

    def test_load_tle_l5_mtl1(self):
        acq = acquisitions(LS5_SCENE1).get_acquisitions()[0]
        data = load_tle(acq, TLE_DIR)
        self.assertIsInstance(data, ephem.EarthSatellite)

    def test_load_tle_l7_mtl1(self):
        acq = acquisitions(LS7_SCENE1).get_acquisitions()[0]
        data = load_tle(acq, TLE_DIR)
        self.assertIsInstance(data, ephem.EarthSatellite)

    def test_load_tle_l8_mtl1(self):
        acq = acquisitions(LS8_SCENE1).get_acquisitions()[0]
        data = load_tle(acq, TLE_DIR)
        self.assertIsInstance(data, ephem.EarthSatellite)


if __name__ == '__main__':
    unittest.main()
