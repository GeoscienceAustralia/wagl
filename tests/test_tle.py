from __future__ import absolute_import
from os.path import join as pjoin, abspath, dirname
import unittest
import ephem
from gaip.acquisition import acquisitions
from gaip.tle import load_tle

DATA_DIR = TLE_DIR = pjoin(dirname(abspath(__file__)), 'data')
L5_MTL1 = pjoin(DATA_DIR, 'LANDSAT5', 'L5090081_08120090407_MTL.txt')
L5_MTL2 = pjoin(DATA_DIR, 'LANDSAT5', 'LT05_L1TP_095066_20100601_20170222_01_T1_MTL.txt')
L7_MTL1 = pjoin(DATA_DIR, 'LANDSAT7', 'L71090081_08120090415_MTL.txt')
L7_MTL2 = pjoin(DATA_DIR, 'LANDSAT7', 'LE07_L1TP_112066_20020218_20170221_01_T1_MTL.txt')
L8_MTL1 = pjoin(DATA_DIR, 'LANDSAT8', 'LO80900842013284ASA00_MTL.txt')
L8_MTL2 = pjoin(DATA_DIR, 'LANDSAT8', 'LC80990842016277LGN00_MTL.txt')


class TLELoadingTest(unittest.TestCase):

    def test_load_tle_l5_mtl1(self):
        acq = acquisitions(L5_MTL1).get_acquisitions()[0]
        data = load_tle(acq, TLE_DIR)
        self.assertIsInstance(data, ephem.EarthSatellite)

    def test_load_tle_l5_mtl2(self):
        acq = acquisitions(L5_MTL2).get_acquisitions()[0]
        data = load_tle(acq, TLE_DIR)
        self.assertIsInstance(data, ephem.EarthSatellite)

    def test_load_tle_l7_mtl1(self):
        acq = acquisitions(L7_MTL1).get_acquisitions()[0]
        data = load_tle(acq, TLE_DIR)
        self.assertIsInstance(data, ephem.EarthSatellite)

    def test_load_tle_l7_mtl2(self):
        acq = acquisitions(L7_MTL2).get_acquisitions()[0]
        data = load_tle(acq, TLE_DIR)
        self.assertIsInstance(data, ephem.EarthSatellite)

    def test_load_tle_l8_mtl1(self):
        acq = acquisitions(L8_MTL1).get_acquisitions()[0]
        data = load_tle(acq, TLE_DIR)
        self.assertIsInstance(data, ephem.EarthSatellite)

    def test_load_tle_l8_mtl2(self):
        acq = acquisitions(L8_MTL2).get_acquisitions()[0]
        data = load_tle(acq, TLE_DIR)
        self.assertIsInstance(data, ephem.EarthSatellite)


if __name__ == '__main__':
    unittest.main()
