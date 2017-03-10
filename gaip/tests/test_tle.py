from __future__ import absolute_import
import unittest
import gaip
import ephem
import os

DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')
L5_MTL1 = os.path.join(DATA_DIR, 'L5090081_08120090407_MTL.txt')
L5_MTL2 = os.path.join(DATA_DIR, 'LT05_L1TP_095066_20100601_20170222_01_T1_MTL.txt')
L7_MTL1 = os.path.join(DATA_DIR, 'L71090081_08120090415_MTL.txt')
L7_MTL2 = os.path.join(DATA_DIR, 'LE07_L1TP_112066_20020218_20170221_01_T1_MTL.txt')
L8_MTL1 = os.path.join(DATA_DIR, 'LO80900842013284ASA00_MTL.txt')
L8_MTL2 = os.path.join(DATA_DIR, 'LO80900842013284ASA00_MTL.txt')

TLE_DIR = '/g/data/v10/eoancillarydata/sensor-specific'

class TLELoadingTest(unittest.TestCase):

    def test_load_tle_l5_mtl1(self):
        acq = gaip.acquisitions(L5_MTL1).get_acquisitions()[0]
        data = gaip.load_tle(acq, TLE_DIR)
        self.assertIsInstance(data, ephem.EarthSatellite)

    def test_load_tle_l5_mtl2(self):
        acq = gaip.acquisitions(L5_MTL2).get_acquisitions()[0]
        data = gaip.load_tle(acq, TLE_DIR)
        self.assertIsInstance(data, ephem.EarthSatellite)

    def test_load_tle_l7_mtl1(self):
        acq = gaip.acquisitions(L7_MTL1).get_acquisitions()[0]
        data = gaip.load_tle(acq, TLE_DIR)
        self.assertIsInstance(data, ephem.EarthSatellite)

    def test_load_tle_l7_mtl2(self):
        acq = gaip.acquisitions(L7_MTL2).get_acquisitions()[0]
        data = gaip.load_tle(acq, TLE_DIR)
        self.assertIsInstance(data, ephem.EarthSatellite)

    def test_load_tle_l8_mtl1(self):
        acq = gaip.acquisitions(L8_MTL1).get_acquisitions()[0]
        data = gaip.load_tle(acq, TLE_DIR)
        self.assertIsInstance(data, ephem.EarthSatellite)

    def test_load_tle_l8_mtl2(self):
        acq = gaip.acquisitions(L8_MTL2).get_acquisitions()[0]
        data = gaip.load_tle(acq, TLE_DIR)
        self.assertIsInstance(data, ephem.EarthSatellite)


if __name__ == '__main__':
    unittest.main()
