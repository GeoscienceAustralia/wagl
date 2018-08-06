from __future__ import absolute_import
from os.path import join as pjoin, abspath, dirname
import unittest
import ephem

from wagl.acquisition import acquisitions
from wagl.tle import load_tle

from .data import LS5_SCENE1, LS7_SCENE1, LS8_SCENE1, TLE_DIR


class TLELoadingTest(unittest.TestCase):

    def test_load_tle_l5_mtl1(self):
        acq = acquisitions(LS5_SCENE1).get_acquisitions()[0]
        data = load_tle(acq, TLE_DIR)
        self.assertIsInstance(data, ephem.EarthSatellite)

    def test_load_tle_l7_mtl1(self):
        acq = acquisitions(LS7_SCENE1).get_acquisitions(group='RES-GROUP-1')[0]
        data = load_tle(acq, TLE_DIR)
        self.assertIsInstance(data, ephem.EarthSatellite)

    def test_load_tle_l8_mtl1(self):
        acq = acquisitions(LS8_SCENE1).get_acquisitions(group='RES-GROUP-1')[0]
        data = load_tle(acq, TLE_DIR)
        self.assertIsInstance(data, ephem.EarthSatellite)
