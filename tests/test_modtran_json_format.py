#!/usr/bin/env python

"""
Test several different configurations of the json input  file format.

"""

from os.path import join as pjoin, abspath, dirname
import unittest
from wagl import modtran_profile_json as  mpjson
import json
from wagl.modtran import JsonEncoder


DATA_DIR = pjoin(dirname(abspath(__file__)), 'data')
FNAME1 = pjoin(DATA_DIR, 'midlat_summber_albedo.json')
FNAME2 = pjoin(DATA_DIR, 'tropical_albedo.json')
FNAME3 = pjoin(DATA_DIR, 'thermal_channel.json')


class JsonTest(unittest.TestCase):

    """
    Test different configurations of the json file format,
    and that the input/output round trips.
    """

    def test_midlat_summer_albedo(self):

        """
        Test the mid latitude summer albedo configuration.
        """
        test = mpjson.midlat_summer_albedo(name="POINT-0-ALBEDO-0",
                                           water=1.85,
                                           ozone=0.025,
                                           visibility=-0.04,
                                           doy =97,
                                           lats=-29.0,
                                           lons=209.0,
                                           time=23.0,
                                           sat_azimuth=279.5,
                                           elevation=0.292,
                                           sat_height=705.0,
                                           sat_view=171.06,
                                           albedo=0.0,
                                           filter_function="landsat5_vsir.flt",
                                           binary=False)

        d = json.dumps(test, cls=JsonEncoder, indent=4)

        with open(FNAME1, 'r') as src:
            data = ''.join(src.readlines())

        self.assertTrue(d == data)

    def test_tropical_albedo(self):

        """
        Test the tropical  albedo configuration.
        """

        test = mpjson.tropical_albedo(name="POINT-0-ALBEDO-0",
                                           water=1.85,
                                           ozone=0.025,
                                           visibility=-0.04,
                                           doy=97,
                                           lats=-29.0,
                                           lons=209.0,
                                           time=23.0,
                                           sat_azimuth=279.5,
                                           elevation=0.292,
                                           sat_height=705.0,
                                           sat_view=171.06,
                                           albedo=0.0,
                                           filter_function="landsat5_vsir.flt",
                                           binary=False)

        d = json.dumps(test, cls=JsonEncoder, indent=4)

        with open(FNAME2, 'r') as src:
            data = ''.join(src.readlines())

        self.assertTrue(d == data)

    def test_thermal_channel(self):

        """
        Test surface brightness  configuration.
        """
        pass


if __name__ == '__main__':

    unittest.main()
