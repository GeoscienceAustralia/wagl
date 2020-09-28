#!/usr/bin/env python

"""
Test configurations of the tp5 file format for MODTRAN 5.4.
"""

from pathlib import Path
import unittest
from wagl import modtran_profiles as mp

FNAME1 = Path(__file__).parent.joinpath("data", "MIDLAT-SUMMER-ALBEDO.tp5")
FNAME2 = Path(__file__).parent.joinpath("data", "TROPICAL-ALBEDO.tp5")


class Tp5Test(unittest.TestCase):

    """
    Test configurations of the tp5 file format for MODTRAN 5.4 Point-Spread-Function Run.
    """

    def test_midlat_summer_albedo(self):
        """
        Test the mid latitude summer albedo configuration.
        """
        test = mp.MIDLAT_SUMMER_ALBEDO.format(
            aerosol_type=3,
            albedo=0.0,
            water=1.07000122070313,
            ozone=0.28499999642372,
            filter_function="landsat7_vsir.flt",
            visibility=-0.02264800000000,
            elevation=0.70900000000000,
            sat_height=705.0,
            sat_view=171.000748,
            doy=212,
            lat=-29.33856209871443,
            lon=209.88857485506449,
            time=23.73920805027778,
            sat_azimuth=279.408417,
        )

        with open(FNAME1, "r") as src:
            data = "".join(src.readlines())

        self.assertTrue(test == data)

    def test_tropical_albedo(self):
        """
        Test the tropical albedo configuration.
        """
        test = mp.TROPICAL_ALBEDO.format(
            aerosol_type=3,
            albedo=0.0,
            water=1.3500000000000001,
            ozone=0.25600001,
            filter_function="landsat8_vsir.flt",
            visibility=-0.043157435953617096,
            elevation=0.378,
            sat_height=705.0,
            sat_view=171.00043,
            doy=200,
            lat=-20.247597626228778,
            lon=229.23617402910139,
            time=1.2133087877777777,
            sat_azimuth=278.77069,
        )

        with open(FNAME2, "r") as src:
            data = "".join(src.readlines())

        self.assertTrue(test == data)


if __name__ == "__main__":
    unittest.main()
