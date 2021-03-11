from __future__ import absolute_import
from os.path import join as pjoin, abspath, dirname
import unittest
import datetime
from wagl.mtl import load_mtl, parse_type

DATA_DIR = pjoin(dirname(abspath(__file__)), "data")

L5_MTL1 = pjoin(DATA_DIR, "LANDSAT5", "L5090081_08120090407_MTL.txt")
L5_MTL2 = pjoin(DATA_DIR, "LANDSAT5", "LT05_L1TP_095066_20100601_20170222_01_T1_MTL.txt")
L7_MTL1 = pjoin(DATA_DIR, "LANDSAT7", "L71090081_08120090415_MTL.txt")
L7_MTL2 = pjoin(DATA_DIR, "LANDSAT7", "LE07_L1TP_112066_20020218_20170221_01_T1_MTL.txt")
L8_MTL1 = pjoin(DATA_DIR, "LANDSAT8", "LO80900842013284ASA00_MTL.txt")
L8_MTL2 = pjoin(DATA_DIR, "LANDSAT8", "LO80900842013284ASA00_MTL.txt")
L8_MTL1C2 = pjoin(
    DATA_DIR, "LANDSAT8", "LC08_L1TP_092084_20201029_20201106_02_T1_MTL.txt"
)
L8_MTRTC2 = pjoin(
    DATA_DIR,
    "LANDSAT8",
    "LC81060632021051LGN00__C2_RT",
    "LC08_L1TP_106063_20210220_20210220_02_RT_MTL.txt",
)


class TypeParserTest(unittest.TestCase):
    def test_integer(self):
        num = parse_type("1")
        self.assertEqual(num, 1)

    def test_float(self):
        num = parse_type("1.0")
        self.assertEqual(num, 1.0)

    def test_datetime(self):
        dt0 = parse_type("2013-11-07T01:42:41Z")
        dt1 = datetime.datetime(2013, 11, 7, 1, 42, 41)
        self.assertEqual(dt0, dt1)

    def test_quoted_datetime(self):
        dt0 = parse_type('"2013-11-07T01:42:41Z"')
        dt1 = datetime.datetime(2013, 11, 7, 1, 42, 41)
        self.assertEqual(dt0, dt1)

    def test_date(self):
        dt0 = parse_type("2013-11-07")
        dt1 = datetime.date(2013, 11, 7)
        self.assertEqual(dt0, dt1)

    def test_time(self):
        dt0 = parse_type("23:46:09.1442826Z")
        dt1 = datetime.time(23, 46, 9, 144282)
        self.assertEqual(dt0, dt1)

    def test_quoted_time(self):
        dt0 = parse_type('"23:46:09.1442826Z"')
        dt1 = datetime.time(23, 46, 9, 144282)
        self.assertEqual(dt0, dt1)

    def test_yes(self):
        resp = parse_type("Y")
        self.assertTrue(resp is True)

    def test_no(self):
        resp = parse_type("N")
        self.assertTrue(resp is False)

    def test_none(self):
        val = parse_type("NONE")
        self.assertIsNone(val)

    def test_str(self):
        s = parse_type("1adsd")
        self.assertEqual(s, "1adsd")


class Landsat5MTL1ParserTest(unittest.TestCase):
    def test_load(self):
        tree = load_mtl(L5_MTL1)
        self.assertEqual(len(tree), 9)
        self.assertTrue("METADATA_FILE_INFO" in tree)
        self.assertTrue("PRODUCT_METADATA" in tree)
        self.assertTrue("MIN_MAX_RADIANCE" in tree)
        self.assertTrue("MIN_MAX_PIXEL_VALUE" in tree)
        self.assertTrue("PRODUCT_PARAMETERS" in tree)
        self.assertTrue("CORRECTIONS_APPLIED" in tree)
        self.assertTrue("PROJECTION_PARAMETERS" in tree)
        self.assertTrue("UTM_PARAMETERS" in tree)


class Landsat5MTL2ParserTest(unittest.TestCase):
    def test_load(self):
        tree = load_mtl(L5_MTL2)
        self.assertEqual(len(tree), 10)
        self.assertTrue("METADATA_FILE_INFO" in tree)
        self.assertTrue("PRODUCT_METADATA" in tree)
        self.assertTrue("MIN_MAX_RADIANCE" in tree)
        self.assertTrue("MIN_MAX_REFLECTANCE" in tree)
        self.assertTrue("MIN_MAX_PIXEL_VALUE" in tree)
        self.assertTrue("PRODUCT_PARAMETERS" in tree)
        self.assertTrue("PROJECTION_PARAMETERS" in tree)
        self.assertTrue("IMAGE_ATTRIBUTES" in tree)
        self.assertTrue("THERMAL_CONSTANTS" in tree)


class Landsat7MTL1ParserTest(unittest.TestCase):
    def test_load(self):
        tree = load_mtl(L7_MTL1)
        self.assertEqual(len(tree), 8)
        self.assertTrue("METADATA_FILE_INFO" in tree)
        self.assertTrue("PRODUCT_METADATA" in tree)
        self.assertTrue("MIN_MAX_RADIANCE" in tree)
        self.assertTrue("MIN_MAX_PIXEL_VALUE" in tree)
        self.assertTrue("PRODUCT_PARAMETERS" in tree)
        self.assertTrue("CORRECTIONS_APPLIED" in tree)
        self.assertTrue("PROJECTION_PARAMETERS" in tree)
        self.assertTrue("UTM_PARAMETERS" in tree)


class Landsat7MTL2ParserTest(unittest.TestCase):
    def test_load(self):
        tree = load_mtl(L7_MTL2)
        self.assertEqual(len(tree), 10)
        self.assertTrue("METADATA_FILE_INFO" in tree)
        self.assertTrue("PRODUCT_METADATA" in tree)
        self.assertTrue("MIN_MAX_RADIANCE" in tree)
        self.assertTrue("MIN_MAX_REFLECTANCE" in tree)
        self.assertTrue("MIN_MAX_PIXEL_VALUE" in tree)
        self.assertTrue("PRODUCT_PARAMETERS" in tree)
        self.assertTrue("PROJECTION_PARAMETERS" in tree)
        self.assertTrue("IMAGE_ATTRIBUTES" in tree)
        self.assertTrue("THERMAL_CONSTANTS" in tree)


class Landsat8MTL1ParserTest(unittest.TestCase):
    def test_load(self):
        tree = load_mtl(L8_MTL1)
        self.assertEqual(len(tree), 9)
        self.assertTrue("METADATA_FILE_INFO" in tree)
        self.assertTrue("PRODUCT_METADATA" in tree)
        self.assertTrue("IMAGE_ATTRIBUTES" in tree)
        self.assertTrue("MIN_MAX_RADIANCE" in tree)
        self.assertTrue("MIN_MAX_REFLECTANCE" in tree)
        self.assertTrue("MIN_MAX_PIXEL_VALUE" in tree)
        self.assertTrue("RADIOMETRIC_RESCALING" in tree)
        self.assertTrue("TIRS_THERMAL_CONSTANTS" in tree)
        self.assertTrue("PROJECTION_PARAMETERS" in tree)


class Landsat8MTL1C2ParserTest(unittest.TestCase):
    def test_load(self):
        tree = load_mtl(L8_MTL1C2)
        self.assertEqual(len(tree), 10)
        self.assertTrue("PRODUCT_CONTENTS" in tree)  # was METADATA_FILE_INFO
        # self.assertTrue("PRODUCT_METADATA" in tree)
        self.assertTrue("IMAGE_ATTRIBUTES" in tree)  # PRODUCT_METADATA info added here
        self.assertTrue("LEVEL1_MIN_MAX_RADIANCE" in tree)
        self.assertTrue("LEVEL1_MIN_MAX_REFLECTANCE" in tree)
        self.assertTrue("LEVEL1_MIN_MAX_PIXEL_VALUE" in tree)
        self.assertTrue("LEVEL1_RADIOMETRIC_RESCALING" in tree)
        self.assertTrue("LEVEL1_THERMAL_CONSTANTS" in tree)
        self.assertTrue("LEVEL1_PROJECTION_PARAMETERS" in tree)


class Landsat8MTL1C2ParserTest(unittest.TestCase):
    def test_load(self):
        tree = load_mtl(L8_MTRTC2)
        self.assertEqual(len(tree), 10)
        self.assertTrue("PRODUCT_CONTENTS" in tree)
        self.assertTrue("IMAGE_ATTRIBUTES" in tree)
        self.assertTrue("LEVEL1_MIN_MAX_RADIANCE" in tree)
        self.assertTrue("LEVEL1_MIN_MAX_REFLECTANCE" in tree)
        self.assertTrue("LEVEL1_MIN_MAX_PIXEL_VALUE" in tree)
        self.assertTrue("LEVEL1_RADIOMETRIC_RESCALING" in tree)
        self.assertTrue("LEVEL1_THERMAL_CONSTANTS" in tree)
        self.assertTrue("LEVEL1_PROJECTION_PARAMETERS" in tree)


class Landsat8MTL2ParserTest(unittest.TestCase):
    def test_load(self):
        tree = load_mtl(L8_MTL2)
        self.assertEqual(len(tree), 9)
        self.assertTrue("METADATA_FILE_INFO" in tree)
        self.assertTrue("PRODUCT_METADATA" in tree)
        self.assertTrue("IMAGE_ATTRIBUTES" in tree)
        self.assertTrue("MIN_MAX_RADIANCE" in tree)
        self.assertTrue("MIN_MAX_REFLECTANCE" in tree)
        self.assertTrue("MIN_MAX_PIXEL_VALUE" in tree)
        self.assertTrue("RADIOMETRIC_RESCALING" in tree)
        self.assertTrue("TIRS_THERMAL_CONSTANTS" in tree)
        self.assertTrue("PROJECTION_PARAMETERS" in tree)


if __name__ == "__main__":
    unittest.main()
