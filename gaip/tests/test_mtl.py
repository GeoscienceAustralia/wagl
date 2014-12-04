import unittest
import datetime
import gaip
import os

DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')

L5_MTL = os.path.join(DATA_DIR, 'L5090081_08120090407_MTL.txt')
L5_DIR = os.path.join(DATA_DIR, 'L1T', 'LS5_90-84_1996-08-25', 'UTM',
                      'LS5_TM_OTH_P51_GALPGS01-002_090_084_19960825')

L7_MTL = os.path.join(DATA_DIR, 'L71090081_08120090415_MTL.txt')
L7_DIR = os.path.join(DATA_DIR, 'L1T', 'LS7_90-84_2000-09-13', 'UTM',
                      'LS7_ETM_OTH_P51_GALPGS01-002_090_084_20000913')

L8_MTL = os.path.join(DATA_DIR, 'LO80900842013284ASA00_MTL.txt')
L8_DIR = os.path.join(DATA_DIR, 'L1T', 'LS8_90_84_2013-10-11', 'UTM',
                      'LS8_OLITIRS_OTH_P51_GALPGS01-002_090_084_20131011')


class TypeParserTest(unittest.TestCase):

    def test_integer(self):
        num = gaip.parse_type('1')
        self.assertEqual(num, 1)

    def test_float(self):
        num = gaip.parse_type('1.0')
        self.assertEqual(num, 1.0)

    def test_datetime(self):
        dt0 = gaip.parse_type('2013-11-07T01:42:41Z')
        dt1 = datetime.datetime(2013, 11, 7, 1, 42, 41)
        self.assertEqual(dt0, dt1)

    def test_date(self):
        dt0 = gaip.parse_type('2013-11-07')
        dt1 = datetime.date(2013, 11, 7)
        self.assertEqual(dt0, dt1)

    def test_time(self):
        dt0 = gaip.parse_type('23:46:09.1442826Z')
        dt1 = datetime.time(23, 46, 9, 144282)
        self.assertEqual(dt0, dt1)

    def test_yes(self):
        resp = gaip.parse_type('Y')
        self.assertTrue(resp is True)

    def test_no(self):
        resp = gaip.parse_type('N')
        self.assertTrue(resp is False)

    def test_none(self):
        val = gaip.parse_type('NONE')
        self.assertIsNone(val)

    def test_str(self):
        s = gaip.parse_type('1adsd')
        self.assertEqual(s, '1adsd')


class MTLParserTest(unittest.TestCase):

    def test_load(self):
        tree = gaip.load_mtl(L7_MTL)
        self.assertEqual(len(tree), 8)
        self.assertTrue(tree.has_key('METADATA_FILE_INFO'))
        self.assertTrue(tree.has_key('PRODUCT_METADATA'))
        self.assertTrue(tree.has_key('MIN_MAX_RADIANCE'))
        self.assertTrue(tree.has_key('MIN_MAX_PIXEL_VALUE'))
        self.assertTrue(tree.has_key('PRODUCT_PARAMETERS'))
        self.assertTrue(tree.has_key('CORRECTIONS_APPLIED'))
        self.assertTrue(tree.has_key('PROJECTION_PARAMETERS'))
        self.assertTrue(tree.has_key('UTM_PARAMETERS'))

if __name__ == '__main__':
    unittest.main()
