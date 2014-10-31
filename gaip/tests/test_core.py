import unittest
import gaip.core as core
import datetime
import os

MTL_FILE_DATA = os.path.join(os.path.dirname(__file__), 'data', 'mtl.txt')

class TypeParserTest(unittest.TestCase):

    def test_integer(self):
        num = core.parse_type('1')
        self.assertEqual(num, 1)

    def test_float(self):
        num = core.parse_type('1.0')
        self.assertEqual(num, 1.0)

    def test_datetime(self):
        dt0 = core.parse_type('2013-11-07T01:42:41Z')
        dt1 = datetime.datetime(2013, 11, 7, 1, 42, 41)
        self.assertEqual(dt0, dt1)

    def test_date(self):
        dt0 = core.parse_type('2013-11-07')
        dt1 = datetime.date(2013, 11, 7)
        self.assertEqual(dt0, dt1)

    def test_time(self):
        dt0 = core.parse_type('23:46:09.1442826Z')
        dt1 = datetime.time(23, 46, 9, 144282)
        self.assertEqual(dt0, dt1)

    def test_yes(self):
        resp = core.parse_type('Y')
        self.assertTrue(resp is True)

    def test_no(self):
        resp = core.parse_type('N')
        self.assertTrue(resp is False)

    def test_none(self):
        val = core.parse_type('NONE')
        self.assertIsNone(val)

    def test_str(self):
        s = core.parse_type('1adsd')
        self.assertEqual(s, '1adsd')



class MTLParserTest(unittest.TestCase):

    def test_load(self):
        tree = core.load_mtl(MTL_FILE_DATA)
        self.assertEqual(len(tree), 8)
        self.assertTrue(tree.has_key('METADATA_FILE_INFO'))
        self.assertTrue(tree.has_key('PRODUCT_METADATA'))
        self.assertTrue(tree.has_key('MIN_MAX_RADIANCE'))
        self.assertTrue(tree.has_key('MIN_MAX_PIXEL_VALUE'))
        self.assertTrue(tree.has_key('PRODUCT_PARAMETERS'))
        self.assertTrue(tree.has_key('CORRECTIONS_APPLIED'))
        self.assertTrue(tree.has_key('PROJECTION_PARAMETERS'))
        self.assertTrue(tree.has_key('UTM_PARAMETERS'))
        

class AcquisitionTest(unittest.TestCase):

    def test_load_acquisitions(self):
        acq = core.acquisitions(MTL_FILE_DATA)
        self.assertEqual(len(acq), 9)

    def test_acquisition(self):
        acq = core.acquisitions(MTL_FILE_DATA)[0]
        self.assertEqual(acq.band_name, 'band1')
        self.assertEqual(acq.band_num, 1)
        self.assertEqual(acq.file_name, 'L71090084_08420131003_B10.TIF')
        self.assertEqual(acq.lmin, -6.2)
        self.assertEqual(acq.lmax, 191.60)
        self.assertEqual(acq.qcalmin, 1.0)
        self.assertEqual(acq.qcalmax, 255.0)


if __name__ == '__main__':
    unittest.main()
