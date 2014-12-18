import unittest
import gaip
import datetime
import rasterio
import os
from osgeo import osr

L8_DIR = './data/luigi_test/LS8_OLITIRS_OTH_P51_GALPGS01-002_115_075_20141014'

class Landsat8AcquisitionTest(unittest.TestCase):

    def setUp(self):
        self.acqs = gaip.acquisitions(L8_DIR)

    def test_discovery(self):
        acqs = gaip.acquisitions(L8_DIR)
        self.assertEqual(len(acqs), 12)

    def test_type(self):
        for acq in self.acqs:
            self.assertTrue(isinstance(acq, gaip.Landsat8Acquisition))

    def test_lines(self):
        for acq in self.acqs:
            filename = os.path.join(acq.dir_name, acq.file_name)
            with rasterio.open(filename) as a:
                self.assertEqual(acq.lines, a.height)

    def test_samples(self):
        for acq in self.acqs:
            filename = os.path.join(acq.dir_name, acq.file_name)
            with rasterio.open(filename) as a:
                self.assertEqual(acq.samples, a.width)

    def test_scene_center_time(self):
        for acq in self.acqs:
            self.assertEqual(acq.scene_center_time,
                datetime.time(2, 21, 25, 902494))

    def test_read_with_all_selected(self):
        try:
            acqs_subset, bands, geo_box = gaip.stack_data(self.acqs)
            print acqs_subset
            print bands.shape
            self.fail("Should have got ValueError exception "
                "because Panchromatic band is wrong size" )
        except ValueError:
            pass

    def test_band_type(self):
        self.assertEqual(self.acqs[0].band_type, gaip.REF)
        self.assertEqual(self.acqs[1].band_type, gaip.REF)
        self.assertEqual(self.acqs[2].band_type, gaip.REF)
        self.assertEqual(self.acqs[3].band_type, gaip.REF)
        self.assertEqual(self.acqs[4].band_type, gaip.REF)
        self.assertEqual(self.acqs[5].band_type, gaip.REF)
        self.assertEqual(self.acqs[6].band_type, gaip.REF)
        self.assertEqual(self.acqs[7].band_type, gaip.PAN)
        self.assertEqual(self.acqs[8].band_type, gaip.ATM)
        self.assertEqual(self.acqs[9].band_type, gaip.THM)
        self.assertEqual(self.acqs[10].band_type, gaip.THM)
        self.assertEqual(self.acqs[11].band_type, gaip.BQA)

    def test_grid_cell_size(self):
        self.assertEqual(self.acqs[0].grid_cell_size, 25.0)
        self.assertEqual(self.acqs[1].grid_cell_size, 25.0)
        self.assertEqual(self.acqs[2].grid_cell_size, 25.0)
        self.assertEqual(self.acqs[3].grid_cell_size, 25.0)
        self.assertEqual(self.acqs[4].grid_cell_size, 25.0)
        self.assertEqual(self.acqs[5].grid_cell_size, 25.0)
        self.assertEqual(self.acqs[6].grid_cell_size, 25.0)
        self.assertEqual(self.acqs[7].grid_cell_size, 12.5)
        self.assertEqual(self.acqs[8].grid_cell_size, 25.0)
        self.assertEqual(self.acqs[9].grid_cell_size, 25.0)
        self.assertEqual(self.acqs[10].grid_cell_size, 25.0)
        self.assertEqual(self.acqs[11].grid_cell_size, 25.0)

    def test_scene_center_datetime(self):
        for acq in self.acqs:
            self.assertEqual(acq.scene_center_datetime,
                datetime.datetime(2014, 10, 14, 2, 21, 25, 902494))

    def no_test_min_max_radiance_band1(self):
        self.assertEqual(self.acqs[0].min_radiance, -64.75256)
        self.assertEqual(self.acqs[0].max_radiance, 784.11609)

    def no_test_min_max_radiance_band2(self):
        self.assertEqual(self.acqs[0].min_radiance, -64.75256)
        self.assertEqual(self.acqs[0].max_radiance, 784.11609)

    def no_test_min_max_radiance_band3(self):
        self.assertEqual(self.acqs[0].min_radiance, -64.75256)
        self.assertEqual(self.acqs[0].max_radiance, 784.11609)

    def no_test_lmin(self):
        self.assertEqual(self.acqs[0].lmin, -64.75256)

    def no_test_lmax(self):
        self.assertEqual(self.acqs[0].lmax, 784.11609)

    def no_test_qcalmin(self):
        self.assertEqual(self.acqs[0].qcalmin, 1)

    def no_test_qcalmax(self):
        self.assertEqual(self.acqs[0].qcalmax, 65535)

    def no_test_zone_number(self):
        self.assertEqual(self.acqs[0].zone_number, -55)

    def no_test_sun_azimuth(self):
        self.assertEqual(self.acqs[0].sun_azimuth, 50.86088724)

    def no_test_sun_elevation(self):
        self.assertEqual(self.acqs[0].sun_elevation, 52.25003864)

    def no_test_gain(self):
        self.assertAlmostEquals(self.acqs[0].gain, 0.012953)

    def no_test_bias(self):
        self.assertAlmostEquals(self.acqs[0].bias, -64.76551)



if __name__ == '__main__':
    unittest.main()
