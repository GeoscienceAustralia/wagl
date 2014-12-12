import unittest
import gaip
import datetime
import rasterio
import os
from osgeo import osr

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


class AcquisitionTest(unittest.TestCase):

    def test_load_acquisitions(self):
        acq = gaip.acquisitions(L7_MTL)
        self.assertEqual(len(acq), 9)

    def test_acquisition(self):
        acq = gaip.acquisitions(L7_MTL)[0]

        self.assertEqual(acq.band_name, 'band1')
        self.assertEqual(acq.band_num, 1)

        self.assertFalse(hasattr(acq, 'band1_sl_gain_change'))


class Landsat5AcquisitionTest(unittest.TestCase):

    def setUp(self):
        self.acqs = gaip.acquisitions(L5_DIR)

    def test_discovery(self):
        acqs = gaip.acquisitions(L5_DIR)
        self.assertEqual(len(acqs), 7)

    def test_type(self):
        for acq in self.acqs:
            self.assertTrue(isinstance(acq, gaip.LandsatAcquisition))

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

    def test_band_type(self):
        self.assertEqual(self.acqs[0].band_type, gaip.REF)
        self.assertEqual(self.acqs[1].band_type, gaip.REF)
        self.assertEqual(self.acqs[2].band_type, gaip.REF)
        self.assertEqual(self.acqs[3].band_type, gaip.REF)
        self.assertEqual(self.acqs[4].band_type, gaip.REF)
        self.assertEqual(self.acqs[5].band_type, gaip.THM)
        self.assertEqual(self.acqs[6].band_type, gaip.REF)

    def test_grid_cell_size(self):
        self.assertEqual(self.acqs[0].grid_cell_size, 25.0)
        self.assertEqual(self.acqs[1].grid_cell_size, 25.0)
        self.assertEqual(self.acqs[2].grid_cell_size, 25.0)
        self.assertEqual(self.acqs[3].grid_cell_size, 25.0)
        self.assertEqual(self.acqs[4].grid_cell_size, 25.0)
        self.assertEqual(self.acqs[5].grid_cell_size, 25.0)
        self.assertEqual(self.acqs[6].grid_cell_size, 25.0)

    def test_scene_center_time(self):
        for acq in self.acqs:
            self.assertEqual(acq.scene_center_time,
                             datetime.time(23, 6, 53, 583075))

    def test_scene_center_datetime(self):
        for acq in self.acqs:
            self.assertEqual(acq.scene_center_datetime,
                             datetime.datetime(1996, 8, 25, 23, 6, 53, 583075))

    def test_min_max_radiance_band1(self):
        self.assertEqual(self.acqs[0].min_radiance, -1.520)
        self.assertEqual(self.acqs[0].max_radiance, 193.0)

    def test_min_max_radiance_band2(self):
        self.assertEqual(self.acqs[0].min_radiance, -1.52)
        self.assertEqual(self.acqs[0].max_radiance, 193.0)

    def test_min_max_radiance_band3(self):
        self.assertEqual(self.acqs[0].min_radiance, -1.52)
        self.assertEqual(self.acqs[0].max_radiance, 193.0)

    def test_lmin(self):
        self.assertEqual(self.acqs[0].lmin, -1.520)

    def test_lmax(self):
        self.assertEqual(self.acqs[0].lmax, 193.0)

    def test_qcalmin(self):
        self.assertEqual(self.acqs[0].qcalmin, 1.0)

    def test_qcalmax(self):
        self.assertEqual(self.acqs[0].qcalmax, 255.0)

    def test_zone_number(self):
        self.assertEqual(self.acqs[0].zone_number, -55)

    def test_sun_azimuth(self):
        self.assertEqual(self.acqs[0].sun_azimuth, 51.1073225)

    def test_sun_elevation(self):
        self.assertEqual(self.acqs[0].sun_elevation, 28.704069)

    def test_gain(self):
        self.assertAlmostEquals(self.acqs[0].gain, 0.7658267716535433)

    def test_bias(self):
        self.assertAlmostEquals(self.acqs[0].bias, -2.2858267716535465)


class Landsat7AcquisitionTest(unittest.TestCase):

    def setUp(self):
        self.acqs = gaip.acquisitions(L7_DIR)

    def test_discovery(self):
        acqs = gaip.acquisitions(L7_DIR)
        self.assertEqual(len(acqs), 9)

    def test_type(self):
        for acq in self.acqs:
            self.assertTrue(isinstance(acq, gaip.LandsatAcquisition))

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

    def test_band_type(self):
        self.assertEqual(self.acqs[0].band_type, gaip.REF)
        self.assertEqual(self.acqs[1].band_type, gaip.REF)
        self.assertEqual(self.acqs[2].band_type, gaip.REF)
        self.assertEqual(self.acqs[3].band_type, gaip.REF)
        self.assertEqual(self.acqs[4].band_type, gaip.REF)
        self.assertEqual(self.acqs[5].band_type, gaip.THM)
        self.assertEqual(self.acqs[6].band_type, gaip.THM)
        self.assertEqual(self.acqs[7].band_type, gaip.REF)
        self.assertEqual(self.acqs[8].band_type, gaip.PAN)

    def test_grid_cell_size(self):
        self.assertEqual(self.acqs[0].grid_cell_size, 25.0)
        self.assertEqual(self.acqs[1].grid_cell_size, 25.0)
        self.assertEqual(self.acqs[2].grid_cell_size, 25.0)
        self.assertEqual(self.acqs[3].grid_cell_size, 25.0)
        self.assertEqual(self.acqs[4].grid_cell_size, 25.0)
        self.assertEqual(self.acqs[5].grid_cell_size, 25.0)
        self.assertEqual(self.acqs[6].grid_cell_size, 25.0)
        self.assertEqual(self.acqs[7].grid_cell_size, 25.0)
        self.assertEqual(self.acqs[8].grid_cell_size, 12.5)

    def test_scene_center_time(self):
        for acq in self.acqs:
            self.assertEqual(acq.scene_center_time,
                             datetime.time(23, 40, 55, 160927))

    def test_scene_center_datetime(self):
        for acq in self.acqs:
            self.assertEqual(acq.scene_center_datetime,
                             datetime.datetime(2000, 9, 13, 23, 40, 55, 160927))

    def test_min_max_radiance_band1(self):
        self.assertEqual(self.acqs[0].min_radiance, -6.2)
        self.assertEqual(self.acqs[0].max_radiance, 191.6)

    def test_min_max_radiance_band2(self):
        self.assertEqual(self.acqs[0].min_radiance, -6.2)
        self.assertEqual(self.acqs[0].max_radiance, 191.6)

    def test_min_max_radiance_band3(self):
        self.assertEqual(self.acqs[0].min_radiance, -6.2)
        self.assertEqual(self.acqs[0].max_radiance, 191.6)

    def test_lmin(self):
        self.assertEqual(self.acqs[0].lmin, -6.2)

    def test_lmax(self):
        self.assertEqual(self.acqs[0].lmax, 191.6)

    def test_qcalmin(self):
        self.assertEqual(self.acqs[0].qcalmin, 1.0)

    def test_qcalmax(self):
        self.assertEqual(self.acqs[0].qcalmax, 255.0)

    def test_zone_number(self):
        self.assertEqual(self.acqs[0].zone_number, -55)

    def test_sun_azimuth(self):
        self.assertEqual(self.acqs[0].sun_azimuth, 46.9270726)

    def test_sun_elevation(self):
        self.assertEqual(self.acqs[0].sun_elevation, 40.4612407)

    def test_gain(self):
        self.assertAlmostEquals(self.acqs[0].gain, 0.7787401574803149)

    def test_bias(self):
        self.assertAlmostEquals(self.acqs[0].bias, -6.978740157480303)




class Landsat8AcquisitionTest(unittest.TestCase):

    def setUp(self):
        self.acqs = gaip.acquisitions(L8_DIR)

    def test_discovery(self):
        acqs = gaip.acquisitions(L8_DIR)
        self.assertEqual(len(acqs), 10)

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
                             datetime.time(23, 52, 10, 108347))

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
        self.assertEqual(self.acqs[9].band_type, gaip.BQA)

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

    def test_scene_center_datetime(self):
        for acq in self.acqs:
            self.assertEqual(acq.scene_center_datetime,
                             datetime.datetime(2013, 10, 11, 23, 52, 10,
                                               108347))

    def test_min_max_radiance_band1(self):
        self.assertEqual(self.acqs[0].min_radiance, -64.75256)
        self.assertEqual(self.acqs[0].max_radiance, 784.11609)

    def test_min_max_radiance_band2(self):
        self.assertEqual(self.acqs[0].min_radiance, -64.75256)
        self.assertEqual(self.acqs[0].max_radiance, 784.11609)

    def test_min_max_radiance_band3(self):
        self.assertEqual(self.acqs[0].min_radiance, -64.75256)
        self.assertEqual(self.acqs[0].max_radiance, 784.11609)

    def test_lmin(self):
        self.assertEqual(self.acqs[0].lmin, -64.75256)

    def test_lmax(self):
        self.assertEqual(self.acqs[0].lmax, 784.11609)

    def test_qcalmin(self):
        self.assertEqual(self.acqs[0].qcalmin, 1)

    def test_qcalmax(self):
        self.assertEqual(self.acqs[0].qcalmax, 65535)

    def test_zone_number(self):
        self.assertEqual(self.acqs[0].zone_number, -55)

    def test_sun_azimuth(self):
        self.assertEqual(self.acqs[0].sun_azimuth, 50.86088724)

    def test_sun_elevation(self):
        self.assertEqual(self.acqs[0].sun_elevation, 52.25003864)

    def test_gain(self):
        self.assertAlmostEquals(self.acqs[0].gain, 0.012953)

    def test_bias(self):
        self.assertAlmostEquals(self.acqs[0].bias, -64.76551)


class L1TDataTest(unittest.TestCase):

    def test_acquisition(self):
        self.assertTrue(os.path.exists(DATA_DIR))


if __name__ == '__main__':
    unittest.main()
