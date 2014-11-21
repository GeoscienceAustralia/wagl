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

    def test_min_max_reflectance_band1(self):
        self.assertEqual(self.acqs[0].min_reflectance, 1.0)
        self.assertEqual(self.acqs[0].max_reflectance, 255.0)

    def test_min_max_reflectance_band2(self):
        self.assertEqual(self.acqs[0].min_reflectance, 1.0)
        self.assertEqual(self.acqs[0].max_reflectance, 255.0)

    def test_min_max_reflectance_band3(self):
        self.assertEqual(self.acqs[0].min_reflectance, 1.0)
        self.assertEqual(self.acqs[0].max_reflectance, 255.0)

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

    def test_min_max_reflectance_band1(self):
        self.assertEqual(self.acqs[0].min_reflectance, 1.0)
        self.assertEqual(self.acqs[0].max_reflectance, 255.0)

    def test_min_max_reflectance_band2(self):
        self.assertEqual(self.acqs[0].min_reflectance, 1.0)
        self.assertEqual(self.acqs[0].max_reflectance, 255.0)

    def test_min_max_reflectance_band3(self):
        self.assertEqual(self.acqs[0].min_reflectance, 1.0)
        self.assertEqual(self.acqs[0].max_reflectance, 255.0)

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

    def test_min_max_reflectance_band1(self):
        self.assertEqual(self.acqs[0].max_reflectance, 1.210700)
        self.assertEqual(self.acqs[0].min_reflectance, -0.099980)

    def test_min_max_reflectance_band2(self):
        self.assertEqual(self.acqs[0].max_reflectance, 1.210700)
        self.assertEqual(self.acqs[0].min_reflectance, -0.099980)

    def test_min_max_reflectance_band3(self):
        self.assertEqual(self.acqs[0].max_reflectance, 1.210700)
        self.assertEqual(self.acqs[0].min_reflectance, -0.099980)

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





    def test_gridded_geo_box(self):
        box = self.acqs[0].gridded_geo_box()
        self.assertEqual(type(box), gaip.GriddedGeoBox)
        self.assertEqual(box.origin, (644000.0, 6283000.0))
        self.assertEqual(box.corner, (879025.0, 6055975.0))
        self.assertEqual(box.shape, (9081, 9401))
        self.assertEqual(box.pixelsize, (25.0, 25.0))

    def test_gridded_geo_box_crs(self):
        box = self.acqs[0].gridded_geo_box()
        crs = box.crs
        self.assertEqual(type(crs), osr.SpatialReference)
        self.assertEqual(crs.IsProjected(), 1)
        self.assertEqual(crs.IsGeographic(), 0)
        self.assertEqual(crs.IsGeocentric(), 0)
        self.assertEqual(crs.GetUTMZone(), -55)
        self.assertEqual(crs.GetLinearUnits(), 1.0)
        self.assertEqual(crs.GetLinearUnitsName(), 'metre')
        self.assertAlmostEqual(crs.GetAngularUnits(), 0.0174532925199433)
        self.assertEqual(crs.ExportToProj4(), "+proj=utm +zone=55 +south +ellps=GRS80 " \
             "+towgs84=0,0,0,0,0,0,0 +units=m +no_defs " )

    def test_single_band_read(self):
        acq = self.acqs[0]
        band_data = acq.data()
        self.assertEqual(band_data.shape, (9081, 9401))

    def test_single_band_read_with_gridded_geo_box(self):
        acq = self.acqs[0]
        band_data, box = acq.data_and_box()
        self.assertEqual(band_data.shape, (9081, 9401))
        self.assertEqual(type(box), gaip.GriddedGeoBox)
        self.assertEqual(box.origin, (644000.0, 6283000.0))
        self.assertEqual(box.corner, (879025.0, 6055975.0))
        self.assertEqual(box.shape, (9081, 9401))
        self.assertEqual(box.pixelsize, (25.0, 25.0))

#    def test_multi_band_read(self):
#        bands = gaip.vstack_data(self.acqs)
#        self.assertEqual(bands.shape, (7, 9081, 9401))



class L1TDataTest(unittest.TestCase):

    def test_acquisition(self):
        self.assertTrue(os.path.exists(DATA_DIR))


if __name__ == '__main__':
    unittest.main()
