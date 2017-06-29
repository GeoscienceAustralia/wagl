from __future__ import absolute_import
import unittest
import datetime
import rasterio
import os
from osgeo import osr
from gaip.acquisition import acquisitions
from gaip.acquisition import Landsat8Acquisition, LandsatAcquisition
from gaip.constants import BandType

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')

L5_MTL1 = os.path.join(DATA_DIR, 'L5090081_08120090407_MTL.txt')
L5_MTL2 = os.path.join(DATA_DIR, 'LT05_L1TP_095066_20100601_20170222_01_T1_MTL.txt')
L7_MTL1 = os.path.join(DATA_DIR, 'L71090081_08120090415_MTL.txt')
L7_MTL2 = os.path.join(DATA_DIR, 'LE07_L1TP_112066_20020218_20170221_01_T1_MTL.txt')
L8_MTL1 = os.path.join(DATA_DIR, 'LO80900842013284ASA00_MTL.txt')
L8_MTL2 = os.path.join(DATA_DIR, 'LC80990842016277LGN00_MTL.txt')


class AcquisitionLoadMtlTest(unittest.TestCase):

    def test_load_acquisitions_l5_mtl1(self):
        acq = acquisitions(L5_MTL1).get_acquisitions()
        self.assertEqual(len(acq), 7)

    def test_load_acquisitions_l5_mtl2(self):
        acq = acquisitions(L5_MTL2).get_acquisitions()
        self.assertEqual(len(acq), 8)

    def test_load_acquisitions_l7_mtl1(self):
        acq = acquisitions(L7_MTL1).get_acquisitions()
        self.assertEqual(len(acq), 9)

    def test_load_acquisitions_l7_mtl2(self):
        acq = acquisitions(L7_MTL2).get_acquisitions()
        self.assertEqual(len(acq), 10)

    def test_load_acquisitions_l8_mtl1(self):
        acq = acquisitions(L8_MTL1).get_acquisitions()
        self.assertEqual(len(acq), 10)

    def test_load_acquisitions_l8_mtl2(self):
        acq = acquisitions(L8_MTL2).get_acquisitions()
        self.assertEqual(len(acq), 12)


class AcquisitionsContainerTest(unittest.TestCase):

    def test_groups_l5_mtl1(self):
        scene = acquisitions(L5_MTL1)
        self.assertEqual(len(scene.groups), 1)

    def test_groups_l5_mtl2(self):
        scene = acquisitions(L5_MTL2)
        self.assertEqual(len(scene.groups), 1)

    def test_groups_l7_mtl1(self):
        scene = acquisitions(L7_MTL1)
        self.assertEqual(len(scene.groups), 1)

    def test_groups_l7_mtl2(self):
        scene = acquisitions(L7_MTL2)
        self.assertEqual(len(scene.groups), 1)

    def test_groups_l8_mtl1(self):
        scene = acquisitions(L8_MTL1)
        self.assertEqual(len(scene.groups), 1)

    def test_granules_ls5_mtl1(self):
        scene = acquisitions(L5_MTL1)
        self.assertIsNone(scene.granules[0])

    def test_granules_ls5_mtl2(self):
        scene = acquisitions(L5_MTL2)
        self.assertIsNone(scene.granules[0])

    def test_granules_ls7_mtl1(self):
        scene = acquisitions(L7_MTL1)
        self.assertIsNone(scene.granules[0])

    def test_granules_ls7_mtl2(self):
        scene = acquisitions(L7_MTL2)
        self.assertIsNone(scene.granules[0])

    def test_granules_ls8_mtl1(self):
        scene = acquisitions(L8_MTL1)
        self.assertIsNone(scene.granules[0])

    def test_granules_ls8_mtl2(self):
        scene = acquisitions(L8_MTL2)
        self.assertIsNone(scene.granules[0])


class Landsat5Mtl1AcquisitionTest(unittest.TestCase):

    def setUp(self):
        self.acqs = acquisitions(L5_MTL1).get_acquisitions()

    def test_type(self):
        for acq in self.acqs:
            self.assertTrue(isinstance(acq, LandsatAcquisition))

    def test_band_type(self):
        self.assertEqual(self.acqs[0].band_type, BandType.Reflective)
        self.assertEqual(self.acqs[1].band_type, BandType.Reflective)
        self.assertEqual(self.acqs[2].band_type, BandType.Reflective)
        self.assertEqual(self.acqs[3].band_type, BandType.Reflective)
        self.assertEqual(self.acqs[4].band_type, BandType.Reflective)
        self.assertEqual(self.acqs[5].band_type, BandType.Thermal)
        self.assertEqual(self.acqs[6].band_type, BandType.Reflective)

    # TODO: need extra name mappings to get the different names across versions
    # def test_grid_cell_size(self):
    #     self.assertEqual(self.acqs[0].grid_cell_size, 25.0)
    #     self.assertEqual(self.acqs[1].grid_cell_size, 25.0)
    #     self.assertEqual(self.acqs[2].grid_cell_size, 25.0)
    #     self.assertEqual(self.acqs[3].grid_cell_size, 25.0)
    #     self.assertEqual(self.acqs[4].grid_cell_size, 25.0)
    #     self.assertEqual(self.acqs[5].grid_cell_size, 25.0)
    #     self.assertEqual(self.acqs[6].grid_cell_size, 25.0)

    def test_scene_center_time(self):
        for acq in self.acqs:
            self.assertEqual(acq.scene_center_time,
                             datetime.time(23, 36, 9, 88050))

    def test_scene_center_datetime(self):
        for acq in self.acqs:
            self.assertEqual(acq.scene_center_datetime,
                             datetime.datetime(2009, 4, 7, 23, 36, 9, 88050))

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
        self.assertEqual(self.acqs[0].zone_number, -56)

    def test_sun_azimuth(self):
        self.assertEqual(self.acqs[0].sun_azimuth, 48.1772887)

    def test_sun_elevation(self):
        self.assertEqual(self.acqs[0].sun_elevation, 39.4014194)

    def test_gain(self):
        self.assertAlmostEqual(self.acqs[0].gain, 0.7658267716535433)

    def test_bias(self):
        self.assertAlmostEqual(self.acqs[0].bias, -2.2858267716535465)

    def test_sensor(self):
        for acq in self.acqs:
            self.assertEqual(acq.sensor_id, 'TM')

    def test_satellite(self):
        for acq in self.acqs:
            self.assertEqual(acq.spacecraft_id, 'LANDSAT_5')


class Landsat5Mtl2AcquisitionTest(unittest.TestCase):

    def setUp(self):
        self.acqs = acquisitions(L5_MTL2).get_acquisitions()

    def test_type(self):
        for acq in self.acqs:
            self.assertTrue(isinstance(acq, LandsatAcquisition))

    def test_band_type(self):
        self.assertEqual(self.acqs[0].band_type, BandType.Reflective)
        self.assertEqual(self.acqs[1].band_type, BandType.Reflective)
        self.assertEqual(self.acqs[2].band_type, BandType.Reflective)
        self.assertEqual(self.acqs[3].band_type, BandType.Reflective)
        self.assertEqual(self.acqs[4].band_type, BandType.Reflective)
        self.assertEqual(self.acqs[5].band_type, BandType.Thermal)
        self.assertEqual(self.acqs[6].band_type, BandType.Reflective)
        self.assertEqual(self.acqs[7].band_type, BandType.Quality)

    # TODO: need extra name mappings to get the different names across versions
    # def test_grid_cell_size(self):
    #     self.assertEqual(self.acqs[0].grid_cell_size, 25.0)
    #     self.assertEqual(self.acqs[1].grid_cell_size, 25.0)
    #     self.assertEqual(self.acqs[2].grid_cell_size, 25.0)
    #     self.assertEqual(self.acqs[3].grid_cell_size, 25.0)
    #     self.assertEqual(self.acqs[4].grid_cell_size, 25.0)
    #     self.assertEqual(self.acqs[5].grid_cell_size, 25.0)
    #     self.assertEqual(self.acqs[6].grid_cell_size, 25.0)

    def test_scene_center_time(self):
        for acq in self.acqs:
            self.assertEqual(acq.scene_center_time,
                             datetime.time(0, 4, 43, 174081))

    def test_scene_center_datetime(self):
        for acq in self.acqs:
            self.assertEqual(acq.scene_center_datetime,
                             datetime.datetime(2010, 6, 1, 0, 4, 43, 174081))

    def test_min_max_radiance_band1(self):
        self.assertEqual(self.acqs[0].min_radiance, -1.520)
        self.assertEqual(self.acqs[0].max_radiance, 193.0)

    def test_min_max_radiance_band2(self):
        self.assertEqual(self.acqs[0].min_radiance, -1.52)
        self.assertEqual(self.acqs[0].max_radiance, 193.0)

    def test_min_max_radiance_band3(self):
        self.assertEqual(self.acqs[0].min_radiance, -1.52)
        self.assertEqual(self.acqs[0].max_radiance, 193.0)

    def test_zone_number(self):
        self.assertEqual(self.acqs[0].utm_zone, -55)

    def test_sun_azimuth(self):
        self.assertEqual(self.acqs[0].sun_azimuth, 43.24285506)

    def test_sun_elevation(self):
        self.assertEqual(self.acqs[0].sun_elevation, 47.53234255)

    def test_gain(self):
        self.assertAlmostEqual(self.acqs[0].gain, 0.7658267716535433)

    def test_bias(self):
        self.assertAlmostEqual(self.acqs[0].bias, -2.2858267716535465)

    def test_sensor(self):
        for acq in self.acqs:
            self.assertEqual(acq.sensor_id, 'TM')

    def test_satellite(self):
        for acq in self.acqs:
            self.assertEqual(acq.spacecraft_id, 'LANDSAT_5')


class Landsat7Mtl1AcquisitionTest(unittest.TestCase):

    def setUp(self):
        self.acqs = acquisitions(L7_MTL1).get_acquisitions()

    def test_type(self):
        for acq in self.acqs:
            self.assertTrue(isinstance(acq, LandsatAcquisition))

    def test_band_type(self):
        self.assertEqual(self.acqs[0].band_type, BandType.Reflective)
        self.assertEqual(self.acqs[1].band_type, BandType.Reflective)
        self.assertEqual(self.acqs[2].band_type, BandType.Reflective)
        self.assertEqual(self.acqs[3].band_type, BandType.Reflective)
        self.assertEqual(self.acqs[4].band_type, BandType.Reflective)
        self.assertEqual(self.acqs[5].band_type, BandType.Thermal)
        self.assertEqual(self.acqs[6].band_type, BandType.Thermal)
        self.assertEqual(self.acqs[7].band_type, BandType.Reflective)
        self.assertEqual(self.acqs[8].band_type, BandType.Panchromatic)

    # TODO: need extra name mappings to get the different names across versions
    # def test_grid_cell_size(self):
    #     self.assertEqual(self.acqs[0].grid_cell_size, 25.0)
    #     self.assertEqual(self.acqs[1].grid_cell_size, 25.0)
    #     self.assertEqual(self.acqs[2].grid_cell_size, 25.0)
    #     self.assertEqual(self.acqs[3].grid_cell_size, 25.0)
    #     self.assertEqual(self.acqs[4].grid_cell_size, 25.0)
    #     self.assertEqual(self.acqs[5].grid_cell_size, 25.0)
    #     self.assertEqual(self.acqs[6].grid_cell_size, 25.0)
    #     self.assertEqual(self.acqs[7].grid_cell_size, 25.0)
    #     self.assertEqual(self.acqs[8].grid_cell_size, 12.5)

    def test_scene_center_time(self):
        for acq in self.acqs:
            self.assertEqual(acq.scene_center_time,
                             datetime.time(23, 39, 26, 931462))

    def test_scene_center_datetime(self):
        for acq in self.acqs:
            self.assertEqual(acq.scene_center_datetime,
                             datetime.datetime(2009, 4, 15, 23, 39, 26, 931462))

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
        self.assertEqual(self.acqs[0].zone_number, -56)

    def test_sun_azimuth(self):
        self.assertEqual(self.acqs[0].sun_azimuth, 44.5023798)

    def test_sun_elevation(self):
        self.assertEqual(self.acqs[0].sun_elevation, 37.9491813)

    def test_gain(self):
        self.assertAlmostEqual(self.acqs[0].gain, 0.7787401574803149)

    def test_bias(self):
        self.assertAlmostEqual(self.acqs[0].bias, -6.978740157480303)

    def test_sensor(self):
        for acq in self.acqs:
            self.assertEqual(acq.sensor_id, 'ETM+')

    def test_satellite(self):
        for acq in self.acqs:
            self.assertEqual(acq.spacecraft_id, 'LANDSAT_7')


class Landsat7Mtl2AcquisitionTest(unittest.TestCase):

    def setUp(self):
        self.acqs = acquisitions(L7_MTL2).get_acquisitions()

    def test_type(self):
        for acq in self.acqs:
            self.assertTrue(isinstance(acq, LandsatAcquisition))

    def test_band_type(self):
        self.assertEqual(self.acqs[0].band_type, BandType.Reflective)
        self.assertEqual(self.acqs[1].band_type, BandType.Reflective)
        self.assertEqual(self.acqs[2].band_type, BandType.Reflective)
        self.assertEqual(self.acqs[3].band_type, BandType.Reflective)
        self.assertEqual(self.acqs[4].band_type, BandType.Reflective)
        self.assertEqual(self.acqs[5].band_type, BandType.Thermal)
        self.assertEqual(self.acqs[6].band_type, BandType.Thermal)
        self.assertEqual(self.acqs[7].band_type, BandType.Reflective)
        self.assertEqual(self.acqs[8].band_type, BandType.Panchromatic)
        self.assertEqual(self.acqs[9].band_type, BandType.Quality)

    # TODO: need extra name mappings to get the different names across versions
    # def test_grid_cell_size(self):
    #     self.assertEqual(self.acqs[0].grid_cell_size, 25.0)
    #     self.assertEqual(self.acqs[1].grid_cell_size, 25.0)
    #     self.assertEqual(self.acqs[2].grid_cell_size, 25.0)
    #     self.assertEqual(self.acqs[3].grid_cell_size, 25.0)
    #     self.assertEqual(self.acqs[4].grid_cell_size, 25.0)
    #     self.assertEqual(self.acqs[5].grid_cell_size, 25.0)
    #     self.assertEqual(self.acqs[6].grid_cell_size, 25.0)
    #     self.assertEqual(self.acqs[7].grid_cell_size, 25.0)
    #     self.assertEqual(self.acqs[8].grid_cell_size, 12.5)

    def test_scene_center_time(self):
        for acq in self.acqs:
            self.assertEqual(acq.scene_center_time,
                             datetime.time(1, 47, 55, 878250))

    def test_scene_center_datetime(self):
        for acq in self.acqs:
            self.assertEqual(acq.scene_center_datetime,
                             datetime.datetime(2002, 2, 18, 1, 47, 55, 878250))

    def test_min_max_radiance_band1(self):
        self.assertEqual(self.acqs[0].min_radiance, -6.2)
        self.assertEqual(self.acqs[0].max_radiance, 191.6)

    def test_min_max_radiance_band2(self):
        self.assertEqual(self.acqs[0].min_radiance, -6.2)
        self.assertEqual(self.acqs[0].max_radiance, 191.6)

    def test_min_max_radiance_band3(self):
        self.assertEqual(self.acqs[0].min_radiance, -6.2)
        self.assertEqual(self.acqs[0].max_radiance, 191.6)

    def test_zone_number(self):
        self.assertEqual(self.acqs[0].utm_zone, -51)

    def test_sun_azimuth(self):
        self.assertEqual(self.acqs[0].sun_azimuth, 98.1470638)

    def test_sun_elevation(self):
        self.assertEqual(self.acqs[0].sun_elevation, 55.95447861)

    def test_gain(self):
        self.assertAlmostEqual(self.acqs[0].gain, 0.7787401574803149)

    def test_bias(self):
        self.assertAlmostEqual(self.acqs[0].bias, -6.978740157480303)

    def test_sensor(self):
        for acq in self.acqs:
            self.assertEqual(acq.sensor_id, 'ETM+')

    def test_satellite(self):
        for acq in self.acqs:
            self.assertEqual(acq.spacecraft_id, 'LANDSAT_7')


class Landsat8Mtl1AcquisitionTest(unittest.TestCase):

    def setUp(self):
        self.acqs = acquisitions(L8_MTL1).get_acquisitions()

    def test_type(self):
        for acq in self.acqs:
            self.assertTrue(isinstance(acq, Landsat8Acquisition))

    def test_scene_center_time(self):
        for acq in self.acqs:
            self.assertEqual(acq.scene_center_time,
                             datetime.time(23, 52, 10, 108347))

    def test_band_type(self):
        self.assertEqual(self.acqs[0].band_type, BandType.Reflective)
        self.assertEqual(self.acqs[1].band_type, BandType.Reflective)
        self.assertEqual(self.acqs[2].band_type, BandType.Reflective)
        self.assertEqual(self.acqs[3].band_type, BandType.Reflective)
        self.assertEqual(self.acqs[4].band_type, BandType.Reflective)
        self.assertEqual(self.acqs[5].band_type, BandType.Reflective)
        self.assertEqual(self.acqs[6].band_type, BandType.Reflective)
        self.assertEqual(self.acqs[7].band_type, BandType.Panchromatic)
        self.assertEqual(self.acqs[8].band_type, BandType.Atmosphere)
        self.assertEqual(self.acqs[9].band_type, BandType.Quality)

    # TODO: need extra name mappings to get the different names across versions
    # def test_grid_cell_size(self):
    #     self.assertEqual(self.acqs[0].grid_cell_size, 25.0)
    #     self.assertEqual(self.acqs[1].grid_cell_size, 25.0)
    #     self.assertEqual(self.acqs[2].grid_cell_size, 25.0)
    #     self.assertEqual(self.acqs[3].grid_cell_size, 25.0)
    #     self.assertEqual(self.acqs[4].grid_cell_size, 25.0)
    #     self.assertEqual(self.acqs[5].grid_cell_size, 25.0)
    #     self.assertEqual(self.acqs[6].grid_cell_size, 25.0)
    #     self.assertEqual(self.acqs[7].grid_cell_size, 12.5)
    #     self.assertEqual(self.acqs[8].grid_cell_size, 25.0)
    #     self.assertEqual(self.acqs[9].grid_cell_size, 25.0)

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
        self.assertEqual(self.acqs[0].utm_zone, -55)

    def test_sun_azimuth(self):
        self.assertEqual(self.acqs[0].sun_azimuth, 50.86088724)

    def test_sun_elevation(self):
        self.assertEqual(self.acqs[0].sun_elevation, 52.25003864)

    def test_gain(self):
        self.assertAlmostEqual(self.acqs[0].gain, 0.012953)

    def test_bias(self):
        self.assertAlmostEqual(self.acqs[0].bias, -64.76551)

    def test_sensor(self):
        for acq in self.acqs:
            self.assertEqual(acq.sensor_id, 'OLI')

    def test_satellite(self):
        for acq in self.acqs:
            self.assertEqual(acq.spacecraft_id, 'LANDSAT_8')


class Landsat8Mtl2AcquisitionTest(unittest.TestCase):

    def setUp(self):
        self.acqs = acquisitions(L8_MTL2).get_acquisitions()

    def test_type(self):
        for acq in self.acqs:
            self.assertTrue(isinstance(acq, Landsat8Acquisition))

    def test_scene_center_time(self):
        for acq in self.acqs:
            self.assertEqual(acq.scene_center_time,
                             datetime.time(0, 46, 10, 530409))

    def test_band_type(self):
        self.assertEqual(self.acqs[0].band_type, BandType.Reflective)
        self.assertEqual(self.acqs[1].band_type, BandType.Thermal)
        self.assertEqual(self.acqs[2].band_type, BandType.Thermal)
        self.assertEqual(self.acqs[3].band_type, BandType.Reflective)
        self.assertEqual(self.acqs[4].band_type, BandType.Reflective)
        self.assertEqual(self.acqs[5].band_type, BandType.Reflective)
        self.assertEqual(self.acqs[6].band_type, BandType.Reflective)
        self.assertEqual(self.acqs[7].band_type, BandType.Reflective)
        self.assertEqual(self.acqs[8].band_type, BandType.Reflective)
        self.assertEqual(self.acqs[9].band_type, BandType.Panchromatic)
        self.assertEqual(self.acqs[10].band_type, BandType.Atmosphere)
        self.assertEqual(self.acqs[11].band_type, BandType.Quality)

    # TODO: need extra name mappings to get the different names across versions
    # def test_grid_cell_size(self):
    #     self.assertEqual(self.acqs[0].grid_cell_size, 25.0)
    #     self.assertEqual(self.acqs[1].grid_cell_size, 25.0)
    #     self.assertEqual(self.acqs[2].grid_cell_size, 25.0)
    #     self.assertEqual(self.acqs[3].grid_cell_size, 25.0)
    #     self.assertEqual(self.acqs[4].grid_cell_size, 25.0)
    #     self.assertEqual(self.acqs[5].grid_cell_size, 25.0)
    #     self.assertEqual(self.acqs[6].grid_cell_size, 25.0)
    #     self.assertEqual(self.acqs[7].grid_cell_size, 12.5)
    #     self.assertEqual(self.acqs[8].grid_cell_size, 25.0)
    #     self.assertEqual(self.acqs[9].grid_cell_size, 25.0)

    def test_scene_center_datetime(self):
        for acq in self.acqs:
            self.assertEqual(acq.scene_center_datetime,
                             datetime.datetime(2016, 10, 3, 0, 46, 10, 530409))

    def test_min_max_radiance_band1(self):
        self.assertEqual(self.acqs[0].min_radiance, -62.69242)
        self.assertEqual(self.acqs[0].max_radiance, 759.16895)

    def test_min_max_radiance_band2(self):
        self.assertEqual(self.acqs[3].min_radiance, -64.19780)
        self.assertEqual(self.acqs[3].max_radiance, 777.39825)

    def test_min_max_radiance_band3(self):
        self.assertEqual(self.acqs[4].min_radiance, -59.15772)
        self.assertEqual(self.acqs[4].max_radiance, 716.36584)

    def test_lmin(self):
        self.assertEqual(self.acqs[0].lmin, -62.69242)

    def test_lmax(self):
        self.assertEqual(self.acqs[0].lmax, 759.16895)

    def test_qcalmin(self):
        self.assertEqual(self.acqs[0].qcalmin, 1)

    def test_qcalmax(self):
        self.assertEqual(self.acqs[0].qcalmax, 65535)

    def test_zone_number(self):
        self.assertEqual(self.acqs[0].utm_zone, -53)

    def test_sun_azimuth(self):
        self.assertEqual(self.acqs[0].sun_azimuth, 48.79660801)

    def test_sun_elevation(self):
        self.assertEqual(self.acqs[0].sun_elevation, 48.83189159)

    def test_gain(self):
        self.assertAlmostEqual(self.acqs[0].gain, 0.012541)

    def test_bias(self):
        self.assertAlmostEqual(self.acqs[0].bias, -62.70496)

    def test_sensor(self):
        for acq in self.acqs:
            self.assertEqual(acq.sensor_id, 'OLI_TIRS')

    def test_satellite(self):
        for acq in self.acqs:
            self.assertEqual(acq.spacecraft_id, 'LANDSAT_8')


if __name__ == '__main__':
    unittest.main()
