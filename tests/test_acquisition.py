from __future__ import absolute_import
import unittest
import datetime
from wagl.acquisition import acquisitions
from wagl.acquisition.landsat import Landsat8Acquisition, LandsatAcquisition
from wagl.constants import BandType
from wagl.temperature import temperature_at_sensor

from .data import LS5_SCENE1, LS7_SCENE1, LS8_SCENE1, LS8_SCENE1C2, LS8_SCENERTC2


class AcquisitionLoadMtlTest(unittest.TestCase):
    def test_load_acquisitions_ls5_scene1(self):
        acq_cont = acquisitions(LS5_SCENE1)
        self.assertEqual(len(acq_cont.get_acquisitions()), 7)
        self.assertEqual(len(acq_cont.get_acquisitions(only_supported_bands=False)), 7)

    def test_load_acquisitions_ls7_scene1(self):
        acq_cont = acquisitions(LS7_SCENE1)
        self.assertEqual(len(acq_cont.get_acquisitions()), 1)

    def test_highest_resolution_ls7_scene1(self):
        acq_cont = acquisitions(LS7_SCENE1)
        self.assertEqual(len(acq_cont.get_highest_resolution()[0]), 1)

    def test_res_group_1_ls7_scene1(self):
        acq_cont = acquisitions(LS7_SCENE1)
        self.assertEqual(len(acq_cont.get_acquisitions(group="RES-GROUP-1")), 8)

    def test_load_acquisitions_ls8_scene1(self):
        acq_cont = acquisitions(LS8_SCENE1)
        self.assertEqual(len(acq_cont.get_acquisitions()), 1)

    def test_highest_resolution_ls8_scene1(self):
        acq_cont = acquisitions(LS8_SCENE1)
        self.assertEqual(len(acq_cont.get_highest_resolution()[0]), 1)

    def test_res_group_1_ls8_scene1(self):
        acq_cont = acquisitions(LS8_SCENE1)
        self.assertEqual(len(acq_cont.get_acquisitions(group="RES-GROUP-1")), 9)

    def test_load_acquisitions_ls8c2_scene1(self):
        acq_cont = acquisitions(LS8_SCENE1C2)
        self.assertEqual(len(acq_cont.get_acquisitions()), 1)

    def test_highest_resolution_ls8c2_scene1(self):
        acq_cont = acquisitions(LS8_SCENE1C2)
        self.assertEqual(len(acq_cont.get_highest_resolution()[0]), 1)

    def test_res_group_1_ls8_scene1c2(self):
        acq_cont = acquisitions(LS8_SCENE1C2)
        # The data for this test is bad though
        # it is a copy of the C1 test data
        # This is used to keep the file size small
        self.assertEqual(len(acq_cont.get_acquisitions(group="RES-GROUP-1")), 9)

    def test_load_acquisitions_ls8c2_RT(self):
        acq_cont = acquisitions(LS8_SCENERTC2)
        self.assertEqual(len(acq_cont.get_acquisitions()), 1)

    def test_highest_resolution_ls8c2_RT(self):
        acq_cont = acquisitions(LS8_SCENERTC2)
        self.assertEqual(len(acq_cont.get_highest_resolution()[0]), 1)


class AcquisitionsContainerTest(unittest.TestCase):
    def test_groups_ls5_scene1(self):
        scene = acquisitions(LS5_SCENE1)
        self.assertEqual(len(scene.groups), 1)

    def test_groups_ls7_scene1(self):
        scene = acquisitions(LS7_SCENE1)
        self.assertEqual(len(scene.groups), 2)

    def test_groups_ls8_scene1(self):
        scene = acquisitions(LS8_SCENE1)
        self.assertEqual(len(scene.groups), 2)

    def test_groups_ls8_scene1rt(self):
        scene = acquisitions(LS8_SCENE1C2)
        self.assertEqual(len(scene.groups), 2)

    def test_granules_ls5_scene1(self):
        scene = acquisitions(LS5_SCENE1)
        self.assertEqual(scene.granules[0], "LT50900812009097ASA00")

    def test_granules_ls7_scene1(self):
        scene = acquisitions(LS7_SCENE1)
        self.assertEqual(scene.granules[0], "LE70900812009105ASA00")

    def test_granules_ls8_scene1(self):
        scene = acquisitions(LS8_SCENE1)
        self.assertEqual(scene.granules[0], "LC80900842013284LGN00")

    def test_granules_ls8_scene1C2(self):
        scene = acquisitions(LS8_SCENE1C2)
        self.assertEqual(scene.granules[0], "LC80920842020303LGN00")


class Landsat5Scene1AcquisitionTest(unittest.TestCase):
    def setUp(self):
        self.acqs = acquisitions(LS5_SCENE1).get_acquisitions()

    def test_type(self):
        for acq in self.acqs:
            self.assertTrue(isinstance(acq, LandsatAcquisition))

    def test_band_type(self):
        self.assertEqual(self.acqs[0].band_type, BandType.REFLECTIVE)
        self.assertEqual(self.acqs[1].band_type, BandType.REFLECTIVE)
        self.assertEqual(self.acqs[2].band_type, BandType.REFLECTIVE)
        self.assertEqual(self.acqs[3].band_type, BandType.REFLECTIVE)
        self.assertEqual(self.acqs[4].band_type, BandType.REFLECTIVE)
        self.assertEqual(self.acqs[5].band_type, BandType.THERMAL)
        self.assertEqual(self.acqs[6].band_type, BandType.REFLECTIVE)

    def test_acquisition_datetime(self):
        for acq in self.acqs:
            self.assertEqual(
                acq.acquisition_datetime, datetime.datetime(2009, 4, 7, 23, 36, 9, 88050)
            )

    def test_min_radiance_band1(self):
        self.assertEqual(self.acqs[0].min_radiance, -1.520)

    def test_max_radiance_band1(self):
        self.assertEqual(self.acqs[0].max_radiance, 193.0)

    def test_min_quantize_band1(self):
        self.assertEqual(self.acqs[0].min_quantize, 1.0)

    def test_max_quantize_band1(self):
        self.assertEqual(self.acqs[0].max_quantize, 255.0)

    def test_solar_azimuth(self):
        self.assertEqual(self.acqs[0].solar_azimuth, 48.17689881)

    def test_solar_elevation(self):
        self.assertEqual(self.acqs[0].solar_elevation, 39.40143058)

    def test_gain(self):
        self.assertAlmostEqual(self.acqs[0].gain, 0.7658267716535433)

    def test_bias(self):
        self.assertAlmostEqual(self.acqs[0].bias, -2.2858267716535465)

    def test_sensor_id(self):
        for acq in self.acqs:
            self.assertEqual(acq.sensor_id, "TM")

    def test_platform_id(self):
        for acq in self.acqs:
            self.assertEqual(acq.platform_id, "LANDSAT_5")

    def test_samples(self):
        self.assertEqual(self.acqs[0].samples, 74)

    def test_lines(self):
        self.assertEqual(self.acqs[0].lines, 65)

    def test_read(self):
        self.assertEqual(self.acqs[0].data()[40, 30], 53)

    def test_spectral_filter_cfg_vsir(self):
        self.assertEqual(self.acqs[0].spectral_filter_name, "landsat5_vsir.flt")

    def test_spectral_filter_cfg_thermal(self):
        self.assertEqual(self.acqs[5].spectral_filter_name, "landsat5_thermal.flt")

    def test_temperature(self):
        result = temperature_at_sensor(self.acqs[5], window=((40, 41), (40, 41)))
        self.assertAlmostEqual(result[0, 0], 293.76944042)

    def test_tzinfo(self):
        for acq in self.acqs:
            self.assertTrue(acq.acquisition_datetime, None)


class Landsat7Mtl1AcquisitionTest(unittest.TestCase):
    def setUp(self):
        self.acqs = acquisitions(LS7_SCENE1).get_acquisitions(group="RES-GROUP-1")

    def test_type(self):
        for acq in self.acqs:
            self.assertTrue(isinstance(acq, LandsatAcquisition))

    def test_band_type(self):
        self.assertEqual(self.acqs[0].band_type, BandType.REFLECTIVE)
        self.assertEqual(self.acqs[1].band_type, BandType.REFLECTIVE)
        self.assertEqual(self.acqs[2].band_type, BandType.REFLECTIVE)
        self.assertEqual(self.acqs[3].band_type, BandType.REFLECTIVE)
        self.assertEqual(self.acqs[4].band_type, BandType.REFLECTIVE)
        self.assertEqual(self.acqs[5].band_type, BandType.THERMAL)
        self.assertEqual(self.acqs[6].band_type, BandType.THERMAL)
        self.assertEqual(self.acqs[7].band_type, BandType.REFLECTIVE)

    def test_acquisition_datetime(self):
        for acq in self.acqs:
            self.assertEqual(
                acq.acquisition_datetime,
                datetime.datetime(2009, 4, 15, 23, 39, 26, 931462),
            )

    def test_min_radiance_band1(self):
        self.assertEqual(self.acqs[0].min_radiance, -6.2)

    def test_max_radiance_band1(self):
        self.assertEqual(self.acqs[0].max_radiance, 191.6)

    def test_min_quantize_band1(self):
        self.assertEqual(self.acqs[0].min_quantize, 1.0)

    def test_max_quantize_band1(self):
        self.assertEqual(self.acqs[0].max_quantize, 255.0)

    def test_solar_azimuth(self):
        self.assertEqual(self.acqs[0].solar_azimuth, 44.50200305)

    def test_solar_elevation(self):
        self.assertEqual(self.acqs[0].solar_elevation, 37.94917208)

    def test_gain(self):
        self.assertAlmostEqual(self.acqs[0].gain, 0.7787401574803149)

    def test_bias(self):
        self.assertAlmostEqual(self.acqs[0].bias, -6.978740157480303)

    def test_sensor_id(self):
        for acq in self.acqs:
            self.assertEqual(acq.sensor_id, "ETM+")

    def test_platform_id(self):
        for acq in self.acqs:
            self.assertEqual(acq.platform_id, "LANDSAT_7")

    def test_samples(self):
        self.assertEqual(self.acqs[0].samples, 75)

    def test_lines(self):
        self.assertEqual(self.acqs[0].lines, 65)

    def test_read(self):
        self.assertEqual(self.acqs[0].data()[40, 30], 50)

    def test_spectral_filter_cfg_vsir(self):
        self.assertEqual(self.acqs[0].spectral_filter_name, "landsat7_vsir.flt")

    def test_spectral_filter_cfg_thermal(self):
        self.assertEqual(self.acqs[5].spectral_filter_name, "landsat7_thermal.flt")

    def test_temperature61(self):
        result = temperature_at_sensor(self.acqs[5], window=((41, 42), (41, 42)))
        self.assertAlmostEqual(result[0, 0], 293.93158705)

    def test_temperature62(self):
        result = temperature_at_sensor(self.acqs[6], window=((41, 42), (41, 42)))
        self.assertAlmostEqual(result[0, 0], 293.990413299)

    def test_tzinfo(self):
        for acq in self.acqs:
            self.assertTrue(acq.acquisition_datetime, None)


class Landsat7PanAcquisitionTest(unittest.TestCase):
    def setUp(self):
        self.acqs = acquisitions(LS7_SCENE1).get_acquisitions(group="RES-GROUP-0")

    def test_type(self):
        for acq in self.acqs:
            self.assertTrue(isinstance(acq, LandsatAcquisition))

    def test_band_type(self):
        self.assertEqual(self.acqs[0].band_type, BandType.REFLECTIVE)

    def test_tzinfo(self):
        for acq in self.acqs:
            self.assertTrue(acq.acquisition_datetime, None)


class Landsat8Mtl1AcquisitionTest(unittest.TestCase):
    def setUp(self):
        self.acqs = acquisitions(LS8_SCENE1).get_acquisitions(group="RES-GROUP-1")

    def test_type(self):
        for acq in self.acqs:
            self.assertTrue(isinstance(acq, Landsat8Acquisition))

    def test_band_type(self):
        self.assertEqual(self.acqs[0].band_type, BandType.REFLECTIVE)
        self.assertEqual(self.acqs[1].band_type, BandType.THERMAL)
        self.assertEqual(self.acqs[2].band_type, BandType.THERMAL)
        self.assertEqual(self.acqs[3].band_type, BandType.REFLECTIVE)
        self.assertEqual(self.acqs[4].band_type, BandType.REFLECTIVE)
        self.assertEqual(self.acqs[5].band_type, BandType.REFLECTIVE)
        self.assertEqual(self.acqs[6].band_type, BandType.REFLECTIVE)
        self.assertEqual(self.acqs[6].band_type, BandType.REFLECTIVE)
        self.assertEqual(self.acqs[6].band_type, BandType.REFLECTIVE)

    def test_acquisition_datetime(self):
        for acq in self.acqs:
            self.assertEqual(
                acq.acquisition_datetime,
                datetime.datetime(2013, 10, 11, 23, 52, 10, 570334),
            )

    def test_min_radiance_band1(self):
        self.assertEqual(self.acqs[0].min_radiance, -63.00884)

    def test_max_radiance_band1(self):
        self.assertEqual(self.acqs[0].max_radiance, 763.00067)

    def test_min_quantize_band1(self):
        self.assertEqual(self.acqs[0].min_quantize, 1.0)

    def test_max_quantize_band1(self):
        self.assertEqual(self.acqs[0].max_quantize, 65535)

    def test_solar_azimuth(self):
        self.assertEqual(self.acqs[0].solar_azimuth, 50.86391564)

    def test_solar_elevation(self):
        self.assertEqual(self.acqs[0].solar_elevation, 52.04105874)

    def test_gain(self):
        self.assertAlmostEqual(self.acqs[0].gain, 0.0126, 4)

    def test_bias(self):
        self.assertAlmostEqual(self.acqs[0].bias, -63.0214, 4)

    def test_sensor_id(self):
        for acq in self.acqs:
            self.assertEqual(acq.sensor_id, "OLI")

    def test_platform_id(self):
        for acq in self.acqs:
            self.assertEqual(acq.platform_id, "LANDSAT_8")

    def test_samples(self):
        self.assertEqual(self.acqs[0].samples, 74)

    def test_lines(self):
        self.assertEqual(self.acqs[0].lines, 75)

    def test_read(self):
        self.assertEqual(self.acqs[0].data()[70, 30], 0)

    def test_spectral_filter_cfg_vsir(self):
        self.assertEqual(self.acqs[0].spectral_filter_name, "landsat8_vsir.flt")

    def test_spectral_filter_cfg_thermal(self):
        self.assertEqual(self.acqs[1].spectral_filter_name, "landsat8_thermal.flt")

    def test_temperature10(self):
        result = temperature_at_sensor(self.acqs[1], window=((41, 42), (41, 42)))
        self.assertAlmostEqual(result[0, 0], 299.91454310)

    def test_temperature11(self):
        result = temperature_at_sensor(self.acqs[2], window=((41, 42), (41, 42)))
        self.assertAlmostEqual(result[0, 0], 298.049253923)

    def test_tzinfo(self):
        for acq in self.acqs:
            self.assertTrue(acq.acquisition_datetime, None)


class Landsat8PanAcquisitionTest(unittest.TestCase):
    def setUp(self):
        self.acqs = acquisitions(LS8_SCENE1).get_acquisitions(group="RES-GROUP-0")

    def test_type(self):
        for acq in self.acqs:
            self.assertTrue(isinstance(acq, LandsatAcquisition))

    def test_band_type(self):
        self.assertEqual(self.acqs[0].band_type, BandType.REFLECTIVE)

    def test_tzinfo(self):
        for acq in self.acqs:
            self.assertTrue(acq.acquisition_datetime, None)


class Landsat8Mtl1c2AcquisitionTest(unittest.TestCase):
    def setUp(self):
        self.acqs = acquisitions(LS8_SCENE1C2).get_acquisitions(group="RES-GROUP-1")

    def test_type(self):
        for acq in self.acqs:
            self.assertTrue(isinstance(acq, Landsat8Acquisition))

    def test_band_type(self):
        self.assertEqual(self.acqs[0].band_type, BandType.REFLECTIVE)
        self.assertEqual(self.acqs[1].band_type, BandType.THERMAL)
        self.assertEqual(self.acqs[2].band_type, BandType.THERMAL)
        self.assertEqual(self.acqs[3].band_type, BandType.REFLECTIVE)
        self.assertEqual(self.acqs[4].band_type, BandType.REFLECTIVE)
        self.assertEqual(self.acqs[5].band_type, BandType.REFLECTIVE)
        self.assertEqual(self.acqs[6].band_type, BandType.REFLECTIVE)
        self.assertEqual(self.acqs[6].band_type, BandType.REFLECTIVE)
        self.assertEqual(self.acqs[6].band_type, BandType.REFLECTIVE)

    def test_acquisition_datetime(self):
        for acq in self.acqs:
            self.assertEqual(
                acq.acquisition_datetime,
                datetime.datetime(2020, 10, 29, 0, 2, 59, 26835),
            )

    def test_min_radiance_band1(self):
        self.assertEqual(self.acqs[0].min_radiance, -63.61862)

    def test_max_radiance_band1(self):
        self.assertEqual(self.acqs[0].max_radiance, 770.38470)

    def test_min_quantize_band1(self):
        self.assertEqual(self.acqs[0].min_quantize, 1.0)

    def test_max_quantize_band1(self):
        self.assertEqual(self.acqs[0].max_quantize, 65535)

    def test_solar_azimuth(self):
        self.assertEqual(self.acqs[0].solar_azimuth, 57.65543514)

    def test_solar_elevation(self):
        self.assertEqual(self.acqs[0].solar_elevation, 56.77807119)

    def test_gain(self):
        self.assertAlmostEqual(self.acqs[0].gain, 0.012726, 4)

    def test_bias(self):
        self.assertAlmostEqual(self.acqs[0].bias, -63.63135, 4)

    def test_sensor_id(self):
        for acq in self.acqs:
            self.assertEqual(acq.sensor_id, "OLI")

    def test_platform_id(self):
        for acq in self.acqs:
            self.assertEqual(acq.platform_id, "LANDSAT_8")

    def test_samples(self):
        self.assertEqual(self.acqs[0].samples, 74)

    def test_lines(self):
        self.assertEqual(self.acqs[0].lines, 75)

    def test_read(self):
        self.assertEqual(self.acqs[0].data()[70, 30], 0)

    def test_spectral_filter_cfg_vsir(self):
        self.assertEqual(self.acqs[0].spectral_filter_name, "landsat8_vsir.flt")

    def test_spectral_filter_cfg_thermal(self):
        self.assertEqual(self.acqs[1].spectral_filter_name, "landsat8_thermal.flt")

    def test_temperature10(self):
        result = temperature_at_sensor(self.acqs[1], window=((41, 42), (41, 42)))
        self.assertAlmostEqual(result[0, 0], 299.91454310)

    def test_temperature11(self):
        result = temperature_at_sensor(self.acqs[2], window=((41, 42), (41, 42)))
        self.assertAlmostEqual(result[0, 0], 298.049253923)

    def test_tzinfo(self):
        for acq in self.acqs:
            self.assertTrue(acq.acquisition_datetime, None)


class Landsat8C2PanAcquisitionTest(unittest.TestCase):
    def setUp(self):
        self.acqs = acquisitions(LS8_SCENE1C2).get_acquisitions(group="RES-GROUP-0")

    def test_type(self):
        for acq in self.acqs:
            self.assertTrue(isinstance(acq, LandsatAcquisition))

    def test_band_type(self):
        self.assertEqual(self.acqs[0].band_type, BandType.REFLECTIVE)

    def test_tzinfo(self):
        for acq in self.acqs:
            self.assertTrue(acq.acquisition_datetime, None)


if __name__ == "__main__":
    unittest.main()
