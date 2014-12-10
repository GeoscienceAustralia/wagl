import luigi
import unittest
import workflow
import logging
import gaip
import os

from os.path import join as pjoin, dirname, exists


L1T_PATH = pjoin(dirname(gaip.__file__), 'tests', 'data', 'L1T')
LS7_PATH = pjoin(L1T_PATH, 'LS7_90-81_2009-04-15', 'UTM',
                 'LS7_ETM_OTH_P51_GALPGS01-002_090_081_20090415')


logging.disable(logging.INFO)


class ElevationAncillaryDataTest(unittest.TestCase):

    def setUp(self):
        self.task = workflow.GetElevationAncillaryDataTask(LS7_PATH)

    def test_invoke(self):
        luigi.build([self.task], local_scheduler=True)
        self.assertTrue(self.task.output().exists())

    def test_value(self):
        luigi.build([self.task], local_scheduler=True)
        self.assertTrue(self.task.complete())
        value = workflow.load(self.task.output())
        self.assertAlmostEqual(value['value'], 0.292)


class OzoneAncillaryDataTest(unittest.TestCase):

    def setUp(self):
        self.task = workflow.GetOzoneAncillaryDataTask(LS7_PATH)

    def test_invoke(self):
        luigi.build([self.task], local_scheduler=True)
        self.assertTrue(self.task.output().exists())

    def test_value(self):
        luigi.build([self.task], local_scheduler=True)
        self.assertTrue(self.task.complete())
        value = workflow.load(self.task.output())
        self.assertAlmostEqual(value['value'], 0.26800001)


class SolarIrradianceAncillaryDataTest(unittest.TestCase):

    def setUp(self):
        self.task = workflow.GetSolarIrradianceAncillaryDataTask(LS7_PATH)

    def test_invoke(self):
        luigi.build([self.task], local_scheduler=True)
        self.assertTrue(self.task.output().exists())

    def test_value(self):
        luigi.build([self.task], local_scheduler=True)
        self.assertTrue(self.task.complete())
        value = workflow.load(self.task.output())
        self.assertAlmostEqual(value[1], 1997.0)


class SolarDistanceAncillaryDataTest(unittest.TestCase):

    def setUp(self):
        self.task = workflow.GetSolarDistanceAncillaryDataTask(LS7_PATH)

    def test_invoke(self):
        luigi.build([self.task], local_scheduler=True)
        self.assertTrue(self.task.output().exists())

    def test_value(self):
        luigi.build([self.task], local_scheduler=True)
        self.assertTrue(self.task.complete())
        value = workflow.load(self.task.output())
        self.assertAlmostEqual(value, 1.00325)


class WaterVapourAncillaryDataTest(unittest.TestCase):

    def setUp(self):
        self.task = workflow.GetWaterVapourAncillaryDataTask(LS7_PATH)

    def test_invoke(self):
        luigi.build([self.task], local_scheduler=True)
        self.assertTrue(self.task.output().exists())

    def test_value(self):
        luigi.build([self.task], local_scheduler=True)
        self.assertTrue(self.task.complete())
        value = workflow.load(self.task.output())
        self.assertAlmostEqual(value['value'], 1.5209991455078127)


class AerosolAncillaryDataTest(unittest.TestCase):

    def setUp(self):
        self.task = workflow.GetAerosolAncillaryDataTask(LS7_PATH)

    def test_invoke(self):
        luigi.build([self.task], local_scheduler=True)
        self.assertTrue(self.task.output().exists())

    def test_value(self):
        luigi.build([self.task], local_scheduler=True)
        self.assertTrue(self.task.complete())
        value = workflow.load(self.task.output())
        self.assertAlmostEqual(value['value'], 0.041519)


class BrdfAncillaryDataTest(unittest.TestCase):

    def setUp(self):
        self.task = workflow.GetBrdfAncillaryDataTask(LS7_PATH)

    def test_invoke(self):
        luigi.build([self.task], local_scheduler=True)
        self.assertTrue(self.task.output().exists())

    def test_value(self):
        luigi.build([self.task], local_scheduler=True)
        self.assertTrue(self.task.complete())
        value = workflow.load(self.task.output())

        self.assertAlmostEqual(value[(1, 'iso')]['value'], 0.0450349829466)
        self.assertAlmostEqual(value[(3, 'iso')]['value'], 0.0898902284207)
        self.assertAlmostEqual(value[(4, 'iso')]['value'], 0.2509117168502)
        self.assertAlmostEqual(value[(5, 'iso')]['value'], 0.2685744810199)
        #self.assertAlmostEqual(value[(6,'iso')]['value'], 0.2685744810199)
        self.assertAlmostEqual(value[(7, 'iso')]['value'], 0.1589079075588)

        self.assertAlmostEqual(value[(1, 'geo')]['value'], 0.0097146516233)
        self.assertAlmostEqual(value[(2, 'geo')]['value'], 0.0208474632336)
        self.assertAlmostEqual(value[(3, 'geo')]['value'], 0.0208628134638)
        self.assertAlmostEqual(value[(4, 'geo')]['value'], 0.0227648405092)
        self.assertAlmostEqual(value[(5, 'geo')]['value'], 0.0385139467520)
        self.assertAlmostEqual(value[(7, 'geo')]['value'], 0.0284719055920)

        self.assertAlmostEqual(value[(1, 'vol')]['value'], 0.0169988198954)
        self.assertAlmostEqual(value[(2, 'vol')]['value'], 0.0342593849179)
        self.assertAlmostEqual(value[(3, 'vol')]['value'], 0.0321581645880)
        self.assertAlmostEqual(value[(4, 'vol')]['value'], 0.1914701312401)
        self.assertAlmostEqual(value[(5, 'vol')]['value'], 0.1136554213528)
        self.assertAlmostEqual(value[(7, 'vol')]['value'], 0.0491882392382)


class AncillaryDataTest(unittest.TestCase):

    def setUp(self):
        self.task = workflow.GetAncillaryData(LS7_PATH)

    def test_invoke(self):
        luigi.build([self.task], local_scheduler=True)
        self.assertTrue(self.task.complete())


class CalculateLatGridTest(unittest.TestCase):

    def setUp(self):
        self.task = workflow.CalculateLatGrid(LS7_PATH)

    def test_invoke(self):
        luigi.build([self.task], local_scheduler=True)
        self.assertTrue(self.task.complete())


class CalculateLonGridTest(unittest.TestCase):

    def setUp(self):
        self.task = workflow.CalculateLonGrid(LS7_PATH)

    def test_invoke(self):
        luigi.build([self.task], local_scheduler=True)
        self.assertTrue(self.task.complete())


class CalculateLonLatGridsTest(unittest.TestCase):

    def setUp(self):
        self.task = workflow.CalculateLatLonGrids(LS7_PATH)

    def test_invoke(self):
        luigi.build([self.task], local_scheduler=True)
        self.assertTrue(self.task.complete())


if __name__ == '__main__':
    unittest.main()
