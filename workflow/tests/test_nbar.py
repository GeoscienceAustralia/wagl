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
        self.assertAlmostEqual(value['value'], 1.5209991455078127)

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
        self.assertAlmostEqual(value['value'], 1.5209991455078127)

if __name__ == '__main__':
    unittest.main()
