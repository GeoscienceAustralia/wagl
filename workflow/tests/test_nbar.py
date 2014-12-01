import luigi
import unittest
import workflow
import os

class GetElevationAncillaryDataTest(unittest.TestCase):

    def setUp(self):
        l1t_path = '../../gaip/tests/data/L1T/LS7_90-81_2009-04-15/UTM/LS7_ETM_OTH_P51_GALPGS01-002_090_081_20090415'
        self.task = workflow.GetElevationAncillaryDataTask(l1t_path)

    def test_invoke(self):
        luigi.build([self.task], local_scheduler=True)
        self.assertTrue(os.path.exists(self.task.output().path))

if __name__ == '__main__':
    unittest.main()
