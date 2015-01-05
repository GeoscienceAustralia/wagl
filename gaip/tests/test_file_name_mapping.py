import unittest
from gaip.pq import nbar_name_from_l1t, pqa_name_from_l1t

class TestDatasetNameMappings(unittest.TestCase):

    def test_l1t_to_nbar(self):
        l1t_name = 'LS5_TM_OTH_P51_GALPGS01-002_100_081_20100228'
        nbar_name = 'LS5_TM_NBAR_P54_GANBAR01-002_100_081_20100228'
        self.assertEqual(nbar_name, nbar_name_from_l1t(l1t_name))

    def test_LS8_l1t_to_nbar(self):
        l1t_name = 'LS8_OLITIRS_OTH_P51_GALPGS01-015_101_078_20141012'
        nbar_name = 'LS8_OLI_TIRS_NBAR_P54_GANBAR01-015_101_078_20141012'
        self.assertEqual(nbar_name, nbar_name_from_l1t(l1t_name))

    def test_l1t_to_pqa(self):
        l1t_name = 'LS5_TM_OTH_P51_GALPGS01-002_100_081_20100228'
        pqa_name = 'LS5_TM_PQ_P55_GAPQ01-002_100_081_20100228'
        self.assertEqual(pqa_name, pqa_name_from_l1t(l1t_name))

    def test_LS8_l1t_to_pqa(self):
        l1t_name = 'LS8_OLITIRS_OTH_P51_GALPGS01-015_101_078_20141012'
        pqa_name = 'LS8_OLI_TIRS_PQ_P55_GAPQ01-015_101_078_20141012'
        self.assertEqual(pqa_name, pqa_name_from_l1t(l1t_name))



if __name__ == '__main__':
    unittest.main()
