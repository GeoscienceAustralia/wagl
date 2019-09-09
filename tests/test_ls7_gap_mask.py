#!/usr/bin/env python


import unittest
from wagl.acquisition import acquisitions


from .data import LS7_GAP_MASK


class GapMaskRadianceTest(unittest.TestCase):
    """
    Test that the SLC gap mask loads correctly.
    The values to check against were derived by manually selecting
    the band and the corresponding gap mask, and creating a null mask.
    """

    def setUp(self):
        self.acqs = acquisitions(LS7_GAP_MASK).get_all_acquisitions()

    def test_band8(self):
        acq = self.acqs[0]
        mask = acq.radiance_data() == -999
        count = mask.sum()
        self.assertEqual(count, 259512)

    def test_band1(self):
        acq = self.acqs[1]
        mask = acq.radiance_data() == -999
        count = mask.sum()
        self.assertEqual(count, 64746)

    def test_band2(self):
        acq = self.acqs[2]
        mask = acq.radiance_data() == -999
        count = mask.sum()
        self.assertEqual(count, 64766)

    def test_band3(self):
        acq = self.acqs[3]
        mask = acq.radiance_data() == -999
        count = mask.sum()
        self.assertEqual(count, 64761)

    def test_band4(self):
        acq = self.acqs[4]
        mask = acq.radiance_data() == -999
        count = mask.sum()
        self.assertEqual(count, 64770)

    def test_band5(self):
        acq = self.acqs[5]
        mask = acq.radiance_data() == -999
        count = mask.sum()
        self.assertEqual(count, 64769)

    def test_band61(self):
        acq = self.acqs[6]
        mask = acq.radiance_data() == -999
        count = mask.sum()
        self.assertEqual(count, 64862)

    def test_band62(self):
        acq = self.acqs[7]
        mask = acq.radiance_data() == -999
        count = mask.sum()
        self.assertEqual(count, 64898)

    def test_band7(self):
        acq = self.acqs[8]
        mask = acq.radiance_data() == -999
        count = mask.sum()
        self.assertEqual(count, 64747)

