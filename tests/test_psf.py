import unittest

from pathlib import Path

import numpy as np
import numpy.testing as npt

from wagl.psf import read_tp7, compute_fwhm, compute_filter_matrix

DATA_DIR = Path(__file__).parent.joinpath("data", "WATER_ATCOR")


class TestPointSpreadFunction(unittest.TestCase):
    def setUp(self):
        tp7_file = DATA_DIR.joinpath("MM_alb_0.tp7")
        self.psf_data = read_tp7(tp7_file, 8)

    def test_readtp7(self):
        npt.assert_almost_equal(np.sum(self.psf_data['Band_8']), 13803.3752)
        npt.assert_almost_equal(np.sum(self.psf_data['Band_7']), 9282.04350)
        npt.assert_almost_equal(np.sum(self.psf_data['Band_6']), 7131.43700)
        npt.assert_almost_equal(np.sum(self.psf_data['Band_5']), 9144.7286)
        npt.assert_almost_equal(np.sum(self.psf_data['Band_4']), 13116.7809)
        npt.assert_almost_equal(np.sum(self.psf_data['Band_3']), 13957.8518)
        npt.assert_almost_equal(np.sum(self.psf_data['Band_2']), 15481.3534)
        npt.assert_almost_equal(np.sum(self.psf_data['Band_1']), 16452.4749)

    def test_computefwhm(self):
        npt.assert_almost_equal(compute_fwhm(self.psf_data['Band_8'], np.arange(251) *  10.0, 10.0), 467.444982995936)
        npt.assert_almost_equal(compute_fwhm(self.psf_data['Band_7'], np.arange(251) *  10.0, 10.0), 746.22377909361)
        npt.assert_almost_equal(compute_fwhm(self.psf_data['Band_6'], np.arange(251) *  10.0, 10.0), 742.7550338823073)
        npt.assert_almost_equal(compute_fwhm(self.psf_data['Band_5'], np.arange(251) *  10.0, 10.0), 544.4609683735308)
        npt.assert_almost_equal(compute_fwhm(self.psf_data['Band_4'], np.arange(251) *  10.0, 10.0), 465.13088437672735)
        npt.assert_almost_equal(compute_fwhm(self.psf_data['Band_3'], np.arange(251) *  10.0, 10.0), 468.07414661500195)
        npt.assert_almost_equal(compute_fwhm(self.psf_data['Band_2'], np.arange(251) *  10.0, 10.0), 469.7774722852594)
        npt.assert_almost_equal(compute_fwhm(self.psf_data['Band_1'], np.arange(251) *  10.0, 10.0), 473.56584102060157)


    def test_computefiltermatrix(self):
        kernel_matrix = compute_filter_matrix(self.psf_data['Band_1'], 30.0, 30.0, 40, 10.0)
        with open(DATA_DIR.joinpath("psf_b1.txt"), 'r') as fid:
            data = np.asarray([np.asarray([float(item) for item in val.rstrip().split()]) for val in fid.readlines()[9:]])
            self.assertEqual(kernel_matrix.shape, data.shape)
            npt.assert_almost_equal(np.sum(data), np.sum(kernel_matrix), decimal=2)


if __name__ == '__main__':
    unittest.main()
