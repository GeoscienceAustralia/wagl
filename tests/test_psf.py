#! /usr/bin/env python

import unittest
from pathlib import Path
import tempfile
import json
import numpy as np
import numpy.testing as npt
from wagl.constants import Albedos
from wagl.acquisition import acquisitions
from wagl.psf import (
    compute_filter_matrix,
    compute_fwhm,
    read_tp7,
    _max_filter_size,
    prepare_modtran54,
    format_tp5,
)

DATA_DIR = Path(__file__).parent.joinpath("data", "WATER_ATCOR")
L8_SCENE = Path(__file__).parent.joinpath(
    "data", "LANDSAT8", "LS8_OLITIRS_OTH_P51_GALPGS01-032_090_084_20131011"
)


class TestPointSpreadFunction(unittest.TestCase):
    def setUp(self):
        tp7_file = DATA_DIR.joinpath("MM_alb_0.tp7")
        self.psf_data = read_tp7(tp7_file)

    def test_readtp7(self):
        npt.assert_almost_equal(np.sum(self.psf_data["BAND-8"]), 13803.3752)
        npt.assert_almost_equal(np.sum(self.psf_data["BAND-7"]), 9282.04350)
        npt.assert_almost_equal(np.sum(self.psf_data["BAND-6"]), 7131.43700)
        npt.assert_almost_equal(np.sum(self.psf_data["BAND-5"]), 9144.7286)
        npt.assert_almost_equal(np.sum(self.psf_data["BAND-4"]), 13116.7809)
        npt.assert_almost_equal(np.sum(self.psf_data["BAND-3"]), 13957.8518)
        npt.assert_almost_equal(np.sum(self.psf_data["BAND-2"]), 15481.3534)
        npt.assert_almost_equal(np.sum(self.psf_data["BAND-1"]), 16452.4749)

    def test_computefwhm(self):
        npt.assert_almost_equal(
            compute_fwhm(self.psf_data["BAND-8"], np.arange(251) * 10.0, 10.0),
            467.444982995936,
        )
        npt.assert_almost_equal(
            compute_fwhm(self.psf_data["BAND-7"], np.arange(251) * 10.0, 10.0),
            746.22377909361,
        )
        npt.assert_almost_equal(
            compute_fwhm(self.psf_data["BAND-6"], np.arange(251) * 10.0, 10.0),
            742.7550338823073,
        )
        npt.assert_almost_equal(
            compute_fwhm(self.psf_data["BAND-5"], np.arange(251) * 10.0, 10.0),
            544.4609683735308,
        )
        npt.assert_almost_equal(
            compute_fwhm(self.psf_data["BAND-4"], np.arange(251) * 10.0, 10.0),
            465.13088437672735,
        )
        npt.assert_almost_equal(
            compute_fwhm(self.psf_data["BAND-3"], np.arange(251) * 10.0, 10.0),
            468.07414661500195,
        )
        npt.assert_almost_equal(
            compute_fwhm(self.psf_data["BAND-2"], np.arange(251) * 10.0, 10.0),
            469.7774722852594,
        )
        npt.assert_almost_equal(
            compute_fwhm(self.psf_data["BAND-1"], np.arange(251) * 10.0, 10.0),
            473.56584102060157,
        )

    def test_maxfiltersize(self):
        self.assertEqual(_max_filter_size(30.0, 30.0, 40), 67)
        self.assertEqual(_max_filter_size(15.0, 15.0, 40), 133)

    def test_preparemodtran54(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_dir = Path(tmp_dir)
            modtran_exe = tmp_dir.joinpath(tmp_dir, "MODTRAN", "modtran_exe")
            acqs = acquisitions(L8_SCENE).get_acquisitions(group="RES-GROUP-1")

            with self.assertRaises(OSError):
                prepare_modtran54(acqs, 4, Albedos.ALBEDO_0, tmp_dir, modtran_exe)

            modtran_exe.parent.joinpath("DATA").mkdir(parents=True)
            prepare_modtran54(acqs[0], 4, Albedos.ALBEDO_0, tmp_dir, modtran_exe)
            modtran_work = tmp_dir.joinpath("POINT-4", "ALBEDO-0")

            self.assertTrue(modtran_work.joinpath("DATA").exists())
            self.assertTrue(
                modtran_work.joinpath(
                    Path(acqs[0].spectral_filter_filepath).name
                ).exists()
            )
            self.assertTrue(modtran_work.joinpath("mod5root.in").exists())

    def test_formattp5(self):
        with open(DATA_DIR.joinpath("MODTRAN-INPUT-DATA.json"), "r") as src:
            data = json.load(src)
            data["aerosol_type"] = 3
            data1_tp5, _ = format_tp5(data)

        with open(DATA_DIR.joinpath("POINT-4-ALBEDO-0.tp5"), "r") as src:
            data2_tp5 = "".join(src.readlines())

        self.assertEqual(data1_tp5, data2_tp5)

    def test_computefiltermatrix(self):
        kernel_matrix = compute_filter_matrix(
            self.psf_data["BAND-1"], 30.0, 30.0, 40, 10.0
        )
        with open(DATA_DIR.joinpath("psf_b1.txt"), "r") as fid:
            data = np.asarray(
                [
                    np.asarray([float(item) for item in val.rstrip().split()])
                    for val in fid.readlines()[9:]
                ]
            )
            npt.assert_array_almost_equal(kernel_matrix, kernel_matrix.T, decimal=10)
            npt.assert_array_almost_equal(data, kernel_matrix, decimal=2)
            npt.assert_array_almost_equal(kernel_matrix, data.T, decimal=2)


if __name__ == "__main__":
    unittest.main()
