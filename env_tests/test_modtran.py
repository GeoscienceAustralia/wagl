#!/usr/bin/env python

"""
Test several different configurations of the json input file format.

"""
import os
import shutil
import json
from os.path import join as pjoin, abspath, dirname
from posixpath import join as ppjoin
import unittest
import subprocess
import tempfile
import glob
from datetime import datetime
import pkg_resources
from functools import partial

import pandas as pd


from wagl.modtran import coefficients, _get_solar_angles, run_modtran, read_spectral_response

from wagl.hdf5 import read_h5_table

from wagl.constants import (
    Albedos, BandType, GroupName, DatasetName, Workflow, POINT_FMT,
    ALBEDO_FMT, POINT_ALBEDO_FMT
)

import mock


DATA_DIR = pjoin(dirname(abspath(__file__)), 'data')
INPUT_JSON = pjoin(DATA_DIR, 'TL_alb_0.json')
EXPECTED_CSV = pjoin(DATA_DIR, 'TL_alb_0.expected.csv')
MODTRAN_EXE = 'mod6c_cons'
SPECTRAL_RESPONSE_LS8 = pkg_resources.resource_filename(
    'wagl',
    'spectral_response/landsat8_vsir.flt'
)


def mock_spectral_response():
    """ Quick callable to return the spectral response from the repository """
    with pkg_resources.resource_stream(
            'wagl', 'spectral_response/landsat8_vsir.flt') as rsc_stream:
        return read_spectral_response(
            SPECTRAL_RESPONSE_LS8,
            False,
            range(2600, 349, -1)
        )


class ModtranTest(unittest.TestCase):

    def test_modtran_run(self):
        """
        Tests that the interface to modtran (run_modtran)
        works for known inputs.
        Used to validate environment configuration/setup
        """

        band_names = [
            'BAND-1', 'BAND-2', 'BAND-3', 'BAND-4',
            'BAND-5', 'BAND-6', 'BAND-7', 'BAND-8'
        ]
        point = 0
        albedo = Albedos.ALBEDO_0

        # setup mock acquistions object
        acquisitions = []
        for bandn in band_names:
            acq = mock.MagicMock()
            acq.acquisition_datetime = datetime(2001, 1, 1)
            acq.band_type = BandType.REFLECTIVE
            acq.spectral_response = mock_spectral_response
            acquisitions.append(acq)

        # setup mock atmospherics group
        attrs = { 'lonlat': 'TEST' }
        atmospherics = mock.MagicMock()
        atmospherics.attrs = attrs
        atmospherics_group = {
            POINT_FMT.format(p=point): atmospherics
        }

        # Compute base path -- prefix for hdf5 file
        base_path = ppjoin(GroupName.ATMOSPHERIC_RESULTS_GRP.value,
                           POINT_FMT.format(p=point))

        with tempfile.TemporaryDirectory() as workdir:
            run_dir = pjoin(workdir, POINT_FMT.format(p=point), ALBEDO_FMT.format(a=albedo.value))
            os.makedirs(run_dir)

            # TODO replace json_input copy with json input generation
            with open(INPUT_JSON, 'r') as fd:
                json_data = json.load(fd)
                for mod_input in json_data['MODTRAN']:
                    mod_input['MODTRANINPUT']['SPECTRAL']['FILTNM'] = SPECTRAL_RESPONSE_LS8

            with open(pjoin(run_dir, POINT_ALBEDO_FMT.format(p=point, a=albedo.value)) + ".json", 'w') as fd:
                json.dump(json_data, fd)

            fid = run_modtran(
                acquisitions,
                atmospherics_group,
                Workflow.STANDARD,
                npoints=12,  # number of track points
                point=point,
                albedos=[albedo],
                modtran_exe=MODTRAN_EXE,
                basedir=workdir,
                out_group=None
            )
            assert fid

            # Test base attrs
            assert fid[base_path].attrs['lonlat'] == 'TEST'
            assert fid[base_path].attrs['datetime'] == datetime(2001, 1, 1).isoformat()
            # test albedo headers?
            # Summarise modtran results to surface reflectance coefficients
            test_grp = fid[base_path][ALBEDO_FMT.format(a=albedo.value)]
            nbar_coefficients, _ = coefficients(
                read_h5_table(fid, pjoin(base_path, ALBEDO_FMT.format(a=albedo.value), DatasetName.CHANNEL.value)),
                read_h5_table(fid, pjoin(base_path, ALBEDO_FMT.format(a=albedo.value), DatasetName.SOLAR_ZENITH_CHANNEL.value))
            )

            expected = pd.read_csv(EXPECTED_CSV, index_col='band_name')
            pd.testing.assert_frame_equal(nbar_coefficients, expected, check_less_precise=True)


if __name__ == '__main__':
    unittest.main()
