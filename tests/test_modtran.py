#!/usr/bin/env python

"""
Test several different configurations of the json input  file format.

"""
from os.path import join as pjoin, abspath, dirname
import unittest
import subprocess
import tempfile
import glob
from wagl.modtran import coefficients
from wagl.modtran import _get_solar_angles
import pandas as pd

DATA_DIR = pjoin(dirname(abspath(__file__)), 'data')
input_json = pjoin(DATA_DIR, 'TL_alb_0.json')
trans_file = pjoin(DATA_DIR, 'TL_alb.txt')
modtran_exe = 'mod6c_cons'


class ModtranTest(unittest.TestCase):

    """
    Test the modtran 6.0.1 program to see if runs completes
    successfully and outputs the expected '.chn' file and
    generates required atmospheric correction parameters

    """
    def test_modtran_run(self):

        """
        Test to run MODTRAN 6.0.1 and check if the transmittance
        result generated is as expected as in the 'trans_file'
        """
        with tempfile.TemporaryDirectory() as tmpdir:

            nbands = 8

            data = pd.read_csv(trans_file, nrows=nbands, header=None, delim_whitespace=True)

            columns = ['FS', 'FV', 'A', 'B', 'S', 'DIR', 'DIF', 'TS']
            df = pd.DataFrame(columns=columns, index=data.index)
            df['band_name'] = ['BAND-1', 'BAND-2', 'BAND-3', 'BAND-4', 'BAND-5', 'BAND-6', 'BAND-7', 'BAND-8']
            df['FS'] = data[1]
            df['FV'] = data[2]
            df['A'] = data[3]
            df['B'] = data[4]
            df['S'] = data[5]
            df['DIR'] = data[6]
            df['DIF'] = data[7]
            df['TS'] = data[8]
            df.set_index('band_name', inplace=True)

            subprocess.check_call([modtran_exe, input_json], cwd=tmpdir)

            chn_fname = glob.glob(pjoin(tmpdir, '*.chn'))[0]
            tp6_fname = glob.glob(pjoin(tmpdir, '*.tp6'))[0]

            solar_zenith = _get_solar_angles(tp6_fname)

            df_sz_angle = pd.DataFrame()
            df_sz_angle['solar_zenith'] = solar_zenith

            chn_data = pd.read_csv(chn_fname, skiprows=5, header=None, nrows=nbands,
                                   delim_whitespace=True)
            chn_data['band_name'] = chn_data[26]
            chn_data.drop(26, inplace=True, axis=1)

            chn_data.set_index('band_name', inplace=True)
            chn_data.columns = chn_data.columns.astype(str)

            nbar, _ = coefficients(chn_data, df_sz_angle, upward_radiation=None, downward_radiation=None)

            pd.testing.assert_frame_equal(df, nbar, check_less_precise=True)


if __name__ == '__main__':
    unittest.main()
