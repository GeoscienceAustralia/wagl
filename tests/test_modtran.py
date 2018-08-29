#!/usr/bin/env python

"""
Test several different configurations of the json input  file format.

"""

from os.path import join as pjoin, abspath, dirname
import unittest
import subprocess
import tempfile
import glob


DATA_DIR = pjoin(dirname(abspath(__file__)), 'data')
input_json = pjoin(DATA_DIR, 'POINT-0-ALBEDO-0.json')
chn_file = pjoin(DATA_DIR, 'POINT-0-ALBEDO-0.chn')
modtran_exe = '/g/data/v10/private/modules/MODTRAN/MODTRAN-6.0.1/MODTRAN6.0/bin/linux/mod6c_cons'
data_dir = pjoin('/g/data1a/u46/users/pd1813/test-run/', 'DATA')

class ModtranTest(unittest.TestCase):

    """
    Test the modtran 6.0.1 program to see if runs completes
    successfully and outputs the expected '.chn' file

    """
    def test_modtran_run(self):
        """
        Test to run modtran 6.0.1 program and check '*.chn' file
        output

        """
        with tempfile.TemporaryDirectory() as tmpdir:


            subprocess.check_call([modtran_exe, input_json, data_dir], cwd=tmpdir)

            chn_fname = glob.glob(pjoin(tmpdir, '*.chn'))[0]

            with open(chn_fname,'r') as f:
                data1 = ''.join(f.readlines())

            with open(chn_file, 'r') as src:
                data = ''.join(src.readlines())


            self.assertTrue(data1 == data)

if __name__ == '__main__':
    unittest.main()
