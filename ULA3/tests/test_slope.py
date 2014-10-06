#!/usr/bin/env python

import sys
import os
import argparse
import gc

import numpy
import numpy.testing as npt
from osgeo import gdal
from osgeo import osr

from EOtools.DatasetDrivers import SceneDataset

from ULA3.tc import run_slope
from ULA3._gdal_tools.py import Buffers

from general_tools import write_img, read_img, find_file

def calculate_slope_angles(scene_dataset, ref_dir, outdir, pixel_buffer=250):
    """
    
    """

    # Define the pixel buffering object
    pixel_buf = Buffers(pixel_buffer)

    # Check and load the required files from disk
    fname_sat_v = 'SAT_V.bin'
    fname_sat_az = 'SAT_AZ.bin'
    fname_sol_z = 'SOL_Z.bin'
    fname_sol_az = 'SOL_AZ.bin'
    fname_dsm = 'smoothed_dsm'

    view_angle = read_img(find_file(ref_dir, fname_sat_v))
    azi_angle = read_img(find_file(ref_dir, fname_sat_az))
    solar_angle = read_img(find_file(ref_dir, fname_sol_z))
    sazi_angle = read_img(find_file(ref_dir, fname_sol_az))
    dsm_data = read_img(find_file(ref_dir, fname_dsm))

    is_utm =  not scene_dataset.IsGeographic()

    slope_results = run_slope(scene_dataset, dsm_data, solar_angle, view_angle, sazi_angle, azi_angle, pixel_buf, is_utm)

    slope_results.dump_arrays(dump_path, scene_dataset, "ENVI", ".bin")


class test_slope_filenames():
    """
    
    """

    def __init__(self):
        """
        
        """

        self.reference_dir = reference_dir
        self.test_dir = test_dir

        self.fname_mask_self = 'mask_self'
        self.fname_slope = 'slope'
        self.fname_aspect = 'aspect'
        self.fname_incident = 'incident'
        self.fname_exiting = 'exiting'
        self.fname_azi_incident = 'azi_incident'
        self.fname_azi_exiting = 'azi_exiting'
        self.fname_rela_slope = 'rela_slope'

    def test_mask_self_ref(self):
        """
        Check that the mask self reference file exists.
        """
        fname = os.path.join(self.reference_dir, self.fname_mask_self)
        self.assertIs(os.path.exists(fname), True,
                      'Reference file does not exist: %s'%fname)

    def test_mask_self_tst(self):
        """
        Check that the mask self test file exists.
        """
        fname = os.path.join(self.test_dir, self.fname_mask_self)
        self.assertIs(os.path.exists(fname), True,
                      'Test file does not exist: %s'%fname)

    def test_slope_ref(self):
        """
        Check that the slope reference file exists.
        """
        fname = os.path.join(self.reference_dir, self.fname_slope)
        self.assertIs(os.path.exists(fname), True,
                      'Reference file does not exist: %s'%fname)

    def test_slope_tst(self):
        """
        Check that the slope test file exists.
        """
        fname = os.path.join(self.test_dir, self.fname_slope)
        self.assertIs(os.path.exists(fname), True,
                      'Test file does not exist: %s'%fname)

    def test_aspect_ref(self):
        """
        Check that the aspect reference file exists.
        """
        fname = os.path.join(self.reference_dir, self.fname_aspect)
        self.assertIs(os.path.exists(fname), True,
                      'Reference file does not exist: %s'%fname)

    def test_aspect_tst(self):
        """
        Check that the aspect test file exists.
        """
        fname = os.path.join(self.test_dir, self.fname_aspect)
        self.assertIs(os.path.exists(fname), True,
                      'Test file does not exist: %s'%fname)

    def test_incident_ref(self):
        """
        Check that the incident reference file exists.
        """
        fname = os.path.join(self.reference_dir, self.fname_incident)
        self.assertIs(os.path.exists(fname), True,
                      'Reference file does not exist: %s'%fname)

    def test_incident_tst(self):
        """
        Check that the mask self test file exists.
        """
        fname = os.path.join(self.test_dir, self.fname_incident)
        self.assertIs(os.path.exists(fname), True,
                      'Test file does not exist: %s'%fname)

    def test_exiting_ref(self):
        """
        Check that the exiting reference file exists.
        """
        fname = os.path.join(self.reference_dir, self.fname_exiting)
        self.assertIs(os.path.exists(fname), True,
                      'Reference file does not exist: %s'%fname)

    def test_exiting_tst(self):
        """
        Check that the exiting self test file exists.
        """
        fname = os.path.join(self.test_dir, self.fname_exiting)
        self.assertIs(os.path.exists(fname), True,
                      'Test file does not exist: %s'%fname)

    def test_azi_incident_ref(self):
        """
        Check that the azimuth incident reference file exists.
        """
        fname = os.path.join(self.reference_dir, self.fname_azi_incident)
        self.assertIs(os.path.exists(fname), True,
                      'Reference file does not exist: %s'%fname)

    def test_azi_incident_tst(self):
        """
        Check that the azimuth incident test file exists.
        """
        fname = os.path.join(self.test_dir, self.fname_azi_incident)
        self.assertIs(os.path.exists(fname), True,
                      'Test file does not exist: %s'%fname)

    def test_azi_exiting_ref(self):
        """
        Check that the azimuth exiting reference file exists.
        """
        fname = os.path.join(self.reference_dir, self.fname_azi_exiting)
        self.assertIs(os.path.exists(fname), True,
                      'Reference file does not exist: %s'%fname)

    def test_azi_exiting_tst(self):
        """
        Check that the azimuth exiting test file exists.
        """
        fname = os.path.join(self.test_dir, self.fname_azi_exiting)
        self.assertIs(os.path.exists(fname), True,
                      'Test file does not exist: %s'%fname)

    def test_rela_slope_ref(self):
        """
        Check that the relative slope reference file exists.
        """
        fname = os.path.join(self.reference_dir, self.fname_rela_slope)
        self.assertIs(os.path.exists(fname), True,
                      'Reference file does not exist: %s'%fname)

    def test_rela_slope_tst(self):
        """
        Check that the relative slope test file exists.
        """
        fname = os.path.join(self.test_dir, self.fname_rela_slope)
        self.assertIs(os.path.exists(fname), True,
                      'Test file does not exist: %s'%fname)


class test_slope_angles():

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(description='Calculates slope, aspect, incident, exiting & relative slope.')

    parser.add_argument('--L1T_dir', required=True, help='A directory path of a L1T scene.')
    parser.add_argument('--nbar_work_dir', required=True, help='A directory path to the associated NBAR working directory.')
    parser.add_argument('--outdir', required=True, help='A directory path that will contain the output files.')
    parser.add_argument('--dec_precision', default=4, help='The decimal precision used for array comparison')
    parser.add_argument('--int_precision', default=1, help='The integer precision used for array comparison')

    parsed_args = parser.parse_args()

