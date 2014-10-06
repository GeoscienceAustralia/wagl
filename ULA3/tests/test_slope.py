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

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(description='Calculates slope, aspect, incident, exiting & relative slope.')

    parser.add_argument('--L1T_dir', required=True, help='A directory path of a L1T scene.')
    parser.add_argument('--nbar_work_dir', required=True, help='A directory path to the associated NBAR working directory.')
    parser.add_argument('--outdir', required=True, help='A directory path that will contain the output files.')
    parser.add_argument('--dec_precision', default=4, help='The decimal precision used for array comparison')
    parser.add_argument('--int_precision', default=1, help='The integer precision used for array comparison')

    parsed_args = parser.parse_args()

