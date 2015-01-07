import sys

from ancillary import *
from acquisition import *
from data import *
from mtl import *
from water import *
from GriddedGeoBox import GriddedGeoBox #FIXME: PEP8
from fast_cache import *
from tle import *
from brdf import *
from get_brdf import *
from calculate_lon_lat_arrays import *
from land_sea_masking import *
from land_sea import get_land_sea_mask 
from modtran import *
from margins import *
from constants import PQAConstants
from pqa_result import PQAResult


from saturation_masking import setSaturationBits
from contiguity_masking import setContiguityBit
from thermal_conversion import get_landsat_temperature
from acca_cloud_masking import calc_acca_cloud_mask
from fmask_cloud_masking_wrapper import FMaskCloudMask
from cloud_shadow_masking import Cloud_Shadow

try:
    from set_satmod import set_satmod # F2Py
    from set_times import set_times # F2Py
    from angle_all import angle # F2Py
    from _shade_main_landsat_pixel import shade_main_landsat_pixel # F2Py
    from _slope_pixelsize_newpole import slope_pixelsize_newpole # F2Py
    from _brdf_terrain_newdiff_all import terrain_correction # F2Py
    from calculate_angles import *
    from tc import *
    from dsm import get_dsm
except ImportError:
    msg = ('Run Makefile to build the Fortran modules.\n'
           'Some functionality in library is disabled')
    print sys.stderr, msg
