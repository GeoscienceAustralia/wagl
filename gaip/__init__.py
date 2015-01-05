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
from modtran import *
from tc import *
from dsm import get_dsm

try:
    from calculate_angles import *
    from set_satmod import set_satmod # F2Py
    from set_times import set_times # F2Py
    from angle_all import angle # F2Py
    from filtering import filter_float as filter # F2Py
    from _shade_main_landsat_pixel import shade_main_landsat_pixel # F2Py
    from _slope_pixelsize_newpole import slope_pixelsize_newpole # F2Py
    from _brdf_terrain_newdiff_all import terrain_correction # F2Py
except ImportError:
    print 'You need to run `python setup.py build_ext -i` to compile modules'
    print 'Some functionality in library is disabled'
