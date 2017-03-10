from __future__ import print_function
from __future__ import absolute_import

from .hdf5 import *
from .data import *
from .metadata import *
from .ancillary import *
from .acquisition import *
from .mtl import *
from .geobox import GriddedGeoBox
from .tle import *
from .brdf import *
from .interpolation import *
from .calculate_lon_lat_arrays import *
from .land_sea_masking import *
from .land_sea import get_land_sea_mask
from .modtran_profiles import *
from .modtran import *
from .margins import *
from .constants import PQAConstants
from .pqa_result import PQAResult
from .saturation_masking import set_saturation_bits
from .contiguity_masking import set_contiguity_bit
from .fc_utils import *
from .endmembers import *
from .thermal_conversion import get_landsat_temperature
from .acca_cloud_masking import calc_acca_cloud_mask
from .acca_cloud_masking import majority_filter
from .fmask_cloud_masking_wrapper import fmask_cloud_mask
from .cloud_shadow_masking import cloud_shadow
from .dsm import *

try:
    from __satellite_model import set_satmod # F2Py
    from __track_time_info import set_times # F2Py
    from __sat_sol_angles import angle # F2Py
    from __cast_shadow_mask import cast_shadow_main # F2Py
    from __exiting_angle import exiting_angle # F2Py
    from __incident_angle import incident_angle # F2Py
    from __slope_aspect import slope_aspect # F2Py
    from __surface_reflectance import reflectance # F2Py
    from __bilinear_interpolation import bilinear_interpolation #F2Py
    from .calculate_angles import *
    from .calculate_shadow_masks import *
    from .calculate_reflectance import calculate_reflectance
    from .calculate_incident_exiting_angles import *
    from .calculate_slope_aspect import *
except ImportError:
    msg = ('FORTRAN modules have not been built.\n'
           'Some functionality in library is disabled')
    print(msg)

from ._version import get_versions

__version__ = get_versions()['version']
del get_versions
