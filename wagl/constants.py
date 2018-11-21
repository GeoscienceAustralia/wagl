"""
Constants
---------
"""
# pylint: disable=attribute-defined-outside-init

from __future__ import absolute_import, print_function
from os.path import join as pjoin
import re
from enum import Enum


POINT_FMT = 'POINT-{p}'
ALBEDO_FMT = 'ALBEDO-{a}'
POINT_ALBEDO_FMT = ''.join([POINT_FMT, '-', ALBEDO_FMT])


class Workflow(Enum):
    """
    Represents the different workflow that wagl can run.

    *standard* Indicates both NBAR and SBT workflows will run
    *nbar* Indicates NBAR only
    *sbt* Indicates SBT only
    """

    STANDARD = 1
    NBAR = 2
    SBT = 3

    @property
    def atmos_coefficients(self):
        """
        Returns the atmospheric coefficients names used for interpolation
        for a given Workflow.<option>.
        """
        atmos_var = list(AtmosphericCoefficients)
        fmap = {Workflow.STANDARD: atmos_var,
                Workflow.NBAR: atmos_var[0:8],
                Workflow.SBT: atmos_var[8:]}
        return fmap.get(self)

    @property
    def albedos(self):
        """
        Returns the albedo names used for specific Atmospheric
        evaluations for a given Workflow.<option>.
        """
        albs = list(Albedos)
        amap = {Workflow.STANDARD: albs,
                Workflow.NBAR: albs[0:-1],
                Workflow.SBT: [albs[-1]]}
        return amap.get(self)

    @property
    def ard_products(self):
        """
        Returns the ARD products available for a given
        Workflow.<option>.
        """
        products = list(ArdProducts)
        amap = {Workflow.STANDARD: products,
                Workflow.NBAR: products[0:-1],
                Workflow.SBT: [products[-1]]}
        return amap.get(self)


class BandType(Enum):
    """
    Represents the Band Type a given acquisition falls under.
    """

    REFLECTIVE = 0
    THERMAL = 1
    PANCHROMATIC = 2
    ATMOSPHERE = 3
    QUALITY = 4


class DatasetName(Enum):
    """
    Defines the dataset names or format descriptors, that are used
    for creating and accessing throughout the code base.
    """

    # wagl.ancillary
    COORDINATOR = 'COORDINATOR'
    DEWPOINT_TEMPERATURE = 'DEWPOINT-TEMPERATURE'
    TEMPERATURE_2M = 'TEMPERATURE-2METRE'
    SURFACE_PRESSURE = 'SURFACE-PRESSURE'
    SURFACE_GEOPOTENTIAL = 'SURFACE-GEOPOTENTIAL-HEIGHT'
    SURFACE_RELATIVE_HUMIDITY = 'SURFACE-RELATIVE-HUMIDITY'
    GEOPOTENTIAL = 'GEO-POTENTIAL'
    TEMPERATURE = 'TEMPERATURE'
    RELATIVE_HUMIDITY = 'RELATIVE-HUMIDITY'
    ATMOSPHERIC_PROFILE = 'ATMOSPHERIC-PROFILE'
    AEROSOL = 'AEROSOL'
    WATER_VAPOUR = 'WATER-VAPOUR'
    OZONE = 'OZONE'
    ELEVATION = 'ELEVATION'
    BRDF_FMT = "BRDF-{parameter}-{band_name}"
    ECMWF_PATH_FMT = pjoin('{product}', '{year}', 'tif', '{product}_*.tif')

    # wagl.longitude_latitude_arrays
    LON = 'LONGITUDE'
    LAT = 'LATITUDE'

    # wagl.satellite_solar_angles
    SATELLITE_VIEW = 'SATELLITE-VIEW'
    SATELLITE_AZIMUTH = 'SATELLITE-AZIMUTH'
    SOLAR_ZENITH = 'SOLAR-ZENITH'
    SOLAR_ZENITH_CHANNEL = 'SOLAR-ZENITH-CHANNEL'
    SOLAR_AZIMUTH = 'SOLAR-AZIMUTH'
    RELATIVE_AZIMUTH = 'RELATIVE-AZIMUTH'
    TIME = 'TIMEDELTA'
    CENTRELINE = 'CENTRELINE'
    BOXLINE = 'BOXLINE'
    SPHEROID = 'SPHEROID'
    ORBITAL_ELEMENTS = 'ORBITAL-ELEMENTS'
    SATELLITE_MODEL = 'SATELLITE-MODEL'
    SATELLITE_TRACK = 'SATELLITE-TRACK'
    GENERIC = 'GENERIC'

    # wagl.incident_exiting_angles
    INCIDENT = 'INCIDENT'
    AZIMUTHAL_INCIDENT = 'AZIMUTHAL-INCIDENT'
    EXITING = 'EXITING'
    AZIMUTHAL_EXITING = 'AZIMUTHAL-EXITING'
    RELATIVE_SLOPE = 'RELATIVE-SLOPE'

    # wagl.reflectance
    REFLECTANCE_FMT = 'REFLECTANCE/{product}/{band_name}'

    # wagl.temperature
    TEMPERATURE_FMT = 'THERMAL/{product}/{band_name}'

    # wagl.terrain_shadow_masks
    SELF_SHADOW = 'SELF-SHADOW'
    CAST_SHADOW_FMT = 'CAST-SHADOW-{source}'
    COMBINED_SHADOW = 'COMBINED-TERRAIN-SHADOW'

    # wagl.slope_aspect
    SLOPE = 'SLOPE'
    ASPECT = 'ASPECT'

    # wagl.dsm
    DSM = 'DSM'
    DSM_SMOOTHED = 'DSM-SMOOTHED'

    # wagl.interpolation
    INTERPOLATION_FMT = '{coefficient}/{band_name}'

    # wagl.modtran
    MODTRAN_INPUT = 'MODTRAN-INPUT-DATA'
    FLUX = 'FLUX'
    ALTITUDES = 'ALTITUDES'
    SOLAR_IRRADIANCE = 'SOLAR-IRRADIANCE'
    UPWARD_RADIATION_CHANNEL = 'UPWARD-RADIATION-CHANNEL'
    DOWNWARD_RADIATION_CHANNEL = 'DOWNWARD-RADIATION-CHANNEL'
    CHANNEL = 'CHANNEL'
    NBAR_COEFFICIENTS = 'NBAR-COEFFICIENTS'
    SBT_COEFFICIENTS = 'SBT-COEFFICIENTS'

    # wagl.pq
    PQ_FMT = 'PIXEL-QUALITY/{produt}/PIXEL-QUALITY'

    # metadata
    NBAR_YAML = 'METADATA/NBAR-METADATA'
    PQ_YAML = 'METADATA/PQ-METADATA'
    SBT_YAML = 'METADATA/SBT-METADATA'


class GroupName(Enum):
    """
    Defines the group names or format descriptors, that are used
    for creating and accessing throughout the code base.
    """

    LON_LAT_GROUP = 'LONGITUDE-LATITUDE'
    SAT_SOL_GROUP = 'SATELLITE-SOLAR'
    ANCILLARY_GROUP = 'ANCILLARY'
    ANCILLARY_AVG_GROUP = 'AVERAGED-ANCILLARY'
    ATMOSPHERIC_INPUTS_GRP = 'ATMOSPHERIC-INPUTS'
    ATMOSPHERIC_RESULTS_GRP = 'ATMOSPHERIC-RESULTS'
    COEFFICIENTS_GROUP = 'ATMOSPHERIC-COEFFICIENTS'
    INTERP_GROUP = 'INTERPOLATED-ATMOSPHERIC-COEFFICIENTS'
    ELEVATION_GROUP = 'ELEVATION'
    SLP_ASP_GROUP = 'SLOPE-ASPECT'
    INCIDENT_GROUP = 'INCIDENT-ANGLES'
    EXITING_GROUP = 'EXITING-ANGLES'
    REL_SLP_GROUP = 'RELATIVE-SLOPE'
    SHADOW_GROUP = 'SHADOW-MASKS'
    STANDARD_GROUP = 'STANDARDISED-PRODUCTS'


class Method(Enum):
    """
    Defines the Interpolation method used for interpolating the
    atmospheric coefficients.
    """

    BILINEAR = 0
    FBILINEAR = 1
    SHEAR = 2
    SHEARB = 3
    RBF = 4


class BrdfParameters(Enum):
    """
    Defines the BRDF Parameters used in BRDF correction.
    """

    ISO = 'ISO'
    VOL = 'VOL'
    GEO = 'GEO'


class ArdProducts(Enum):
    """
    Defines the output ARD products that wagl produces.
    """

    NBAR = 'NBAR'
    NBART = 'NBART'
    LAMBERTIAN = 'LAMBERTIAN'
    SBT = 'SBT'


class Albedos(Enum):
    """
    Defines the albedo labels that wagl uses.
    """

    ALBEDO_0 = '0'
    ALBEDO_TH = 'TH'


class AtmosphericCoefficients(Enum):# param, coeff, vari... what to use
    """
    Defines the atmospheric coefficient names that wagl uses.
    """

    FS = 'FS'
    FV = 'FV'
    A = 'A'
    B = 'B'
    S = 'S'
    DIR = 'DIR'
    DIF = 'DIF'
    TS = 'TS'
    PATH_UP = 'PATH-UP'
    PATH_DOWN = 'PATH-DOWN'
    TRANSMITTANCE_UP = 'TRANSMITTANCE-UP'


class TrackIntersection(Enum):
    """
    Defines the type of track intersection an acquisition
    will have.
    """

    FULL = 0
    PARTIAL = 1
    EMPTY = 2


class PQbits(Enum):
    BAND_1_SATURATED = 0
    BAND_2_SATURATED = 1
    BAND_3_SATURATED = 2
    BAND_4_SATURATED = 3
    BAND_5_SATURATED = 4
    BAND_6_SATURATED = 5
    BAND_7_SATURATED = 6
    CONTIGUITY = 7
    LAND_OBS = 8
    CLOUD_ACCA = 9
    CLOUD_FMASK = 10
    CLOUD_SHADOW_ACCA = 11
    CLOUD_SHADOW_FMASK = 12


# TODO: get rid or redo this class
class PQAConstants(object):

    """
    A class object that contains the majority of constants used throughout
    the PQA process.  Such constants include bands for specific tests, bit
    positions for various tests and thresholds used within various tests.
    """

    def __init__(self, sensor):
        assert sensor is not None
        self.sensor = sensor
        # Initialise everything for immediate access
        self.set_saturation_bands()
        self.set_saturation_bits()
        self.set_acca()
        self.set_fmask()
        self.set_cloud_shadow()
        self.set_test_bits()
        self.set_available_bands()
        self.set_run_cloud_shadow()
        self.set_run_cloud()
        self.set_olitirs()
        self.set_thermal_band()

    def set_saturation_bands(self):
        """
        Get the band numbers associated with saturation tests for a given
        sensor. The band numbers are (to some degree) the band names. This
        may change to be an ordered list, ie 1-n_bands.
        """
        saturation = {
            'TM': ['1', '2', '3', '4', '5', '6', '7'],
            'ETM+': ['1', '2', '3', '4', '5', '61', '62', '7'],
            'OLI_TIRS': ['2', '3', '4', '5', '6', '7', '10', '11'],
            'OLI': ['2', '3', '4', '5', '6', '7'],
            'TIRS': ['10', '11']
        }

        self.saturation_bands = saturation[self.sensor]

    def set_saturation_bits(self):
        """
        Get the relevant bit positions for setting the saturation tests.
        The order should be the same as that returned by the
        set_saturation_bands() function.
        """
        bits = {
            'TM': [0, 1, 2, 3, 4, 5, 7],
            'ETM+': [0, 1, 2, 3, 4, 5, 6, 7],
            'OLI_TIRS': [0, 1, 2, 3, 4, 7, 5, 6],
            'OLI': [0, 1, 2, 3, 4, 7],
            'TIRS': [5, 6]
        }

        self.saturation_bits = bits[self.sensor]

    def set_acca(self):
        """
        Set the threshold constants for the ACCA test.
        """
        # Potentially can configure thresholds per sensor
        self.acca_thresh_f1 = 0.08
        self.acca_thresh_f2 = 0.7
        self.acca_thresh_f3 = 300
        self.acca_thresh_f4 = 225
        self.acca_thresh_f5 = 2
        self.acca_thresh_f6 = 2
        self.acca_thresh_f7 = 1
        self.acca_thresh_f8 = 210
        self.acca_desert_index = 0.5
        self.acca_cold_cloud_pop = 0.4
        self.acca_cold_cloud_mean = 295
        self.acca_thermal_effect = 40.0
        self.acca_snow_threshold = 1

    def set_fmask(self):
        """
        Set the threshold constants for the Fmask test.

        Note: Most of the thresholds are still defined in the python function
        and not here due to licencing constraints.
        """
        # Potentially can configure thresholds per sensor
        self.fmask_cloudprob = 22.5
        # Threshold for water.
        # NB: This seems to miss some clouds over water (which end up having
        # about 35-40% probability, not >50%)
        self.fmask_wclr_max = 50

    def set_cloud_shadow(self):
        """
        Set the threshold constants for the Cloud shadow test.
        """
        # Potentially can configure thresholds per sensor
        self.cshadow_wt_ndvi = 0.1
        self.cshadow_wt_b4 = 0.04
        self.cshadow_wt_b5 = 0.05
        self.cshadow_vrat_th = 0.08
        self.cshadow_btt_th = 293
        self.cshadow_rt_b3 = 0.4
        self.cshadow_rt_b4 = 0.6
        self.cshadow_srt_low = 0.9
        self.cshadow_srt_hi = 1.3
        self.cshadow_lapse_wet = 4.8
        self.cshadow_lapse_standard = 6.4
        self.cshadow_lapse_dry = 9.8
        self.cshadow_stdv_native_bush = 0.04
        self.cshadow_stdv_spectral_flat_water = 0.008
        self.cshadow_mndwi_thresh = 0.1
        self.cshadow_dense_veg = 0.5
        self.cshadow_slope_b34 = 0.11
        self.cshadow_slope_b45 = 0.005
        self.cshadow_slope_b47a = 0.01
        self.cshadow_slope_b47b = 0.05
        self.cshadow_stdv_multiplier = 2.5

    def set_test_bits(self):
        """
        Set the bit positions for each Pixel Quality test.
        """
        self.contiguity = 8
        self.land_sea = 9
        self.acca = 10
        self.fmask = 11
        self.acca_shadow = 12
        self.fmask_shadow = 13
        self.topo_shadow = 14
        self.reserved = 15

    def set_available_bands(self):
        """
        Set the availble bands for a given sensor.
        """
        band_numbers = {
            'TM': ['1', '2', '3', '4', '5', '6', '7'],
            'ETM+': ['1', '2', '3', '4', '5', '61', '62', '7'],
            'OLI_TIRS': ['1', '2', '3', '4', '5', '6', '7', '9', '10', '11'],
            'OLI': ['1', '2', '3', '4', '5', '6', '7', '9'],
            'TIRS': ['10', '11']
        }

        self.available_bands = band_numbers[self.sensor]

    def get_array_band_lookup(self, band_numbers):
        """
        Get the correspoding array indices for a given list of band number
        identifiers. This is only meant to be used wherever
        dataset.ReadAsArray() is used, otherwise the array index lookup
        could be incorrect.
        """
        idx = [self.available_bands.index(bn) for bn in band_numbers]
        return idx

    def set_run_cloud_shadow(self):
        """
        Determine and set (True/False) as to whether or not the cloud shadow
        algorithm will be run. This is so due to the algorithm needing both
        spectral and temperature arrays.
        """
        sensor_list = ['TM', 'ETM+', 'OLI_TIRS']
        self.run_cloud_shadow = bool(self.sensor in sensor_list)

    def set_run_cloud(self):
        """
        Determine and set (True/False) as to whether or not the cloud
        algorithm will be run. This is so due to the algorithm needing
        both spectral and temperature arrays.
        """
        sensor_list = ['TM', 'ETM+', 'OLI_TIRS']
        self.run_cloud = bool(self.sensor in sensor_list)

    def set_olitirs(self):
        """
        Determine and set (True/False) as to whether or not the sensor in
        question is OLI_TIRS. This will be used for both ACCA and the cloud
        shadow algorithm, where the argument input "image_stack" doesn't use
        the coastal aerosol band, but is automatcally read by the
        ReadAsArray() method.
        """
        self.oli_tirs = bool(self.sensor == 'OLI_TIRS')

    def set_thermal_band(self):
        """
        Set the relevant thermal band used for the cloud and cloud shadow
        algorithms. The thermal_band variable will be set to an integer
        corresponding to the band number for a given sensors thermal band.
        If no band is found, then a string is returned.
        """
        self.thermal_band = {
            'TM': '6',
            'ETM+': '61',
            'OLI_TIRS': '10'
        }.get(self.sensor, 'Error! No Thermal Band Found.')

def combine_satellite_sensor(satellite, sensor):
    """
    A small utility to deal with GA's and USGS's various naming
    conventions.
    This joins the two strings into an ugly looking string.
    """
    # NOTE: GA Landsat products use both '-' and '_' as a seperator
    # Remove any occurences of - and _ then convert to lowercase
    satellite_name = re.sub('[-_]', '', satellite).lower()
    sensor_name = re.sub('[-_]', '', sensor).lower()
    return ''.join((satellite_name, sensor_name))


def sbt_bands(satellite, sensor):
    """
    Retrieve the thermal bands to be processed through to SBT for
    a given satellite sensor.
    """
    combined = combine_satellite_sensor(satellite, sensor)

    lookup = {'landsat5tm': ['6'],
              'landsat7etm+': ['61', '62'],
              'landsat8olitirs': ['10']} # band 11 is not stable

    return lookup.get(combined, [])
