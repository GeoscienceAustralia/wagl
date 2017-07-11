"""
Constants
---------
"""
# pylint: disable=attribute-defined-outside-init

from __future__ import absolute_import, print_function
from os.path import join as pjoin
import re
from enum import Enum


ALL_FACTORS = ['fs',
               'fv',
               'a',
               'b',
               's',
               'dir',
               'dif',
               'ts',
               'path-up',
               'path-down',
               'transmittance-up']

ALL_ALBEDOS = ['0', '1', 't', 'th']

POINT_FMT = 'point-{p}'
ALBEDO_FMT = 'albedo-{a}'
POINT_ALBEDO_FMT = ''.join([POINT_FMT, '-', ALBEDO_FMT])

ARD_PRODUCTS = ['brdf', 'terrain', 'lambertian', 'thermal']


class Model(Enum):
    """
    Represents the model workflow that gaip can run.

    *standard* Indicates both NBAR and SBT workflows will run
    *nbar* Indicates NBAR only
    *sbt* Indicates SBT only
    """

    standard = 1
    nbar = 2
    sbt = 3

    @property
    def factors(self):
        """
        Returns the factor names used for interpolation for a given
        Model.<option>.
        """
        fmap = {Model.standard: ALL_FACTORS,
                Model.nbar: ALL_FACTORS[0:8],
                Model.sbt: ALL_FACTORS[8:]}
        return fmap.get(self)

    @property
    def albedos(self):
        """
        Returns the albedo names used for specific Atmospheric
        evaluations for a given Model.<option>.
        """
        amap = {Model.standard: ALL_ALBEDOS,
                Model.nbar: ALL_ALBEDOS[0:-1],
                Model.sbt: [ALL_ALBEDOS[-1]]}
        return amap.get(self)

    @property
    def ard_products(self):
        """
        Returns the ARD products available for a given
        Model.<option>.
        """
        amap = {Model.standard: ARD_PRODUCTS,
                Model.nbar: ARD_PRODUCTS[0:-1],
                Model.sbt: [ARD_PRODUCTS[-1]]}
        return amap.get(self)


class BandType(Enum):
    """
    Represents the Band Type a given acquisition falls under.
    """

    Reflective = 0
    Thermal = 1
    Panchromatic = 2
    Atmosphere = 3
    Quality = 4


class DatasetName(Enum):
    """
    Defines the dataset names or format descriptors, that are used
    for creating and accessing throughout the code base.
    """

    # gaip.ancillary
    ancillary_group = 'ancillary'
    coordinator = 'coordinator'
    dewpoint_temperature = 'dewpoint-temperature'
    temperature_2m = 'temperature-2metre'
    surface_pressure = 'surface-pressure'
    surface_geopotential = 'surface-geopotential-height'
    surface_relative_humidity = 'surface-relative-humidity'
    geopotential = 'geo-potential'
    temperature = 'temperature'
    relative_humidity = 'relative-humidity'
    atmospheric_profile = 'atmospheric-profile'
    aerosol = 'aerosol'
    water_vapour = 'water-vapour'
    ozone = 'ozone'
    elevation = 'elevation'
    brdf_fmt = "BRDF-Band-{band}-{factor}"
    ecmwf_path_fmt = pjoin('{product}', '{year}', 'tif', '{product}_*.tif')

    # gaip.longitude_latitude_arrays
    lon_lat_group = 'longitude-latitude'
    lon = 'longitude'
    lat = 'latitude'

    # gaip.satellite_solar_angles
    sat_sol_group = 'satellite-solar'
    satellite_view = 'satellite-view'
    satellite_azimuth = 'satellite-azimuth'
    solar_zenith = 'solar-zenith'
    solar_azimuth = 'solar-azimuth'
    relative_azimuth = 'relative-azimuth'
    acquisition_time = 'acquisition-time'
    centreline = 'centreline'
    boxline = 'boxline'
    spheroid = 'spheroid'
    orbital_elements = 'orbital-elements'
    satellite_model = 'satellite-model'
    satellite_track = 'satellite-track'

    # gaip.incident_exiting_angles
    inc_exi_group = 'incident-exiting-angles'
    incident = 'incident'
    azimuthal_incident = 'azimuthal-incident'
    exiting = 'exiting'
    azimuthal_exiting = 'azimuthal-exiting'
    relative_slope = 'relative-slope'

    # gaip.reflectance
    standard_group = 'standard-products'
    reflectance_fmt = '{product}/reflectance-band-{band}'

    # gaip.temperature
    temperature_fmt = 'thermal/surface-brightness-temperature-band-{band}'

    # gaip.terrain_shadow_masks
    shadow_group = 'shadow-masks'
    self_shadow = 'self-shadow'
    cast_shadow_fmt = 'cast-shadow-{source}'
    combined_shadow = 'combined-shadow'

    # gaip.slope_aspect
    slp_asp_group = 'slope-aspect'
    slope = 'slope'
    aspect = 'aspect'

    # gaip.dsm
    elevation_group = 'elevation'
    dsm = 'dsm'
    dsm_smoothed = 'dsm-smoothed'

    # gaip.interpolation
    interp_group = 'interpolated-coefficients'
    interpolation_fmt = '{factor}-band-{band}'

    # gaip.modtran
    atmos_inputs_group = 'atmospheric-inputs'
    atmos_results_group = 'atmospheric-results'
    tp5 = 'tp5-data'
    flux = 'flux'
    altitudes = 'altitudes'
    solar_irradiance = 'solar-irradiance'
    upward_radiation_channel = 'upward-radiation-channel'
    downward_radiation_channel = 'downward-radiation-channel'
    channel = 'channel'
    nbar_coefficients = 'nbar-coefficients'
    sbt_coefficients = 'sbt-coefficients'

    # gaip.pq
    pixel_quality = 'pixel-quality/pixel-quality'

    # metadata
    nbar_yaml = 'metadata/nbar-metadata'
    pq_yaml = 'metadata/pq-metadata'
    sbt_yaml = 'metadata/sbt-metadata'


class PQbits(Enum):
    band_1_saturated = 0
    band_2_saturated = 1
    band_3_saturated = 2
    band_4_saturated = 3
    band_5_saturated = 4
    band_6_saturated = 5
    band_7_saturated = 6
    contiguity = 7
    land_obs = 8
    cloud_acca = 9
    cloud_fmask = 10
    cloud_shadow_acca = 11
    cloud_shadow_fmask = 12


# TODO: Re-work this entire file, and the class structures
#       A pain but required as this file contains neccessary
#       band exclusions, and hardwired BRDF wavelength matchups.

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
        if self.sensor in sensor_list:
            self.run_cloud_shadow = True
        else:
            self.run_cloud_shadow = False

    def set_run_cloud(self):
        """
        Determine and set (True/False) as to whether or not the cloud
        algorithm will be run. This is so due to the algorithm needing
        both spectral and temperature arrays.
        """
        sensor_list = ['TM', 'ETM+', 'OLI_TIRS']
        if self.sensor in sensor_list:
            self.run_cloud = True
        else:
            self.run_cloud = False

    def set_olitirs(self):
        """
        Determine and set (True/False) as to whether or not the sensor in
        question is OLI_TIRS. This will be used for both ACCA and the cloud
        shadow algorithm, where the argument input "image_stack" doesn't use
        the coastal aerosol band, but is automatcally read by the
        ReadAsArray() method.
        """
        if self.sensor == 'OLI_TIRS':
            self.oli_tirs = True
        else:
            self.oli_tirs = False

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

def brdf_wavelength_lut(satellite_sensor):
    """
    Retrieves the BRDF wavelengths for a given satellite-sensor.

    :param satellite_sensor:
        A string containing a valid satellite-sensor combination.
        Valid combinations are:
        landsat5tm
        landsat7etm
        landsat8oli
        landsat8olitirs
        sentinel2amsi

    :return:
        A dictionary containing the Band numbers of a sensor as the
        keys, and the BRDF wavelengths as the values.
    """

    input_str = str(satellite_sensor)

    brdf_lut = {'landsat5tm': {'1': '0459_0479nm',
                               '2': '0545_0565nm',
                               '3': '0620_0670nm',
                               '4': '0841_0876nm',
                               '5': '1628_1652nm',
                               '7': '2105_2155nm'},
                'landsat7etm+': {'1': '0459_0479nm',
                                 '2': '0545_0565nm',
                                 '3': '0620_0670nm',
                                 '4': '0841_0876nm',
                                 '5': '1628_1652nm',
                                 '7': '2105_2155nm'},
                'landsat8oli': {'1': '0459_0479nm',
                                '2': '0459_0479nm',
                                '3': '0545_0565nm',
                                '4': '0620_0670nm',
                                '5': '0841_0876nm',
                                '6': '1628_1652nm',
                                '7': '2105_2155nm'},
                'landsat8olitirs': {'1': '0459_0479nm',
                                    '2': '0459_0479nm',
                                    '3': '0545_0565nm',
                                    '4': '0620_0670nm',
                                    '5': '0841_0876nm',
                                    '6': '1628_1652nm',
                                    '7': '2105_2155nm'},
                'sentinel2amsi': {'1': '0459_0479nm',
                                  '2': '0459_0479nm',
                                  '3': '0545_0565nm',
                                  '4': '0620_0670nm',
                                  '5': '0620_0670nm',
                                  '6': '0620_0670nm',
                                  '7': '0841_0876nm',
                                  '8': '0841_0876nm',
                                  '8a': '0841_0876nm',
                                  '9': '0841_0876nm',
                                  '11': '1628_1652nm',
                                  '12': '2105_2155nm'}}.get(input_str, 'Error')

    return brdf_lut


def nbar_bands_lut(satellite_sensor):
    """
    Given a satellite_sensor string, retrieve a list bands to
    process through the NBAR algorithm.

    :param satellite_sensor:
        A string containing a valid satellite-sensor combination.
        Valid combinations are:
        landsat5tm
        landsat7etm
        landsat8oli
        landsat8olitirs
        sentinel2amsi

    :return:
    """

    input_str = str(satellite_sensor)

    nbar_lut = {'landsat5tm': ['1', '2', '3', '4', '5', '7'],
                'landsat7etm+': ['1', '2', '3', '4', '5', '7'],
                'landsat8oli': ['1', '2', '3', '4', '5', '6', '7'],
                'landsat8olitirs': ['1', '2', '3', '4', '5', '6', '7'],
                'sentinel2amsi': ['1',
                                  '2',
                                  '3',
                                  '4',
                                  '5',
                                  '6',
                                  '7',
                                  '8',
                                  '8a',
                                  '9',
                                  '11',
                                  '12']}

    return nbar_lut.get(input_str, 'Error')


def avg_reflectance_lut(satellite_sensor):
    """
    Retrieves the average reflectance values for a given
    satellite-sensor.
    Only those bands processed through NBAR will be returned.

    :param satellite_sensor:
        A string containing a valid satellite-sensor combination.
        Valid combinations are:
        landsat5tm
        landsat7etm
        landsat8oli
        landsat8olitirs
        sentinel2amsi

    :return:
        A dictionary containing the Band numbers of a sensor as the
        keys, and the average reflectance as the values.

    :notes:
        These were copied from the files files in
        /g/data/v10/ULA3-TESTDATA/brdf_modis_band%i.txt. They were
        contained in the last line of those files.
    """

    input_str = str(satellite_sensor)

    avg_reflectance_values = {'landsat5tm': {'1': 0.0365,
                                             '2': 0.0667,
                                             '3': 0.0880,
                                             '4': 0.2231,
                                             '5': 0.2512,
                                             '7': 0.1648},
                              'landsat7etm+': {'1': 0.0365,
                                               '2': 0.0667,
                                               '3': 0.0880,
                                               '4': 0.2231,
                                               '5': 0.2512,
                                               '7': 0.1648},
                              'landsat8oli': {'1': 0.0365,
                                              '2': 0.0365,
                                              '3': 0.0667,
                                              '4': 0.0880,
                                              '5': 0.2231,
                                              '6': 0.2512,
                                              '7': 0.1648},
                              'landsat8olitirs': {'1': 0.0365,
                                                  '2': 0.0365,
                                                  '3': 0.0667,
                                                  '4': 0.0880,
                                                  '5': 0.2231,
                                                  '6': 0.2512,
                                                  '7': 0.1648},
                              'sentinel2amsi': {'1': 0.0365,
                                                '2': 0.0365,
                                                '3': 0.0667,
                                                '4': 0.088,
                                                '5': 0.088,
                                                '6': 0.088,
                                                '7': 0.2231,
                                                '8': 0.2231,
                                                '8a': 0.2231,
                                                '9': 0.2231,
                                                '11': 0.2512,
                                                '12': 0.1648}
                             }

    return avg_reflectance_values.get(input_str, 'Error')


class NBARConstants(object):

    """
    NBAR Constants
    """

    def __init__(self, satellite, sensor):
        """
        Initialise.
        """
        self.satellite = satellite
        self.sensor = sensor
        # NOTE: GA Landsat products use both '-' and '_' as a seperator
        # Remove any occurences of - and _ then convert to lowercase
        satellite_name = re.sub('[-_]', '', self.satellite).lower()
        sensor_name = re.sub('[-_]', '', self.sensor).lower()

        self.sat_sensor = ''.join((satellite_name, sensor_name))

    def get_brdf_lut(self):
        """
        Get the BRDF lookup table.
        """

        brdf_wavelengths = brdf_wavelength_lut(self.sat_sensor)

        return brdf_wavelengths

    def get_brdf_factors(self):
        """
        Get the BRDF factors.
        """
        factors = ['geo', 'iso', 'vol']

        return factors

    def get_nbar_lut(self):
        """
        Get the NBAR lookup table.
        """

        nbar_lut = nbar_bands_lut(self.sat_sensor)

        return nbar_lut

    def get_avg_ref_lut(self):
        """
        Get the average reflectance lookup table.
        """
        avg_ref_lut = avg_reflectance_lut(self.sat_sensor)

        return avg_ref_lut


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
