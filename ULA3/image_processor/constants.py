#!/usr/bin/env python

import re

class pqaContants:
    """
    A Class object that contains the majority of constants used throughout the PQA process.
    Such constants include bands for specific tests, bit positions for various tests and thresholds used
    within various tests.
    """
    def __init__(self, sensor):
        self.sensor = sensor
        # Initialise everything for immediate access
        self.setSaturationBands()
        self.setSaturationBits()
        self.setACCA()
        self.setFmask()
        self.setCloudShadow()
        self.setTestBits()
        self.setAvailableBands()
        #self.setBandIndexLookup()
        self.setRunCloudShadow()
        self.setRunCloud()
        self.setOLITIRS()
        self.setThermalBand()
        self.setBandNumSequence()

    def setSaturationBands(self):
        """
        Get the band numbers associated with saturation tests for a given sensor.
        The band numbers are (to some degree) the band names. This may change to
        be an ordered list, ie 1-n_bands.
        """
        saturation = {
                     'TM' : [1,2,3,4,5,6,7],
                     'ETM+' : [1,2,3,4,5,61,62,7],
                     'OLI_TIRS' : [2,3,4,5,6,7,10,11],
                     'OLI' : [2,3,4,5,6,7],
                     'TIRS' : [10,11]
                     }

        self.saturation_bands = saturation[self.sensor]

    def setSaturationBits(self):
        """
        Get the relevant bit positions for setting the saturation tests.
        The order should be the same as that returned by the setSaturationBands() function.
        """
        bits = {
               'TM' : [0,1,2,3,4,5,7],
               'ETM+' : [0,1,2,3,4,5,6,7],
               'OLI_TIRS' : [0,1,2,3,4,7,5,6],
               'OLI' : [0,1,2,3,4,7],
               'TIRS' : [6,6]
               }

        self.saturation_bits = bits[self.sensor]

    def setACCA(self):
       """
       Set the threshold constants for the ACCA test.
       """
       # Potentially can configure thresholds per sensor
       self.acca_thresh_f1      = 0.08
       self.acca_thresh_f2      = 0.7
       self.acca_thresh_f3      = 300
       self.acca_thresh_f4      = 225
       self.acca_thresh_f5      = 2
       self.acca_thresh_f6      = 2
       self.acca_thresh_f7      = 1
       self.acca_thresh_f8      = 210
       self.acca_desertIndex    = 0.5
       self.acca_coldCloud_pop  = 0.4
       self.acca_coldCloud_mean = 295
       self.acca_thermal_effect = 40.0
       self.acca_snow_threshold = 1

    def setFmask(self):
       """
       Set the threshold constants for the Fmask test.

       Note: Most of the thresholds are still defined in the python function and not here due to licencing constraints.
       """
       # Potentially can configure thresholds per sensor
       self.fmask_cloudprob = 22.5
       # Threshold for water.
       # NB: This seems to miss some clouds over water (which end up having about 35-40% probability, not >50%)
       self.fmask_wclr_max = 50

    def setCloudShadow(self):
       """
       Set the threshold constants for the Cloud shadow test.
       """
       # Potentially can configure thresholds per sensor
       self.cshadow_wt_ndvi                  = 0.1
       self.cshadow_wt_b4                    = 0.04
       self.cshadow_wt_b5                    = 0.05
       self.cshadow_vrat_th                  = 0.08
       self.cshadow_btt_th                   = 293
       self.cshadow_rt_b3                    = 0.4
       self.cshadow_rt_b4                    = 0.6
       self.cshadow_srt_low                  = 0.9
       self.cshadow_srt_hi                   = 1.3
       self.cshadow_lapse_wet                = 4.8
       self.cshadow_lapse_standard           = 6.4
       self.cshadow_lapse_dry                = 9.8
       self.cshadow_stdv_native_bush         = 0.04
       self.cshadow_stdv_spectral_flat_water = 0.008
       self.cshadow_mndwi_thresh             = 0.1
       self.cshadow_dense_veg                = 0.5
       self.cshadow_slope_b34                = 0.11
       self.cshadow_slope_b45                = 0.005
       self.cshadow_slope_b47a               = 0.01
       self.cshadow_slope_b47b               = 0.05
       self.cshadow_stdv_multiplier          = 2.5

    def setTestBits(self):
        """
        Set the bit positions for each Pixel Quality test.
        """
        self.contiguity   = 8
        self.land_sea     = 9
        self.acca         = 10
        self.fmask        = 11
        self.acca_shadow  = 12
        self.fmask_shadow = 13
        self.topo_shadow  = 14
        self.reserved     = 15

    def setAvailableBands(self):
        """
        Set the availble bands for a given sensor.
        """
        band_numbers = {
                       'TM' : [1,2,3,4,5,6,7],
                       'ETM+' : [1,2,3,4,5,61,62,7],
                       'OLI_TIRS' : [1,2,3,4,5,6,7,9,10,11],
                       'OLI' : [1,2,3,4,5,6,7,9],
                       'TIRS' : [10,11]
                       }

        self.available_bands = band_numbers[self.sensor]

    def getArrayBandLookup(self, band_numbers):
        """
        Get the correspoding array indices for a given list of band number identifiers.
        This is only meant to be used wherever dataset.ReadAsArray() is used, otherwise the 
        array index lookup could be incorrect.
        """
        idx = [self.available_bands.index(bn) for bn in band_numbers]
        return idx

    def setRunCloudShadow(self):
        """
        Determine and set (True/False) as to whether or not the cloud shadow algorithm will be run.
        This is so due to the algorithm needing both spectral and temperature arrays.
        """
        sensor_list = ['TM','ETM+','OLI_TIRS']
        if self.sensor in sensor_list:
            self.run_cloud_shadow = True
        else:
            self.run_cloud_shadow = False

    def setRunCloud(self):
        """
        Determine and set (True/False) as to whether or not the cloud algorithm will be run.
        This is so due to the algorithm needing both spectral and temperature arrays.
        """
        sensor_list = ['TM','ETM+','OLI_TIRS']
        if self.sensor in sensor_list:
            self.run_cloud = True
        else:
            self.run_cloud = False

    def setOLITIRS(self):
        """
        Determine and set (True/False) as to whether or not the sensor in question is OLI_TIRS.
        This will be used for both ACCA and the cloud shadow algorithm, where the argument input
        "image_stack" doesn't use the coastal aerosol band, but is automatcally read by the
        ReadAsArray() method.
        """
        if self.sensor == 'OLI_TIRS':
            self.oli_tirs = True
        else:
            self.oli_tirs = False

    def setThermalBand(self):
        """
        Set the relevant thermal band used for the cloud and cloud shadow algorithms.
        The thermal_band variable will be set to an integer corresponding to the band number
        for a given sensors thermal band.  If no band is found, then a string is returned.
        """
        self.thermal_band = {
                            'TM' : 6,
                            'ETM+' : 61,
                            'OLI_TIRS' : 10
                            }.get(self.sensor, 'Error! No Thermal Band Found.')

    def setBandNumSequence(self):
        """
        Set the band numbering sequence 1:n for a given sensor.
        This is intended to be used with the SceneDataset class's
        gain and bias method.
        """
        sequence = {
                   'TM' : {
                           1 : 1,
                           2 : 2,
                           3 : 3,
                           4 : 4,
                           5 : 5,
                           6 : 6,
                           7 : 7
                          },
                   'ETM+' : {
                             1 : 1,
                             2 : 2,
                             3 : 3,
                             4 : 4,
                             5 : 5,
                             61 : 6,
                             62 : 7,
                             7 : 8
                            },
                   'OLI_TIRS' : {
                                 1 : 1,
                                 2 : 2,
                                 3 : 3,
                                 4 : 4,
                                 5 : 5,
                                 6 : 6,
                                 7 : 7,
                                 8 : 8,
                                 9 : 9,
                                 10 : 10,
                                 11 : 11
                                },
                   'OLI' : {
                            1 : 1,
                            2 : 2,
                            3 : 3,
                            4 : 4,
                            5 : 5,
                            6 : 6,
                            7 : 7,
                            8 : 8,
                            9 : 9
                           },
                   'TIRS' : {
                             10 : 10,
                             11 : 11
                            }
                   }
        self.band_num_sequence = sequence[self.sensor]


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

    :return:
        A dictionary containing the Band numbers of a sensor as the
        keys, and the BRDF wavelengths as the values.
    """

    input_str = str(satellite_sensor)

    BRDF_LUT = {
        'landsat5tm' : { 1 : '0459_0479nm',
                         2 : '0545_0565nm',
                         3 : '0620_0670nm',
                         4 : '0841_0876nm',
                         5 : '1628_1652nm',
                         7 : '2105_2155nm'
                       },
        'landsat7etm+' : { 1 : '0459_0479nm',
                          2 : '0545_0565nm',
                          3 : '0620_0670nm',
                          4 : '0841_0876nm',
                          5 : '1628_1652nm',
                          7 : '2105_2155nm'
                        },
        'landsat8oli' : { 1 : '0459_0479nm',
                          2 : '0459_0479nm',
                          3 : '0545_0565nm',
                          4 : '0620_0670nm',
                          5 : '0841_0876nm',
                          6 : '1628_1652nm',
                          7 : '2105_2155nm'
                        },
        'landsat8olitirs' : { 1 : '0459_0479nm',
                               2 : '0459_0479nm',
                               3 : '0545_0565nm',
                               4 : '0620_0670nm',
                               5 : '0841_0876nm',
                               6 : '1628_1652nm',
                               7 : '2105_2155nm'
                            }
               }.get(input_str, 'Error')

    return BRDF_LUT


class NBARConstants:
    """
    
    """

    def __init__(self, satellite, sensor):
        """
        
        """
        self.satellite = satellite
        self.sensor = sensor

    def getBRDFlut(self):
        """
        
        """
        # NOTE: GA Landsat products use both '-' and '_' as a seperator
        # Remove any occurences of - and _ then convert to lowercase
        satellite_name = re.sub('[-_]', '', self.satellite).lower()
        sensor_name = re.sub('[-_]', '', self.sensor).lower()
        
        sat_sensor = ''.join((satellite_name, sensor_name))
        
        brdf_wavelengths = brdf_wavelength_lut(sat_sensor)

        return brdf_wavelengths

    def getBRDFfactors(self):
        """
        
        """
        factors = ['geo', 'iso', 'vol']

        return factors
