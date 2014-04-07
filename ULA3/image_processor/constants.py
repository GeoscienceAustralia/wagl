#!/usr/bin/env python

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
        self.setBandIndexLookup()

    def setSaturationBands(self):
        """
        Get the band numbers associated with saturation tests for a given sensor.
        The band numbers are sequential ordering, eg ETM+ Band61, Band62, Band7 are 6, 7, 8.
        """
        saturation = {
                     'TM' : [1,2,3,4,5,6,7],
                     'ETM+' : [1,2,3,4,5,6,7,8],
                     'OLI_TIRS' : [2,3,4,5,6,7,10,11]
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
               'OLI_TIRS' : [0,1,2,3,4,7,5,6]
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
       self.cloud_shadow_wt_ndvi        = 0.1
       self.cloud_shadow_wt_b4          = 0.04
       self.cloud_shadow_wt_b5          = 0.05
       self.cloud_shadow_vrat_th        = 0.08
       self.cloud_shadow_btt_th         = 293
       self.cloud_shadow_rt_b3          = 0.4
       self.cloud_shadow_rt_b4          = 0.6
       self.cloud_shadow_srt_low        = 0.9
       self.cloud_shadow_srt_hi         = 1.3
       self.cloud_shadow_lapse_wet      = 4.8
       self.cloud_shadow_lapse_standard = 6.4
       self.cloud_shadow_lapse_dry      = 9.8

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

    def setBandIndexLookup(self):
        """
        Set the array index lookup for the bands from a given sensor.
        """
        band_numbers = {
                       'TM' : [1,2,3,4,5,6,7],
                       'ETM+' : [1,2,3,4,5,61,62,7],
                       'OLI_TIRS' : [1,2,3,4,5,6,7,9,10,11],
                       'OLI' : [1,2,3,4,5,6,7,9]
                       'TIRS' : [10,11]
                       }

        self.band_index_lookup = band_numbers[self.sensor]
