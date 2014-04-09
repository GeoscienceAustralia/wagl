#! /usr/bin/env python

import logging
import numpy, numexpr
from osgeo import gdal
from ULA3.dataset import SceneDataset
from ULA3.image_processor import ProcessorConfig
from ULA3.image_processor import constants
from ULA3 import DataManager, DataGrid

logger = logging.getLogger('root.' + __name__)

def process(subprocess_list=[], resume=False):
    logger.info('%s.process(%s, %s) called', __name__, subprocess_list, resume)

    CONFIG = ProcessorConfig()
    DATA = DataManager()

    l1t_input_dataset = DATA.get_item(CONFIG.input['l1t']['path'], SceneDataset)
    assert l1t_input_dataset, 'Unable to retrieve SceneDataset object for L1T input scene dataset'
    logger.debug( 'SceneDataset object for %s retrieved', l1t_input_dataset.pathname)

    l1t_stack = DATA.get_item('l1t_stack', numpy.ndarray)
    assert l1t_stack is not None, 'Unable to retrieve ndarray object for l1t_stack'
    logger.debug( 'ndarray object for l1t_stack retrieved')


    def LandsatTemperature(image_stack, input_dataset):
        """
        Converts a Landsat TM/ETM+ thermal band into degrees Kelvin.

        Required input is the image to be in byte scaled DN form (0-255).

        :param band_array:
            A 2D Numpy array containing the thermal band.

        :return:
            A 2D Numpy array containing degrees Kelvin.

        :author:
            Josh Sixsmith, joshua.sixsmith@ga.gov.au
        """

        def radiance_conversion(band_array, gain, bias):
            """
            Converts the input image into radiance.

            Two methods could be used; the gain and bias or the spectral
            radiance scaling. Defined to use the gain and bias method.

            Gain and bias method:
            B6_gain = (LMAX_BAND6 - LMIN_BAND6) / (QCALMAX_Band6 - QCALMIN_Band6)
            B6_bias = LMIN_BAND6 - (B6_gain * QCALMIN_Band6)
            rad = gain * image + bias

            Spectral Radiance Scaling method
            rad = ((LMAX - LMIN)/(QCALMAX - QCALMIN)) * (image - QCALMIN) + LMIN

            :param band_array:
                The Thermal band to be converted to radiance, requires DN.

            :param gain:
                Floating point value.

            :param bias:
                Floating point value.

            :return:
                The thermal band converted to at-sensor radiance in
                watts/(meter squared * ster * um) as a 2D Numpy array.
            """

            logger.debug('gain = %f, bias = %f', gain, bias)

            return numexpr.evaluate('gain * band_array + bias')


        def temperature_conversion(band_array, k1, k2):
            """
            Converts the radiance image to degrees Kelvin.

            :param image:
                A 2D Numpy array containing the thermal band converted to
                radiance.

            :param met_data:
                A dictionary containg the metadata of the input image.

            :return:
                A 2D Numpy array of the thermal band coverted to at-sensor
                degrees Kelvin.
            """

            logger.debug('k1 = %f, k2 = %f', k1, k2)

            return k2/(numpy.log(k1/band_array + 1))



        #------------------------Processing Goes Here-----------------------------

        pq_const = constants.pqaContants(l1t_input_dataset.sensor)

        #full_band_list = sorted(set(l1t_input_dataset.bands('REFLECTIVE')) |
        #                    set(l1t_input_dataset.bands('THERMAL')))
        full_band_list = pq_const.available_bands

        # Find index of band 62 for LS7, band 60 for LS5 ***Should be band 61 for LS7
        #thermal_band = max(l1t_input_dataset.bands('THERMAL')) # TODO: May break with LS8
        #thermal_band_index = full_band_list.index(thermal_band)
        # Get Band 6 for TM, Band 61 for ETM+ and Band 10 for OLI_TIRS
        thermal_band = pq_const.thermal_band
        thermal_band_index = pq_const.getArrayBandLookup([thermal_band])

        logger.debug('thermal_band = %d, thermal_band_index = %d', thermal_band, thermal_band_index)

        radiance_array = radiance_conversion(l1t_stack[thermal_band_index],
                                             l1t_input_dataset.gain[thermal_band],
                                             l1t_input_dataset.bias[thermal_band])

        kelvin_array = temperature_conversion(radiance_array,
                                              l1t_input_dataset.satellite.k[0],
                                              l1t_input_dataset.satellite.k[1])

        return kelvin_array.astype('float32')

    # Compute and store temperature value
    kelvin_grid = DataGrid(array=LandsatTemperature(l1t_stack, l1t_input_dataset))
    DATA.set_item('kelvin.tif', kelvin_grid)

#    del l1t_stack; # Remove l1t_stack - not needed downstream
##### Oh yes it is - see (for instance) saturaton_masking.py
#    DATA.free_item('l1t_stack', numpy.ndarray)
#                   save_first=CONFIG.DEBUG)


