#! /usr/bin/env python

import logging
import numpy
import numexpr

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

    logging.debug('gain = %f, bias = %f', gain, bias)

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

    logging.debug('k1 = %f, k2 = %f', k1, k2)

    return k2/(numpy.log(k1/band_array + 1))

def get_landsat_temperature(l1t_stack, acquisitions, pq_const):
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

    full_band_list = pq_const.available_bands

    # Find index of band 62 for LS7, band 60 for LS5 ***Should be band 61 for LS7
    #thermal_band = max(l1t_input_dataset.bands('THERMAL')) # TODO: May break with LS8
    #thermal_band_index = full_band_list.index(thermal_band)
    # Get Band 6 for TM, Band 61 for ETM+ and Band 10 for OLI_TIRS
    thermal_band = pq_const.thermal_band

    if (type(thermal_band) == str):
        logging.debug('No thermal band defined in constants.py for sensor %s.' \
            ' Generating a blank float32 array.', l1t_input_dataset.sensor)
        kelvin_array = numpy.zeros((l1t_stack.shape[1],l1t_stack.shape[2]), dtype='float32')
        return kelvin_array

    # Function returns a list of one item. Take the first item.
    thermal_band_index = pq_const.getArrayBandLookup([thermal_band])[0] 

    logging.debug('thermal_band = %d, thermal_band_index = %d', thermal_band, thermal_band_index)

    radiance_array = radiance_conversion(l1t_stack[thermal_band_index],
        acquisitions[thermal_band_index].gain,
        acquisitions[thermal_band_index].bias)
#        l1t_input_dataset.gain[pq_const.band_num_sequence[thermal_band]],
#        l1t_input_dataset.bias[pq_const.band_num_sequence[thermal_band]])

    kelvin_array = temperature_conversion(radiance_array,
        acquisitions[thermal_band_index].K1,
        acquisitions[thermal_band_index].K2)
#        l1t_input_dataset.satellite.k[0],
#        l1t_input_dataset.satellite.k[1])

    return kelvin_array.astype('float32')
