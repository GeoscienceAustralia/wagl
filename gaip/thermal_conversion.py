"""
Conversion Routines
-------------------
"""
from __future__ import absolute_import, print_function
import logging
import numpy
import numexpr


def radiance_conversion(band_array, gain, bias):
    """
    Converts the input image into radiance using the gain and bias
    method.

    :param band_array:
        A `NumPy` array containing the scaled DN to be converted
        to radiance at sensor.

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

    :param k1:
        Conversion constant 1.

    :param k2:
        Conversion constant 2.

    :return:
        A 2D Numpy array of the thermal band coverted to at-sensor
        degrees Kelvin.
    """

    logging.debug('k1 = %f, k2 = %f', k1, k2)

    return k2 / (numpy.log(k1 / band_array + 1))


def get_landsat_temperature(l1t_stack, acquisitions, pq_const):
    """
    Converts a Landsat TM/ETM+ thermal band into degrees Kelvin.
    Required input is the image to be in byte scaled DN form (0-255).

    :param l1t_stack:
        A 3D `numpy.ndarray` containing the thermal band.

    :param acquisitions:
        A list of acquisition instances.

    :param pq_const:
        An instance of the PQ constants.

    :return:
        A 2D Numpy array containing degrees Kelvin.
    """
    acqs = acquisitions
    thermal_band = pq_const.thermal_band

    if type(thermal_band) == str:
        kelvin_array = numpy.zeros((l1t_stack.shape[1],
                                    l1t_stack.shape[2]), dtype='float32')
        return kelvin_array

    # Function returns a list of one item. Take the first item.
    thermal_band_index = pq_const.get_array_band_lookup([thermal_band])[0]

    logging.debug('thermal_band = %d, thermal_band_index = %d',
                  thermal_band, thermal_band_index)

    radiance_array = radiance_conversion(l1t_stack[thermal_band_index],
                                         acqs[thermal_band_index].gain,
                                         acqs[thermal_band_index].bias)

    kelvin_array = temperature_conversion(radiance_array,
                                          acqs[thermal_band_index].K1,
                                          acqs[thermal_band_index].K2)

    return kelvin_array.astype('float32')
