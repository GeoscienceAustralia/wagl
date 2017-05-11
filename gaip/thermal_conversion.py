#!/usr/bin/env python

"""
Various routines for converting radiance to temperature
-------------------------------------------------------
"""

from __future__ import absolute_import, print_function
import logging
import numpy
import numexpr
import h5py

from gaip.constants import DatasetName
from gaip.hdf5 import dataset_compression_kwargs
from gaip.hdf5 import attach_image_attributes
from gaip.metadata import create_ard_yaml
from gaip.tiling import generate_tiles


def _surface_brightness_temperature(acquisition, bilinear_fname,
                                    ancillary_fname, out_fname, compression,
                                    y_tile):
    """
    A private wrapper for dealing with the internal custom workings of the
    NBAR workflow.
    """
    band_num = acquisition.band_num
    dname_fmt = DatasetName.interpolation_fmt.value
    with h5py.File(bilinear_fname, 'r') as fid:
        dname = dname_fmt.format(factor='path-up', band=band_num)
        upwelling_dset = fid[dname]
        dname = dname_fmt.format(factor='transmittance-up', band=band_num)
        transmittance_dset = fid[dname]

        kwargs = {'acquisition': acquisition,
                  'upwelling_radiation': upwelling_dset,
                  'transmittance': transmittance_dset,
                  'out_fname': out_fname,
                  'compression': compression,
                  'y_tile': y_tile}

        rfid = surface_brightness_temperature(**kwargs)

    create_ard_yaml(acquisition, ancillary_fname, rfid, True)

    rfid.close()
    return


def surface_brightness_temperature(acquisition, upwelling_radiation,
                                   transmittance, out_fname=None,
                                   compression='lzf', y_tile=100):
    """
    Convert Thermal acquisition to Surface Brightness Temperature.

    T[Kelvin] = k2 / ln( 1 + (k1 / I[0]) )

    where T is the surface brightness temperature (the surface temperature
    if the surface is assumed to be an ideal black body i.e. unit emissivity),
    k1 & k2 are calibration constants specific to the platform/sensor/band,
    and I[0] is the surface radiance (the integrated band radiance, in
    Watts per square metre per steradian per thousand nanometres).

    I = t I[0] + d

    where I is the radiance at the sensor, t is the transmittance (through
    the atmosphere), and d is radiance from the atmosphere itself.

    :param acquisition:
        An instance of an acquisition object.

    :param upwelling_radiation:
        A `NumPy` or `NumPy` like dataset that allows indexing
        and returns a `NumPy` dataset containing the MODTRAN
        factor `upwelling_radiation` data values when indexed/sliced.
        
    :param transmittance:
        A `NumPy` or `NumPy` like dataset that allows indexing
        and returns a `NumPy` dataset containing the MODTRAN
        factor `upwelling_transmittance` data values when indexed/sliced.

    :param out_fname:
        If set to None (default) then the results will be returned
        as an in-memory hdf5 file, i.e. the `core` driver.
        Otherwise it should be a string containing the full file path
        name to a writeable location on disk in which to save the HDF5
        file.

        The dataset names will be as follows:

        * surface-brightness-temperature-band-{number}

    :param compression:
        The compression filter to use. Default is 'lzf'.
        Options include:

        * 'lzf' (Default)
        * 'lz4'
        * 'mafisc'
        * An integer [1-9] (Deflate/gzip)

    :param y_tile:
        Defines the tile size along the y-axis. Default is 100.

    :return:
        An opened `h5py.File` object, that is either in-memory using the
        `core` driver, or on disk.
    """
    acq = acquisition
    geobox = acq.gridded_geo_box()

    # tiling scheme
    tiles = generate_tiles(acq.samples, acq.lines, acq.samples, y_tile)

    # Initialise the output file
    if out_fname is None:
        fid = h5py.File('surface-temperature.h5', driver='core',
                        backing_store=False)
    else:
        fid = h5py.File(out_fname, 'w')

    kwargs = dataset_compression_kwargs(compression=compression,
                                        chunks=(1, acq.samples))
    kwargs['shape'] = (acq.lines, acq.samples)
    kwargs['fillvalue'] = -999
    kwargs['dtype'] = 'float32'

    # attach some attributes to the image datasets
    attrs = {'crs_wkt': geobox.crs.ExportToWkt(),
             'geotransform': geobox.transform.to_gdal(),
             'no_data_value': kwargs['fillvalue'],
             'sattelite': acq.spacecraft_id,
             'sensor': acq.sensor_id,
             'band number': acq.band_num}

    name_fmt = DatasetName.temperature_fmt.value
    dataset_name = name_fmt.format(band=acq.band_num)
    out_dset = fid.create_dataset(dataset_name, **kwargs)

    desc = "Surface Brightness Temperature in Kelvin."
    attrs['Description'] = desc
    attach_image_attributes(out_dset, attrs)

    # constants
    k1 = acq.K1
    k2 = acq.K2

    # process each tile
    for tile in tiles:
        idx = (slice(tile[0][0], tile[0][1]), slice(tile[1][0], tile[1][1]))

        acq_args = {'window': tile,
                    'masked': False,
                    'apply_gain_offset': acq.scaled_radiance,
                    'out_no_data': kwargs['fillvalue']}

        radiance = acq.data(**acq_args)
        path_up = upwelling_radiation[idx]
        trans = transmittance[idx]
        expr = "(radiance-path_up) / trans"
        corrected_radiance = numexpr.evaluate(expr)
        mask = corrected_radiance <= 0
        expr = "k2 / log(k1 / corrected_radiance + 1)"
        brightness_temp = numexpr.evaluate(expr)
        brightness_temp[mask] = kwargs['fillvalue']

        out_dset[idx] = brightness_temp

    return fid


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


def get_landsat_temperature(acquisitions, pq_const):
    """
    Converts a Landsat TM/ETM+ thermal band into degrees Kelvin.
    Required input is the image to be in byte scaled DN form (0-255).

    :param acquisitions:
        A list of acquisition instances.

    :param pq_const:
        An instance of the PQ constants.

    :return:
        A 2D Numpy array containing degrees Kelvin.
    """
    acqs = acquisitions
    thermal_band = pq_const.thermal_band

    # Function returns a list of one item. Take the first item.
    acq = acqs[pq_const.get_array_band_lookup([thermal_band])[0]]
    radiance = acq.data(apply_gain_offset=True)

    kelvin_array = temperature_conversion(radiance, acq.K1, acq.K2)

    return kelvin_array.astype('float32')
