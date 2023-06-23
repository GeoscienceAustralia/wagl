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

from wagl.constants import DatasetName, GroupName, ArdProducts
from wagl.constants import AtmosphericCoefficients as AC
from wagl.hdf5 import H5CompressionFilter, attach_image_attributes
from wagl.metadata import create_ard_yaml

NO_DATA_VALUE = -999


def _surface_brightness_temperature(
    acquisition,
    acquisitions,
    bilinear_fname,
    ancillary_fname,
    out_fname,
    normalized_solar_zenith,
    compression=H5CompressionFilter.LZF,
    filter_opts=None,
):
    """
    A private wrapper for dealing with the internal custom workings of the
    NBAR workflow.
    """
    with h5py.File(bilinear_fname, "r") as interp_fid, h5py.File(
        ancillary_fname, "r"
    ) as fid_anc, h5py.File(out_fname, "w") as fid:
        grp1 = interp_fid[GroupName.INTERP_GROUP.value]
        surface_brightness_temperature(acquisition, grp1, fid, compression, filter_opts)

        grp2 = fid_anc[GroupName.ANCILLARY_GROUP.value]
        create_ard_yaml(acquisitions, grp2, fid, normalized_solar_zenith, True)


def surface_brightness_temperature(
    acquisition,
    interpolation_group,
    out_group=None,
    compression=H5CompressionFilter.LZF,
    filter_opts=None,
):
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

    :param interpolation_group:
        The root HDF5 `Group` that contains the interpolated
        atmospheric coefficients.
        The dataset pathnames are given by the following string format:

        * DatasetName.INTERPOLATION_FMT

    :param out_group:
        If set to None (default) then the results will be returned
        as an in-memory hdf5 file, i.e. the `core` driver. Otherwise,
        a writeable HDF5 `Group` object.

        The dataset names will be given by the format string detailed
        by:

        * DatasetName.TEMPERATURE_FMT

    :param compression:
        The compression filter to use.
        Default is H5CompressionFilter.LZF

    :filter_opts:
        A dict of key value pairs available to the given configuration
        instance of H5CompressionFilter. For example
        H5CompressionFilter.LZF has the keywords *chunks* and *shuffle*
        available.
        Default is None, which will use the default settings for the
        chosen H5CompressionFilter instance.

    :return:
        An opened `h5py.File` object, that is either in-memory using the
        `core` driver, or on disk.

    :notes:
        This function used to accept `NumPy` like datasets as inputs,
        but as this functionality was never used, it was simpler to
        parse through the H5 Group object, which in most cases
        reduced the number or parameters being parsed through.
        Thereby simplifying the overall workflow, and making it
        consistant with other functions within the overall workflow.
    """
    acq = acquisition
    geobox = acq.gridded_geo_box()
    bn = acq.band_name

    # retrieve the upwelling radiation and transmittance datasets
    dname_fmt = DatasetName.INTERPOLATION_FMT.value
    dname = dname_fmt.format(coefficient=AC.PATH_UP.value, band_name=bn)
    upwelling_radiation = interpolation_group[dname]
    dname = dname_fmt.format(coefficient=AC.TRANSMITTANCE_UP.value, band_name=bn)
    transmittance = interpolation_group[dname]

    # Initialise the output file
    if out_group is None:
        fid = h5py.File("surface-temperature.h5", "w", driver="core", backing_store=False)
    else:
        fid = out_group

    if GroupName.STANDARD_GROUP.value not in fid:
        fid.create_group(GroupName.STANDARD_GROUP.value)

    if filter_opts is None:
        filter_opts = {}
    else:
        filter_opts = filter_opts.copy()
    filter_opts["chunks"] = acq.tile_size

    group = fid[GroupName.STANDARD_GROUP.value]
    kwargs = compression.config(**filter_opts).dataset_compression_kwargs()
    kwargs["shape"] = (acq.lines, acq.samples)
    kwargs["fillvalue"] = NO_DATA_VALUE
    kwargs["dtype"] = "float32"

    # attach some attributes to the image datasets
    attrs = {
        "crs_wkt": geobox.crs.ExportToWkt(),
        "geotransform": geobox.transform.to_gdal(),
        "no_data_value": kwargs["fillvalue"],
        "platform_id": acq.platform_id,
        "sensor_id": acq.sensor_id,
        "band_id": acq.band_id,
        "band_name": bn,
        "alias": acq.alias,
    }

    name_fmt = DatasetName.TEMPERATURE_FMT.value
    dataset_name = name_fmt.format(product=ArdProducts.SBT.value, band_name=acq.band_name)
    out_dset = group.create_dataset(dataset_name, **kwargs)

    desc = "Surface Brightness Temperature in Kelvin."
    attrs["description"] = desc
    attach_image_attributes(out_dset, attrs)

    # pylint: disable=unused-variable
    # constants
    k1 = acq.K1  # noqa: F841
    k2 = acq.K2  # noqa: F841

    # process each tile
    for tile in acq.tiles():
        idx = (slice(tile[0][0], tile[0][1]), slice(tile[1][0], tile[1][1]))

        radiance = acq.radiance_data(window=tile, out_no_data=NO_DATA_VALUE)  # noqa: F841
        path_up = upwelling_radiation[idx]  # noqa: F841
        trans = transmittance[idx]
        mask = ~numpy.isfinite(trans)
        expr = "(radiance - path_up) / trans"
        corrected_radiance = numexpr.evaluate(expr)
        mask |= corrected_radiance <= 0
        expr = "k2 / log(k1 / corrected_radiance + 1)"
        brightness_temp = numexpr.evaluate(expr)
        brightness_temp[mask] = kwargs["fillvalue"]

        out_dset[idx] = brightness_temp
    acq.close()  # If dataset is cached; clear it

    if out_group is None:
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

    logging.debug("gain = %f, bias = %f", gain, bias)

    return numexpr.evaluate("gain * band_array + bias")


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

    logging.debug("k1 = %f, k2 = %f", k1, k2)

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
    acq = [a for a in acqs if a.band_id == thermal_band][0]
    radiance = acq.radiance_data()

    kelvin_array = temperature_conversion(radiance, acq.K1, acq.K2)

    return kelvin_array.astype("float32")


def temperature_at_sensor(thermal_acquisition, window=None):
    """
    Given a thermal acquisition, convert to at sensor temperature
    in Kelivn.

    :param thermal_acquisition:
        An acquisition with a band_type of BandType.Thermal.

    :param window:
        Defines a subset ((ystart, yend), (xstart, xend)) in array
        co-ordinates. Default is None.

    :return:
        A `NumPy` array of whose shape is given by:
        (thermal_acquisition.lines, thermal_acquisition.samples)
        or the dimensions given by the `window` parameter.
    """
    k1 = thermal_acquisition.K1  # pylint: disable=unused-variable # noqa: F841
    k2 = thermal_acquisition.K2  # pylint: disable=unused-variable # noqa: F841

    # pylint: disable=unused-variable
    data = thermal_acquisition.radiance_data(window=window)  # noqa: F841
    result = numexpr.evaluate("k2 / (log(k1 / data + 1))")

    return result
