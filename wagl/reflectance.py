#!/usr/bin/env python

"""
Calculates the Lambertian, BRDF corrected and BRDF + Terrain corrected
----------------------------------------------------------------------

reflectance
-----------
"""

from __future__ import absolute_import, print_function
import numpy
import numexpr
import h5py
from scipy import ndimage

from wagl.constants import DatasetName, GroupName, BrdfDirectionalParameters
from wagl.constants import AtmosphericCoefficients as AC
from wagl.constants import ArdProducts as AP
from wagl.data import as_array
from wagl.hdf5 import H5CompressionFilter, attach_image_attributes
from wagl.hdf5 import create_external_link, find
from wagl.metadata import create_ard_yaml
from wagl.__surface_reflectance import reflectance

NO_DATA_VALUE = -999
NAN = numpy.nan


def _sequential_valid_rows(mask):
    """
    Check that the mask of null data contains sequential rows.
    i.e. non-sequential means an invalid row is surrounded by valid
    rows.
    """
    nrows = mask.shape[0]
    row_ids = numpy.arange(nrows)

    # identify any rows that are completely null
    invalid_rows_mask = numpy.all(mask, axis=1)

    valid_rows = row_ids[~invalid_rows_mask]
    start_idx = valid_rows[0]
    end_idx = valid_rows[-1] + 1

    # this section is probably not required anymore
    if (start_idx == 0) and (end_idx == nrows):
        all_valid = True
    else:
        all_valid = False

    # check that we are ascending by 1
    # i.e. an invalid row has valid rows before and after
    # TODO;
    # don't raise, simply return and we use another method to fill nulls
    sequential = numpy.all(numpy.diff(valid_rows) == 1)
    if not sequential:
        msg = "Rows with valid data are non-sequential."
        raise Exception(msg)

    return all_valid, start_idx, end_idx


def _fill_nulls(data, mask):
    """
    Calculate run-length averages and insert (inplace) at null pixels.
    In future this could also be an averaging kernel.
    We are assuming 2D arrays only (y, x).
    """
    # TODO: how to account for a row of data that is all null?
    # potentially use alternate methods if row length average doesn't satisfy
    for row in range(data.shape[0]):
        data[row][mask[row]] = numpy.mean(data[row][~mask[row]])


def scale_reflectance(data, clip_range=(1, 10000), clip=True):
    """
    Scale reflectance data to the range [1, 10000], with a null
    value of -999, and a datatype of int16.
    Scaling follows the formulae used in the original f90 code:

        'data * 10000 + 0.5'
    
    The data is also clipped to the range [1, 10000].

    :param data:
        A NumPy array to be scaled by 10000 and clipped.

    :param clip_range:
        A tuple containing the range to clip the data limits to.
        Default is (1, 10000).

    :return:
    """
    expr = "data * 10000 + p5"
    p5 = numpy.float32(0.5)  # noqa # pylint: disable
    data_mask = ~numpy.isfinite(data)

    # evaluate
    result = numexpr.evaluate(expr).astype('int16')

    # clip the data range, if desired
    if clip:
        minv, maxv = clip_range
        numpy.clip(result, minv, maxv, out=result)

    # re-insert the null data pixels
    result[data_mask] = NO_DATA_VALUE

    return result


def lambertian_block(radiance, a, b, s):
    """
    A convenience function for evaluating lambertian reflectance
    over a single block of data.

    :param radiance:
        Top of Atmosphere Radiance.

    :param a:
        An atmpsheric coefficient derived from the radiative transfer.
        Calculated as:
            (DIR + DIF) / pi * (TV + TDV)
        Where:
            DIR = Direct irradiance at the surface
            DIF = Diffuse irradiance at the surface
            TV = Dirrect transmittance in the view direction
            TDV = Diffuse transmittance in the view direction

    :param b:
        Path radiance due to atmospheric scattering.

    :param s:
        Atmospheric albedo.

    :return:
        An numpy.ndarray.
    """
    expr = "(radiance - b) / (a + s * (radiance - b))"
    result = numexpr.evaluate(expr)

    return result


def lambertian_tiled(acquisition, a, b, s, esun=None, outds=None):
    """
    Calculate lambertian reflectance coupled with a correction
    for atmospheric adjacency (if a psf kernel is supplied).

    :param a:
        An atmpsheric coefficient derived from the radiative transfer.
        Calculated as:
            (DIR + DIF) / pi * (TV + TDV)
        Where:
            DIR = Direct irradiance at the surface
            DIF = Diffuse irradiance at the surface
            TV = Dirrect transmittance in the view direction
            TDV = Diffuse transmittance in the view direction

    :param b:
        Path radiance due to atmospheric scattering.

    :param s:
        Atmospheric albedo.

    :param esun:
        Exoatmospheric solar irradiance. Default is None.

    :param outds:
        A HDF5 Dataset to contain the float32 output result.

    :return:
        A float32 2D NumPy array of type float32 with NaN's populating
        invalid elements.

    :notes:
        Internally processed as tiles to conserve memory as
        calculating TOARadiance will output float64.
    """
    dims = (acquisition.lines, acquisition.samples)
    data_mask = numpy.zeros(dims, dtype='bool')

    # the processing is biased towards receiving a HDF5 Dataset object
    # as such, if None received we'll be allocating additional memory for
    # the data subsets which then get copied onto the full dimensional array
    if outds is None:
        outds = numpy.zeros(dims, dtype='float32')

    # process by tile
    for tile in acquisition.tiles():
        # read the data corresponding to the current tile for all dataset
        # the original f90 routine specified single precision
        rad = acquisition.radiance_data(window=tile, out_no_data=NAN,
                                        esun=esun).astype('float32')

        # define the data mask (identify NaN pixels)
        data_mask = ~numpy.isfinite(rad)

        # lambertian
        lambt = lambertian_block(rad, a[tile], b[tile], s[tile])

        # account for original nulls and any evaluated nan's
        lambt[data_mask] = NAN

        outds[tile] = lambt

    if isinstance(outds, numpy.ndarray):
        return outds


def average_lambertian(acquisition, a, b, s, psf_kernel, esun=None,
                       normalise=True):
    """
    Calculate the lambertian reflectance by taking into account the
    the surrounding pixels.
    Essentially we convolve the point spread function as a kernel
    across the lambertian reflectance to get the surrounding
    contribution.
    Prior to convolution, the data is first smoothed by filling
    null values with an average. Currently this null filling process
    is evaluated by a run length average.
    After convolution and before returning, nulls are re-inserted back
    into the result.

    TODO:
        * Include an option to run convolution via fourier, to enable
          faster processing for large kernels.
        * Incorporate astropy for convolution to enable handling of
          null data as well as NaN's.

    :param acquisition:
        This process was designed to be applied to lambertian
        reflectance, but in practice it could be any 2D NumPy array.
        Currently, the process is assuming null data to be -999.

    :param psf_kernel:
        A 2D kernel/filter representing the point spread function of
        the atmosphere to resolve the atmospheric adjacency.

    :param normalise:
        A boolean indicating whether or not to normalise the
        psf_kernel. Default is True.

    :return:
        A 2D NumPy array of type float32, with NaN's populating
        invalid elements.
    """
    # normalise the kernel or not
    if normalise:
        psf_kernel = psf_kernel / psf_kernel.sum()

    data = lambertian_tiled(acquisition, a, b, s, esun)
    data_mask = ~numpy.isfinite(data)  # locate NaN's

    # can we correctly apply row-length averages?
    _, start_idx, end_idx = _sequential_valid_rows(data_mask)  # ignore all_valid for time being

    # fill nulls with run-length averages
    # _fill_nulls(data[start_idx:end_idx], null_mask[start_idx:end_idx])
    _fill_nulls(data, data_mask)

    # apply convolution
    result = numpy.full(data.shape, fill_value=numpy.nan, dtype='float32')
    ndimage.convolve(data[start_idx:end_idx], psf_kernel,
                     output=result[start_idx:end_idx])

    # insert nulls back into the array
    result[data_mask] = NAN

    return result


def adjacency_correction(lambertian, lambertian_average, fv):
    """
    Calculate lambertian reflectance with adjacency correction.

    :param lambertian:
        The lambertian reflectance.

    :param lambertian_average:
        Average lambertian reflectance of the surrounding pixels.
        Lambertian reflectance convolved with the point spread
        function describing the atmospheric adjacency effect.
        which is computed via the average_lambertian funtion.

    :param fv:
        Direct fraction of radiation in the view direction.

    :return:
        A numpy.ndarray of lambertian reflectance with atmospheric
        adjacency correction.
    """
    expr = "lambertian + ((1 - fv) / fv) * (lambertian - lambertian_average)"
    result = numexpr.evaluate(expr)

    return result


def lambertian_adjacency(acquisition, a, b, s, fv, psf_kernel, esun=None,
                         normalise=True, outds=None):
    """
    Workflow to calculate lambertian reflectance with atmospheric
    adjacency correction.

    :param acquisition:
        An instance of an acquisition object.

    :param a:
        An atmpsheric coefficient derived from the radiative transfer.
        Calculated as:
            (DIR + DIF) / pi * (TV + TDV)
        Where:
            DIR = Direct irradiance at the surface
            DIF = Diffuse irradiance at the surface
            TV = Dirrect transmittance in the view direction
            TDV = Diffuse transmittance in the view direction

    :param b:
        Path radiance due to atmospheric scattering.

    :param s:
        Atmospheric albedo.

    :param fv:
        Direct fraction of radiation in the view direction.

    :param psf_kernel:
        A 2D kernel/filter representing the point spread function of
        the atmosphere to resolve the atmospheric adjacency.

    :param esun:
        Exoatmospheric solar irradiance. Default is None.

    :param normalise:
        A boolean indicating whether or not to normalise the
        psf_kernel. Default is True.

    :param outds:
        A HDF5 Dataset to contain the float32 output result.

    :return:
        A float32 2D NumPy array of type float32 with NaN's populating
        invalid elements.
    """
    dims = (acquisition.lines, acquisition.samples)

    if outds is None:
        outds = numpy.zeros(dims, dtype='float32')

    avg_lambt = average_lambertian(acquisition, a, b, s, psf_kernel, esun,
                                   normalise)

    # process by tile
    for tile in acquisition.tiles():
        # read the data corresponding to the current tile for all dataset
        # the original f90 routine specified single precision
        rad = acquisition.radiance_data(window=tile, out_no_data=NAN,
                                        esun=esun).astype('float32')

        # define the data mask (identify NaN pixels)
        data_mask = ~numpy.isfinite(rad)

        # lambertian
        lambt = lambertian_block(rad, a[tile], b[tile], s[tile])

        # correct for atmospheric adjacency
        adj_cor = adjacency_correction(lambt, avg_lambt[tile], fv[tile])

        # account for original nulls and any evaluated nan's
        adj_cor[data_mask] = NAN

        # copy to final destination (numpy.ndarray, or h5py.Dataset)
        outds[tile] = adj_cor

    if isinstance(outds, numpy.ndarray):
        return outds


def sky_glint(satellite_view, refractive_index=1.34, outds=None, tiles=None):
    """
    Calculate sky glint; Fresnel reflectance of a flat water body.
    Based on Fresnel reflectance at the sensor zenith view angle.
    Calculations are done in single precision (same as the original
    FORTRAN code that this was derived from).

    :param satellite_view:
        The satellite zenith view angle in degrees.

    :param refractive_index:
        The refractive index of water. Default is 1.34 and internally
        recast as a float32.

    :param outds:
        A HDF5 Dataset to contain the float32 output result.

    :tiles:
        Optionally a list of tile or windows of data subsets on which to
        iterate over. Default is None, in which case the entire array is
        worked on at once. ((ystart, yend), (xstart, xend))

    :return:
        A 2D NumPy array of type float32.
    """
    if outds is None:
        outds = numpy.zeros(satellite_view.shape, dtype='float32')

    # if no tiles defined, then one big tile covering the entire array
    if tiles is None:
        tiles = [(slice(None, None, None), slice(None, None, None))]

    # force constants to float32 (reduce memory for array computations)
    rw = numpy.float32(refractive_index)  # noqa # pylint: disable
    p5 = numpy.float32(0.5)  # noqa # pylint: disable

    # value used for pixels that are flagged by the tolerence test
    value = numpy.abs((rw - 1) / (rw + 1))**2

    for tile in tiles:
        theta = numpy.deg2rad(satellite_view[tile], dtype='float32')  # noqa # pylint: disable

        expr = "arcsin(sin(theta) / rw)"
        theta_prime = numexpr.evaluate(expr)  # noqa # pylint: disable

        # tolerance mask; suitable for very low angles, and those at nadir
        expr = "abs(theta + theta_prime) < 1.0e-5"
        tolerance_mask = numexpr.evaluate(expr)

        expr = ("p5 * ((sin(theta-theta_prime) / sin(theta+theta_prime))**2 "
                "+ (tan(theta-theta_prime) / tan(theta+theta_prime))**2)")

        # sky glint (taken as fresnel reflectance at theta)
        # this part is now less optimal if dealing with no tiles
        # as it is now the full array
        sky_g = numexpr.evaluate(expr)

        # insert value for any pixels that are flagged by the tolerance test
        sky_g[tolerance_mask] = value

        outds[tile] = sky_g

    if isinstance(outds, numpy.ndarray):
        return outds


def scattering(average_lambertian, s):
    """
    Needs documentation;
    This function is not detailed in the paper, but is defined in the
    sample code that was provided.

    :param average_lambertian:
        Average lambertian reflectance of the surrounding pixels.
        Lambertian reflectance convolved with the point spread
        function describing the atmospheric adjacency effect.
        which is computed via the average_lambertian funtion.

    :param s:
        Atmospheric albedo.

    :return:
        A 2D NumPy array of type float32.
    """
    expr = "s * average_lambertian / (1 - s * average_lambertian)"
    result = numexpr.evaluate(expr)

    return result


def sky_glint_correction(lambertian, fs, s, lambertian_average, satellite_view,
                         tiles=None, refractive_index=1.34):
    """
    Apply sky glint correction.

    :param lambertian:
        Lambertian reflectance with atmospheric adjacency correction.

    :param fs:
        Direct fraction in the sun direction.

    :param s:
        Atmospheric albedo.

    :param average_lambertian:
        Average lambertian reflectance of the surrounding pixels.
        Lambertian reflectance convolved with the point spread
        function describing the atmospheric adjacency effect.
        which is computed via the average_lambertian funtion.

    :param satellite_view:
        The satellite zenith view angle in degrees.

    :tiles:
        Optionally a list of tile or windows of data subsets on which to
        iterate over. Default is None, in which case the entire array is
        worked on at once. ((ystart, yend), (xstart, xend))

    :param refractive_index:
        The refractive index of water. Default is 1.34 and internally
        recast as a float32.

    :return:
        A 2D NumPy array of type float32.
    """
    result = numpy.zeros(lambertian.shape, dtype='float32')

    # if no tiles defined, then one big tile covering the entire array
    if tiles is None:
        tiles = [(slice(None, None, None), slice(None, None, None))]

    for tile in tiles:
        lambt = lambertian[tile]
        data_mask = ~numpy.isfinite(lambt)
        scat = scattering(lambertian_average[tile], s[tile])  # noqa # pylint: disable

        expr = "lambt - ((1 - fs) + scat) * sky_g"
        sky_g = sky_glint(satellite_view[tile], refractive_index)  # noqa # pylint: disable

        # evaluate; insert nulls back into array
        numexpr.evaluate(expr, out=result[tile])
        result[tile][data_mask] = NAN

    return result


def _calculate_reflectance(acquisition, acquisitions, interpolation_fname,
                           satellite_solar_angles_fname, slope_aspect_fname,
                           relative_slope_fname, incident_angles_fname,
                           exiting_angles_fname, shadow_masks_fname,
                           ancillary_fname, rori, out_fname, compression,
                           filter_opts, normalized_solar_zenith):
    """
    A private wrapper for dealing with the internal custom workings of the
    NBAR workflow.
    """
    with h5py.File(interpolation_fname, 'r') as fid_interp,\
        h5py.File(satellite_solar_angles_fname, 'r') as fid_sat_sol,\
        h5py.File(slope_aspect_fname, 'r') as fid_slp_asp,\
        h5py.File(relative_slope_fname, 'r') as fid_rel_slp,\
        h5py.File(incident_angles_fname, 'r') as fid_inc,\
        h5py.File(exiting_angles_fname, 'r') as fid_exi,\
        h5py.File(shadow_masks_fname, 'r') as fid_shadow,\
        h5py.File(ancillary_fname, 'r') as fid_anc,\
        h5py.File(out_fname, 'w') as fid:

        grp1 = fid_interp[GroupName.INTERP_GROUP.value]
        grp2 = fid_sat_sol[GroupName.SAT_SOL_GROUP.value]
        grp3 = fid_slp_asp[GroupName.SLP_ASP_GROUP.value]
        grp4 = fid_rel_slp[GroupName.REL_SLP_GROUP.value]
        grp5 = fid_inc[GroupName.INCIDENT_GROUP.value]
        grp6 = fid_exi[GroupName.EXITING_GROUP.value]
        grp7 = fid_shadow[GroupName.SHADOW_GROUP.value]
        grp8 = fid_anc[GroupName.ANCILLARY_GROUP.value]
        calculate_reflectance(acquisition, grp1, grp2, grp3, grp4, grp5, grp6,
                              grp7, grp8, rori, fid, compression, filter_opts,
                              normalized_solar_zenith)

        create_ard_yaml(acquisitions, grp8, fid, normalized_solar_zenith)


def calculate_reflectance(acquisition, interpolation_group,
                          satellite_solar_group, slope_aspect_group,
                          relative_slope_group, incident_angles_group,
                          exiting_angles_group, shadow_masks_group,
                          ancillary_group, rori, out_group=None,
                          compression=H5CompressionFilter.LZF,
                          filter_opts=None, normalized_solar_zenith=45.0,
                          esun=None, psf_kernel=None):
    """
    Calculates Lambertian, BRDF corrected and BRDF + terrain
    illumination corrected surface reflectance.

    :param acquisition:
        An instance of an acquisition object.

    :param interpolation_group:
        The root HDF5 `Group` that contains the interpolated
        atmospheric coefficients.
        The dataset pathnames are given by:

        * DatasetName.INTERPOLATION_FMT

    :param satellite_solar_group:
        The root HDF5 `Group` that contains the solar zenith and
        solar azimuth datasets specified by the pathnames given by:

        * DatasetName.SOLAR_ZENITH
        * DatasetName.SOLAR_AZIMUTH
        * DatasetName.SATELLITE_VIEW
        * DatasetName.SATELLITE_AZIMUTH
        * DatasetName.RELATIVE_AZIMUTH
        
    :param slope_aspect_group:
        The root HDF5 `Group` that contains the slope and aspect
        datasets specified by the pathnames given by:

        * DatasetName.SLOPE
        * DatasetName.ASPECT

    :param relative_slope_group:
        The root HDF5 `Group` that contains the relative slope dataset
        specified by the pathname given by:

        * DatasetName.RELATIVE_SLOPE

    :param incident_angles_group:
        The root HDF5 `Group` that contains the incident
        angle dataset specified by the pathname given by:

        * DatasetName.INCIDENT

    :param exiting_angles_group:
        The root HDF5 `Group` that contains the exiting
        angle dataset specified by the pathname given by:

        * DatasetName.EXITING

    :param shadow_masks_group:
        The root HDF5 `Group` that contains the combined shadow
        masks; self shadow, cast shadow (solar),
        cast shadow (satellite), dataset specified by the pathname
        given by:

        * DatasetName.COMBINED_SHADOW

    :param ancillary_group:
        The root HDF5 `Group` that contains the Isotropic (iso),
        RossThick (vol), and LiSparseR (geo) BRDF scalar parameters.
        The dataset pathnames are given by:

        * DatasetName.BRDF_FMT

    :param rori:
        Threshold for terrain correction. Fuqin to document.

    :param out_group:
        If set to None (default) then the results will be returned
        as an in-memory hdf5 file, i.e. the `core` driver. Otherwise,
        a writeable HDF5 `Group` object.

        The dataset names will be given by the format string detailed
        by:

        * DatasetName.REFLECTANCE_FMT

        The reflectance products are:

        * lambertian
        * nbar (BRDF corrected reflectance)
        * nbart (BRDF + terrain illumination corrected reflectance)

    :param compression:
        The compression filter to use.
        Default is H5CompressionFilter.LZF 

    :param filter_opts:
        A dict of key value pairs available to the given configuration
        instance of H5CompressionFilter. For example
        H5CompressionFilter.LZF has the keywords *chunks* and *shuffle*
        available.
        Default is None, which will use the default settings for the
        chosen H5CompressionFilter instance.

    :param normalized_solar_zenith:
        A float value type to normalize reflectance to a particular angle.

    :param esun
        A float value type. A solar solar irradiance normal to atmosphere
        in unit of W/sq cm/sr/nm.

    :param psf_kernel:
        A 2D kernel/filter representing the point spread function of
        the atmosphere to resolve the atmospheric adjacency.
        If not supplied then atmospheric adjacency and sky glint
        contributions will not be calculated.

    :return:
        An opened `h5py.File` object, that is either in-memory using the
        `core` driver, or on disk.
    """
    geobox = acquisition.gridded_geo_box()
    bn = acquisition.band_name

    # NOTES:
    #    Lots of unused variables, but that is ok for the time being.
    #    I've closed off the section where they'll be used. But for now,
    #    we only need to get functionality enabled for processing
    #    atmospheric correction over water
    dname_fmt = DatasetName.INTERPOLATION_FMT.value
    fv_dataset = interpolation_group[dname_fmt.format(coefficient=AC.FV.value,
                                                      band_name=bn)]
    fs_dataset = interpolation_group[dname_fmt.format(coefficient=AC.FS.value,
                                                      band_name=bn)]
    b_dataset = interpolation_group[dname_fmt.format(coefficient=AC.B.value,
                                                     band_name=bn)]
    s_dataset = interpolation_group[dname_fmt.format(coefficient=AC.S.value,
                                                     band_name=bn)]
    a_dataset = interpolation_group[dname_fmt.format(coefficient=AC.A.value,
                                                     band_name=bn)]
    dir_dataset = interpolation_group[dname_fmt.format(coefficient=AC.DIR.value,
                                                       band_name=bn)]
    dif_dataset = interpolation_group[dname_fmt.format(coefficient=AC.DIF.value,
                                                       band_name=bn)]
    ts_dataset = interpolation_group[dname_fmt.format(coefficient=AC.TS.value,
                                                      band_name=bn)]
    solar_zenith_dset = satellite_solar_group[DatasetName.SOLAR_ZENITH.value]
    solar_azimuth_dset = satellite_solar_group[DatasetName.SOLAR_AZIMUTH.value]
    satellite_v_dset = satellite_solar_group[DatasetName.SATELLITE_VIEW.value]
    relative_a_dset = satellite_solar_group[DatasetName.RELATIVE_AZIMUTH.value]
    slope_dataset = slope_aspect_group[DatasetName.SLOPE.value]
    aspect_dataset = slope_aspect_group[DatasetName.ASPECT.value]
    relative_s_dset = relative_slope_group[DatasetName.RELATIVE_SLOPE.value]
    incident_angle_dataset = incident_angles_group[DatasetName.INCIDENT.value]
    exiting_angle_dataset = exiting_angles_group[DatasetName.EXITING.value]
    shadow_dataset = shadow_masks_group[DatasetName.COMBINED_SHADOW.value]

    dname_fmt = DatasetName.BRDF_FMT.value
    dname = dname_fmt.format(band_name=bn, parameter=BrdfDirectionalParameters.ALPHA_1.value)
    brdf_alpha1 = ancillary_group[dname][()]

    dname = dname_fmt.format(band_name=bn, parameter=BrdfDirectionalParameters.ALPHA_2.value)
    brdf_alpha2 = ancillary_group[dname][()]

    # Initialise the output file
    if out_group is None:
        fid = h5py.File('surface-reflectance.h5', driver='core',
                        backing_store=False)
    else:
        fid = out_group

    if GroupName.STANDARD_GROUP.value not in fid:
        fid.create_group(GroupName.STANDARD_GROUP.value)

    if filter_opts is None:
        filter_opts = {}
    else:
        filter_opts = filter_opts.copy()
    filter_opts['chunks'] = acquisition.tile_size

    kwargs = compression.config(**filter_opts).dataset_compression_kwargs()
    grp = fid[GroupName.STANDARD_GROUP.value]
    kwargs['shape'] = (acquisition.lines, acquisition.samples)
    kwargs['fillvalue'] = NO_DATA_VALUE
    kwargs['dtype'] = 'int16'

    # create the integer datasets
    # lambertian
    dname_fmt = DatasetName.REFLECTANCE_FMT.value
    dname = dname_fmt.format(product=AP.LAMBERTIAN.value, band_name=bn)
    lmbrt_dset = grp.create_dataset(dname, **kwargs)

    # nbar
    dname = dname_fmt.format(product=AP.NBAR.value, band_name=bn)
    nbar_dset = grp.create_dataset(dname, **kwargs)

    # nbart
    dname = dname_fmt.format(product=AP.NBART.value, band_name=bn)
    nbart_dset = grp.create_dataset(dname, **kwargs)

    # attach some attributes to the image datasets
    attrs = {'crs_wkt': geobox.crs.ExportToWkt(),
             'geotransform': geobox.transform.to_gdal(),
             'no_data_value': kwargs['fillvalue'],
             'rori_threshold_setting': rori,
             'platform_id': acquisition.platform_id,
             'sensor_id': acquisition.sensor_id,
             'band_id': acquisition.band_id,
             'band_name': bn,
             'alias': acquisition.alias}

    # lambertian with no additional corrections applied
    desc = "Contains the lambertian reflectance data scaled by 10000."
    attrs['description'] = desc
    attach_image_attributes(lmbrt_dset, attrs)

    # nbar
    desc = "Contains the brdf corrected reflectance data scaled by 10000."
    attrs['description'] = desc
    attach_image_attributes(nbar_dset, attrs)

    # nbart
    desc = ("Contains the brdf and terrain corrected reflectance data scaled "
            "by 10000.")
    attrs['description'] = desc
    attach_image_attributes(nbart_dset, attrs)

    # lambertian with atmospheric adjacency correction
    # TODO:
    #    more utilisation of HDF5 datasets to cut down the memory usage
    #    i.e. parse through the output H5 dataset to store the lambertian
    #    adjacency correction result.
    #    sky_glint can be produced using a tiling routine, and output direct
    #    to a HDF5 dataset.
    if psf_kernel is not None:
        # NOTES:
        #    for the time being we're outputting float32 while water atcor details
        #    are figured out on how sun and sky glint is to be applied.
        #    it is unknown that if we truncate to int16 now, will this have an
        #    impact applying glint correction later
        kwargs['fillvalue'] = numpy.nan
        kwargs['dtype'] = 'float32'

        # update the attrs for the float32 datasets
        attrs['no_data_value'] = numpy.nan

        # *** lambertian with adjacency correction; may not be required later ***
        # TODO:
        #    confirm what is required and how to be delivered once production
        #    details have been determined.
        dname = dname_fmt.format(product=AP.ADJ.value, band_name=bn)
        adj_dset = grp.create_dataset(dname, **kwargs)

        # *** sky glint dataset; may not be required later ***
        # TODO:
        #    confirm what is required and how to be delivered once production
        #    details have been determined.
        dname = dname_fmt.format(product=AP.SKY.value, band_name=bn)
        skyg_dset = grp.create_dataset(dname, **kwargs)

        # TODO:
        #    change the description once final product details are defined
        #    currently output unscaled float32 data
        desc = "Contains the lambertian reflectace corrected for atmospheric adjacency."
        attrs['description'] = desc
        attach_image_attributes(adj_dset, attrs)

        # *** sky glint dataset; may not be required later ***
        # TODO:
        #    confirm what is required and how to be delivered once production
        #    details have been determined.
        desc = "Contains the sky glint coefficient."
        attrs['description'] = desc
        attach_image_attributes(skyg_dset, attrs)

        # calculate
        lambertian_adjacency(acquisition, a_dataset, b_dataset, s_dataset,
                             fv_dataset, psf_kernel, esun, adj_dset)
        sky_glint(satellite_v_dset, outds=skyg_dset, tiles=acquisition.tiles())

    # NOTES:
    #    for the time being, output the atmospheric adjacency corrected
    #    lambertian reflectance alongside the sky glint coefficients.
    #    we won't apply sky glint correction, nor sun glint correction
    #    in order for users to apply sun + sky glint correction, they'll also
    #    need the fs dataset (Direct fraction in the sun direction).

    # process by tile
    for tile in acquisition.tiles():
        # define some static arguments
        acq_args = {'window': tile,
                    'out_no_data': NAN,
                    'esun': esun}
        f32_args = {'dtype': numpy.float32, 'transpose': True}

        # nested function call ... not the most pretty way of doing something
        # TODO
        #    create a func that takes acq, window, a, b, s and calculates
        #    lambertian; we need float32 and will save as int16
        ref_lm = lambertian_block(
            acquisition.radiance_data(**acq_args).astype('float32'),
            a_dataset[tile],
            b_dataset[tile],
            s_dataset[tile])

        # NOTES
        #    for this prototype, use lambertian with adjacency correction
        #    as input into nbar and nbart correction

        shadow = as_array(shadow_dataset[tile], numpy.int8, transpose=True)
        solar_zenith = as_array(solar_zenith_dset[tile], **f32_args)
        solar_azimuth = as_array(solar_azimuth_dset[tile], **f32_args)
        satellite_view = as_array(satellite_v_dset[tile], **f32_args)
        relative_angle = as_array(relative_a_dset[tile], **f32_args)
        slope = as_array(slope_dataset[tile], **f32_args)
        aspect = as_array(aspect_dataset[tile], **f32_args)
        incident_angle = as_array(incident_angle_dataset[tile], **f32_args)
        exiting_angle = as_array(exiting_angle_dataset[tile], **f32_args)
        relative_slope = as_array(relative_s_dset[tile], **f32_args)
        s_mod = as_array(s_dataset[tile], **f32_args)
        fs = as_array(fs_dataset[tile], **f32_args)
        fv = as_array(fv_dataset[tile], **f32_args)
        ts = as_array(ts_dataset[tile], **f32_args)
        direct = as_array(dir_dataset[tile], **f32_args)
        diffuse = as_array(dif_dataset[tile], **f32_args)

        # TODO
        #    use a more elegant way so we don't have to do an
        #    *if* check on every loop
        if psf_kernel is not None:
            input_lambertian = adj_dset[tile]
        else:
            input_lambertian = ref_lm

        # Allocate the output arrays
        ysize, xsize = ref_lm.shape
        ref_brdf = numpy.full((ysize, xsize), numpy.nan, dtype='float32')
        ref_terrain = numpy.full((ysize, xsize), numpy.nan, dtype='float32')

        # Run terrain correction
        reflectance(xsize, ysize, rori, brdf_alpha1, brdf_alpha2,
                    acquisition.reflectance_adjustment, numpy.float32(NAN),
                    normalized_solar_zenith, shadow, solar_zenith,
                    solar_azimuth, satellite_view, relative_angle,
                    slope, aspect, incident_angle, exiting_angle,
                    relative_slope, s_mod, fs, fv, ts, direct, diffuse,
                    input_lambertian.transpose(),
                    ref_brdf.transpose(), ref_terrain.transpose())

        # Write the current tile to disk
        # TODO
        #    use h5_dataset.write_direct()
        lmbrt_dset[tile] = scale_reflectance(ref_lm)
        nbar_dset[tile] = scale_reflectance(ref_brdf)
        nbart_dset[tile] = scale_reflectance(ref_terrain)

    # close any still opened files, arrays etc associated with the acquisition
    acquisition.close()

    if out_group is None:
        return fid


def link_standard_data(input_fnames, out_fname):
    # TODO: incorporate linking for multi-granule and multi-group
    #       datasets
    """
    Links the individual reflectance and surface temperature
    results into a single file for easier access.
    """
    for fname in input_fnames:
        with h5py.File(fname, 'r') as fid:
            dataset_names = find(fid, dataset_class='IMAGE')

        for dname in dataset_names:
            create_external_link(fname, dname, out_fname, dname)

        # metadata
        with h5py.File(fname, 'r') as fid:
            with h5py.File(out_fname) as out_fid:
                yaml_dname = DatasetName.NBAR_YAML.value
                if yaml_dname in fid and yaml_dname not in out_fid:
                    fid.copy(yaml_dname, out_fid, name=yaml_dname)

                yaml_dname = DatasetName.SBT_YAML.value
                if yaml_dname in fid and yaml_dname not in out_fid:
                    fid.copy(yaml_dname, out_fid, name=yaml_dname)
