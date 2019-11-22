#!/usr/bin/env python

"""
Calculates the Lambertian, BRDF corrected and BRDF + Terrain corrected
----------------------------------------------------------------------

reflectance
-----------
"""

from __future__ import absolute_import, print_function
from os.path import join as pjoin
import tempfile
import numpy
import numexpr
import h5py

from wagl.constants import DatasetName, GroupName, BrdfDirectionalParameters
from wagl.constants import AtmosphericCoefficients as AC
from wagl.constants import ArdProducts as AP
from wagl.convolution import convolve
from wagl.data import as_array
from wagl.hdf5 import H5CompressionFilter, attach_image_attributes
from wagl.hdf5 import create_external_link, find
from wagl.metadata import create_ard_yaml
from wagl.__surface_reflectance import reflectance

NO_DATA_VALUE = -999
NAN = numpy.nan
PI = numpy.pi


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
    result = numexpr.evaluate(expr).astype("int16")

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

    # the processing is biased towards receiving a HDF5 Dataset object
    # as such, if None received we'll be allocating additional memory for
    # the data subsets which then get copied onto the full dimensional array
    if outds is None:
        outds = numpy.zeros(dims, dtype="float32")

    # process by tile
    for tile in acquisition.tiles():
        # read the data corresponding to the current tile for all dataset
        # the original f90 routine specified single precision
        rad = acquisition.radiance_data(window=tile, out_no_data=NAN, esun=esun).astype(
            "float32"
        )

        # define the data mask (identify NaN pixels)
        data_mask = ~numpy.isfinite(rad)

        # lambertian
        lambt = lambertian_block(rad, a[tile], b[tile], s[tile])

        # account for original nulls and any evaluated nan's
        lambt[data_mask] = NAN

        outds[tile] = lambt

    if isinstance(outds, numpy.ndarray):
        return outds


def average_lambertian(lambertian, psf_kernel, outds=None):
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

    :return:
        A 2D NumPy array of type float32, with NaN's populating
        invalid elements.
    """

    # normalise the kernel
    psf_kernel = psf_kernel / psf_kernel.sum()

    data_mask = ~numpy.isfinite(lambertian)  # locate NaN's

    # apply convolution
    result = convolve(lambertian, psf_kernel, data_mask, True)

    # insert nulls back into the array
    result[data_mask] = NAN

    if outds is not None:
        outds.write_direct(result)
        return

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


def sky_glint(satellite_view, refractive_index=1.34):
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

    :return:
        A 2D NumPy array of type float32.
    """
    # force constants to float32 (reduce memory for array computations)
    rw = numpy.float32(refractive_index)  # noqa # pylint: disable
    p5 = numpy.float32(0.5)  # noqa # pylint: disable

    # value used for pixels that are flagged by the tolerence test
    value = numpy.abs((rw - 1) / (rw + 1)) ** 2

    theta = numpy.deg2rad(satellite_view, dtype="float32")  # noqa # pylint: disable

    expr = "arcsin(sin(theta) / rw)"
    theta_prime = numexpr.evaluate(expr)  # noqa # pylint: disable

    # tolerance mask; suitable for very low angles, and those at nadir
    expr = "abs(theta + theta_prime) < 1.0e-5"
    tolerance_mask = numexpr.evaluate(expr)

    expr = (
        "p5 * ((sin(theta-theta_prime) / sin(theta+theta_prime))**2 "
        "+ (tan(theta-theta_prime) / tan(theta+theta_prime))**2)"
    )

    # sky glint (taken as fresnel reflectance at theta)
    # this part is now less optimal if dealing with no tiles
    # as it is now the full array
    sky_g = numexpr.evaluate(expr)

    # insert value for any pixels that are flagged by the tolerance test
    sky_g[tolerance_mask] = value

    return sky_g


def sun_glint(
    satellite_view, solar_zenith, relative_azimuth, wind_speed, refractive_index=1.34
):
    """
    Calculate sun glint based on the Cox and Munk (1954) model.

    Calculations are done in single precision (same as the original
    FORTRAN code that this was derived from).

    :param satellite_view:
        The satellite zenith view angle in degrees.

    :param solar_zenith:
        The solar zenith angle in degrees.

    :param relative_azimuth:
        The relative azimuth angle between sun and view direction.

    :param wind_speed:
        The wind speed in m/s.

    :param refractive_index:
        The refractive index of water. Default is 1.34 and internally
        recast as a float32.

    :return:
        A 2D NumPy array of type float32.
    """

    expr = "relative_azimuth > 180.0"
    angle_mask = numexpr.evaluate(expr)
    relative_azimuth[angle_mask] = relative_azimuth[angle_mask] - 360.0

    # force constants to float32 (reduce memory for array computations)
    rw = numpy.float32(refractive_index)  # noqa # pylint: disable
    p5 = numpy.float32(0.5)  # noqa # pylint: disable

    theta_view = numpy.deg2rad(
        satellite_view, dtype="float32"
    )  # noqa # pylint: disable
    theta_sun = numpy.deg2rad(solar_zenith, dtype="float32")  # noqa # pylint: disable
    theta_phi = numpy.deg2rad(
        relative_azimuth, dtype="float32"
    )  # noqa # pylint: disable

    expr = "cos(theta_sun) * cos(theta_view) + sin(theta_view) * sin(theta_sun) * cos(theta_phi)"
    cos_psi = numexpr.evaluate(expr)  # noqa # pylint: disable

    expr = "sqrt((1.0 + cos_psi) / 2.0)"
    cos_omega = numexpr.evaluate(expr)

    expr = "(cos(theta_view) + cos(theta_sun)) / (2.0 * cos_omega)"
    cos_beta = numexpr.evaluate(expr)

    expr = "0.003 + 0.00512 * wind_speed"
    sigma2 = numexpr.evaluate(expr)

    expr = "((1.0 / cos_beta**2) - 1.0) / sigma2"
    fac = numexpr.evaluate(expr)

    expr = "exp(-fac) / (sigma2 * PI)"
    pval = numexpr.evaluate(expr)

    # tolerance mask
    expr = "fac > -log(1.0e-5)"
    tolerance_mask = numexpr.evaluate(expr)

    # insert 0.0 for any pixels that are flagged by the tolerance test
    pval[tolerance_mask] = 0.0

    expr = "pval * cos_omega / (4.0 * cos_beta**3)"
    gval = numexpr.evaluate(expr)

    # value used for pixels that are flagged by the tolerance test
    value = numpy.abs((rw - 1) / (rw + 1)) ** 2

    expr = "arccos(cos_omega)"
    omega = numexpr.evaluate(expr)

    expr = "arcsin(sin(omega) / rw)"
    omega_prime = numexpr.evaluate(expr)  # noqa # pylint: disable

    # tolerance mask
    expr = "abs(omega + omega_prime) < 1.0e-5"
    tolerance_mask = numexpr.evaluate(expr)

    expr = (
        "p5 * ((sin(omega-omega_prime) / sin(omega+omega_prime))**2 "
        "+ (tan(omega-omega_prime) / tan(omega+omega_prime))**2)"
    )

    pf_omega = numexpr.evaluate(expr)

    # insert value for any pixels that are flagged by the tolerance test
    pf_omega[tolerance_mask] = value

    expr = "pf_omega * gval * PI"
    frsun = numexpr.evaluate(expr)

    return frsun


def scattering(mean_lambertian, s):
    """
    Needs documentation;
    This function is not detailed in the paper, but is defined in the
    sample code that was provided.

    :param mean_lambertian:
        Average lambertian reflectance of the surrounding pixels.
        Lambertian reflectance convolved with the point spread
        function describing the atmospheric adjacency effect.
        which is computed via the average_lambertian funtion.

    :param s:
        Atmospheric albedo.

    :return:
        A 2D NumPy array of type float32.
    """
    expr = "s * mean_lambertian / (1 - s * mean_lambertian)"
    result = numexpr.evaluate(expr)

    return result


def sky_glint_correction(adjacency_corrected, fs, scattering, sky_glint):
    """
    Apply skyglint correction.

    :param adjacency_corrected:
        Lambertian reflectance corrected for atmospheric adjacency.

    :param fs:
        Direct fraction in the sun direction.

    :param scattering:
        TODO; Document. Was not detailed in the paper, but is defined
        in the sample code that was provided.

    :param sky_glint:
        Fresnel reflectance of a flat water body.

    :return:
        A 2D numpy.ndarray of type float32.
    """
    expr = "adjacency_corrected - ((1 - fs) + scattering) * sky_glint"
    sky_glint_corrected = numexpr.evaluate(expr)

    return sky_glint_corrected


def sun_glint_correction(sky_glint_corrected, fs, sun_glint):
    """
    Apply sun glint correction.

    :param sky_glint_corrected:
       surface reflectance corrected for sky glint.

    :param fs:
        Direct fraction in the sun direction.

    :param sun_glint:
        sun glint.

    :return:
        A 2D numpy.ndarray of type float32.
    """
    expr = "sky_glint_corrected - fs * sun_glint"
    sun_glint_corrected = numexpr.evaluate(expr)

    return sun_glint_corrected


def lambertian_corrections(
    lambertian,
    fs,
    fv,
    s,
    satellite_view,
    psf_kernel,
    refractive_index=1.34,
    tiles=None,
    out_adjacency=None,
    out_skyglint=None,
):
    """
    Workflow to apply atmospheric adjacency and skyglint correction
    to the lambertian reflectance.

    :param lambertian:
        A 2D NumPy array containing lambertian reflectance.

    :param fs:
        Direct fraction in the sun direction.

    :param fv:
        Direct fraction of radiation in the view direction.

    :param s:
        Atmospheric albedo.

    :param satellite_view:
        The satellite zenith view angle in degrees.

    :param psf_kernel:
        A 2D kernel/filter representing the point spread function of
        the atmosphere to resolve the atmospheric adjacency.

    :param refractive_index:
        The refractive index of water. Default is 1.34 and internally
        recast as a float32.

    :param tiles:
        Optionally a list of tile or windows of data subsets on which to
        iterate over. Default is None, in which case the entire array is
        worked on at once. ((ystart, yend), (xstart, xend))

    :param out_adjacency:
        A HDF5 Dataset to save the output of the atmospheric adjacency
        correction.

    :param out_skyglint:
        A HDF5 Dataset to save the output of the sky glint correction.

    :return:
        If either the out_adjacency param or the out_skyglint param are
        not provided, then a tuple of NumPy float32 arrays will be returned.
        (out_adjacency, out_skyglint).

    :notes:
        Temporary files will be created in order to minimise the memory
        footprint. They will be cleaned up automatically.
    """
    dims = lambertian.shape
    if isinstance(lambertian, h5py.Dataset):
        chunks = lambertian.chunks
    else:
        chunks = None

    # if no tiles defined, then one big tile covering the entire array
    if tiles is None:
        tiles = [(slice(None, None, None), slice(None, None, None))]

    if out_adjacency is None or out_skyglint is None:
        out_adjacency = numpy.full(dims, fill_value=numpy.nan, dtype="float32")
        out_skyglint = numpy.full(dims, fill_value=numpy.nan, dtype="float32")

    # create a temp workspace for datasets that we don't need to keep
    with tempfile.TemporaryDirectory(".tmp", "corrections-") as tmpd:
        with h5py.File(pjoin(tmpd, "lambertian-corrections"), "w") as fid:
            # temp file
            avg_ds = fid.create_dataset(
                "lambertian-average",
                shape=dims,
                compression="lzf",
                shuffle=True,
                dtype=lambertian.dtype,
                chunks=chunks,
            )

            # average lambertian of surrounding pixels via the psf
            average_lambertian(lambertian, psf_kernel, avg_ds)

            # process each block of data
            for tile in tiles:
                lambt = lambertian[tile]
                data_mask = ~numpy.isfinite(lambt)

                # atmospheric adjacency correction
                adj_cor = adjacency_correction(lambt, avg_ds[tile], fv[tile])
                adj_cor[data_mask] = NAN

                scat = scattering(avg_ds[tile], s[tile])  # noqa # pylint: disable

                sky_g = sky_glint(
                    satellite_view[tile], refractive_index
                )  # noqa # pylint: disable

                # sky glint correction
                skyg_c = sky_glint_correction(adj_cor, fs[tile], scat, sky_g)
                skyg_c[data_mask] = NAN

                # scale to skyg_c to int16 if written direct to h5 dataset
                # retain adj_cor still in float32 to use in brdf correction
                if not isinstance(out_adjacency, numpy.ndarray):
                    skyg_c = scale_reflectance(skyg_c, clip=False)

                out_adjacency[tile] = adj_cor
                out_skyglint[tile] = skyg_c

    if isinstance(out_adjacency, numpy.ndarray):
        return out_adjacency, out_skyglint


def _calculate_reflectance(
    acquisition,
    acquisitions,
    interpolation_fname,
    satellite_solar_angles_fname,
    slope_aspect_fname,
    relative_slope_fname,
    incident_angles_fname,
    exiting_angles_fname,
    shadow_masks_fname,
    ancillary_fname,
    rori,
    out_fname,
    compression,
    filter_opts,
    normalized_solar_zenith,
):
    """
    A private wrapper for dealing with the internal custom workings of the
    NBAR workflow.
    """
    with h5py.File(interpolation_fname, "r") as fid_interp, h5py.File(
        satellite_solar_angles_fname, "r"
    ) as fid_sat_sol, h5py.File(slope_aspect_fname, "r") as fid_slp_asp, h5py.File(
        relative_slope_fname, "r"
    ) as fid_rel_slp, h5py.File(
        incident_angles_fname, "r"
    ) as fid_inc, h5py.File(
        exiting_angles_fname, "r"
    ) as fid_exi, h5py.File(
        shadow_masks_fname, "r"
    ) as fid_shadow, h5py.File(
        ancillary_fname, "r"
    ) as fid_anc, h5py.File(
        out_fname, "w"
    ) as fid:

        grp1 = fid_interp[GroupName.INTERP_GROUP.value]
        grp2 = fid_sat_sol[GroupName.SAT_SOL_GROUP.value]
        grp3 = fid_slp_asp[GroupName.SLP_ASP_GROUP.value]
        grp4 = fid_rel_slp[GroupName.REL_SLP_GROUP.value]
        grp5 = fid_inc[GroupName.INCIDENT_GROUP.value]
        grp6 = fid_exi[GroupName.EXITING_GROUP.value]
        grp7 = fid_shadow[GroupName.SHADOW_GROUP.value]
        grp8 = fid_anc[GroupName.ANCILLARY_GROUP.value]
        calculate_reflectance(
            acquisition,
            grp1,
            grp2,
            grp3,
            grp4,
            grp5,
            grp6,
            grp7,
            grp8,
            rori,
            fid,
            compression,
            filter_opts,
            normalized_solar_zenith,
        )

        create_ard_yaml(acquisitions, grp8, fid, normalized_solar_zenith)


def calculate_reflectance(
    acquisition,
    interpolation_group,
    satellite_solar_group,
    slope_aspect_group,
    relative_slope_group,
    incident_angles_group,
    exiting_angles_group,
    shadow_masks_group,
    ancillary_group,
    rori,
    out_group=None,
    compression=H5CompressionFilter.LZF,
    filter_opts=None,
    normalized_solar_zenith=45.0,
    esun=None,
    psf_kernel=None,
):
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
    fv_dataset = interpolation_group[
        dname_fmt.format(coefficient=AC.FV.value, band_name=bn)
    ]
    fs_dataset = interpolation_group[
        dname_fmt.format(coefficient=AC.FS.value, band_name=bn)
    ]
    b_dataset = interpolation_group[
        dname_fmt.format(coefficient=AC.B.value, band_name=bn)
    ]
    s_dataset = interpolation_group[
        dname_fmt.format(coefficient=AC.S.value, band_name=bn)
    ]
    a_dataset = interpolation_group[
        dname_fmt.format(coefficient=AC.A.value, band_name=bn)
    ]
    dir_dataset = interpolation_group[
        dname_fmt.format(coefficient=AC.DIR.value, band_name=bn)
    ]
    dif_dataset = interpolation_group[
        dname_fmt.format(coefficient=AC.DIF.value, band_name=bn)
    ]
    ts_dataset = interpolation_group[
        dname_fmt.format(coefficient=AC.TS.value, band_name=bn)
    ]
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
    dname = dname_fmt.format(
        band_name=bn, parameter=BrdfDirectionalParameters.ALPHA_1.value
    )
    brdf_alpha1 = ancillary_group[dname][()]

    dname = dname_fmt.format(
        band_name=bn, parameter=BrdfDirectionalParameters.ALPHA_2.value
    )
    brdf_alpha2 = ancillary_group[dname][()]

    # Initialise the output file
    if out_group is None:
        fid = h5py.File("surface-reflectance.h5", driver="core", backing_store=False)
    else:
        fid = out_group

    if GroupName.STANDARD_GROUP.value not in fid:
        fid.create_group(GroupName.STANDARD_GROUP.value)

    if filter_opts is None:
        filter_opts = {}
    else:
        filter_opts = filter_opts.copy()
    filter_opts["chunks"] = acquisition.tile_size

    kwargs = compression.config(**filter_opts).dataset_compression_kwargs()
    grp = fid[GroupName.STANDARD_GROUP.value]
    kwargs["shape"] = (acquisition.lines, acquisition.samples)
    kwargs["fillvalue"] = NO_DATA_VALUE
    kwargs["dtype"] = "int16"

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
    attrs = {
        "crs_wkt": geobox.crs.ExportToWkt(),
        "geotransform": geobox.transform.to_gdal(),
        "no_data_value": kwargs["fillvalue"],
        "rori_threshold_setting": rori,
        "platform_id": acquisition.platform_id,
        "sensor_id": acquisition.sensor_id,
        "band_id": acquisition.band_id,
        "band_name": bn,
        "alias": acquisition.alias,
    }

    # lambertian with no additional corrections applied
    desc = "Contains the lambertian reflectance data scaled by 10000."
    attrs["description"] = desc
    attach_image_attributes(lmbrt_dset, attrs)

    # nbar
    desc = "Contains the brdf corrected reflectance data scaled by 10000."
    attrs["description"] = desc
    attach_image_attributes(nbar_dset, attrs)

    # nbart
    desc = (
        "Contains the brdf and terrain corrected reflectance data scaled " "by 10000."
    )
    attrs["description"] = desc
    attach_image_attributes(nbart_dset, attrs)

    # a HDF5 workspace for holding various temp datasets
    tmpdir = tempfile.TemporaryDirectory(suffix=".tmp", prefix="lambertian-")
    tmp_fid = h5py.File(pjoin(tmpdir.name, "lambertian-workspace.h5"), "w")

    kwargs["fillvalue"] = numpy.nan
    kwargs["dtype"] = "float32"

    # temporary file to hold float32 lambertian
    lamb_f32 = tmp_fid.create_dataset("lambertian", **kwargs)
    lambertian_tiled(acquisition, a_dataset, b_dataset, s_dataset, esun, lamb_f32)

    # lambertian with atmospheric adjacency correction
    if psf_kernel is not None:
        # temporary file to hold float32 lambertian adjacent corrected reflectance
        adj_dset_f32 = tmp_fid.create_dataset("lmbadj_corrected", **kwargs)

        # update values to reflect int16
        kwargs["fillvalue"] = NO_DATA_VALUE
        kwargs["dtype"] = "int16"

        # *** lambertian with adjacency correction; may not be required later ***
        # TODO:
        #    confirm what is required and how to be delivered once production
        #    details have been determined.
        dname = dname_fmt.format(product=AP.ADJ.value, band_name=bn)
        adj_dset = grp.create_dataset(dname, **kwargs)

        # *** sky glint correction dataset; may not be required later ***
        # TODO:
        #    confirm what is required and how to be delivered once production
        #    details have been determined.
        dname = dname_fmt.format(product=AP.SKY.value, band_name=bn)
        skygc_dset = grp.create_dataset(dname, **kwargs)

        # TODO:
        #    change the description once final product details are defined
        #    currently output unscaled float32 data

        desc = (
            "Contains the lambertian reflectace corrected for " "atmospheric adjacency."
        )
        attrs["description"] = desc
        attach_image_attributes(adj_dset, attrs)

        # *** sky glint dataset; may not be required later ***
        # TODO:
        #    confirm what is required and how to be delivered once production
        #    details have been determined.
        desc = "Contains the sky glint coefficient."
        attrs["description"] = desc

        attach_image_attributes(skygc_dset, attrs)

        # calculate
        lambertian_corrections(
            lamb_f32,
            fs_dataset,
            fv_dataset,
            s_dataset,
            satellite_v_dset,
            psf_kernel,
            1.34,
            acquisition.tiles(),
            adj_dset_f32,
            skygc_dset,
        )

    # NOTES:
    #    for the time being, output the atmospheric adjacency corrected
    #    lambertian reflectance alongside the sky glint coefficients.
    #    we won't apply sky glint correction, nor sun glint correction
    #    in order for users to apply sun + sky glint correction, they'll also
    #    need the fs dataset (Direct fraction in the sun direction).

    # process by tile
    for tile in acquisition.tiles():
        # define some static arguments

        f32_args = {"dtype": numpy.float32, "transpose": True}

        # load standard lambertian
        ref_lm = lamb_f32[tile]

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
            input_lambertian = adj_dset_f32[tile]
        else:
            input_lambertian = lamb_f32[tile]

        # Allocate the output arrays
        ysize, xsize = ref_lm.shape
        ref_brdf = numpy.full((ysize, xsize), numpy.nan, dtype="float32")
        ref_terrain = numpy.full((ysize, xsize), numpy.nan, dtype="float32")

        # Run terrain correction
        reflectance(
            xsize,
            ysize,
            rori,
            brdf_alpha1,
            brdf_alpha2,
            acquisition.reflectance_adjustment,
            numpy.float32(NAN),
            normalized_solar_zenith,
            shadow,
            solar_zenith,
            solar_azimuth,
            satellite_view,
            relative_angle,
            slope,
            aspect,
            incident_angle,
            exiting_angle,
            relative_slope,
            s_mod,
            fs,
            fv,
            ts,
            direct,
            diffuse,
            input_lambertian.transpose(),
            ref_brdf.transpose(),
            ref_terrain.transpose(),
        )

        # Write the current tile to disk
        # we will not clip data to for testing phase
        if psf_kernel is not None:
            adj_data = adj_dset_f32[tile]
            adj_dset.write_direct(
                scale_reflectance(adj_data, clip=False), dest_sel=tile
            )
        lmbrt_dset.write_direct(scale_reflectance(ref_lm, clip=False), dest_sel=tile)
        nbar_dset.write_direct(scale_reflectance(ref_brdf, clip=False), dest_sel=tile)
        nbart_dset.write_direct(
            scale_reflectance(ref_terrain, clip=False), dest_sel=tile
        )

    # close any still opened files, arrays etc associated with the acquisition
    acquisition.close()
    tmpdir.cleanup()
    tmp_fid.close()

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
        with h5py.File(fname, "r") as fid:
            dataset_names = find(fid, dataset_class="IMAGE")

        for dname in dataset_names:
            create_external_link(fname, dname, out_fname, dname)

        # metadata
        with h5py.File(fname, "r") as fid:
            with h5py.File(out_fname) as out_fid:
                yaml_dname = DatasetName.NBAR_YAML.value
                if yaml_dname in fid and yaml_dname not in out_fid:
                    fid.copy(yaml_dname, out_fid, name=yaml_dname)

                yaml_dname = DatasetName.SBT_YAML.value
                if yaml_dname in fid and yaml_dname not in out_fid:
                    fid.copy(yaml_dname, out_fid, name=yaml_dname)
