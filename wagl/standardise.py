#!/usr/bin/env python

from os.path import join as pjoin
import tempfile
import json

from posixpath import join as ppjoin
import h5py

from wagl.acquisition import acquisitions
from wagl.ancillary import collect_ancillary
from wagl.constants import ArdProducts as AP, GroupName, Workflow, BandType
from wagl.constants import ALBEDO_FMT, POINT_FMT, POINT_ALBEDO_FMT, Albedos
from wagl.constants import DatasetName
from wagl.constants import AtmosphericCoefficients
from wagl.dsm import get_dsm
from wagl.hdf5 import H5CompressionFilter, read_h5_table
from wagl.incident_exiting_angles import incident_angles, exiting_angles
from wagl.incident_exiting_angles import relative_azimuth_slope
from wagl.interpolation import interpolate
from wagl.longitude_latitude_arrays import create_lon_lat_grids
from wagl.metadata import create_ard_yaml
from wagl.modtran import format_json, prepare_modtran, run_modtran
from wagl.modtran import calculate_coefficients
from wagl.reflectance import calculate_reflectance
from wagl.satellite_solar_angles import calculate_angles
from wagl.terrain_shadow_masks import self_shadow, calculate_cast_shadow
from wagl.terrain_shadow_masks import combine_shadow_masks
from wagl.slope_aspect import slope_aspect_arrays
from wagl.temperature import surface_brightness_temperature
from wagl.pq import can_pq, run_pq
from wagl.modtran import JsonEncoder
from wagl.psf import run_modtran54
from wagl.logs import STATUS_LOGGER


# pylint disable=too-many-arguments
def card4l(level1, granule, workflow, vertices, method, pixel_quality, landsea,
           tle_path, aerosol, brdf, ozone_path,
           water_vapour, dem_path, dsm_fname, invariant_fname, modtran_exe, modtran54_exe,
           out_fname, ecmwf_path=None, rori=0.52, buffer_distance=8000,
           compression=H5CompressionFilter.LZF, filter_opts=None,
           h5_driver=None, acq_parser_hint=None, normalized_solar_zenith=45.):
    """
    CEOS Analysis Ready Data for Land.
    A workflow for producing standardised products that meet the
    CARD4L specification.

    :param level1:
        A string containing the full file pathname to the level1
        dataset.

    :param granule:
        A string containing the granule id to process.

    :param workflow:
        An enum from wagl.constants.Workflow representing which
        workflow workflow to run.

    :param vertices:
        An integer 2-tuple indicating the number of rows and columns
        of sample-locations ("coordinator") to produce.
        The vertex columns should be an odd number.

    :param method:
        An enum from wagl.constants.Method representing the
        interpolation method to use during the interpolation
        of the atmospheric coefficients.

    :param pixel_quality:
        A bool indicating whether or not to run pixel quality.

    :param landsea:
        A string containing the full file pathname to the directory
        containing the land/sea mask datasets.

    :param tle_path:
        A string containing the full file pathname to the directory
        containing the two line element datasets.

    :param aerosol:
        A string containing the full file pathname to the HDF5 file
        containing the aerosol data.

    :param brdf:
        A dict containing either user-supplied BRDF values, or the
        full file pathname to the directory containing the BRDF data
        and the decadal averaged BRDF data used for acquisitions
        prior to TERRA/AQUA satellite operations.

    :param ozone_path:
        A string containing the full file pathname to the directory
        containing the ozone datasets.

    :param water_vapour:
        A string containing the full file pathname to the directory
        containing the water vapour datasets.

    :param dem_path:
        A string containing the full file pathname to the directory
        containing the reduced resolution DEM.

    :param dsm_path:
        A string containing the full file pathname to the directory
        containing the Digital Surface Workflow for use in terrain
        illumination correction.

    :param invariant_fname:
        A string containing the full file pathname to the image file
        containing the invariant geo-potential data for use within
        the SBT process.

    :param modtran_exe:
        A string containing the full file pathname to the MODTRAN
        executable.

    :param modtran54_exe:
        A string containing the full file pathname to MODTRAN5.4
        executable.

    :param out_fname:
        A string containing the full file pathname that will contain
        the output data from the data standardisation process.
        executable.

    :param ecmwf_path:
        A string containing the full file pathname to the directory
        containing the data from the European Centre for Medium Weather
        Forcast, for use within the SBT process.

    :param rori:
        A floating point value for surface reflectance adjustment.
        TODO Fuqin to add additional documentation for this parameter.
        Default is 0.52.

    :param buffer_distance:
        A number representing the desired distance (in the same
        units as the acquisition) in which to calculate the extra
        number of pixels required to buffer an image.
        Default is 8000, which for an acquisition using metres would
        equate to 8000 metres.

    :param compression:
        An enum from hdf5.compression.H5CompressionFilter representing
        the desired compression filter to use for writing H5 IMAGE and
        TABLE class datasets to disk.
        Default is H5CompressionFilter.LZF.

    :param filter_opts:
        A dict containing any additional keyword arguments when
        generating the configuration for the given compression Filter.
        Default is None.

    :param h5_driver:
        The specific HDF5 file driver to use when creating the output
        HDF5 file.
        See http://docs.h5py.org/en/latest/high/file.html#file-drivers
        for more details.
        Default is None; which writes direct to disk using the
        appropriate driver for the underlying OS.

    :param acq_parser_hint:
        A string containing any hints to provide the acquisitions
        loader with.

    :param normalized_solar_zenith:
        Solar zenith angle to normalize for (in degrees). Default is 45 degrees.
    """
    json_fmt = pjoin(POINT_FMT, ALBEDO_FMT, ''.join([POINT_ALBEDO_FMT, '.json']))
    nvertices = vertices[0] * vertices[1]
    container = acquisitions(level1, hint=acq_parser_hint)

    # TODO: pass through an acquisitions container rather than pathname
    with h5py.File(out_fname, 'w', driver=h5_driver) as fid:
        fid.attrs['level1_uri'] = level1

        for grp_name in container.supported_groups:
            log = STATUS_LOGGER.bind(level1=container.label, granule=granule,
                                     granule_group=grp_name)

            # root group for a given granule and resolution group
            root = fid.create_group(ppjoin(granule, grp_name))
            acqs = container.get_acquisitions(granule=granule, group=grp_name)

            # include the resolution as a group attribute
            root.attrs['resolution'] = acqs[0].resolution

            # longitude and latitude
            log.info('Latitude-Longitude')
            create_lon_lat_grids(acqs[0], root, compression, filter_opts)

            # satellite and solar angles
            log.info('Satellite-Solar-Angles')
            calculate_angles(acqs[0], root[GroupName.LON_LAT_GROUP.value],
                             root, compression, filter_opts, tle_path)

            if workflow in (Workflow.STANDARD, Workflow.NBAR):

                # DEM
                log.info('DEM-retriveal')
                get_dsm(acqs[0], dsm_fname, buffer_distance, root, compression,
                        filter_opts)

                # slope & aspect
                log.info('Slope-Aspect')
                slope_aspect_arrays(acqs[0],
                                    root[GroupName.ELEVATION_GROUP.value],
                                    buffer_distance, root, compression,
                                    filter_opts)

                # incident angles
                log.info('Incident-Angles')
                incident_angles(root[GroupName.SAT_SOL_GROUP.value],
                                root[GroupName.SLP_ASP_GROUP.value],
                                root, compression, filter_opts)

                # exiting angles
                log.info('Exiting-Angles')
                exiting_angles(root[GroupName.SAT_SOL_GROUP.value],
                               root[GroupName.SLP_ASP_GROUP.value],
                               root, compression, filter_opts)

                # relative azimuth slope
                log.info('Relative-Azimuth-Angles')
                incident_group_name = GroupName.INCIDENT_GROUP.value
                exiting_group_name = GroupName.EXITING_GROUP.value
                relative_azimuth_slope(root[incident_group_name],
                                       root[exiting_group_name],
                                       root, compression, filter_opts)

                # self shadow
                log.info('Self-Shadow')
                self_shadow(root[incident_group_name],
                            root[exiting_group_name], root, compression,
                            filter_opts)

                # cast shadow solar source direction
                log.info('Cast-Shadow-Solar-Direction')
                dsm_group_name = GroupName.ELEVATION_GROUP.value
                calculate_cast_shadow(acqs[0], root[dsm_group_name],
                                      root[GroupName.SAT_SOL_GROUP.value],
                                      buffer_distance, root, compression,
                                      filter_opts)

                # cast shadow satellite source direction
                log.info('Cast-Shadow-Satellite-Direction')
                calculate_cast_shadow(acqs[0], root[dsm_group_name],
                                      root[GroupName.SAT_SOL_GROUP.value],
                                      buffer_distance, root, compression,
                                      filter_opts, False)

                # combined shadow masks
                log.info('Combined-Shadow')
                combine_shadow_masks(root[GroupName.SHADOW_GROUP.value],
                                     root[GroupName.SHADOW_GROUP.value],
                                     root[GroupName.SHADOW_GROUP.value],
                                     root, compression, filter_opts)

        # nbar and sbt ancillary
        log = STATUS_LOGGER.bind(level1=container.label, granule=granule,
                                 granule_group=None)

        # granule root group
        root = fid[granule]

        # get the highest resolution group containing supported bands
        acqs, grp_name = container.get_highest_resolution(granule=granule)

        grn_con = container.get_granule(granule=granule, container=True)
        res_group = root[grp_name]

        log.info('Ancillary-Retrieval')
        nbar_paths = {'aerosol_dict': aerosol,
                      'water_vapour_dict': water_vapour,
                      'ozone_path': ozone_path,
                      'dem_path': dem_path,
                      'brdf_dict': brdf}
        collect_ancillary(grn_con, res_group[GroupName.SAT_SOL_GROUP.value],
                          nbar_paths, ecmwf_path, invariant_fname,
                          vertices, root, compression, filter_opts)

        # atmospherics
        log.info('Atmospherics')

        ancillary_group = root[GroupName.ANCILLARY_GROUP.value]

        # satellite/solar angles and lon/lat for a resolution group
        sat_sol_grp = res_group[GroupName.SAT_SOL_GROUP.value]
        lon_lat_grp = res_group[GroupName.LON_LAT_GROUP.value]

        # TODO: supported acqs in different groups pointing to different response funcs
        json_data, _ = format_json(acqs, ancillary_group, sat_sol_grp,
                                   lon_lat_grp, workflow, root)
        # do we run psf here or at the end?
        run_modtran54(acqs, json_data, nvertices, modtran54_exe)
        # atmospheric inputs group
        inputs_grp = root[GroupName.ATMOSPHERIC_INPUTS_GRP.value]

        # radiative transfer for each point and albedo
        for key in json_data:
            point, albedo = key

            log.info('Radiative-Transfer', point=point, albedo=albedo.value)

            with tempfile.TemporaryDirectory() as tmpdir:

                prepare_modtran(acqs, point, [albedo], tmpdir)

                point_dir = pjoin(tmpdir, POINT_FMT.format(p=point))
                workdir = pjoin(point_dir, ALBEDO_FMT.format(a=albedo.value))

                json_mod_infile = pjoin(tmpdir, json_fmt.format(p=point, a=albedo.value))

                with open(json_mod_infile, 'w') as src:
                    json_dict = json_data[key]

                    if albedo == Albedos.ALBEDO_TH:

                        json_dict["MODTRAN"][0]["MODTRANINPUT"]["SPECTRAL"]["FILTNM"] = \
                            "%s/%s" % (workdir, json_dict["MODTRAN"][0]["MODTRANINPUT"]["SPECTRAL"]["FILTNM"])
                        json_dict["MODTRAN"][1]["MODTRANINPUT"]["SPECTRAL"]["FILTNM"] = \
                            "%s/%s" % (workdir, json_dict["MODTRAN"][1]["MODTRANINPUT"]["SPECTRAL"]["FILTNM"])

                    else:

                        json_dict["MODTRAN"][0]["MODTRANINPUT"]["SPECTRAL"]["FILTNM"] = \
                            "%s/%s" % (workdir, json_dict["MODTRAN"][0]["MODTRANINPUT"]["SPECTRAL"]["FILTNM"])

                    json.dump(json_dict, src, cls=JsonEncoder, indent=4)

                run_modtran(acqs, inputs_grp, workflow, nvertices, point, [albedo],
                            modtran_exe, tmpdir, root, compression, filter_opts)

        # atmospheric coefficients
        log.info('Coefficients')
        results_group = root[GroupName.ATMOSPHERIC_RESULTS_GRP.value]
        calculate_coefficients(results_group, root, compression, filter_opts)
        esun_values = {}
        # interpolate coefficients
        for grp_name in container.supported_groups:
            log = STATUS_LOGGER.bind(level1=container.label, granule=granule,
                                     granule_group=grp_name)
            log.info('Interpolation')

            # acquisitions and available bands for the current group level
            acqs = container.get_acquisitions(granule=granule, group=grp_name)
            nbar_acqs = [acq for acq in acqs if
                         acq.band_type == BandType.REFLECTIVE]
            sbt_acqs = [acq for acq in acqs if
                        acq.band_type == BandType.THERMAL]

            res_group = root[grp_name]
            sat_sol_grp = res_group[GroupName.SAT_SOL_GROUP.value]
            comp_grp = root[GroupName.COEFFICIENTS_GROUP.value]

            for coefficient in workflow.atmos_coefficients:
                if coefficient is AtmosphericCoefficients.ESUN:
                    continue
                if coefficient in Workflow.NBAR.atmos_coefficients:
                    band_acqs = nbar_acqs
                else:
                    band_acqs = sbt_acqs

                for acq in band_acqs:
                    log.info('Interpolate', band_id=acq.band_id,
                             coefficient=coefficient.value)
                    interpolate(acq, coefficient, ancillary_group, sat_sol_grp,
                                comp_grp, res_group, compression, filter_opts,
                                method)

            # standardised products
            band_acqs = []

            if workflow in (Workflow.STANDARD, Workflow.NBAR):
                band_acqs.extend(nbar_acqs)

            if workflow in (Workflow.STANDARD, Workflow.SBT):
                band_acqs.extend(sbt_acqs)

            for acq in band_acqs:
                interp_grp = res_group[GroupName.INTERP_GROUP.value]

                if acq.band_type == BandType.THERMAL:
                    log.info('SBT', band_id=acq.band_id)
                    surface_brightness_temperature(acq, interp_grp, res_group,
                                                   compression, filter_opts)
                else:
                    atmos_coefs = read_h5_table(comp_grp, DatasetName.NBAR_COEFFICIENTS.value)
                    esun_values[acq.band_name] = (
                        atmos_coefs
                        [atmos_coefs.band_name == acq.band_name]
                        [AtmosphericCoefficients.ESUN.value]
                    ).values[0]

                    slp_asp_grp = res_group[GroupName.SLP_ASP_GROUP.value]
                    rel_slp_asp = res_group[GroupName.REL_SLP_GROUP.value]
                    incident_grp = res_group[GroupName.INCIDENT_GROUP.value]
                    exiting_grp = res_group[GroupName.EXITING_GROUP.value]
                    shadow_grp = res_group[GroupName.SHADOW_GROUP.value]

                    log.info('Surface-Reflectance', band_id=acq.band_id)
                    calculate_reflectance(acq, interp_grp, sat_sol_grp,
                                          slp_asp_grp, rel_slp_asp,
                                          incident_grp, exiting_grp,
                                          shadow_grp, ancillary_group,
                                          rori, res_group, compression,
                                          filter_opts, normalized_solar_zenith,
                                          esun_values[acq.band_name])

            # pixel quality
            sbt_only = workflow == Workflow.SBT
            if pixel_quality and can_pq(level1, acq_parser_hint) and not sbt_only:
                run_pq(level1, res_group, landsea, res_group, compression, filter_opts, AP.NBAR, acq_parser_hint)
                run_pq(level1, res_group, landsea, res_group, compression, filter_opts, AP.NBART, acq_parser_hint)

        def get_band_acqs(grp_name):
            acqs = container.get_acquisitions(granule=granule, group=grp_name)
            nbar_acqs = [acq for acq in acqs if acq.band_type == BandType.REFLECTIVE]
            sbt_acqs = [acq for acq in acqs if acq.band_type == BandType.THERMAL]

            band_acqs = []
            if workflow in (Workflow.STANDARD, Workflow.NBAR):
                band_acqs.extend(nbar_acqs)

            if workflow in (Workflow.STANDARD, Workflow.SBT):
                band_acqs.extend(sbt_acqs)

            return band_acqs

        # wagl parameters
        parameters = {'vertices': list(vertices),
                      'method': method.value,
                      'rori': rori,
                      'buffer_distance': buffer_distance,
                      'normalized_solar_zenith': normalized_solar_zenith,
                      'esun': esun_values}

        # metadata yaml's
        metadata = root.create_group(DatasetName.METADATA.value)
        create_ard_yaml({grp_name: get_band_acqs(grp_name) for grp_name in container.supported_groups},
                        ancillary_group, metadata, parameters, workflow)
