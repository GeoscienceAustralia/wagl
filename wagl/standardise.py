#!/usr/bin/env python

from os.path import join as pjoin
import logging
import tempfile

from posixpath import join as ppjoin
from structlog import wrap_logger
from structlog.processors import JSONRenderer
import h5py

from wagl.acquisition import acquisitions
from wagl.ancillary import collect_ancillary
from wagl.constants import ArdProducts as AP, GroupName, Model, BandType
from wagl.constants import ALBEDO_FMT, POINT_FMT, POINT_ALBEDO_FMT
from wagl.dsm import get_dsm
from wagl.incident_exiting_angles import incident_angles, exiting_angles
from wagl.incident_exiting_angles import relative_azimuth_slope
from wagl.interpolation import interpolate
from wagl.longitude_latitude_arrays import create_lon_lat_grids
from wagl.metadata import create_ard_yaml
from wagl.modtran import format_tp5, prepare_modtran, run_modtran
from wagl.modtran import calculate_coefficients
from wagl.reflectance import calculate_reflectance
from wagl.satellite_solar_angles import calculate_angles
from wagl.terrain_shadow_masks import self_shadow, calculate_cast_shadow
from wagl.terrain_shadow_masks import combine_shadow_masks
from wagl.slope_aspect import slope_aspect_arrays
from wagl.temperature import surface_brightness_temperature
from wagl.pq import can_pq, run_pq

LOG = wrap_logger(logging.getLogger('wagl-status'),
                  processors=[JSONRenderer(indent=1, sort_keys=True)])


# pylint disable=too-many-arguments
def card4l(level1, granule, model, vertices, method, pixel_quality, landsea,
           tle_path, aerosol, brdf_path, brdf_premodis_path, ozone_path,
           water_vapour, dem_path, dsm_fname, invariant_fname, modtran_exe,
           out_fname, ecmwf_path=None, rori=0.52, buffer_distance=8000,
           compression='lzf', acq_parser_hint=None):
    """
    CEOS Analysis Ready Data for Land.
    A workflow for producing standardised products that meet the
    CARD4L specification.
    """
    tp5_fmt = pjoin(POINT_FMT, ALBEDO_FMT, ''.join([POINT_ALBEDO_FMT, '.tp5']))
    nvertices = vertices[0] * vertices[1]

    container = acquisitions(level1, hint=acq_parser_hint)

    # TODO: pass through an acquisitions container rather than pathname
    with h5py.File(out_fname, 'w') as fid:
        fid.attrs['level1_uri'] = level1

        for grp_name in container.supported_groups:
            log = LOG.bind(level1=container.label, granule=granule,
                           granule_group=grp_name)

            # root group for a given granule and resolution group
            root = fid.create_group(ppjoin(granule, grp_name))
            acqs = container.get_acquisitions(granule=granule, group=grp_name)

            # longitude and latitude
            log.info('Latitude-Longitude')
            create_lon_lat_grids(acqs[0], root, compression)

            # satellite and solar angles
            log.info('Satellite-Solar-Angles')
            calculate_angles(acqs[0], root[GroupName.lon_lat_group.value],
                             root, compression, tle_path)

            if model == Model.standard or model == model.nbar:

                # DEM
                log.info('DEM-retriveal')
                get_dsm(acqs[0], dsm_fname, buffer_distance, root, compression)

                # slope & aspect
                log.info('Slope-Aspect')
                slope_aspect_arrays(acqs[0],
                                    root[GroupName.elevation_group.value],
                                    buffer_distance, root, compression)

                # incident angles
                log.info('Incident-Angles')
                incident_angles(root[GroupName.sat_sol_group.value],
                                root[GroupName.slp_asp_group.value],
                                root, compression)

                # exiting angles
                log.info('Exiting-Angles')
                exiting_angles(root[GroupName.sat_sol_group.value],
                               root[GroupName.slp_asp_group.value],
                               root, compression)

                # relative azimuth slope
                log.info('Relative-Azimuth-Angles')
                incident_group_name = GroupName.incident_group.value
                exiting_group_name = GroupName.exiting_group.value
                relative_azimuth_slope(root[incident_group_name],
                                       root[exiting_group_name],
                                       root, compression)

                # self shadow
                log.info('Self-Shadow')
                self_shadow(root[incident_group_name],
                            root[exiting_group_name], root, compression)

                # cast shadow solar source direction
                log.info('Cast-Shadow-Solar-Direction')
                dsm_group_name = GroupName.elevation_group.value
                calculate_cast_shadow(acqs[0], root[dsm_group_name],
                                      root[GroupName.sat_sol_group.value],
                                      buffer_distance, root, compression)

                # cast shadow satellite source direction
                log.info('Cast-Shadow-Satellite-Direction')
                calculate_cast_shadow(acqs[0], root[dsm_group_name],
                                      root[GroupName.sat_sol_group.value],
                                      buffer_distance, root, compression, False)

                # combined shadow masks
                log.info('Combined-Shadow')
                combine_shadow_masks(root[GroupName.shadow_group.value],
                                     root[GroupName.shadow_group.value],
                                     root[GroupName.shadow_group.value],
                                     root, compression)

        # nbar and sbt ancillary
        log = LOG.bind(level1=container.label, granule=granule,
                       granule_group=None)

        # granule root group
        root = fid[granule]

        # get the highest resoltion group cotaining supported bands
        acqs, grp_name = container.get_highest_resolution(granule=granule)

        grn_con = container.get_granule(granule=granule, container=True)
        res_group = root[grp_name]

        log.info('Ancillary-Retrieval')
        nbar_paths = {'aerosol_dict': aerosol,
                      'water_vapour_dict': water_vapour,
                      'ozone_path': ozone_path,
                      'dem_path': dem_path,
                      'brdf_path': brdf_path,
                      'brdf_premodis_path': brdf_premodis_path}
        collect_ancillary(grn_con, res_group[GroupName.sat_sol_group.value],
                          nbar_paths, ecmwf_path, invariant_fname,
                          vertices, root, compression)

        # atmospherics
        log.info('Atmospherics')

        ancillary_group = root[GroupName.ancillary_group.value]

        # satellite/solar angles and lon/lat for a resolution group
        sat_sol_grp = res_group[GroupName.sat_sol_group.value]
        lon_lat_grp = res_group[GroupName.lon_lat_group.value]

        # TODO: supported acqs in different groups pointing to different response funcs
        # tp5 files
        tp5_data, _ = format_tp5(acqs, ancillary_group, sat_sol_grp,
                                 lon_lat_grp, model, root)

        # atmospheric inputs group
        inputs_grp = root[GroupName.atmospheric_inputs_grp.value]

        # radiative transfer for each point and albedo
        for key in tp5_data:
            point, albedo = key

            log.info('Radiative-Transfer', point=point, albedo=albedo.value)
            with tempfile.TemporaryDirectory() as tmpdir:

                prepare_modtran(acqs, point, [albedo], tmpdir, modtran_exe)

                # tp5 data
                fname = pjoin(tmpdir,
                              tp5_fmt.format(p=point, a=albedo.value))
                with open(fname, 'w') as src:
                    src.writelines(tp5_data[key])

                run_modtran(acqs, inputs_grp, model, nvertices, point,
                            [albedo], modtran_exe, tmpdir, root, compression)

        # atmospheric coefficients
        log.info('Coefficients')
        results_group = root[GroupName.atmospheric_results_grp.value]
        calculate_coefficients(results_group, root, compression)

        # interpolate coefficients
        for grp_name in container.supported_groups:
            log = LOG.bind(level1=container.label, granule=granule,
                           granule_group=grp_name)
            log.info('Interpolation')

            # acquisitions and available bands for the current group level
            acqs = container.get_acquisitions(granule=granule, group=grp_name)
            nbar_acqs = [acq for acq in acqs if
                         acq.band_type == BandType.Reflective]
            sbt_acqs = [acq for acq in acqs if
                        acq.band_type == BandType.Thermal]

            res_group = root[grp_name]
            sat_sol_grp = res_group[GroupName.sat_sol_group.value]
            comp_grp = root[GroupName.coefficients_group.value]

            for coefficient in model.atmos_coefficients:
                if coefficient in Model.nbar.atmos_coefficients:
                    band_acqs = nbar_acqs
                else:
                    band_acqs = sbt_acqs

                for acq in band_acqs:
                    log.info('Interpolate', band_id=acq.band_id,
                             coefficient=coefficient.value)
                    interpolate(acq, coefficient, ancillary_group, sat_sol_grp,
                                comp_grp, res_group, compression, method)

            # standardised products
            band_acqs = []
            if model == Model.standard or model == model.nbar:
                band_acqs.extend(nbar_acqs)

            if model == Model.standard or model == model.sbt:
                band_acqs.extend(sbt_acqs)

            for acq in band_acqs:
                interp_grp = res_group[GroupName.interp_group.value]

                if acq.band_type == BandType.Thermal:
                    log.info('SBT', band_id=acq.band_id)
                    surface_brightness_temperature(acq, interp_grp, res_group,
                                                   compression)
                else:
                    slp_asp_grp = res_group[GroupName.slp_asp_group.value]
                    rel_slp_asp = res_group[GroupName.rel_slp_group.value]
                    incident_grp = res_group[GroupName.incident_group.value]
                    exiting_grp = res_group[GroupName.exiting_group.value]
                    shadow_grp = res_group[GroupName.shadow_group.value]

                    log.info('Surface-Reflectance', band_id=acq.band_id)
                    calculate_reflectance(acq, interp_grp, sat_sol_grp,
                                          slp_asp_grp, rel_slp_asp,
                                          incident_grp, exiting_grp,
                                          shadow_grp, ancillary_group,
                                          rori, res_group, compression)

            # metadata yaml's
            if model == Model.standard or model == Model.nbar:
                create_ard_yaml(band_acqs, ancillary_group, res_group)

            if model == Model.standard or model == Model.sbt:
                create_ard_yaml(band_acqs, ancillary_group, res_group, True)

            # pixel quality
            sbt_only = model == Model.sbt
            if pixel_quality and can_pq(level1, acq_parser_hint) and not sbt_only:
                run_pq(level1, res_group, landsea, res_group, compression, AP.nbar, acq_parser_hint)
                run_pq(level1, res_group, landsea, res_group, compression, AP.nbart, acq_parser_hint)
