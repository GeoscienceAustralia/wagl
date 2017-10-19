#!/usr/bin/env python

from os.path import join as pjoin
from posixpath import join as ppjoin
import logging
import tempfile
from structlog import wrap_logger
from structlog.processors import JSONRenderer
import h5py

from gaip.acquisition import acquisitions
from gaip.ancillary import collect_ancillary, aggregate_ancillary
from gaip import constants
from gaip.constants import GroupName, Model, BandType
from gaip.constants import ALBEDO_FMT, POINT_FMT, POINT_ALBEDO_FMT
from gaip.dsm import get_dsm
from gaip.incident_exiting_angles import incident_angles, exiting_angles
from gaip.incident_exiting_angles import relative_azimuth_slope
from gaip.interpolation import interpolate
from gaip.longitude_latitude_arrays import create_lon_lat_grids
from gaip.metadata import create_ard_yaml
from gaip.modtran import format_tp5, prepare_modtran, run_modtran
from gaip.modtran import calculate_coefficients
from gaip.reflectance import calculate_reflectance
from gaip.satellite_solar_angles import calculate_angles
from gaip.terrain_shadow_masks import self_shadow, calculate_cast_shadow
from gaip.terrain_shadow_masks import combine_shadow_masks
from gaip.slope_aspect import slope_aspect_arrays
from gaip.temperature import surface_brightness_temperature


LOG = wrap_logger(logging.getLogger('gaip-status'),
                  processors=[JSONRenderer(indent=1, sort_keys=True)])


def get_buffer(group):
    buf = {'product': 250,
           'R10m': 700,
           'R20m': 350,
           'R60m': 120}
    return buf[group]


def card4l(level1, model, vertices, method, pixel_quality, landsea, tle_path,
           aerosol_fname, brdf_path, brdf_premodis_path, ozone_path,
           water_vapour_path, dem_path, dsm_fname, invariant_fname,
           modtran_exe, out_fname, ecmwf_path=None, rori=0.52,
           compression='lzf', y_tile=100):
    """
    CEOS Analysis Ready Data for Land.
    A workflow for producing standardised products that meet the
    CARD4L specification.
    """
    tp5_fmt = pjoin(POINT_FMT, ALBEDO_FMT, ''.join([POINT_ALBEDO_FMT, '.tp5']))
    nvertices = vertices[0] * vertices[1]

    scene = acquisitions(level1)

    with h5py.File(out_fname, 'w') as fid:
        for grn_name in scene.granules:
            if grn_name is None:
                granule_group = fid['/']
            else:
                granule_group = fid.create_group(grn_name)

            for grp_name in scene.groups:
                log = LOG.bind(scene=scene.label, granule=grn_name,
                               granule_group=grp_name)
                group = granule_group.create_group(grp_name)
                acqs = scene.get_acquisitions(granule=grn_name, group=grp_name)

                # longitude and latitude
                log.info('Latitude-Longitude')
                create_lon_lat_grids(acqs[0].gridded_geo_box(), group,
                                     compression=compression, y_tile=y_tile)

                # satellite and solar angles
                log.info('Satellite-Solar-Angles')
                calculate_angles(acqs[0],
                                 group[GroupName.lon_lat_group.value], group,
                                 compression, tle_path, y_tile)

                if model == Model.standard or model == model.nbar:

                    # DEM
                    log.info('DEM-retriveal')
                    get_dsm(acqs[0], dsm_fname, get_buffer(grp_name), group,
                            compression, y_tile)

                    # slope & aspect
                    log.info('Slope-Aspect')
                    slope_aspect_arrays(acqs[0],
                                        group[GroupName.elevation_group.value],
                                        get_buffer(grp_name), group,
                                        compression, y_tile)

                    # incident angles
                    log.info('Incident-Angles')
                    incident_angles(group[GroupName.sat_sol_group.value],
                                    group[GroupName.slp_asp_group.value],
                                    group, compression, y_tile)

                    # exiting angles
                    log.info('Exiting-Angles')
                    exiting_angles(group[GroupName.sat_sol_group.value],
                                   group[GroupName.slp_asp_group.value],
                                   group, compression, y_tile)

                    # relative azimuth slope
                    log.info('Relative-Azimuth-Angles')
                    incident_group_name = GroupName.incident_group.value
                    exiting_group_name = GroupName.exiting_group.value
                    relative_azimuth_slope(group[incident_group_name],
                                           group[exiting_group_name],
                                           group, compression, y_tile)

                    # self shadow
                    log.info('Self-Shadow')
                    self_shadow(group[incident_group_name],
                                group[exiting_group_name], group, compression,
                                y_tile)

                    # cast shadow solar source direction
                    log.info('Cast-Shadow-Solar-Direction')
                    dsm_group_name = GroupName.elevation_group.value
                    calculate_cast_shadow(acqs[0], group[dsm_group_name],
                                          group[GroupName.sat_sol_group.value],
                                          get_buffer(grp_name), 500, 500, 
                                          group, compression, y_tile) 

                    # cast shadow satellite source direction
                    log.info('Cast-Shadow-Satellite-Direction')
                    calculate_cast_shadow(acqs[0], group[dsm_group_name],
                                          group[GroupName.sat_sol_group.value],
                                          get_buffer(grp_name), 500, 500, 
                                          group, compression, y_tile, False) 

                    # combined shadow masks
                    log.info('Combined-Shadow')
                    combine_shadow_masks(group[GroupName.shadow_group.value],
                                         group[GroupName.shadow_group.value],
                                         group[GroupName.shadow_group.value],
                                         group, compression, y_tile)

            # nbar and sbt ancillary
            LOG.info('Ancillary-Retrieval', scene=scene.label,
                     granule=grn_name, granule_group=None)
            nbar_paths = {'aerosol_fname': aerosol_fname,
                          'water_vapour_path': water_vapour_path,
                          'ozone_path': ozone_path,
                          'dem_path': dem_path,
                          'brdf_path': brdf_path,
                          'brdf_premodis_path': brdf_premodis_path}
            collect_ancillary(acqs[0], group[GroupName.sat_sol_group.value], 
                              nbar_paths, ecmwf_path, invariant_fname,
                              vertices, granule_group, compression)

        if scene.tiled:
            LOG.info('Aggregate-Ancillary', scene=scene.label,
                     granule='All Granules', granule_group=None)
            granule_groups = [fid[granule] for granule in scene.granules]
            aggregate_ancillary(granule_groups)

        # atmospherics
        for grn_name in scene.granules:
            log = LOG.bind(scene=scene.label, granule=grn_name,
                           granule_group=None)
            log.info('Atmospherics')

            granule_group = fid[scene.get_root(granule=grn_name)]

            # any resolution group is fine
            grp_name = scene.groups[0]
            acqs = scene.get_acquisitions(granule=grn_name, group=grp_name)
            root_path = ppjoin(scene.get_root(granule=grn_name), grp_name)

            # TODO check that the average ancilary group can be parsed to reflectance and other functions
            if scene.tiled:
                ancillary_group = granule_group[GroupName.ancillary_group.value]
            else:
                ancillary_group = fid[GroupName.ancillary_group.value]

            # satellite/solar angles and lon/lat for a resolution group
            pth = ppjoin(root_path, GroupName.sat_sol_group.value)
            sat_sol_grp = granule_group[pth]
            pth = ppjoin(root_path, GroupName.lon_lat_group.value)
            lon_lat_grp = granule_group[pth]

            # tp5 files
            tp5_data, _ = format_tp5(acqs, ancillary_group, sat_sol_grp,
                                     lon_lat_grp, model, granule_group)

            # atmospheric inputs group
            inputs_grp = granule_group[GroupName.atmospheric_inputs_grp.value]

            # radiative transfer for each point and albedo
            for key in tp5_data:
                point, albedo = key

                log.info('Radiative-Transfer', point=point, albedo=albedo)
                with tempfile.TemporaryDirectory() as tmpdir:

                    prepare_modtran(acqs, point, [albedo], tmpdir, modtran_exe)

                    # tp5 data
                    fname = pjoin(tmpdir, tp5_fmt.format(p=point, a=albedo))
                    with open(fname, 'w') as src:
                        src.writelines(tp5_data[key])


                    run_modtran(acqs, inputs_grp, model, nvertices, point,
                                [albedo], modtran_exe, tmpdir,
                                granule_group, compression)

            # coefficients
            log.info('Coefficients')
            pth = GroupName.atmospheric_results_grp.value
            results_group = granule_group[pth]
            calculate_coefficients(results_group, granule_group, compression)

            # interpolate coefficients
            for grp_name in scene.groups:
                log = LOG.bind(scene=scene.label, granule=grn_name,
                               granule_group=grp_name)
                log.info('Interpolation')

                # acquisitions and available bands for the current group level
                acqs = scene.get_acquisitions(granule=grn_name, group=grp_name)
                nbar_acqs = [acq for acq in acqs if
                             acq.band_type == BandType.Reflective]
                sbt_acqs = [acq for acq in acqs if
                            acq.band_type == BandType.Thermal]


                group = granule_group[grp_name]
                sat_sol_grp = group[GroupName.sat_sol_group.value]
                coef_grp = granule_group[GroupName.coefficients_group.value]

                for factor in model.factors:
                    if factor in Model.nbar.factors:
                        band_acqs = nbar_acqs
                    else:
                        band_acqs = sbt_acqs

                    for acq in band_acqs:
                        log.info('Interpolate', band_id=acq.band_id,
                                 factor=factor)
                        interpolate(acq, factor, ancillary_group, sat_sol_grp,
                                    coef_grp, group, compression, y_tile,
                                    method)

                # standardised products
                band_acqs = []
                if model == Model.standard or model == model.nbar:
                    band_acqs.extend(nbar_acqs)

                if model == Model.standard or model == model.sbt:
                    band_acqs.extend(sbt_acqs)

                for acq in band_acqs:
                    interp_grp = group[GroupName.interp_group.value]

                    if acq.band_type == BandType.Thermal:
                        log.info('SBT', band_id=acq.band_id)
                        surface_brightness_temperature(acq, interp_grp, group,
                                                       compression, y_tile)
                    else:
                        slp_asp_grp = group[GroupName.slp_asp_group.value]
                        rel_slp_asp = group[GroupName.rel_slp_group.value]
                        incident_grp = group[GroupName.incident_group.value]
                        exiting_grp = group[GroupName.exiting_group.value]
                        shadow_grp = group[GroupName.shadow_group.value]
                        log.info('Surface-Reflectance',
                                 band_id=acq.band_id)
                        calculate_reflectance(acq, interp_grp, sat_sol_grp,
                                              slp_asp_grp, rel_slp_asp,
                                              incident_grp, exiting_grp,
                                              shadow_grp, ancillary_group,
                                              rori, group, compression, y_tile)

                # metadata yaml's
                if model == Model.standard or model == model.nbar:
                    create_ard_yaml(band_acqs, ancillary_group, group)

                if model == Model.standard or model == model.sbt:
                    create_ard_yaml(band_acqs, ancillary_group, group, True)
