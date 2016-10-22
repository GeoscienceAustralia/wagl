#!/usr/bin/env python
"""
NBAR Workflow
-------------

Workflow settings can be configured in `nbar.cfg` file.

"""
# pylint: disable=missing-docstring,no-init,too-many-function-args
# pylint: disable=too-many-locals

from datetime import datetime as dt
import cPickle as pickle
import os
from os.path import join as pjoin, basename, dirname, exists
import subprocess
import logging
import shutil
import tempfile
import argparse
import yaml
from yaml.representer import Representer
import luigi
from pathlib import Path
import numpy
from eodatasets.run import package_newly_processed_data_folder
from eodatasets.drivers import PACKAGE_DRIVERS
from eodatasets import type as ptype
import gaip


def save(target, value):
    """Save `value` to `target` where `target` is a `luigi.Target` object. If
    the target filename ends with `pkl` then pickle the data. Otherwise, save
    as text."""
    with target.open('w') as outfile:
        if target.fn.endswith('pkl'):
            pickle.dump(value, outfile)
        else:
            print >>outfile, value


def load(target):
    """Load data from `target` where `target` is a `luigi.Target`."""
    if not target.fn.endswith('pkl'):
        raise IOError('Cannot load non-pickled object')
    with target.open('r') as infile:
        return pickle.load(infile)


def load_value(target):
    """Load the value from `target`."""
    if isinstance(target, str):
        target = luigi.LocalTarget(target)
    data = load(target)
    try:
        return data['value']
    except KeyError:
        return data


class CreateWorkingDirectoryTree(luigi.Task):

    """Creates the output working directory tree."""

    level1 = luigi.Parameter()
    nbar_root = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        out_path = pjoin(self.nbar_root, CONFIG.get('work', 'targets_root'))
        target = pjoin(out_path, 'CreateWorkingDirectoryTree.task')
        return luigi.LocalTarget(target)

    def run(self):
        container = gaip.acquisitions(self.level1)

        for granule in container.granules:
            out_path = container.get_root(self.nbar_root, granule=granule)

            # base/top level output targets
            targets_dir = CONFIG.get('work', 'targets_root')
            base_targets = pjoin(out_path, targets_dir)
            if not exists(base_targets):
                os.makedirs(base_targets)

            # modtran directory tree
            modtran_exe_root = CONFIG.get('modtran', 'root')
            modtran_root = pjoin(out_path, CONFIG.get('work', 'modtran_root'))
            input_format = CONFIG.get('modtran', 'input_format')
            workpath_format = CONFIG.get('modtran', 'workpath_format')
            coords = CONFIG.get('modtran', 'coords').split(',')
            albedos = CONFIG.get('modtran', 'albedos').split(',')

            gaip.create_modtran_dirs(coords, albedos, modtran_root,
                                     modtran_exe_root,
                                     workpath_format,
                                     input_format)

            # brdf intermediates
            brdf_path = pjoin(out_path, CONFIG.get('work', 'brdf_root'))
            if not exists(brdf_path):
                os.makedirs(brdf_path)


            for group in container.groups:
                grp_path = pjoin(out_path, group)

                # group level output targets
                grp_targets = pjoin(grp_path, targets_dir)

                # terrain correction intermediates
                tc_path = pjoin(grp_path, CONFIG.get('work', 'tc_root'))

                # reflectance outputs
                rfl_path = pjoin(grp_path,
                                 CONFIG.get('work', 'reflectance_root'))


                # bilinear
                bil_path = pjoin(grp_path, CONFIG.get('work', 'bilinear_root'))

                if not exists(grp_targets):
                    os.makedirs(grp_targets)

                if not exists(tc_path):
                    os.makedirs(tc_path)

                if not exists(rfl_path):
                    os.makedirs(rfl_path)

                if not exists(bil_path):
                    os.makedirs(bil_path)

        save(self.output(), 'completed')


class GetElevationAncillaryData(luigi.Task):

    """Get ancillary elevation data."""

    level1 = luigi.Parameter()
    nbar_root = luigi.Parameter()
    granule = luigi.Parameter()

    def requires(self):
        return [CreateWorkingDirectoryTree(self.level1, self.nbar_root)]

    def output(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, granule=self.granule)
        out_path = pjoin(out_path, CONFIG.get('work', 'targets_root'))
        target = pjoin(out_path, 'GetElevationAncillaryData.task')
        return luigi.LocalTarget(target)

    def run(self):
        container = gaip.acquisitions(self.level1)
        acqs = container.get_acquisitions(granule=self.granule)
        geobox = acqs[0].gridded_geo_box()

        dem_path = CONFIG.get('ancillary', 'dem_path')
        out_path = container.get_root(self.nbar_root, granule=self.granule)
        value = gaip.get_elevation_data(geobox.centre_lonlat, dem_path)

        out_fname = pjoin(out_path, CONFIG.get('work', 'dem_fname'))
        save(luigi.LocalTarget(out_fname), value)
        save(self.output(), 'completed')


class GetOzoneAncillaryData(luigi.Task):

    """Get ancillary ozone data."""

    level1 = luigi.Parameter()
    nbar_root = luigi.Parameter()
    granule = luigi.Parameter()

    def requires(self):
        return [CreateWorkingDirectoryTree(self.level1, self.nbar_root)]

    def output(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, granule=self.granule)
        out_path = pjoin(out_path, CONFIG.get('work', 'targets_root'))
        target = pjoin(out_path, 'GetOzoneAncillaryData.task')
        return luigi.LocalTarget(target)

    def run(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, granule=self.granule)
        acqs = container.get_acquisitions(granule=self.granule)
        geobox = acqs[0].gridded_geo_box()

        ozone_path = CONFIG.get('ancillary', 'ozone_path')
        centre = geobox.centre_lonlat
        centre_datetime = acqs[0].scene_center_datetime
        value = gaip.get_ozone_data(ozone_path, centre, centre_datetime)

        out_fname = pjoin(out_path, CONFIG.get('work', 'ozone_fname'))
        save(luigi.LocalTarget(out_fname), value)
        save(self.output(), 'completed')


class GetWaterVapourAncillaryData(luigi.Task):

    """Get ancillary water vapour data."""

    level1 = luigi.Parameter()
    nbar_root = luigi.Parameter()
    granule = luigi.Parameter()

    def requires(self):
        return [CreateWorkingDirectoryTree(self.level1, self.nbar_root)]

    def output(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, granule=self.granule)
        out_path = pjoin(out_path, CONFIG.get('work', 'targets_root'))
        target = pjoin(out_path, 'GetWaterVapourAncillaryData.task')
        return luigi.LocalTarget(target)

    def run(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, granule=self.granule)
        acqs = container.get_acquisitions(granule=self.granule)

        vapour_path = CONFIG.get('ancillary', 'vapour_path')
        value = gaip.get_water_vapour(acqs[0], vapour_path)

        out_fname = pjoin(out_path, CONFIG.get('work', 'vapour_fname'))
        save(luigi.LocalTarget(out_fname), value)
        save(self.output(), 'completed')


class GetAerosolAncillaryData(luigi.Task):

    """Get ancillary aerosol data."""

    level1 = luigi.Parameter()
    nbar_root = luigi.Parameter()
    granule = luigi.Parameter()

    def requires(self):
        return [CreateWorkingDirectoryTree(self.level1, self.nbar_root)]

    def output(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, granule=self.granule)
        out_path = pjoin(out_path, CONFIG.get('work', 'targets_root'))
        target = pjoin(out_path, 'GetAerosolAncillaryData.task')
        return luigi.LocalTarget(target)

    def run(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, granule=self.granule)
        acqs = container.get_acquisitions(granule=self.granule)

        aerosol_path = CONFIG.get('ancillary', 'aerosol_path')
        value = gaip.get_aerosol_data(acqs[0], aerosol_path)
        # aerosol_path = CONFIG.get('ancillary', 'aerosol_fname') # version 2
        # value = gaip.get_aerosol_data_v2(acqs[0], aerosol_path) # version 2

        out_fname = pjoin(out_path, CONFIG.get('work', 'aerosol_fname'))
        save(luigi.LocalTarget(out_fname), value)
        save(self.output(), 'completed')


class GetBrdfAncillaryData(luigi.Task):

    """Get ancillary BRDF data."""

    level1 = luigi.Parameter()
    nbar_root = luigi.Parameter()
    granule = luigi.Parameter()

    def requires(self):
        return [CreateWorkingDirectoryTree(self.level1, self.nbar_root)]

    def output(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, granule=self.granule)
        out_path = pjoin(out_path, CONFIG.get('work', 'targets_root'))
        target = pjoin(out_path, 'GetBrdfAncillaryData.task')
        return luigi.LocalTarget(target)

    def run(self):
        container = gaip.acquisitions(self.level1)
        acqs = container.get_acquisitions(granule=self.granule)
        out_path = container.get_root(self.nbar_root, granule=self.granule)

        work_path = pjoin(out_path, CONFIG.get('work', 'brdf_root'))
        brdf_path = CONFIG.get('ancillary', 'brdf_path')
        brdf_premodis_path = CONFIG.get('ancillary', 'brdf_premodis_path')
        value = gaip.get_brdf_data(acqs[0], brdf_path, brdf_premodis_path,
                                   work_path)

        out_fname = pjoin(out_path, CONFIG.get('work', 'brdf_fname'))
        save(luigi.LocalTarget(out_fname), value)
        save(self.output(), 'completed')


class GetAncillaryData(luigi.Task):

    """Get all ancillary data. This a helper task."""

    level1 = luigi.Parameter()
    nbar_root = luigi.Parameter()
    granule = luigi.Parameter()

    def requires(self):
        args = [self.level1, self.nbar_root, self.granule]
        return [GetElevationAncillaryData(*args),
                GetOzoneAncillaryData(*args),
                GetWaterVapourAncillaryData(*args),
                GetAerosolAncillaryData(*args),
                GetBrdfAncillaryData(*args)]

    def complete(self):
        return all([t.complete() for t in self.requires()])


class CalculateLonGrid(luigi.Task):

    """Calculate the longitude grid."""

    level1 = luigi.Parameter()
    nbar_root = luigi.Parameter()
    granule = luigi.Parameter()
    group = luigi.Parameter()

    def requires(self):
        return [CreateWorkingDirectoryTree(self.level1, self.nbar_root)]

    def output(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, group=self.group,
                                      granule=self.granule)
        out_path = pjoin(out_path, CONFIG.get('work', 'targets_root'))
        target = pjoin(out_path, 'CalculateLonGrid.task')
        return luigi.LocalTarget(target)

    def run(self):
        container = gaip.acquisitions(self.level1)
        acqs = container.get_acquisitions(group=self.group,
                                          granule=self.granule)
        out_path = container.get_root(self.nbar_root, group=self.group,
                                      granule=self.granule)

        out_fname = pjoin(out_path, CONFIG.get('work', 'lon_grid_fname'))
        gaip.create_lon_grid(acqs[0], out_fname)

        save(self.output(), 'completed')


class CalculateLatGrid(luigi.Task):

    """Calculate the latitude grid."""

    level1 = luigi.Parameter()
    nbar_root = luigi.Parameter()
    granule = luigi.Parameter()
    group = luigi.Parameter()

    def requires(self):
        return [CreateWorkingDirectoryTree(self.level1, self.nbar_root)]

    def output(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, group=self.group,
                                      granule=self.granule)
        out_path = pjoin(out_path, CONFIG.get('work', 'targets_root'))
        target = pjoin(out_path, 'CalculateLatGrid.task')
        return luigi.LocalTarget(target)

    def run(self):
        container = gaip.acquisitions(self.level1)
        acqs = container.get_acquisitions(group=self.group,
                                          granule=self.granule)
        out_path = container.get_root(self.nbar_root, group=self.group,
                                      granule=self.granule)

        out_fname = pjoin(out_path, CONFIG.get('work', 'lat_grid_fname'))
        gaip.create_lat_grid(acqs[0], out_fname)

        save(self.output(), 'completed')


class CalculateSatelliteAndSolarGrids(luigi.Task):

    """Calculate the satellite and solar grids."""

    level1 = luigi.Parameter()
    nbar_root = luigi.Parameter()
    granule = luigi.Parameter()
    group = luigi.Parameter()

    def requires(self):
        args = [self.level1, self.nbar_root, self.granule, self.group]
        return [CalculateLatGrid(*args),
                CalculateLonGrid(*args)]

    def output(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, group=self.group,
                                      granule=self.granule)
        out_path = pjoin(out_path, CONFIG.get('work', 'targets_root'))
        target = pjoin(out_path, 'CalculateSatelliteAndSolarGrids.task')
        return luigi.LocalTarget(target)

    def run(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, group=self.group,
                                      granule=self.granule)
        fnames = [CONFIG.get('work', 'sat_view_fname'),
                  CONFIG.get('work', 'sat_azimuth_fname'),
                  CONFIG.get('work', 'solar_zenith_fname'),
                  CONFIG.get('work', 'solar_azimuth_fname'),
                  CONFIG.get('work', 'relative_azimuth_fname'),
                  CONFIG.get('work', 'time_fname')]
        out_fnames = [pjoin(out_path, f) for f in fnames]
        centreline_fname = pjoin(out_path,
                                 CONFIG.get('work', 'centreline_fname'))
        coordinator_fname = pjoin(out_path,
                                  CONFIG.get('work', 'coordinator_fname'))
        boxline_fname = pjoin(out_path, CONFIG.get('work', 'boxline_fname'))
        lon_fname = pjoin(out_path, CONFIG.get('work', 'lon_grid_fname'))
        lat_fname = pjoin(out_path, CONFIG.get('work', 'lat_grid_fname'))

        acqs = container.get_acquisitions(group=self.group,
                                          granule=self.granule)

        geobox = acqs[0].gridded_geo_box()
        cols = acqs[0].samples
        view_max = acqs[0].maximum_view_angle

        (satellite_zenith, _, _, _, _, _, y_cent, x_cent, n_cent) = \
            gaip.calculate_angles(acqs[0], lon_fname, lat_fname,
                                  npoints=12, out_fnames=out_fnames)

        gaip.create_centreline_file(geobox, y_cent, x_cent, n_cent, cols,
                                    view_max=view_max,
                                    outfname=centreline_fname)

        gaip.create_boxline_file(satellite_zenith, y_cent, x_cent,
                                 boxline_fname=boxline_fname,
                                 max_angle=view_max,
                                 coordinator_fname=coordinator_fname)

        save(self.output(), 'completed')


class AggregateAncillary(luigi.Task):

    """Aggregates the ancillary data from each granule."""

    level1 = luigi.Parameter()
    nbar_root = luigi.Parameter()

    def requires(self):
        container = gaip.acquisitions(self.level1)
        tasks = []

        for granule in container.granules:
            args1 = [self.level1, self.nbar_root, granule]
            tasks.append(GetAncillaryData(args1))
            for group in container.groups:
                args2 = [self.level1, self.nbar_root, granule, group]
                tasks.append(CalculateSatelliteAndSolarGrids(*args2))
                tasks.append(CalculateLatGrid(*args2))
                tasks.append(CalculateLonGrid(*args2))

        return tasks

    def output(self):
        out_path = self.nbar_root
        out_path = pjoin(out_path, CONFIG.get('work', 'targets_root'))
        target = pjoin(out_path, 'AggregateAncillary.task')
        return luigi.LocalTarget(target)

    def run(self):
        container = gaip.acquisitions(self.level1)
        n_tiles = len(container.granules)
        out_path = self.nbar_root

        # initialise the mean result
        ozone = vapour = aerosol = elevation = 0.0

        # loop over each granule and retrieve the acnillary
        for granule in container.granules:
            grn_path = container.get_root(out_path, granule)

            # load the ancillary point values
            ozone_fname = pjoin(grn_path, CONFIG.get('work', 'ozone_fname'))
            vapour_fname = pjoin(grn_path, CONFIG.get('work', 'vapour_fname'))
            aerosol_fname = pjoin(grn_path, CONFIG.get('work',
                                                       'aerosol_fname'))
            elevation_fname = pjoin(grn_path, CONFIG.get('work', 'dem_fname'))

            ozone += load_value(ozone_fname)
            vapour += load_value(vapour_fname)
            aerosol += load_value(aerosol_fname)
            elevation += load_value(elevation_fname)

        ozone /= n_tiles
        vapour /= n_tiles
        aerosol /= n_tiles
        elevation /= n_tiles

        # out filenames
        ozone_outfname = pjoin(out_path, CONFIG.get('work', 'ozone_fname'))
        vapour_outfname = pjoin(out_path, CONFIG.get('work', 'vapour_fname'))
        aerosol_outfname = pjoin(out_path, CONFIG.get('work', 'aerosol_fname'))
        elevation_fname = pjoin(out_path, CONFIG.get('work', 'dem_fname'))

        # write the mean ancillary values
        data = {'data_source': 'granule_average'}
        data['value'] = ozone
        save(luigi.LocalTarget(ozone_outfname), data)
        data['value'] = vapour
        save(luigi.LocalTarget(vapour_outfname), data)
        data['value'] = aerosol
        save(luigi.LocalTarget(aerosol_outfname), data)
        data['value'] = elevation
        save(luigi.LocalTarget(elevation_fname), data)

        save(self.output(), 'completed')


class WriteTp5(luigi.Task):

    """Output the `tp5` formatted files."""

    level1 = luigi.Parameter()
    nbar_root = luigi.Parameter()
    granule = luigi.Parameter()

    def requires(self):
        # for consistancy, we'll wait for dependencies on all granules and
        # groups of acquisitions
        # current method requires to compute an average from all granules
        # if the scene is tiled up that way
        container = gaip.acquisitions(self.level1)
        tasks = []

        # are we dealing with tiles/granules?
        if container.tiled:
            tasks.append(AggregateAncillary(self.level1, self.nbar_root))

        for granule in container.granules:
            args1 = [self.level1, self.nbar_root, granule]
            tasks.append(GetAncillaryData(args1))
            for group in container.groups:
                args2 = [self.level1, self.nbar_root, granule, group]
                tasks.append(CalculateSatelliteAndSolarGrids(*args2))
                tasks.append(CalculateLatGrid(*args2))
                tasks.append(CalculateLonGrid(*args2))

        return tasks

    def output(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, granule=self.granule)
        out_path = pjoin(out_path, CONFIG.get('work', 'targets_root'))
        target = pjoin(out_path, 'WriteTp5.task')
        return luigi.LocalTarget(target)

    def run(self):
        container = gaip.acquisitions(self.level1)
        grn_path = container.get_root(self.nbar_root, granule=self.granule)

        coords = CONFIG.get('modtran', 'coords').split(',')
        albedos = CONFIG.get('modtran', 'albedos').split(',')
        output_format = CONFIG.get('write_tp5', 'output_format')
        workdir = pjoin(grn_path, CONFIG.get('work', 'modtran_root'))
        out_fname_format = pjoin(workdir, output_format)

        # get the first group name
        group = container.groups[0]
        grp_path = container.get_root(grn_path, group=group)

        # get the filenames for the coordinator,
        # satellite view zenith, azimuth, and latitude/longitude arrays
        coord_fname = pjoin(grp_path, CONFIG.get('work', 'coordinator_fname'))
        sat_view_fname = pjoin(grp_path, CONFIG.get('work', 'sat_view_fname'))
        sat_azi_fname = pjoin(grp_path, CONFIG.get('work',
                                                   'sat_azimuth_fname'))
        lon_fname = pjoin(grp_path, CONFIG.get('work', 'lon_grid_fname'))
        lat_fname = pjoin(grp_path, CONFIG.get('work', 'lat_grid_fname'))

        # load the ancillary point values
        out_path = self.nbar_root
        ozone_fname = pjoin(out_path, CONFIG.get('work', 'ozone_fname'))
        vapour_fname = pjoin(out_path, CONFIG.get('work', 'vapour_fname'))
        aerosol_fname = pjoin(out_path, CONFIG.get('work', 'aerosol_fname'))
        elevation_fname = pjoin(out_path, CONFIG.get('work', 'dem_fname'))
        ozone = load_value(ozone_fname)
        vapour = load_value(vapour_fname)
        aerosol = load_value(aerosol_fname)
        elevation = load_value(elevation_fname)

        # load an acquisition
        acq = container.get_acquisitions(group=group, granule=self.granule)[0]

        # run
        gaip.write_tp5(acq, coord_fname, sat_view_fname, sat_azi_fname,
                       lat_fname, lon_fname, ozone, vapour, aerosol, elevation,
                       coords, albedos, out_fname_format)

        save(self.output(), 'completed')


class LegacyOutputs(luigi.Task):

    """
    A task that isn't required for production, but can be useful for
    users undergoing validation. The files produced by this
    task are neccessary for her version of NBAR to run, whereas
    the production version skips producing file `a` and goes direct
    to file `b`, or simply is not used by the system and hence has
    been removed from production.
    """

    level1 = luigi.Parameter()
    nbar_root = luigi.Parameter()

    def requires(self):
        container = gaip.acquisitions(self.level1)
        tasks = []

        # are we dealing with tiles/granules?
        if container.tiled:
            tasks.append(AggregateAncillary(self.level1, self.nbar_root))

        for granule in container.granules:
            args1 = [self.level1, self.nbar_root, granule]
            tasks.append(GetAncillaryData(args1))
            for group in container.groups:
                args2 = [self.level1, self.nbar_root, granule, group]
                tasks.append(CalculateSatelliteAndSolarGrids(*args2))
                tasks.append(CalculateLatGrid(*args2))
                tasks.append(CalculateLonGrid(*args2))

        return tasks

    def output(self):
        out_path = self.nbar_root
        out_path = pjoin(out_path, CONFIG.get('work', 'targets_root'))
        target = pjoin(out_path, 'LegacyOutputs.task')
        return luigi.LocalTarget(target)

    def run(self):
        out_path = self.nbar_root
        container = gaip.acquisitions(self.level1)
        acqs = container.get_acqisitions()

        # satellite filter file
        satfilterpath = CONFIG.get('ancillary', 'satfilter_path')
        filter_fname = CONFIG.get('legacy_outputs', 'sat_filter_fname')
        out_fname = pjoin(out_path, filter_fname)
        gaip.create_satellite_filter_file(acqs, satfilterpath, out_fname)

        # modtran input file
        ozone_fname = pjoin(out_path, CONFIG.get('work', 'ozone_fname'))
        vapour_fname = pjoin(out_path, CONFIG.get('work', 'vapour_fname'))
        aerosol_fname = pjoin(out_path, CONFIG.get('work', 'aerosol_fname'))
        elevation_fname = pjoin(out_path, CONFIG.get('work', 'dem_fname'))
        mod_input_fname = CONFIG.get('legacy_outputs', 'modtran_input_fname')
        out_fname = pjoin(out_path, mod_input_fname)
        ozone = load_value(ozone_fname)
        vapour = load_value(vapour_fname)
        aerosol = load_value(aerosol_fname)
        elevation = load_value(elevation_fname)
        gaip.write_modtran_input(acqs, out_fname, ozone, vapour, aerosol,
                                 elevation)

        for granule in container.granules:
            for group in container.groups:
                acqs = container.get_acquisitions(group=group, granule=granule)
                grp_path = container.get_root(group=group, granule=granule)

                # modtran input files
                coordinator_fname = CONFIG.get('work', 'coordinator_fname')
                coordinator_fname = pjoin(grp_path, coordinator_fname)
                sat_view_zenith_fname = CONFIG.get('work', 'sat_view_fname')
                sat_view_zenith_fname = pjoin(grp_path, sat_view_zenith_fname)
                sat_azimuth_fname = CONFIG.get('work', 'sat_azimuth_fname')
                sat_azimuth_fname = pjoin(grp_path, sat_azimuth_fname)
                lon_grid_fname = CONFIG.get('work', 'lon_grid_fname')
                lon_grid_fname = pjoin(grp_path, lon_grid_fname)
                lat_grid_fname = CONFIG.get('work', 'lat_grid_fname')
                lat_grid_fname = pjoin(grp_path, lat_grid_fname)

                coords = CONFIG.get('modtran', 'coords').split(',')
                albedos = CONFIG.get('legacy_outputs', 'albedos').split(',')
                fname_format = CONFIG.get('legacy_outputs', 'output_format')

                out_fname_fmt = pjoin(grp_path, fname_format)
                gaip.write_modtran_inputs(acqs[0], coordinator_fname,
                                          sat_view_zenith_fname,
                                          sat_azimuth_fname,
                                          lat_grid_fname, lon_grid_fname,
                                          ozone, vapour, aerosol, elevation,
                                          coords, albedos, out_fname_fmt)

                # header angle file
                view_max = acqs[0].maximum_view_angle
                hdr_angle_fname = CONFIG.get('work', 'header_angle_fname')
                header_angle_fname = pjoin(grp_path, hdr_angle_fname)
                gaip.create_header_angle_file(acqs[0], view_max=view_max,
                                              outfname=header_angle_fname)

        save(self.output(), 'completed')


class PrepareModtranInput(luigi.Task):

    """Prepare MODTRAN inputs. This is a helper task."""

    level1 = luigi.Parameter()
    nbar_root = luigi.Parameter()

    def requires(self):
        args = [self.level1, self.nbar_root]
        # do we need to output the legacy files?
        if CONFIG.getboolean('legacy_outputs', 'required'):
            tasks = [LegacyOutputs(*args),
                     WriteTp5(*args)]
        else:
            tasks = [WriteTp5(*args)]

        return tasks

    def complete(self):
        return all([t.complete() for t in self.requires()])


class RunModtranCase(luigi.Task):

    """Run MODTRAN for a specific `coord` and `albedo`. This task is
    parameterised this way to allow parallel instances of MODTRAN to run."""

    level1 = luigi.Parameter()
    nbar_root = luigi.Parameter()
    granule = luigi.Parameter()
    coord = luigi.Parameter()
    albedo = luigi.Parameter()

    def requires(self):
        args = [self.level1, self.nbar_root, self.granule]
        return [PrepareModtranInput(*args)]

    def output(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, granule=self.granule)
        out_path = pjoin(out_path, CONFIG.get('work', 'targets_root'))
        target = pjoin(out_path, 'RunModtranCase_{coord}_{albedo}.task')
        target = target.format(coord=self.coord, albedo=self.albedo)
        return luigi.LocalTarget(target)

    def run(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, granule=self.granule)

        modtran_exe = CONFIG.get('modtran', 'exe')
        workpath_format = CONFIG.get('modtran', 'workpath_format')
        modtran_root = pjoin(out_path, CONFIG.get('work', 'modtran_root'))
        workpath = workpath_format.format(coord=self.coord, albedo=self.albedo)
        gaip.run_modtran(modtran_exe, pjoin(modtran_root, workpath))

        save(self.output(), 'completed')


class RunModtran(luigi.Task):

    """Run MODTRAN for all coords and albedos. This is a helper task."""

    level1 = luigi.Parameter()
    nbar_root = luigi.Parameter()
    granule = luigi.Parameter()

    def requires(self):
        coords = CONFIG.get('modtran', 'coords').split(',')
        albedos = CONFIG.get('modtran', 'albedos').split(',')
        reqs = [PrepareModtranInput(self.level1, self.nbar_root)]
        for coord in coords:
            for albedo in albedos:
                reqs.append(RunModtranCase(self.level1, self.nbar_root, 
                                           self.granule, coord, albedo))
        return reqs

    def complete(self):
        return all([t.complete() for t in self.requires()])


class RunAccumulateSolarIrradianceCase(luigi.Task):

    """
    Run calculate_solar_radiation for a given case, with case being
    a given albedo and coordinate.
    """

    level1 = luigi.Parameter()
    nbar_root = luigi.Parameter()
    granule = luigi.Parameter()
    coord = luigi.Parameter()
    albedo = luigi.Parameter()

    def requires(self):
        return [RunModtran(self.level1, self.nbar_root, self.granule)]

    def output(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, granule=self.granule)
        out_path = pjoin(out_path, CONFIG.get('work', 'targets_root'))
        task = 'RunAccumulateSolarIrradianceCase_{coord}_{albedo}.task'
        target = pjoin(out_path, task)
        target = target.format(coord=self.coord, albedo=self.albedo)
        return luigi.LocalTarget(target)

    def run(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, granule=self.granule)
        acqs = container.get_acquisitions(granule=self.granule)

        modtran_root = pjoin(out_path, CONFIG.get('work', 'modtran_root'))
        input_format = CONFIG.get('extract_flux', 'input_format')
        input_format = pjoin(modtran_root, input_format)
        flux_fname = input_format.format(coord=self.coord, albedo=self.albedo)

        satfilterpath = CONFIG.get('ancillary', 'satfilter_path')
        response_fname = pjoin(satfilterpath, acqs[0].spectral_filter_file)

        output_format = CONFIG.get('extract_flux', 'output_format')
        output_format = pjoin(modtran_root, output_format)
        out_fname = output_format.format(coord=self.coord, albedo=self.albedo)

        transmittance = True if self.albedo == 't' else False
        result = gaip.calculate_solar_radiation(flux_fname, response_fname,
                                                transmittance)
        result.to_csv(out_fname, index=False, sep='\t')

        save(self.output(), 'completed')


class AccumulateSolarIrradiance(luigi.Task):

    """
    Extract the flux data from the MODTRAN outputs, and calculate
    the accumulative solar irradiance for a given spectral
    response function.

    This is a helper class that kicks off individual coordinate albedo
    RunAccumulateSolarIrradianceCase tasks.
    """

    level1 = luigi.Parameter()
    nbar_root = luigi.Parameter()
    granule = luigi.Parameter()

    def requires(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, granule=self.granule)
        coords = CONFIG.get('modtran', 'coords').split(',')
        albedos = CONFIG.get('modtran', 'albedos').split(',')
        modtran_root = pjoin(out_path, CONFIG.get('work', 'modtran_root'))
        output_format = CONFIG.get('extract_flux', 'output_format')
        output_format = pjoin(modtran_root, output_format)

        reqs = []
        for coord in coords:
            for albedo in albedos:
                reqs.append(RunAccumulateSolarIrradianceCase(self.level1,
                                                             self.nbar_root, 
                                                             self.granule,
                                                             coord, albedo))
        return reqs

    def complete(self):
        return all([t.complete() for t in self.requires()])


class CalculateCoefficients(luigi.Task):

    """
    Calculate the atmospheric parameters needed by BRDF and atmospheric
    correction model.
    """

    level1 = luigi.Parameter()
    nbar_root = luigi.Parameter()
    granule = luigi.Parameter()

    def requires(self):
        args = [self.level1, self.nbar_root, self.granule]
        return [AccumulateSolarIrradiance(*args)]

    def output(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, granule=self.granule)
        out_path = pjoin(out_path, CONFIG.get('work', 'targets_root'))
        target = pjoin(out_path, 'CalculateCoefficients.task')
        return luigi.LocalTarget(target)

    def run(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, granule=self.granule)

        coords = CONFIG.get('modtran', 'coords').split(',')
        chn_input_format = CONFIG.get('coefficients', 'chn_input_format')
        dir_input_format = CONFIG.get('extract_flux', 'output_format')
        output_format1 = CONFIG.get('coefficients', 'output_format1')
        workpath = pjoin(out_path, CONFIG.get('work', 'modtran_root'))
        output_format2 = CONFIG.get('coefficients', 'output_format2')

        gaip.calculate_coefficients(coords, chn_input_format, dir_input_format,
                                    output_format1, output_format2,
                                    workpath)

        save(self.output(), 'completed')


class BilinearInterpolationBand(luigi.Task):
    """
    Runs the bilinear interpolation function for a given band.
    """

    level1 = luigi.Parameter()
    nbar_root = luigi.Parameter()
    granule = luigi.Parameter()
    group = luigi.Parameter()
    band_num = luigi.IntParameter()
    factor = luigi.Parameter()

    def requires(self):
        args1 = [self.level1, self.nbar_root, self.granule]
        args2 = [self.level1, self.nbar_root, self.granule, self.group]
        return [CalculateCoefficients(*args1),
                CalculateSatelliteAndSolarGrids(*args2)]

    def output(self):
        band_num = self.band_num
        factor = self.factor
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, group=self.group,
                                      granule=self.granule)
        out_path = pjoin(out_path, CONFIG.get('work', 'targets_root'))
        task = 'BilinearInterpolationBand_{band}_{factor}.task'
        target = pjoin(out_path, task.format(band=band_num, factor=factor))
        return luigi.LocalTarget(target)

    def run(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, group=self.group,
                                      granule=self.granule)
        coordinator = pjoin(out_path,
                            CONFIG.get('work', 'coordinator_fname'))
        boxline = pjoin(out_path,
                        CONFIG.get('work', 'boxline_fname'))
        centreline = pjoin(out_path,
                           CONFIG.get('work', 'centreline_fname'))
        input_format = CONFIG.get('coefficients', 'output_format2')
        output_format = CONFIG.get('bilinear', 'output_format')
        workpath = pjoin(out_path,
                         CONFIG.get('work', 'bilinear_root'))

        mod_work_path = pjoin(container.get_root(self.nbar_root,
                                                 granule=self.granule),
                              CONFIG.get('work', 'modtran_root'))
        input_format = pjoin(mod_work_path, input_format)

        acqs = container.get_acquisitions(group=self.group,
                                          granule=self.granule)

        # get the acquisition we wish to process
        acqs = [acq for acq in acqs if acq.band_num == self.band_num]

        _ = gaip.bilinear_interpolate(acqs, [self.factor], coordinator,
                                      boxline, centreline, input_format,
                                      output_format, workpath)

        save(self.output(), 'completed')


class BilinearInterpolation(luigi.Task):
    """
    Issues BilinearInterpolationBand tasks.
    This is a helper task.
    """

    level1 = luigi.Parameter()
    nbar_root = luigi.Parameter()
    granule = luigi.Parameter()
    group = luigi.Parameter()

    def requires(self):
        factors = CONFIG.get('bilinear', 'factors').split(',')
        container = gaip.acquisitions(self.level1)
        acqs = container.get_acquisitions(group=self.group,
                                          granule=self.granule)

        # Retrieve the satellite and sensor for the acquisition
        satellite = acqs[0].spacecraft_id
        sensor = acqs[0].sensor_id

        # Get the required nbar bands list for processing
        nbar_constants = gaip.constants.NBARConstants(satellite, sensor)
        bands_to_process = nbar_constants.get_nbar_lut()

        bands = [a.band_num for a in acqs if a.band_num in bands_to_process]
        tasks = []
        for factor in factors:
            for band in bands:
                tasks.append(BilinearInterpolationBand(self.level1,
                                                       self.nbar_root,
                                                       self.granule,
                                                       self.group,
                                                       band, factor))
        return tasks

    def output(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, group=self.group,
                                      granule=self.granule)
        out_path = pjoin(out_path, CONFIG.get('work', 'targets_root'))
        target = pjoin(out_path, 'BilinearInterpolation.task')
        return luigi.LocalTarget(target)

    def run(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, group=self.group,
                                      granule=self.granule)
        factors = CONFIG.get('bilinear', 'factors').split(',')
        input_format = CONFIG.get('coefficients', 'output_format2')
        output_format = CONFIG.get('bilinear', 'output_format')
        workpath = pjoin(out_path, CONFIG.get('work', 'bilinear_root'))
        input_format = pjoin(workpath, input_format)

        acqs = container.get_acquisitions(group=self.group,
                                          granule=self.granule)
        bands = [a.band_num for a in acqs]

        # Retrieve the satellite and sensor for the acquisition
        satellite = acqs[0].spacecraft_id
        sensor = acqs[0].sensor_id

        # Get the required nbar bands list for processing
        nbar_constants = gaip.constants.NBARConstants(satellite, sensor)
        bands_to_process = nbar_constants.get_nbar_lut()

        # Initialise the dict to store the locations of the bilinear outputs
        bilinear_fnames = {}

        for band in bands:
            if band not in bands_to_process:
                continue
            for factor in factors:
                fname = output_format.format(factor=factor, band=band)
                fname = pjoin(workpath, fname)
                bilinear_fnames[(band, factor)] = fname

        out_fname = pjoin(out_path,
                          CONFIG.get('work', 'bilinear_outputs_fname'))
        save(luigi.LocalTarget(out_fname), bilinear_fnames)

        save(self.output(), 'completed')


class DEMExctraction(luigi.Task):

    """
    Extract the DEM covering the acquisition extents plus an
    arbitrary buffer. The subset is then smoothed with a gaussian
    filter.
    """

    level1 = luigi.Parameter()
    nbar_root = luigi.Parameter()
    granule = luigi.Parameter()
    group = luigi.Parameter()

    def requires(self):
        return [CreateWorkingDirectoryTree(self.level1, self.nbar_root)]

    def output(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, group=self.group,
                                      granule=self.granule)
        out_path = pjoin(out_path, CONFIG.get('work', 'targets_root'))
        target = pjoin(out_path, 'DEMExctraction.task')
        return luigi.LocalTarget(target)

    def run(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, group=self.group,
                                      granule=self.granule)
        acqs = container.get_acquisitions(group=self.group,
                                          granule=self.granule)

        work_path = pjoin(out_path, CONFIG.get('work', 'tc_root'))
        national_dsm = CONFIG.get('ancillary', 'dem_tc')
        subset_fname = CONFIG.get('extract_dsm', 'dsm_subset')
        smoothed_fname = CONFIG.get('extract_dsm', 'dsm_smooth_subset')
        buffer = int(CONFIG.get('extract_dsm', 'dsm_buffer_width'))
        dsm_subset_fname = pjoin(work_path, subset_fname)
        dsm_subset_smooth_fname = pjoin(work_path, smoothed_fname)

        gaip.get_dsm(acqs[0], national_dsm, buffer, dsm_subset_fname,
                     dsm_subset_smooth_fname)

        save(self.output(), 'completed')


class SlopeAndAspect(luigi.Task):

    """
    Compute the slope and aspect images.
    """

    level1 = luigi.Parameter()
    nbar_root = luigi.Parameter()
    granule = luigi.Parameter()
    group = luigi.Parameter()

    def requires(self):
        args = [self.level1, self.nbar_root, self.granule, self.group]
        return [DEMExctraction(*args)]
        # return [DEMExctraction(self.level1, self.nbar_root, self.group)]

    def output(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, group=self.group,
                                      granule=self.granule)
        out_path = pjoin(out_path, CONFIG.get('work', 'targets_root'))
        target = pjoin(out_path, 'SlopeAndAspect.task')
        return luigi.LocalTarget(target)

    def run(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, group=self.group,
                                      granule=self.granule)
        work_path = pjoin(out_path, CONFIG.get('work', 'tc_root'))

        acqs = container.get_acquisitions(group=self.group,
                                          granule=self.granule)

        # Input filename
        smoothed_dsm_fname = pjoin(work_path, CONFIG.get('extract_dsm',
                                                         'dsm_smooth_subset'))
        margins = int(CONFIG.get('extract_dsm', 'dsm_buffer_width'))

        # Output filenames
        slope_fname = pjoin(work_path,
                            CONFIG.get('self_shadow', 'slope_fname'))
        aspect_fname = pjoin(work_path,
                             CONFIG.get('self_shadow', 'aspect_fname'))
        header_slope_fname = pjoin(work_path,
                                   CONFIG.get('work', 'header_slope_fname'))

        gaip.slope_aspect_arrays(acqs[0], smoothed_dsm_fname, margins,
                                 slope_fname, aspect_fname,
                                 header_slope_fname)

        save(self.output(), 'completed')


class IncidentAngles(luigi.Task):

    """
    Compute the incident angles.
    """

    level1 = luigi.Parameter()
    nbar_root = luigi.Parameter()
    granule = luigi.Parameter()
    group = luigi.Parameter()

    def requires(self):
        args = [self.level1, self.nbar_root, self.granule, self.group]
        return [CalculateSatelliteAndSolarGrids(*args),
                SlopeAndAspect(*args)]

    def output(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, group=self.group,
                                      granule=self.granule)
        out_path = pjoin(out_path, CONFIG.get('work', 'targets_root'))
        target = pjoin(out_path, 'IncidentAngles.task')
        return luigi.LocalTarget(target)

    def run(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, group=self.group,
                                      granule=self.granule)
        work_path = pjoin(out_path, CONFIG.get('work', 'tc_root'))

        # Input filenames
        solar_zenith_fname = pjoin(out_path,
                                   CONFIG.get('work', 'solar_zenith_fname'))
        solar_azimuth_fname = pjoin(out_path,
                                    CONFIG.get('work',
                                               'solar_azimuth_fname'))
        slope_fname = pjoin(work_path,
                            CONFIG.get('self_shadow', 'slope_fname'))
        aspect_fname = pjoin(work_path,
                             CONFIG.get('self_shadow', 'aspect_fname'))

        # Get the processing tile sizes
        x_tile = int(CONFIG.get('work', 'x_tile_size'))
        y_tile = int(CONFIG.get('work', 'y_tile_size'))
        x_tile = None if x_tile <= 0 else x_tile
        y_tile = None if y_tile <= 0 else y_tile

        # Output filenames
        incident_fname = pjoin(work_path,
                               CONFIG.get('self_shadow', 'incident_fname'))
        azi_incident_fname = pjoin(work_path,
                                   CONFIG.get('self_shadow',
                                              'azimuth_incident_fname'))

        gaip.incident_angles(solar_zenith_fname, solar_azimuth_fname,
                             slope_fname, aspect_fname,
                             incident_fname, azi_incident_fname,
                             x_tile, y_tile)

        save(self.output(), 'completed')


class ExitingAngles(luigi.Task):

    """
    Compute the exiting angles.
    """

    level1 = luigi.Parameter()
    nbar_root = luigi.Parameter()
    granule = luigi.Parameter()
    group = luigi.Parameter()

    def requires(self):
        args = [self.level1, self.nbar_root, self.granule, self.group]
        return [CalculateSatelliteAndSolarGrids(*args),
                SlopeAndAspect(*args)]

    def output(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, group=self.group,
                                      granule=self.granule)
        out_path = pjoin(out_path, CONFIG.get('work', 'targets_root'))
        target = pjoin(out_path, 'ExitingAngles.task')
        return luigi.LocalTarget(target)

    def run(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, group=self.group,
                                      granule=self.granule)
        work_path = pjoin(out_path, CONFIG.get('work', 'tc_root'))

        # Input filenames
        satellite_view_fname = pjoin(out_path,
                                     CONFIG.get('work', 'sat_view_fname'))
        satellite_azimuth_fname = pjoin(out_path,
                                        CONFIG.get('work',
                                                   'sat_azimuth_fname'))
        slope_fname = pjoin(work_path,
                            CONFIG.get('self_shadow', 'slope_fname'))
        aspect_fname = pjoin(work_path,
                             CONFIG.get('self_shadow', 'aspect_fname'))

        # Get the processing tile sizes
        x_tile = int(CONFIG.get('work', 'x_tile_size'))
        y_tile = int(CONFIG.get('work', 'y_tile_size'))
        x_tile = None if x_tile <= 0 else x_tile
        y_tile = None if y_tile <= 0 else y_tile

        # Output filenames
        exiting_fname = pjoin(work_path,
                              CONFIG.get('self_shadow',
                                         'exiting_fname'))
        azi_exiting_fname = pjoin(work_path,
                                  CONFIG.get('self_shadow',
                                             'azimuth_exiting_fname'))

        gaip.exiting_angles(satellite_view_fname, satellite_azimuth_fname,
                            slope_fname, aspect_fname,
                            exiting_fname, azi_exiting_fname, x_tile, y_tile)

        save(self.output(), 'completed')


class RelativeAzimuthSlope(luigi.Task):

    """
    Compute the relative azimuth angle on the slope surface.
    """

    level1 = luigi.Parameter()
    nbar_root = luigi.Parameter()
    granule = luigi.Parameter()
    group = luigi.Parameter()

    def requires(self):
        args = [self.level1, self.nbar_root, self.granule, self.group]
        return [IncidentAngles(*args),
                ExitingAngles(*args)]

    def output(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, group=self.group,
                                      granule=self.granule)
        out_path = pjoin(out_path, CONFIG.get('work', 'targets_root'))
        target = pjoin(out_path, 'RelativeAzimuthSlope.task')
        return luigi.LocalTarget(target)

    def run(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, group=self.group,
                                      granule=self.granule)
        work_path = pjoin(out_path, CONFIG.get('work', 'tc_root'))

        # Input filenames
        azi_incident_fname = pjoin(work_path,
                                   CONFIG.get('self_shadow',
                                              'azimuth_incident_fname'))
        azi_exiting_fname = pjoin(work_path,
                                  CONFIG.get('self_shadow',
                                             'azimuth_exiting_fname'))

        # Get the processing tile sizes
        x_tile = int(CONFIG.get('work', 'x_tile_size'))
        y_tile = int(CONFIG.get('work', 'y_tile_size'))
        x_tile = None if x_tile <= 0 else x_tile
        y_tile = None if y_tile <= 0 else y_tile

        # Output filename
        out_fname = pjoin(work_path, CONFIG.get('self_shadow',
                                                'relative_slope_fname'))

        gaip.relative_azimuth_slope(azi_incident_fname, azi_exiting_fname,
                                    out_fname, x_tile, y_tile)

        save(self.output(), 'completed')


class SelfShadow(luigi.Task):

    """
    Calculate the self shadow mask.
    """

    level1 = luigi.Parameter()
    nbar_root = luigi.Parameter()
    granule = luigi.Parameter()
    group = luigi.Parameter()

    def requires(self):
        args = [self.level1, self.nbar_root, self.granule, self.group]
        return [IncidentAngles(*args),
                ExitingAngles(*args)]

    def output(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, group=self.group,
                                      granule=self.granule)
        out_path = pjoin(out_path, CONFIG.get('work', 'targets_root'))
        target = pjoin(out_path, 'SelfShadow.task')
        return luigi.LocalTarget(target)

    def run(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, group=self.group,
                                      granule=self.granule)
        work_path = pjoin(out_path, CONFIG.get('work', 'tc_root'))

        # Input filenames
        incident_fname = pjoin(work_path,
                               CONFIG.get('self_shadow', 'incident_fname'))
        exiting_fname = pjoin(work_path,
                              CONFIG.get('self_shadow', 'exiting_fname'))

        # Get the processing tile sizes
        x_tile = int(CONFIG.get('work', 'x_tile_size'))
        y_tile = int(CONFIG.get('work', 'y_tile_size'))
        x_tile = None if x_tile <= 0 else x_tile
        y_tile = None if y_tile <= 0 else y_tile

        # Output filename
        out_fname = pjoin(work_path, CONFIG.get('self_shadow',
                                                'self_shadow_fname'))

        gaip.self_shadow(incident_fname, exiting_fname, out_fname,
                         x_tile, y_tile)

        save(self.output(), 'completed')


class CalculateCastShadow(luigi.Task):

    """Calculate cast shadow masks. This is a helper task."""

    level1 = luigi.Parameter()
    nbar_root = luigi.Parameter()
    granule = luigi.Parameter()
    group = luigi.Parameter()

    def requires(self):
        args = [self.level1, self.nbar_root, self.granule, self.group]
        return [CalculateCastShadowSun(*args),
                CalculateCastShadowSatellite(*args)]

    def complete(self):
        return all([t.complete() for t in self.requires()])


class CalculateCastShadowSun(luigi.Task):

    """
    Calculates the Cast shadow mask in the direction back to the
    sun.
    """

    level1 = luigi.Parameter()
    nbar_root = luigi.Parameter()
    granule = luigi.Parameter()
    group = luigi.Parameter()

    def requires(self):
        args = [self.level1, self.nbar_root, self.granule, self.group]
        return [CalculateSatelliteAndSolarGrids(*args),
                DEMExctraction(*args)]

    def output(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, group=self.group,
                                      granule=self.granule)
        out_path = pjoin(out_path, CONFIG.get('work', 'targets_root'))
        target = pjoin(out_path, 'CalculateCastShadowSun.task')
        return luigi.LocalTarget(target)

    def run(self):
        container = gaip.acquisitions(self.level1)
        acqs = container.get_acquisitions(group=self.group,
                                          granule=self.granule)
        out_path = container.get_root(self.nbar_root, group=self.group,
                                      granule=self.granule)
        tc_work_path = pjoin(out_path, CONFIG.get('work', 'tc_root'))

        # Input filenames
        smoothed_dsm_fname = pjoin(tc_work_path,
                                   CONFIG.get('extract_dsm',
                                              'dsm_smooth_subset'))
        solar_zenith_fname = pjoin(out_path,
                                   CONFIG.get('work', 'solar_zenith_fname'))
        solar_azimuth_fname = pjoin(out_path,
                                    CONFIG.get('work', 'solar_azimuth_fname'))
        buffer = int(CONFIG.get('extract_dsm', 'dsm_buffer_width'))
        window_height = int(CONFIG.get('terrain_correction',
                                       'shadow_sub_matrix_height'))
        window_width = int(CONFIG.get('terrain_correction',
                                      'shadow_sub_matrix_width'))

        # Output filename
        out_fname = pjoin(tc_work_path,
                          CONFIG.get('cast_shadow', 'sun_direction_fname'))

        gaip.calculate_cast_shadow(acqs[0], smoothed_dsm_fname, buffer,
                                   window_height, window_width,
                                   solar_zenith_fname, solar_azimuth_fname,
                                   out_fname)

        save(self.output(), 'completed')


class CalculateCastShadowSatellite(luigi.Task):

    """
    Calculates the Cast shadow mask in the direction back to the
    sun.
    """

    level1 = luigi.Parameter()
    nbar_root = luigi.Parameter()
    granule = luigi.Parameter()
    group = luigi.Parameter()

    def requires(self):
        args = [self.level1, self.nbar_root, self.granule, self.group]
        return [CalculateSatelliteAndSolarGrids(*args),
                DEMExctraction(*args)]

    def output(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, group=self.group,
                                      granule=self.granule)
        out_path = pjoin(out_path, CONFIG.get('work', 'targets_root'))
        target = pjoin(out_path, 'CalculateCastShadowSatellite.task')
        return luigi.LocalTarget(target)

    def run(self):
        container = gaip.acquisitions(self.level1)
        acqs = container.get_acquisitions(group=self.group,
                                          granule=self.granule)
        out_path = container.get_root(self.nbar_root, group=self.group,
                                      granule=self.granule)
        tc_work_path = pjoin(out_path, CONFIG.get('work', 'tc_root'))

        # Input filenames
        smoothed_dsm_fname = pjoin(tc_work_path,
                                   CONFIG.get('extract_dsm',
                                              'dsm_smooth_subset'))
        satellite_view_fname = pjoin(out_path,
                                     CONFIG.get('work', 'sat_view_fname'))
        satellite_azimuth_fname = pjoin(out_path,
                                        CONFIG.get('work',
                                                   'sat_azimuth_fname'))
        buffer = int(CONFIG.get('extract_dsm', 'dsm_buffer_width'))
        window_height = int(CONFIG.get('terrain_correction',
                                       'shadow_sub_matrix_height'))
        window_width = int(CONFIG.get('terrain_correction',
                                      'shadow_sub_matrix_width'))

        # Output fnames
        satellite_fname = pjoin(tc_work_path,
                                CONFIG.get('cast_shadow',
                                           'satellite_direction_fname'))

        gaip.calculate_cast_shadow(acqs[0], smoothed_dsm_fname, buffer,
                                   window_height, window_width,
                                   satellite_view_fname,
                                   satellite_azimuth_fname, satellite_fname)

        save(self.output(), 'completed')


class RunTCBand(luigi.Task):

    """Run the terrain correction over a given band."""

    level1 = luigi.Parameter()
    nbar_root = luigi.Parameter()
    granule = luigi.Parameter()
    group = luigi.Parameter()
    band_num = luigi.IntParameter()

    def requires(self):
        args = [self.level1, self.nbar_root, self.granule, self.group]
        return [BilinearInterpolation(*args),
                DEMExctraction(*args),
                RelativeAzimuthSlope(*args),
                SelfShadow(*args),
                CalculateCastShadow(*args)]

    def output(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, group=self.group,
                                      granule=self.granule)
        out_path = pjoin(out_path, CONFIG.get('work', 'targets_root'))
        target = pjoin(out_path, 'RunTCBand{band}.task')
        return luigi.LocalTarget(target.format(band=self.band_num))

    def run(self):
        container = gaip.acquisitions(self.level1)
        acqs = container.get_acquisitions(group=self.group,
                                          granule=self.granule)
        out_path = container.get_root(self.nbar_root, group=self.group,
                                      granule=self.granule)
        granule_path = container.get_root(self.nbar_root,
                                          granule=self.granule)

        # Get the necessary config params
        tc_path = pjoin(out_path, CONFIG.get('work', 'tc_root'))
        outdir = pjoin(out_path, CONFIG.get('work', 'reflectance_root'))
        bilinear_fname = pjoin(out_path, 
                               CONFIG.get('work', 'bilinear_outputs_fname'))
        bilinear_data = load_value(bilinear_fname)
        rori = float(CONFIG.get('terrain_correction', 'rori'))

        brdf_fname = pjoin(granule_path, CONFIG.get('work', 'brdf_fname'))
        brdf_data = load_value(brdf_fname)

        # Get the reflectance levels and base output format
        rfl_levels = CONFIG.get('terrain_correction', 'rfl_levels').split(',')
        output_format = CONFIG.get('terrain_correction', 'output_format')

        # Input filenames (images)
        self_shadow_fname = pjoin(tc_path,
                                  CONFIG.get('self_shadow',
                                             'self_shadow_fname'))
        slope_fname = pjoin(tc_path,
                            CONFIG.get('self_shadow', 'slope_fname'))
        aspect_fname = pjoin(tc_path,
                             CONFIG.get('self_shadow', 'aspect_fname'))
        incident_fname = pjoin(tc_path,
                               CONFIG.get('self_shadow', 'incident_fname'))
        exiting_fname = pjoin(tc_path,
                              CONFIG.get('self_shadow', 'exiting_fname'))
        relative_slope_fname = pjoin(tc_path,
                                     CONFIG.get('self_shadow',
                                                'relative_slope_fname'))
        sun_fname = pjoin(tc_path,
                          CONFIG.get('cast_shadow', 'sun_direction_fname'))
        satellite_fname = pjoin(tc_path,
                                CONFIG.get('cast_shadow',
                                           'satellite_direction_fname'))
        solar_zenith_fname = pjoin(out_path,
                                   CONFIG.get('work', 'solar_zenith_fname'))
        solar_azimuth_fname = pjoin(out_path,
                                    CONFIG.get('work', 'solar_azimuth_fname'))
        satellite_view_fname = pjoin(out_path,
                                     CONFIG.get('work', 'sat_view_fname'))
        relative_angle_fname = pjoin(out_path,
                                     CONFIG.get('work',
                                                'relative_azimuth_fname'))

        # Get the processing tile sizes
        x_tile = int(CONFIG.get('work', 'x_tile_size'))
        y_tile = int(CONFIG.get('work', 'y_tile_size'))
        x_tile = None if x_tile <= 0 else x_tile
        y_tile = None if y_tile <= 0 else y_tile

        # get the acquisition we wish to process
        acqs = [acq for acq in acqs if acq.band_num == self.band_num]

        # Output filenames
        # Create a dict of filenames per reflectance level per band
        rfl_lvl_fnames = {}
        for level in rfl_levels:
            outfname = output_format.format(level=level, band=self.band_num)
            rfl_lvl_fnames[(self.band_num, level)] = pjoin(outdir, outfname)

        # calculate reflectance for lambertian, brdf, and terrain correction 
        gaip.calculate_reflectance(acqs, bilinear_data, rori, brdf_data,
                                   self_shadow_fname, sun_fname,
                                   satellite_fname, solar_zenith_fname,
                                   solar_azimuth_fname, satellite_view_fname,
                                   relative_angle_fname, slope_fname,
                                   aspect_fname, incident_fname,
                                   exiting_fname, relative_slope_fname,
                                   rfl_lvl_fnames, x_tile, y_tile)

        save(self.output(), 'completed')


class TerrainCorrection(luigi.Task):

    """Perform the terrain correction."""

    level1 = luigi.Parameter()
    nbar_root = luigi.Parameter()
    granule = luigi.Parameter()
    group = luigi.Parameter()

    def requires(self):
        container = gaip.acquisitions(self.level1)
        acqs = container.get_acquisitions(group=self.group,
                                          granule=self.granule)

        # Retrieve the satellite and sensor for the acquisition
        satellite = acqs[0].spacecraft_id
        sensor = acqs[0].sensor_id
        bands = [acq.band_num for acq in acqs]

        # Get the required nbar bands list for processing
        nbar_constants = gaip.constants.NBARConstants(satellite, sensor)
        avail_bands = nbar_constants.get_nbar_lut()
        bands_to_process = [bn for bn in bands if bn in avail_bands]

        # define the bands to compute reflectance for
        tc_bands = []
        for band in bands_to_process:
            tc_bands.append(RunTCBand(self.level1, self.nbar_root,
                                      self.granule, self.group, band))

        return tc_bands

    def complete(self):
        return all([t.complete() for t in self.requires()])


class WriteMetadata(luigi.Task):

    """Write metadata."""

    level1 = luigi.Parameter()
    nbar_root = luigi.Parameter()

    def requires(self):
        container = gaip.acquisitions(self.level1)
        tasks = []
        for granule in container.granules:
            for group in container.groups:
                tasks.append(TerrainCorrection(self.level1, self.nbar_root,
                                               granule, group))
        return tasks

    def output(self):
        out_path = self.nbar_root
        out_path = pjoin(out_path, CONFIG.get('work', 'targets_root'))
        target = pjoin(out_path, 'WriteMetadata.task')
        return luigi.LocalTarget(target)

    def run(self):
        # TODO: retrieve single acquisition
        # TODO: a single acquisition doesn't account for all granules
        # TODO: do we write metadata per granule???
        container = gaip.acquisitions(self.level1)
        acq = container.get_acquisitions()[0]
        out_path = self.nbar_root

        source_info = {}
        source_info['source_scene'] = self.level1
        source_info['scene_centre_datetime'] = acq.scene_centre_datetime
        source_info['platform'] = acq.spacecraft_id
        source_info['sensor'] = acq.sensor_id
        source_info['path'] = acq.path
        source_info['row'] = acq.row

        # ancillary metadata tracking
        md = gaip.extract_ancillary_metadata(self.level1)
        for key in md:
            source_info[key] = md[key]

        targets = [pjoin(out_path, CONFIG.get('work', 'aerosol_fname')),
                   pjoin(out_path, CONFIG.get('work', 'vapour_fname')),
                   pjoin(out_path, CONFIG.get('work', 'ozone_fname')),
                   pjoin(out_path, CONFIG.get('work', 'dem_fname'))]

        sources = ['aerosol',
                   'water_vapour',
                   'ozone',
                   'elevation']

        sources_targets = zip(sources, targets)

        ancillary = {}
        for source, target in sources_targets:
            with open(target, 'rb') as src:
                ancillary[source] = pickle.load(src)

        # Get the required BRDF LUT & factors list
        nbar_constants = gaip.constants.NBARConstants(acq.spacecraft_id,
                                                      acq.sensor_id)

        bands = nbar_constants.get_brdf_lut()
        brdf_factors = nbar_constants.get_brdf_factors()

        target = pjoin(out_path, CONFIG.get('work', 'brdf_fname'))
        with open(target, 'rb') as src:
            brdf_data = pickle.load(src)

        brdf = {}
        band_fmt = 'band_{}'
        for band in bands:
            data = {}
            for factor in brdf_factors:
                data[factor] = brdf_data[(band, factor)]
            brdf[band_fmt.format(band)] = data

        ancillary['brdf'] = brdf

        # TODO (a) retrieve software version from git once deployed
        algorithm = {}
        algorithm['algorithm_version'] = 2.0 # hardcode for now see TODO (a)
        algorithm['software_repository'] = ('https://github.com/'
                                            'GeoscienceAustralia/'
                                            'ga-neo-landsat-processor.git')
        algorithm['arg25_doi'] = 'http://dx.doi.org/10.4225/25/5487CC0D4F40B'
        algorithm['nbar_doi'] = 'http://dx.doi.org/10.1109/JSTARS.2010.2042281'
        algorithm['nbar_terrain_corrected_doi'] = ('http://dx.doi.org/10.1016/'
                                                   'j.rse.2012.06.018')

        system_info = {}
        proc = subprocess.Popen(['uname', '-a'], stdout=subprocess.PIPE)
        system_info['node'] = proc.stdout.read()
        system_info['time_processed'] = dt.utcnow()

        metadata = {}
        metadata['system_information'] = system_info
        metadata['source_data'] = source_info
        metadata['ancillary_data'] = ancillary
        metadata['algorithm_information'] = algorithm
        
        # Account for NumPy dtypes
        yaml.add_representer(numpy.float, Representer.represent_float)
        yaml.add_representer(numpy.float32, Representer.represent_float)
        yaml.add_representer(numpy.float64, Representer.represent_float)

        # output
        out_fname = pjoin(out_path, CONFIG.get('work', 'metadata_fname'))
        with open(out_fname, 'w') as src:
            yaml.dump(metadata, src, default_flow_style=False)

        save(self.output(), 'completed')


class Packager(luigi.Task):

    """Packages an nbar or nbart product."""

    level1 = luigi.Parameter()
    work_root = luigi.Parameter()
    nbar_root = luigi.Parameter()
    product = luigi.Parameter()

    def requires(self):
        return [WriteMetadata(self.level1, self.nbar_root)]

    def output(self):
        out_path = pjoin(self.nbar_root, CONFIG.get('work', 'targets_root'))
        target = pjoin(out_path, 'Packager_{}.task')
        return luigi.LocalTarget(target.format(self.product))

    def run(self):
        # run the packager
        kwargs = {'driver': PACKAGE_DRIVERS[self.product],
                  'input_data_paths': [Path(self.nbar_root)],
                  'destination_path': Path(pjoin(self.work_root, self.product)),
                  'parent_dataset_paths': [Path(self.level1)],
                  'metadata_expand_fn': lambda dataset: dataset.lineage.machine.note_current_system_software(),
                  'hard_link': False}
        package_newly_processed_data_folder(**kwargs)

        save(self.output(), 'completed')


class PackageTC(luigi.Task):

    """Issues nbar &/or nbart packaging depending on the config."""

    level1 = luigi.Parameter()
    work_root = luigi.Parameter()
    nbar_root = luigi.Parameter()

    def requires(self):
        products = CONFIG.get('packaging', 'products').split(',')
        tasks = []
        for product in products:
            tasks.append(Packager(self.level1, self.work_root, self.nbar_root,
                                  product))
        return tasks

    def output(self):
        out_format = '{}.completed'
        base_dir = dirname(self.nbar_root)
        fname = out_format.format(basename(self.nbar_root))
        out_fname = pjoin(base_dir, fname)
        return luigi.LocalTarget(out_fname)

    def run(self):
        with self.output().open('w') as src:
            src.write('Task completed')

        # cleanup the entire nbar scene working directory
        cleanup = CONFIG.getboolean('cleanup', 'cleanup')
        if cleanup:
            shutil.rmtree(self.nbar_root)


def is_valid_path(parser, arg):
    """Used by argparse"""
    if not exists(arg):
        parser.error("{path} does not exist".format(path=arg))
    else:
        return arg


def scatter(iterable, P=1, p=1):
    """
    Scatter an iterator across `P` processors where `p` is the index
    of the current processor. This partitions the work evenly across
    processors.
    """
    import itertools
    return itertools.islice(iterable, p-1, None, P)
 
 
def main(level1_list, work_root, nnodes=1, nodenum=1):
    # Setup Software Versions for Packaging
    ptype.register_software_version(
        software_code='gaip',
        version=gaip.get_version(),
        repo_url='https://github.com/GeoscienceAustralia/ga-neo-landsat-processor.git'
    )
    ptype.register_software_version(
        software_code='modtran',
        version=CONFIG.get('modtran', 'version'),
        repo_url='http://www.ontar.com/software/productdetails.aspx?item=modtran'
    )

    # create product output dirs
    products = CONFIG.get('packaging', 'products').split(',')
    for product in products:
        product_dir = pjoin(work_root, product)
        if not exists(product_dir): os.makedirs(product_dir)

    tasks = []
    for level1 in open(level1_list).readlines():
        level1 = level1.strip()
        if level1 == '': continue
        bf = basename(level1)

        # create a workpath for the given level1 dataset
        nbar_root = pjoin(work_root, bf + ".nbar-work")
        tasks.append(PackageTC(level1, work_root, nbar_root))

    tasks = [f for f in scatter(tasks, nnodes, nodenum)]
    ncpus = int(os.getenv('PBS_NCPUS', '1'))
    luigi.build(tasks, local_scheduler=True, workers=ncpus / nnodes)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--level1_list",
                        help="A full file path listing of level1 datasets",
                        required=True,
                        type=lambda x: is_valid_path(parser, x))
    parser.add_argument("--work_root", help=("Path to directory where NBAR "
                        "datasets are to be written"), required=True)
    parser.add_argument('--cfg_file',
                        help='Path to a user defined configuration file.')
    parser.add_argument("--log_path", help=("Path to directory where log"
                        " files will be written"), default='.',
                        type=lambda x: is_valid_path(parser, x))
    parser.add_argument("--debug", help=("Selects more detail logging (default"
                        " is INFO)"), default=False, action='store_true')

    args = parser.parse_args()

    # TODO: save the cfg file in a HDF5(???) dataset somewhere
    cfg = args.cfg_file

    # Setup the config file
    global CONFIG
    if cfg is None:
        CONFIG = luigi.configuration.get_config()
        CONFIG.add_config_path(pjoin(dirname(__file__), 'nbar.cfg'))
    else:
        CONFIG = luigi.configuration.get_config()
        CONFIG.add_config_path(cfg)


    # setup logging
    logfile = "{log_path}/run_nbar_{uname}_{pid}.log"
    logfile = logfile.format(log_path=args.log_path, uname=os.uname()[1],
                             pid=os.getpid())
    logging_level = logging.INFO
    if args.debug:
        logging_level = logging.DEBUG
    logging.basicConfig(filename=logfile, level=logging_level,
                        format=("%(asctime)s: [%(name)s] (%(levelname)s) "
                                "%(message)s "), datefmt='%H:%M:%S')

    # TODO: save the input list into a HDF5(???) dataset somewhere
    logging.info("nbar.py started")
    logging.info('out_path={path}'.format(path=args.out_path))
    logging.info('log_path={path}'.format(path=args.log_path))

    size = int(os.getenv('PBS_NNODES', '1'))
    rank = int(os.getenv('PBS_VNODENUM', '1'))
    main(args.level1_list, args.work_root, size, rank)
