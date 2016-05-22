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
import glob
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


class GetElevationAncillaryData(luigi.Task):

    """Get ancillary elevation data."""

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        out_path = self.out_path
        target = pjoin(out_path, CONFIG.get('work', 'dem_target'))
        return luigi.LocalTarget(target)

    def run(self):
        acqs = gaip.acquisitions(self.l1t_path)
        geobox = acqs[0].gridded_geo_box()
        dem_path = CONFIG.get('ancillary', 'dem_path')
        value = gaip.get_elevation_data(geobox.centre_lonlat, dem_path)
        save(self.output(), value)


class GetOzoneAncillaryData(luigi.Task):

    """Get ancillary ozone data."""

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        out_path = self.out_path
        target = pjoin(out_path, CONFIG.get('work', 'ozone_target'))
        return luigi.LocalTarget(target)

    def run(self):
        acqs = gaip.acquisitions(self.l1t_path)
        geobox = acqs[0].gridded_geo_box()
        ozone_path = CONFIG.get('ancillary', 'ozone_path')
        centre = geobox.centre_lonlat
        dt = acqs[0].scene_center_datetime
        value = gaip.get_ozone_data(ozone_path, centre, dt)
        save(self.output(), value)


class GetSolarIrradianceAncillaryData(luigi.Task):

    """Get ancillary solar irradiance data."""

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        out_path = self.out_path
        target = pjoin(out_path, CONFIG.get('work', 'irrad_target'))
        return luigi.LocalTarget(target)

    def run(self):
        acqs = gaip.acquisitions(self.l1t_path)
        solar_path = CONFIG.get('ancillary', 'solarirrad_path')
        value = gaip.get_solar_irrad(acqs, solar_path)
        save(self.output(), value)


class GetSolarDistanceAncillaryData(luigi.Task):

    """Get ancillary solar distance data."""

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        out_path = self.out_path
        target = pjoin(out_path, CONFIG.get('work', 'sundist_target'))
        return luigi.LocalTarget(target)

    def run(self):
        acqs = gaip.acquisitions(self.l1t_path)
        sundist_path = CONFIG.get('ancillary', 'sundist_path')
        value = gaip.get_solar_dist(acqs[0], sundist_path)
        save(self.output(), value)


class GetWaterVapourAncillaryData(luigi.Task):

    """Get ancillary water vapour data."""

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        out_path = self.out_path
        target = pjoin(out_path, CONFIG.get('work', 'vapour_target'))
        return luigi.LocalTarget(target)

    def run(self):
        acqs = gaip.acquisitions(self.l1t_path)
        vapour_path = CONFIG.get('ancillary', 'vapour_path')
        value = gaip.get_water_vapour(acqs[0], vapour_path)
        save(self.output(), value)


class GetAerosolAncillaryData(luigi.Task):

    """Get ancillary aerosol data."""

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        out_path = self.out_path
        target = pjoin(out_path, CONFIG.get('work', 'aerosol_target'))
        return luigi.LocalTarget(target)

    def run(self):
        acqs = gaip.acquisitions(self.l1t_path)
        aerosol_path = CONFIG.get('ancillary', 'aerosol_path')
        value = gaip.get_aerosol_data(acqs[0], aerosol_path)
        # aerosol_path = CONFIG.get('ancillary', 'aerosol_fname') # version 2
        # value = gaip.get_aerosol_data_v2(acqs[0], aerosol_path) # version 2
        save(self.output(), value)


class GetBrdfAncillaryData(luigi.Task):

    """Get ancillary BRDF data."""

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        out_path = self.out_path
        target = pjoin(out_path, CONFIG.get('work', 'brdf_target'))
        return luigi.LocalTarget(target)

    def run(self):
        acqs = gaip.acquisitions(self.l1t_path)
        out_path = self.out_path
        brdf_path = CONFIG.get('ancillary', 'brdf_path')
        brdf_premodis_path = CONFIG.get('ancillary', 'brdf_premodis_path')
        value = gaip.get_brdf_data(acqs[0], brdf_path, brdf_premodis_path,
                                   out_path)
        save(self.output(), value)


class GetAncillaryData(luigi.Task):

    """Get all ancillary data. This a helper task."""

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return [GetElevationAncillaryData(self.l1t_path, self.out_path),
                GetOzoneAncillaryData(self.l1t_path, self.out_path),
                GetSolarDistanceAncillaryData(self.l1t_path, self.out_path),
                GetSolarIrradianceAncillaryData(self.l1t_path, self.out_path),
                GetWaterVapourAncillaryData(self.l1t_path, self.out_path),
                GetAerosolAncillaryData(self.l1t_path, self.out_path),
                GetBrdfAncillaryData(self.l1t_path, self.out_path)]

    def complete(self):
        return all([t.complete() for t in self.requires()])


class CalculateLonGrid(luigi.Task):

    """Calculate the longitude grid."""

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        out_path = self.out_path
        target = pjoin(out_path, CONFIG.get('work', 'lon_grid_target'))
        return luigi.LocalTarget(target)

    def run(self):
        acqs = gaip.acquisitions(self.l1t_path)
        target = self.output()
        gaip.create_lon_grid(acqs[0], target.fn)


class CalculateLatGrid(luigi.Task):

    """Calculate the latitude grid."""

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        out_path = self.out_path
        target = pjoin(out_path, CONFIG.get('work', 'lat_grid_target'))
        return luigi.LocalTarget(target)

    def run(self):
        acqs = gaip.acquisitions(self.l1t_path)
        target = self.output()
        gaip.create_lat_grid(acqs[0], target.fn)


class CalculateLatLonGrids(luigi.Task):

    """Calculate the longitude and latitude grids. This is a helper task."""

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return [CalculateLatGrid(self.l1t_path, self.out_path),
                CalculateLonGrid(self.l1t_path, self.out_path)]


class CalculateSatelliteAndSolarGrids(luigi.Task):

    """Calculate the satellite and solar grids."""

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return [CalculateLatGrid(self.l1t_path, self.out_path),
                CalculateLonGrid(self.l1t_path, self.out_path)]

    def output(self):
        out_path = self.out_path
        targets = [CONFIG.get('work', 'sat_view_target'),
                   CONFIG.get('work', 'sat_azimuth_target'),
                   CONFIG.get('work', 'solar_zenith_target'),
                   CONFIG.get('work', 'solar_azimuth_target'),
                   CONFIG.get('work', 'relative_azimuth_target'),
                   CONFIG.get('work', 'time_target'),
                   CONFIG.get('work', 'centreline_target'),
                   CONFIG.get('work', 'header_angle_target'),
                   CONFIG.get('work', 'coordinator_target'),
                   CONFIG.get('work', 'boxline_target')]
        return [luigi.LocalTarget(pjoin(out_path, t)) for t in targets]

    def run(self):
        out_path = self.out_path
        targets = [CONFIG.get('work', 'sat_view_target'),
                   CONFIG.get('work', 'sat_azimuth_target'),
                   CONFIG.get('work', 'solar_zenith_target'),
                   CONFIG.get('work', 'solar_azimuth_target'),
                   CONFIG.get('work', 'relative_azimuth_target'),
                   CONFIG.get('work', 'time_target')]
        targets = [pjoin(out_path, t) for t in targets]
        centreline_target = pjoin(out_path,
                                  CONFIG.get('work', 'centreline_target'))
        header_angle_target = pjoin(out_path,
                                    CONFIG.get('work', 'header_angle_target'))
        coordinator_target = pjoin(out_path,
                                   CONFIG.get('work', 'coordinator_target'))
        boxline_target = pjoin(out_path, CONFIG.get('work', 'boxline_target'))
        lon_target = pjoin(out_path, CONFIG.get('work', 'lon_grid_target'))
        lat_target = pjoin(out_path, CONFIG.get('work', 'lat_grid_target'))

        acqs = gaip.acquisitions(self.l1t_path)

        geobox = acqs[0].gridded_geo_box()
        cols = acqs[0].samples

        (satellite_zenith, satellite_azimuth, solar_zenith, solar_azimuth,
         relative_azimuth, time, y_cent, x_cent, n_cent) = \
            gaip.calculate_angles(acqs[0], lon_target, lat_target,
                                  npoints=12, out_fnames=targets)

        gaip.create_centreline_file(geobox, y_cent, x_cent, n_cent, cols,
                                    view_max=9.0, outfname=centreline_target)

        gaip.create_header_angle_file(acqs[0], view_max=9.0,
                                      outfname=header_angle_target)

        gaip.create_boxline_file(satellite_zenith, y_cent, x_cent,
                                 boxline_fname=boxline_target,
                                 coordinator_fname=coordinator_target)


class CalculateGridsTask(luigi.Task):

    """Calculate all the grids. This is a helper task."""

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return [CalculateLatLonGrids(self.l1t_path, self.out_path),
                CalculateSatelliteAndSolarGrids(self.l1t_path, self.out_path)]


class CreateModtranDirectories(luigi.Task):

    """Create the MODTRAN work directories and input driver files."""

    out_path = luigi.Parameter()

    def output(self):
        out_path = self.out_path
        input_format = CONFIG.get('modtran', 'input_format')
        coords = CONFIG.get('modtran', 'coords').split(',')
        albedos = CONFIG.get('modtran', 'albedos').split(',')
        modtran_root = pjoin(out_path, CONFIG.get('work', 'modtran_root'))
        targets = []
        for coord in coords:
            for albedo in albedos:
                targets.append(input_format.format(coord=coord,
                                                   albedo=albedo))
        return [luigi.LocalTarget(pjoin(modtran_root, t)) for t in targets]

    def run(self):
        out_path = self.out_path
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


class CreateSatelliteFilterFile(luigi.Task):

    """Create the satellite filter file."""

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def output(self):
        out_path = self.out_path
        target = pjoin(out_path, CONFIG.get('work', 'sat_filter_target'))
        return luigi.LocalTarget(target)

    def run(self):
        out_path = self.out_path
        acqs = gaip.acquisitions(self.l1t_path)
        satfilterpath = CONFIG.get('ancillary', 'satfilter_path')
        target = pjoin(out_path, CONFIG.get('work', 'sat_filter_target'))
        gaip.create_satellite_filter_file(acqs, satfilterpath,
                                          target)


# Keep around for testing; for the time being
class CreateModtranInputFile(luigi.Task):

    """Create the MODTRAN input file."""

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return [GetAncillaryData(self.l1t_path, self.out_path)]

    def output(self):
        out_path = self.out_path
        target = pjoin(out_path, CONFIG.get('work', 'modtran_input_target'))
        return luigi.LocalTarget(target)

    def run(self):
        out_path = self.out_path
        ozone_target = pjoin(out_path, CONFIG.get('work', 'ozone_target'))
        vapour_target = pjoin(out_path, CONFIG.get('work', 'vapour_target'))
        aerosol_target = pjoin(out_path, CONFIG.get('work', 'aerosol_target'))
        elevation_target = pjoin(out_path, CONFIG.get('work', 'dem_target'))
        acqs = gaip.acquisitions(self.l1t_path)
        target = self.output().fn
        ozone = load_value(ozone_target)
        vapour = load_value(vapour_target)
        aerosol = load_value(aerosol_target)
        elevation = load_value(elevation_target)
        gaip.write_modtran_input(acqs, target, ozone, vapour, aerosol,
                                 elevation)


class CreateModisBrdfFiles(luigi.Task):

    """Create the Modis BRDF files."""

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return [GetAncillaryData(self.l1t_path, self.out_path)]

    def output(self):
        acqs = gaip.acquisitions(self.l1t_path)
        out_path = self.out_path
        modis_brdf_format = pjoin(out_path,
                                  CONFIG.get('brdf', 'modis_brdf_format'))

        # Retrieve the satellite and sensor for the acquisition
        satellite = acqs[0].spacecraft_id
        sensor = acqs[0].sensor_id

        # Get the required nbar bands list for processing
        nbar_constants = gaip.constants.NBARConstants(satellite, sensor)
        bands_to_process = nbar_constants.get_nbar_lut()

        targets = []
        for acq in acqs:
            band = acq.band_num
            if band not in bands_to_process:
                continue
            modis_brdf_filename = modis_brdf_format.format(band_num=band)
            target = pjoin(out_path, modis_brdf_filename)
            targets.append(luigi.LocalTarget(target))
        return targets

    def run(self):
        acqs = gaip.acquisitions(self.l1t_path)
        outdir = self.out_path
        modis_brdf_format = pjoin(outdir,
                                  CONFIG.get('brdf', 'modis_brdf_format'))
        brdf_target = pjoin(outdir, CONFIG.get('work', 'brdf_target'))
        brdf_data = load_value(brdf_target)
        irrad_target = pjoin(outdir, CONFIG.get('work', 'irrad_target'))
        solar_irrad_data = load_value(irrad_target)
        solar_dist_target = pjoin(outdir, CONFIG.get('work', 'sundist_target'))
        solar_dist_data = load_value(solar_dist_target)
        gaip.write_modis_brdf_files(acqs, modis_brdf_format, brdf_data,
                                    solar_irrad_data, solar_dist_data)


# Keep around for testing; for the time being
class GenerateModtranInputFiles(luigi.Task):

    """Generate the MODTRAN input files."""

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return [CreateModtranDirectories(self.out_path),
                GetAncillaryData(self.l1t_path, self.out_path),
                CalculateSatelliteAndSolarGrids(self.l1t_path, self.out_path),
                CalculateLatGrid(self.l1t_path, self.out_path),
                CalculateLonGrid(self.l1t_path, self.out_path)]

    def output(self):
        out_path = self.out_path
        coords = CONFIG.get('input_modtran', 'coords').split(',')
        albedos = CONFIG.get('input_modtran', 'albedos').split(',')
        output_format = CONFIG.get('input_modtran', 'output_format')
        workdir = pjoin(out_path, CONFIG.get('work', 'input_modtran_cwd'))
        output_format = pjoin(workdir, output_format)

        targets = []
        for coord in coords:
            for albedo in albedos:
                targets.append(output_format.format(coord=coord,
                                                    albedo=albedo))
        return [luigi.LocalTarget(t) for t in targets]

    def run(self):
        out_path = self.out_path
        # sources
        modtran_input_target = pjoin(out_path,
                                     CONFIG.get('work',
                                                'modtran_input_target'))
        coordinator_target = pjoin(out_path,
                                   CONFIG.get('work', 'coordinator_target'))
        sat_view_zenith_target = pjoin(out_path,
                                       CONFIG.get('work', 'sat_view_target'))
        sat_azimuth_target = pjoin(out_path,
                                   CONFIG.get('work', 'sat_azimuth_target'))
        lon_grid_target = pjoin(out_path,
                                CONFIG.get('work', 'lon_grid_target'))
        lat_grid_target = pjoin(out_path,
                                CONFIG.get('work', 'lat_grid_target'))

        coords = CONFIG.get('input_modtran', 'coords').split(',')
        albedos = CONFIG.get('input_modtran', 'albedos').split(',')
        fname_format = CONFIG.get('input_modtran', 'output_format')
        workdir = pjoin(out_path, CONFIG.get('work', 'input_modtran_cwd'))

        acqs = gaip.acquisitions(self.l1t_path)
        ozone_target = pjoin(out_path, CONFIG.get('work', 'ozone_target'))
        vapour_target = pjoin(out_path, CONFIG.get('work', 'vapour_target'))
        aerosol_target = pjoin(out_path, CONFIG.get('work', 'aerosol_target'))
        elevation_target = pjoin(out_path, CONFIG.get('work', 'dem_target'))
        ozone = load_value(ozone_target)
        vapour = load_value(vapour_target)
        aerosol = load_value(aerosol_target)
        elevation = load_value(elevation_target)
        out_fname_fmt = pjoin(workdir, fname_format)
        gaip.write_modtran_inputs(acqs[0], coordinator_target,
                                  sat_view_zenith_target, sat_azimuth_target,
                                  lat_grid_target, lon_grid_target, ozone,
                                  vapour, aerosol, elevation, coords, albedos,
                                  out_fname_fmt)


class WriteTp5(luigi.Task):

    """Output the `tp5` formatted files."""

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return [CreateModtranDirectories(self.out_path),
                GetAncillaryData(self.l1t_path, self.out_path),
                CalculateSatelliteAndSolarGrids(self.l1t_path, self.out_path),
                CalculateLatGrid(self.l1t_path, self.out_path),
                CalculateLonGrid(self.l1t_path, self.out_path)]

    def output(self):
        out_path = self.out_path
        coords = CONFIG.get('write_tp5', 'coords').split(',')
        albedos = CONFIG.get('write_tp5', 'albedos').split(',')
        output_format = CONFIG.get('write_tp5', 'output_format')
        workdir = pjoin(out_path, CONFIG.get('work', 'modtran_root'))
        targets = []
        for coord in coords:
            for albedo in albedos:
                targets.append(output_format.format(coord=coord,
                                                    albedo=albedo))
        return [luigi.LocalTarget(pjoin(workdir, t)) for t in targets]

    def run(self):
        out_path = self.out_path
        coords = CONFIG.get('write_tp5', 'coords').split(',')
        albedos = CONFIG.get('write_tp5', 'albedos').split(',')
        output_format = CONFIG.get('write_tp5', 'output_format')
        workdir = pjoin(out_path, CONFIG.get('work', 'modtran_root'))
        out_fname_format = pjoin(workdir, output_format)

        # get the filenames for the coordinator,
        # satellite view zenith, azimuth, and latitude/longitude arrays
        coord_fname = pjoin(out_path, CONFIG.get('work', 'coordinator_target'))
        sat_view_fname = pjoin(out_path, CONFIG.get('work', 'sat_view_target'))
        sat_azi_fname = pjoin(out_path, CONFIG.get('work',
                                                   'sat_azimuth_target'))
        lon_fname = pjoin(out_path, CONFIG.get('work', 'lon_grid_target'))
        lat_fname = pjoin(out_path, CONFIG.get('work', 'lat_grid_target'))

        # load the ancillary point values
        ozone_fname = pjoin(out_path, CONFIG.get('work', 'ozone_target'))
        vapour_fname = pjoin(out_path, CONFIG.get('work', 'vapour_target'))
        aerosol_fname = pjoin(out_path, CONFIG.get('work', 'aerosol_target'))
        elevation_fname = pjoin(out_path, CONFIG.get('work', 'dem_target'))
        ozone = load_value(ozone_fname)
        vapour = load_value(vapour_fname)
        aerosol = load_value(aerosol_fname)
        elevation = load_value(elevation_fname)

        # load an acquisition
        acq = gaip.acquisitions(self.l1t_path)[0]

        # run
        gaip.write_tp5(acq, coord_fname, sat_view_fname, sat_azi_fname,
                       lat_fname, lon_fname, ozone, vapour, aerosol, elevation,
                       coords, albedos, out_fname_format)


class PrepareModtranInput(luigi.Task):

    """Prepare MODTRAN inputs. This is a helper task."""

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
                # retain for testing; for the time being.
                # CreateModtranInputFile(self.l1t_path, self.out_path),
                # GenerateModtranInputFiles(self.l1t_path, self.out_path),
        return [CreateModtranDirectories(self.out_path),
                CreateSatelliteFilterFile(self.l1t_path, self.out_path),
                WriteTp5(self.l1t_path, self.out_path)]

    def complete(self):
        return all([t.complete() for t in self.requires()])


class RunModtranCase(luigi.Task):

    """Run MODTRAN for a specific `coord` and `albedo`. This task is
    parameterised this way to allow parallel instances of MODTRAN to run."""

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()
    coord = luigi.Parameter()
    albedo = luigi.Parameter()

    def requires(self):
        return [PrepareModtranInput(self.l1t_path, self.out_path)]

    def output(self):
        out_path = self.out_path
        modtran_root = pjoin(out_path, CONFIG.get('work', 'modtran_root'))
        flux_format = CONFIG.get('modtran', 'flx_output_format')
        flux_format = pjoin(modtran_root, flux_format)
        coef_format = CONFIG.get('modtran', 'chn_output_format')
        coef_format = pjoin(modtran_root, coef_format)
        flx_target = flux_format.format(coord=self.coord, albedo=self.albedo)
        chn_target = coef_format.format(coord=self.coord, albedo=self.albedo)
        return [luigi.LocalTarget(flx_target),
                luigi.LocalTarget(chn_target)]

    def run(self):
        out_path = self.out_path
        modtran_exe = CONFIG.get('modtran', 'exe')
        workpath_format = CONFIG.get('modtran', 'workpath_format')
        modtran_root = pjoin(out_path, CONFIG.get('work', 'modtran_root'))
        workpath = workpath_format.format(coord=self.coord, albedo=self.albedo)
        gaip.run_modtran(modtran_exe, pjoin(modtran_root, workpath))


class RunModtran(luigi.Task):

    """Run MODTRAN for all coords and albedos. This is a helper task."""

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        coords = CONFIG.get('modtran', 'coords').split(',')
        albedos = CONFIG.get('modtran', 'albedos').split(',')
        reqs = [PrepareModtranInput(self.l1t_path, self.out_path)]
        for coord in coords:
            for albedo in albedos:
                reqs.append(RunModtranCase(self.l1t_path, self.out_path, 
                                           coord, albedo))
        return reqs

    def complete(self):
        return all([t.complete() for t in self.requires()])


class ExtractFlux(luigi.Task):

    """Extract the flux data from the MODTRAN outputs. This runs the
       Fortran binary `read_flux_albedo`."""

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return [RunModtran(self.l1t_path, self.out_path)]

    def output(self):
        out_path = self.out_path
        coords = CONFIG.get('extract_flux', 'coords').split(',')
        albedos = CONFIG.get('extract_flux', 'albedos').split(',')
        modtran_root = pjoin(out_path, CONFIG.get('work', 'modtran_root'))
        output_format = CONFIG.get('extract_flux', 'output_format')
        output_format = pjoin(modtran_root, output_format)
        targets = []
        for coord in coords:
            for albedo in albedos:
                target = output_format.format(coord=coord, albedo=albedo)
                targets.append(luigi.LocalTarget(target))
        return targets

    def run(self):
        out_path = self.out_path
        coords = CONFIG.get('extract_flux', 'coords').split(',')
        albedos = CONFIG.get('extract_flux', 'albedos').split(',')
        modtran_root = pjoin(out_path, CONFIG.get('work', 'modtran_root'))
        input_format = CONFIG.get('extract_flux', 'input_format')
        input_format = pjoin(modtran_root, input_format)
        output_format = CONFIG.get('extract_flux', 'output_format')
        output_format = pjoin(modtran_root, output_format)
        satfilter = pjoin(out_path, CONFIG.get('work', 'sat_filter_target'))

        gaip.extract_flux(coords, albedos, input_format, output_format,
                          satfilter)


class ExtractFluxTrans(luigi.Task):

    """Extract the flux data from the MODTRAN output in the transmissive
    case. This runs the Fortran binary `read_flux_transmittance`."""

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return [RunModtran(self.l1t_path, self.out_path)]

    def output(self):
        out_path = self.out_path
        coords = CONFIG.get('extract_flux_trans', 'coords').split(',')
        modtran_root = pjoin(out_path, CONFIG.get('work', 'modtran_root'))
        output_format = CONFIG.get('extract_flux_trans', 'output_format')
        output_format = pjoin(modtran_root, output_format)
        targets = []
        for coord in coords:
            target = output_format.format(coord=coord)
            targets.append(luigi.LocalTarget(target))
        return targets

    def run(self):
        out_path = self.out_path
        coords = CONFIG.get('extract_flux_trans', 'coords').split(',')
        modtran_root = pjoin(out_path, CONFIG.get('work', 'modtran_root'))
        input_format = CONFIG.get('extract_flux_trans', 'input_format')
        input_format = pjoin(modtran_root, input_format)
        output_format = CONFIG.get('extract_flux_trans', 'output_format')
        output_format = pjoin(modtran_root, output_format)
        satfilter = pjoin(out_path, CONFIG.get('work', 'sat_filter_target'))

        gaip.extract_flux_trans(coords, input_format, output_format,
                                satfilter)


class CalculateCoefficients(luigi.Task):

    """Calculate the atmospheric parameters needed by BRDF and atmospheric
    correction model. This runs the Fortran binary `calculate_coefficients`."""

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return [ExtractFlux(self.l1t_path, self.out_path),
                ExtractFluxTrans(self.l1t_path, self.out_path)]

    def output(self):
        out_path = self.out_path
        coords = CONFIG.get('coefficients', 'coords').split(',')
        modtran_root = pjoin(out_path, CONFIG.get('work', 'modtran_root'))
        output_format = CONFIG.get('coefficients', 'output_format')
        output_format = pjoin(modtran_root, output_format)
        targets = []
        for coord in coords:
            target = output_format.format(coord=coord)
            targets.append(luigi.LocalTarget(target))
        return targets

    def run(self):
        out_path = self.out_path
        coords = CONFIG.get('coefficients', 'coords').split(',')
        chn_input_format = CONFIG.get('coefficients', 'chn_input_format')
        dir_input_format = CONFIG.get('coefficients', 'dir_input_format')
        output_format = CONFIG.get('coefficients', 'output_format')
        satfilter = pjoin(out_path, CONFIG.get('work', 'sat_filter_target'))
        workpath = pjoin(out_path, CONFIG.get('work', 'modtran_root'))

        gaip.calc_coefficients(coords, chn_input_format, dir_input_format,
                               output_format, satfilter, workpath)


class ReformatAtmosphericParameters(luigi.Task):

    """Reformat the atmospheric parameters produced by MODTRAN for four boxes.
    These are needed to conduct bilinear interpolation. This runs the binary
    `reformat_modtran_output`. """

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return [CalculateCoefficients(self.l1t_path, self.out_path),
                CreateSatelliteFilterFile(self.l1t_path, self.out_path)]

    def output(self):
        out_path = self.out_path
        factors = CONFIG.get('read_modtran', 'factors').split(',')
        modtran_root = pjoin(out_path, CONFIG.get('work', 'modtran_root'))
        output_format = CONFIG.get('read_modtran', 'output_format')
        output_format = pjoin(modtran_root, output_format)
        acqs = gaip.acquisitions(self.l1t_path)

        # Retrieve the satellite and sensor for the acquisition
        satellite = acqs[0].spacecraft_id
        sensor = acqs[0].sensor_id

        # Get the required nbar bands list for processing
        nbar_constants = gaip.constants.NBARConstants(satellite, sensor)
        bands_to_process = nbar_constants.get_nbar_lut()

        bands = [a.band_num for a in acqs]
        targets = []
        for factor in factors:
            for band in bands:
                if band not in bands_to_process:
                    # Skip
                    continue
                target = output_format.format(factor=factor, band=band)
                targets.append(luigi.LocalTarget(target))
        return targets

    def run(self):
        out_path = self.out_path
        coords = CONFIG.get('read_modtran', 'coords').split(',')
        factors = CONFIG.get('read_modtran', 'factors').split(',')
        workpath = pjoin(out_path, CONFIG.get('work', 'modtran_root'))
        input_format = CONFIG.get('read_modtran', 'input_format')
        input_format = pjoin(workpath, input_format)
        output_format = CONFIG.get('read_modtran', 'output_format')
        output_format = pjoin(workpath, output_format)
        satfilter = pjoin(out_path, CONFIG.get('work', 'sat_filter_target'))

        acqs = gaip.acquisitions(self.l1t_path)

        # Retrieve the satellite and sensor for the acquisition
        satellite = acqs[0].spacecraft_id
        sensor = acqs[0].sensor_id

        # Get the required nbar bands list for processing
        nbar_constants = gaip.constants.NBARConstants(satellite, sensor)
        bands_to_process = nbar_constants.get_nbar_lut()

        # Initialise the list to contain the acquisitions we wish to process
        acqs_to_process = []
        for acq in acqs:
            band_number = acq.band_num
            if band_number in bands_to_process:
                acqs_to_process.append(acq)

        gaip.reformat_atmo_params(acqs_to_process, coords, satfilter, factors,
                                  input_format, output_format, workpath)


class BilinearInterpolationBand(luigi.Task):
    """
    Runs the bilinear interpolation function for a given band.
    """

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()
    band_num = luigi.IntParameter()
    factor = luigi.Parameter()

    def requires(self):
        return [ReformatAtmosphericParameters(self.l1t_path, self.out_path),
                CalculateSatelliteAndSolarGrids(self.l1t_path, self.out_path)]

    def output(self):
        out_path = self.out_path
        modtran_root = pjoin(out_path, CONFIG.get('work', 'modtran_root'))
        output_format = CONFIG.get('bilinear', 'output_format')
        output_format = pjoin(modtran_root, output_format)

        out_fname = output_format.format(factor=self.factor,
                                         band=self.band_num)
        return luigi.LocalTarget(out_fname)

    def run(self):
        out_path = self.out_path
        coordinator = pjoin(out_path,
                            CONFIG.get('work', 'coordinator_target'))
        boxline = pjoin(out_path,
                        CONFIG.get('work', 'boxline_target'))
        centreline = pjoin(out_path,
                           CONFIG.get('work', 'centreline_target'))
        input_format = CONFIG.get('bilinear', 'input_format')
        output_format = CONFIG.get('bilinear', 'output_format')
        workpath = pjoin(out_path,
                         CONFIG.get('work', 'modtran_root'))
        input_format = pjoin(workpath, input_format)

        acqs = gaip.acquisitions(self.l1t_path)

        # get the acquisition we wish to process
        acqs = [acq for acq in acqs if acq.band_num == self.band_num]

        # Retrieve the satellite and sensor for the acquisition
        satellite = acqs[0].spacecraft_id
        sensor = acqs[0].sensor_id

        bilinear_fnames = gaip.bilinear_interpolate(acqs, [self.factor],
                                                    coordinator, boxline,
                                                    centreline, input_format,
                                                    output_format, workpath)


class BilinearInterpolation(luigi.Task):
    """
    Issues BilinearInterpolationBand tasks.
    This is a helper task.
    """

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        out_path = self.out_path
        modtran_root = pjoin(out_path, CONFIG.get('work', 'modtran_root'))
        factors = CONFIG.get('bilinear', 'factors').split(',')
        output_format = CONFIG.get('bilinear', 'output_format')
        output_format = pjoin(modtran_root, output_format)
        acqs = gaip.acquisitions(self.l1t_path)

        # Retrieve the satellite and sensor for the acquisition
        satellite = acqs[0].spacecraft_id
        sensor = acqs[0].sensor_id

        # Get the required nbar bands list for processing
        nbar_constants = gaip.constants.NBARConstants(satellite, sensor)
        bands_to_process = nbar_constants.get_nbar_lut()

        bands = [a.band_num for a in acqs]
        tasks = []
        for factor in factors:
            for band in bands:
                if band not in bands_to_process:
                    # Skip
                    continue
                tasks.append(BilinearInterpolationBand(self.l1t_path,
                                                       self.out_path,
                                                       band, factor))
        return tasks

    def output(self):
        out_path = self.out_path
        target = pjoin(out_path,
                       CONFIG.get('work', 'bilinear_outputs_target'))
        return luigi.LocalTarget(target)

    def run(self):
        out_path = self.out_path
        factors = CONFIG.get('bilinear', 'factors').split(',')
        input_format = CONFIG.get('bilinear', 'input_format')
        output_format = CONFIG.get('bilinear', 'output_format')
        workpath = pjoin(out_path,
                         CONFIG.get('work', 'modtran_root'))
        input_format = pjoin(workpath, input_format)

        acqs = gaip.acquisitions(self.l1t_path)
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

        save(self.output(), bilinear_fnames)


class CreateTCRflDirs(luigi.Task):
    """
    Setup the directories to contain the Intermediate files
    produced for terrain corection.
    """

    out_path = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        out_path = self.out_path
        tc_path = pjoin(out_path, CONFIG.get('work', 'tc_intermediates'))
        rfl_path = pjoin(out_path, CONFIG.get('work', 'rfl_output_dir'))

        targets = [luigi.LocalTarget(tc_path), luigi.LocalTarget(rfl_path)]
        return targets

    def run(self):
        out_path = self.out_path
        tc_path = pjoin(out_path, CONFIG.get('work', 'tc_intermediates'))
        rfl_path = pjoin(out_path, CONFIG.get('work', 'rfl_output_dir'))
        if not exists(tc_path):
            os.makedirs(tc_path)
        if not exists(rfl_path):
            os.makedirs(rfl_path)

class DEMExctraction(luigi.Task):

    """
    Extract the DEM covering the acquisition extents plus an
    arbitrary buffer. The subset is then smoothed with a gaussian
    filter.
    """

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return [CreateTCRflDirs(self.out_path)]

    def output(self):
        out_path = self.out_path
        work_path = pjoin(out_path, CONFIG.get('work', 'tc_intermediates'))
        subset_target = pjoin(work_path,
                              CONFIG.get('extract_dsm', 'dsm_subset'))
        smoothed_target = pjoin(work_path,
                                CONFIG.get('extract_dsm', 'dsm_smooth_subset'))
        targets = [luigi.LocalTarget(subset_target),
                   luigi.LocalTarget(smoothed_target)]
        return targets

    def run(self):
        acqs = gaip.acquisitions(self.l1t_path)
        out_path = self.out_path
        work_path = pjoin(out_path, CONFIG.get('work', 'tc_intermediates'))
        national_dsm = CONFIG.get('ancillary', 'dem_tc')
        subset_target = CONFIG.get('extract_dsm', 'dsm_subset')
        smoothed_target = CONFIG.get('extract_dsm', 'dsm_smooth_subset')
        buffer = int(CONFIG.get('extract_dsm', 'dsm_buffer_width'))
        dsm_subset_fname = pjoin(work_path, subset_target)
        dsm_subset_smooth_fname = pjoin(work_path, smoothed_target)

        gaip.get_dsm(acqs[0], national_dsm, buffer, dsm_subset_fname,
                     dsm_subset_smooth_fname)

class SlopeAndAspect(luigi.Task):

    """
    Compute the slope and aspect images.
    """

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return [DEMExctraction(self.l1t_path, self.out_path)]

    def output(self):
        out_path = self.out_path
        work_path = pjoin(out_path, CONFIG.get('work', 'tc_intermediates'))

        slope_target = pjoin(work_path,
                             CONFIG.get('self_shadow', 'slope_target'))
        aspect_target = pjoin(work_path,
                              CONFIG.get('self_shadow', 'aspect_target'))
        header_slope_target = pjoin(work_path,
                                    CONFIG.get('work', 'header_slope_target'))

        targets = [luigi.LocalTarget(slope_target),
                   luigi.LocalTarget(aspect_target),
                   luigi.LocalTarget(header_slope_target)]

        return targets

    def run(self):
        out_path = self.out_path
        work_path = pjoin(out_path, CONFIG.get('work', 'tc_intermediates'))

        acqs = gaip.acquisitions(self.l1t_path)

        # Input targets
        smoothed_dsm_fname = pjoin(work_path, CONFIG.get('extract_dsm',
                                                         'dsm_smooth_subset'))
        margins = int(CONFIG.get('extract_dsm', 'dsm_buffer_width'))

        # Output targets
        slope_target = pjoin(work_path,
                             CONFIG.get('self_shadow', 'slope_target'))
        aspect_target = pjoin(work_path,
                              CONFIG.get('self_shadow', 'aspect_target'))
        header_slope_target = pjoin(work_path,
                                    CONFIG.get('work', 'header_slope_target'))

        gaip.slope_aspect_arrays(acqs[0], smoothed_dsm_fname, margins,
                                 slope_target, aspect_target,
                                 header_slope_target)


class IncidentAngles(luigi.Task):

    """
    Compute the incident angles.
    """

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return [CalculateSatelliteAndSolarGrids(self.l1t_path, self.out_path),
                SlopeAndAspect(self.l1t_path, self.out_path)]

    def output(self):
        out_path = self.out_path
        work_path = pjoin(out_path, CONFIG.get('work', 'tc_intermediates'))

        incident_target = pjoin(work_path,
                                CONFIG.get('self_shadow', 'incident_target'))
        azi_incident_target = pjoin(work_path,
                                    CONFIG.get('self_shadow',
                                               'azimuth_incident_target'))

        targets = [luigi.LocalTarget(incident_target),
                   luigi.LocalTarget(azi_incident_target)]

        return targets

    def run(self):
        out_path = self.out_path
        work_path = pjoin(out_path, CONFIG.get('work', 'tc_intermediates'))

        # Input targets
        solar_zenith_fname = pjoin(out_path,
                                   CONFIG.get('work', 'solar_zenith_target'))
        solar_azimuth_fname = pjoin(out_path,
                                    CONFIG.get('work',
                                               'solar_azimuth_target'))
        slope_target = pjoin(work_path,
                             CONFIG.get('self_shadow', 'slope_target'))
        aspect_target = pjoin(work_path,
                              CONFIG.get('self_shadow', 'aspect_target'))

        # Get the processing tile sizes
        x_tile = int(CONFIG.get('work', 'x_tile_size'))
        y_tile = int(CONFIG.get('work', 'y_tile_size'))
        x_tile = None if x_tile <= 0 else x_tile
        y_tile = None if y_tile <= 0 else y_tile

        # Output targets
        incident_target = pjoin(work_path,
                                CONFIG.get('self_shadow', 'incident_target'))
        azi_incident_target = pjoin(work_path,
                                    CONFIG.get('self_shadow',
                                               'azimuth_incident_target'))

        gaip.incident_angles(solar_zenith_fname, solar_azimuth_fname,
                             slope_target, aspect_target,
                             incident_target, azi_incident_target,
                             x_tile, y_tile)


class ExitingAngles(luigi.Task):

    """
    Compute the exiting angles.
    """

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return [CalculateSatelliteAndSolarGrids(self.l1t_path, self.out_path),
                SlopeAndAspect(self.l1t_path, self.out_path)]

    def output(self):
        out_path = self.out_path
        work_path = pjoin(out_path, CONFIG.get('work', 'tc_intermediates'))

        exiting_target = pjoin(work_path,
                               CONFIG.get('self_shadow',
                                          'exiting_target'))
        azi_exiting_target = pjoin(work_path,
                                   CONFIG.get('self_shadow',
                                              'azimuth_exiting_target'))

        targets = [luigi.LocalTarget(exiting_target),
                   luigi.LocalTarget(azi_exiting_target)]

        return targets

    def run(self):
        out_path = self.out_path
        work_path = pjoin(out_path, CONFIG.get('work', 'tc_intermediates'))

        # Input targets
        satellite_view_fname = pjoin(out_path,
                                     CONFIG.get('work', 'sat_view_target'))
        satellite_azimuth_fname = pjoin(out_path,
                                        CONFIG.get('work',
                                                   'sat_azimuth_target'))
        slope_target = pjoin(work_path,
                             CONFIG.get('self_shadow', 'slope_target'))
        aspect_target = pjoin(work_path,
                              CONFIG.get('self_shadow', 'aspect_target'))

        # Get the processing tile sizes
        x_tile = int(CONFIG.get('work', 'x_tile_size'))
        y_tile = int(CONFIG.get('work', 'y_tile_size'))
        x_tile = None if x_tile <= 0 else x_tile
        y_tile = None if y_tile <= 0 else y_tile

        # Output targets
        exiting_target = pjoin(work_path,
                               CONFIG.get('self_shadow',
                                          'exiting_target'))
        azi_exiting_target = pjoin(work_path,
                                   CONFIG.get('self_shadow',
                                              'azimuth_exiting_target'))

        gaip.exiting_angles(satellite_view_fname, satellite_azimuth_fname,
                            slope_target, aspect_target,
                            exiting_target, azi_exiting_target, x_tile, y_tile)


class RelativeAzimuthSlope(luigi.Task):

    """
    Compute the relative azimuth angle on the slope surface.
    """

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return [IncidentAngles(self.l1t_path, self.out_path),
                ExitingAngles(self.l1t_path, self.out_path)]

    def output(self):
        out_path = self.out_path
        work_path = pjoin(out_path, CONFIG.get('work', 'tc_intermediates'))

        relative_azimuth_slope_target = pjoin(work_path,
                                              CONFIG.get('self_shadow',
                                                      'relative_slope_target'))

        return luigi.LocalTarget(relative_azimuth_slope_target)

    def run(self):
        out_path = self.out_path
        work_path = pjoin(out_path, CONFIG.get('work', 'tc_intermediates'))

        # Input targets
        azi_incident_target = pjoin(work_path,
                                    CONFIG.get('self_shadow',
                                               'azimuth_incident_target'))
        azi_exiting_target = pjoin(work_path,
                                   CONFIG.get('self_shadow',
                                              'azimuth_exiting_target'))

        # Get the processing tile sizes
        x_tile = int(CONFIG.get('work', 'x_tile_size'))
        y_tile = int(CONFIG.get('work', 'y_tile_size'))
        x_tile = None if x_tile <= 0 else x_tile
        y_tile = None if y_tile <= 0 else y_tile

        # Output target
        relative_azimuth_slope_target = pjoin(work_path,
                                              CONFIG.get('self_shadow',
                                                      'relative_slope_target'))

        gaip.relative_azimuth_slope(azi_incident_target, azi_exiting_target,
                                    relative_azimuth_slope_target,
                                    x_tile, y_tile)

class SelfShadow(luigi.Task):

    """
    Calculate the self shadow mask.
    """

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return [IncidentAngles(self.l1t_path, self.out_path),
                ExitingAngles(self.l1t_path, self.out_path)]

    def output(self):
        out_path = self.out_path
        work_path = pjoin(out_path, CONFIG.get('work', 'tc_intermediates'))

        self_shadow_target = pjoin(work_path,
                                   CONFIG.get('self_shadow',
                                              'self_shadow_target'))

        return luigi.LocalTarget(self_shadow_target)

    def run(self):
        out_path = self.out_path
        work_path = pjoin(out_path, CONFIG.get('work', 'tc_intermediates'))

        # Input targets
        incident_target = pjoin(work_path,
                                CONFIG.get('self_shadow', 'incident_target'))
        exiting_target = pjoin(work_path,
                               CONFIG.get('self_shadow', 'exiting_target'))

        # Get the processing tile sizes
        x_tile = int(CONFIG.get('work', 'x_tile_size'))
        y_tile = int(CONFIG.get('work', 'y_tile_size'))
        x_tile = None if x_tile <= 0 else x_tile
        y_tile = None if y_tile <= 0 else y_tile

        # Output target
        self_shadow_target = pjoin(work_path,
                                   CONFIG.get('self_shadow',
                                              'self_shadow_target'))

        gaip.self_shadow(incident_target, exiting_target, self_shadow_target,
                         x_tile, y_tile)


class CalculateCastShadow(luigi.Task):

    """Calculate cast shadow masks. This is a helper task."""

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return [CalculateCastShadowSun(self.l1t_path, self.out_path),
                CalculateCastShadowSatellite(self.l1t_path, self.out_path)]

    def complete(self):
        return all([t.complete() for t in self.requires()])


class CalculateCastShadowSun(luigi.Task):

    """
    Calculates the Cast shadow mask in the direction back to the
    sun.
    """

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return [CalculateSatelliteAndSolarGrids(self.l1t_path, self.out_path),
                DEMExctraction(self.l1t_path, self.out_path)]

    def output(self):
        out_path = pjoin(self.out_path,
                         CONFIG.get('work', 'tc_intermediates'))
        sun_target = pjoin(out_path,
                           CONFIG.get('cast_shadow', 'sun_direction_target'))

        target = luigi.LocalTarget(sun_target)

        return target

    def run(self):
        acqs = gaip.acquisitions(self.l1t_path)
        out_path = self.out_path
        tc_work_path = pjoin(out_path, CONFIG.get('work', 'tc_intermediates'))

        # Input targets
        smoothed_dsm_fname = pjoin(tc_work_path,
                                   CONFIG.get('extract_dsm',
                                              'dsm_smooth_subset'))
        solar_zenith_target = pjoin(out_path,
                                    CONFIG.get('work', 'solar_zenith_target'))
        solar_azimuth_target = pjoin(out_path,
                                     CONFIG.get('work',
                                                'solar_azimuth_target'))
        buffer = int(CONFIG.get('extract_dsm', 'dsm_buffer_width'))
        window_height = int(CONFIG.get('terrain_correction',
                                       'shadow_sub_matrix_height'))
        window_width = int(CONFIG.get('terrain_correction',
                                      'shadow_sub_matrix_width'))

        # Output targets
        sun_target = pjoin(tc_work_path,
                           CONFIG.get('cast_shadow', 'sun_direction_target'))

        gaip.calculate_cast_shadow(acqs[0], smoothed_dsm_fname, buffer,
                                   window_height, window_width,
                                   solar_zenith_target, solar_azimuth_target,
                                   sun_target)


class CalculateCastShadowSatellite(luigi.Task):

    """
    Calculates the Cast shadow mask in the direction back to the
    sun.
    """

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return [CalculateSatelliteAndSolarGrids(self.l1t_path, self.out_path),
                DEMExctraction(self.l1t_path, self.out_path)]

    def output(self):
        out_path = pjoin(self.out_path,
                         CONFIG.get('work', 'tc_intermediates'))
        satellite_target = pjoin(out_path,
                                 CONFIG.get('cast_shadow',
                                            'satellite_direction_target'))
        target = luigi.LocalTarget(satellite_target)

        return target

    def run(self):
        acqs = gaip.acquisitions(self.l1t_path)
        out_path = self.out_path
        tc_work_path = pjoin(out_path, CONFIG.get('work', 'tc_intermediates'))

        # Input targets
        smoothed_dsm_fname = pjoin(tc_work_path,
                                   CONFIG.get('extract_dsm',
                                              'dsm_smooth_subset'))
        satellite_view_target = pjoin(out_path,
                                      CONFIG.get('work', 'sat_view_target'))
        satellite_azimuth_target = pjoin(out_path,
                                         CONFIG.get('work',
                                                    'sat_azimuth_target'))
        buffer = int(CONFIG.get('extract_dsm', 'dsm_buffer_width'))
        window_height = int(CONFIG.get('terrain_correction',
                                       'shadow_sub_matrix_height'))
        window_width = int(CONFIG.get('terrain_correction',
                                      'shadow_sub_matrix_width'))

        # Output targets
        satellite_target = pjoin(tc_work_path,
                                 CONFIG.get('cast_shadow',
                                            'satellite_direction_target'))

        gaip.calculate_cast_shadow(acqs[0], smoothed_dsm_fname, buffer,
                                   window_height, window_width,
                                   satellite_view_target,
                                   satellite_azimuth_target, satellite_target)


class RunTCBand(luigi.Task):

    """Run the terrain correction over a given band."""

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()
    band_num = luigi.IntParameter()

    def requires(self):
        return [BilinearInterpolation(self.l1t_path, self.out_path),
                DEMExctraction(self.l1t_path, self.out_path),
                RelativeAzimuthSlope(self.l1t_path, self.out_path),
                SelfShadow(self.l1t_path, self.out_path),
                CalculateCastShadow(self.l1t_path, self.out_path),
                CreateModisBrdfFiles(self.l1t_path, self.out_path)]

    def output(self):
        acqs = gaip.acquisitions(self.l1t_path)

        # Get the reflectance levels and base output format
        rfl_levels = CONFIG.get('terrain_correction', 'rfl_levels').split(',')
        output_format = CONFIG.get('terrain_correction', 'output_format')

        # Output directory
        out_path = pjoin(self.out_path,
                         CONFIG.get('work', 'rfl_output_dir'))

        # get the acquisition we wish to process
        acqs = [acq for acq in acqs if acq.band_num == self.band_num]

        # Create the targets
        targets = []
        for level in rfl_levels:
            target = pjoin(out_path,
                           output_format.format(level=level,
                                                band=self.band_num))
            targets.append(luigi.LocalTarget(target))
        return targets

    def run(self):
        acqs = gaip.acquisitions(self.l1t_path)
        out_path = self.out_path

        # Get the necessary config params
        tc_path = pjoin(out_path, CONFIG.get('work', 'tc_intermediates'))
        outdir = pjoin(out_path, CONFIG.get('work', 'rfl_output_dir'))
        bilinear_target = pjoin(out_path, 
                                CONFIG.get('work', 'bilinear_outputs_target'))
        bilinear_target = load_value(bilinear_target)
        rori = float(CONFIG.get('terrain_correction', 'rori'))
        modis_brdf_format = pjoin(out_path,
                                  CONFIG.get('brdf', 'modis_brdf_format'))
        new_modis_brdf_format = pjoin(tc_path,
                                      CONFIG.get('brdf',
                                                 'new_modis_brdf_format'))

        # Get the reflectance levels and base output format
        rfl_levels = CONFIG.get('terrain_correction', 'rfl_levels').split(',')
        output_format = CONFIG.get('terrain_correction', 'output_format')

        # Input targets (images)
        self_shadow_target = pjoin(tc_path,
                                   CONFIG.get('self_shadow',
                                              'self_shadow_target'))
        slope_target = pjoin(tc_path,
                             CONFIG.get('self_shadow', 'slope_target'))
        aspect_target = pjoin(tc_path,
                              CONFIG.get('self_shadow', 'aspect_target'))
        incident_target = pjoin(tc_path,
                                CONFIG.get('self_shadow', 'incident_target'))
        exiting_target = pjoin(tc_path,
                               CONFIG.get('self_shadow', 'exiting_target'))
        relative_slope_target = pjoin(tc_path,
                                      CONFIG.get('self_shadow',
                                                 'relative_slope_target'))
        sun_target = pjoin(tc_path,
                           CONFIG.get('cast_shadow', 'sun_direction_target'))
        satellite_target = pjoin(tc_path,
                                 CONFIG.get('cast_shadow',
                                            'satellite_direction_target'))
        solar_zenith_target = pjoin(out_path,
                                    CONFIG.get('work', 'solar_zenith_target'))
        solar_azimuth_target = pjoin(out_path,
                                     CONFIG.get('work',
                                                'solar_azimuth_target'))
        satellite_view_target = pjoin(out_path,
                                      CONFIG.get('work', 'sat_view_target'))
        relative_angle_target = pjoin(out_path,
                                      CONFIG.get('work',
                                                 'relative_azimuth_target'))

        # Get the processing tile sizes
        x_tile = int(CONFIG.get('work', 'x_tile_size'))
        y_tile = int(CONFIG.get('work', 'y_tile_size'))
        x_tile = None if x_tile <= 0 else x_tile
        y_tile = None if y_tile <= 0 else y_tile

        # get the acquisition we wish to process
        acqs = [acq for acq in acqs if acq.band_num == self.band_num]

        # Output targets
        # Create a dict of filenames per reflectance level per band
        rfl_lvl_fnames = {}
        for level in rfl_levels:
            outfname = output_format.format(level=level, band=self.band_num)
            rfl_lvl_fnames[(self.band_num, level)] = pjoin(outdir, outfname)

        # calculate reflectance for lambertian, brdf, and terrain correction 
        gaip.calculate_reflectance(acqs, bilinear_target, rori,
                                   self_shadow_target, sun_target,
                                   satellite_target, solar_zenith_target,
                                   solar_azimuth_target, satellite_view_target,
                                   relative_angle_target, slope_target,
                                   aspect_target, incident_target,
                                   exiting_target, relative_slope_target,
                                   rfl_lvl_fnames, modis_brdf_format,
                                   new_modis_brdf_format, x_tile, y_tile)


class TerrainCorrection(luigi.Task):

    """Perform the terrain correction."""

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        acqs = gaip.acquisitions(self.l1t_path)

        # Retrieve the satellite and sensor for the acquisition
        satellite = acqs[0].spacecraft_id
        sensor = acqs[0].sensor_id

        # Get the required nbar bands list for processing
        nbar_constants = gaip.constants.NBARConstants(satellite, sensor)
        bands_to_process = nbar_constants.get_nbar_lut()

        # define the bands to compute reflectance for
        tc_bands = []
        for band in bands_to_process:
            tc_bands.append(RunTCBand(self.l1t_path, self.out_path, band))

        return tc_bands

    def complete(self):
        return all([t.complete() for t in self.requires()])


class WriteMetadata(luigi.Task):

    """Write metadata."""

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        return [TerrainCorrection(self.l1t_path, self.out_path)]

    def output(self):
        out_path = self.out_path
        out_fname = pjoin(out_path, CONFIG.get('work', 'metadata_target'))
        return luigi.LocalTarget(out_fname)

    def run(self):
        acqs = gaip.acquisitions(self.l1t_path)
        acq = acqs[0]
        out_path = self.out_path

        source_info = {}
        source_info['source_scene'] = self.l1t_path
        source_info['scene_centre_datetime'] = acq.scene_centre_datetime
        source_info['platform'] = acq.spacecraft_id
        source_info['sensor'] = acq.sensor_id
        source_info['path'] = acq.path
        source_info['row'] = acq.row

        # ancillary metadata tracking
        md = gaip.extract_ancillary_metadata(self.l1t_path)
        for key in md:
            source_info[key] = md[key]

        targets = [pjoin(out_path, CONFIG.get('work', 'aerosol_target')),
                   pjoin(out_path, CONFIG.get('work', 'sundist_target')),
                   pjoin(out_path, CONFIG.get('work', 'vapour_target')),
                   pjoin(out_path, CONFIG.get('work', 'ozone_target')),
                   pjoin(out_path, CONFIG.get('work', 'dem_target'))]

        sources = ['aerosol',
                   'solar_distance',
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

        target = pjoin(out_path, CONFIG.get('work', 'brdf_target'))
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
        algorithm['software_version'] = gaip.get_version()
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
        with self.output().open('w') as src:
            yaml.dump(metadata, src, default_flow_style=False)


class Packager(luigi.Task):

    """Packages an nbar or nbart product."""

    l1t_path = luigi.Parameter()
    out_path = luigi.Parameter()
    work_path = luigi.Parameter()
    product = luigi.Parameter()

    def requires(self):
        return [WriteMetadata(self.l1t_path, self.work_path)]

    def output(self):
        out_format = '{}-packaging.completed'
        out_fname = pjoin(self.work_path, out_format.format(self.product))
        return luigi.LocalTarget(out_fname)

    def run(self):
        # run the packager
        kwargs = {'driver': PACKAGE_DRIVERS[self.product],
                  'input_data_paths': [Path(self.work_path)],
                  'destination_path': Path(pjoin(self.out_path, self.product)),
                  'parent_dataset_paths': [Path(self.l1t_path)],
                  'hard_link': False}
        package_newly_processed_data_folder(**kwargs)

        # output a checkpoint
        with self.output().open('w') as src:
            src.write('{} packaging completed'.format(self.product))


class PackageTC(luigi.Task):

    """Issues nbar &/or nbart packaging depending on the config."""

    l1t_path = luigi.Parameter()
    work_path = luigi.Parameter()
    out_path = luigi.Parameter()

    def requires(self):
        products = CONFIG.get('packaging', 'products').split(',')
        tasks = []
        for product in products:
            tasks.append(Packager(self.l1t_path, self.out_path,
                                  self.work_path, product))
        return tasks

    def output(self):
        out_format = '{}.completed'
        base_dir = dirname(self.work_path)
        fname = out_format.format(basename(self.work_path))
        out_fname = pjoin(base_dir, fname)
        return luigi.LocalTarget(out_fname)

    def run(self):
        with self.output().open('w') as src:
            src.write('Task completed')

        # cleanup the entire nbar scene working directory
        cleanup = bool(int(CONFIG.get('cleanup', 'cleanup')))
        if cleanup:
            shutil.rmtree(self.work_path)


def is_valid_directory(parser, arg):
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
 
 
def main(inpath, outpath, workpath, nnodes=1, nodenum=1):
    l1t_files = sorted([pjoin(inpath, f) for f in os.listdir(inpath) if
                        '_OTH_' in f])
    filtered_l1t = []
    for l1t in l1t_files:
        acq = gaip.acquisitions(l1t)[0]
        if ((87 <= acq.path <= 116) & (67 <= acq.row <= 91)):
            filtered_l1t.append(l1t)
        else:
            msg = "Skipping {}".format(acq.dir_name)
            print msg

    # create product output dirs
    products = CONFIG.get('packaging', 'products').split(',')
    for product in products:
        product_dir = pjoin(outpath, product)
        if not exists(product_dir):
            os.makedirs(product_dir)

    l1t_files = [f for f in scatter(filtered_l1t, nnodes, nodenum)]
    nbar_files = [pjoin(workpath, os.path.basename(f).replace('OTH', 'NBAR'))
                  for f in l1t_files]
    # tasks = [TerrainCorrection(l1t, nbar) for l1t, nbar in
    # tasks = [WriteMetadata(l1t, nbar) for l1t, nbar in
    tasks = [PackageTC(l1t, nbar, outpath) for l1t, nbar in
             zip(l1t_files, nbar_files)]
    ncpus = int(os.getenv('PBS_NCPUS', '1'))
    luigi.build(tasks, local_scheduler=True, workers=ncpus / nnodes)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--l1t_path", help=("Path to directory containing L1T "
                        "datasets"), required=True,
                        type=lambda x: is_valid_directory(parser, x))
    parser.add_argument("--out_path", help=("Path to directory where NBAR "
                        "dataset are to be written"), required=True,
                        type=lambda x: is_valid_directory(parser, x))
    parser.add_argument('--cfg',
                        help='Path to a user defined configuration file.')
    parser.add_argument("--log_path", help=("Path to directory where where log"
                        " files will be written"), default='.',
                        type=lambda x: is_valid_directory(parser, x))
    parser.add_argument("--debug", help=("Selects more detail logging (default"
                        " is INFO)"), default=False, action='store_true')
    parser.add_argument("--work_path", help=("Path to a directory where the "
                        "intermediate files will be written."), required=False,
                        type=lambda x: is_valid_directory(parser, x))

    args = parser.parse_args()

    cfg = args.cfg

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

    # use the disk of the local node if we can
    # working directly off the lustre drive seems to flaky
    if args.work_path is None:
        work_path = tempfile.mkdtemp()
    else:
        work_path = args.out_path

    logging.info("nbar.py started")
    logging.info('l1t_path={path}'.format(path=args.l1t_path))
    logging.info('out_path={path}'.format(path=args.out_path))
    logging.info('log_path={path}'.format(path=args.log_path))

    size = int(os.getenv('PBS_NNODES', '1'))
    rank = int(os.getenv('PBS_VNODENUM', '1'))
    main(args.l1t_path, args.out_path, work_path, size, rank)
