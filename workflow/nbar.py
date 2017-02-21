#!/usr/bin/env python
"""
NBAR Workflow
-------------

Workflow settings can be configured in `nbar.cfg` file.

"""
# pylint: disable=missing-docstring,no-init,too-many-function-args
# pylint: disable=too-many-locals

import os
from os.path import join as pjoin, basename, dirname, exists, splitext
import logging
import shutil
import argparse
import luigi
from pathlib import Path
from eodatasets.run import package_newly_processed_data_folder
from eodatasets.drivers import PACKAGE_DRIVERS
from eodatasets import type as ptype
import gaip


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


def get_buffer(group):
    buf = {'product': 250,
           'R10m': 700,
           'R20m': 350,
           'R60m': 120}
    return buf[group]


def get_tile_sizes():
    x_tile = CONFIG.getint('work', 'x_tile_size')
    y_tile = CONFIG.getint('work', 'y_tile_size')
    x_tile = None if x_tile <= 0 else x_tile
    y_tile = None if y_tile <= 0 else y_tile

    return x_tile, y_tile


class GetAncillaryData(luigi.Task):

    """Get all ancillary data."""

    level1 = luigi.Parameter()
    nbar_root = luigi.Parameter()
    granule = luigi.Parameter()

    def output(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, granule=self.granule)
        out_fname = pjoin(out_path, CONFIG.get('work', 'ancillary_fname'))
        return luigi.LocalTarget(out_fname)

    def run(self):
        container = gaip.acquisitions(self.level1)
        acqs = container.get_acquisitions(granule=self.granule)

        aerosol_path = CONFIG.get('ancillary', 'aerosol_path')
        water_vapour_path = CONFIG.get('ancillary', 'vapour_path')
        ozone_path = CONFIG.get('ancillary', 'ozone_path')
        dem_path = CONFIG.get('ancillary', 'dem_path')
        brdf_path = CONFIG.get('ancillary', 'brdf_path')
        brdf_premodis_path = CONFIG.get('ancillary', 'brdf_premodis_path')

        with self.output().temporary_path() as out_fname:
            gaip._collect_ancillary_data(acqs[0], aerosol_path,
                                         water_vapour_path, ozone_path,
                                         dem_path, brdf_path,
                                         brdf_premodis_path, out_fname)


class CalculateLonGrid(luigi.Task):

    """Calculate the longitude grid."""

    level1 = luigi.Parameter()
    nbar_root = luigi.Parameter()
    granule = luigi.Parameter()
    group = luigi.Parameter()

    def output(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, group=self.group,
                                      granule=self.granule)
        out_fname = pjoin(out_path, CONFIG.get('work', 'lon_grid_fname'))
        return luigi.LocalTarget(out_fname)

    def run(self):
        container = gaip.acquisitions(self.level1)
        acqs = container.get_acquisitions(group=self.group,
                                          granule=self.granule)
        compression = CONFIG.get('work', 'compression')

        with self.output().temporary_path() as out_fname:
            gaip.create_lon_grid(acqs[0], out_fname, compression=compression)


class CalculateLatGrid(luigi.Task):

    """Calculate the latitude grid."""

    level1 = luigi.Parameter()
    nbar_root = luigi.Parameter()
    granule = luigi.Parameter()
    group = luigi.Parameter()

    def output(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, group=self.group,
                                      granule=self.granule)
        out_fname = pjoin(out_path, CONFIG.get('work', 'lat_grid_fname'))
        return luigi.LocalTarget(out_fname)

    def run(self):
        container = gaip.acquisitions(self.level1)
        acqs = container.get_acquisitions(group=self.group,
                                          granule=self.granule)
        compression = CONFIG.get('work', 'compression')

        with self.output().temporary_path() as out_fname:
            gaip.create_lat_grid(acqs[0], out_fname, compression=compression)


class CalculateSatelliteAndSolarGrids(luigi.Task):

    """Calculate the satellite and solar grids."""

    level1 = luigi.Parameter()
    nbar_root = luigi.Parameter()
    granule = luigi.Parameter()
    group = luigi.Parameter()

    def requires(self):
        args = [self.level1, self.nbar_root, self.granule, self.group]
        return {'lat': CalculateLatGrid(*args),
                'lon': CalculateLonGrid(*args)}

    def output(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, group=self.group,
                                      granule=self.granule)
        out_fname = CONFIG.get('work', 'satellite_solar_fname')

        return luigi.LocalTarget(pjoin(out_path, out_fname))

    def run(self):
        container = gaip.acquisitions(self.level1)
        lat_fname = self.input()['lat'].path
        lon_fname = self.input()['lon'].path

        acqs = container.get_acquisitions(group=self.group,
                                          granule=self.granule)

        view_max = acqs[0].maximum_view_angle
        tle_path = CONFIG.get('work', 'tle_path')
        compression = CONFIG.get('work', 'compression')

        with self.output().temporary_path() as out_fname:
            gaip._calculate_angles(acqs[0], lon_fname, lat_fname, out_fname,
                                   npoints=12, compression=compression,
                                   max_angle=view_max, tle_path=tle_path)


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
        tasks = {}

        for granule in container.granules:
            key1 = (granule, 'ancillary')
            args1 = [self.level1, self.nbar_root, granule]
            tasks[key1] = GetAncillaryData(*args1)
            for group in container.groups:
                key2 = (granule, group)
                args2 = [self.level1, self.nbar_root, granule, group]
                tsks = {'sat_sol': CalculateSatelliteAndSolarGrids(*args2),
                        'lat': CalculateLatGrid(*args2),
                        'lon': CalculateLonGrid(*args2)}
                tasks[key2] = tsks

        return tasks

    def output(self):
        container = gaip.acquisitions(self.level1)
        grn_path = container.get_root(self.nbar_root, granule=self.granule)
        modtran_root = pjoin(grn_path, CONFIG.get('work', 'modtran_root'))

        coords = CONFIG.get('modtran', 'coords').split(',')
        albedos = CONFIG.get('modtran', 'albedos').split(',')
        out_format = CONFIG.get('write_tp5', 'output_format')

        targets = {}
        for coord in coords:
            for albedo in albedos:
                fname = out_format.format(coord=coord, albedo=albedo)
                out_fname = pjoin(modtran_root, fname)
                targets[(coord, albedo)] = luigi.LocalTarget(out_fname)
        return targets

    def run(self):
        container = gaip.acquisitions(self.level1)

        # as we have an all granules groups dependency, it doesn't matter which
        # group, so just get the first
        group = container.groups[0]

        coords = CONFIG.get('modtran', 'coords').split(',')
        albedos = CONFIG.get('modtran', 'albedos').split(',')

        # input file
        inputs = self.input()

        # all ancillary filenames from each granule
        fnames = [inputs[key].path for key in inputs if 'ancillary' in key]

        if container.tiled:
            ancillary_fname = pjoin(self.nbar_root, 'averaged-ancillary.h5')
            gaip.aggregate_ancillary(fnames, ancillary_fname)
        else:
            ancillary_fname = fnames[0]

        sat_sol_fname = inputs[(self.granule, group)]['sat_sol'].path
        lon_fname = inputs[(self.granule, group)]['lon'].path
        lat_fname = inputs[(self.granule, group)]['lat'].path

        # load an acquisition
        acq = container.get_acquisitions(group=group, granule=self.granule)[0]

        tp5_data = gaip._format_tp5(acq, sat_sol_fname, lon_fname, lat_fname,
                                    ancillary_fname, coords, albedos)

        for key, target in self.output().items():
            with target.temporary_path() as out_fname:
                with open(out_fname, 'w') as src:
                    src.writelines(tp5_data[key])


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
        tasks = {'tp5': WriteTp5(*args)}

        return tasks

    def output(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, granule=self.granule)
        modtran_root = pjoin(out_path, CONFIG.get('work', 'modtran_root'))
        output_format = CONFIG.get('extract_flux', 'input_format')
        out_fname = output_format.format(coord=self.coord, albedo=self.albedo)

        return luigi.LocalTarget(pjoin(modtran_root, out_fname))

    def run(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, granule=self.granule)

        modtran_exe = CONFIG.get('modtran', 'exe')
        workpath_format = CONFIG.get('modtran', 'workpath_format')
        modtran_root = pjoin(out_path, CONFIG.get('work', 'modtran_root'))

        workpath = workpath_format.format(coord=self.coord, albedo=self.albedo)
        modtran_work = pjoin(modtran_root, workpath)

        gaip.prepare_modtran(self.coord, self.albedo, modtran_work,
                             modtran_exe)
        gaip.run_modtran(modtran_exe, modtran_work)


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
        args = [self.level1, self.nbar_root, self.granule,
                self.coord, self.albedo]
        return RunModtranCase(*args)

    def output(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, granule=self.granule)
        modtran_root = pjoin(out_path, CONFIG.get('work', 'modtran_root'))
        out_fmt = CONFIG.get('extract_flux', 'output_format')
        out_fname = pjoin(modtran_root, out_fmt.format(coord=self.coord,
                                                       albedo=self.albedo))
        return luigi.LocalTarget(out_fname)

    def run(self):
        container = gaip.acquisitions(self.level1)
        acqs = container.get_acquisitions(granule=self.granule)

        flux_fname = self.input().path

        satfilterpath = CONFIG.get('ancillary', 'satfilter_path')
        response_fname = pjoin(satfilterpath, acqs[0].spectral_filter_file)

        transmittance = True if self.albedo == 't' else False
        compression = CONFIG.get('work', 'compression')

        with self.output().temporary_path() as out_fname:
            gaip._calculate_solar_radiation(flux_fname, response_fname,
                                            transmittance, out_fname,
                                            compression)


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

        for coord in coords:
            for albedo in albedos:
                args = [self.level1, self.nbar_root, self.granule, coord,
                        albedo]
                yield RunAccumulateSolarIrradianceCase(*args)

    def output(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, granule=self.granule)
        modtran_root = pjoin(out_path, CONFIG.get('work', 'modtran_root'))
        out_fname = CONFIG.get('extract_flux', 'out_fname')
        out_fname = pjoin(modtran_root, out_fname)
        return luigi.LocalTarget(out_fname)

    def run(self):
        compression = CONFIG.get('work', 'compression')
        with self.output().temporary_path() as out_fname:
            for target in self.input():
                gaip.create_solar_irradiance_tables(target.path, out_fname,
                                                    compression=compression)
                target.remove()


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
        return AccumulateSolarIrradiance(*args)

    def output(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, granule=self.granule)
        out_path = pjoin(out_path, CONFIG.get('work', 'modtran_root'))
        out_fname = pjoin(out_path, CONFIG.get('coefficients', 'out_fname'))
        return luigi.LocalTarget(out_fname)
               

    def run(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, granule=self.granule)
        workpath = pjoin(out_path, CONFIG.get('work', 'modtran_root'))

        coords = CONFIG.get('modtran', 'coords').split(',')
        chn_input_fmt = CONFIG.get('coefficients', 'chn_input_format')
        dir_input_fmt = CONFIG.get('extract_flux', 'output_format')
        compression = CONFIG.get('work', 'compression')

        with self.output().temporary_path() as out_fname:
            gaip._calculate_coefficients(coords, chn_input_fmt, dir_input_fmt,
                                         workpath, out_fname, compression)


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
        return {'coef': CalculateCoefficients(*args1),
                'satsol': CalculateSatelliteAndSolarGrids(*args2)}

    def output(self):
        band_num = self.band_num
        factor = self.factor
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, group=self.group,
                                      granule=self.granule)
        workpath = pjoin(out_path,
                         CONFIG.get('work', 'bilinear_root'))
        output_format = CONFIG.get('bilinear', 'output_format')
        out_fname = pjoin(workpath, output_format.format(band=band_num,
                                                         factor=factor))
        return luigi.LocalTarget(out_fname)

    def run(self):
        container = gaip.acquisitions(self.level1)
        sat_sol_angles_fname = self.input()['satsol'].path
        coefficients_fname = self.input()['coef'].path
        compression = CONFIG.get('work', 'compression')

        acqs = container.get_acquisitions(group=self.group,
                                          granule=self.granule)
        acq = [acq for acq in acqs if acq.band_num == self.band_num][0]

        with self.output().temporary_path() as out_fname:
            gaip._bilinear_interpolate(acq, self.factor, sat_sol_angles_fname,
                                       coefficients_fname, out_fname,
                                       compression)


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
        tasks = {}
        for factor in factors:
            for band in bands:
                key = (band, factor)
                args = [self.level1, self.nbar_root, self.granule, self.group,
                        band, factor]
                tasks[key] = BilinearInterpolationBand(*args)
        return tasks

    def output(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, group=self.group,
                                      granule=self.granule)
        out_fname = pjoin(out_path, CONFIG.get('bilinear', 'out_fname'))
        return luigi.LocalTarget(out_fname)

    def run(self):
        bilinear_fnames = {}
        for key, value in self.input().items():
            bilinear_fnames[key] = value.path

        with self.output().temporary_path() as out_fname:
            gaip.link_bilinear_data(bilinear_fnames, out_fname)


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
        return

    def output(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, group=self.group,
                                      granule=self.granule)
        out_path = pjoin(out_path, CONFIG.get('work', 'tc_root'))
        out_fname = CONFIG.get('extract_dsm', 'dsm_fname')
        return luigi.LocalTarget(pjoin(out_path, out_fname))

    def run(self):
        container = gaip.acquisitions(self.level1)
        acqs = container.get_acquisitions(group=self.group,
                                          granule=self.granule)

        national_dsm = CONFIG.get('ancillary', 'dem_tc')
        margins = get_buffer(self.group)
        compression = CONFIG.get('work', 'compression')

        with self.output().temporary_path() as out_fname:
            _ = gaip.get_dsm(acqs[0], national_dsm, margins, out_fname,
                             compression)


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
        return DEMExctraction(*args)

    def output(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, group=self.group,
                                      granule=self.granule)
        out_path = pjoin(out_path, CONFIG.get('work', 'tc_root'))
        out_fname = pjoin(out_path,
                          CONFIG.get('self_shadow', 'slope_aspect_fname'))

        return luigi.LocalTarget(out_fname)

    def run(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, group=self.group,
                                      granule=self.granule)

        acqs = container.get_acquisitions(group=self.group,
                                          granule=self.granule)

        dsm_fname = self.input().path
        margins = get_buffer(self.group)
        compression = CONFIG.get('work', 'compression')

        with self.output().temporary_path() as out_fname:
            gaip._slope_aspect_arrays_wrapper(acqs[0], dsm_fname, margins,
                                              out_fname, compression)


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
        return {'sat_sol': CalculateSatelliteAndSolarGrids(*args),
                'slp_asp': SlopeAndAspect(*args)}

    def output(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, group=self.group,
                                      granule=self.granule)
        out_path = pjoin(out_path, CONFIG.get('work', 'tc_root'))
        out_fname = pjoin(out_path,
                          CONFIG.get('self_shadow', 'incident_fname'))
        return luigi.LocalTarget(out_fname)

    def run(self):
        # input filenames
        sat_sol_fname = self.input()['sat_sol'].path
        slope_aspect_fname = self.input()['slp_asp'].path

        # get the processing tile sizes
        x_tile, y_tile = get_tile_sizes()
        compression = CONFIG.get('work', 'compression')

        with self.output().temporary_path() as out_fname:
            gaip._incident_angles(sat_sol_fname, slope_aspect_fname, out_fname,
                                  compression, x_tile, y_tile)


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
        return {'sat_sol': CalculateSatelliteAndSolarGrids(*args),
                'slp_asp': SlopeAndAspect(*args)}

    def output(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, group=self.group,
                                      granule=self.granule)
        out_path = pjoin(out_path, CONFIG.get('work', 'tc_root'))
        out_fname = pjoin(out_path, CONFIG.get('self_shadow', 'exiting_fname'))
        return luigi.LocalTarget(out_fname)

    def run(self):
        # input filenames
        sat_sol_fname = self.input()['sat_sol'].path
        slope_aspect_fname = self.input()['slp_asp'].path

        # get the processing tile sizes
        x_tile, y_tile = get_tile_sizes()
        compression = CONFIG.get('work', 'compression')

        with self.output().temporary_path() as out_fname:
            gaip._exiting_angles(sat_sol_fname, slope_aspect_fname, out_fname,
                                 compression, x_tile, y_tile)


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
        return {'incident': IncidentAngles(*args),
                'exiting': ExitingAngles(*args)}

    def output(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, group=self.group,
                                      granule=self.granule)
        out_path = pjoin(out_path, CONFIG.get('work', 'tc_root'))
        out_fname = pjoin(out_path,
                          CONFIG.get('self_shadow', 'relative_slope_fname'))
        return luigi.LocalTarget(out_fname)

    def run(self):
        # input filenames
        incident_fname = self.input()['incident'].path
        exiting_fname = self.input()['exiting'].path

        # get the processing tile sizes
        x_tile, y_tile = get_tile_sizes()
        compression = CONFIG.get('work', 'compression')

        with self.output().temporary_path() as out_fname:
            gaip._relative_azimuth_slope(incident_fname, exiting_fname,
                                         out_fname, compression, x_tile,
                                         y_tile)


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
        return {'incident': IncidentAngles(*args),
                'exiting': ExitingAngles(*args)}

    def output(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, group=self.group,
                                      granule=self.granule)
        out_path = pjoin(out_path, CONFIG.get('work', 'tc_root'))
        out_fname = pjoin(out_path, CONFIG.get('self_shadow',
                                               'self_shadow_fname'))
        return luigi.LocalTarget(out_fname)

    def run(self):
        # input filenames
        incident_fname = self.input()['incident'].path
        exiting_fname = self.input()['exiting'].path

        # get the processing tile sizes
        x_tile, y_tile = get_tile_sizes()
        compression = CONFIG.get('work', 'compression')

        with self.output().temporary_path() as out_fname:
            gaip._self_shadow(incident_fname, exiting_fname, out_fname,
                              compression, x_tile, y_tile)


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
        return {'sat_sol': CalculateSatelliteAndSolarGrids(*args),
                'dsm': DEMExctraction(*args)}

    def output(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, group=self.group,
                                      granule=self.granule)
        out_path = pjoin(out_path, CONFIG.get('work', 'tc_root'))
        out_fname = pjoin(out_path, CONFIG.get('cast_shadow',
                                               'sun_direction_fname'))
        return luigi.LocalTarget(out_fname)

    def run(self):
        container = gaip.acquisitions(self.level1)
        acqs = container.get_acquisitions(group=self.group,
                                          granule=self.granule)

        # input filenames
        dsm_fname = self.input()['dsm'].path
        sat_sol_fname = self.input()['sat_sol'].path

        margins = get_buffer(self.group)
        compression = CONFIG.get('work', 'compression')
        window_height = CONFIG.getint('terrain_correction',
                                      'shadow_sub_matrix_height')
        window_width = CONFIG.getint('terrain_correction',
                                     'shadow_sub_matrix_width')

        with self.output().temporary_path() as out_fname:
            gaip._calculate_cast_shadow(acqs[0], dsm_fname, margins,
                                        window_height, window_width,
                                        sat_sol_fname, out_fname, compression)


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
        return {'sat_sol': CalculateSatelliteAndSolarGrids(*args),
                'dsm': DEMExctraction(*args)}

    def output(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, group=self.group,
                                      granule=self.granule)
        out_path = pjoin(out_path, CONFIG.get('work', 'tc_root'))
        out_fname = pjoin(out_path, CONFIG.get('cast_shadow',
                                               'satellite_direction_fname'))
        return luigi.LocalTarget(out_fname)

    def run(self):
        container = gaip.acquisitions(self.level1)
        acqs = container.get_acquisitions(group=self.group,
                                          granule=self.granule)

        # input filenames
        dsm_fname = self.input()['dsm'].path
        sat_sol_fname = self.input()['sat_sol'].path

        margins = get_buffer(self.group)
        compression = CONFIG.get('work', 'compression')
        window_height = CONFIG.getint('terrain_correction',
                                      'shadow_sub_matrix_height')
        window_width = CONFIG.getint('terrain_correction',
                                     'shadow_sub_matrix_width')

        with self.output().temporary_path() as out_fname:
            gaip._calculate_cast_shadow(acqs[0], dsm_fname, margins,
                                        window_height, window_width,
                                        sat_sol_fname, out_fname, compression,
                                        False)


class CalculateShadowMasks(luigi.Task):

    """
    Issues self and cast shadow tasks for two direction sources;
    the sun and the satellite. Acts as a helper task,
    but combines the results into a single file.
    """

    level1 = luigi.Parameter()
    nbar_root = luigi.Parameter()
    granule = luigi.Parameter()
    group = luigi.Parameter()

    def requires(self):
        args = [self.level1, self.nbar_root, self.granule, self.group]
        return {'sun': CalculateCastShadowSun(*args),
                'sat': CalculateCastShadowSatellite(*args),
                'self': SelfShadow(*args)}

    def output(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, group=self.group,
                                      granule=self.granule)
        out_path = pjoin(out_path, CONFIG.get('work', 'tc_root'))
        out_fname = pjoin(out_path, CONFIG.get('cast_shadow', 'out_fname'))
        return luigi.LocalTarget(out_fname)

    def run(self):
        with self.output().temporary_path() as out_fname:
            for key in self.input():
                fname = self.input()[key].path
                dname = splitext(basename(fname))[0]
                gaip.create_external_link(fname, dname, out_fname, dname)


class RunTCBand(luigi.Task):

    """Run the terrain correction over a given band."""

    level1 = luigi.Parameter()
    nbar_root = luigi.Parameter()
    granule = luigi.Parameter()
    group = luigi.Parameter()
    band_num = luigi.IntParameter()

    def requires(self):
        args = [self.level1, self.nbar_root, self.granule, self.group]
        reqs = {'bilinear': BilinearInterpolation(*args),
                'ancillary': GetAncillaryData(*args[:-1]),
                'dsm': DEMExctraction(*args),
                'rel_slope': RelativeAzimuthSlope(*args),
                'shadow': CalculateShadowMasks(*args),
                'slp_asp': SlopeAndAspect(*args),
                'incident': IncidentAngles(*args),
                'exiting': ExitingAngles(*args),
                'sat_sol': CalculateSatelliteAndSolarGrids(*args)}

        return reqs

    def output(self):
        band = self.band_num
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, group=self.group,
                                      granule=self.granule)
        reflectance_root = CONFIG.get('work', 'reflectance_root')
        out_path = pjoin(out_path, reflectance_root)
        out_fname = CONFIG.get('terrain_correction', 'out_fname')

        return luigi.LocalTarget(pjoin(out_path, out_fname.format(band=band)))

    def run(self):
        container = gaip.acquisitions(self.level1)
        acqs = container.get_acquisitions(group=self.group,
                                          granule=self.granule)

        # TODO: what is rori???
        rori = CONFIG.getfloat('terrain_correction', 'rori')
        compression = CONFIG.get('work', 'compression')

        # inputs
        inputs = self.input()
        bilinear_fname = inputs['bilinear'].path
        slp_asp_fname = inputs['slp_asp'].path
        incident_fname = inputs['incident'].path
        exiting_fname = inputs['exiting'].path
        relative_slope_fname = inputs['rel_slope'].path
        shadow_fname = inputs['shadow'].path
        sat_sol_fname = inputs['sat_sol'].path

        # get the processing tile sizes
        x_tile, y_tile = get_tile_sizes()

        # get the acquisition we wish to process
        acq = [acq for acq in acqs if acq.band_num == self.band_num][0]

        if container.tiled:
            ancillary_fname = pjoin(self.nbar_root, 'averaged-ancillary.h5')
        else:
            ancillary_fname = inputs['ancillary'].path

        with self.output().temporary_path() as out_fname:
            gaip._calculate_reflectance(acq, bilinear_fname, sat_sol_fname,
                                        slp_asp_fname, relative_slope_fname,
                                        incident_fname, exiting_fname,
                                        shadow_fname, ancillary_fname, rori,
                                        out_fname, compression, x_tile, y_tile)


class TerrainCorrection(luigi.WrapperTask):

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
        for band in bands_to_process:
            yield RunTCBand(self.level1, self.nbar_root, self.granule,
                            self.group, band)


# TODO: re-work with modified I/O
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
        out_fname = pjoin(out_path, CONFIG.get('work', 'metadata_fname'))
        return luigi.LocalTarget(out_fname)

    def run(self):
        # TODO: retrieve single acquisition
        # TODO: a single acquisition doesn't account for all granules
        # TODO: do we write metadata per granule???
        # TODO: rework to use self.input()
        container = gaip.acquisitions(self.level1)
        acq = container.get_acquisitions()[0]
        out_path = self.nbar_root

        targets = [pjoin(out_path, CONFIG.get('work', 'aerosol_fname')),
                   pjoin(out_path, CONFIG.get('work', 'vapour_fname')),
                   pjoin(out_path, CONFIG.get('work', 'ozone_fname')),
                   pjoin(out_path, CONFIG.get('work', 'dem_fname')),
                   pjoin(out_path, CONFIG.get('work', 'brdf_fname'))]

        aerosol_data = load(luigi.LocalTarget(targets[0]))
        water_vapour_data = load(luigi.LocalTarget(targets[1]))
        ozone_data = load(luigi.LocalTarget(targets[2]))
        elevation_data = load(luigi.LocalTarget(targets[3]))
        brdf_data = load(luigi.LocalTarget(targets[4]))

        # output
        with self.output().temporary_path() as out_fname:
            gaip.write_nbar_yaml(acq, self.level1, ozone_data, aerosol_data,
                                 water_vapour_data, elevation_data, brdf_data,
                                 out_fname)


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


class RunGaip(luigi.WrapperTask):

    level1_list = luigi.Parameter()
    work_root = luigi.Parameter()
    nbar_root = luigi.Parameter()

    def requires(self):
        with open(self.level1_list) as src:
           scenes = src.readlines()
        for level1 in scenes:
            nbar_root = pjoin(work_root, basename(level1) + ".nbar-work")
            yield PackageTC(level1, self.work_root, nbar_root)


def scatter(iterable, P=1, p=1):
    """
    Scatter an iterator across `P` processors where `p` is the index
    of the current processor. This partitions the work evenly across
    processors.
    """
    import itertools
    return itertools.islice(iterable, p-1, None, P)
 
 
if __name__ == '__main__':
    luigi.run()
