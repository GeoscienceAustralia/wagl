#!/usr/bin/env python
"""
NBAR Workflow
-------------

Workflow settings can be configured in `nbar.cfg` file.

"""
# pylint: disable=missing-docstring,no-init,too-many-function-args
# pylint: disable=too-many-locals

import cPickle as pickle
import os
from os.path import join as pjoin, basename, dirname, exists
import subprocess
import logging
import shutil
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


class GetAncillarData(luigi.Task):

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
        out_path = container.get_root(self.nbar_root, granule=self.granule)

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
        return {'lon': luigi.LocalTarget(out_fname)}

    def run(self):
        container = gaip.acquisitions(self.level1)
        acqs = container.get_acquisitions(group=self.group,
                                          granule=self.granule)
        compression = CONFIG.get('work', 'compression')

        with self.output()['lon'].temporary_path() as out_fname:
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
        return {'lat': luigi.LocalTarget(out_fname)}

    def run(self):
        container = gaip.acquisitions(self.level1)
        acqs = container.get_acquisitions(group=self.group,
                                          granule=self.granule)
        compression = CONFIG.get('work', 'compression')

        with self.output()['lat'].temporary_path() as out_fname:
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
        lat_fname = self.input()['lat']['lat'].path
        lon_fname = self.input()['lon']['lon'].path

        acqs = container.get_acquisitions(group=self.group,
                                          granule=self.granule)

        view_max = acqs[0].maximum_view_angle
        tle_path = CONFIG.get('work', 'tle_path')
        compression = CONFIG.get('work', 'compression')

        with self.output().temporary_path() as out_fname:
            _calculate_angles(acqs[0], lon_fname, lat_fname, out_fname,
                              npoints=12, compression=compression,
                              max_angle=view_max, tle_path=tle_path)


# TODO: re-work with hdf5; re-writing into the ancillary file removes atomicity
class AggregateAncillary(luigi.Task):

    """Aggregates the ancillary data from each granule."""

    level1 = luigi.Parameter()
    nbar_root = luigi.Parameter()

    def requires(self):
        container = gaip.acquisitions(self.level1)
        tasks = {}

        for granule in container.granules:
            args = [self.level1, self.nbar_root, granule]
            tasks[granule] = GetAncillaryData(*args)

        return tasks

    def output(self):
        out_path = self.nbar_root
        ozone_fname = pjoin(out_path, CONFIG.get('work', 'ozone_fname'))
        vapour_fname = pjoin(out_path, CONFIG.get('work', 'vapour_fname'))
        aerosol_fname = pjoin(out_path, CONFIG.get('work', 'aerosol_fname'))
        elevation_fname = pjoin(out_path, CONFIG.get('work', 'dem_fname'))
        targets = {'ozone': luigi.LocalTarget(ozone_fname),
                   'vapour': luigi.LocalTarget(vapour_fname),
                   'aerosol': luigi.LocalTarget(aerosol_fname),
                   'elevation': luigi.LocalTarget(elevation_fname)}
        return targets

    def run(self):
        container = gaip.acquisitions(self.level1)
        n_tiles = len(container.granules)

        # initialise the mean result
        ozone = vapour = aerosol = elevation = 0.0

        for _, value in self.requires().inputs().items():
            ozone_fname = value['ozone'].path
            vapour_fname = value['vapour'].path
            aerosol_fname = value['aerosol'].path
            elevation_fname = value['elevation'].path

            ozone += load_value(ozone_fname)
            vapour += load_value(vapour_fname)
            aerosol += load_value(aerosol_fname)
            elevation += load_value(elevation_fname)

        ozone /= n_tiles
        vapour /= n_tiles
        aerosol /= n_tiles
        elevation /= n_tiles

        # write the mean ancillary values
        outputs = self.output()
        data = {'data_source': 'granule_average'}

        with outputs['ozone'].temporary_path() as out_fname1,\
            outputs['vapour'].temporary_path() as out_fname2,\
            outputs['aerosol'].temporary_path() as out_fname3,\
            outputs['elevation'].temporary_path() as out_fname4:

            data['value'] = ozone
            save(out_fname1, data)

            data['value'] = vapour
            save(out_fname2, data)

            data['value'] = aerosol
            save(out_fname3, data)

            data['value'] = elevation
            save(out_fname4, data)


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

        # are we dealing with tiles/granules? aggregate the ancillary values
        if container.tiled:
            args = [self.level1, self.nbar_root]
            tasks['ancillary'] = AggregateAncillary(*args)

        for granule in container.granules:
            # grab non-aggregated ancillary
            if not container.tiled:
                key = (granule, 'ancillary')
                args = [self.level1, self.nbar_root, granule]
                tasks[key] = GetAncillaryData(*args)
            for group in container.groups:
                key = (granule, group)
                args = [self.level1, self.nbar_root, granule, group]
                tsks = {'sat_sol': CalculateSatelliteAndSolarGrids(*args),
                        'lat': CalculateLatGrid(*args),
                        'lon': CalculateLonGrid(*args)}
                tasks[key] = tsks

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
        # inputs = self.requires().input()
        inputs = self.input()

        if container.tiled:
            ancillary = inputs['ancillary']
        else:
            ancillary = self.requires()[(self.granule, 'ancillary')].input()
            # ancillary = inputs[(self.granule, 'ancillary')]

        sat_sol_fname = inputs[(self.granule, group)]['sat_sol'].path
        lon_fname = inputs[(self.granule, group)]['lon']['lon'].path
        lat_fname = inputs[(self.granule, group)]['lat']['lat'].path

        ancillary_fname = ancillary.path

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

        # do we need to output the legacy files?
        if CONFIG.getboolean('legacy_outputs', 'required'):
            tasks = {'legacy': LegacyOutputs(*args),
                     'tp5': WriteTp5(*args)}
        else:
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
            _calculate_solar_radiation(flux_fname, response_fname,
                                       transmittance, out_fname, compression)


# TODO: do we need wrapper tasks???
class AccumulateSolarIrradiance(luigi.WrapperTask):

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

    def output(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, granule=self.granule)
        modtran_root = pjoin(out_path, CONFIG.get('work', 'modtran_root'))
        out_fname = CONFIG.get('extract_flux', 'out_fname')
        out_fname = pjoin(modtran_root, out_fname)
        return luigi.LocalTarget(out_fname)

    def run(self):
        with self.output().temporary_path() as out_fname:
            for target in self.input():
                create_solar_irradiance_tables(target.path, out_fname,
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
        return [AccumulateSolarIrradiance(*args)]

    def output(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, granule=self.granule)
        out_path = pjoin(out_path, CONFIG.get('work', 'modtran_root'))
        out_fname1 = pjoin(out_path, CONFIG.get('coefficients', 'out_fname1'))
        out_fname2 = pjoin(out_path, CONFIG.get('coefficients', 'out_fname2'))
        targets = {'coef1': luigi.LocalTarget(out_fname1),
                   'coef2': luigi.LocalTarget(out_fname2)}
        return targets
               

    def run(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, granule=self.granule)
        workpath = pjoin(out_path, CONFIG.get('work', 'modtran_root'))

        coords = CONFIG.get('modtran', 'coords').split(',')
        chn_input_format = CONFIG.get('coefficients', 'chn_input_format')
        dir_input_format = CONFIG.get('extract_flux', 'output_format')

        coef1, coef2 = gaip.calculate_coefficients(coords, chn_input_format,
                                                   dir_input_format,
                                                   workpath)

        with self.output()['coef1'].temporary_path() as out_fname1,\
            self.output()['coef2'].temporary_path() as out_fname2:

            print "---out_fname1---", out_fname1
            print "---out_fname2---", out_fname2
            save(out_fname1, coef1)
            save(out_fname2, coef2)


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
        coord_fname = self.input()['satsol']['coordinator'].path
        box_fname = self.input()['satsol']['boxline'].path
        centre_fname = self.input()['satsol']['centreline'].path

        band_num = self.band_num
        factor = self.factor
        coefficients = load(self.input()['coef']['coef2'])[(band_num, factor)]

        acqs = container.get_acquisitions(group=self.group,
                                          granule=self.granule)

        acq = [acq for acq in acqs if acq.band_num == band_num][0]
        with self.output().temporary_path() as out_fname:
            gaip.bilinear_interpolate(acq, coord_fname, box_fname,
                                      centre_fname, coefficients, out_fname)


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
        out_fname = pjoin(out_path,
                          CONFIG.get('work', 'bilinear_outputs_fname'))
        return luigi.LocalTarget(out_fname)

    def run(self):
        # Initialise the dict to store the locations of the bilinear outputs
        bilinear_fnames = {}

        for key, value in self.input().items():
            bilinear_fnames[key] = value.path

        with self.output().temporary_path() as out_fname:
            save(out_fname, bilinear_fnames)


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
        subset_fname = CONFIG.get('extract_dsm', 'dsm_subset')
        smoothed_fname = CONFIG.get('extract_dsm', 'dsm_smooth_subset')
        return {'subset': luigi.LocalTarget(pjoin(out_path, subset_fname)),
                'smoothed': luigi.LocalTarget(pjoin(out_path, smoothed_fname))}

    def run(self):
        container = gaip.acquisitions(self.level1)
        acqs = container.get_acquisitions(group=self.group,
                                          granule=self.granule)

        national_dsm = CONFIG.get('ancillary', 'dem_tc')
        #buffer = int(CONFIG.get('extract_dsm', 'dsm_buffer_width'))
        buffer = get_buffer(self.group)

        with self.output()['subset'].temporary_path() as subset_fname,\
            self.output()['smoothed'].temporary_path() as smoothed_fname:

            gaip.get_dsm(acqs[0], national_dsm, buffer, subset_fname,
                         smoothed_fname)


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
        slope_fname = pjoin(out_path,
                            CONFIG.get('self_shadow', 'slope_fname'))
        aspect_fname = pjoin(out_path,
                             CONFIG.get('self_shadow', 'aspect_fname'))
        targets = {'slope': luigi.LocalTarget(slope_fname),
                   'aspect': luigi.LocalTarget(aspect_fname)}

        return targets

    def run(self):
        container = gaip.acquisitions(self.level1)
        out_path = container.get_root(self.nbar_root, group=self.group,
                                      granule=self.granule)
        work_path = pjoin(out_path, CONFIG.get('work', 'tc_root'))

        acqs = container.get_acquisitions(group=self.group,
                                          granule=self.granule)

        dsm_fname = self.input()['smoothed'].path

        #margins = int(CONFIG.get('extract_dsm', 'dsm_buffer_width'))
        margins = get_buffer(self.group)

        # TODO: remove header_slope file
        # Output filenames
        header_slope_fname = pjoin(work_path,
                                   CONFIG.get('work', 'header_slope_fname'))

        with self.output()['slope'].temporary_path() as slope_fname,\
            self.output()['aspect'].temporary_path() as aspect_fname:

            gaip.slope_aspect_arrays(acqs[0], dsm_fname, margins, slope_fname,
                                     aspect_fname, header_slope_fname)


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
        incident_fname = pjoin(out_path,
                               CONFIG.get('self_shadow', 'incident_fname'))
        azi_incident_fname = pjoin(out_path,
                                   CONFIG.get('self_shadow',
                                              'azimuth_incident_fname'))
        targets = {'incident': luigi.LocalTarget(incident_fname),
                   'azi_incident': luigi.LocalTarget(azi_incident_fname)}
        return targets

    def run(self):
        # input filenames
        sol_zen_fname = self.input()['sat_sol']['sol_zenith'].path
        sol_azi_fname = self.input()['sat_sol']['sol_azimuth'].path
        slope_fname = self.input()['slp_asp']['slope'].path
        aspect_fname = self.input()['slp_asp']['aspect'].path

        # get the processing tile sizes
        x_tile, y_tile = get_tile_sizes()

        with self.output()['incident'].temporary_path() as incident_fname,\
            self.output()['azi_incident'].temporary_path() as azi_i_fname:

            gaip.incident_angles(sol_zen_fname, sol_azi_fname, slope_fname,
                                 aspect_fname, incident_fname, azi_i_fname,
                                 x_tile, y_tile)


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
        exiting_fname = pjoin(out_path,
                              CONFIG.get('self_shadow',
                                         'exiting_fname'))
        azi_exiting_fname = pjoin(out_path,
                                  CONFIG.get('self_shadow',
                                             'azimuth_exiting_fname'))
        targets = {'exiting': luigi.LocalTarget(exiting_fname),
                   'azi_exiting': luigi.LocalTarget(azi_exiting_fname)}
        return targets

    def run(self):
        # input filenames
        sat_view_fname = self.input()['sat_sol']['sat_view'].path
        sat_azi_fname = self.input()['sat_sol']['sat_azimuth'].path
        slope_fname = self.input()['slp_asp']['slope'].path
        aspect_fname = self.input()['slp_asp']['aspect'].path

        # get the processing tile sizes
        x_tile, y_tile = get_tile_sizes()

        with self.output()['exiting'].temporary_path() as exiting_fname,\
            self.output()['azi_exiting'].temporary_path() as azi_e_fname:

            gaip.exiting_angles(sat_view_fname, sat_azi_fname,
                                slope_fname, aspect_fname, exiting_fname,
                                azi_e_fname, x_tile, y_tile)


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
        out_fname = pjoin(out_path, CONFIG.get('self_shadow',
                                               'relative_slope_fname'))
        return luigi.LocalTarget(out_fname)

    def run(self):
        # input filenames
        azi_incident_fname = self.input()['incident']['azi_incident'].path
        azi_exiting_fname = self.input()['exiting']['azi_exiting'].path

        # get the processing tile sizes
        x_tile, y_tile = get_tile_sizes()

        with self.output().temporary_path() as out_fname:
            gaip.relative_azimuth_slope(azi_incident_fname, azi_exiting_fname,
                                        out_fname, x_tile, y_tile)


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
        incident_fname = self.input()['incident']['incident'].path
        exiting_fname = self.input()['exiting']['exiting'].path

        # get the processing tile sizes
        x_tile, y_tile = get_tile_sizes()

        with self.output().temporary_path() as out_fname:
            gaip.self_shadow(incident_fname, exiting_fname, out_fname,
                             x_tile, y_tile)


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
        dsm_fname = self.input()['dsm']['smoothed'].path
        sol_zen_fname = self.input()['sat_sol']['sol_zenith'].path
        sol_azi_fname = self.input()['sat_sol']['sol_azimuth'].path

        #buffer = int(CONFIG.get('extract_dsm', 'dsm_buffer_width'))
        buffer = get_buffer(self.group)
        window_height = CONFIG.getint('terrain_correction',
                                      'shadow_sub_matrix_height')
        window_width = CONFIG.getint('terrain_correction',
                                     'shadow_sub_matrix_width')

        with self.output().temporary_path() as out_fname:
            gaip.calculate_cast_shadow(acqs[0], dsm_fname, buffer,
                                       window_height, window_width,
                                       sol_zen_fname, sol_azi_fname,
                                       out_fname)


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
        dsm_fname = self.input()['dsm']['smoothed'].path
        sat_view_fname = self.input()['sat_sol']['sat_view'].path
        sat_azi_fname = self.input()['sat_sol']['sat_azimuth'].path

        #buffer = int(CONFIG.get('extract_dsm', 'dsm_buffer_width'))
        buffer = get_buffer(self.group)
        window_height = CONFIG.getint('terrain_correction',
                                      'shadow_sub_matrix_height')
        window_width = CONFIG.getint('terrain_correction',
                                     'shadow_sub_matrix_width')

        with self.output().temporary_path() as out_fname:
            gaip.calculate_cast_shadow(acqs[0], dsm_fname, buffer,
                                       window_height, window_width,
                                       sat_view_fname, sat_azi_fname,
                                       out_fname)


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
                'brdf': GetBrdfAncillaryData(*args[:-1]),
                'dsm': DEMExctraction(*args),
                'rel_slope': RelativeAzimuthSlope(*args),
                'self_shadow': SelfShadow(*args),
                'cast_shadow_sun': CalculateCastShadowSun(*args),
                'cast_shadow_sat': CalculateCastShadowSatellite(*args),
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

        # get the reflectance levels and base output format
        lambertian_fmt = CONFIG.get('terrain_correction', 'lambertian_format')
        brdf_fmt = CONFIG.get('terrain_correction', 'brdf_format')
        terrain_fmt = CONFIG.get('terrain_correction', 'terrain_format')

        lambertian_fname = pjoin(out_path, lambertian_fmt.format(band=band))
        brdf_fname = pjoin(out_path, brdf_fmt.format(band=band))
        terrain_fname = pjoin(out_path, terrain_fmt.format(band=band))

        targets = {'lambertian': luigi.LocalTarget(lambertian_fname),
                   'brdf': luigi.LocalTarget(brdf_fname),
                   'terrain': luigi.LocalTarget(terrain_fname)}

        return targets

    def run(self):
        container = gaip.acquisitions(self.level1)
        acqs = container.get_acquisitions(group=self.group,
                                          granule=self.granule)

        # TODO: what is rori???
        rori = CONFIG.getfloat('terrain_correction', 'rori')

        # inputs
        inputs = self.input()
        bilinear_fname = inputs['bilinear'].path
        brdf_fname = inputs['brdf'].path
        self_shadow_fname = inputs['self_shadow'].path
        slope_fname = inputs['slp_asp']['slope'].path
        aspect_fname = inputs['slp_asp']['aspect'].path
        incident_fname = inputs['incident']['incident'].path
        exiting_fname = inputs['exiting']['exiting'].path
        relative_slope_fname = inputs['rel_slope'].path
        cast_shadow_sun_fname = inputs['cast_shadow_sun'].path
        cast_shadow_sat_fname = inputs['cast_shadow_sat'].path
        solar_zenith_fname = inputs['sat_sol']['sol_zenith'].path
        solar_azimuth_fname = inputs['sat_sol']['sol_azimuth'].path
        satellite_view_fname = inputs['sat_sol']['sat_view'].path
        relative_angle_fname = inputs['sat_sol']['rel_azimuth'].path

        bilinear_data = load_value(bilinear_fname)
        brdf_data = load_value(brdf_fname)

        # get the processing tile sizes
        x_tile, y_tile = get_tile_sizes()

        # get the acquisition we wish to process
        acq = [acq for acq in acqs if acq.band_num == self.band_num][0]


        outputs = self.output()
        with outputs['lambertian'].temporary_path() as lambertian_fname,\
            outputs['brdf'].temporary_path() as brdf_fname,\
            outputs['terrain'].temporary_path() as terrain_fname:

            key1 = (self.band_num, 'lambertian')
            key2 = (self.band_num, 'brdf')
            key3 = (self.band_num, 'terrain')
            reflectance_fnames = {key1: lambertian_fname,
                                  key2: brdf_fname,
                                  key3: terrain_fname}

            # TODO: update reflectance function:
            # accept single acquisition, no string formats, brdf data,
            # no dicts of band numbers etc
            # arguments
            args = [[acq],
                    bilinear_data,
                    rori,
                    brdf_data,
                    self_shadow_fname,
                    cast_shadow_sun_fname,
                    cast_shadow_sat_fname,
                    solar_zenith_fname,
                    solar_azimuth_fname,
                    satellite_view_fname,
                    relative_angle_fname,
                    slope_fname,
                    aspect_fname,
                    incident_fname,
                    exiting_fname,
                    relative_slope_fname,
                    reflectance_fnames,
                    x_tile,
                    y_tile]

            # calculate reflectance
            gaip.calculate_reflectance(*args)


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
        tc_bands = []
        for band in bands_to_process:
            tc_bands.append(RunTCBand(self.level1, self.nbar_root,
                                      self.granule, self.group, band))

        return tc_bands


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
    logging.info('out_path={path}'.format(path=args.work_root))
    logging.info('log_path={path}'.format(path=args.log_path))

    size = int(os.getenv('PBS_NNODES', '1'))
    rank = int(os.getenv('PBS_VNODENUM', '1'))
    main(args.level1_list, args.work_root, size, rank)
