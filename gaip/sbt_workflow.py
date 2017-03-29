#!/usr/bin/env python
"""
NBAR Workflow
-------------

Workflow settings can be configured in `luigi.cfg` file.

"""
# pylint: disable=missing-docstring,no-init,too-many-function-args
# pylint: disable=too-many-locals
# pylint: disable=protected-access

from __future__ import absolute_import, print_function
from os.path import join as pjoin, basename, splitext, dirname
import luigi
from luigi.local_target import LocalFileSystem
from luigi.util import inherits, requires
from gaip.acquisition import acquisitions
from gaip.ancillary import collect_ancillary_data, _collect_sbt_ancillary
from gaip.calculate_angles import _calculate_angles
from gaip.calculate_incident_exiting_angles import _incident_angles
from gaip.calculate_incident_exiting_angles import _exiting_angles
from gaip.calculate_incident_exiting_angles import _relative_azimuth_slope
from gaip.calculate_lon_lat_arrays import create_lon_grid, create_lat_grid
from gaip import constants
from gaip.hdf5 import create_external_link
from gaip.modtran import _format_tp5, _run_modtran, _calculate_solar_radiation
from gaip.modtran import _calculate_coefficients, prepare_modtran
from gaip.modtran import POINT_FMT, ALBEDO_FMT, POINT_ALBEDO_FMT
from gaip.interpolation import _bilinear_interpolate, link_bilinear_data
from gaip.nbar_workflow import CalculateLonGrid, CalculateLatGrid, WorkRoot
from gaip.nbar_workflow import CalculateSatelliteAndSolarGrids
from gaip.nbar_workflow import RunModtranCase, WriteTp5, GetAncillaryData


# TODO: Remove the dpendency on the satelite solar grids.
#       This means moving the coordinator calculator out of sat sol grids.
#       It'll clean up the workflow and streamline it a little

@requires(CalculateSatelliteAndSolarGrids)
class SBTAncillary(luigi.Task):

    """Collect the ancillary data required for SBT."""

    vertices = luigi.TupleParameter(default=(5, 5), significant=False)
    dewpoint_path = luigi.Parameter(significant=False)
    temp_2m_path = luigi.Parameter(significant=False)
    surface_pressure_path = luigi.Parameter(significant=False)
    geopotential_path = luigi.Parameter(significant=False)
    temperature_path = luigi.Parameter(significant=False)
    relative_humidity_path = luigi.Parameter(significant=False)
    invariant_height_fname = luigi.Parameter(significant=False)

    def output(self):
        out_path = acquisitions(self.level1).get_root(self.work_root,
                                                      granule=self.granule)
        return luigi.LocalTarget(pjoin(out_path, 'sbt-ancillary.h5'))

    def run(self):
        container = acquisitions(self.level1)
        acqs = container.get_acquisitions(granule=self.granule)
        work_root =  container.get_root(self.work_root, granule=self.granule)

        sat_sol_fname = self.input().path

        with self.output().temporary_path() as out_fname:
            _collect_sbt_ancillary(acqs[0], sat_sol_fname, self.dewpoint_path,
                                   self.temp_2m_path,
                                   self.surface_pressure_path,
                                   self.geopotential_path,
                                   self.temperature_path,
                                   self.relative_humidity_path,
                                   self.invariant_height_fname,
                                   out_fname, self.compression)


# TODO: Remove albedos as a command option.
#       We should default to calculate all.
#       We could still put a switch to say 3 or 4 modtran calls per point??

class ThermalTp5(WriteTp5):

    """Output the `tp5` formatted files."""

    vertices = luigi.TupleParameter(default=(5, 5), significant=False)
    albedos = luigi.ListParameter(default=['th'], significant=False)

    def requires(self):
        # for consistancy, we'll wait for dependencies on all granules and
        # groups of acquisitions
        # current method requires to compute an average from all granules
        # if the scene is tiled up that way
        container = acquisitions(self.level1)
        tasks = {}

        for granule in container.granules:
            args1 = [self.level1, self.work_root, granule]
            kwargs = {'level1': self.level1,
                      'work_root': self.work_root,
                      'granule': granule,
                      'group': container.groups[0],
                      'vertices': self.vertices}
            tasks[(granule, 'sbt-ancillary')] = SBTAncillary(**kwargs)
            tasks[(granule, 'ancillary')] = GetAncillaryData(*args1)
            for group in container.groups:
                args2 = [self.level1, self.work_root, granule, group]
                kwargs['group'] = group
                tsks = {'sat_sol': CalculateSatelliteAndSolarGrids(**kwargs),
                        'lat': CalculateLatGrid(*args2),
                        'lon': CalculateLonGrid(*args2)}
                tasks[(granule, group)] = tsks

        return tasks


class SBTAccumulateSolarIrradiance(AccumulateSolarIrradiance):

    """
    Accumulate the solar irradiance specifically for the SBT workflow.
    """

    albedos = luigi.ListParameter(default=['th'], significant=False)


@requires(SBTAccumulateSolarIrradiance)
class SBTCoefficients(CalculateCoefficients):


class SBTBilinearInterpolation(BilinearInterpolation):

    # TODO: which factors to use

    factors = luigi.ListParameter(default=[], significant=False)

    def requires(self):
        container = acquisitions(self.level1)
        acqs = container.get_acquisitions(group=self.group,
                                          granule=self.granule)

        # Retrieve the satellite and sensor for the acquisition
        satellite = acqs[0].spacecraft_id
        sensor = acqs[0].sensor_id

        # Get the required thermal bands list for processing
        bands_to_process = constants.sbt_bands(satellite, sensor)
        bands = [a.band_num for a in acqs if a.band_num in bands_to_process]

        tasks = {}
        for factor in self.factors:
            for band in bands:
                key = (band, factor)
                kwargs = {'level1': self.level1, 'work_root': self.work_root,
                          'granule': self.granule, 'group': self.group,
                          'band_num': band, 'factor': factor}
                tasks[key] = BilinearInterpolationBand(**kwargs)
        return tasks


class SurfaceTemperature(luigi.Task):

    """Calculates surface brightness temperature."""

    # TODO: Re-write the class contents (output, run)

    def requires(self):

        # TODO: add other upstream tasks as needed

        args = [self.level1, self.work_root, self.granule, self.group]
        reqs = {'bilinear': SBTBilinearInterpolation(*args),
                'sat_sol': self.clone(CalculateSatelliteAndSolarGrids)}
        return reqs


class SBT(luigi.WrapperTask):

    """Kicks off SurfaceTemperature tasks for each level1 entry."""

    level1_csv = luigi.Parameter()
    output_directory = luigi.Parameter()
    work_extension = luigi.Parameter(default='.gaip-work', significant=False)

    def requires(self):
        with open(self.level1_csv) as src:
            level1_scenes = [scene.strip() for scene in src.readlines()]

        for scene in level1_scenes:
            work_name = basename(scene) + self.work_extension
            work_root = pjoin(self.output_directory, work_name)
            container = acquisitions(scene)
            for granule in container.granules:
                for group in container.groups:
                    yield SurfaceTemperature(scene, work_root, granule, group)


if __name__ == '__main__':
    luigi.run()
