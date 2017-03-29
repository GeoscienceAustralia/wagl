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

    # def run(self):
    #     container = acquisitions(self.level1)
    #     # as we have an all granules groups dependency, it doesn't matter which
    #     # group, so just get the first and use it to retrieve the angles
    #     group = container.groups[0]
    #     acq = container.get_acquisitions(group, granule=self.granule)[0]

    #     # input data files, and the output format
    #     inputs = self.input()
    #     output_fmt = pjoin(POINT_FMT, ALBEDO_FMT,
    #                        ''.join([POINT_ALBEDO_FMT, '.tp5']))

    #     # all ancillary filenames from each granule
    #     fnames = [inputs[key].path for key in inputs if 'ancillary' in key]

    #     if container.tiled:
    #         ancillary_fname = pjoin(self.work_root, 'averaged-ancillary.h5')
    #         aggregate_ancillary(fnames, ancillary_fname)
    #     else:
    #         ancillary_fname = fnames[0]

    #     sat_sol_fname = inputs[(self.granule, group)]['sat_sol'].path
    #     lon_fname = inputs[(self.granule, group)]['lon'].path
    #     lat_fname = inputs[(self.granule, group)]['lat'].path
    #     sbt_ancillary_fname = inputs[(self.granule, 'sbt-ancillary')].path

    #     with self.output().temporary_path() as out_fname:
    #         tp5_data = _format_tp5(acq, sat_sol_fname, lon_fname, lat_fname,
    #                                ancillary_fname, out_fname, self.albedos,
    #                                sbt_ancillary_fname)

    #         # keep this as an indented block, that way the target will remain
    #         # atomic and be moved upon closing
    #         for key in tp5_data:
    #             point, albedo = key
    #             tp5_out_fname = output_fmt.format(p=point, a=albedo)
    #             target = pjoin(dirname(out_fname), tp5_out_fname)
    #             with luigi.LocalTarget(target).open('w') as src:
    #                 src.writelines(tp5_data[key])


@requires(ThermalTp5)
class SBTModtranCase(RunModtranCase):
    pass


class SurfaceTemperature(luigi.Task):
    """Calculates surface brightness temperature."""


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
