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
from gaip.ancillary import _collect_ancillary
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


class SBTAncillary(GetAncillaryData):

    """Collect the ancillary data required for SBT."""

    vertices = luigi.TupleParameter(default=(5, 5), significant=False)
    dewpoint_path = luigi.Parameter(significant=False)
    temp_2m_path = luigi.Parameter(significant=False)
    surface_pressure_path = luigi.Parameter(significant=False)
    geopotential_path = luigi.Parameter(significant=False)
    temperature_path = luigi.Parameter(significant=False)
    relative_humidity_path = luigi.Parameter(significant=False)
    invariant_height_fname = luigi.Parameter(significant=False)

    def run(self):
        container = acquisitions(self.level1)
        acqs = container.get_acquisitions(granule=self.granule)
        work_root =  container.get_root(self.work_root, granule=self.granule)

        nbar_paths = {'aerosol_fname': self.aerosol_fname,
                      'water_vapour_path': self.water_vapour_path,
                      'ozone_path': self.ozone_path,
                      'dem_path': self.dem_path,
                      'brdf_path': self.brdf_path,
                      'brdf_premodis_path': self.brdf_premodis_path}

        sbt_paths = {'dewpoint_path': self.dewpoint_path,
                     'temperature_2m_path': self.temp_2m_path,
                     'surface_pressure_path': self.surface_pressure_path,
                     'geopotential_path': self.geopotential_path,
                     'temperature_path': self.temperature_path,
                     'relative_humidity_path': self.relative_humidity_path,
                     'invariant_fname': self.invariant_height_fname}

        with self.output().temporary_path() as out_fname:
            _collect_ancillary(acq, self.input().path, nbar_paths, sbt_paths,
                               vertices=self.vertices, out_fname=out_fname,
                               work_path=work_root,
                               compression=self.compression)


class ThermalTp5(WriteTp5):

    """Output the `tp5` formatted files."""

    def requires(self):
        container = acquisitions(self.level1)
        tasks = {}

        for granule in container.granules:
            args1 = [self.level1, self.work_root, granule]
            tasks[(granule, 'ancillary')] = SBTAncillary(*args1)
            for group in container.groups:
                args2 = [self.level1, self.work_root, granule, group]
                tsks = {'sat_sol': CalculateSatelliteAndSolarGrids(*args2),
                        'lat': CalculateLatGrid(*args2),
                        'lon': CalculateLonGrid(*args2)}
                tasks[(granule, group)] = tsks

        return tasks


class SBTAccumulateSolarIrradiance(AccumulateSolarIrradiance):

    """
    Accumulate the solar irradiance specifically for the SBT workflow.
    """


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
