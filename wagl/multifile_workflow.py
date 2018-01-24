#!/usr/bin/env python
"""
Multifile workflow for producing NBAR and SBT
---------------------------------------------

This workflow is geared around a Multiple Independent File workflow, thus
allowing a form a parallelism. HDF5 Linking via a post task then allows
the workflow to appear as if the IO is through a single file.

The multifile workflow approach does allow more freedom of control in
accessing individual coefficients of the entire workflow, and easier for
a user to rapidly test new features.

This workflow can also pick-up exactly where it left off, if the files
generated persist on disk.
However, this method could flood the scheduler if thousands of level1
datasets are submitted at once in a single call.

Workflow settings can be configured in `luigi.cfg` file.
"""
# pylint: disable=missing-docstring,no-init,too-many-function-args
# pylint: disable=too-many-locals
# pylint: disable=protected-access

from __future__ import absolute_import, print_function
import os
from os.path import join as pjoin, basename, dirname, normpath, splitext
from posixpath import join as ppjoin
import logging
import traceback

from structlog import wrap_logger
from structlog.processors import JSONRenderer
import h5py
import luigi
from luigi.local_target import LocalFileSystem
from luigi.util import inherits, requires
from wagl.acquisition import acquisitions
from wagl.ancillary import _collect_ancillary
from wagl.satellite_solar_angles import _calculate_angles
from wagl.incident_exiting_angles import _incident_exiting_angles
from wagl.incident_exiting_angles import _relative_azimuth_slope
from wagl.longitude_latitude_arrays import _create_lon_lat_grids
from wagl.reflectance import _calculate_reflectance, link_standard_data
from wagl.terrain_shadow_masks import _self_shadow, _calculate_cast_shadow
from wagl.terrain_shadow_masks import _combine_shadow
from wagl.slope_aspect import _slope_aspect_arrays
from wagl.constants import Model, BandType, Method, AtmosphericCoefficients
from wagl.constants import POINT_FMT, ALBEDO_FMT, POINT_ALBEDO_FMT, Albedos
from wagl.dsm import _get_dsm
from wagl.modtran import _format_tp5, _run_modtran
from wagl.modtran import _calculate_coefficients, prepare_modtran
from wagl.modtran import link_atmospheric_results
from wagl.interpolation import _interpolate, link_interpolated_data
from wagl.temperature import _surface_brightness_temperature
from wagl.pq import can_pq, _run_pq
from wagl.hdf5 import create_external_link


ERROR_LOGGER = wrap_logger(logging.getLogger('wagl-error'),
                           processors=[JSONRenderer(indent=1, sort_keys=True)])


def get_buffer(group):
    buf = {'product': 250,
           'R10m': 700,
           'R20m': 350,
           'R60m': 120}
    return buf[group]


@luigi.Task.event_handler(luigi.Event.FAILURE)
def on_failure(task, exception):
    """Capture any Task Failure here."""
    ERROR_LOGGER.error(task=task.get_task_family(),
                       params=task.to_str_params(),
                       level1=task.level1,
                       exception=exception.__str__(),
                       traceback=traceback.format_exc().splitlines())


class WorkRoot(luigi.Task):

    """
    Create the work root directory space, and sub directories that
    could compete later in a race condition of creation.
    """

    level1 = luigi.Parameter()
    work_root = luigi.Parameter(significant=False)
    acq_parser_hint = luigi.Parameter(default=None)
    reflectance_dir = '_standardised'
    shadow_dir = '_shadow'
    interpolation_dir = '_interpolation'

    def output(self):
        out_dirs = [self.reflectance_dir, self.shadow_dir,
                    self.interpolation_dir]
        container = acquisitions(self.level1, self.acq_parser_hint)
        for granule in container.granules:
            for group in container.groups:
                pth = container.get_root(self.work_root, group, granule)
                for out_dir in out_dirs:
                    yield luigi.LocalTarget(pjoin(pth, out_dir))

    def run(self):
        local_fs = LocalFileSystem()
        for target in self.output():
            local_fs.mkdir(target.path)


class CalculateLonLatGrids(luigi.Task):

    """Calculates the longitude and latitude grids."""

    level1 = luigi.Parameter()
    work_root = luigi.Parameter(significant=False)
    granule = luigi.Parameter(default=None)
    group = luigi.Parameter()
    acq_parser_hint = luigi.Parameter(default=None)
    compression = luigi.Parameter(default='lzf', significant=False)

    def requires(self):
        # we want to pass the level1 root not the granule root
        root = dirname(self.work_root) if self.granule else self.work_root
        return WorkRoot(self.level1, root)

    def output(self):
        out_path = pjoin(self.work_root, self.group)
        return luigi.LocalTarget(pjoin(out_path, 'longitude-latitude.h5'))

    def run(self):
        acq = (
            acquisitions(self.level1, self.acq_parser_hint)
            .get_acquisitions(self.group, self.granule)
        )[0]

        with self.output().temporary_path() as out_fname:
            _create_lon_lat_grids(acq, out_fname, self.compression)


@inherits(CalculateLonLatGrids)
class CalculateSatelliteAndSolarGrids(luigi.Task):

    """Calculate the satellite and solar grids."""

    tle_path = luigi.Parameter(significant=False)

    def requires(self):
        args = [self.level1, self.work_root, self.granule, self.group]
        return CalculateLonLatGrids(*args)

    def output(self):
        out_path = pjoin(self.work_root, self.group)
        return luigi.LocalTarget(pjoin(out_path, 'satellite-solar.h5'))

    def run(self):
        acqs = (
            acquisitions(self.level1, self.acq_parser_hint)
            .get_acquisitions(self.group, self.granule)
        )

        with self.output().temporary_path() as out_fname:
            _calculate_angles(acqs[0], self.input().path, out_fname,
                              self.compression, self.tle_path)


class AncillaryData(luigi.Task):

    """Get all ancillary data."""

    level1 = luigi.Parameter()
    work_root = luigi.Parameter(significant=False)
    granule = luigi.Parameter(default=None)
    vertices = luigi.TupleParameter()
    model = luigi.EnumParameter(enum=Model)
    acq_parser_hint = luigi.Parameter(default=None)
    aerosol = luigi.DictParameter({'user': 0.05}, significant=False)
    brdf_path = luigi.Parameter(significant=False)
    brdf_premodis_path = luigi.Parameter(significant=False)
    ozone_path = luigi.Parameter(significant=False)
    water_vapour = luigi.DictParameter({'user': 1.5}, significant=False)
    dem_path = luigi.Parameter(significant=False)
    ecmwf_path = luigi.Parameter(significant=False)
    invariant_height_fname = luigi.Parameter(significant=False)
    compression = luigi.Parameter(default='lzf', significant=False)

    def requires(self):
        group = acquisitions(self.level1, self.acq_parser_hint).groups[0]
        args = [self.level1, self.work_root, self.granule, group]
        return CalculateSatelliteAndSolarGrids(*args)

    def output(self):
        return luigi.LocalTarget(pjoin(self.work_root, 'ancillary.h5'))

    def run(self):
        container = acquisitions(self.level1, self.acq_parser_hint))
        grn = container.get_granule(granule=self.granule, container=True)
        sbt_path = None

        nbar_paths = {'aerosol_dict': self.aerosol,
                      'water_vapour_dict': self.water_vapour,
                      'ozone_path': self.ozone_path,
                      'dem_path': self.dem_path,
                      'brdf_path': self.brdf_path,
                      'brdf_premodis_path': self.brdf_premodis_path}

        if self.model == Model.standard or self.model == Model.sbt:
            sbt_path = self.ecmwf_path

        with self.output().temporary_path() as out_fname:
            _collect_ancillary(grn, self.input().path, nbar_paths, sbt_path,
                               self.invariant_height_fname, self.vertices,
                               out_fname, self.compression)


class WriteTp5(luigi.Task):

    """Output the `tp5` formatted files."""

    level1 = luigi.Parameter()
    work_root = luigi.Parameter(significant=False)
    granule = luigi.Parameter(default=None)
    vertices = luigi.TupleParameter()
    acq_parser_hint = luigi.Parameter(default=None)
    model = luigi.EnumParameter(enum=Model)
    base_dir = luigi.Parameter(default='_atmospherics', significant=False)
    compression = luigi.Parameter(default='lzf', significant=False)

    def requires(self):
        container = acquisitions(self.level1, self.acq_parser_hint)
        tasks = {}

        tasks['ancillary'] = AncillaryData(self.level1, self.work_root,
                                           self.granule, self.vertices,
                                           self.model)

        for group in container.groups:
            args = [self.level1, self.work_root, self.granule, group]
            tsks = {'sat_sol': CalculateSatelliteAndSolarGrids(*args),
                    'lon_lat': CalculateLonLatGrids(*args)}
            tasks[group] = tsks

        return tasks

    def output(self):
        out_fname = pjoin(self.work_root, 'atmospheric-inputs.h5')
        return luigi.LocalTarget(out_fname)

    def run(self):
        container = acquisitions(self.level1, self.acq_parser_hint)
        group = container.groups[0]
        acqs = container.get_acquisitions(group, granule=self.granule)

        # output filename format
        output_fmt = pjoin(POINT_FMT, ALBEDO_FMT,
                           ''.join([POINT_ALBEDO_FMT, '.tp5']))

        # input filenames
        ancillary_fname = self.input()['ancillary'].path
        sat_sol_fname = self.input()[group]['sat_sol'].path
        lon_lat_fname = self.input()[group]['lon_lat'].path

        with self.output().temporary_path() as out_fname:
            tp5_data = _format_tp5(acqs, sat_sol_fname, lon_lat_fname,
                                   ancillary_fname, out_fname, self.model)

            # keep this as an indented block, that way the target will remain
            # atomic and be moved upon closing
            for key in tp5_data:
                point, albedo = key
                tp5_fname = output_fmt.format(p=point, a=albedo.value)
                target = pjoin(dirname(out_fname), self.base_dir, tp5_fname)
                with luigi.LocalTarget(target).open('w') as src:
                    src.writelines(tp5_data[key])


@requires(WriteTp5)
class AtmosphericsCase(luigi.Task):

    """
    Run MODTRAN for a specific point (vertex) and albedo.
    This task is parameterised this wat to allow parallel instances
    of MODTRAN to run.
    """

    point = luigi.Parameter()
    albedos = luigi.ListParameter()
    modtran_exe = luigi.Parameter(significant=False)

    def output(self):
        out_path = pjoin(self.work_root, self.base_dir)
        albedos = '-'.join([a for a in self.albedos])
        out_fname = ''.join([POINT_ALBEDO_FMT.format(p=self.point,
                                                     a=albedos), '.h5'])
        return luigi.LocalTarget(pjoin(out_path, out_fname))

    def run(self):
        container = acquisitions(self.level1, self.acq_parser_hint)
        # out_path = container.get_root(self.work_root, granule=self.granule)
        acqs = container.get_acquisitions(granule=self.granule)
        atmospheric_inputs_fname = self.input().path
        base_dir = pjoin(self.work_root, self.base_dir)
        albedos = [Albedos(a) for a in self.albedos]

        prepare_modtran(acqs, self.point, albedos, base_dir, self.modtran_exe)

        with self.output().temporary_path() as out_fname:
            nvertices = self.vertices[0] * self.vertices[1]
            _run_modtran(acqs, self.modtran_exe, base_dir, self.point, albedos,
                         self.model, nvertices, atmospheric_inputs_fname,
                         out_fname, self.compression)


@inherits(WriteTp5)
class Atmospherics(luigi.Task):

    """
    Kicks off MODTRAN calculations for all points and albedos.
    """

    model = luigi.EnumParameter(enum=Model)
    separate = luigi.BoolParameter()

    def requires(self):
        args = [self.level1, self.work_root, self.granule, self.vertices]
        for point in range(self.vertices[0] * self.vertices[1]):
            kwargs = {'point': point, 'model': self.model}
            if self.separate:
                for albedo in self.model.albedos:
                    kwargs['albedos'] = [albedo.value]
                    yield AtmosphericsCase(*args, **kwargs)
            else:
                kwargs['albedos'] = [a.value for a in self.model.albedos]
                yield AtmosphericsCase(*args, **kwargs)

    def output(self):
        out_fname = pjoin(self.work_root, 'atmospheric-results.h5')
        return luigi.LocalTarget(out_fname)

    def run(self):
        nvertices = self.vertices[0] * self.vertices[1]
        with self.output().temporary_path() as out_fname:
            link_atmospheric_results(self.input(), out_fname, nvertices,
                                     self.model)


@requires(Atmospherics)
class CalculateCoefficients(luigi.Task):

    """
    Calculate the atmospheric coefficients needed by BRDF and atmospheric
    correction model.
    """

    def output(self):
        out_fname = pjoin(self.work_root, 'atmospheric-coefficients.h5')
        return luigi.LocalTarget(out_fname)

    def run(self):
        with self.output().temporary_path() as out_fname:
            _calculate_coefficients(self.input().path, out_fname,
                                  self.compression)


@inherits(CalculateLonLatGrids)
class InterpolateCoefficient(luigi.Task):
    """
    Runs the interpolation function for a given band for a
    given atmospheric coefficient.
    """

    vertices = luigi.TupleParameter()
    band_id = luigi.Parameter()
    coefficient = luigi.EnumParameter(enum=AtmosphericCoefficients)
    base_dir = luigi.Parameter(default='_interpolation', significant=False)
    model = luigi.EnumParameter(enum=Model)
    method = luigi.EnumParameter(enum=Method, default=Method.shear)

    def requires(self):
        args = [self.level1, self.work_root, self.granule, self.vertices]
        return {'comp': CalculateCoefficients(*args, model=self.model),
                'satsol': self.clone(CalculateSatelliteAndSolarGrids),
                'ancillary': AncillaryData(*args, model=self.model)}

    def output(self):
        out_path = pjoin(self.work_root, self.group, self.base_dir)
        out_fname = '{}-BAND-{}.h5'.format(self.coefficient.value, self.band_id)
        return luigi.LocalTarget(pjoin(out_path, out_fname))

    def run(self):
        acqs = (
            acquisitions(self.level1, self.acq_parser_hint)
            .get_acquisitions(self.group, self.granule)
        )

        sat_sol_angles_fname = self.input()['satsol'].path
        coefficients_fname = self.input()['comp'].path
        ancillary_fname = self.input()['ancillary'].path

        acq = [acq for acq in acqs if acq.band_id == self.band_id][0]

        with self.output().temporary_path() as out_fname:
            _interpolate(acq, self.coefficient, sat_sol_angles_fname,
                         coefficients_fname, ancillary_fname, out_fname,
                         self.compression, self.method)


@inherits(CalculateLonLatGrids)
class InterpolateCoefficients(luigi.Task):

    """
    Issues InterpolateCoefficient tasks.
    This acts as a helper task, and links the results from each
    InterpolateCoefficient task single HDF5 file.
    """

    vertices = luigi.TupleParameter()
    model = luigi.EnumParameter(enum=Model)
    method = luigi.EnumParameter(enum=Method, default=Method.shear)

    def requires(self):
        container = acquisitions(self.level1, self.acq_parser_hint)
        acqs = container.get_acquisitions(group=self.group,
                                          granule=self.granule)

        # NBAR & SBT acquisitions
        nbar_acqs = [a for a in acqs if a.band_type == BandType.Reflective]
        sbt_acqs = [a for a in acqs if a.band_type == BandType.Thermal]

        tasks = {}
        for coefficient in self.model.atmos_coefficients:
            if coefficient in Model.nbar.atmos_coefficients:
                band_acqs = nbar_acqs
            else:
                band_acqs = sbt_acqs

            for acq in band_acqs:
                key = (acq.band_id, coefficient)
                kwargs = {'level1': self.level1, 'work_root': self.work_root,
                          'granule': self.granule, 'group': self.group,
                          'band_id': acq.band_id, 'coefficient': coefficient,
                          'model': self.model, 'vertices': self.vertices,
                          'method': self.method}
                tasks[key] = InterpolateCoefficient(**kwargs)
        return tasks

    def output(self):
        out_path = pjoin(self.work_root, self.group)
        return luigi.LocalTarget(pjoin(out_path, 'interpolated-coefficients.h5'))

    def run(self):
        fnames = {}
        for key, value in self.input().items():
            fnames[key] = value.path

        with self.output().temporary_path() as out_fname:
            link_interpolated_data(fnames, out_fname)


@inherits(CalculateLonLatGrids)
class DEMExtraction(luigi.Task):

    """
    Extract the DEM covering the acquisition extents plus an
    arbitrary buffer. The subset is then smoothed with a gaussian
    filter.
    """

    dsm_fname = luigi.Parameter(significant=False)

    def requires(self):
        # we want to pass the level1 root not the granule root
        root = dirname(self.work_root) if self.granule else self.work_root
        return WorkRoot(self.level1, root)

    def output(self):
        out_path = pjoin(self.work_root, self.group)
        return luigi.LocalTarget(pjoin(out_path, 'dsm-subset.h5'))

    def run(self):
        acqs = (
            acquisitions(self.level1, self.acq_parser_hint)
            .get_acquisitions(self.group, self.granule)
        )
        margins = get_buffer(self.group)

        with self.output().temporary_path() as out_fname:
            _get_dsm(acqs[0], self.dsm_fname, margins, out_fname,
                     self.compression)


@requires(DEMExtraction)
class SlopeAndAspect(luigi.Task):

    """
    Compute the slope and aspect images.
    """

    def output(self):
        out_path = pjoin(self.work_root, self.group)
        return luigi.LocalTarget(pjoin(out_path, 'slope-aspect.h5'))

    def run(self):
        acqs = (
            acquisitions(self.level1, self.acq_parser_hint)
            .get_acquisitions(self.group, self.granule)
        )
        dsm_fname = self.input().path
        margins = get_buffer(self.group)

        with self.output().temporary_path() as out_fname:
            _slope_aspect_arrays(acqs[0], dsm_fname, margins, out_fname,
                                 self.compression)


@inherits(CalculateLonLatGrids)
class IncidentAngles(luigi.Task):

    """
    Compute the incident angles.
    """

    dsm_fname = luigi.Parameter(significant=False)

    def requires(self):
        args = [self.level1, self.work_root, self.granule, self.group]
        return {'sat_sol': self.clone(CalculateSatelliteAndSolarGrids),
                'slp_asp': SlopeAndAspect(*args, dsm_fname=self.dsm_fname)}

    def output(self):
        out_path = pjoin(self.work_root, self.group)
        return luigi.LocalTarget(pjoin(out_path, 'incident-angles.h5'))

    def run(self):
        # input filenames
        sat_sol_fname = self.input()['sat_sol'].path
        slope_aspect_fname = self.input()['slp_asp'].path

        with self.output().temporary_path() as out_fname:
            _incident_exiting_angles(sat_sol_fname, slope_aspect_fname,
                                     out_fname, self.compression)


@inherits(IncidentAngles)
class ExitingAngles(luigi.Task):

    """
    Compute the exiting angles.
    """

    def requires(self):
        args = [self.level1, self.work_root, self.granule, self.group]
        return {'sat_sol': self.clone(CalculateSatelliteAndSolarGrids),
                'slp_asp': SlopeAndAspect(*args, dsm_fname=self.dsm_fname)}

    def output(self):
        out_path = pjoin(self.work_root, self.group)
        return luigi.LocalTarget(pjoin(out_path, 'exiting-angles.h5'))

    def run(self):
        # input filenames
        sat_sol_fname = self.input()['sat_sol'].path
        slope_aspect_fname = self.input()['slp_asp'].path

        with self.output().temporary_path() as out_fname:
            _incident_exiting_angles(sat_sol_fname, slope_aspect_fname,
                                     out_fname, self.compression, False)


@inherits(IncidentAngles)
class RelativeAzimuthSlope(luigi.Task):

    """
    Compute the relative azimuth angle on the slope surface.
    """

    def requires(self):
        return {'incident': self.clone(IncidentAngles),
                'exiting': self.clone(ExitingAngles)}

    def output(self):
        out_path = (
            acquisitions(self.level1, self.acq_parser_hint)
            .get_root(self.work_root, self.group, self.granule)
        )
        return luigi.LocalTarget(pjoin(out_path, 'relative-slope.h5'))

    def run(self):
        # input filenames
        incident_fname = self.input()['incident'].path
        exiting_fname = self.input()['exiting'].path

        with self.output().temporary_path() as out_fname:
            _relative_azimuth_slope(incident_fname, exiting_fname,
                                    out_fname, self.compression)


@inherits(IncidentAngles)
class SelfShadow(luigi.Task):

    """
    Calculate the self shadow mask.
    """
    base_dir = luigi.Parameter(default='_shadow', significant=False)

    def requires(self):
        return {'incident': self.clone(IncidentAngles),
                'exiting': self.clone(ExitingAngles)}

    def output(self):
        out_path = pjoin(self.work_root, self.group, self.base_dir)
        return luigi.LocalTarget(pjoin(out_path, 'self-shadow.h5'))

    def run(self):
        # input filenames
        incident_fname = self.input()['incident'].path
        exiting_fname = self.input()['exiting'].path

        with self.output().temporary_path() as out_fname:
            _self_shadow(incident_fname, exiting_fname, out_fname,
                         self.compression)


@inherits(SelfShadow)
class CalculateCastShadowSun(luigi.Task):

    """
    Calculates the Cast shadow mask in the direction back to the
    sun.
    """

    def requires(self):
        args = [self.level1, self.work_root, self.granule, self.group]
        return {'sat_sol': self.clone(CalculateSatelliteAndSolarGrids),
                'dsm': DEMExtraction(*args, dsm_fname=self.dsm_fname)}

    def output(self):
        out_path = pjoin(self.work_root, self.group, self.base_dir)
        return luigi.LocalTarget(pjoin(out_path, 'cast-shadow-sun.h5'))

    def run(self):
        acqs = (
            acquisitions(self.level1, self.acq_parser_hint)
            .get_acquisitions(self.group, self.granule)
        )

        # input filenames
        dsm_fname = self.input()['dsm'].path
        sat_sol_fname = self.input()['sat_sol'].path

        # TODO: convert to a func of distance and resolution
        margins = get_buffer(self.group)
        window_height = 500
        window_width = 500

        with self.output().temporary_path() as out_fname:
            _calculate_cast_shadow(acqs[0], dsm_fname, margins, window_height,
                                   window_width, sat_sol_fname, out_fname,
                                   self.compression)


@inherits(SelfShadow)
class CalculateCastShadowSatellite(luigi.Task):

    """
    Calculates the Cast shadow mask in the direction back to the
    sun.
    """

    def requires(self):
        args = [self.level1, self.work_root, self.granule, self.group]
        return {'sat_sol': self.clone(CalculateSatelliteAndSolarGrids),
                'dsm': DEMExtraction(*args, dsm_fname=self.dsm_fname)}

    def output(self):
        out_path = pjoin(self.work_root, self.group, self.base_dir)
        return luigi.LocalTarget(pjoin(out_path, 'cast-shadow-satellite.h5'))

    def run(self):
        acqs = (
            acquisitions(self.level1, self.acq_parser_hint)
            .get_acquisitions(self.group, self.granule)
        )

        # input filenames
        dsm_fname = self.input()['dsm'].path
        sat_sol_fname = self.input()['sat_sol'].path

        # TODO: convert to a func of distance and resolution
        margins = get_buffer(self.group)
        window_height = 500
        window_width = 500

        with self.output().temporary_path() as out_fname:
            _calculate_cast_shadow(acqs[0], dsm_fname, margins, window_height,
                                   window_width, sat_sol_fname, out_fname,
                                   self.compression, False)


@inherits(IncidentAngles)
class CalculateShadowMasks(luigi.Task):

    """
    Issues self and cast shadow tasks for two direction sources;
    the sun and the satellite. Acts as a helper task,
    but combines the results into a single file.
    """

    def requires(self):
        return {'sun': self.clone(CalculateCastShadowSun),
                'sat': self.clone(CalculateCastShadowSatellite),
                'self': self.clone(SelfShadow)}

    def output(self):
        out_path = pjoin(self.work_root, self.group)
        return luigi.LocalTarget(pjoin(out_path, 'shadow-masks.h5'))

    def run(self):
        with self.output().temporary_path() as out_fname:
            inputs = self.input()
            _combine_shadow(inputs['self'].path, inputs['sun'].path,
                            inputs['sat'].path, out_fname, self.compression)


@inherits(InterpolateCoefficients)
class SurfaceReflectance(luigi.Task):

    """Run the terrain correction over a given band."""

    band_id = luigi.Parameter()
    rori = luigi.FloatParameter(default=0.52, significant=False)
    base_dir = luigi.Parameter(default='_standardised', significant=False)
    dsm_fname = luigi.Parameter(significant=False)

    def requires(self):
        reqs = {'interpolation': self.clone(InterpolateCoefficients),
                'ancillary': self.clone(AncillaryData),
                'rel_slope': self.clone(RelativeAzimuthSlope),
                'shadow': self.clone(CalculateShadowMasks),
                'slp_asp': self.clone(SlopeAndAspect),
                'incident': self.clone(IncidentAngles),
                'exiting': self.clone(ExitingAngles),
                'sat_sol': self.clone(CalculateSatelliteAndSolarGrids)}

        return reqs

    def output(self):
        out_path = pjoin(self.work_root, self.group, self.base_dir)
        fname = 'reflectance-{band}.h5'.format(band=self.band_id)
        return luigi.LocalTarget(pjoin(out_path, fname))

    def run(self):
        container = acquisitions(self.level1, self.acq_parser_hint)
        acqs = container.get_acquisitions(self.group, self.granule)

        # inputs
        inputs = self.input()
        interpolation_fname = inputs['interpolation'].path
        slp_asp_fname = inputs['slp_asp'].path
        incident_fname = inputs['incident'].path
        exiting_fname = inputs['exiting'].path
        relative_slope_fname = inputs['rel_slope'].path
        shadow_fname = inputs['shadow'].path
        sat_sol_fname = inputs['sat_sol'].path
        ancillary_fname = inputs['ancillary'].path

        # get the acquisition we wish to process
        acq = [acq for acq in acqs if acq.band_id == self.band_id][0]

        with self.output().temporary_path() as out_fname:
            _calculate_reflectance(acq, acqs, interpolation_fname,
                                   sat_sol_fname, slp_asp_fname,
                                   relative_slope_fname, incident_fname,
                                   exiting_fname, shadow_fname,
                                   ancillary_fname, self.rori, out_fname,
                                   self.compression)


@inherits(SurfaceReflectance)
class SurfaceTemperature(luigi.Task):

    """
    Calculates surface brightness temperature for a given band.
    """

    def requires(self):
        reqs = {'interpolation': self.clone(InterpolateCoefficients),
                'ancillary': self.clone(AncillaryData)}
        return reqs

    def output(self):
        out_path = pjoin(self.work_root, self.group, self.base_dir)
        fname = 'temperature-{band}.h5'.format(band=self.band_id)
        return luigi.LocalTarget(pjoin(out_path, fname))

    def run(self):
        container = acquisitions(self.level1, self.acq_parser_hint)
        acqs = container.get_acquisitions(self.group, self.granule)
        acq = [acq for acq in acqs if acq.band_id == self.band_id][0]

        with self.output().temporary_path() as out_fname:
            interpolation_fname = self.input()['interpolation'].path
            ancillary_fname = self.input()['ancillary'].path
            _surface_brightness_temperature(acq, acqs, interpolation_fname,
                                            ancillary_fname, out_fname,
                                            self.compression)


@inherits(InterpolateCoefficients)
class DataStandardisation(luigi.Task):

    """
    Issues standardisation (analysis ready) tasks for both
    SurfaceReflectance and SurfaceTemperature.
    """
    land_sea_path = luigi.Parameter()
    pixel_quality = luigi.BoolParameter()
    dsm_fname = luigi.Parameter(significant=False)

    def requires(self):
        band_acqs = []
        container = acquisitions(self.level1, self.acq_parser_hint)
        acqs = container.get_acquisitions(group=self.group,
                                          granule=self.granule)

        # NBAR acquisitions
        if self.model == Model.standard or self.model == Model.nbar:
            band_acqs.extend([a for a in acqs if
                              a.band_type == BandType.Reflective])

        # SBT acquisitions
        if self.model == Model.standard or self.model == Model.sbt:
            band_acqs.extend([a for a in acqs if
                              a.band_type == BandType.Thermal])

        tasks = []
        for acq in band_acqs:
            kwargs = {'level1': self.level1, 'work_root': self.work_root,
                      'granule': self.granule, 'group': self.group,
                      'band_id': acq.band_id, 'model': self.model,
                      'vertices': self.vertices, 'method': self.method}
            if acq.band_type == BandType.Thermal:
                tasks.append(SurfaceTemperature(**kwargs))
            else:
                kwargs['dsm_fname'] = self.dsm_fname
                tasks.append(SurfaceReflectance(**kwargs))

        return tasks

    def output(self):
        out_path = pjoin(self.work_root, self.group)
        return luigi.LocalTarget(pjoin(out_path, 'standardised-products.h5'))

    def run(self):
        with self.output().temporary_path() as out_fname:
            fnames = [target.path for target in self.input()]
            link_standard_data(fnames, out_fname)
            sbt_only = self.model == Model.sbt
            if self.pixel_quality and can_pq(self.level1, self.acq_parser_hint) and not sbt_only:
                _run_pq(self.level1, out_fname, self.group, self.land_sea_path,
                        self.compression, self.acq_parser_hint)


class LinkwaglOutputs(luigi.Task):

    """
    Link all the multifile outputs from wagl into a single file.
    """

    level1 = luigi.Parameter()
    work_root = luigi.Parameter()
    granule = luigi.Parameter(default=None)
    acq_parser_hint = luigi.Parameter(default=None)
    model = luigi.EnumParameter(enum=Model)
    vertices = luigi.TupleParameter(default=(5, 5))
    pixel_quality = luigi.BoolParameter()
    method = luigi.EnumParameter(enum=Method, default=Method.shear)
    acq_parser_hint = luigi.Parameter(default=None)
    dsm_fname = luigi.Parameter(significant=False)

    def requires(self):
        container = acquisitions(self.level1, self.acq_parser_hint)
        for group in container.groups:
            kwargs = {'level1': self.level1, 'work_root': self.work_root,
                      'granule': self.granule, 'group': group,
                      'model': self.model, 'vertices': self.vertices,
                      'pixel_quality': self.pixel_quality,
                      'method': self.method, 'dsm_fname': self.dsm_fname}
            yield DataStandardisation(**kwargs)

    def output(self):
        return luigi.LocalTarget('{}.h5'.format(normpath(self.work_root)))

    def run(self):
        with self.output().temporary_path() as out_fname:
            for root, _, files in os.walk(self.work_root):
                # skip any private files
                if basename(root)[0] == '_':
                    continue

                for file_ in files:
                    if splitext(file_)[1] == '.h5':
                        fname = pjoin(root, file_)
                        grp_name = dirname(fname.replace(self.work_root, ''))

                        with h5py.File(fname, 'r') as fid:
                            groups = [g for g in fid]

                        for pth in groups:
                            new_path = ppjoin(grp_name, pth)
                            create_external_link(fname, pth, out_fname,
                                                 new_path)

            with h5py.File(out_fname) as fid:
                container = acquisitions(self.level1, self.acq_parser_hint)
                fid.attrs['level1_uri'] = self.level1
                fid.attrs['tiled'] = container.tiled


class ARD(luigi.WrapperTask):

    """Kicks off ARD tasks for each level1 entry."""

    level1_list = luigi.Parameter()
    outdir = luigi.Parameter()
    model = luigi.EnumParameter(enum=Model)
    vertices = luigi.TupleParameter(default=(5, 5))
    pixel_quality = luigi.BoolParameter()
    method = luigi.EnumParameter(enum=Method, default=Method.shear)
    dsm_fname = luigi.Parameter(significant=False)
    acq_parser_hint = luigi.Parameter(default=None)

    def requires(self):
        with open(self.level1_list) as src:
            level1_list = [level1.strip() for level1 in src.readlines()]

        for level1 in level1_list:
            container = acquisitions(level1, self.acq_parser_hint)
            work_name = '{}.wagl'.format(container.label)
            for granule in container.granules:
                sub_path = granule if granule else ''
                work_root = pjoin(self.outdir, work_name, sub_path)
                kwargs = {'level1': level1, 'work_root': normpath(work_root),
                          'granule': granule, 'model': self.model,
                          'vertices': self.vertices,
                          'pixel_quality': self.pixel_quality,
                          'method': self.method, 'dsm_fname': self.dsm_fname}

                yield LinkwaglOutputs(**kwargs)


class CallTask(luigi.WrapperTask):

    """An entry point for calling most tasks defined in the above
       workflow. Useful for submitting a list of datasets to process
       a given task that could be the entire workflow, or only to
       the desired task.
    """

    level1_list = luigi.Parameter()
    acq_parser_hint = luigi.Parameter(default=None)
    outdir = luigi.Parameter()
    task = luigi.TaskParameter()

    def requires(self):
        with open(self.level1_list) as src:
            level1_list = [level1.strip() for level1 in src.readlines()]

        for level1 in level1_list:
            work_name = '{}.wagl'.format(basename(level1))
            work_root = pjoin(self.outdir, work_name)
            container = acquisitions(level1, self.acq_parser_hint)
            for granule in container.granules:
                if 'group' in self.task.get_param_names():
                    for group in container.groups:
                        yield self.task(level1, work_root, granule, group)
                else:
                    yield self.task(level1, work_root, granule)


if __name__ == '__main__':
    luigi.run()
