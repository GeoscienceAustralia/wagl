#!/usr/bin/env python
"""
Single file workflow for producing NBAR and SBT
-----------------------------------------------

This workflow is geared to minimise the number of files on disk
and provide a kind of direct to archive compute, and retain all
the necessary intermediate files, which comprise a mixture of
imagery, tables, and point/scalar datasets.

It also provides a consistant logical structure allowing an easier
comparison between 'archives' from different production runs, or
versions of wagl.

This workflow is more suited to full production runs, where testing
has ensured that the workflow is sound, and more easilt allows
thousands of scenes to be submitted to the scheduler at once.

Workflow settings can be configured in `luigi.cfg` file.
"""
# pylint: disable=missing-docstring,no-init,too-many-function-args
# pylint: disable=too-many-locals
# pylint: disable=protected-access

from os.path import join as pjoin, basename
import logging
import traceback
from structlog import wrap_logger
from structlog.processors import JSONRenderer
import luigi

from wagl.constants import Model, Method
from wagl.standardise import card4l


ERROR_LOGGER = wrap_logger(logging.getLogger('wagl-error'),
                           processors=[JSONRenderer(indent=1, sort_keys=True)])
INTERFACE_LOGGER = logging.getLogger('luigi-interface')


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
                       scene=task.level1,
                       exception=exception.__str__(),
                       traceback=traceback.format_exc().splitlines())


class DataStandardisation(luigi.Task):

    """
    Runs the standardised product workflow.
    """
    level1 = luigi.Parameter()
    outdir = luigi.Parameter()
    model = luigi.EnumParameter(enum=Model)
    vertices = luigi.TupleParameter(default=(5, 5))
    method = luigi.EnumParameter(enum=Method, default=Method.shear)
    pixel_quality = luigi.BoolParameter()
    land_sea_path = luigi.Parameter()
    aerosol = luigi.DictParameter({'user': 0.05}, significant=False)
    brdf_path = luigi.Parameter(significant=False)
    brdf_premodis_path = luigi.Parameter(significant=False)
    ozone_path = luigi.Parameter(significant=False)
    water_vapour = luigi.DictParameter({'user': 1.5}, significant=False)
    dem_path = luigi.Parameter(significant=False)
    ecmwf_path = luigi.Parameter(significant=False)
    invariant_height_fname = luigi.Parameter(significant=False)
    dsm_fname = luigi.Parameter(significant=False)
    modtran_exe = luigi.Parameter(significant=False)
    tle_path = luigi.Parameter(significant=False)
    rori = luigi.FloatParameter(default=0.52, significant=False)
    compression = luigi.Parameter(default='lzf', significant=False)
    acq_parser_hint = luigi.Parameter(default=None)

    def output(self):
        fmt = '{scene}.wagl.h5'
        scene = basename(self.level1)
        out_fname = fmt.format(scene=scene, model=self.model.name)
        return luigi.LocalTarget(pjoin(self.outdir, out_fname))

    def run(self):
        if self.model == Model.standard or self.model == Model.sbt:
            ecmwf_path = self.ecmwf_path
        else:
            ecmwf_path = None

        with self.output().temporary_path() as out_fname:
            card4l(self.level1, self.model, self.vertices, self.method,
                   self.pixel_quality, self.land_sea_path, self.tle_path,
                   self.aerosol, self.brdf_path, self.brdf_premodis_path,
                   self.ozone_path, self.water_vapour, self.dem_path,
                   self.dsm_fname, self.invariant_height_fname,
                   self.modtran_exe, out_fname, ecmwf_path, self.rori,
                   self.compression, self.acq_parser_hint)


@inherits(DataStandardisation)
class ARD(luigi.WrapperTask):

    """Kicks off ARD tasks for each level1 entry."""

    level1_list = luigi.Parameter()

    # override so it's not required at the command line
    level1 = luigi.Parameter(default=None, significant=False)

    def requires(self):
        with open(self.level1_list) as src:
            level1_scenes = [scene.strip() for scene in src.readlines()]

        for scene in level1_scenes:
            kwargs = {'level1': scene,
                      'model': self.model,
                      'vertices': self.vertices,
                      'pixel_quality': self.pixel_quality,
                      'method': self.method,
                      'modtran_exe': self.modtran_exe,
                      'outdir': self.outdir,
                      'land_sea_path': self.land_sea_path,
                      'aerosol': self.aerosol,
                      'brdf_path': self.brdf_path,
                      'brdf_premodis_path': self.brdf_premodis_path,
                      'ozone_path': self.ozone_path,
                      'water_vapour': self.water_vapour,
                      'dem_path': self.dem_path,
                      'ecmwf_path': ecmwf_path,
                      'invariant_height_fname': self.invariant_height_fname,
                      'dsm_fname': self.dsm_fname,
                      'tle_path': self.tle_path,
                      'rori': self.rori}
            yield DataStandardisation(**kwargs)

        
if __name__ == '__main__':
    luigi.run()
