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
thousands of level1 datasets to be submitted to the scheduler at once.

Workflow settings can be configured in `luigi.cfg` file.
"""
# pylint: disable=missing-docstring,no-init,too-many-function-args
# pylint: disable=too-many-locals
# pylint: disable=protected-access

from os.path import join as pjoin, basename
import traceback
import luigi
from luigi.util import inherits

from wagl.acquisition import acquisitions
from wagl.constants import Workflow, Method
from wagl.hdf5 import H5CompressionFilter
from wagl.standardise import card4l

from wagl.logs import TASK_LOGGER


@luigi.Task.event_handler(luigi.Event.FAILURE)
def on_failure(task, exception):
    """Capture any Task Failure here."""
    TASK_LOGGER.exception(task=task.get_task_family(),
                          params=task.to_str_params(),
                          level1=getattr(task, 'level1', ''),
                          stack_info=True,
                          status='failure',
                          exception=exception.__str__(),
                          traceback=traceback.format_exc().splitlines())


@luigi.Task.event_handler(luigi.Event.SUCCESS)
def on_success(task):
    """Capture any Task Success here."""
    TASK_LOGGER.info(task=task.get_task_family(),
                     params=task.to_str_params(),
                     level1=getattr(task, 'level1', ''),
                     status='success')


class DataStandardisation(luigi.Task):

    """
    Runs the standardised product workflow.
    """
    level1 = luigi.Parameter()
    outdir = luigi.Parameter()
    granule = luigi.OptionalParameter(default='')
    workflow = luigi.EnumParameter(enum=Workflow, default=Workflow.STANDARD)
    vertices = luigi.TupleParameter(default=(5, 5))
    method = luigi.EnumParameter(enum=Method, default=Method.SHEAR)
    pixel_quality = luigi.BoolParameter()
    land_sea_path = luigi.Parameter()
    aerosol = luigi.DictParameter(default={'user': 0.05}, significant=False)
    brdf_path = luigi.Parameter(significant=False)
    brdf_premodis_path = luigi.Parameter(significant=False)
    ozone_path = luigi.Parameter(significant=False)
    water_vapour = luigi.DictParameter(default={'user': 1.5},
                                       significant=False)
    dem_path = luigi.Parameter(significant=False)
    ecmwf_path = luigi.Parameter(significant=False)
    invariant_height_fname = luigi.Parameter(significant=False)
    dsm_fname = luigi.Parameter(significant=False)
    modtran_exe = luigi.Parameter(significant=False)
    tle_path = luigi.Parameter(significant=False)
    rori = luigi.FloatParameter(default=0.52, significant=False)
    compression = luigi.EnumParameter(enum=H5CompressionFilter,
                                      default=H5CompressionFilter.LZF,
                                      significant=False)
    filter_opts = luigi.DictParameter(default=None, significant=False)
    acq_parser_hint = luigi.OptionalParameter(default='')
    buffer_distance = luigi.FloatParameter(default=8000, significant=False)
    h5_driver = luigi.OptionalParameter(default='', significant=False)

    def output(self):
        fmt = '{label}.wagl.h5'
        label = self.granule if self.granule else basename(self.level1)
        out_fname = fmt.format(label=label)
         
        return luigi.LocalTarget(pjoin(self.outdir, out_fname))

    def run(self):
        if self.workflow == Workflow.STANDARD or self.workflow == Workflow.SBT:
            ecmwf_path = self.ecmwf_path
        else:
            ecmwf_path = None

        with self.output().temporary_path() as out_fname:
            card4l(self.level1, self.granule, self.workflow, self.vertices,
                   self.method, self.pixel_quality, self.land_sea_path,
                   self.tle_path, self.aerosol, self.brdf_path,
                   self.brdf_premodis_path, self.ozone_path, self.water_vapour,
                   self.dem_path, self.dsm_fname, self.invariant_height_fname,
                   self.modtran_exe, out_fname, ecmwf_path, self.rori,
                   self.buffer_distance, self.compression, self.filter_opts,
                   self.h5_driver, self.acq_parser_hint)


@inherits(DataStandardisation)
class ARD(luigi.WrapperTask):

    """Kicks off ARD tasks for each level1 entry."""

    level1_list = luigi.Parameter()

    # override here so it's not required at the command line or config
    level1 = luigi.OptionalParameter(default='', significant=False)

    def requires(self):
        with open(self.level1_list) as src:
            level1_list = [level1.strip() for level1 in src.readlines()]

        for level1 in level1_list:
            container = acquisitions(level1)
            outdir = pjoin(self.outdir, '{}.wagl'.format(container.label))
            for granule in container.granules:
                kwargs = {'level1': level1,
                          'granule': granule,
                          'workflow': self.workflow,
                          'vertices': self.vertices,
                          'pixel_quality': self.pixel_quality,
                          'method': self.method,
                          'modtran_exe': self.modtran_exe,
                          'outdir': outdir,
                          'land_sea_path': self.land_sea_path,
                          'aerosol': self.aerosol,
                          'brdf_path': self.brdf_path,
                          'brdf_premodis_path': self.brdf_premodis_path,
                          'ozone_path': self.ozone_path,
                          'water_vapour': self.water_vapour,
                          'dem_path': self.dem_path,
                          'ecmwf_path': self.ecmwf_path,
                          'invariant_height_fname': self.invariant_height_fname,
                          'dsm_fname': self.dsm_fname,
                          'tle_path': self.tle_path,
                          'rori': self.rori,
                          'compression': self.compression,
                          'filter_opts': self.filter_opts,
                          'buffer_distance': self.buffer_distance,
                          'h5_driver': self.h5_driver}
                yield DataStandardisation(**kwargs)

        
if __name__ == '__main__':
    luigi.run()
