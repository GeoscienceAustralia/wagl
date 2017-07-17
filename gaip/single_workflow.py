#!/usr/bin/env python
"""
Single file workflow for producing NBAR and SBT
-----------------------------------------------

Workflow settings can be configured in `luigi.cfg` file.
"""
# pylint: disable=missing-docstring,no-init,too-many-function-args
# pylint: disable=too-many-locals
# pylint: disable=protected-access

from os.path import join as pjoin, basename, dirname
import logging
import traceback
import luigi
from luigi.util import inherits, requires

from gaip.constants import Model
from gaip.standardise import card4l

ERROR_LOGGER = logging.getLogger('luigi-error')


def get_buffer(group):
    buf = {'product': 250,
           'R10m': 700,
           'R20m': 350,
           'R60m': 120}
    return buf[group]


@luigi.Task.event_handler(luigi.Event.FAILURE)
def on_failure(task, exception):
    """Capture any Task Failure here."""
    fmt = "Error processing scene:\n{}\npath:\n{}"
    msg = fmt.format(basename(task.level1), task.level1)
    excp_msg = exception.__str__()
    traceback_msg = traceback.format_exc()
    ERROR_LOGGER.error(msg)
    ERROR_LOGGER.error(excp_msg)
    ERROR_LOGGER.error(traceback_msg)


class Standard(luigi.Task):

    """
    Runs the standardised product workflow.
    """
    level1 = luigi.Parameter()
    outdir = luigi.Parameter()
    model = luigi.EnumParameter(enum=Model)
    vertices = luigi.TupleParameter(default=(5, 5))
    method = luigi.Parameter(default='shear')
    pixel_quality = luigi.BoolParameter()
    land_sea_path = luigi.Parameter()
    aerosol_fname = luigi.Parameter(significant=False)
    brdf_path = luigi.Parameter(significant=False)
    brdf_premodis_path = luigi.Parameter(significant=False)
    ozone_path = luigi.Parameter(significant=False)
    water_vapour_path = luigi.Parameter(significant=False)
    dem_path = luigi.Parameter(significant=False)
    ecmwf_path = luigi.Parameter(significant=False)
    invariant_height_fname = luigi.Parameter(significant=False)

    def output(self):
        fmt = '{scene}_{model}.h5'
        scene = basename(self.level1)
        out_fname = fmt.format(scene=scene, model=self.model.name)
        return luigi.LocalTarget(pjoin(self.outdir, out_fname))

    def run(self):
        with self.output().temporary_path() as out_fname:
            card4l(self.level1, self.model, self.vertices, self.method,
                   self.pixel_quality, self.land_sea_path, self.ecmwf_path, 
                   self.tle_path, self.aerosol_fname, self.brdf_path,
                   self.brdf_premodis_path, self.ozone_path,
                   self.water_vapour_path, self.dsm_path,
                   self.invariant_height_fname, self.modtran_exe, out_fname,
                   self.rori, self.compression, self.y_tile)


@inherits(Standard)
class ARD(luigi.WrapperTask):

    """Kicks off ARD tasks for each level1 entry."""

    level1_list = luigi.Parameter()

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
                      'aerosol_fname': self.aerosol_fname,
                      'brdf_path': self.brdf_path,
                      'brdf_premodis_path': self.brdf_premodis_path,
                      'ozone_path': self.ozone_path,
                      'water_vapour_path': self.water_vapour_path,
                      'dem_path': self.dem_path,
                      'ecmwf_path': self.ecmwf_path,
                      'invariant_height_fname': self.invariant_height_fname}
            yield Standard(**kwargs)

        
if __name__ == '__main__':
    luigi.run()
