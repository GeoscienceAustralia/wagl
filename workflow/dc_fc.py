#!/bin/env python

import luigi.contrib.mpi as mpi
import luigi
import argparse
import os
from os.path import join as pjoin, dirname, exists, splitext, basename
import logging
import gaip

from datacube.api.model import DatasetType, Satellite
from datacube.api.query import list_tiles_as_list
from datacube.config import Config
from datetime import date

CONFIG = luigi.configuration.get_config()
CONFIG.add_config_path(pjoin(dirname(__file__), 'fc.cfg'))

class CreateDirs(luigi.Task):
    """
    Create the output directory.
    """

    out_path = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        out_path = self.out_path
        return luigi.LocalTarget(out_path)

    def run(self):
        out_path = self.out_path
        if not exists(out_path):
            os.makedirs(out_path)


class FractionalCoverTask(luigi.Task):
    nbar_ds = luigi.Parameter()
    out_dir = luigi.Parameter()

    def requires(self):
        #return [CreateDirs(self.out_dir)]
        return []

    def output(self):
        base_name = splitext(basename(self.nbar_ds.path))[0]
        base_name = base_name.replace('NBAR', 'FC') + '.tif'
        out_fname = pjoin(self.out_dir, base_name)
        return luigi.LocalTarget(out_fname)

    def run(self):
        base_name = splitext(basename(self.nbar_ds.path))[0]
        base_name = base_name.replace('NBAR', 'FC') + '.tif'
        out_fname = pjoin(self.out_dir, base_name)

        # Get the processing tile sizes
        x_tile = int(CONFIG.get('work', 'x_tile_size'))
        y_tile = int(CONFIG.get('work', 'y_tile_size'))
        x_tile = None if x_tile <= 0 else x_tile
        y_tile = None if y_tile <= 0 else y_tile

        # Run
        gaip.fractional_cover(self.nbar_ds, x_tile, y_tile, out_fname)


if __name__ == '__main__':
    # command line arguments

    parser = argparse.ArgumentParser()
    parser.add_argument("--debug", help="selects more detail logging (default is INFO)", \
        default=False, action='store_true')

    args = parser.parse_args()

    # setup logging
    
    logfile = "run_fc_{}_{}.log".format(os.uname()[1], os.getpid())
    logfile = os.path.join(args.log_path, logfile)
    logging_level = logging.INFO
    if args.debug:
        logging_level = logging.DEBUG
    logging.basicConfig(filename=logfile, level=logging_level, \
        format= '[%(asctime)s] {%(pathname)s:%(lineno)d} %(levelname)s - %(message)s',
        datefmt='%H:%M:%S')
    logging.info("dc_fc.py started")

    # Define the database config, dataset type and satellites
    config = Config()
    ds_type = DatasetType.ARG25
    satellites = [Satellite.LS7, Satellite.LS5, Satellite.LS8]
    tiles = list_tiles_as_list(x=[140], y=[-35], acq_min=date(2000, 1, 1),
                               acq_max=date(2010, 12, 31),
                               satellites=satellites, datasets=ds_type,
                               database=config.get_db_database(),
                               user=config.get_db_username(),
                               password=config.get_db_password(),
                               host=config.get_db_host(),
                               port=config.get_db_port())

    cell_out_dir = '{cell_x}_{cell_y}'.format(cell_x=tiles[0].x,
                                              cell_y=tiles[0].y)

    out_dir = CONFIG.get('work', 'output_directory')
    out_dir = pjoin(out_dir, cell_out_dir)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Create the task list
    tasks = []
    for ds in tiles:
        reflectance_ds = ds[ds_type]
        tasks.append(FractionalCoverTask(reflectance_ds, out_dir))

    mpi.run(tasks)
