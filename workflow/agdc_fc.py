#!/bin/env python

import luigi.contrib.mpi as mpi

import luigi
import os
from os.path import join as pjoin, dirname, exists, splitext, basename
import gaip
import cPickle as pickle

from datacube.api.model import DatasetType, Satellite
from datacube.api.query import list_tiles_as_list
from datacube.config import Config
from datetime import date
from EOtools.DatasetDrivers.stacked_dataset import StackedDataset

CONFIG = luigi.configuration.get_config()
CONFIG.add_config_path(pjoin(dirname(__file__), 'fc.cfg'))

class TileQuery(luigi.Task):
    """
    Queries the database based on the input parameters.
    """
    out_path = luigi.Parameter()
    satellites = luigi.Parameter()
    min_date = luigi.DateParameter()
    max_date = luigi.DateParameter()

    config = Config()
    ds_type = DatasetType.ARG25

    def requires(self):
        return []

    def output(self):
        out_file = pjoin(self.out_path, CONFIG.get('agdc', 'query_result'))
        return luigi.LocalTarget(out_file)

    def run(self):
        print "Executing DB query!"
        satellites = [Satellite(i) for i in self.satellites.split(',')]
        tiles = list_tiles_as_list(x=[140], y=[-35], acq_min=self.min_date,
                                   acq_max=self.max_date,
                                   satellites=satellites,
                                   datasets=self.ds_type,
                                   database=self.config.get_db_database(),
                                   user=self.config.get_db_username(),
                                   password=self.config.get_db_password(),
                                   host=self.config.get_db_host(),
                                   port=self.config.get_db_port())

        with self.output().open('w') as outf:
            pickle.dump(tiles, outf)


class FractionalCoverTask(luigi.Task):
    """
    Computes fractional cover for a given input.
    """
    fname = luigi.Parameter()
    out_dir = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        base_name = splitext(basename(self.fname))[0]
        base_name = base_name.replace('NBAR', 'FC') + '.tif'

        # QDERM TEST
        #base_name = base_name.replace('dbgm4', 'FC') + '_byte.tif'

        out_fname = pjoin(self.out_dir, base_name)
        return luigi.LocalTarget(out_fname)

    def run(self):
        out_fname = self.output().path

        # Get the processing tile sizes
        x_tile = int(CONFIG.get('work', 'x_tile_size'))
        y_tile = int(CONFIG.get('work', 'y_tile_size'))
        x_tile = None if x_tile <= 0 else x_tile
        y_tile = None if y_tile <= 0 else y_tile

        # Run
        sd = StackedDataset(self.fname)
        gaip.fractional_cover(sd, x_tile, y_tile, out_fname)


class ProcessFC(luigi.Task):
    """
    Creates FractionalCoverTask's based on the results determined from
    the database query.
    """
    out_path = luigi.Parameter()

    ds_type = DatasetType.ARG25

    def requires(self):
        query_result = pjoin(self.out_path, CONFIG.get('agdc', 'query_result'))
        with open(query_result, 'r') as f:
            query = pickle.load(f)

        tasks = []
        for ds in query:
            fname = ds.datasets[self.ds_type].path

            cell_out_dir = '{cell_x}_{cell_y}'.format(cell_x=ds.x,
                                                      cell_y=ds.y)

            out_dir = pjoin(self.out_path, cell_out_dir)
            if not exists(out_dir):
                os.makedirs(out_dir)

            tasks.append(FractionalCoverTask(fname, out_dir))

        return tasks

    def output(self):
        out_fname = pjoin(self.out_path, 'ProcessFC.completed')
        return luigi.LocalTarget(out_fname)

    def run(self):
        with self.output().open('w') as outf:
            outf.write('Completed')


if __name__ == '__main__':
    out_dir = CONFIG.get('agdc', 'output_directory')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    satellites = CONFIG.get('agdc', 'satellites')

    min_date = CONFIG.get('agdc', 'min_date')
    min_date = [int(i) for i in min_date.split('_')]
    min_date = date(min_date[0], min_date[1], min_date[2])

    max_date = CONFIG.get('agdc', 'max_date')
    max_date = [int(i) for i in max_date.split('_')]
    max_date = date(max_date[0], max_date[1], max_date[2])

    tasks = [TileQuery(out_dir, satellites, min_date, max_date)]
    mpi.run(tasks)
    tasks = [ProcessFC(out_dir)]
    mpi.run(tasks)




    # *********************** QDERM TEST ******************************
    #task = [FractionalCoverTask('/g/data1/v10/testing_ground/jps547/FC/QDERM_test_data/l5tmre_p095r084_20090411/l5tmre_p095r084_20090411_dbgm4.img',
    #                            '/g/data1/v10/testing_ground/jps547/FC/QDERM_test_data/l5tmre_p095r084_20090411')]
    #mpi.run(task)
