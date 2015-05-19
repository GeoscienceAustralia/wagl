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

        # Comment out the next line if QDERM test is to run
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
    idx1 = luigi.IntParameter()
    idx2 = luigi.IntParameter()

    ds_type = DatasetType.ARG25

    def requires(self):
        query_result = pjoin(self.out_path, CONFIG.get('agdc', 'query_result'))
        with open(query_result, 'r') as f:
            query = pickle.load(f)

        tasks = []
        for ds in query[self.idx1:self.idx2]:
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
    # Setup command-line arguments
    desc = "Processes fractional cover.."
    hlp = ("The tile/chunk index to retieve from the tiles list. "
           "(Needs to have been previously computed to a file named chunks.pkl")
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('--tile', type=int, help=hlp)

    parsed_args = parser.parse_args()
    tile_idx = parsed_args.tile

    # setup logging
    log_dir = CONFIG.get('work', 'logs_directory')
    if not exists(log_dir):
        os.makedirs(log_dir)

    logfile = "{log_path}/run_wc_{uname}_{pid}.log"
    logfile = logfile.format(log_path=log_dir, uname=os.uname()[1],
                             pid=os.getpid())
    logging_level = logging.INFO
    logging.basicConfig(filename=logfile, level=logging_level,
                        format=("%(asctime)s: [%(name)s] (%(levelname)s) "
                                "%(message)s "))

    out_dir = CONFIG.get('agdc', 'output_directory')

    # Get the chunks list (each node will process a chunk determined by the
    # chunks.pkl, and the index to the chunk is determined at run time
    chunks_list_fname = pjoin(out_dir, CONFIG.get('agdc', 'chunks_list'))
    with open(chunks_list_fname, 'r') as f:
        chunks = pickle.load(f)

    # Initialise the job to luigi
    chunk = chunks[tile_idx]
    tasks = [ProcessFC(out_dir, chunk[0], chunk[1])]
    luigi.build(tasks, local_scheduler=True, workers=16)
    luigi.run()


    # *********************** QDERM TEST ******************************
    # For the QDERM test to run, a few lines need to be commented and
    # uncommented both here in agdc_fc.py and fc_utils.py

    #task = [FractionalCoverTask('/g/data1/v10/testing_ground/jps547/FC/QDERM_test_data/l5tmre_p095r084_20090411/l5tmre_p095r084_20090411_dbgm4.img',
    #                            '/g/data1/v10/testing_ground/jps547/FC/QDERM_test_data/l5tmre_p095r084_20090411')]
    #mpi.run(task)
