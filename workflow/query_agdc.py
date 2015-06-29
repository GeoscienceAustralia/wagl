#!/usr/bin/env python

import luigi
from os.path import join as pjoin, dirname
import cPickle as pickle
import subprocess
from datetime import date
import argparse

from datacube.api.model import DatasetType, Satellite
from datacube.api.query import list_tiles_as_list
from datacube.config import Config


PBS_DSH = (
"""#!/bin/bash

#PBS -P {project}
#PBS -q {queue}
#PBS -l walltime={walltime},ncpus={ncpus},mem={mem}GB,jobfs=350GB
#PBS -l wd
#PBS -me
#PBS -M {email}

NNODES={nnodes}

for i in $(seq 1 $NNODES); do
   pbsdsh -n $((16 *$i)) -- bash -l -c "{modules} python {pyfile} \
   --tile $[$i - 1]{config}" &
done;
wait
""")


def create_tiles(array_size, tile_size=25):
    """
    A minor function to tile a 1D array or list into smaller manageable
    portions.
    """
    idx_start = range(0, array_size, tile_size)
    idx_tiles = []
    for idx_st in idx_start:
        if idx_st + tile_size < array_size:
            idx_end = idx_st + tile_size
        else:
            idx_end = array_size
        idx_tiles.append((idx_st, idx_end))

    return idx_tiles


def query(output_path, cfg=None):
    """
    Queries the database based on the parameters given by the
    config file (default is fc.cfg).
    """
    out_fname = pjoin(output_path, CONFIG.get('agdc', 'query_result'))

    config = Config()
    ds_type = DatasetType.ARG25

    # Get the satellites we wish to query
    satellites = CONFIG.get('agdc', 'satellites')
    satellites = [Satellite(i) for i in satellites.split(',')]

    # Get the min/max date range to query
    min_date = CONFIG.get('agdc', 'min_date')
    min_date = [int(i) for i in min_date.split('_')]
    min_date = date(min_date[0], min_date[1], min_date[2])
    max_date = CONFIG.get('agdc', 'max_date')
    max_date = [int(i) for i in max_date.split('_')]
    max_date = date(max_date[0], max_date[1], max_date[2])

    # Get the min/max cell range to query
    cell_xmin = int(CONFIG.get('agdc', 'cell_xmin'))
    cell_ymin = int(CONFIG.get('agdc', 'cell_ymin'))
    cell_xmax = int(CONFIG.get('agdc', 'cell_xmax'))
    cell_ymax = int(CONFIG.get('agdc', 'cell_ymax'))
    xcells = range(cell_xmin, cell_xmax, 1)
    ycells = range(cell_ymin, cell_ymax, 1)

    # Query the DB
    tiles = list_tiles_as_list(x=xcells, y=ycells, acq_min=min_date,
                               acq_max=max_date,
                               satellites=satellites,
                               datasets=ds_type,
                               database=config.get_db_database(),
                               user=config.get_db_username(),
                               password=config.get_db_password(),
                               host=config.get_db_host(),
                               port=config.get_db_port()) 

    # Output the result to disk
    with open(out_fname, 'w') as outf:
        pickle.dump(tiles, outf)

    # Setup the modules to use for the job
    modules = CONFIG.get('pbs', 'modules').split(',')

    modules = ['module load {}; '.format(module) for module in modules]
    modules.insert(0, 'module use /projects/u46/opt/modules/modulefiles;')
    modules = ''.join(modules)

    # Tile/chunk up the number of items we have to process
    chunk_size = int(CONFIG.get('agdc', 'chunks_per_node'))
    chunks = create_tiles(len(tiles), tile_size=chunk_size)
    chunks_out_fname = pjoin(output_path, CONFIG.get('agdc', 'chunks_list'))
    with open(chunks_out_fname, 'w') as outf:
        pickle.dump(chunks, outf)

    # Calculate the job node and cpu requirements
    nnodes = len(chunks)
    ncpus = nnodes * 16
    mem = nnodes * 32

    # Populate the PBS shell script
    project = CONFIG.get('pbs', 'project')
    queue = CONFIG.get('pbs', 'queue')
    walltime = CONFIG.get('pbs', 'walltime')
    email = CONFIG.get('pbs', 'email')
    py_file = pjoin(dirname(__file__), 'agdc_fc.py')
    if cfg is None:
        cfg = ''
    else:
        cfg = ' --cfg {}'.format(cfg)
    pbs_job = PBS_DSH.format(project=project, queue=queue,
                             walltime=walltime, ncpus=ncpus, mem=mem,
                             email=email, nnodes=nnodes,
                             modules=modules, pyfile=py_file, config=cfg)

    # Output the shell script to disk
    pbs_fname = pjoin(output_path, CONFIG.get('pbs', 'pbs_filename'))
    with open(pbs_fname, 'w') as out_file:
        out_file.write(pbs_job)

    # Execute the job
    subprocess.check_call(['qsub', pbs_fname])

if __name__ == '__main__':
    desc = ("The database query mechanism used for processing "
            "NBAR to Fractional Cover.")
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('--cfg',
                       help='Path to a user defined configuration file.')

    parsed_args = parser.parse_args()
    cfg = parsed_args.cfg

    # Setup the config file
    global CONFIG
    if cfg is None:
        CONFIG = luigi.configuration.get_config()
        CONFIG.add_config_path(pjoin(dirname(__file__), 'fc.cfg'))
    else:
        CONFIG = luigi.configuration.get_config()
        CONFIG.add_config_path(cfg)

    out_dir = CONFIG.get('agdc', 'output_directory')
    query(out_dir, cfg)
