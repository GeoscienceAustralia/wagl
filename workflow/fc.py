#!/bin/env python

from __future__ import absolute_import
import luigi
import argparse
import os
from os.path import join as pjoin, dirname
import logging
import gaip
import numpy

CONFIG = luigi.configuration.get_config()
CONFIG.add_config_path(pjoin(dirname(__file__), 'fc.cfg'))

class FractionalCoverTask(luigi.Task):
    nbar_path = luigi.Parameter()
    fc_path = luigi.Parameter()

    def output(self):
        return FCDataset(self.fc_path)

    def requires(self):
        return NBARTask(self.nbar_path)

    def run(self):
        logging.info("In FractionalCoverTask.run method, NBAR={}, "
                     "output={}".format(self.input().nbar_path,
                                          self.output().path))

        # create the output directory
        logging.debug("creating output directory {}".format(self.fc_path))
        gaip.create_dir(self.fc_path)

        # Get the processing tile sizes
        x_tile = int(CONFIG.get('work', 'x_tile_size'))
        y_tile = int(CONFIG.get('work', 'y_tile_size'))
        x_tile = None if x_tile <= 0 else x_tile
        y_tile = None if y_tile <= 0 else y_tile

        # Get the fractional component short names and base output format
        fraction_names = CONFIG.get('scene', 'fractions').split(',')
        output_format = CONFIG.get('scene', 'output_format')

        # Get the wavelengths to filter by
        min_lambda = float(CONFIG.get('scene', 'min_lambda'))
        max_lambda = float(CONFIG.get('scene', 'max_lambda'))

        # Define the output targets
        out_fnames = []
        for component in fraction_names:
            out_fname = output_format.format(fraction=component)
            out_fnames.append(pjoin(self.fc_path, out_fname))

        # Get the acquisitions and filter by wavelength
        # Using the centre of wavelength range would be more ideal
        acqs = gaip.acquisitions(self.nbar_path)
        acqs = [acq for acq in acqs if (acq.band_type == gaip.REF and
                                        acq.wavelength[1] > min_lambda and 
                                        acq.wavelength[1] <= max_lambda)]

        # Run
        gaip.fractional_cover(acqs, x_tile, y_tile, out_fnames)

        logging.info("Done processing")

class FCDataset(luigi.Target):

    def __init__(self, path):
        self.path = path

    def exists(self):
        return os.path.exists(self.path)

class NBARTask(luigi.ExternalTask):
    nbar_path = luigi.Parameter()

    def output(self):
        return NBARdataset(self.nbar_path)

class NBARdataset(luigi.Target):

    def __init__(self, nbar_path):
        self.nbar_path = nbar_path
        self.acquisitions = gaip.acquisitions(self.nbar_path)

    def exists(self):
        return os.path.exists(self.nbar_path)

def is_valid_directory(parser, arg):
    """Used by argparse"""
    if not os.path.exists(arg):
        parser.error("{} does not exist".format(arg))
    else:
        return arg

def fc_name_from_nbar(nbar_fname):
    """
    Return an NBAR file name given a L1T file name
    """
    return nbar_fname.replace('NBAR', 'FC')

def scatter(iterable, P=1, p=1):
    """
    Scatter an iterator across `P` processors where `p` is the index
    of the current processor. This partitions the work evenly across
    processors.
    """
    import itertools
    return itertools.islice(iterable, p-1, None, P)
 
 
def main(nbar_path, outpath, nnodes=1, nodenum=1):
    nbar_files = sorted([pjoin(nbar_path, f) for f in os.listdir(nbar_path) if
                        '_NBAR_' in f])
    nbar_files = [basename(f) for f in scatter(l1t_files, nnodes, nodenum)]
    ncpus = int(os.getenv('PBS_NCPUS', '1'))
    tasks = []
    for nbar_file in nbar_files:
        nbar_dataset_path = pjoin(nbar_path, nbar_file)
        fc_dataset_path = os.path.join(outpath, fc_name_from_nbar(nbar_file))
        tasks.append(FractionalCoverTask(nbar_dataset_path,
                                         fc_dataset_path))

    luigi.build(tasks, local_scheduler=True, workers=16)


if __name__ == '__main__':
    # command line arguments

    parser = argparse.ArgumentParser()
    parser.add_argument("--nbar_path",
                        help="path to directory containing NBAR datasets",
                        required=True,
                        type=lambda x: is_valid_directory(parser, x))
    parser.add_argument("--out_path",
                        help=("path to directory where FC dataset is to "
                              "be written"), required=True,
                        type=lambda x: is_valid_directory(parser, x))
    parser.add_argument("--log_path",
                        help=("path to directory where where log files will "
                              "be written"), default='.',
                        type=lambda x: is_valid_directory(parser, x))
    parser.add_argument("--debug",
                        help="selects more detail logging (default is INFO)",
                        default=False, action='store_true')

    args = parser.parse_args()

    # setup logging
    fmt = '[%(asctime)s] {%(pathname)s:%(lineno)d} %(levelname)s - %(message)s'
    
    logfile = "run_fc_{}_{}.log".format(os.uname()[1], os.getpid())
    logfile = os.path.join(args.log_path, logfile)
    logging_level = logging.INFO
    if args.debug:
        logging_level = logging.DEBUG
    logging.basicConfig(filename=logfile, level=logging_level,
                        format=fmt, datefmt='%H:%M:%S')
    logging.info("fc.py started")


    logging.info('nbar_path={}'.format(args.nbar_path))
    logging.info('out_path={}'.format(args.out_path))
    logging.info('log_path={}'.format(args.log_path))


    size = int(os.getenv('PBS_NNODES', '1'))
    rank = int(os.getenv('PBS_VNODENUM', '1'))
    main(args.l1t_path, args.nbar_path, args.land_sea_path, args.out_path,
         size, rank)
