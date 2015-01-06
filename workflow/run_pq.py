#!/bin/env python
# 
# Runs Pixel Quality workflow against specified input data
#
#

import argparse
import luigi.contrib.mpi as mpi
import logging
import os
import sys
from pq import PixelQualityTask, nbar_name_from_l1t, pqa_name_from_l1t

def is_valid_directory(parser, arg):
    """Used by argparse"""
    if not os.path.exists(arg):
        parser.error("%s does not exist" % (arg, ))
    else:
        return arg

def nbar_name_from(l1t_fname):
    """
    Return an NBAR file name given a L1T file name
    """

if __name__ == '__main__':

    # command line arguments

    parser = argparse.ArgumentParser()
    parser.add_argument("--l1t_path", help="path to directory containing OTH datasets", \
        required=True, type=lambda x: is_valid_directory(parser, x))
    parser.add_argument("--nbar_path", help="path to directory containing NBAR datasets", \
        required=True, type=lambda x: is_valid_directory(parser, x))
    parser.add_argument("--land_sea_path", help="path to directory containing Land/Sea datasets", \
        default='/g/data/v10/eoancillarydata/Land_Sea_Rasters', \
        type=lambda x: is_valid_directory(parser, x))
    parser.add_argument("--out_path", help="path to directory where PQA dataset is to be written", \
        required=True, type=lambda x: is_valid_directory(parser, x))
    parser.add_argument("--log_path", help="path to directory where where log files will be written", \
        default='.', type=lambda x: is_valid_directory(parser, x))
    parser.add_argument("--debug", help="selects more detail logging (default is INFO)", \
        default=False, action='store_true')

    args = parser.parse_args()

    # setup logging
    
    logfile = "%s/run_pq_%s_%d.log" % (args.log_path, os.uname()[1], os.getpid())
    logging_level = logging.INFO
    if args.debug:
        logging_level = logging.DEBUG
    logging.basicConfig(filename=logfile, level=logging_level, \
        format= '[%(asctime)s] {%(pathname)s:%(lineno)d} %(levelname)s - %(message)s',
        datefmt='%H:%M:%S')
    logging.info("run_pq started")


    logging.info('l1t_path=%s' % (args.l1t_path, ))
    logging.info('nbar_path=%s' % (args.nbar_path, ))
    logging.info('land_sea_path=%s' % (args.land_sea_path, ))
    logging.info('out_path=%s' % (args.out_path, ))
    logging.info('log_path=%s' % (args.log_path, ))

    # create the task list based on L1T files to processa

    tasks = []
    for l1t_file in [f for f in os.listdir(args.l1t_path) if '_OTH_' in f]:
        l1t_dataset_path = os.path.join(args.l1t_path, l1t_file)
        nbar_dataset_path = os.path.join(args.nbar_path, nbar_name_from_l1t(l1t_file))
        pqa_dataset_path = os.path.join(args.out_path, pqa_name_from_l1t(l1t_file))

        tasks.append( \
            PixelQualityTask( \
                l1t_dataset_path,
                nbar_dataset_path,
                args.land_sea_path,
                pqa_dataset_path \
            ) \
        )
        
        print nbar_dataset_path

    mpi.run(tasks)
