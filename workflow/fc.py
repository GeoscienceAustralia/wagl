#!/bin/env python

<<<<<<< HEAD
import luigi.contrib.mpi as mpi
import luigi
import argparse
import os
import logging
import gaip
import numpy
import numexpr
from memuseFilter import MemuseFilter


class FractionalCoverTask(luigi.Task):
    nbar_path = luigi.Parameter()
    fc_path = luigi.Parameter()

    def output(self):
        return FCDataset(self.fc_path)

    def requires(self):
        return NBARTask(self.nbar_path)

    def run(self):
        logging.info("In FractionalCoverTask.run method, NBAR=%s, output=%s" % (self.input().nbar_path, self.output().path))

        # get acquisition stack 

        logging.debug("reading acquistions")
        acqs = gaip.acquisitions(self.nbar_path)
        logging.debug("stacking data")
        (acqs, stack, geo_box) = gaip.stack_data(acqs, \
            fn=(lambda acq: acq.band_type == gaip.REF and \
            acq.wavelength[1] > 0.52))

        logging.debug("read %d NBAR bands, coverage: %s" % (len(acqs), str(geo_box)))

        # set no data values to zero

        logging.debug("converting no_data to zeros")
        zero = numpy.int16(0)
        no_data = acqs[0].no_data 
        if no_data is None:
            no_data = -999
        numpy.maximum(stack, zero, out=stack)

        # call the fortran routiine

        logging.debug("calling the unmix() function")
        (green, dead1, dead2, bare, err) = gaip.unmix(stack)

        # change zero values back to no_data value

        logging.debug("unmix() done, setting no_data_values")
        wh_any = numpy.any(numexpr.evaluate("stack == zero"), axis=0)

        # scale factors

        sf2 = numpy.float32(0.01)
        sf3 = 10000

        # compute results

        logging.debug("computing the result bands")
        green = numexpr.evaluate("green * sf3").astype('int16')
        green[wh_any] = no_data
        dead = numexpr.evaluate("(dead1 + dead2) * sf3").astype('int16')
        dead[wh_any] = no_data
        bare = numexpr.evaluate("bare * sf3").astype('int16')
        bare[wh_any] = no_data
        unmix_err = numexpr.evaluate("err * sf2 * sf3").astype('int16')
        unmix_err[wh_any] = no_data

        # create the output directory

        logging.debug("creating output directory %s" % (self.fc_path, ))
        gaip.create_dir(self.fc_path)

        # output bands (one per GeoTiff file)

        ids = ['PV', 'NPV', 'BS', 'UE']
        new_stack = [green, dead, bare, unmix_err]

        for band_no in  range(len(new_stack)):
            logging.debug("writing band %s" % (ids[band_no], ))
            file_path = "%s/fc_%s.tif" % (self.fc_path, ids[band_no])
            gaip.write_img(new_stack[band_no], file_path, fmt="GTiff", \
                nodata=no_data)

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
        parser.error("%s does not exist" % (arg, ))
    else:
        return arg

def fc_name_from_nbar(nbar_fname):
    """
    Return an NBAR file name given a L1T file name
    """
    return nbar_fname.replace('NBAR', 'FC')

if __name__ == '__main__':
#     logging.info("FC started")
#     luigi.run()
#     import sys
#     sys.exit(0)

    # command line arguments

    parser = argparse.ArgumentParser()
    parser.add_argument("--nbar_path", help="path to directory containing NBAR datasets", \
        required=True, type=lambda x: is_valid_directory(parser, x))
    parser.add_argument("--out_path", help="path to directory where FC dataset is to be written", \
        required=True, type=lambda x: is_valid_directory(parser, x))
    parser.add_argument("--log_path", help="path to directory where where log files will be written", \
        default='.', type=lambda x: is_valid_directory(parser, x))
    parser.add_argument("--debug", help="selects more detail logging (default is INFO)", \
        default=False, action='store_true')

    args = parser.parse_args()

    # setup logging
    
    logfile = "%s/run_fc_%s_%d.log" % (args.log_path, os.uname()[1], os.getpid())
    logging_level = logging.INFO
    if args.debug:
        logging_level = logging.DEBUG
    logging.basicConfig(filename=logfile, level=logging_level, \
        format= '[%(asctime)s] {%(pathname)s:%(lineno)d} %(levelname)s - %(message)s',
        datefmt='%H:%M:%S')
    logging.info("fc.py started")


    logging.info('nbar_path=%s' % (args.nbar_path, ))
    logging.info('out_path=%s' % (args.out_path, ))
    logging.info('log_path=%s' % (args.log_path, ))

    # create the task list based on L1T files to process

    tasks = []
    for nbar_file in [f for f in os.listdir(args.nbar_path) if '_NBAR_' in f]:
        nbar_dataset_path = os.path.join(args.nbar_path, nbar_file)
        fc_dataset_path = os.path.join(args.out_path, fc_name_from_nbar(nbar_file))

        tasks.append( \
            FractionalCoverTask( \
                nbar_dataset_path,
                fc_dataset_path \
            ) \
        )
        
        print nbar_dataset_path

    mpi.run(tasks)


# if __name__ == '__main__':
#    logging.config.fileConfig('logging.conf') # Get basic config
#     log = logging.getLogger('')               # Get root logger
#     f = MemuseFilter()                        # Create filter
#     log.handlers[0].addFilter(f)         # The ugly part:adding filter to handler
#     log.info("FC started")
#     luigi.run()
