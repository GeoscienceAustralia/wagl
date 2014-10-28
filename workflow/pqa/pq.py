#!/bin/env python

import luigi
import os
from EOtools.DatasetDrivers import SceneDataset
import logging
from memuseFilter import MemuseFilter
from constants import PQAConstants
from pqa_result import PQAResult
from saturation_masking import setSaturationBits
from contiguity_masking import setContiguityBit


class PixelQualityTask(luigi.Task):
    l1t_path = luigi.Parameter()
    nbar_path = luigi.Parameter()
    pq_path = luigi.Parameter()

    def output(self):
        return PQDataset(self.pq_path)

    def requires(self):
        return NBARTask(self.nbar_path)

    def run(self):
        logging.info("In PixelQualityTask.run method, L1T=%s NBAR=%s, output=%s" %\
           (self.l1t_path, self.input().nbar_path, self.output().path))

        # read L1T data
        logging.debug("Creating L1T SceneDataset")
        l1t_sd = SceneDataset(self.l1t_path)
        logging.debug("Reading L1T bands")
        l1t_data = l1t_sd.ReadAsArray()
        logging.debug("l1t_data shape=%s" % (str(l1t_data.shape)))

        # read NBAR data
        logging.debug("Creating NBAR SceneDataset")
        nbar_sd = SceneDataset(self.nbar_path)
        logging.debug("Reading NBAR bands")
        nbar_data = nbar_sd.ReadAsArray()
        logging.debug("nbar_data shape=%s" % (str(nbar_data.shape)))

        # constants to be use for this PQA computation 

        sensor = l1t_sd.sensor
        logging.debug("setting constants for sensor=%s" % (sensor, ))
        pq_const = PQAConstants(sensor)

        # the PQAResult object for this run

        pqaResult = PQAResult(l1t_data[0].shape)

        # Saturation

        logging.debug("setting saturation bits")
        setSaturationBits(l1t_data, pq_const, pqaResult)
        logging.debug("done setting saturation bits")

        # contiguity

        logging.debug("setting contiguity bit")
        setContiguityBit(l1t_data, l1t_sd.satellite, pq_const, pqaResult)
        logging.debug("done setting contiguity bit")


        


class PQDataset(luigi.Target):

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
        self.dataset = SceneDataset(pathname=self.nbar_path)

    def exists(self):
        return os.path.exists(self.nbar_path)


if __name__ == '__main__':
    logging.config.fileConfig('logging.conf') # Get basic config
    log = logging.getLogger('')               # Get root logger
    f = MemuseFilter()                        # Create filter
    log.handlers[0].addFilter(f)         # The ugly part:adding filter to handler
    logging.info("PQA started")
    luigi.run()
    logging.info("PQA done")

