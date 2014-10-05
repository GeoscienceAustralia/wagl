#!/bin/env python

import luigi
import os
from ULA3.fc import fractional_cover
from EOtools.DatasetDrivers import SceneDataset
import logging
from memuseFilter import MemuseFilter


class FractionalCoverTask(luigi.Task):
    nbar_path = luigi.Parameter()
    fc_path = luigi.Parameter()

    def output(self):
        return FCDataset(self.fc_path)

    def requires(self):
        return NBARTask(self.nbar_path)

    def run(self):
        logging.info("In FractionalCoverTask.run method")
        result = fractional_cover(self.input().nbar_path,  asfloat32=False,
            fc_data_path=self.output().path, single_tif=False)


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
        self.dataset = SceneDataset(pathname=self.nbar_path)

    def exists(self):
        return os.path.exists(self.nbar_path)


if __name__ == '__main__':
    logging.config.fileConfig('logging.conf') # Get basic config
    log = logging.getLogger('')               # Get root logger
    f = MemuseFilter()                        # Create filter
    log.handlers[0].addFilter(f)         # The ugly part:adding filter to handler
    log.info("FC started")
    luigi.run()
