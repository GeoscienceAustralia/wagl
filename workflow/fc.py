#!/bin/env python

import luigi
import os
from EOtools.DatasetDrivers import SceneDataset


class FractionalCoverGenerator(luigi.Task):
    nbar_path = luigi.Parameter()
    output_base_dir = luigi.Parameter()

    def output(self):
        output_path = "%s/dummy_fc_output.txt" % (self.output_base_dir,)     
        return FractionalCoverDataset(output_path)

    def requires(self):
        return NBAR_dataset_generator(self.nbar_path)

    def run(self):
        print "Hello world - processing %s" % (self.nbar_path,)
        print "Input is %s" % (self.input())
        print "Metadata is %s" % (self.input().dataset.GetMetadata())

        # some code here to write data such that self.output().exists() returns True


class FractionalCoverDataset(luigi.Target):

    def __init__(self, output_path):
        self.output_path = output_path

    def exists(self):
        return os.path.exists(self.output_path)

class NBAR_dataset_generator(luigi.ExternalTask):
    nbar_path = luigi.Parameter()

    def output(self):
        return NBAR_dataset(self.nbar_path)


class NBAR_dataset(luigi.Target):

    def __init__(self, nbar_path):
        self.nbar_path = nbar_path
        self.dataset = SceneDataset(pathname=self.nbar_path)

    def exists(self):
        return os.path.exists(self.nbar_path)


if __name__ == '__main__':
    luigi.run()
   




