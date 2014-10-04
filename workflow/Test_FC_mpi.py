#!/bin/env python

#import luigi
import luigi.contrib.mpi as mpi
from fc import FractionalCoverGenerator
import os

if __name__ == '__main__':
    outputRoot = os.getenv('OUTPUT_ROOT', '/tmp/smr')

    tasks = []
    with open('NBAR_test_scenes_64.txt', 'r') as infile:
        for line in infile:
            nbarPath = line.strip()
            fcName = os.path.basename(nbarPath).replace('NBAR', 'FC')
            tasks.append(FractionalCoverTask(nbarPath, '%s/output/%s' % (outputRoot, fcName)))

    mpi.run(tasks)
