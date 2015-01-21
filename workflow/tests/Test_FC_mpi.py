#!/bin/env python

#import luigi
import luigi.contrib.mpi as mpi
from fc import FractionalCoverTask
import os
import logging
from memuseFilter import MemuseFilter


if __name__ == '__main__':
    logging.config.fileConfig('logging.conf') # Get basic config
    log = logging.getLogger('')               # Get root logger
    f = MemuseFilter()                        # Create filter
    log.handlers[0].addFilter(f)         # The ugly part:adding filter to handler
    log.info("FC started")

    outputRoot = os.getenv('OUTPUT_ROOT', '/tmp/smr')

    tasks = []
    with open('NBAR_test_scenes_64.txt', 'r') as infile:
        for line in infile:
            nbarPath = line.strip()
            fcName = os.path.basename(nbarPath).replace('NBAR', 'FC')
            tasks.append(FractionalCoverTask(nbarPath, '%s/output/%s' % (outputRoot, fcName)))

    mpi.run(tasks)
