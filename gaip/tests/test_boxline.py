#!/usr/bin/env python

from os.path import join as pjoin
import numpy as np
import pandas
import rasterio


def read_centreline(fname):
    df = pandas.read_csv(fname, skiprows=2, header=None, sep=r'\s+\s+',
                         engine='python',
                         names=['line', 'centre', 'npoints', 'lat', 'lon'])
    return df


def read_coordinator(fname):
    df = pandas.read_csv(fname, header=None, sep=r'\s+\s+', engine='python',
                         names=['row', 'col'])
    return df

def read_boxline(fname):
    df = pandas.read_csv(fname, header=None, sep=r'\s+\s+', engine='python',
                         names=['line', 'cstart', 'cend'])
    return df


# with open('nbar-dirs.txt', 'r') as src:
with open('scenes.txt', 'r') as src:
    dirs = src.readlines()
dirs = [d.strip() for d in dirs]

centreline_base = 'CENTRELINE'
coordinator_base = 'COORDINATOR'
boxline_base = 'BOXLINE'
ref_dir = '/g/data/v10/agdc/jez/galpgs/candidate7/out/ls5/nbar/2008/05/output'
test_dir = '/g/data/v10/testing_ground/jps547/test-boxline/output'

for d in dirs:
    print "Scene: {}".format(d)
    coordinator_ref_fname = pjoin(pjoin(ref_dir, d), coordinator_base)
    coordinator_test_fname = pjoin(pjoin(test_dir, d), coordinator_base)
    boxline_ref_fname = pjoin(pjoin(ref_dir, d), boxline_base)
    boxline_test_fname = pjoin(pjoin(test_dir, d), boxline_base)


    coordinator_ref = read_coordinator(coordinator_ref_fname)
    coordinator_test = read_coordinator(coordinator_test_fname)
    boxline_ref = read_boxline(boxline_ref_fname)
    boxline_test = read_boxline(boxline_test_fname)

    print "boxline check"
    print "line check: {}".format((boxline_ref['line'] -
                                   boxline_test['line']).sum())
    print "cstart check {}".format((boxline_ref['cstart'] -
                                    boxline_test['cstart']).sum())
    print "cend check: {}".format((boxline_ref['cend'] -
                                   boxline_test['cend']).sum())
    print "coordinator check"
    print "row check: {}".format((coordinator_ref['row'] -
                                  coordinator_test['row']).sum())
    print "col check: {}".format((coordinator_ref['col'] -
                                  coordinator_test['col']).sum())
