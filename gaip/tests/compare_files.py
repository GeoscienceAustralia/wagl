#!/usr/bin/env python

import rasterio
from os.path import join as pjoin, exists

with open('image-files.txt', 'r') as src:
    files = src.readlines()

with open('scenes.txt', 'r') as src:
    scenes = src.readlines()

files = [f.strip() for f in files]
scenes = [s.strip() for s in scenes]

s_dir = '/g/data2/v10/testing_ground/jps547/test-bilinear/single/output'
m_dir = '/g/data2/v10/testing_ground/jps547/test-bilinear/multi/output'
f_dir = '/g/data2/v10/testing_ground/jps547/test-bilinear/f2py/output'

ref_dir = '/g/data/v10/agdc/jez/galpgs/candidate7/out/ls5/nbar/2008/05/output'
test_dir = '/g/data/v10/testing_ground/jps547/test-boxline/output'

for scene in scenes:
    ref_scene = pjoin(ref_dir, scene)
    test_scene = pjoin(test_dir, scene)
    for f in files:
        ref_fname = pjoin(ref_scene, f)
        test_fname = pjoin(test_scene, f)
        if not exists(ref_fname):
            continue
        with rasterio.open(ref_fname) as ref_ds, rasterio.open(test_fname) as test_ds:
            print "Testing\nScene: {}\n File: {}".format(scene, f)
            ref_data = ref_ds.read(1)
            test_data = test_ds.read(1)
            diff = (ref_data - test_data).sum()
            if diff != 0:
                msg = "Mismatch:\nRef: {}\nTest: {}\nDifference: {}\n"
                print msg.format(ref_fname, test_fname, diff)
