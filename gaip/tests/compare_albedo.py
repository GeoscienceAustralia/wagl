#!/usr/bin/env python

from os.path import join as pjoin, exists
import numpy


with open('scenes.txt', 'r') as src:
    scenes = src.readlines()

with open('albedo-files.txt', 'r') as src:
    files = src.readlines()

scenes = [s.strip() for s in scenes]
files = [f.strip() for f in files]

ref_dir = '/g/data/v10/testing_ground/jps547/test-mod-base/output'
test_dir = '/g/data/v10/testing_ground/jps547/test-mod-input2/output'

for scene in scenes:
    ref_scene = pjoin(ref_dir, scene)
    test_scene = pjoin(test_dir, scene)
    for f in files:
        ref_fname = pjoin(pjoin(ref_scene, 'mod'), f)
        test_fname = pjoin(pjoin(test_scene, 'mod'), f)
        if not exists(ref_fname):
            continue
        print "Testing\nScene: {}\n File: {}".format(scene, f)
        with open(ref_fname, 'r') as ref, open(test_fname) as test:
            ref_data = ref.readlines()
            ref_data = [rd.strip() for rd in ref_data]
            test_data = test.readlines()
            test_data = [td.strip() for td in test_data]
            for i in range(len(ref_data)):
                try:
                    float(ref_data[i])
                    if not numpy.isclose(float(ref_data[i]),
                                         float(test_data[i])):
                        msg = "Line {} is not equivilent.".format(i)
                        print msg
                except ValueError:
                    if ref_data[i] != test_data[i]:
                        msg = "Line {} is not equivilent.".format(i)
                        print msg
