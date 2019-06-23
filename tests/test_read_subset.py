#!/usr/bin/env python

from __future__ import absolute_import
import os
import shutil
import tempfile
import unittest

import numpy
import h5py

from wagl.data import read_subset
from wagl.data import write_img
from wagl import unittesting_tools as ut


class TestReadSubset(unittest.TestCase):

    img, geobox = ut.create_test_image((1200, 1500))
    img[:] = 1

    fid = h5py.File('test-subset.h5', backing_store=False, driver='core')
    ds = fid.create_dataset('data', data=img)
    ds.attrs['geotransform'] = geobox.transform.to_gdal()
    ds.attrs['crs_wkt'] = geobox.crs.ExportToWkt()
    ds.attrs['fillvalue'] = 0

    subs_shape = (200, 300)

    @unittest.skip("Refactor DSM subsetting logic; TODO update test")
    def testWestBounds(self):
        """
        Test that a co-ordinate west of the image domain returns an
        index error.
        The subset attempts to read a 20 by 20 block with half contained
        within the image bounds and half contained outside the image
        bounds.
        """
        img, geobox = ut.create_test_image()
    
        # Temporarily write the image to disk
        temp_dir = tempfile.mkdtemp()
        fname = os.path.join(temp_dir, 'testWestBounds')
        write_img(img, fname, geobox=geobox)
    
        # Create box to read 10 pixels left of the image bounds
        UL = geobox.convert_coordinates((-9, 0))
        UR = geobox.convert_coordinates((9, 0))
        LR = geobox.convert_coordinates((9, 10))
        LL = geobox.convert_coordinates((-9, 10))
    
        kwds = {'fname': fname,
                 'ul_xy': UL,
                 'ur_xy': UR,
                 'lr_xy': LR,
                 'll_xy': LL}

        self.assertRaises(IndexError, read_subset, **kwds)
    
        # Cleanup
        shutil.rmtree(temp_dir)

    
    @unittest.skip("Refactor DSM subsetting logic; TODO update test")
    def testEastBounds(self):
        """
        Test that a co-ordinate east of the image domain returns an
        index error.
        The subset attempts to read a 20 by 20 block with half contained
        within the image bounds and half contained outside the image
        """
        img, geobox = ut.create_test_image()
    
        cols, rows = geobox.get_shape_xy()
    
        # Temporarily write the image to disk
        temp_dir = tempfile.mkdtemp()
        fname = os.path.join(temp_dir, 'testEastBounds')
        write_img(img, fname, geobox=geobox)
    
        # Create box to read 10 pixels right of the image bounds
        UL = geobox.convert_coordinates((cols-9, 0))
        UR = geobox.convert_coordinates((cols+10, 0))
        LR = geobox.convert_coordinates((cols+10, 10))
        LL = geobox.convert_coordinates((cols-9, 10))
    
        kwds = {'fname': fname,
                 'ul_xy': UL,
                 'ur_xy': UR,
                 'lr_xy': LR,
                 'll_xy': LL}

        self.assertRaises(IndexError, read_subset, **kwds)
    
        # Cleanup
        shutil.rmtree(temp_dir)


    @unittest.skip("Refactor DSM subsetting logic; TODO update test")
    def testNorthBounds(self):
        """
        Test that a co-ordinate north of the image domain returns an
        index error.
        The subset attempts to read a 20 by 20 block with half contained
        within the image bounds and half contained outside the image
        """
        img, geobox = ut.create_test_image()
    
        # Temporarily write the image to disk
        temp_dir = tempfile.mkdtemp()
        fname = os.path.join(temp_dir, 'testNorthBounds')
        write_img(img, fname, geobox=geobox)
    
        # Create box to read 10 pixels above the image bounds
        UL = geobox.convert_coordinates((0, -9))
        UR = geobox.convert_coordinates((10, -9))
        LR = geobox.convert_coordinates((10, 10))
        LL = geobox.convert_coordinates((0, 10))
    
        kwds = {'fname': fname,
                 'ul_xy': UL,
                 'ur_xy': UR,
                 'lr_xy': LR,
                 'll_xy': LL}

        self.assertRaises(IndexError, read_subset, **kwds)
    
        # Cleanup
        shutil.rmtree(temp_dir)


    @unittest.skip("Refactor DSM subsetting logic; TODO update test")
    def testSouthBounds(self):
        """
        Test that a co-ordinate south of the image domain returns an
        index error.
        The subset attempts to read a 20 by 20 block with half contained
        within the image bounds and half contained outside the image
        """
        img, geobox = ut.create_test_image()
    
        cols, rows = geobox.get_shape_xy()
    
        # Temporarily write the image to disk
        temp_dir = tempfile.mkdtemp()
        fname = os.path.join(temp_dir, 'testSouthBounds')
        write_img(img, fname, geobox=geobox)
    
        # Create box to read 10 pixels below the image bounds
        UL = geobox.convert_coordinates((0, rows-9))
        UR = geobox.convert_coordinates((10, rows-9))
        LR = geobox.convert_coordinates((10, rows+10))
        LL = geobox.convert_coordinates((0, rows+10))
    
        kwds = {'fname': fname,
                 'ul_xy': UL,
                 'ur_xy': UR,
                 'lr_xy': LR,
                 'll_xy': LL}

        self.assertRaises(IndexError, read_subset, **kwds)
    
        # Cleanup
        shutil.rmtree(temp_dir)


    @unittest.skip('Requires refactoring')
    def test_correct_subset(self):
        """
        Test that the subset is what we expect.
        Read a 10 by 10 starting at the UL corner.
        """
        img, geobox = ut.create_test_image()
    
        cols, rows = geobox.get_shape_xy()
    
        # Temporarily write the image to disk
        temp_dir = tempfile.mkdtemp()
        fname = os.path.join(temp_dir, 'test_image')
        write_img(img, fname, geobox=geobox)

        # Create box to read 10 pixels below the image bounds
        UL = geobox.convert_coordinates((0, 0))
        UR = geobox.convert_coordinates((9, 0))
        LR = geobox.convert_coordinates((9, 9))
        LL = geobox.convert_coordinates((0, 9))
    
        kwds = {'fname': fname,
                 'ul_xy': UL,
                 'ur_xy': UR,
                 'lr_xy': LR,
                 'll_xy': LL}

        subs, geobox = read_subset(**kwds)

        base = img[0:10,0:10]

        result = numpy.sum(base - subs)

        self.assertTrue(result == 0)
    
        # Cleanup
        shutil.rmtree(temp_dir)

    def test_case_a(self):
        """
        Origin = (-150, -50); `O`

           O+----+
            -    -
            -  +-----------+
            -  - -         -
            +----+         -
               -           -
               -           -
               -           -
               -           -
               +-----------+
        """
        # indices based on full array
        ul = (-150, -50)
        ur = (ul[0], ul[1] + self.subs_shape[1])
        lr = (ul[0] + self.subs_shape[0], ul[1] + self.subs_shape[1])
        ll = (ul[0] + self.subs_shape[0], ul[1])

        # real world coords (note reversing (y, x) to (x, y)
        ul_xy_map = ul[::-1] * self.geobox.transform
        ur_xy_map = ur[::-1] * self.geobox.transform
        lr_xy_map = lr[::-1] * self.geobox.transform
        ll_xy_map = ll[::-1] * self.geobox.transform

        # read subset
        data, gb = read_subset(self.ds, ul_xy_map, ur_xy_map, lr_xy_map, ll_xy_map)

        count = 50 * 250
        self.assertTrue(data.sum() == count)

    def test_case_b(self):
        """
        Origin = (-150, 600); `O`

           O+---+
            -   -
        +----   ----+
        -   -   -   -
        -   +---+   -
        -           -
        -           -
        -           -
        -           -
        +-----------+
        """
        # indices based on full array
        ul = (-150, 600)
        ur = (ul[0], ul[1] + self.subs_shape[1])
        lr = (ul[0] + self.subs_shape[0], ul[1] + self.subs_shape[1])
        ll = (ul[0] + self.subs_shape[0], ul[1])

        # real world coords (note reversing (y, x) to (x, y)
        ul_xy_map = ul[::-1] * self.geobox.transform
        ur_xy_map = ur[::-1] * self.geobox.transform
        lr_xy_map = lr[::-1] * self.geobox.transform
        ll_xy_map = ll[::-1] * self.geobox.transform

        # read subset
        data, gb = read_subset(self.ds, ul_xy_map, ur_xy_map, lr_xy_map, ll_xy_map)

        count = 50 * 300
        self.assertTrue(data.sum() == count)

    def test_case_c(self):
        """
        Origin = (-150, 1400); `O`

                  O+---+
                   -   -
        +------------+ -
        -          - - -
        -          +---+
        -            -
        -            -
        -            -
        -            -
        +------------+
        """
        # indices based on full array
        ul = (-150, 1400)
        ur = (ul[0], ul[1] + self.subs_shape[1])
        lr = (ul[0] + self.subs_shape[0], ul[1] + self.subs_shape[1])
        ll = (ul[0] + self.subs_shape[0], ul[1])

        # real world coords (note reversing (y, x) to (x, y)
        ul_xy_map = ul[::-1] * self.geobox.transform
        ur_xy_map = ur[::-1] * self.geobox.transform
        lr_xy_map = lr[::-1] * self.geobox.transform
        ll_xy_map = ll[::-1] * self.geobox.transform

        # read subset
        data, gb = read_subset(self.ds, ul_xy_map, ur_xy_map, lr_xy_map, ll_xy_map)

        count = 50 * 100
        self.assertTrue(data.sum() == count)

    def test_case_d(self):
        """
        Origin = (600, -50); `O`

              +-----------+
              -           -
          O+-----+        -
           -  -  -        -
           -  -  -        -
           +-----+        -
              -           -
              +-----------+
        """
        # indices based on full array
        ul = (600, -50)
        ur = (ul[0], ul[1] + self.subs_shape[1])
        lr = (ul[0] + self.subs_shape[0], ul[1] + self.subs_shape[1])
        ll = (ul[0] + self.subs_shape[0], ul[1])

        # real world coords (note reversing (y, x) to (x, y)
        ul_xy_map = ul[::-1] * self.geobox.transform
        ur_xy_map = ur[::-1] * self.geobox.transform
        lr_xy_map = lr[::-1] * self.geobox.transform
        ll_xy_map = ll[::-1] * self.geobox.transform

        # read subset
        data, gb = read_subset(self.ds, ul_xy_map, ur_xy_map, lr_xy_map, ll_xy_map)

        count = 200 * 250
        self.assertTrue(data.sum() == count)

    def test_case_e(self):
        """
        Origin = (600, 600); `O`

        +-----------+
        -           -
        -  O+---+   -
        -   -   -   -
        -   -   -   -
        -   +---+   -
        -           -
        +-----------+
        """
        # indices based on full array
        ul = (600, 600)
        ur = (ul[0], ul[1] + self.subs_shape[1])
        lr = (ul[0] + self.subs_shape[0], ul[1] + self.subs_shape[1])
        ll = (ul[0] + self.subs_shape[0], ul[1])

        # real world coords (note reversing (y, x) to (x, y)
        ul_xy_map = ul[::-1] * self.geobox.transform
        ur_xy_map = ur[::-1] * self.geobox.transform
        lr_xy_map = lr[::-1] * self.geobox.transform
        ll_xy_map = ll[::-1] * self.geobox.transform

        # read subset
        data, gb = read_subset(self.ds, ul_xy_map, ur_xy_map, lr_xy_map, ll_xy_map)

        count = 200 * 300
        self.assertTrue(data.sum() == count)

    def test_case_f(self):
        """
        Origin = (600, 1400); `O`

        +-----------+
        -           -
        -       O+-----+
        -        -  -  -
        -        -  -  -
        -        +-----+
        -           -
        +-----------+
        """
        # indices based on full array
        ul = (600, 1400)
        ur = (ul[0], ul[1] + self.subs_shape[1])
        lr = (ul[0] + self.subs_shape[0], ul[1] + self.subs_shape[1])
        ll = (ul[0] + self.subs_shape[0], ul[1])

        # real world coords (note reversing (y, x) to (x, y)
        ul_xy_map = ul[::-1] * self.geobox.transform
        ur_xy_map = ur[::-1] * self.geobox.transform
        lr_xy_map = lr[::-1] * self.geobox.transform
        ll_xy_map = ll[::-1] * self.geobox.transform

        # read subset
        data, gb = read_subset(self.ds, ul_xy_map, ur_xy_map, lr_xy_map, ll_xy_map)

        count = 200 * 100
        self.assertTrue(data.sum() == count)

    def test_case_g(self):
        """
        Origin = (1100, -50); `O`

              +-----------+
              -           -
              -           -
              -           -
              -           -
          O+-----+        -
           -  -  -        -
           -  +-----------+
           -     -
           -     -
           +-----+
        """
        # indices based on full array
        ul = (1100, -50)
        ur = (ul[0], ul[1] + self.subs_shape[1])
        lr = (ul[0] + self.subs_shape[0], ul[1] + self.subs_shape[1])
        ll = (ul[0] + self.subs_shape[0], ul[1])

        # real world coords (note reversing (y, x) to (x, y)
        ul_xy_map = ul[::-1] * self.geobox.transform
        ur_xy_map = ur[::-1] * self.geobox.transform
        lr_xy_map = lr[::-1] * self.geobox.transform
        ll_xy_map = ll[::-1] * self.geobox.transform

        # read subset
        data, gb = read_subset(self.ds, ul_xy_map, ur_xy_map, lr_xy_map, ll_xy_map)

        count = 100 * 250
        self.assertTrue(data.sum() == count)

    def test_case_h(self):
        """
        Origin = (1100, 600); `O`

        +-----------+
        -           -
        -           -
        -           -
        -           -
        -  O+----+  -
        -   -    -  -
        +-----------+
            -    -
            -    -
            +----+
        """
        # indices based on full array
        ul = (1100, 600)
        ur = (ul[0], ul[1] + self.subs_shape[1])
        lr = (ul[0] + self.subs_shape[0], ul[1] + self.subs_shape[1])
        ll = (ul[0] + self.subs_shape[0], ul[1])

        # real world coords (note reversing (y, x) to (x, y)
        ul_xy_map = ul[::-1] * self.geobox.transform
        ur_xy_map = ur[::-1] * self.geobox.transform
        lr_xy_map = lr[::-1] * self.geobox.transform
        ll_xy_map = ll[::-1] * self.geobox.transform

        # read subset
        data, gb = read_subset(self.ds, ul_xy_map, ur_xy_map, lr_xy_map, ll_xy_map)

        count = 100 * 300
        self.assertTrue(data.sum() == count)

    def test_case_i(self):
        """
        Origin = (1100, 1400); `O`

        +-----------+
        -           -
        -           -
        -           -
        -           -
        -       O+-----+
        -        -     -
        +-----------+---
                 -     -
                 -     -
                 +-----+
        """
        # indices based on full array
        ul = (1100, 1400)
        ur = (ul[0], ul[1] + self.subs_shape[1])
        lr = (ul[0] + self.subs_shape[0], ul[1] + self.subs_shape[1])
        ll = (ul[0] + self.subs_shape[0], ul[1])

        # real world coords (note reversing (y, x) to (x, y)
        ul_xy_map = ul[::-1] * self.geobox.transform
        ur_xy_map = ur[::-1] * self.geobox.transform
        lr_xy_map = lr[::-1] * self.geobox.transform
        ll_xy_map = ll[::-1] * self.geobox.transform

        # read subset
        data, gb = read_subset(self.ds, ul_xy_map, ur_xy_map, lr_xy_map, ll_xy_map)

        count = 100 * 100
        self.assertTrue(data.sum() == count)

if __name__ == '__main__':
    unittest.main()
