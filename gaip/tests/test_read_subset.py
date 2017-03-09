#!/usr/bin/env python

from __future__ import absolute_import
import os
import shutil
import tempfile
import unittest

import numpy

from gaip import read_subset
from gaip import write_img
from gaip.tests import unittesting_tools as ut


class TestReadSubset(unittest.TestCase):

    def testWestBounds(self):
        """
        Test that a co-ordinate west of the image domain returns an
        index error.
        The subset attempts to read a 20 by 20 block with half contained
        within the image bounds and half contained outside the image
        bounds.
        """
        img, geobox = ut.createTestImage()
    
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
                 'ULxy': UL,
                 'URxy': UR,
                 'LRxy': LR,
                 'LLxy': LL}

        self.assertRaises(IndexError, read_subset, **kwds)
    
        # Cleanup
        shutil.rmtree(temp_dir)

    
    def testEastBounds(self):
        """
        Test that a co-ordinate east of the image domain returns an
        index error.
        The subset attempts to read a 20 by 20 block with half contained
        within the image bounds and half contained outside the image
        """
        img, geobox = ut.createTestImage()
    
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
                 'ULxy': UL,
                 'URxy': UR,
                 'LRxy': LR,
                 'LLxy': LL}

        self.assertRaises(IndexError, read_subset, **kwds)
    
        # Cleanup
        shutil.rmtree(temp_dir)


    def testNorthBounds(self):
        """
        Test that a co-ordinate north of the image domain returns an
        index error.
        The subset attempts to read a 20 by 20 block with half contained
        within the image bounds and half contained outside the image
        """
        img, geobox = ut.createTestImage()
    
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
                 'ULxy': UL,
                 'URxy': UR,
                 'LRxy': LR,
                 'LLxy': LL}

        self.assertRaises(IndexError, read_subset, **kwds)
    
        # Cleanup
        shutil.rmtree(temp_dir)


    def testSouthBounds(self):
        """
        Test that a co-ordinate south of the image domain returns an
        index error.
        The subset attempts to read a 20 by 20 block with half contained
        within the image bounds and half contained outside the image
        """
        img, geobox = ut.createTestImage()
    
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
                 'ULxy': UL,
                 'URxy': UR,
                 'LRxy': LR,
                 'LLxy': LL}

        self.assertRaises(IndexError, read_subset, **kwds)
    
        # Cleanup
        shutil.rmtree(temp_dir)


    def test_correct_subset(self):
        """
        Test that the subset is what we expect.
        Read a 10 by 10 starting at the UL corner.
        """
        img, geobox = ut.createTestImage()
    
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
                 'ULxy': UL,
                 'URxy': UR,
                 'LRxy': LR,
                 'LLxy': LL}

        subs, geobox = read_subset(**kwds)

        base = img[0:10,0:10]

        result = numpy.sum(base - subs)

        self.assertTrue(result == 0)
    
        # Cleanup
        shutil.rmtree(temp_dir)


if __name__ == '__main__':
    unittest.main()
