#!/usr/bin/env python

import os
import shutil
import tempfile
import unittest

import numpy

from gaip import load_2D_bin_file
from gaip import write_img
from gaip.tests import unittesting_tools as ut


class TestLoad2DBinFile(unittest.TestCase):

    def test_read_uint8(self):
        """
        Test that an UInt8 file is read correctly.
        """
        # Generate some test data
        img, geobox = ut.createTestImage()

        # Temporarily write the image to disk
        temp_dir = tempfile.mkdtemp()
        fname = os.path.join(temp_dir, 'testUInt8.dat')

        # Write to disk using an ENVI format (flat binary file)
        write_img(img, fname, geobox=geobox)

        # Read back into memory
        test = load_2D_bin_file(fname, img.shape[0], img.shape[1],
            img.dtype.name)

        # Compare
        total_diff = numpy.sum(img - test)

        self.assertEqual(total_diff, 0)

        # Cleanup
        shutil.rmtree(temp_dir)


    def test_read_int16(self):
        """
        Test that an Int16 file is read correctly.
        """
        # Generate some test data
        img, geobox = ut.createTestImage(dtype='Int16')

        # Temporarily write the image to disk
        temp_dir = tempfile.mkdtemp()
        fname = os.path.join(temp_dir, 'testInt16.dat')

        # Write to disk using an ENVI format (flat binary file)
        write_img(img, fname, geobox=geobox)

        # Read back into memory
        test = load_2D_bin_file(fname, img.shape[0], img.shape[1],
            img.dtype.name)

        # Compare
        total_diff = numpy.sum(img - test)

        self.assertEqual(total_diff, 0)

        # Cleanup
        shutil.rmtree(temp_dir)


    def test_read_uint16(self):
        """
        Test that an UInt16 file is read correctly.
        """
        # Generate some test data
        img, geobox = ut.createTestImage(dtype='UInt16')

        # Temporarily write the image to disk
        temp_dir = tempfile.mkdtemp()
        fname = os.path.join(temp_dir, 'testUInt16.dat')

        # Write to disk using an ENVI format (flat binary file)
        write_img(img, fname, geobox=geobox)

        # Read back into memory
        test = load_2D_bin_file(fname, img.shape[0], img.shape[1],
            img.dtype.name)

        # Compare
        total_diff = numpy.sum(img - test)

        self.assertEqual(total_diff, 0)

        # Cleanup
        shutil.rmtree(temp_dir)


    def test_read_int32(self):
        """
        Test that an Int32 file is read correctly.
        """
        # Generate some test data
        img, geobox = ut.createTestImage(dtype='Int32')

        # Temporarily write the image to disk
        temp_dir = tempfile.mkdtemp()
        fname = os.path.join(temp_dir, 'testInt32.dat')

        # Write to disk using an ENVI format (flat binary file)
        write_img(img, fname, geobox=geobox)

        # Read back into memory
        test = load_2D_bin_file(fname, img.shape[0], img.shape[1],
            img.dtype.name)

        # Compare
        total_diff = numpy.sum(img - test)

        self.assertEqual(total_diff, 0)

        # Cleanup
        shutil.rmtree(temp_dir)


    def test_read_uint32(self):
        """
        Test that an UInt32 file is read correctly.
        """
        # Generate some test data
        img, geobox = ut.createTestImage(dtype='UInt32')

        # Temporarily write the image to disk
        temp_dir = tempfile.mkdtemp()
        fname = os.path.join(temp_dir, 'testUInt32.dat')

        # Write to disk using an ENVI format (flat binary file)
        write_img(img, fname, geobox=geobox)

        # Read back into memory
        test = load_2D_bin_file(fname, img.shape[0], img.shape[1],
            img.dtype.name)

        # Compare
        total_diff = numpy.sum(img - test)

        self.assertEqual(total_diff, 0)

        # Cleanup
        shutil.rmtree(temp_dir)


    def test_read_float32(self):
        """
        Test that an Float32 file is read correctly.
        """
        # Generate some test data
        img, geobox = ut.createTestImage(dtype='Float32')

        # Temporarily write the image to disk
        temp_dir = tempfile.mkdtemp()
        fname = os.path.join(temp_dir, 'testFloat32.dat')

        # Write to disk using an ENVI format (flat binary file)
        write_img(img, fname, geobox=geobox)

        # Read back into memory
        test = load_2D_bin_file(fname, img.shape[0], img.shape[1],
            img.dtype.name)

        # Compare
        total_diff = numpy.sum(img - test)

        self.assertEqual(total_diff, 0)

        # Cleanup
        shutil.rmtree(temp_dir)


    def test_read_float64(self):
        """
        Test that an Float64 file is read correctly.
        """
        # Generate some test data
        img, geobox = ut.createTestImage(dtype='Float64')

        # Temporarily write the image to disk
        temp_dir = tempfile.mkdtemp()
        fname = os.path.join(temp_dir, 'testFloat64.dat')

        # Write to disk using an ENVI format (flat binary file)
        write_img(img, fname, geobox=geobox)

        # Read back into memory
        test = load_2D_bin_file(fname, img.shape[0], img.shape[1],
            img.dtype.name)

        # Compare
        total_diff = numpy.sum(img - test)

        self.assertEqual(total_diff, 0)

        # Cleanup
        shutil.rmtree(temp_dir)


if __name__ == '__main__':
    unittest.main()
