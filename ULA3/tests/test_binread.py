import numpy as np
import unittest
import tempfile
import os

from os.path import join as pjoin
from numpy.testing import assert_array_equal
from ULA3.utils import load_bin_file

class TestBinaryRead(unittest.TestCase):

    def test_int8_tofile(self):
        fn = tempfile.NamedTemporaryFile().name
        orig = np.arange(1, 41, dtype=np.int8).reshape((10,4))
        orig.tofile(fn)
        data = load_bin_file(fn, 10, 4, np.int8)
        assert_array_equal(orig, data)
        os.remove(fn)

    def test_int8_fromfile(self):
        fn = tempfile.NamedTemporaryFile().name
        orig = np.arange(1, 41, dtype=np.int8).reshape((10,4))
        orig.tofile(fn)
        data = np.fromfile(fn, np.int8).reshape((10,4))
        assert_array_equal(orig, data)
        os.remove(fn)

    def test_int8_fortran(self):
        fn = pjoin('data', 'int8.bin')
        orig = np.arange(1, 41, dtype=np.int8).reshape((10,4))
        data = load_bin_file(fn, 10, 4, np.int8)
        assert_array_equal(orig, data)

    def test_int16_tofile(self):
        fn = tempfile.NamedTemporaryFile().name
        orig = np.arange(1, 41, dtype=np.int16).reshape((10,4))
        orig.tofile(fn)
        data = load_bin_file(fn, 10, 4, np.int16)
        assert_array_equal(orig, data)
        os.remove(fn)

    def test_int16_fromfile(self):
        fn = tempfile.NamedTemporaryFile().name
        orig = np.arange(1, 41, dtype=np.int16).reshape((10,4))
        orig.tofile(fn)
        data = np.fromfile(fn, np.int16).reshape((10,4))
        assert_array_equal(orig, data)
        os.remove(fn)

    def test_int16_fortran(self):
        fn = pjoin('data', 'int16.bin')
        orig = np.arange(1, 41, dtype=np.int16).reshape((10,4))
        data = load_bin_file(fn, 10, 4, np.int16)
        assert_array_equal(orig, data)

    def test_int32_tofile(self):
        fn = tempfile.NamedTemporaryFile().name
        orig = np.arange(1, 41, dtype=np.int32).reshape((10,4))
        orig.tofile(fn)
        data = load_bin_file(fn, 10, 4, np.int32)
        assert_array_equal(orig, data)
        os.remove(fn)

    def test_int32_fromfile(self):
        fn = tempfile.NamedTemporaryFile().name
        orig = np.arange(1, 41, dtype=np.int32).reshape((10,4))
        orig.tofile(fn)
        data = np.fromfile(fn, np.int32).reshape((10,4))
        assert_array_equal(orig, data)
        os.remove(fn)

    def test_int32_fortran(self):
        fn = pjoin('data', 'int32.bin')
        orig = np.arange(1, 41, dtype=np.int32).reshape((10,4))
        data = load_bin_file(fn, 10, 4, np.int32)
        assert_array_equal(orig, data)

    def test_float32_tofile(self):
        fn = tempfile.NamedTemporaryFile().name
        orig = np.arange(1, 41, dtype=np.float32).reshape((10,4))
        orig.tofile(fn)
        data = load_bin_file(fn, 10, 4, np.float32)
        assert_array_equal(orig, data)
        os.remove(fn)

    def test_float32_fromfile(self):
        fn = tempfile.NamedTemporaryFile().name
        orig = np.arange(1, 41, dtype=np.float32).reshape((10,4))
        orig.tofile(fn)
        data = np.fromfile(fn, np.float32).reshape((10,4))
        assert_array_equal(orig, data)
        os.remove(fn)

    def test_float32_fortran(self):
        fn = pjoin('data', 'float32.bin')
        orig = np.arange(1, 41, dtype=np.float32).reshape((10,4))
        data = load_bin_file(fn, 10, 4, np.float32)
        assert_array_equal(orig, data)

if __name__ == '__main__':
    unittest.main()
