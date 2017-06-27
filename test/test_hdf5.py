#!/usr/bin/env python

import datetime
import unittest
import numpy
import h5py
import pandas
from gaip import hdf5


class HDF5Test(unittest.TestCase):

    default_kwargs = {'compression': 'lzf',
                      'shuffle': True,
                      'chunks': True,
                      'compression_opts': None}

    mafisc_kwargs = {'compression': 32002,
                     'shuffle': False,
                     'chunks': True,
                     'compression_opts': (1, 0)}

    bitshuffle_kwargs = {'compression': 32008,
                         'shuffle': False,
                         'chunks': True,
                         'compression_opts': (0, 2)}

    memory_kwargs = {'driver': 'core', 'backing_store': False}

    def test_default_filter(self):
        """
        Test the default compression keyword settings.
        """
        kwargs = hdf5.dataset_compression_kwargs()
        self.assertDictEqual(kwargs, self.default_kwargs)

    def test_mafisc_filter(self):
        """
        Test the mafisc compression keyword settings.
        """
        kwargs = hdf5.dataset_compression_kwargs(compression='mafisc')
        self.assertDictEqual(kwargs, self.mafisc_kwargs)

    def test_bitshuffle_filter(self):
        """
        Test the bitshuffle compression keyword settings.
        """
        kwargs = hdf5.dataset_compression_kwargs(compression='bitshuffle')
        self.assertEqual(kwargs, self.bitshuffle_kwargs)

    def test_scalar(self):
        """
        Test the read and write functionality for scalar datasets.
        """
        attrs = {'test_attribute': 'this is a scalar'}
        data = {'value': 66,
                'CLASS': 'SCALAR',
                'VERSION': '0.1'}

        # insert the attribute into the data dict
        for k, v in data.items():
            data[k] = v

        with h5py.File('test_scalar.h5', **self.memory_kwargs) as fid:
            hdf5.write_scalar(data['value'], 'test-scalar', fid, attrs=attrs)

            self.assertDictEqual(hdf5.read_scalar(fid, 'test-scalar'), data)

    def test_datetime_attrs(self):
        """
        Test that datetime objects will be converted to iso format
        when writing attributes.
        """
        attrs = {'timestamp': datetime.datetime.now()}

        with h5py.File('datetime_attrs.h5', **self.memory_kwargs) as fid:
            hdf5.write_scalar(66.0, 'scalar', fid, attrs=attrs)

            data = hdf5.read_scalar(fid, 'scalar')
            self.assertEqual(data['timestamp'], attrs['timestamp'].isoformat())

    def test_attach_attributes(self):
        """
        Test the attach_attributes function.
        """
        data = numpy.random.randint(0, 256, (10, 10))
        attrs = {'alpha': 1, 'beta': 2}
        with h5py.File('attach_attributes.h5', **self.memory_kwargs) as fid:
            dset = fid.create_dataset('data', data=data)
            hdf5.attach_attributes(dset, attrs)
            test = {k: v for k, v in dset.attrs.items()}
            self.assertEqual(test, attrs)


if __name__ == '__main__':
    unittest.main()
