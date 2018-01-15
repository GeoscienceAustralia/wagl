#!/usr/bin/env python

"""
Test the various utilites contained in the wagl.hdf5 module.
"""

import datetime
import unittest
import numpy
import h5py
import pandas
from wagl import hdf5


class HDF5Test(unittest.TestCase):

    """
    Test the various utilites contained in the wagl.hdf5 module.
    """

    scalar_data = 66
    image_data = numpy.random.randint(0, 256, (10, 10))
    table_dtype = numpy.dtype([('float_data', 'float64'),
                               ('integer_data', 'int64')])
    table_data = numpy.zeros((10), dtype=table_dtype)
    table_data['float_data'] = numpy.random.ranf(10)
    table_data['integer_data'] = numpy.random.randint(0, 10001, (10))

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
        self.assertDictEqual(kwargs, self.bitshuffle_kwargs)

    def test_write_scalar(self):
        """
        Test the write_scalar function.
        """
        data = self.scalar_data
        fname = 'test_write_scalar.h5'
        with h5py.File(fname, **self.memory_kwargs) as fid:
            self.assertIsNone(hdf5.write_scalar(data, 'scalar', fid))

    def test_scalar_attributes(self):
        """
        Test the scalar attributes.
        """
        attrs = {'test_attribute': 'this is a scalar'}
        data = {'value': self.scalar_data,
                'CLASS': 'SCALAR',
                'VERSION': '0.1'}

        # insert the attribute into the data dict
        for k, v in attrs.items():
            data[k] = v

        fname = 'test_scalar_dataset.h5'
        with h5py.File(fname, **self.memory_kwargs) as fid:
            hdf5.write_scalar(data['value'], 'test-scalar', fid, attrs=attrs)

            self.assertDictEqual(hdf5.read_scalar(fid, 'test-scalar'), data)

    def test_datetime_attrs(self):
        """
        Test that datetime objects will be converted to iso format
        when writing attributes.
        """
        attrs = {'timestamp': datetime.datetime.now()}

        fname = 'test_datetime_attrs.h5'
        with h5py.File(fname, **self.memory_kwargs) as fid:
            hdf5.write_scalar(self.scalar_data, 'scalar', fid, attrs=attrs)

            data = hdf5.read_scalar(fid, 'scalar')
            self.assertEqual(data['timestamp'], attrs['timestamp'])

    def test_attach_attributes(self):
        """
        Test the attach_attributes function.
        """
        attrs = {'alpha': 1, 'beta': 2}

        fname = 'test_attach_attributes.h5'
        with h5py.File(fname, **self.memory_kwargs) as fid:
            dset = fid.create_dataset('data', data=self.image_data)
            hdf5.attach_attributes(dset, attrs)
            test = {k: v for k, v in dset.attrs.items()}
            self.assertDictEqual(test, attrs)

    def test_write_h5_image(self):
        """
        Test the write_h5_image function.
        """
        data = self.image_data
        fname = 'test_write_h5_image.h5'
        with h5py.File(fname, **self.memory_kwargs) as fid:
            self.assertIsNone(hdf5.write_h5_image(data, 'image', fid))

    def test_write_h5_table(self):
        """
        Test the write_h5_table function.
        """
        data = self.table_data
        fname = 'test_write_h5_table.h5'
        with h5py.File(fname, **self.memory_kwargs) as fid:
            self.assertIsNone(hdf5.write_h5_table(data, 'table', fid))

    def test_attach_image_attributes(self):
        """
        Test the attach_image_attributes function.
        """
        attrs = {'CLASS': 'IMAGE',
                 'IMAGE_VERSION': '1.2',
                 'DISPLAY_ORIGIN': 'UL'}

        fname = 'test_attach_image_attributes.h5'
        with h5py.File(fname, **self.memory_kwargs) as fid:
            dset = fid.create_dataset('data', data=self.image_data)
            hdf5.attach_image_attributes(dset, attrs)
            test = {k: v for k, v in dset.attrs.items()}
            self.assertDictEqual(test, attrs)

    def test_write_h5_image_attributes(self):
        """
        Test the image attributes of the write_h5_image function.
        """
        attrs = {'CLASS': 'IMAGE',
                 'IMAGE_VERSION': '1.2',
                 'DISPLAY_ORIGIN': 'UL'}

        fname = 'test_write_h5_image_attributes.h5'
        with h5py.File(fname, **self.memory_kwargs) as fid:
            hdf5.write_h5_image(self.image_data, 'image', fid)
            test = {k: v for k, v in fid['image'].attrs.items()}

            # assertDictEqual can't compare a numpy array, so test elsewhere
            del test['IMAGE_MINMAXRANGE']

            self.assertDictEqual(test, attrs)

    def test_write_h5_image_minmax(self):
        """
        Test the IMAGE_MINMAXRANGE attribute is correct.
        """
        minmax = numpy.array([self.image_data.min(), self.image_data.max()])

        fname = 'test_write_h5_image_minmax.h5'
        with h5py.File(fname, **self.memory_kwargs) as fid:
            hdf5.write_h5_image(self.image_data, 'image', fid)

            test = fid['image'].attrs['IMAGE_MINMAXRANGE']

            self.assertTrue((minmax == test).all())

    def test_attach_table_attributes(self):
        """
        Test the attach_table_attributes function.
        """
        attrs = {'CLASS': 'TABLE',
                 'VERSION': '0.2',
                 'TITLE': 'Table',
                 'FIELD_0_NAME': 'float_data',
                 'FIELD_1_NAME': 'integer_data'}

        fname = 'test_attach_table_attributes.h5'
        with h5py.File(fname, **self.memory_kwargs) as fid:
            dset = fid.create_dataset('data', data=self.table_data)
            hdf5.attach_table_attributes(dset, attrs=attrs)
            test = {k: v for k, v in dset.attrs.items()}
            self.assertDictEqual(test, attrs)

    def test_write_dataframe(self):
        """
        Test the write_dataframe function.
        """
        df = pandas.DataFrame(self.table_data)
        fname = 'test_write_dataframe.h5'
        with h5py.File(fname, **self.memory_kwargs) as fid:
            self.assertIsNone(hdf5.write_dataframe(df, 'dataframe', fid))

    def test_dataframe_attributes(self):
        """
        Test the attributes that get created for a dataframe.
        """
        attrs = {'CLASS': 'TABLE',
                 'FIELD_0_NAME': 'index',
                 'FIELD_1_NAME': 'float_data',
                 'FIELD_2_NAME': 'integer_data',
                 'TITLE': 'Table',
                 'VERSION': '0.2',
                 'float_data_dtype': 'float64',
                 'index_dtype': 'int64',
                 'index_names': numpy.array(['index'], dtype=object),
                 'integer_data_dtype': 'int64',
                 'metadata': '`Pandas.DataFrame` converted to HDF5 compound datatype.', # pylint: disable=line-too-long
                 'nrows': 10,
                 'python_type': '`Pandas.DataFrame`'}

        df = pandas.DataFrame(self.table_data)

        fname = 'test_dataframe_attributes.h5'
        with h5py.File(fname, **self.memory_kwargs) as fid:
            hdf5.write_dataframe(df, 'dataframe', fid)

            test = {k: v for k, v in fid['dataframe'].attrs.items()}
            self.assertDictEqual(test, attrs)

    def test_dataframe_roundtrip(self):
        """
        Test that the pandas dataframe roundtrips, i.e. save to HDF5
        and is read back into a dataframe seamlessly.
        Float, integer, datetime and string datatypes will be
        tested.
        """
        df = pandas.DataFrame(self.table_data)
        df['timestamps'] = pandas.date_range('1/1/2000', periods=10, freq='D')
        df['string_data'] = ['period {}'.format(i) for i in range(10)]

        fname = 'test_dataframe_roundtrip.h5'
        with h5py.File(fname, **self.memory_kwargs) as fid:
            hdf5.write_dataframe(df, 'dataframe', fid)
            self.assertTrue(df.equals(hdf5.read_h5_table(fid, 'dataframe')))


if __name__ == '__main__':
    unittest.main()
