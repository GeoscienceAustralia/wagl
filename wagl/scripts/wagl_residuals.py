#!/usr/bin/env python

"""
A recursive utility that compares and evaluates datasets within a
reference file, with the equivalent dataset from a test file.
"""

from __future__ import print_function

import argparse
from functools import partial

from posixpath import join as ppjoin
from posixpath import basename as pbasename
import numpy
import h5py
import pandas

from idl_functions import histogram
from wagl.hdf5 import read_scalar, read_h5_table, write_dataframe
from wagl.hdf5 import write_h5_image, write_h5_table, write_scalar
from wagl.hdf5 import H5CompressionFilter, VLEN_STRING, find
from wagl.geobox import GriddedGeoBox


def distribution(data):
    """
    Evaluates the distribution of a `NumPy` array. Floating point
    arrays will use 256 bins, while integer arrays will use
    a binsize of 1.
    """
    if data.dtype.name in ['float', 'float32', 'float64']:
        try:
            h = histogram(data, omin='omin', omax='omax', locations='loc',
                          nbins=256, nan=True)
        except ValueError:
            h = {}
            h['histogram'] = numpy.zeros((256), dtype='uint32')
            h['omin'] = 0
            h['omax'] = 0
            h['loc'] = numpy.zeros((256), dtype=data.dtype.name)
            h['histogram'][0] = data.size
    else:
        h = histogram(data, omin='omin', omax='omax', locations='loc')

    return h


def image_residual(ref_fid, test_fid, pathname, out_fid,
                   compression=H5CompressionFilter.LZF, save_inputs=False,
                   filter_opts=None):
    """
    Undertake residual analysis for IMAGE CLASS Datasets.
    A histogram and a cumulative histogram of the residuals are
    calculated and recorded as TABLE CLASS Datasets.
    Any NaN's in IMAGE datasets will be handled automatically.

    :param ref_fid:
        A h5py file object (essentially the root Group), containing
        the reference data.

    :param test_fid:
        A h5py file object (essentially the root Group), containing
        the test data.

    :param pathname:
        A `str` containing the pathname to the IMAGE Dataset.

    :param out_fid:
        A h5py file object (essentially the root Group), opened for
        writing the output data.

    :param compression:
        The compression filter to use.
        Default is H5CompressionFilter.LZF 

    :param save_inputs:
        A `bool` indicating whether or not to save the input datasets
        used for evaluating the residuals alongside the results.
        Default is False.

    :filter_opts:
        A dict of key value pairs available to the given configuration
        instance of H5CompressionFilter. For example
        H5CompressionFilter.LZF has the keywords *chunks* and *shuffle*
        available.
        Default is None, which will use the default settings for the
        chosen H5CompressionFilter instance.

    :return:
        None; This routine will only return None or a print statement,
        this is essential for the HDF5 visit routine.
    """
    def evaluate(ref_dset, test_dset):
        """
        Evaluate the image residual.
        Caters for boolean types.
        TODO: geobox intersection if dimensions are different.
        TODO: handle no data values
        TODO: handle classification datasets
        TODO: handle bitwise datasets
        """
        if ref_dset.dtype.name == 'bool':
            result = numpy.logical_xor(ref_dset, test_dset).astype('uint8')
        else:
            result = ref_dset[:] - test_dset
        return result

    class_name = 'IMAGE'
    ref_dset = ref_fid[pathname]
    test_dset = test_fid[pathname]

    # ignore no data values for the time being
    residual = evaluate(ref_dset, test_dset)
    min_residual = numpy.nanmin(residual)
    max_residual = numpy.nanmax(residual)
    pct_difference = (residual != 0).sum() / residual.size * 100

    if filter_opts is None:
        filter_opts = {}
    else:
        fopts = filter_opts.copy()
    fopts['chunks'] = ref_dset.chunks

    geobox = GriddedGeoBox.from_dataset(ref_dset)

    # output residual
    attrs = {
        'crs_wkt': geobox.crs.ExportToWkt(),
        'geotransform': geobox.transform.to_gdal(),
        'description': 'Residual',
        'min_residual': min_residual,
        'max_residual': max_residual,
        'percent_difference': pct_difference
        }

    base_dname = pbasename(pathname)
    group_name = ref_dset.parent.name.strip('/')
    dname = ppjoin('RESULTS', class_name, 'RESIDUALS', group_name, base_dname)
    write_h5_image(residual, dname, out_fid, compression, attrs, fopts)

    # residuals distribution
    h = distribution(residual)
    hist = h['histogram']

    attrs = {
        'description': 'Frequency distribution of the residuals',
        'omin': h['omin'],
        'omax': h['omax']
        }
    dtype = numpy.dtype([('bin_locations', h['loc'].dtype.name),
                         ('residuals_distribution', hist.dtype.name)])
    table = numpy.zeros(hist.shape, dtype=dtype)
    table['bin_locations'] = h['loc']
    table['residuals_distribution'] = hist

    # output
    fopts.drop('chunks')
    dname = ppjoin('RESULTS', class_name, 'FREQUENCY-DISTRIBUTIONS',
                   group_name, base_dname)
    write_h5_table(table, dname, out_fid, compression, attrs=attrs,
                   filter_opts=fopts)

    # cumulative distribution
    h = distribution(numpy.abs(residual))
    hist = h['histogram']
    cdf = numpy.cumsum(hist / hist.sum())

    attrs = {
        'description': 'Cumulative distribution of the residuals',
        'omin': h['omin'],
        'omax': h['omax'],
        '90th_percentile': h['loc'][numpy.searchsorted(cdf, 0.9)],
        '99th_percentile': h['loc'][numpy.searchsorted(cdf, 0.99)]
        }
    dtype = numpy.dtype([('bin_locations', h['loc'].dtype.name),
                         ('cumulative_distribution', cdf.dtype.name)])
    table = numpy.zeros(cdf.shape, dtype=dtype)
    table['bin_locations'] = h['loc']
    table['cumulative_distribution'] = cdf

    # output
    dname = ppjoin('RESULTS', class_name, 'CUMULATIVE-DISTRIBUTIONS',
                   group_name, base_dname)
    write_h5_table(table, dname, out_fid, compression=compression, attrs=attrs,
                   filter_opts=fopts)

    if save_inputs:
        # copy the reference data
        out_grp = out_fid.require_group(ppjoin('REFERENCE-DATA', group_name))
        ref_fid.copy(ref_dset, out_grp)

        # copy the test data
        out_grp = out_fid.require_group(ppjoin('TEST-DATA', group_name))
        test_fid.copy(test_dset, out_grp)


def scalar_residual(ref_fid, test_fid, pathname, out_fid, save_inputs):
    """
    Undertake a simple equivalency test, rather than a numerical
    difference. This allows strings to be compared.

    :param ref_fid:
        A h5py file object (essentially the root Group), containing
        the reference data.

    :param test_fid:
        A h5py file object (essentially the root Group), containing
        the test data.

    :param pathname:
        A `str` containing the pathname to the SCALAR Dataset.

    :param out_fid:
        A h5py file object (essentially the root Group), opened for
        writing the output data.

    :param save_inputs:
        A `bool` indicating whether or not to save the input datasets
        used for evaluating the residuals alongside the results.
        Default is False.

    :return:
        None; This routine will only return None or a print statement,
        this is essential for the HDF5 visit routine.
    """
    class_name = 'SCALAR'
    ref_data = read_scalar(ref_fid, pathname)
    test_data = read_scalar(test_fid, pathname)

    # copy the attrs
    attrs = ref_data.copy()
    attrs.pop('value')
    attrs['description'] = 'Equivalency Test'

    # drop 'file_format' as the conversion tool will try to output that format
    # but currently we're not testing contents, just if it is different
    # so saying we've created a yaml string when it is a simple bool is
    # not correct
    attrs.pop('file_format', None)

    # this'll handle string types, but we won't get a numerical
    # difference value for numerical values, only a bool
    diff = ref_data['value'] == test_data['value']

    # output
    base_dname = pbasename(pathname)
    group_name = ref_fid[pathname].parent.name.strip('/')
    dname = ppjoin('RESULTS', class_name, 'EQUIVALENCY', group_name, base_dname)
    write_scalar(diff, dname, out_fid, attrs)

    if save_inputs:
        # copy the reference data
        out_grp = out_fid.require_group(ppjoin('REFERENCE-DATA', group_name))
        ref_fid.copy(ref_fid[pathname], out_grp)

        # copy the test data
        out_grp = out_fid.require_group(ppjoin('TEST-DATA', group_name))
        test_fid.copy(test_fid[pathname], out_grp)


def table_residual(ref_fid, test_fid, pathname, out_fid,
                   compression=H5CompressionFilter.LZF, save_inputs=False,
                   filter_opts=None):
    """
    Output a residual TABLE of the numerical columns, ignoring
    columns with the dtype `object`.
    An equivalency test using `pandas.DataFrame.equals` is also
    undertaken which if False, requires further investigation to
    determine the column(s) and row(s) that are different.

    :param ref_fid:
        A h5py file object (essentially the root Group), containing
        the reference data.

    :param test_fid:
        A h5py file object (essentially the root Group), containing
        the test data.

    :param pathname:
        A `str` containing the pathname to the TABLE Dataset.

    :param out_fid:
        A h5py file object (essentially the root Group), opened for
        writing the output data.

    :param compression:
        The compression filter to use.
        Default is H5CompressionFilter.LZF 

    :param save_inputs:
        A `bool` indicating whether or not to save the input datasets
        used for evaluating the residuals alongside the results.
        Default is False.

    :filter_opts:
        A dict of key value pairs available to the given configuration
        instance of H5CompressionFilter. For example
        H5CompressionFilter.LZF has the keywords *chunks* and *shuffle*
        available.
        Default is None, which will use the default settings for the
        chosen H5CompressionFilter instance.

    :return:
        None; This routine will only return None or a print statement,
        this is essential for the HDF5 visit routine.
    """
    class_name = 'TABLE'
    ref_df = read_h5_table(ref_fid, pathname)
    test_df = read_h5_table(test_fid, pathname)

    # ignore any `object` dtype columns (mostly just strings)
    cols = [col for col in ref_df.columns if
            ref_df[col].dtype.name != 'object']

    # difference and pandas.DataFrame.equals test
    df = ref_df[cols] - test_df[cols]
    equal = test_df.equals(ref_df)

    # ignored cols
    cols = [col for col in ref_df.columns if
            ref_df[col].dtype.name == 'object']

    # output
    attrs = {'description': 'Residuals of numerical columns only',
             'columns_ignored': numpy.array(cols, VLEN_STRING),
             'equivalent': equal}
    base_dname = pbasename(pathname)
    group_name = ref_fid[pathname].parent.name.strip('/')
    dname = ppjoin('RESULTS', class_name, 'RESIDUALS', group_name, base_dname)
    write_dataframe(df, dname, out_fid, compression, attrs=attrs,
                    filter_opts=filter_opts)

    if save_inputs:
        # copy the reference data
        out_grp = out_fid.require_group(ppjoin('REFERENCE-DATA', group_name))
        ref_fid.copy(ref_fid[pathname], out_grp)

        # copy the test data
        out_grp = out_fid.require_group(ppjoin('TEST-DATA', group_name))
        test_fid.copy(test_fid[pathname], out_grp)


def residuals(ref_fid, test_fid, out_fid, compression=H5CompressionFilter.LZF,
              save_inputs=False, filter_opts=None, pathname=''):
    """
    Undertake residual analysis for each dataset found within the
    reference file.
    Currently only 3 Dataset CLASS types are evaluated:
        * IMAGE
        * TABLE
        * SCALAR

    :param ref_fid:
        A h5py file object (essentially the root Group), containing
        the reference data.

    :param test_fid:
        A h5py file object (essentially the root Group), containing
        the test data.

    :param pathname:
        A `str` containing the pathname to the Dataset.

    :param out_fid:
        A h5py file object (essentially the root Group), opened for
        writing the output data.

    :param compression:
        The compression filter to use.
        Default is H5CompressionFilter.LZF 

    :param save_inputs:
        A `bool` indicating whether or not to save the input datasets
        used for evaluating the residuals alongside the results.
        Default is False.

    :filter_opts:
        A dict of key value pairs available to the given configuration
        instance of H5CompressionFilter. For example
        H5CompressionFilter.LZF has the keywords *chunks* and *shuffle*
        available.
        Default is None, which will use the default settings for the
        chosen H5CompressionFilter instance.

    :return:
        None; This routine will only return None or a print statement,
        this is essential for the HDF5 visit routine.
    """
    pathname = pathname.decode('utf-8')
    if pathname in test_fid:
        obj = ref_fid[pathname]
        if isinstance(obj, h5py.Group):
            # skip any Group objects for now
            return None

        class_name = obj.attrs.get('CLASS')
        if class_name == 'IMAGE':
            image_residual(ref_fid, test_fid, pathname, out_fid, compression,
                           save_inputs, filter_opts)
        elif class_name == 'SCALAR':
            scalar_residual(ref_fid, test_fid, pathname, out_fid, save_inputs)
        elif class_name == 'TABLE':
            table_residual(ref_fid, test_fid, pathname, out_fid, compression,
                           save_inputs, filter_opts)
        else:
            # potentially datasets with no assigned class attribute
            # but good to know what they are so they can be assigned one
            print("Skipping {}\t{}".format(pathname, type(obj)))
    else:
        print('{} not found in test file...Skipping...'.format(pathname))


def image_results(image_group, compression=H5CompressionFilter.LZF,
                  filter_opts=None):
    """
    Combine the residual results of each IMAGE Dataset into a
    single TABLE Dataset.
    """
    # potentially could just use visit...
    img_paths = find(image_group, 'IMAGE')

    min_ = []
    max_ = []
    percent = []
    pct_90 = []
    pct_99 = []
    resid_paths = []
    hist_paths = []
    chist_paths = []
    products = []
    name = []

    for pth in img_paths:
        hist_pth = pth.replace('RESIDUALS', 'FREQUENCY-DISTRIBUTIONS')
        chist_pth = pth.replace('RESIDUALS', 'CUMULATIVE-DISTRIBUTIONS')
        resid_paths.append(ppjoin(image_group.name, pth))
        hist_paths.append(ppjoin(image_group.name, hist_pth))
        chist_paths.append(ppjoin(image_group.name, chist_pth))

        dset = image_group[pth]
        min_.append(dset.attrs['min_residual'])
        max_.append(dset.attrs['max_residual'])
        percent.append(dset.attrs['percent_difference'])
        products.append(pbasename(dset.parent.name))
        name.append(pbasename(dset.name))

        dset = image_group[chist_pth]
        pct_90.append(dset.attrs['90th_percentile'])
        pct_99.append(dset.attrs['99th_percentile'])

    df = pandas.DataFrame({'product': products,
                           'dataset_name': name,
                           'min_residual': min_,
                           'max_residual': max_,
                           'percent_difference': percent,
                           '90th_percentile': pct_90,
                           '99th_percentile': pct_99,
                           'residual_image_pathname': resid_paths,
                           'residual_histogram_pathname': hist_paths,
                           'residual_cumulative_pathname': chist_paths})

    # output
    write_dataframe(df, 'IMAGE-RESIDUALS', image_group, compression,
                    title='RESIDUALS-TABLE', filter_opts=filter_opts)


def scalar_results(scalar_group, compression=H5CompressionFilter.LZF,
                   filter_opts=None):
    """
    Combine the residual results of each SCALAR Dataset into a
    single TABLE Dataset.
    """
    # potentially could just use visit...
    paths = find(scalar_group, 'SCALAR')

    equivalent = []
    products = []
    name = []

    for pth in paths:
        dset = scalar_group[pth]
        equivalent.append(dset[()])
        products.append(pbasename(dset.parent.name))
        name.append(pbasename(dset.name))

    df = pandas.DataFrame({'product': products,
                           'dataset_name': name,
                           'equivalent': equivalent})

    # output
    write_dataframe(df, 'SCALAR-EQUIVALENCY', scalar_group, compression,
                    title='EQUIVALENCY-RESULTS', filter_opts=filter_opts)


def table_results(table_group, compression=H5CompressionFilter.LZF,
                  filter_opts=None):
    """
    Combine the residual results of each TABLE Dataset into a
    single TABLE Dataset.
    """
    # potentially could just use visit...
    paths = find(table_group, 'TABLE')

    equivalent = []
    products = []
    name = []

    for pth in paths:
        dset = table_group[pth]
        equivalent.append(dset.attrs['equal'])
        products.append(pbasename(dset.parent.name))
        name.append(pbasename(dset.name))

    df = pandas.DataFrame({'product': products,
                           'dataset_name': name,
                           'equivalent': equivalent})

    # output
    write_dataframe(df, 'TABLE-EQUIVALENCY', table_group, compression,
                    title='EQUIVALENCY-RESULTS', filter_opts=filter_opts)


def run(reference_fname, test_fname, out_fname, compression, save_inputs,
        filter_opts):
    """ Run the residuals analysis. """
    # note: lower level h5py access is required in order to visit links
    with h5py.File(reference_fname, 'r') as ref_fid:
        with h5py.File(test_fname, 'r') as test_fid:
            with h5py.File(out_fname, 'w') as out_fid:
                root = h5py.h5g.open(ref_fid.id, b'/')
                root.links.visit(partial(residuals, ref_fid, test_fid,
                                         out_fid, compression, save_inputs))

                # create singular TABLES for each Dataset CLASS

                # IMAGE
                grp = out_fid[ppjoin('RESULTS', 'IMAGE')]
                image_results(grp, compression, filter_opts)

                # SCALAR
                grp = out_fid[ppjoin('RESULTS', 'SCALAR')]
                scalar_results(grp)

                # TABLE
                grp = out_fid[ppjoin('RESULTS', 'TABLE')]
                table_results(grp)


def _parser():
    """ Argument parser. """
    description = ("Undertake residual analysis between a reference file "
                   "test file.")
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--test-filename", required=True,
                        help=("The filename of the file containing the test "
                              "datasets."))
    parser.add_argument("--reference-filename", required=True,
                        help=("The filename of the file containing the "
                              "reference datasets."))
    parser.add_argument("--out-filename", required=True,
                        help=("The filename of the file to contain the "
                              "results."))
    parser.add_argument("--compression", default="lzf",
                        help="The comression filter to use.")
    parser.add_argument("--save-inputs", action='store_true',
                        help=("Save the reference and test datasets "
                              "alongside the residual/difference datasets."))

    return parser


def main():
    """ Main execution. """
    parser = _parser()
    args = parser.parse_args()
    run(args.test_filename, args.reference_filename, args.out_filename,
        args.compression, args.save_inputs, args.filter_opts)
