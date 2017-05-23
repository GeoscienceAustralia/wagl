#!/usr/bin/env python

# TODO: move comments into the BLRB code block
"""
Author: Roger Edberg (roger.edberg@ga.gov.au)
Functions for BiLinear Recursive Bisection (BLRB).

All shape references here follow the numpy convention (nrows, ncols), which
makes some of the code harder to follow.
"""

from __future__ import absolute_import
from os.path import basename, splitext
import logging
import numpy
import h5py
from gaip.constants import DatasetName, Model
from gaip.hdf5 import dataset_compression_kwargs
from gaip.hdf5 import write_h5_image
from gaip.hdf5 import read_table
import gaip.interpolate

logger = logging.getLogger(__name__)

DEFAULT_ORIGIN = (0, 0)
DEFAULT_SHAPE = (8, 8)


def bilinear(shape, fUL, fUR, fLR, fLL, dtype=numpy.float64):
    """
    Bilinear interpolation of four scalar values.

    :param shape:
        Shape of interpolated grid (nrows, ncols).

    :param fUL:
        Data value at upper-left (NW) corner.

    :param fUR:
        Data value at upper-right (NE) corner.

    :param fLR:
        Data value at lower-right (SE) corner.

    :param fLL:
        Data value at lower-left (SW) corner.

    :param dtype:
        Data type (numpy I presume?).

    :return:
        Array of data values interpolated between corners.
    """

    s, t = [a.astype(dtype) for a in numpy.ogrid[0:shape[0], 0:shape[1]]]

    s /= (shape[0] - 1.0)
    t /= (shape[1] - 1.0)

    result = (s * (t * fLR + (1.0 - t) * fLL) + (1.0 - s) *
              (t * fUR + (1.0 - t) * fUL))

    return result


def indices(origin=DEFAULT_ORIGIN, shape=DEFAULT_SHAPE):
    """
    Generate corner indices for a grid block.

    :param origin:
        Block origin (2-tuple).

    :param shape:
        Block shape (2-tuple: nrows, ncols).

    :return:
        Corner indices: (xmin, xmax, ymin, ymax).
    """
    return (origin[0], origin[0] + shape[0] - 1,
            origin[1], origin[1] + shape[1] - 1)


def subdivide(origin=DEFAULT_ORIGIN, shape=DEFAULT_SHAPE):
    """
    Generate indices for grid sub-blocks.

    :param origin:
        Block origin (2-tuple).

    :param shape:
        Block shape (nrows, ncols).

    :return:
        Dictionary containing sub-block corner indices:
            { 'UL': <list of 2-tuples>,
              'UR': <list of 2-tuples>,
              'LL': <list of 2-tuples>,
              'LR': <list of 2-tuples> }
    """
    i0, ie, j0, je = indices(origin, shape)
    ic = origin[0] + shape[0] // 2
    jc = origin[1] + shape[1] // 2

    return {
        'UL': [(i0, j0), (i0, jc), (ic, j0), (ic, jc)],
        'LL': [(ic, j0), (ic, jc), (ie, j0), (ie, jc)],
        'UR': [(i0, jc), (i0, je), (ic, jc), (ic, je)],
        'LR': [(ic, jc), (ic, je), (ie, jc), (ie, je)],
    }


def interpolate_block(origin=DEFAULT_ORIGIN, shape=DEFAULT_SHAPE,
                      eval_func=None, grid=None):
    """
    Interpolate a grid block.

    :param origin:
        Block origin (2-tuple).

    :param shape:
        Block shape (nrows, ncols).

    :param eval_func:
        Evaluator function.
    :type eval_func:
        callable; accepts grid indices i, j and returns a scalar value.

    :param grid:
        Grid array.
    :type grid:
        :py:class:`numpy.array`.

    :return:
        Interpolated block array if grid argument is None. If grid argument
        is supplied its elements are modified in place and this function
        does not return a value.
    """
    i0, i1, j0, j1 = indices(origin, shape)

    fUL = eval_func(i0, j0)
    fLL = eval_func(i1, j0)
    fUR = eval_func(i0, j1)
    fLR = eval_func(i1, j1)

    if grid is None:
        return bilinear(shape, fUL, fUR, fLR, fLL)

    grid[i0:i1 + 1, j0:j1 + 1] = bilinear(shape, fUL, fUR, fLR, fLL)


def interpolate_grid(depth=0, origin=DEFAULT_ORIGIN, shape=DEFAULT_SHAPE,
                     eval_func=None, grid=None):
    """
    Interpolate a data grid.

    :param depth:
        Recursive bisection depth.
    :type depth:
        :py:class:`int`

    :param origin:
        Block origin,
    :type origin:
        :py:class:`tuple` of length 2.

    :param shape:
        Block shape.
    :type shape:
        :py:class:`tuple` of length 2 ``(nrows, ncols)``.

    :param eval_func:
        Evaluator function.
    :type eval_func:
        callable; accepts grid indices i, j and returns a scalar value.

    :param grid:
        Grid array.
    :type grid:
        :py:class:`numpy.array`.

    :todo:
        Move arguments ``eval_func`` and ``grid`` to positions 1 and 2, and remove
        defaults (and the check that they are not ``None`` at the top of the function
        body).
    """
    assert eval_func is not None
    assert grid is not None

    if depth == 0:
        interpolate_block(origin, shape, eval_func, grid)
    else:
        blocks = subdivide(origin, shape)
        for (kUL, kUR, kLL, kLR) in blocks.values():
            block_shape = (kLR[0] - kUL[0] + 1, kLR[1] - kUL[1] + 1)
            interpolate_grid(depth - 1, kUL, block_shape, eval_func, grid)


def _bilinear_interpolate(acq, factor, sat_sol_angles_fname,
                          coefficients_fname, ancillary_fname, out_fname,
                          compression, y_tile, method):
    """
    A private wrapper for dealing with the internal custom workings of the
    NBAR workflow.
    """
    with h5py.File(sat_sol_angles_fname, 'r') as sat_sol,\
        h5py.File(coefficients_fname, 'r') as coef,\
        h5py.File(ancillary_fname, 'r') as anc:

        # read the relevant tables into DataFrames
        coord_dset = read_table(anc, DatasetName.coordinator.value)
        centre_dset = read_table(sat_sol, DatasetName.centreline.value)
        box_dset = read_table(sat_sol, DatasetName.boxline.value)

        if factor in Model.nbar.factors:
            dataset_name = DatasetName.nbar_coefficients.value
        elif factor in Model.sbt.factors:
            dataset_name = DatasetName.sbt_coefficients.value
        else:
            msg = "Factor name not found in available factors: {}"
            raise ValueError(msg.format(Model.standard.factors))

        coef_dset = read_table(coef, dataset_name)

        rfid = bilinear_interpolate(acq, factor, coord_dset, box_dset,
                                    centre_dset, coef_dset, out_fname,
                                    compression, y_tile, method)

    rfid.close()
    return


def bilinear_interpolate(acq, factor, coordinator_dataset, boxline_dataset,
                         centreline_dataset, coefficients, out_fname=None,
                         compression='lzf', y_tile=100, method=None):
    # TODO: more docstrings
    """Perform bilinear interpolation."""
    geobox = acq.gridded_geo_box()
    cols, rows = geobox.get_shape_xy()

    coord = numpy.zeros((coordinator_dataset.shape[0], 2), dtype='int')
    map_x = coordinator_dataset.map_x.values
    map_y = coordinator_dataset.map_y.values
    coord[:, 1], coord[:, 0] = (map_x, map_y) * ~geobox.transform
    centre = boxline_dataset.bisection_index.values + 1
    start = boxline_dataset.start_index.values + 1
    end = boxline_dataset.end_index.values + 1

    band = acq.band_num
    band_records = coefficients.band_id == 'BAND {}'.format(band)
    samples = coefficients[factor][band_records].values

    if method is None:
        if len(samples) == 9:
            method = 'linear'
        else:
            method = 'shear'

    func_map = {'linear': gaip.interpolate.fortran_bilinear_interpolate,
                'shear': gaip.interpolate.sheared_bilinear_interpolate,
                'rbf': gaip.interpolate.rbf_interpolate}
    assert method in func_map

    result = func_map[method](cols, rows, coord, samples, start, end, centre)

    # Initialise the output files
    if out_fname is None:
        fid = h5py.File('bilinear.h5', driver='core',
                        backing_store=False)
    else:
        fid = h5py.File(out_fname, 'w')

    # TODO: determine without splitext or basename
    dset_name = splitext(basename(out_fname))[0]
    kwargs = dataset_compression_kwargs(compression=compression,
                                        chunks=(1, geobox.x_size()))
    no_data = -999
    kwargs['fillvalue'] = no_data
    attrs = {'crs_wkt': geobox.crs.ExportToWkt(),
             'geotransform': geobox.transform.to_gdal(),
             'no_data_value': no_data}
    desc = ("Contains the bi-linearly interpolated result of factor {}"
            "for band {} from sensor {}.")
    attrs['Description'] = desc.format(factor, band, acq.satellite_name)
    write_h5_image(result, dset_name, fid, attrs, **kwargs)

    return fid


def link_bilinear_data(data, out_fname):
    """
    Links the individual bilinearly interpolated results into a
    single file for easier access.
    """
    for key in data:
        # band, factor = key
        fname = data[key]
        base_dname = splitext(basename(fname))[0]

        # do we need two group levels?
        # dset_name = ppjoin(band, factor, base_dname)

        with h5py.File(out_fname, 'a') as fid:
            fid[base_dname] = h5py.ExternalLink(fname, base_dname)

    return
