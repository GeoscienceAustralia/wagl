#!/usr/bin/env python

"""
BRDF data extraction utilities
------------------------------

The :ref:`nbar-algorithm-label` and :ref:`tc-algorithm-label` algorithms
require estimates of various atmospheric parameters, which are produced using
`MODTRAN <http://modtran5.com/>`_. MODTRAN, in turn, requires `BRDF
<http://en.wikipedia.org/wiki/Bidirectional_reflectance_distribution_function>`_
estimates. The estimates used in the ULA, are based on `MODIS
<http://modis.gsfc.nasa.gov/>`_ and are produced by CSIRO. For more
information, on how these are used, see :download:`this
<auxiliary/li_etal_2010_05422912.pdf>`.

`MODIS <http://modis.gsfc.nasa.gov/>`_, pre Feb 2001, MODIS data was not
available and an alternative method of deriving `BRDF
<http://en.wikipedia.org/wiki/Bidirectional_reflectance_distribution_function>`_
estimates is required.

"""

from __future__ import absolute_import, print_function
import datetime
import logging
import os
from os.path import join as pjoin

import numpy as np
import rasterio
from rasterio.features import rasterize
from rasterio.crs import CRS
import pyproj
import h5py
from osgeo import ogr
import shapely
import shapely.affinity
import shapely.geometry
from shapely.geometry import box
from shapely import wkt, ops
from wagl.constants import BrdfDirectionalParameters, BrdfModelParameters
from wagl.hdf5 import H5CompressionFilter
from wagl.metadata import current_h5_metadata


log = logging.getLogger('root.' + __name__)


class BRDFLoaderError(Exception):
    """
    BRDF Loader Error
    """


class BRDFLookupError(Exception):
    """
    BRDF Lookup Error
    """


def _date_proximity(cmp_date, date_interpreter=lambda x: x):
    """_date_proximity providers a comparator for an interable
    with an interpreter function. Used to find the closest item
    in a list.

    If two dates are equidistant return the most recent.

    :param cmp_date: date to compare list against
    :param date_interprater: function applied to the list to
        transform items into dates
    """
    def _proximity_comparator(date):
        _date = date_interpreter(date)
        return (
            abs(_date - cmp_date),
            -1 * _date.year,
            -1 * _date.month,
            -1 * _date.day
        )

    return _proximity_comparator


def get_brdf_dirs_modis(brdf_root, scene_date, pattern='%Y.%m.%d'):
    """
    Get list of MODIS BRDF directories for the dataset.

    :param brdf_root:
        BRDF root directory.
    :type brdf_root:
        :py:class:`str`

    :param scene_date:
        Scene Date.
    :type scene_date:
        :py:class:`datetime.date`

    :param pattern:
        A string handed to strptime to interpret directory names into
        observation dates for the brdf ancillary.
    :type pattern:
        :py:class:`str`

    :return:
       A string containing the closest matching BRDF directory.

    """

    dirs = []
    for dname in sorted(os.listdir(brdf_root)):
        try:
            dirs.append(datetime.datetime.strptime(dname, pattern).date())
        except ValueError:
            pass  # Ignore directories that don't match specified pattern

    return min(dirs, key=_date_proximity(scene_date)).strftime(pattern)


def get_brdf_dirs_pre_modis(brdf_root, scene_date):
    """
    Get list of pre-MODIS BRDF directories for the dataset.

    :param brdf_root:
        BRDF root directory.
    :type brdf_root:
        :py:class:`str`

    :param scene_date:
        Scene Date.
    :type scene_date:
        :py:class:`datetime.date`

    :return:
       A string containing the closest matching BRDF directory.

    """
    # Find the N (=n_dirs) BRDF directories with midpoints closest to the
    # scene date.
    # Pre-MODIS BRDF directories are named 'XXX' (day-of-year).
    # Return a list of n_dirs directories to maintain compatibility with
    # the NBAR code, even though we know that the nearest day-of-year
    # database dir will contain usable data.
    # Build list of dates for comparison
    dir_dates = []

    # Standardise names be prepended with leading zeros
    for doy in sorted(os.listdir(brdf_root), key=lambda x: x.zfill(3)):
        dir_dates.append((str(scene_date.year), doy))

    # Add boundary entry for previous year
    dir_dates.insert(0, (str(scene_date.year - 1), dir_dates[-1][1]))
    # Add boundary entry for next year accounting for inserted entry
    dir_dates.append((str(scene_date.year + 1), dir_dates[1][1]))

    # Interpreter function
    doy_intpr = lambda x: datetime.datetime.strptime(' '.join(x), '%Y %j').date()

    # return directory name without year
    return min(dir_dates, key=_date_proximity(scene_date, doy_intpr))[1]


def coord_transformer(src_crs, dst_crs):
    """
    Coordinate transformation function between CRSs.

    :param src_crs:
        Source CRS.
    :type src_crs:
        :py:class:`rasterio.crs.CRS`

    :param dst_crs:
        Destination CRS.
    :type dst_crs:
        :py:class:`rasterio.crs.CRS`

    :return:
        A function that takes a point in the source CRS and returns the same
        point expressed in the destination CRS.
    """
    def crs_to_proj(crs):
        return pyproj.Proj(**crs.to_dict())

    def result(*args, **kwargs):
        return pyproj.transform(crs_to_proj(src_crs), crs_to_proj(dst_crs), *args, **kwargs)

    return result


class BrdfTileSummary:
    """
    A lightweight class to represent the BRDF information gathered from a tile.
    """
    def __init__(self, valid_pixels, brdf_sums, source_files):
        self.valid_pixels = valid_pixels
        self.brdf_sums = brdf_sums
        self.source_files = source_files

    @staticmethod
    def empty():
        """ When the tile is not inside the ROI. """
        return BrdfTileSummary(0, {key: 0.0 for key in BrdfModelParameters}, [])

    def __add__(self, other):
        """ Accumulate information from different tiles. """
        return BrdfTileSummary(self.valid_pixels + other.valid_pixels,
                               {key: self.brdf_sums[key] + other.brdf_sums[key] for key in BrdfModelParameters},
                               self.source_files + other.source_files)

    def mean(self):
        """ Calculate the mean BRDF parameters. """
        if self.valid_pixels == 0:
            raise BRDFLookupError("no brdf datasets found for the ROI")

        # ratio of spatial averages
        averages = {key: self.brdf_sums[key] / self.valid_pixels
                    for key in BrdfModelParameters}

        bands = {BrdfDirectionalParameters.ALPHA_1: BrdfModelParameters.VOL,
                 BrdfDirectionalParameters.ALPHA_2: BrdfModelParameters.GEO}

        return {key: dict(dataset_ids=self.source_files,
                          value=averages[bands[key]] / averages[BrdfModelParameters.ISO])
                for key in BrdfDirectionalParameters}


def valid_region(fname, mask_value=None):
    """
    Return valid data region for input images based on mask value and input image path
    """
    log.info("Valid regions for %s", fname)

    # ensure formats match
    with rasterio.open(str(fname), 'r') as dataset:
        transform = dataset.transform.to_gdal()
        crs = dataset.crs.to_dict()
        img = dataset.read(1)

        if mask_value is not None:
            mask = img & mask_value == mask_value
        else:
            mask = img != 0

    shapes = rasterio.features.shapes(mask.astype('uint8'), mask=mask)
    shape = ops.unary_union([shapely.geometry.shape(shape) for shape, val in shapes if val == 1])

    # convex hull
    geom = shape.convex_hull

    # buffer by 1 pixel
    geom = geom.buffer(1, join_style=3, cap_style=3)

    # simplify with 1 pixel radius
    geom = geom.simplify(1)

    # intersect with image bounding box
    geom = geom.intersection(shapely.geometry.box(0, 0, mask.shape[1], mask.shape[0]))

    # transform from pixel space into CRS space
    geom = shapely.affinity.affine_transform(
        geom, (transform[1], transform[2], transform[4],
               transform[5], transform[0], transform[3])
    )

    return geom, crs


def load_brdf_tile(src_poly, src_crs, fid, dataset_name, is_fallback, fid_mask):
    """
    Summarize BRDF data from a single tile.
    """
    ds = fid[dataset_name]

    def segmentize_src_poly(length_scale):
        src_poly_geom = ogr.CreateGeometryFromWkt(src_poly.wkt)
        src_poly_geom.Segmentize(length_scale)
        return wkt.loads(src_poly_geom.ExportToWkt())

    _, ds_width, ds_height = ds.shape

    dst_geotransform = rasterio.transform.Affine.from_gdal(*ds.attrs['geotransform'])
    dst_crs = CRS.from_wkt(ds.attrs['crs_wkt'])

    # get bounds of BRDF modis tile from h5 dataset 
    left, bottom, right, top = rasterio.transform.array_bounds(ds_height, ds_width, dst_geotransform)
    
    # get a tile window to read from continental coastal mask 
    window = rasterio.windows.from_bounds(left, bottom, right, top, transform=fid_mask.transform)
    
    # read ocean mask file for correspoing tile window
    ocean_data = fid_mask.read(1, window=window)
    # assumes the length scales are the same (m)
    dst_poly = ops.transform(coord_transformer(src_crs, dst_crs),
                             segmentize_src_poly(np.sqrt(np.abs(dst_geotransform.determinant))))

    bound_poly = ops.transform(lambda x, y: dst_geotransform * (x, y), box(0., 0., ds_width, ds_height))
    if not bound_poly.intersects(dst_poly):
        return BrdfTileSummary.empty()

    mask = rasterize([(dst_poly, 1)], fill=0, out_shape=(ds_width, ds_height), transform=dst_geotransform)
    assert ocean_mask.shape == mask.shape    # shape should be same 
    
    # if both ocean_data (1=land, 0=ocean) and dst_poly (1=valid, 0=invalid)
    # then multiplying two masks would result in 1 being if both masks are true
    mask = mask * ocean_data
    
    valid_pixels = np.sum(mask)
    mask = mask.astype(bool)
    assert mask.shape == ds.shape[:-2]

    def layer_sum(i):
        #TODO Discuss with Imam: should we compute individual number of valid pixels for each BRDF parameter? 

        layer = ds[i, :, :]
        fill_value_mask = (layer != ds.attrs['_FillValue'])
        layer = layer.astype('float32')
        layer[~mask] = np.nan
        layer[~fill_value_mask] = np.nan
        layer = ds.attrs['scale_factor'] * (layer - ds.attrs['add_offset'])
        return np.nansum(layer)

    return BrdfTileSummary(valid_pixels,
                           {key: layer_sum(index) for index, key in enumerate(BrdfModelParameters)},
                           [current_h5_metadata(fid)['id']])


def get_brdf_data(acquisition, brdf,
                  compression=H5CompressionFilter.LZF, filter_opts=None):
    """
    Calculates the mean BRDF value for the given acquisition,
    for each BRDF parameter ['geo', 'iso', 'vol'] that covers
    the acquisition's extents.

    :param acquisition:
        An instance of an acquisitions object.

    :param brdf:
        A `dict` defined as either of the following:

        * {'user': {<band-alias>: {'iso': <value>, 'vol': <value>, 'geo': <value>}, ...}}
        * {'brdf_path': <path-to-BRDF>, 'brdf_premodis_path': <path-to-average-BRDF>, 'ocean_mask_path': <path-to-ocean-mask>}
        
        Here <path-to-BRDF> is a string containing the full file system
        path to your directory containing the ource BRDF files
        The BRDF directories are assumed to be yyyy.mm.dd naming convention.

        And <path-to-average-BRDF> is a string containing the full file system
        path to your directory containing the Jupp-Li backup BRDF data.
        To be used for pre-MODIS and potentially post-MODIS acquisitions.

        And <path-to-ocean-mask> is a string containing the full file system path
        to your ocean mask file. To be used for masking ocean pixels from  BRDF data 
        from pre-MODIS and post-MODIS acquisitions. 
  
    :param compression:
        The compression filter to use.
        Default is H5CompressionFilter.LZF

    :filter_opts:
        A dict of key value pairs available to the given configuration
        instance of H5CompressionFilter. For example
        H5CompressionFilter.LZF has the keywords *chunks* and *shuffle*
        available.
        Default is None, which will use the default settings for the
        chosen H5CompressionFilter instance.

    :return:
        A `dict` with the keys:

            * BrdfDirectionalParameters.ALPHA_1
            * BrdfDirectionalParameters.ALPHA_2

        Values for each BRDF Parameter are accessed via the key named
        `value`.

    :notes:
        The keywords compression and filter_opts aren't used as we no
        longer save the BRDF imagery. However, we may need to store
        tables in future, therefore they can remain until we know
        for sure they'll never be used.
    """

    if 'user' in brdf:
        # user-specified override
        return {param: dict(data_source='BRDF', model='user',
                            value=brdf['user'][acquisition.alias][param.value.lower()])
                for param in BrdfDirectionalParameters}

    brdf_primary_path = brdf['brdf_path']
    brdf_secondary_path = brdf['brdf_premodis_path']
    brdf_ocean_mask_path = brdf['ocean_mask_path']
    
    # Get the date of acquisition
    dt = acquisition.acquisition_datetime.date()

    # Compare the scene date and MODIS BRDF start date to select the
    # BRDF data root directory.
    # Scene dates outside the range of the CSIRO mosaic data
    # should use the pre-MODIS, Jupp-Li BRDF.
    brdf_dir_list = sorted(os.listdir(brdf_primary_path))

    try:
        brdf_dir_range = [brdf_dir_list[0], brdf_dir_list[-1]]
        brdf_range = [datetime.date(*[int(x) for x in y.split('.')])
                      for y in brdf_dir_range]

        fallback_brdf = (dt < brdf_range[0] or dt > brdf_range[1])
    except IndexError:
        fallback_brdf = True  # use JuppLi if no primary data available

    if fallback_brdf:
        brdf_base_dir = brdf_secondary_path
        brdf_dirs = get_brdf_dirs_pre_modis(brdf_base_dir, dt)
    else:
        brdf_base_dir = brdf_primary_path
        brdf_dirs = get_brdf_dirs_modis(brdf_base_dir, dt)

    # get all HDF files in the input dir
    dbDir = pjoin(brdf_base_dir, brdf_dirs)
    tile_list = [pjoin(folder, f)
                 for (folder, _, filelist) in os.walk(dbDir) for f in filelist if f.endswith(".h5")]

    src_poly, src_crs = valid_region(acquisition.uri, acquisition.no_data)
    src_crs = rasterio.crs.CRS(**src_crs)

    tally = {}
    with rasterio.open(brdf_ocean_mask_path, 'r') as fid_mask:
        for ds in acquisition.brdf_datasets:
            tally[ds] = BrdfTileSummary.empty(fallback_brdf)
            for tile in tile_list:
                with h5py.File(tile, 'r') as fid:
                    tally[ds] += load_brdf_tile(src_poly, src_crs, fid, ds, fallback_brdf, fid_mask)
            tally[ds] = tally[ds].mean()

    results = {param: dict(data_source='BRDF',
                           dataset_ids=",".join({ds_id
                                                 for ds in acquisition.brdf_datasets
                                                 for ds_id in tally[ds][param]['dataset_ids']}),
                           value=np.mean([tally[ds][param]['value'] for ds in acquisition.brdf_datasets]).item(),
                           model='fallback' if fallback_brdf else 'definitive')
               for param in BrdfDirectionalParameters}

    return results
