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
from skimage.transform import downscale_local_mean

from wagl.constants import BrdfDirectionalParameters, BrdfModelParameters, BrdfTier
from wagl.hdf5 import H5CompressionFilter, VLEN_STRING
from wagl.metadata import current_h5_metadata
from wagl.data import read_subset


_LOG = logging.getLogger(__name__)


# Accurate BRDF requires both Terra and Aqua to be operating
# Aqua launched 2002-05-04, so we'll add a buffer for determining the start
# date for using definitive data.
DEFINITIVE_START_DATE = datetime.datetime(2002, 7, 1).date()


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
            -1 * _date.day,
        )

    return _proximity_comparator


def get_brdf_dirs_viirs(brdf_root, scene_date, pattern="%Y.%m.%d"):
    # our VIIRS collection follows the same folder structure as our MODIS collection
    return get_brdf_dirs_modis(brdf_root, scene_date, pattern=pattern)


def get_brdf_dirs_modis(brdf_root, scene_date, pattern="%Y.%m.%d"):
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

    if dirs == []:
        raise IndexError

    return min(dirs, key=_date_proximity(scene_date)).strftime(pattern)


def get_brdf_dirs_fallback(brdf_root, scene_date):
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

    if dir_dates == []:
        raise IndexError

    # return directory name without year
    return min(
        dir_dates,
        key=_date_proximity(
            scene_date,
            lambda x: datetime.datetime.strptime(" ".join(x), "%Y %j").date(),
        ),
    )[1]


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
        return pyproj.transform(
            crs_to_proj(src_crs), crs_to_proj(dst_crs), *args, **kwargs
        )

    return result


class BrdfTileSummary:
    """
    A lightweight class to represent the BRDF information gathered from a tile.
    """

    def __init__(self, brdf_summaries, source_files):
        self.brdf_summaries = brdf_summaries
        self.source_files = source_files

    @staticmethod
    def empty():
        """When the tile is not inside the ROI."""
        return BrdfTileSummary(
            {key: {"sum": 0.0, "count": 0} for key in BrdfModelParameters}, []
        )

    def is_empty(self):
        return all(self.brdf_summaries[key]["count"] == 0 for key in BrdfModelParameters)

    def __add__(self, other):
        """Accumulate information from different tiles."""

        def add(key):
            this = self.brdf_summaries[key]
            that = other.brdf_summaries[key]
            return {
                "sum": this["sum"] + that["sum"],
                "count": this["count"] + that["count"],
            }

        return BrdfTileSummary(
            {key: add(key) for key in BrdfModelParameters},
            self.source_files + other.source_files,
        )

    def mean(self):
        """Calculate the mean BRDF parameters."""
        if self.is_empty():
            # possibly over the ocean, so lambertian
            return {
                key: dict(id=self.source_files, value=0.0)
                for key in BrdfDirectionalParameters
            }

        # ratio of spatial averages
        averages = {
            key: self.brdf_summaries[key]["sum"] / self.brdf_summaries[key]["count"]
            for key in BrdfModelParameters
        }

        bands = {
            BrdfDirectionalParameters.ALPHA_1: BrdfModelParameters.VOL,
            BrdfDirectionalParameters.ALPHA_2: BrdfModelParameters.GEO,
        }

        return {
            key: dict(
                id=self.source_files,
                value=averages[bands[key]] / averages[BrdfModelParameters.ISO],
            )
            for key in BrdfDirectionalParameters
        }


def valid_region(acquisition, mask_value=None):
    """
    Return valid data region for input images based on mask value and input image path
    """
    img = acquisition.data()
    gbox = acquisition.gridded_geo_box()
    crs = CRS.from_wkt(gbox.crs.ExportToWkt()).to_dict()
    transform = gbox.transform.to_gdal()

    if mask_value is None:
        mask_value = acquisition.no_data

    if mask_value is not None:
        mask = img != mask_value
    else:
        mask = img != 0

    shapes = rasterio.features.shapes(mask.astype("uint8"), mask=mask)
    shape = ops.unary_union(
        [shapely.geometry.shape(shape) for shape, val in shapes if val == 1]
    )

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
        geom,
        (
            transform[1],
            transform[2],
            transform[4],
            transform[5],
            transform[0],
            transform[3],
        ),
    )

    return geom, crs


def extract_VIIRS_geotransform(f):
    """Takes in a VIIRS hdf5 file and returns its geotransform matrix"""

    fileMetadata = f["HDFEOS INFORMATION"]["StructMetadata.0"][
        ()
    ].split()  # Split at newline
    fileMetadata = [m.decode("utf-8") for m in fileMetadata]  # Convert bytes to strings

    for md in fileMetadata:
        if "UpperLeftPointMtrs=(" in md:
            mtrs = md.split("=", 1)[1]
            ulc = np.fromstring(mtrs[1:-1], dtype=float, sep=",")
        elif "LowerRightMtrs=(" in md:
            mtrs = md.split("=", 1)[1]
            lrc = np.fromstring(mtrs[1:-1], dtype=float, sep=",")
        elif "XDim=" in md:
            mtrs = md.split("=", 1)[1]
            width = float(mtrs)
        elif "YDim=" in md:
            mtrs = md.split("=", 1)[1]
            height = float(mtrs)

    xres, yres = np.divide(np.subtract(lrc, ulc), (width, height))

    geoInfo = (ulc[0], xres, 0, ulc[1], 0, yres)
    return geoInfo


def VIIRS_crs():
    """Return VIIRS projection - exact same as MODIS"""
    prj = 'PROJCS["unnamed",\
    GEOGCS["Unknown datum based upon the custom spheroid", \
    DATUM["Not specified (based on custom spheroid)", \
    SPHEROID["Custom spheroid",6371007.181,0]], \
    PRIMEM["Greenwich",0],\
    UNIT["degree",0.0174532925199433]],\
    PROJECTION["Sinusoidal"], \
    PARAMETER["longitude_of_center",0], \
    PARAMETER["false_easting",0], \
    PARAMETER["false_northing",0], \
    UNIT["Meter",1]]'
    return prj


def load_brdf_tile(src_poly, src_crs, fid, dataset_name, fid_mask, satellite):
    """
    Summarize BRDF data from a single tile.
    """
    assert satellite in ["MODIS", "VIIRS"]

    if satellite == "MODIS":
        ds = fid[dataset_name]
        ds_height, ds_width = ds.shape
        dst_geotransform = rasterio.transform.Affine.from_gdal(*ds.attrs["geotransform"])
        dst_crs = CRS.from_wkt(ds.attrs["crs_wkt"])

    else:
        ds = fid["HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/"][dataset_name]
        ds_height, ds_width = ds.shape[:2]
        gt = extract_VIIRS_geotransform(fid)
        dst_geotransform = rasterio.transform.Affine.from_gdal(*gt)
        dst_crs = CRS.from_wkt(VIIRS_crs())

    def segmentize_src_poly(length_scale):
        src_poly_geom = ogr.CreateGeometryFromWkt(src_poly.wkt)
        src_poly_geom.Segmentize(length_scale)
        return wkt.loads(src_poly_geom.ExportToWkt())

    # assumes the length scales are the same (m)
    dst_poly = ops.transform(
        coord_transformer(src_crs, dst_crs),
        segmentize_src_poly(np.sqrt(np.abs(dst_geotransform.determinant))),
    )

    bound_poly = ops.transform(
        lambda x, y: dst_geotransform * (x, y),
        box(0.0, 0.0, ds_width, ds_height, ccw=False),
    )
    if not bound_poly.intersects(dst_poly):
        return BrdfTileSummary.empty()

    ocean_poly = ops.transform(
        lambda x, y: fid_mask.transform * (x, y),
        box(0.0, 0.0, fid_mask.width, fid_mask.height),
    )
    if not ocean_poly.intersects(dst_poly):
        return BrdfTileSummary.empty()

    # read ocean mask file for correspoing tile window
    # land=1, ocean=0
    bound_poly_coords = list(bound_poly.exterior.coords)[:4]
    ocean_mask, _ = read_subset(fid_mask, *bound_poly_coords)
    ocean_mask = ocean_mask.astype(bool)

    # inside=1, outside=0
    roi_mask = rasterize(
        [(dst_poly, 1)],
        fill=0,
        out_shape=(ds_height, ds_width),
        transform=dst_geotransform,
    )
    roi_mask = roi_mask.astype(bool)

    if roi_mask.shape != ocean_mask.shape:
        assert len(roi_mask.shape) == 2 and len(ocean_mask.shape) == 2
        x_ratio = ocean_mask.shape[0] / roi_mask.shape[0]
        assert int(x_ratio) == x_ratio
        y_ratio = ocean_mask.shape[1] / roi_mask.shape[1]
        assert int(y_ratio) == y_ratio
        ocean_mask = downscale_local_mean(ocean_mask, (int(x_ratio), int(y_ratio))) == 1

    # both ocean_mask and mask shape should be same
    if ocean_mask.shape != roi_mask.shape:
        raise ValueError("ocean mask and ROI mask do not have the same shape")
    if roi_mask.shape != (ds_height, ds_width):
        raise ValueError("BRDF dataset and ROI mask do not have the same shape")

    roi_mask = roi_mask & ocean_mask

    def layer_sum(param):
        if satellite == "MODIS":
            layer = ds[param][:, :]
        else:
            map_dict = {
                BrdfModelParameters.ISO: 0,
                BrdfModelParameters.VOL: 1,
                BrdfModelParameters.GEO: 2,
            }
            layer = ds[:, :, map_dict[param]]
        common_mask = roi_mask & (layer != ds.attrs["_FillValue"])
        layer = layer.astype("float32")
        layer[~common_mask] = np.nan
        layer = ds.attrs["scale_factor"] * (layer - ds.attrs["add_offset"])
        return {"sum": np.nansum(layer), "count": np.sum(common_mask)}

    if satellite == "MODIS":
        bts = BrdfTileSummary(
            {param: layer_sum(param.value) for param in BrdfModelParameters},
            [current_h5_metadata(fid)["id"]],
        )
    else:
        bts = BrdfTileSummary(
            {param: layer_sum(param) for param in BrdfModelParameters},
            [fid.attrs["LocalGranuleID"].decode("UTF-8")],
        )

    return bts


def get_tally(mode, brdf_config, brdf_datasets, viirs_datasets, dt, src_poly, src_crs):
    """
    Get all HDF files in the input dir.
    `mode` can be one of `MODIS`, `VIIRS` or `fallback`

    Raises `IndexError` if it can't find the required data.
    """

    def assert_exists(path):
        if not os.path.isdir(path):
            raise IndexError
        return path

    if mode == "fallback":
        brdf_base_dir = assert_exists(brdf_config["brdf_fallback_path"])
        brdf_dirs = get_brdf_dirs_fallback(brdf_base_dir, dt)
        satellite = "MODIS"
        datasets = brdf_datasets

    elif mode == "MODIS":
        brdf_base_dir = assert_exists(brdf_config["brdf_path"])
        brdf_dirs = get_brdf_dirs_modis(brdf_base_dir, dt)
        satellite = "MODIS"
        datasets = brdf_datasets

        # Compare the scene date and MODIS BRDF start date to select the
        # BRDF data root directory.
        # Scene dates outside this range are to use the fallback data
        brdf_dir_list = sorted(os.listdir(brdf_base_dir))
        brdf_dir_range = [brdf_dir_list[0], brdf_dir_list[-1]]
        brdf_range = [
            datetime.date(*[int(x) for x in y.split(".")]) for y in brdf_dir_range
        ]
        if dt < DEFINITIVE_START_DATE or dt > brdf_range[1]:
            raise IndexError

    elif mode == "VIIRS":
        satellite = "VIIRS"
        datasets = list(viirs_datasets.values())[0]

        assert not ("I" in viirs_datasets and "M" in viirs_datasets)
        if "I" in viirs_datasets.keys():
            brdf_base_dir = assert_exists(brdf_config["viirs_i_path"])
            brdf_dirs = get_brdf_dirs_viirs(brdf_base_dir, dt)
        elif "M" in viirs_datasets.keys():
            brdf_base_dir = assert_exists(brdf_config["viirs_m_path"])
            brdf_dirs = get_brdf_dirs_viirs(brdf_base_dir, dt)
        else:
            raise ValueError("No I or M bands in VIIRS band for sensor")

    else:
        raise ValueError(f"Unknown mode {mode}")

    dbDir = pjoin(brdf_base_dir, brdf_dirs)
    tile_list = [
        pjoin(folder, f)
        for (folder, _, filelist) in os.walk(dbDir)
        for f in filelist
        if f.endswith(".h5")
    ]

    tally = {}
    with rasterio.open(brdf_config["ocean_mask_path"], "r") as fid_mask:
        for ds in datasets:
            tally[ds] = BrdfTileSummary.empty()
            for tile in tile_list:
                with h5py.File(tile, "r") as fid:
                    tally[ds] += load_brdf_tile(
                        src_poly, src_crs, fid, ds, fid_mask, satellite
                    )
    return tally


def get_brdf_data(
    acquisition,
    brdf_config,
    mode=None,
    compression=H5CompressionFilter.LZF,
    filter_opts=None,
):
    """
    Calculates the mean BRDF value for the given acquisition,
    for each BRDF parameter ['geo', 'iso', 'vol'] that covers
    the acquisition's extents.

    :param acquisition:
        An instance of an acquisitions object.

    :param brdf_config:
        A `dict` defined as either of the following:
        * {'user': {<band-alias>: {'iso': <value>, 'vol': <value>, 'geo': <value>}, ...}}
        * {'brdf_path': <path-to-BRDF>, 'brdf_fallback_path': <path-to-average-BRDF>,
           'ocean_mask_path': <path-to-ocean-mask>}

        Here <path-to-BRDF> is a string containing the full file system
        path to your directory containing the ource BRDF files
        The BRDF directories are assumed to be yyyy.mm.dd naming convention.

        <path-to-average-BRDF> is a string containing the full file system
        path to your directory containing the fallback BRDF data.
        To be used for pre-MODIS and potentially post-MODIS acquisitions.

        And <path-to-ocean-mask> is a string containing the full file system path
        to your ocean mask file. To be used for masking ocean pixels from  BRDF data
        all acquisitions.

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

    if "user" in brdf_config:
        # user-specified override
        return {
            param: dict(
                data_source="BRDF",
                tier=BrdfTier.USER.name,
                value=brdf_config["user"][acquisition.alias][param.value.lower()],
            )
            for param in BrdfDirectionalParameters
        }

    src_poly, src_crs = valid_region(acquisition)
    src_crs = rasterio.crs.CRS(**src_crs)

    # Get the date of acquisition
    dt = acquisition.acquisition_datetime.date()

    brdf_datasets = acquisition.brdf_datasets
    viirs_datasets = {"I": ["BRDF_Albedo_Parameters_I1"]}

    def get_tally2(mode, dt):
        # brdf_config, brdf_datasets, and viirs datasets are "constants"
        # for the purpose of choosing the data to use (MODIS vs VIIRS vs fallback)
        result = get_tally(
            mode, brdf_config, brdf_datasets, viirs_datasets, dt, src_poly, src_crs
        )

        if any(result[ds].is_empty() for ds in result):
            raise IndexError
        return result

    def back_in_time(from_dt):
        days_back = 1

        while days_back <= 30:
            dt = from_dt - datetime.timedelta(days=days_back)

            try:
                return get_tally2("MODIS", dt)
            except IndexError:
                try:
                    return get_tally2("VIIRS", dt)
                except IndexError:
                    pass

            days_back += 1

        # if we are here, we failed, nothing was found
        raise IndexError

    def select_mode(dt):
        # try MODIS first, then VIIRS, then 30 days back in time, then fallback
        try:
            return get_tally2("MODIS", dt), False
        except IndexError:
            pass

        try:
            return get_tally2("VIIRS", dt), False
        except IndexError:
            pass

        try:
            return back_in_time(dt), False
        except IndexError:
            pass

        try:
            return get_tally2("fallback", dt), True
        except IndexError:
            pass

        raise ValueError(f"No BRDF ancillary found for {dt}")

    if mode is None:
        tally, fallback_brdf = select_mode(dt)
    else:
        tally = get_tally2(mode, dt)
        fallback_brdf = True if mode == "fallback" else False

    spatial_averages = {}
    for ds in tally:
        spatial_averages[ds] = tally[ds].mean()

    results = {
        param: dict(
            data_source="BRDF",
            id=np.array(
                list(
                    {
                        ds_id
                        for ds in spatial_averages
                        for ds_id in spatial_averages[ds][param]["id"]
                    }
                ),
                dtype=VLEN_STRING,
            ),
            value=np.mean(
                [spatial_averages[ds][param]["value"] for ds in spatial_averages]
            ).item(),
            tier=BrdfTier.FALLBACK_DATASET.name
            if fallback_brdf
            else BrdfTier.DEFINITIVE.name,
        )
        for param in BrdfDirectionalParameters
    }

    return results
