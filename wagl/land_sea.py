#!/bin/env python
#
# get Land/Sea mask from UTM dataset
# -----------------------------------
from __future__ import absolute_import, print_function
import rasterio as rio
from osgeo import osr
import numpy

from wagl.geobox import GriddedGeoBox
from wagl.data import write_img

# pylint: disable=invalid-name


def get_utm_zone(pos_longlat):
    """
    Return the UTM zone number corresponding to the supplied position

    Arguments:
        pos_longlat: position as tuple (This in (long, lat)

    Returns:
        The UTM zone number (+ve for North, -ve for South)
    """

    lon, lat = pos_longlat
    z = int(lon / 6) + 31
    if lat < 0:
        z = -z
    return z


def get_land_sea_mask(
    gridded_geo_box, ancillary_path="/g/data/v10/eoancillarydata/Land_Sea_Rasters"
):
    """
    Return a land/sea 2D numpy boolean array in which Land = True, Sea = False
    for the supplied GriddedGeoBox and using the UTM projected data in the
    supplied ancillary_path.

    If the specified gridded_geo_box has a non-UTM CRS or a non-native
    sample frequency, the data will be reprojected/resampled into the the
    gridded_geo_box.
    """

    # get lat/long of geo_box origin

    to_crs = osr.SpatialReference()
    to_crs.SetFromUserInput("EPSG:4326")
    origin_longlat = gridded_geo_box.transform_coordinates(gridded_geo_box.origin, to_crs)

    # get Land/Sea data file for this bounding box
    utmZone = abs(get_utm_zone(origin_longlat))
    utmDataPath = "%s/WORLDzone%d.tif" % (ancillary_path, utmZone)

    # read the land/sea data
    with rio.open(utmDataPath) as ds:
        # get the gridded box for the full dataset extent
        landSeaDataGGB = GriddedGeoBox.from_dataset(ds)

        # read the subset relating to Flinders Islet
        window = landSeaDataGGB.window(gridded_geo_box)
        out = numpy.zeros(gridded_geo_box.shape, dtype=numpy.uint8)
        ds.read(1, window=window, out=out)

        return out


if __name__ == "__main__":
    # need this hack until set_epsilon patch v1.0.5
    # affine.set_epsilon(1e-9)
    import affine

    affine.EPSILON = 1e-9
    affine.EPSILON2 = 1e-18

    # choose Flinders Islet bounding box
    flindersOrigin = (150.927659, -34.453309)
    flindersCorner = (150.931697, -34.457915)

    flindersGGB = GriddedGeoBox.from_corners(flindersOrigin, flindersCorner)
    mask = get_land_sea_mask(flindersGGB)
    print("geobox_shape=%s" % str(flindersGGB.shape))
    print("mask_shape=%s" % str(mask.shape))
    print(mask)

    # same test for AGDC cell around  Darwin area

    scale = 0.00025
    shape = (4000, 4000)
    origin = (130.0, -12.0)

    corner = (shape[1] * scale + origin[0], origin[1] - shape[0] * scale)
    ggb = GriddedGeoBox(shape, origin, pixelsize=(scale, scale))
    print(ggb)

    # now get UTM equilavent

    geo_box = ggb.copy(crs="EPSG:32752")
    print(geo_box)

    # and get the mask

    mask = get_land_sea_mask(geo_box)

    total_pixels = geo_box.shape[1] * geo_box.shape[0]
    land_pixels = sum(sum(mask.astype("uint32")))
    sea_pixels = total_pixels - land_pixels
    sea_pct = 100.0 * sea_pixels / total_pixels
    land_pct = 100.0 * land_pixels / total_pixels

    print("ggb_shape=%s" % str(ggb.shape))
    print("geobox_shape=%s" % str(geo_box.shape))
    print("mask_shape=%s" % str(mask.shape))
    print("total_pixels=%d" % total_pixels)
    print("land_pixels=%d" % land_pixels)
    print("sea_pixels=%d" % sea_pixels)
    # self.assertEqual(land_pixels, 14554858)
    # self.assertEqual(sea_pixels, 1445142)
    # self.assertEqual(total_pixels, 16000000)

    print("land=%f%%, sea=%f%%" % (land_pct, sea_pct))
    write_img(mask, "mask.tif", driver="GTiff", geobox=ggb)
