#!/bin/env python 
#
# get Land/Sea mask from UTM dataset
#-----------------------------------
import rasterio as rio
from rasterio.crs import from_string
import osr
from gaip import GriddedGeoBox
import numpy
from affine import Affine
from rasterio.warp import reproject, RESAMPLING

def getUtmZone(pos_longlat):
    """
    Return the UTM zone number corresponding to the supplied position

    Arguments:
        pos_longlat: position as tuple (This in (long, lat)

    Returns:
        The UTM zone number (+ve for North, -ve for South)
    """

    lon, lat = pos_longlat
    z = int(lon/6) + 31
    if lat < 0:
        z = -z
    return z

def get_land_sea_mask(gridded_geo_box, \
        ancillary_path='/g/data/v10/eoancillarydata/Land_Sea_Rasters'):
    """
    Return a land/sea 2D numpy boolean array in which Land = True, Sea = False
    for the supplied GriddedGeoBox and using the UTM projected data in the 
    supplied ancillary_path.

    If the specified gridded_geo_box has a non-UTM CRS or a non-native
    sample frequency, the data will be reprojected into the the 
    gridded_geo_box grid.
    """

    # get lat/long of geo_box origin

    to_crs = osr.SpatialReference()
    to_crs.SetFromUserInput('EPSG:4326')
    origin_longlat = gridded_geo_box.transform_coordinates(gridded_geo_box.origin, to_crs)
    print "origin_lonlat=%s" % (origin_longlat,)

    # get Land/Sea data file for this bounding box
    utmZone = abs(getUtmZone(origin_longlat))
    utmDataPath = '%s/WORLDzone%d.tif' % (ancillary_path, utmZone)

    # read the data for the Flinders islet region
    with rio.open(utmDataPath) as ds:

        # get the gridded box for the full data extent
        landSeaDataGGB = GriddedGeoBox.from_dataset(ds)
        print "land/sea geo_box=%s" % (str(landSeaDataGGB))
        print "land/ses affine=\n%s" % (str(landSeaDataGGB.affine))

        # read the subset relating to Flinders Islet
        window = landSeaDataGGB.window(gridded_geo_box)
        print "window=%s" % (str(window))
        mask = ds.read(1, window=window)

        # check if we need to resample
        if mask.shape == gridded_geo_box.shape:
            return mask  # no reprojection

        # reprojection required
        print "result shape is %s, expecting shape %s, reprojecting!!!!!!!!!!" \
            % (mask.shape, gridded_geo_box.shape)
        reprojected_mask = numpy.zeros(gridded_geo_box.shape, mask.dtype)
        a = landSeaDataGGB.affine

        (x, y) = ~a * window[0]
        src_affine = Affine(a.a, a.b, x, a.d, a.e, y)
        print "src_geo_box=\n%s" % (str(landSeaDataGGB))
        print "src_affine=\n%s" % (str(src_affine))
        print "dst_geo_box=\n%s" % (str(gridded_geo_box))
        print "dst_affine=\n%s" % (str(gridded_geo_box.affine))
        print "mask="
        print mask
        reproject(
            mask,
            reprojected_mask,
            src_transform=src_affine,
            src_crs=from_string(landSeaDataGGB.crs.ExportToProj4()),
            dst_transform=gridded_geo_box.affine,
            dst_crs=from_string(gridded_geo_box.crs.ExportToProj4()),
            resampling=RESAMPLING.nearest)

        return reprojected_mask


if __name__ == "__main__":
    # need this hack until set_epsilon patch v1.0.5
    # affine.set_epsilon(1e-9)
    import affine
    affine.EPSILON=1e-9
    affine.EPSILON2=1e-18

    # choose Flinders Islet bounding box
    flindersOrigin = (150.927659, -34.453309)
    flindersCorner = (150.931697, -34.457915)

    flindersGGB = GriddedGeoBox.from_corners(flindersOrigin, flindersCorner)
    mask = get_land_sea_mask(flindersGGB)
    print mask
