import os
import numpy
import logging
from osgeo import gdal
from gaip import GriddedGeoBox

def calc_land_sea_mask(gridded_geo_box, utm_zone, ancillary_path):
    """
    Creates a Land/Sea mask.

    :param gridded_geo_box:
        An instance of GriddedGeoBox prescribing the spatial properties
        of the pixel quality mask.

    :param utm_zone:
        The UTM_ZONE for this mask.

    :param ancillary_path:
        Path to land/sea ancillary data directory

    :return:
        A 2D Numpy Boolean array. True = Land, False = Sea.

    :author:
        Josh Sixsmith, joshua.sixsmith@ga.gov.au
    """

    rasfile = os.path.join(ancillary_path, 'WORLDzone%02d.tif' % abs(utm_zone))
    assert os.path.exists(rasfile), 'ERROR: Raster File Not Found (%s)' % rasfile

    # open the land/sea data file

    lsobj   = gdal.Open(rasfile, gdal.gdalconst.GA_ReadOnly)
    land_sea_geo_box = GriddedGeoBox.from_gdal_dataset(lsobj)
    print "shape=%s" % (str(gridded_geo_box.shape), )

    affine = gridded_geo_box.affine
    mUL = (0, 0) * affine
    mLR = gridded_geo_box.shape * affine

    print "mUL=%s" % (str(mUL), )
    print "mLR=%s" % (str(mLR), )

    # Convert the map co-ords into the rasfile image co-ords

    affine = ~ land_sea_geo_box.affine
    iUL = mUL * affine
    iLR = mLR * affine
    print "iUL=%s" % (str(iUL), )
    print "iLR=%s" % (str(iLR), )

    xoff = int(iUL[1])
    yoff = int(iUL[0])
    xsize = int(iLR[1] - xoff)
    ysize = int(iLR[0] - yoff)
    print "xoff=%s" % (str(xoff), )
    print "yoff=%s" % (str(yoff), )
    print "xsize=%s" % (str(xsize), )
    print "ysize=%s" % (str(ysize), )

    # Read in the land/sea array
    ls_arr = lsobj.ReadAsArray(xoff, yoff, xsize, ysize)
    lsobj.FlushCache()
    lsobj = None

    return (ls_arr.astype('bool'))

def setLandSeaBit(gridded_geo_box, utm_zone, ancillary_path, pq_const, pqaResult):
    mask = calc_land_sea_mask(gridded_geo_box, utm_zone, ancillary_path)
    bit_index = pq_const.land_sea
    pqaResult.set_mask(mask, bit_index)
