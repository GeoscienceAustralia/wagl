from __future__ import absolute_import, print_function
import os
import numpy

from osgeo import gdal
from wagl.metadata import extract_ancillary_metadata


def calc_land_sea_mask(
    geo_box, ancillary_path="/g/data/v10/eoancillarydata/Land_Sea_Rasters"
):
    """
    Creates a Land/Sea mask.

    :param geo_box:
        An instance of GriddedGeoBox defining the region for which the
        land/sea mask is required.

            * WARNING: geo_box.crs must be UTM!!!

    :param ancillary_mask:
        The path to the directory containing the land/sea data files.

    :return:
        A 2D Numpy Boolean array. True = Land, False = Sea.

    :note:
        The function does not currently support reprojections. The
        GriddedGeoBox must have CRS and Pixelsize matching the
        ancillary data GeoTiffs.

    :TODO:
        Support reprojection to any arbitrary GriddedGeoBox.
    """

    def img2map(geoTransform, pixel):
        """
        Converts a pixel (image) co-ordinate into a map co-ordinate.
        :param geoTransform:
            The Image co-ordinate information (upper left coords, offset
            and pixel sizes).

        :param pixel:
            A tuple containg the y and x image co-ordinates.

        :return:
            A tuple containg the x and y map co-ordinates.
        """

        mapx = pixel[1] * geoTransform[1] + geoTransform[0]
        mapy = geoTransform[3] - (pixel[0] * (numpy.abs(geoTransform[5])))

        return (mapx, mapy)

    def map2img(geoTransform, location):
        """
        Converts a map co-ordinate into a pixel (image) co-ordinate.

        :param geoTransform:
            The Image co-ordinate information (upper left coords, offset
            and pixel sizes).

        :param location:
            A tuple containg the x and y map co-ordinates.

        :return:
            A tuple containg the y and x image co-ordinates.
        """

        imgx = int(numpy.round((location[0] - geoTransform[0]) / geoTransform[1]))
        imgy = int(
            numpy.round((geoTransform[3] - location[1]) / numpy.abs(geoTransform[5]))
        )
        return (imgy, imgx)

    # get Land/Sea data file for this bounding box
    utm_zone = geo_box.crs.GetUTMZone()

    rasfile = os.path.join(ancillary_path, "WORLDzone%02d.tif" % abs(utm_zone))
    assert os.path.exists(rasfile), "ERROR: Raster File Not Found (%s)" % rasfile

    md = extract_ancillary_metadata(rasfile)
    md["data_source"] = "Rasterised Land/Sea Mask"
    md["data_file"] = rasfile
    metadata = {"land_sea_mask": md}

    geoTransform = geo_box.transform.to_gdal()
    if geoTransform is None:
        raise Exception("Image geotransformation Info is needed")

    dims = geo_box.shape

    lsobj = gdal.Open(rasfile, gdal.gdalconst.GA_ReadOnly)
    ls_geoT = lsobj.GetGeoTransform()

    # Convert the images' image co-ords into map co-ords
    mUL = img2map(geoTransform=geoTransform, pixel=(0, 0))
    mLR = img2map(geoTransform=geoTransform, pixel=(dims[0], dims[1]))

    # Convert the map co-ords into the rasfile image co-ords
    iUL = map2img(geoTransform=ls_geoT, location=mUL)
    iLR = map2img(geoTransform=ls_geoT, location=mLR)

    xoff = iUL[1]
    yoff = iUL[0]
    xsize = iLR[1] - xoff
    ysize = iLR[0] - yoff

    # Read in the land/sea array
    ls_arr = lsobj.ReadAsArray(xoff, yoff, xsize, ysize)

    return (ls_arr.astype("bool"), metadata)


def set_land_sea_bit(
    gridded_geo_box,
    pq_const,
    pqaResult,
    ancillary_path="/g/data/v10/eoancillarydata/Land_Sea_Rasters",
):

    mask, md = calc_land_sea_mask(gridded_geo_box, ancillary_path)
    bit_index = pq_const.land_sea
    pqaResult.set_mask(mask, bit_index)

    return md
