import os
import numpy
import logging
from osgeo import gdal

def calc_land_sea_mask(image_stack, input_dataset, ancillary_path):
    """
    Creates a Land/Sea mask.

    :param image_stack:
        An ordered array of all bands.

    :param input_dataset:
        A 2D Numpy array of the Land/Sea file.

    :param ancillary_path:
        Path to land/sea ancillary data directory

    :return:
        A 2D Numpy Boolean array. True = Land, False = Sea.

    :author:
        Josh Sixsmith, joshua.sixsmith@ga.gov.au
    """

    def img2map(geoTransform, pixel):
        """
        Converts a pixel (image) co-ordinate into a map co-ordinate.

        :param geoTransform:
            The Image co-ordinate information (upper left coords, offset
            and pixel sizes)

        :param pixel:
            A tuple containg the y and x image co-ordinates.

        :return:
            A tuple containg the x and y map co-ordinates.
        """

        mapx = pixel[1] * geoTransform[1] + geoTransform[0]
        mapy = geoTransform[3] - (pixel[0] * (numpy.abs(geoTransform[5])))

        return (mapx,mapy)

    def map2img(geoTransform, location):
        """
        Converts a map co-ordinate into a pixel (image) co-ordinate.

        :param geoTransform:
            The Image co-ordinate information (upper left coords, offset
            and pixel sizes)

        :param location:
            A tuple containg the x and y map co-ordinates.

        :return:
            A tuple containg the y and x image co-ordinates.
        """

        imgx = int(numpy.round((location[0] - geoTransform[0])/geoTransform[1]))
        imgy = int(numpy.round((geoTransform[3] - location[1])/numpy.abs(geoTransform[5])))
        return (imgy,imgx)


    rasfile = os.path.join(ancillary_path, 'WORLDzone%02d.tif' % input_dataset.zone)
    assert os.path.exists(rasfile), 'ERROR: Raster File Not Found (%s)' % rasfile

    prj = input_dataset.GetProjection()
    geoTransform = input_dataset.GetGeoTransform()

    if prj == None: raise Exception('Image projection Infomation is needed')
    if geoTransform == None: raise Exception('Image geotransformation Info is needed')

    dims = image_stack.shape
    if len(dims) >2:
        ncols = dims[2]
        nrows = dims[1]
        dims  = (nrows,ncols)

    lsobj   = gdal.Open(rasfile, gdal.gdalconst.GA_ReadOnly)
    ls_geoT = lsobj.GetGeoTransform()

    # Convert the images' image co-ords into map co-ords
    mUL = img2map(geoTransform=geoTransform, pixel=(0,0))
    mLR = img2map(geoTransform=geoTransform, pixel=(dims[0],dims[1]))

    # Convert the map co-ords into the rasfile image co-ords
    iUL = map2img(geoTransform=ls_geoT, location=mUL)
    iLR = map2img(geoTransform=ls_geoT, location=mLR)

    xoff = iUL[1]
    yoff = iUL[0]
    xsize = iLR[1] - xoff
    ysize = iLR[0] - yoff

    # Read in the land/sea array
    ls_arr = lsobj.ReadAsArray(xoff, yoff, xsize, ysize)

    return (ls_arr.astype('bool'))

def setLandSeaBit(l1t_data, l1t_dataset, ancillary_path, pq_const, pqaResult):
    mask = calc_land_sea_mask(l1t_data, l1t_dataset, ancillary_path)
    bit_index = pq_const.land_sea
    pqaResult.set_mask(mask, bit_index)
