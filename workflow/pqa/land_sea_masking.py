#! /usr/bin/env python

import os, numpy, logging
from osgeo import gdal
from ULA3.dataset import SceneDataset
from ULA3.utils import dump_array
from ULA3.common.pqa_result import PQAResult
from ULA3.image_processor import ProcessorConfig
from ULA3.image_processor import constants
from ULA3 import DataManager

logger = logging.getLogger('root.' + __name__)

def process(subprocess_list=[], resume=False):
    logger.info('%s.process(%s, %s) called', __name__, subprocess_list, resume)

    CONFIG = ProcessorConfig()
    DATA = DataManager()

    l1t_input_dataset = DATA.get_item(CONFIG.input['l1t']['path'], SceneDataset)
    assert l1t_input_dataset, 'Unable to retrieve SceneDataset object for L1T input scene dataset'
    logger.debug( 'SceneDataset object for %s retrieved', l1t_input_dataset.pathname)

    l1t_stack = DATA.get_item('l1t_stack', numpy.ndarray)
    assert l1t_stack is not None, 'Unable to retrieve ndarray object for l1t_stack'
    logger.debug( 'ndarray object for l1t_stack retrieved')

    result = DATA.get_item('result.tif', PQAResult)
    assert result, 'Unable to retrieve PQAResult object for result'
    logger.debug( 'PQAResult object for result retrieved')

    def LandSeaMasking(image_stack, input_dataset):
        """
        Creates a Land/Sea mask.

        :param image_stack:
            An ordered array of all bands.

        :param input_dataset:
            A 2D Numpy array of the Land/Sea file.

        :return:
            A 2D Numpy Boolean array. True = Land, False = Sea.

        :author:
            Josh Sixsmith, joshua.sixsmith@ga.gov.au
        """

        def select_rasterfile(utm_zone):
            """
            Selects a shapefile corresponding to the UTM zone.

            :param utm_zone:
                The UTM zone number as an integer.

            :return:
                A rasterfile name that corresponds to the input zone.
            """

            return os.path.join(CONFIG.DIR_LandSea, 'WORLDzone%02d.tif' % utm_zone)

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


        rasfile = select_rasterfile(input_dataset.zone)
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

    mask = LandSeaMasking(l1t_stack, l1t_input_dataset)

    pq_const = constants.pqaContants(l1t_input_dataset.sensor)
    #bit_index = CONFIG.pqa_test_index['LAND_SEA']
    bit_index = pq_const.land_sea
    result.set_mask(mask, bit_index)
    if CONFIG.debug:
        dump_array(mask,
                   os.path.join(CONFIG.work_path, 'mask_%02d.tif' % bit_index),
                   l1t_input_dataset)
