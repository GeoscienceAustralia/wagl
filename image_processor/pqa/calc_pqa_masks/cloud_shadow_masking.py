#! /usr/bin/env python


import os, datetime, gc, logging
import numpy
import numexpr
from scipy import ndimage
from osgeo import gdal
import osr

from ULA3.image_processor import ProcessorConfig
from ULA3 import DataManager
from ULA3.IDL_functions import IDL_Histogram

CONFIG = ProcessorConfig()
DATA = DataManager()
logger = logging.getLogger('root.' + __name__)

def Cloud_Shadow(image_stack, kelvin_array, cloud_mask, input_dataset,
                 land_sea_mask=None, contiguity_mask=None, cloud_algorithm='ACCA', growregion=False):
    """
    Identifies cloud shadow and creates a mask.

    Uses the sun position to find the shadow direction, and estimated
    cloud height to project a cloud shadows location.

    :param image_stack:
        An ordered array of all bands. The image co-ordinate system needs to be
        in metres.

    :param cloud_mask:
        A 2D Numpy array where 0 = Cloud and 1 = Not Cloud. A boolean of
        False = Cloud and True = Not Cloud is also accepted.

    :param input_dataset:
        SceneDataset object for input dataset.

    :param land_sea_mask:
        A 2D Numpy array where 0 = Sea and 1 = Land.

    :param contiguity_mask:
        A 2D Numpy array where 0 = Non-contiguous and 1 = Contiguous.

    :param cloud_algorithm:
        String indicating which cloud algorithm was used as the input:

    :param growregion:
        Whether to perform region growing on the projected cloud pixels.
        Defaulted to return the cloud projection only (set to False).

    :return:
        An 2D Numpy array mask with 0 for Shadow and the relevant bit specified
        in bitpos for Not Shadow.

    :author:
       Josh Sixsmith, joshua.sixsmith@ga.gov.au
    """

    # Distinguish between potentially concurrent executions for different cloud masks
    global logger
    logger = logging.getLogger('root.' + __name__ + '.' + cloud_algorithm.lower())

    pqa_temp_output = DATA.get_item('pqa_temp_output.dat', str)
    assert pqa_temp_output, 'Unable to retrieve string object for pqa_temp_output'
    logger.debug( 'string object for pqa_temp_output retrieved')

    if len(image_stack) == 0: return None
    if type(image_stack[0]) != numpy.ndarray: raise Exception('Array input is not valid')
    if cloud_mask == None: raise Exception('Cloud Layer input is not valid')

    geoTransform = input_dataset.GetGeoTransform()
    projection = input_dataset.GetProjection()

#-------------------------Filter Thresholds-----------------------------------

    #===========================================================================
    # wt_ndvi = 0.1
    # wt_b4   = 0.04
    # wt_b5   = 0.05
    # vrat_th = 0.08
    # btt_th  = 293
    # rt_b3   = 0.4
    # rt_b4   = 0.6
    # srt_low = 0.9
    # srt_hi  = 1.3
    #===========================================================================
    wt_ndvi = CONFIG.pqa_param['cloud_shadow_wt_ndvi']
    wt_b4   = CONFIG.pqa_param['cloud_shadow_wt_b4']
    wt_b5   = CONFIG.pqa_param['cloud_shadow_wt_b5']
    vrat_th = CONFIG.pqa_param['cloud_shadow_vrat_th']
    btt_th  = CONFIG.pqa_param['cloud_shadow_btt_th']
    rt_b3   = CONFIG.pqa_param['cloud_shadow_rt_b3']
    rt_b4   = CONFIG.pqa_param['cloud_shadow_rt_b4']
    srt_low = CONFIG.pqa_param['cloud_shadow_srt_low']
    srt_hi  = CONFIG.pqa_param['cloud_shadow_srt_hi']

    # Returns the required line from a list of strings
    def linefinder(string_list, string = ""):
        """
        Searches a list for the specified string.

        :param string_list:
            A list containing searchable strings.

        :param string:
            User input containing the string to search.

        :return:
            The line containing the found string.
        """

        for line in string_list:
            if string in str(line):
                return line


#===============================================================================
#    # Reads the metadata file in order to extract the needed parameters
#    def read_metafile(metafile):
#        '''Opens the metadata file and extracs relevant parameters.
#
#           Args:
#               metafile: A full string path name to the metadata file.
#
#           Returns:
#               Dictionary containing the parameters.
#        '''
#
#        f         = open(metafile, 'r')
#        met_array = f.readlines()
#        f.close()
#
#        sfind   = linefinder(met_array, 'SUN_AZIMUTH')
#        s_azi   = float(sfind.split()[2])
#        sfind   = linefinder(met_array, 'SUN_ELEVATION')
#        s_elev  = float(sfind.split()[2])
#
#
#
#        params = {
#                     'Sun_Azimuth'      : s_azi,
#                     'Sun_Elevation'    : s_elev,
#                 }
#
#        print params
#
#        return params
#===============================================================================



    # Convert the origin pixel co-ordinates to map units
    def origin_map(geoTransform, cindex):
        """
        Converts the origin pixel co-ordinates to map units.

        :param geoTransform:
            Image co-ordinate information (upper left coords, offset and pixel
            sizes).

        :param cindex:
            The indices of the cloud locations.

        :return:
            The coordinates of every cloud pixel as a map units, in two Lists
            for the x and y positions.
        """

        omapx = numexpr.evaluate("a*b + c", {'a': cindex[1], 'b': numpy.float32(geoTransform[1]), 'c': numpy.float32(geoTransform[0]) })
        omapy = numexpr.evaluate("a - (b*abs(c))", { 'a': numpy.float32(geoTransform[3]), 'b': cindex[0], 'c': numpy.float32(geoTransform[5]) })
        return omapx, omapy

    # Determine the cloud height
    def cloud_height(ctherm, surface_temp, lapse_rate=numpy.float32(6.4)):
        """
        Determines the height of the cloud.

        Uses a standard environmental lapse rate, the surface temperature,
        and the temperature of the clouds.

        :param ctherm:
            The temperature of the cloud pixels.

        :param surface_temp:
            The surface temperature.

        :param lapse_rate:
            The environmental lapse rate.

        :return:
            An array of cloud heights.
        """

        result = numexpr.evaluate("((surface_temp - ctherm)/lapse_rate)*thou", { 'thou': numpy.float32(1000) }, locals())
        return result

    # Determine the length of the shadow
    def shadow_length(cheight, rad_elev):
        """
        Determines the length of the shadow cast by the cloud.

        :param cheight:
            The height of the cloud.

        :param rad_elev:
            The sun elevation angle in radians.

        :return:
            An array of shadow lengths.
        """

        return cheight/numpy.tan(rad_elev)

    # Retrieve the x and y location in terms of metres from origin
    def rect_xy(shad_length, rad_cor_az):
        """
        Retrieve the x and y distances of projected shadow.

        The distances are in metres from the originating cloud pixel.

        :param shad_length:
            The length of the shadow in metres.

        :param rad_cor_az:
            The corrected azimuth angle in radians.

        :return:
            A complex array containing the x and y distances from the
            originating cloud pixel.
        """
        rectxy = numpy.empty(shad_length.shape, dtype='complex64')
        # Could be done w/ numexpr in a single parallel statement
        # for i in xrange(len(shad_length)):
        #    rectxy[i] = cmath.rect(shad_length[i], rad_cor_az)

        #et = datetime.datetime.now()
        #print 'rect_xy time taken: ', et - st
        #'''OR
        rectxy.real = numexpr.evaluate("shad_length * cos(rad_cor_az)")
        rectxy.imag = numexpr.evaluate("shad_length * sin(rad_cor_az)")
        #'''

        return rectxy

    # Convert the rectangular xy locations into map co-ordinates
    def mapxy(rect_xy, omapx, omapy):
        """
        Convert the x and y locations into map co-ordinates.

        :param rect_xy:
            The complex array containing the x and y distances.

        :param omapx:
            The originating x position of the cloud.

        :param omapy:
            The originating y position of the cloud.

        :return:
            Two arrays, one for the x co-ordinate, the second for the
            y co-ordinate.
        """

        new_mapx = numexpr.evaluate("omapx + real(rect_xy)").astype('float32')
        new_mapy = numexpr.evaluate("omapy + imag(rect_xy)").astype('float32')
        return new_mapx, new_mapy

    # Convert the new mapxy locations to image co-ordinates
    def map2img(new_mapx, new_mapy, geoTransform, dims):
        """
        Converts the x and y map locations in image co-ordinates.

        :param new_mapx:
            The x projected shadow location.

        :param new_mapy:
            The y projected shadow location.

        :param geoTransform:
            The Image co-ordinate information (upper left coords, offset
            and pixel sizes)

        :return:
             A tuple containing the indices of the shadow locations.
        """

        dct = { 'a': numpy.float32(geoTransform[0]), 'b': numpy.float32(geoTransform[1]) }
        imgx = numpy.round(numexpr.evaluate("(new_mapx - a)/b", dct, locals())).astype('int32')

        dct = { 'a': numpy.float32(geoTransform[3]), 'b': numpy.float32(abs(geoTransform[5])) }
        imgy = numpy.round(numexpr.evaluate("(a - new_mapy)/b", dct, locals())).astype('int32')

        mask = numexpr.evaluate("(imgx>=0) & (imgy>=0) & (imgx<d1) & (imgy<d0)", { 'd1': dims[1], 'd0': dims[0] }, locals())

        imgx = imgx[mask]
        imgy = imgy[mask]
        del mask

        return (imgy, imgx)

    def ndvi(red, nir):
        """
        The NDVI function calculates the Normalised Differenced Vegetation Index.

        :param red:
            Band 3 of Landsat TM/ETM+.

        :param nir:
            Band 4 of Landsat TM/ETM+.

        :return:
            An 2D array in the range of -1 to 1.
        """

        ndvi = numexpr.evaluate("(nir - red) / (nir + red)")
        return ndvi

    def WATER_TEST(ndvi, band5):
        """
        The WATER_TEST function is used to identify water pixels.

        :param ndvi:
             Normalised Differenced Vegetation Index.

        :param band5:
             Band 5 of Landsat TM/ETM+

        :return:
             An 2D array of type 'bool'.
        """

        wt = numexpr.evaluate("((ndvi < wt_ndvi) & (band5 < wt_b5))", {'wt_b5' : wt_b5, 'wt_ndvi' : wt_ndvi}, locals())

        return wt

    def mndwi(red, band5):
        """
        Modified Normalised Difference Water Index.

        Used to identify water pixels.

        :param red:
            Band 3 of Landsat TM/ETM+.

        :param band5:
            Band 5 of Landsat TM/ETM+.

        :return:
            An 2D array of type 'bool'.
        """

        wt = numexpr.evaluate("((1 - (band5 / red)) / (1 + (band5 / red))) > mndwi_thresh", {'mndwi_thresh' : mndwi_thresh}, locals())
        #wt = numexpr.evaluate("((1 - (band5 / red)) / (1 + (band5 / red)))")

        return wt

    def stdev(b1, b2, b3, b4, b5, b7):
        """
        Calculates the standard deviation through the bands.

        Creating a bool array that is used to identify water pixels.

        :param b1:
            Band 1 of Landsat TM/ETM+.

        :param b2:
            Band 2 of Landsat TM/ETM+.

        :param b3:
            Band 3 of Landsat TM/ETM+.

        :param b4:
            Band 4 of Landsat TM/ETM+.

        :param b5:
            Band 5 of Landsat TM/ETM+.

        :param b7:
            Band 7 of Landsat TM/ETM+.

        :return:
            An 2D array of type 'bool'.
        """

        xbar = numexpr.evaluate("(b1 + b2 + b3 + b4 + b5 + b7) / 6")
        #stdv = numexpr.evaluate("(sqrt(((b1 - xbar)**2 + (b2 - xbar)**2 + (b3 - xbar)**2 + (b4 - xbar)**2 + (b5 - xbar)**2 + (b7 - xbar)**2) / 5)) < stdev_thresh", {'stdev_thresh' : stdev_thresh}, locals())
        stdv = numexpr.evaluate("(sqrt(((b1 - xbar)**2 + (b2 - xbar)**2 + (b3 - xbar)**2 + (b4 - xbar)**2 + (b5 - xbar)**2 + (b7 - xbar)**2) / 5))")


        return stdv


#---------------------------------------Processing Here------------------------
    global SHADOW_logfile

    SHADOW_logfile = open(os.path.join(pqa_temp_output, 'scene01', '%s_CLOUD_SHADOW_LOGFILE.txt' % cloud_algorithm), 'w')

    start_time = datetime.datetime.now()

    # Get the indices of cloud
    # need the actual indices rather than a boolean array
    cindex = numpy.where(cloud_mask == False)

    # Return mask with all true there is no cloud
    if len(cindex[0]) == 0:
        SHADOW_logfile.write('Cloud Cover = 0% No Shadow to detect!\n')
        SHADOW_logfile.write('Cloud Shadow Percent: 0.0\n')
        cshadow = numpy.ones(cloud_mask.shape, dtype='bool')
        end_time   = datetime.datetime.now()
        time       = end_time - start_time
        SHADOW_logfile.write('Cloud Shadow Process Time: %s\n' %time)
        SHADOW_logfile.close()

        return cshadow

    # Create a spatial reference from the projection string
    sr = osr.SpatialReference()
    sr.ImportFromWkt(projection)

#===========================================================================
#    if (type(metafile) != dict): # Is an actual file in which case read it.
#        parameters = read_metafile(metafile)
#        rad_elev     = numpy.radians(input_dataset.sun_elevation)
#
#        if input_dataset.sun_azimuth >= 180:
#            corrected_az = input_dataset.sun_azimuth - 180
#            if (sr.IsGeographic() == 0):
#                # azimuth is dealt in polar form
#                corrected_az = 360 - corrected_az + 90
#        else:
#            # azimuth < 180
#            corrected_az = input_dataset.sun_azimuth + 180
#            if (sr.IsGeographic() == 0):
#                # azimuth is dealt in polar form
#                corrected_az = 360 - corrected_az + 90
#
#    else: # The metafile is a dictionary
#        rad_elev     = numpy.radians(metafile['Sun_Elevation'])
#
#        if metafile['Sun_Azimuth'] >= 180:
#            corrected_az = metafile['Sun_Azimuth'] - 180
#            if (sr.IsGeographic() == 0):
#                # azimuth is dealt in polar form
#                corrected_az = 360 - corrected_az + 90
#        else:
#            # azimuth < 180
#            corrected_az = metafile['Sun_Azimuth'] + 180
#            if (sr.IsGeographic() == 0):
#                # azimuth is dealt in polar form
#                corrected_az = 360 - corrected_az + 90
#===========================================================================

    rad_elev     = numpy.radians(input_dataset.sun_elevation)
    if input_dataset.sun_azimuth >= 180:
        corrected_az = input_dataset.sun_azimuth - 180
    else: # azimuth < 180
        corrected_az = input_dataset.sun_azimuth + 180

    if (sr.IsGeographic() == 0):
        # azimuth is dealt in polar form
        corrected_az = 360 - corrected_az + 90

    rad_cor_az = numpy.radians(corrected_az)

    # Account for any cloud pixels that are non-contiguous
    if contiguity_mask is not None:
        #cloud_mask[(cloud_mask == 0) & (contiguity_mask == 0)] = 1
        _temp_array = numexpr.evaluate("(cloud_mask | contiguity_mask) == False")
        cloud_mask[_temp_array] = False
        del _temp_array; gc.collect()

    # Expecting surface reflectance with a scale factor of 10000
    #dim = image_stack.shape
    scaling_factor = numpy.float32(0.0001)
    reflectance_stack = image_stack.astype(numpy.float32)
#        reflectance_stack *= (scaling_factor * null_nan_array)
    reflectance_stack = numexpr.evaluate("reflectance_stack * scaling_factor")

    # Get the indices of cloud
    # need the actual indices rather than a boolean array
#        cindex = numpy.where(cloud_mask == False)
    ctherm = kelvin_array[cindex]
    #dims_index   = ctherm.shape
    dims   = kelvin_array.shape

    ndvi      = ndvi(red=reflectance_stack[2], nir=reflectance_stack[3])
    #ni        = ndvi > 0.5 # What if none satisfy?
    ni        = numexpr.evaluate("ndvi > 0.5") # What if none satisfy?
    if (ni.sum() == 0):
        # Then just take non-cloud pixels
        ni = numexpr.evaluate("cloud_mask == True")
    non_cloud = numexpr.evaluate("cloud_mask == True")

    # General water test
    wt   = WATER_TEST(ndvi=ndvi, band5=reflectance_stack[4])
    #wt = numexpr.evaluate("((ndvi < wt_ndvi) & (band5 < wt_b5))")


    if contiguity_mask is not None:
        null = numexpr.evaluate("contiguity_mask == False") # Non-contiguous pixels
        th_copy = kelvin_array.copy()
        th_copy[null] = 0
        gr_zero = numexpr.evaluate("th_copy > 0")
        survivors   = numexpr.evaluate("non_cloud & ni & gr_zero")
        del null, th_copy, gr_zero; gc.collect()
    else:
        survivors   = numexpr.evaluate("non_cloud & ni")

    del ni; gc.collect()


    surfaceTemp = numpy.mean(kelvin_array[survivors], dtype='float64') # What if too low?

    print 'Surface Temperature (Survivors): ', surfaceTemp

    del ndvi, non_cloud

    cshadow = numpy.zeros(dims, dtype='byte')
    # wet, standard and dry
    lapse_rates = numpy.array([CONFIG.pqa_param['cloud_shadow_lapse_wet'],
                               CONFIG.pqa_param['cloud_shadow_lapse_standard'],
                               CONFIG.pqa_param['cloud_shadow_lapse_dry']], dtype='float32')

    if (sr.IsGeographic() == 1):
        R = sr.GetSemiMajor()

        print 'Calculating Cloud Map Origin Coordinates'
        rlon, rlat = origin_map(geoTransform, cindex)
        rlon = numpy.radians(rlon)
        rlat = numpy.radians(rlat)

        i = 1
        for lr in lapse_rates:

            print 'Calculating Cloud Height'
            cheight = cloud_height(ctherm, surface_temp=surfaceTemp, lapse_rate=lr)
            #del ctherm; gc.collect()

            print 'Calculating Cloud Shadow Length'
            d = shadow_length(cheight, rad_elev)
            #del cheight; gc.collect()

            print 'Calculating Cloud Shadow Map Coordinates'
            rlat2 = numpy.arcsin(numpy.sin(rlat)*numpy.cos(d/R) + numpy.cos(rlat)*numpy.sin(d/R)*numpy.cos(rad_cor_az))

            rlon2 = rlon + numpy.arctan2(numpy.sin(rad_cor_az)*numpy.sin(d/R)*numpy.cos(rlat),numpy.cos(d/R)-numpy.sin(rlat)*numpy.sin(rlat))

            rlat2 = numpy.rad2deg(rlat2)
            rlon2 = numpy.rad2deg(rlon2)
            #del rlon, rlat, d; gc.collect()
            #del d; gc.collect()

            print 'Calculating Cloud Shadow Image Coordinates'
            sindex = map2img(new_mapx=rlon2, new_mapy=rlat2, geoTransform=geoTransform, dims=dims)
            #del rlon2, rlat2; gc.collect()

            cshadow[sindex] = i
            i += 1

        del cheight, d, rlon, rlat, rlon2, rlat2, ctherm; gc.collect()
        s_index = numexpr.evaluate("cshadow == True")

        '''
        # Testing projected pixels
        cshadow = numpy.zeros(dims, dtype='byte')
        cshadow[sindex] = 1
        driver = gdal.GetDriverByName("GTiff")
        outfile = driver.Create("geographics_prj_clouds.tif", dims[1], dims[0], 1, gdal.GDT_Byte)
        outband = outfile.GetRasterBand(1)
        outband.WriteArray(cshadow)
        outfile.SetGeoTransform(geoTransform)
        outfile.SetProjection(projection)
        outfile = None
        del cshadow
        '''
    else:
        print 'Calculating Cloud Map Origin Coordinates'
        omapx, omapy = origin_map(geoTransform, cindex)

        i=1
        for lr in lapse_rates:

            print 'Calculating Cloud Height'
            cheight = cloud_height(ctherm, surface_temp=surfaceTemp, lapse_rate=lr)
            #del ctherm; gc.collect()

            print 'Calculating Cloud Shadow Length'
            shad_length = shadow_length(cheight, rad_elev)
            del cheight; gc.collect()

            print 'Converting Polar to Rectangular Coordinates'
            rectxy = rect_xy(shad_length, rad_cor_az)
            del shad_length; gc.collect()

            print 'Calculating Cloud Shadow Map Coordinates'
            new_mapx, new_mapy = mapxy(rectxy, omapx, omapy)
            #del rectxy, omapx, omapy; gc.collect()
            del rectxy; gc.collect()
            #print 'len new_mapx: ', new_mapx.shape

            print 'Calculating Cloud Shadow Image Coordinates'
            sindex = map2img(new_mapx, new_mapy, geoTransform, dims)
            del new_mapx, new_mapy; gc.collect()

            cshadow[sindex] = i
            i += 1

        #sindex = numpy.where(cshadow == 1)
        #sindex = numpy.where(cshadow >= 1)
        #sindex = numexpr.evaluate("where(cshadow >= 1)")
        # or as a bool array
        s_index = numexpr.evaluate("cshadow >= 1")

        # Testing projected pixels
        del cshadow, ctherm, omapx, omapy; gc.collect()



    # Only apply spectral tests when growing a region
    # May need to add or change spectral tests to eliminate some landcovers.
    if growregion is True:
        print 'Applying Region Growing'
        stng = datetime.datetime.now()

        '''
        # Threshold tests
        #q1 = reflectance_stack[3][sindex] < 0.12 # NIR threshold
        #q2 = reflectance_stack[1][sindex] < 0.045 # Green threshold
        q3 = cloud_mask[sindex] >= 1
        if contiguity_mask != None:
            q4 = contiguity_mask[sindex] >= 1
            if land_sea_mask != None:
                q5 = land_sea_mask[sindex] >= 1
                #q6 = q1 & q2 & q3 & q4 & q5
                q6 = q3 & q4 & q5
                #del q1, q2, q3, q4, q5
                del q3, q4, q5
            else:
                #q6 = q1 & q2 & q3 & q4
                q6 = q3 & q4
                #del q1, q2, q3, q4
                del q3, q4
        elif land_sea_mask != None:
            q5 = land_sea_mask[sindex] >= 1
            #q6 = q1 & q2 & q3 & q5
            q6 = q3 & q5
            #del q1, q2, q3, q5
            del q3, q5
        else:
            #q6 = q1 & q2 & q3
            q6 = q2 & q3
            #del q1, q2, q3
            del q3

        qf = numpy.where(q6 == False)
        cy = numpy.delete(sindex[0], qf, 0)
        cx = numpy.delete(sindex[1], qf, 0)
        sindex = (cy,cx)

        del qf, cy, cx; gc.collect()
        '''

        q1 = numexpr.evaluate("(cloud_mask >= 1) & s_index")
        if contiguity_mask != None:
            q2 = numexpr.evaluate("(contiguity_mask >= 1) & s_index")
            if land_sea_mask != None:
                q3 = numexpr.evaluate("(land_sea_mask >= 1) & s_index")
                q4 = numexpr.evaluate("q1 & q2 & q3")
            else:
                q4 = numexpr.evaluate("q1 & q2")
        elif land_sea_mask != None:
            q3 = numexpr.evaluate("(land_sea_mask >= 1) & s_index")
            q4 = numexpr.evaluate("q1 & q3")
        else:
            q4 = q1

        sindex = q4


        s = [[1,1,1],[1,1,1],[1,1,1]]

        # Additional weights generation
        # Will create weights based on the 'slope' of the spectral curve
        # occurring between bands.
        # Using a unit length of 1 for distance between bands
        weights = numpy.ones(dims, dtype='int8')

        # band 3 -> 4 slope
        #slope = reflectance_stack[3] - reflectance_stack[2]
        slope = numexpr.evaluate("(b4 - b3) >= 0.11", {'b4':reflectance_stack[3], 'b3':reflectance_stack[2]})
        #weights[slope >= 0.1] += 1
        #weights[slope >= 0.11] += 1
        weights[slope] += 1

        # band 4 -> 5 slope
        #slope = numpy.abs(reflectance_stack[4] - reflectance_stack[3])
        slope = numexpr.evaluate("abs(b5 - b4) >= 0.005", {'b5':reflectance_stack[4], 'b4':reflectance_stack[3]})
        #weights[slope >= 0.05] += 1
        #weights[slope >= 0.055] += 1
        weights[slope] += 1

        # band 4 -> 7 slope
        #slope = (reflectance_stack[6] - reflectance_stack[3])/2
        slope = numexpr.evaluate("((b7 - b4) / 2) > 0.01", {'b7':reflectance_stack[6], 'b4':reflectance_stack[3]})
        #weights[slope > 0.01] += 1
        weights[slope] += 1
        slope = numexpr.evaluate("abs((b7 - b4) / 2) > 0.05", {'b7':reflectance_stack[6], 'b4':reflectance_stack[3]})
        #weights[numpy.abs(slope) >= 0.05] += 1
        weights[slope] += 1

        # General water test (calculated earlier)
        #ndvi = ndvi(red=reflectance_stack[2], nir=reflectance_stack[3])
        #wt   = WATER_TEST(ndvi=ndvi, band5=reflectance_stack[4])
        # like shadow, water is pretty low, so set weight very high
        weights[wt] += 9

        # Modified normalised difference water index
        wt = mndwi(red=reflectance_stack[2], band5=reflectance_stack[4])
        weights[wt] += 9

        # standard deviation thruogh spectral space
        stdv = stdev(b1=reflectance_stack[0], b2=reflectance_stack[1], b3=reflectance_stack[2], b4=reflectance_stack[3], b5=reflectance_stack[4], b7=reflectance_stack[6])
        # This is for water that is spectrally very flat and near zero
        BOOL = numexpr.evaluate("stdv < 0.008")
        # dilate to get water edges; tends to help with river systems
        ndimage.binary_dilation(BOOL, s, output=wt)
        weights[wt] += 15

        # some aussie native bushland is still picked up. Its stdv is generally
        # higher than shadow though. More testing needed!
        BOOL = numexpr.evaluate("stdv > 0.04")
        weights[BOOL] += 1

        del slope, wt, BOOL, stdv; gc.collect()

        # could implement another weight; avg of bands 1 -> 3 > 450? which
        # would minimise turbid water
        #avg_bands = (reflectance_stack[0] + reflectance_stack[1] + reflectance_stack[2]) / 3
        #weights[avg_bands > 0.045] += 1

        # another for deep water could be avg b4 -> b7 (exclude thermal) <100
        #avg_bands = (reflectance_stack[3] + reflectance_stack[4] + reflectance_stack[6]) / 3
        #weights[avg_bands < 0.01] -= 10

        #weight_sum = ((reflectance_stack[0] + reflectance_stack[1] + reflectance_stack[2] + reflectance_stack[3]) +
        #                2*(reflectance_stack[4] + reflectance_stack[6]))
        #weight_sum = (weights * (reflectance_stack[0] + reflectance_stack[1] + reflectance_stack[2] + reflectance_stack[3]) +
        #                2*(weights * (reflectance_stack[4] + reflectance_stack[6])))

        b1 = reflectance_stack[0]
        b2 = reflectance_stack[1]
        b3 = reflectance_stack[2]
        b4 = reflectance_stack[3]
        b5 = reflectance_stack[4]
        b7 = reflectance_stack[6]
        weight_sum = numexpr.evaluate("weights *(b1 + b2 + b3 + b4) + 2*(weights * (b5 + b7))")

        del b1, b2, b3, b4, b5, b7; gc.collect()

        del weights; gc.collect()


        flatwsum = weight_sum.flatten()


        h = IDL_Histogram(flatwsum, min=0, max=1.0, binsize=0.1, reverse_indices='ri', locations='loc') # may not need loc

        # This might be the only need for the bin locations i.e. loc
        binmean = numpy.zeros(h['loc'].shape[0])
        binstdv = numpy.zeros(h['loc'].shape[0])

        ri = h['ri']
        hist = h['histogram']

        for i in numpy.arange(h['loc'].shape[0]):
            if ri[i+1] > ri[i]:
                temp_array = flatwsum[ri[ri[i]:ri[i+1]]]
                binmean[i] = numpy.mean(temp_array)
                binstdv[i] = numpy.std(temp_array)


        limit = binstdv * 2.5
        upper = binmean + limit
        lower = binmean - limit

        print 'upper: ', upper
        print 'lower: ', lower

        #del binmean, binstdv, weight_sum; gc.collect()
        del binmean, binstdv, flatwsum, temp_array; gc.collect()

        grown_regions = numpy.zeros(dims, dtype='bool').flatten()

        # Global stats/masking method
        lmin = lower.min()
        umax = upper.max()
        #flatmask = (flatwsum >= lower.min()) & (flatwsum <= upper.max())
        #flatmask = numexpr.evaluate("(flatwsum >= lmin) & (flatwsum <= umax)")
        mask = numexpr.evaluate("(weight_sum >= lmin) & (weight_sum <= umax)")
        #mask = flatmask.reshape(dims)
        label_array, num_labels = ndimage.label(mask, structure=s)
        flat_label = label_array.flatten()
        labels = label_array[sindex]
        ulabels = numpy.unique(labels[labels > 0])

        if (ulabels.size > 0): # only apply if labels are identified
            max = numpy.max(ulabels) # this will fail if there are no labels
            h = IDL_Histogram(flat_label, min=0, max=max, reverse_indices='ri')

            # assign other variable names, makes it more readable when indexing
            hist = h['histogram']
            ri = h['ri']

            for i in numpy.arange(ulabels.shape[0]):
                if hist[ulabels[i]] == 0:
                    continue
                grown_regions[ri[ri[ulabels[i]]:ri[ulabels[i]+1]]] = 1

            cshadow = grown_regions.reshape(dims)

        else: # if no labels then output no shadow
            cshadow = numpy.zeros(dims, dtype='byte')


    else:
        cshadow = numpy.zeros(dims, dtype='byte')
        cshadow[sindex] = 1

    # Majority filtering; applying twice, makes a cleaner result
    s = [[1,1,1],[1,1,1],[1,1,1]]
    shadfilt = ndimage.convolve(cshadow, s)
    cshadow = numexpr.evaluate("shadfilt > 4")
    ndimage.convolve(cshadow, s, output=shadfilt)
    cshadow = numexpr.evaluate("shadfilt > 4")

    del shadfilt; gc.collect()

    # Where a shadow pixel is a cloud pixel, change to no shadow
    cshadow[cindex] = False

    # Where a sea pixel is a shadow pixel, change to no shadow
    if land_sea_mask != None:
        sea = numexpr.evaluate("land_sea_mask == 0")
        cshadow[sea] = False
        del sea

    # Where a null pixel is shadow, change to no shadow
    if contiguity_mask != None:
        null = numexpr.evaluate("contiguity_mask == 0")
        cshadow[null] = False
        del null

    Shadow_Percent = (float(cshadow.sum())/cshadow.size) * 100
    SHADOW_logfile.write('Cloud Shadow Percent: %f\n' %Shadow_Percent)


    end_time   = datetime.datetime.now()
    time       = end_time - start_time

    SHADOW_logfile.write('Cloud Shadow Process Time: %s\n' %time)
    SHADOW_logfile.close()

    # Invert array from 1 == "Cloud" to True = "No Cloud"
    return (~cshadow).astype(numpy.bool)

