"""
Cloud Shadow Masking
--------------------
"""
import datetime
import logging
import numpy as np
import numexpr
import gc

from scipy import ndimage
from IDL_functions import histogram


def cloud_shadow(image_stack, kelvin_array, cloud_mask, geo_box, sun_az_deg,
                 sun_elev_deg, pq_const, land_sea_mask=None,
                 contiguity_mask=None, cloud_algorithm='ACCA',
                 growregion=False, aux_data={}):
    """
    Identifies cloud shadow and creates a mask.

    Uses the sun position to find the shadow direction, and estimated
    cloud height to project a cloud shadows location.

    :param image_stack:
        An ordered array of all bands. The image co-ordinate system needs to
        be in metres.

    :param cloud_mask:
        A 2D np array where 0 = Cloud and 1 = Not Cloud. A boolean of
        False = Cloud and True = Not Cloud is also accepted.

    :param geo_box:
        An instance of GriddedGeoBox representing the spatial context of
        the required shadow mask

    :param sun_az_deg:
        the azimutth of the sun in degrees

    :param sun_elev_deg:
        the elevation of the sun in degrees

    :param pq_const:
        An instance of PQAConstants applicable to the reflectance stack
        supplied

    :param land_sea_mask:
        A 2D np array where 0 = Sea and 1 = Land.

    :param contiguity_mask:
        A 2D np array where 0 = Non-contiguous and 1 = Contiguous.

    :param cloud_algorithm:
        String indicating which cloud algorithm was used as the input:

    :param growregion:
        Whether to perform region growing on the projected cloud pixels.
        Defaulted to return the cloud projection only (set to False).

    :param aux_data:
        A dict into which the function will place interesting intermediate
        results and metrics generated during processing - refer to the
        code for details

    :return:
        An 2D np array mask with 0 for Shadow and the relevant bit
        specified in bitpos for Not Shadow.

    :author:
       Josh Sixsmith, joshua.sixsmith@ga.gov.au
    """

    # Distinguish between potentially concurrent executions for different
    # cloud masks

    if len(image_stack) == 0:
        return None
    if type(image_stack[0]) != np.ndarray:
        raise Exception('Array input is not valid')
    if cloud_mask == None:
        raise Exception('Cloud Layer input is not valid')

    geoTransform = geo_box.affine.to_gdal()

    # Filter Thresholds:
    # wt_ndvi = 0.1
    # wt_b4   = 0.04
    # wt_b5   = 0.05
    # vrat_th = 0.08
    # btt_th  = 293
    # rt_b3   = 0.4
    # rt_b4   = 0.6
    # srt_low = 0.9
    # srt_hi  = 1.3

    wt_ndvi = pq_const.cshadow_wt_ndvi
    wt_b4 = pq_const.cshadow_wt_b4
    wt_b5 = pq_const.cshadow_wt_b5
    vrat_th = pq_const.cshadow_vrat_th
    btt_th = pq_const.cshadow_btt_th
    rt_b3 = pq_const.cshadow_rt_b3
    rt_b4 = pq_const.cshadow_rt_b4
    srt_low = pq_const.cshadow_srt_low
    srt_hi = pq_const.cshadow_srt_hi
    lapse_wet = pq_const.cshadow_lapse_wet
    lapse_standard = pq_const.cshadow_lapse_standard
    lapse_dry = pq_const.cshadow_lapse_dry
    stdv_bush = pq_const.cshadow_stdv_native_bush
    stdv_wt = pq_const.cshadow_stdv_spectral_flat_water
    mndwi_thresh = pq_const.cshadow_mndwi_thresh
    dense_veg = pq_const.cshadow_dense_veg
    slope_b34 = pq_const.cshadow_slope_b34
    slope_b45 = pq_const.cshadow_slope_b45
    slope_b47a = pq_const.cshadow_slope_b47a
    slope_b47b = pq_const.cshadow_slope_b47b
    stdv_mltp = pq_const.cshadow_stdv_multiplier

    def linefinder(string_list, string=""):
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

        omapx = numexpr.evaluate("a*b + c",
                                 {'a': cindex[1],
                                  'b': np.float32(geoTransform[1]),
                                  'c': np.float32(geoTransform[0])})

        omapy = numexpr.evaluate("a - (b*abs(c))",
                                 {'a': np.float32(geoTransform[3]),
                                  'b': cindex[0],
                                  'c': np.float32(geoTransform[5])})
        return omapx, omapy

    def cloud_height(ctherm, surface_temp, lapse_rate=np.float32(6.4)):
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

        result = numexpr.evaluate("((surface_temp - ctherm)/lapse_rate)*thou",
                                  {'thou': np.float32(1000)}, locals())
        return result

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

        return cheight / np.tan(rad_elev)

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
        rectxy = np.empty(shad_length.shape, dtype='complex64')
        rectxy.real = numexpr.evaluate("shad_length * cos(rad_cor_az)")
        rectxy.imag = numexpr.evaluate("shad_length * sin(rad_cor_az)")

        return rectxy

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

        dct = {'a': np.float32(geoTransform[0]),
               'b': np.float32(geoTransform[1])}

        imgx = np.round(numexpr.evaluate("(new_mapx - a)/b",
                                         dct, locals())).astype('int32')

        dct = {'a': np.float32(geoTransform[3]),
               'b': np.float32(abs(geoTransform[5]))}
        imgy = np.round(numexpr.evaluate("(a - new_mapy)/b",
                                         dct, locals())).astype('int32')

        mask = numexpr.evaluate("(imgx>=0) & (imgy>=0) &"
                                "(imgx<d1) & (imgy<d0)",
                                {'d1': dims[1], 'd0': dims[0]}, locals())

        imgx = imgx[mask]
        imgy = imgy[mask]
        del mask

        return (imgy, imgx)

    def ndvi(red, nir):
        """
        The NDVI function calculates the Normalised Differenced Vegetation
        Index.

        :param red:
            Band 3 of Landsat TM/ETM+.

        :param nir:
            Band 4 of Landsat TM/ETM+.

        :return:
            An 2D array in the range of -1 to 1.
        """

        ndvi = numexpr.evaluate("(nir - red) / (nir + red)")
        return ndvi

    def water_test(ndvi, band5):
        """
        The water_test function is used to identify water pixels.

        :param ndvi:
             Normalised Differenced Vegetation Index.

        :param band5:
             Band 5 of Landsat TM/ETM+

        :return:
             An 2D array of type 'bool'.
        """

        wt = numexpr.evaluate("((ndvi < wt_ndvi) &"
                              "(band5 < wt_b5))",
                              {'wt_b5': wt_b5,
                               'wt_ndvi': wt_ndvi}, locals())

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

        wt = numexpr.evaluate("((1 - (band5 / red)) / (1 + (band5 / red))) >"
                              "mndwi_thresh",
                              {'mndwi_thresh': mndwi_thresh}, locals())

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
            An 2D array of type 'float32'.
        """

        xbar = numexpr.evaluate("(b1 + b2 + b3 + b4 + b5 + b7) / 6")
        stdv = numexpr.evaluate("(sqrt(((b1 - xbar)**2 + (b2 - xbar)**2 +"
                                "(b3 - xbar)**2 + (b4 - xbar)**2 + "
                                "(b5 - xbar)**2 + (b7 - xbar)**2) / 5))")

        return stdv

    def majority_filter(array, iterations=1):
        """Majority filter."""
        weights_array = [[1, 1, 1], [1, 1, 1], [1, 1, 1]]
        for _ in range(iterations):
            array = ndimage.convolve(array, weights_array)
            array = numexpr.evaluate("array > 4")
        return array


    start_time = datetime.datetime.now()

    # Get the indices of cloud
    # need the actual indices rather than a boolean array
    cindex = np.where(cloud_mask == False)

    # Return mask with all true there is no cloud
    if len(cindex[0]) == 0:
        aux_data['%s_cloud_shadow_percent' % (cloud_algorithm, )] = 0.0
        logging.info('Cloud Shadow Percent: 0.0')
        cshadow = np.ones(cloud_mask.shape, dtype='bool')
        end_time = datetime.datetime.now()
        time = end_time - start_time
        aux_data['%s_cloud_shadow_runtime' % (cloud_algorithm, )] = time
        return cshadow

    # Create a spatial reference from the projection string
    # sr = osr.SpatialReference()
    # sr.ImportFromWkt(projection)
    sr = geo_box.crs

    rad_elev = np.radians(sun_elev_deg)
    if sun_az_deg >= 180:
        corrected_az = sun_az_deg - 180
    else:  # azimuth < 180
        corrected_az = sun_az_deg + 180

    if sr.IsGeographic() == 0:
        # azimuth is dealt in polar form
        corrected_az = 360 - corrected_az + 90

    rad_cor_az = np.radians(corrected_az)

    # Account for any cloud pixels that are non-contiguous
    if contiguity_mask is not None:
        #cloud_mask[(cloud_mask == 0) & (contiguity_mask == 0)] = 1
        _temp_array = numexpr.evaluate(
            "(cloud_mask | contiguity_mask) == False")
        cloud_mask[_temp_array] = False
        del _temp_array
        gc.collect()

    # Expecting surface reflectance with a scale factor of 10000
    #dim = image_stack.shape
    scaling_factor = np.float32(0.0001)
    reflectance_stack = image_stack.astype(np.float32)
    reflectance_stack = numexpr.evaluate("reflectance_stack * scaling_factor")

    # Get the indices of cloud
    # need the actual indices rather than a boolean array
    ctherm = kelvin_array[cindex]
    #dims_index   = ctherm.shape
    dims = kelvin_array.shape

    ndvi = ndvi(red=reflectance_stack[2], nir=reflectance_stack[3])
    ni = numexpr.evaluate("ndvi > dense_veg")  # What if none satisfy?
    if ni.sum() == 0:
        # Then just take non-cloud pixels
        ni = numexpr.evaluate("cloud_mask == True")
    non_cloud = numexpr.evaluate("cloud_mask == True")

    # General water test
    wt = water_test(ndvi=ndvi, band5=reflectance_stack[4])
    #wt = numexpr.evaluate("((ndvi < wt_ndvi) & (band5 < wt_b5))")

    if contiguity_mask is not None:
        # Non-contiguous pixels
        null = numexpr.evaluate("contiguity_mask == False")
        th_copy = kelvin_array.copy()
        th_copy[null] = 0
        gr_zero = numexpr.evaluate("th_copy > 0")
        survivors = numexpr.evaluate("non_cloud & ni & gr_zero")
        del null, th_copy, gr_zero
        gc.collect()
    else:
        survivors = numexpr.evaluate("non_cloud & ni")

    del ni
    gc.collect()

    surfaceTemp = np.mean(
        kelvin_array[survivors], dtype='float64')  # What if too low?

    print 'Surface Temperature (Survivors): ', surfaceTemp

    del ndvi, non_cloud

    cshadow = np.zeros(dims, dtype='byte')
    # wet, standard and dry
    lapse_rates = np.array([lapse_wet,
                            lapse_standard,
                            lapse_dry], dtype='float32')

    if sr.IsGeographic() == 1:
        R = sr.GetSemiMajor()

        print 'Calculating Cloud Map Origin Coordinates'
        rlon, rlat = origin_map(geoTransform, cindex)
        rlon = np.radians(rlon)
        rlat = np.radians(rlat)

        i = 1
        for lr in lapse_rates:

            print 'Calculating Cloud Height'
            cheight = cloud_height(
                ctherm, surface_temp=surfaceTemp, lapse_rate=lr)
            #del ctherm; gc.collect()

            print 'Calculating Cloud Shadow Length'
            d = shadow_length(cheight, rad_elev)
            #del cheight; gc.collect()

            print 'Calculating Cloud Shadow Map Coordinates'
            rlat2 = np.arcsin(np.sin(rlat) * np.cos(d / R) +
                              np.cos(rlat) * np.sin(d / R) *
                              np.cos(rad_cor_az))

            rlon2 = rlon + np.arctan2(np.sin(rad_cor_az) * np.sin(d / R)
                                      * np.cos(rlat), np.cos(d / R) -
                                      np.sin(rlat) * np.sin(rlat))

            rlat2 = np.rad2deg(rlat2)
            rlon2 = np.rad2deg(rlon2)

            print 'Calculating Cloud Shadow Image Coordinates'
            sindex = map2img(new_mapx=rlon2, new_mapy=rlat2,
                             geoTransform=geoTransform, dims=dims)

            cshadow[sindex] = i
            i += 1

        del cheight, d, rlon, rlat, rlon2, rlat2, ctherm
        gc.collect()
        s_index = numexpr.evaluate("cshadow == True")

    else:
        print 'Calculating Cloud Map Origin Coordinates'
        omapx, omapy = origin_map(geoTransform, cindex)

        i = 1
        for lr in lapse_rates:

            print 'Calculating Cloud Height'
            cheight = cloud_height(
                ctherm, surface_temp=surfaceTemp, lapse_rate=lr)
            #del ctherm; gc.collect()

            print 'Calculating Cloud Shadow Length'
            shad_length = shadow_length(cheight, rad_elev)
            del cheight
            gc.collect()

            print 'Converting Polar to Rectangular Coordinates'
            rectxy = rect_xy(shad_length, rad_cor_az)
            del shad_length
            gc.collect()

            print 'Calculating Cloud Shadow Map Coordinates'
            new_mapx, new_mapy = mapxy(rectxy, omapx, omapy)
            #del rectxy, omapx, omapy; gc.collect()
            del rectxy
            gc.collect()
            # print 'len new_mapx: ', new_mapx.shape

            print 'Calculating Cloud Shadow Image Coordinates'
            sindex = map2img(new_mapx, new_mapy, geoTransform, dims)
            del new_mapx, new_mapy
            gc.collect()

            cshadow[sindex] = i
            i += 1

        s_index = numexpr.evaluate("cshadow >= 1")

        del cshadow, ctherm, omapx, omapy
        gc.collect()

    # Only apply spectral tests when growing a region
    # May need to add or change spectral tests to eliminate some landcovers.
    if growregion:
        print 'Applying Region Growing'
        stng = datetime.datetime.now()

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

        s = [[1, 1, 1], [1, 1, 1], [1, 1, 1]]

        # Additional weights generation
        # Will create weights based on the 'slope' of the spectral curve
        # occurring between bands.
        # Using a unit length of 1 for distance between bands
        weights = np.ones(dims, dtype='int8')

        # band 3 -> 4 slope
        slope = numexpr.evaluate("(b4 - b3) >= slope_b34",
                                 {'b4': reflectance_stack[3],
                                  'b3': reflectance_stack[2]}, locals())
        weights[slope] += 1

        # band 4 -> 5 slope
        slope = numexpr.evaluate("abs(b5 - b4) >= slope_b45",
                                 {'b5': reflectance_stack[4],
                                  'b4': reflectance_stack[3]}, locals())
        weights[slope] += 1

        # band 4 -> 7 slope
        slope = numexpr.evaluate("((b7 - b4) / 2) > slope_b47a",
                                 {'b7': reflectance_stack[5],
                                  'b4': reflectance_stack[3]}, locals())
        weights[slope] += 1
        slope = numexpr.evaluate("abs((b7 - b4) / 2) > slope_b47b",
                                 {'b7': reflectance_stack[5],
                                  'b4': reflectance_stack[3]}, locals())
        weights[slope] += 1

        # General water test (calculated earlier)
        weights[wt] += 9

        # Modified normalised difference water index
        wt = mndwi(red=reflectance_stack[2], band5=reflectance_stack[4])
        weights[wt] += 9

        # standard deviation thruogh spectral space
        stdv = stdev(b1=reflectance_stack[0], b2=reflectance_stack[1],
                     b3=reflectance_stack[2], b4=reflectance_stack[3],
                     b5=reflectance_stack[4], b7=reflectance_stack[5])

        # This is for water that is spectrally very flat and near zero
        above = numexpr.evaluate("stdv < stdv_wt")
        # dilate to get water edges; tends to help with river systems
        ndimage.binary_dilation(above, s, output=wt)
        weights[wt] += 15

        # some aussie native bushland is still picked up. Its stdv is
        # generally higher than shadow though. More testing needed!
        above = numexpr.evaluate("stdv > stdv_bush")
        weights[above] += 1

        del slope, wt, above, stdv
        gc.collect()

        b1 = reflectance_stack[0]
        b2 = reflectance_stack[1]
        b3 = reflectance_stack[2]
        b4 = reflectance_stack[3]
        b5 = reflectance_stack[4]
        b7 = reflectance_stack[5]
        weight_sum = numexpr.evaluate("weights *(b1 + b2 + b3 + b4)"
                                      "+ 2*(weights * (b5 + b7))")

        del b1, b2, b3, b4, b5, b7
        gc.collect()

        del weights
        gc.collect()

        flatwsum = weight_sum.flatten()

        h = histogram(flatwsum, min=0, max=1.0, binsize=0.1,
                      reverse_indices='ri', locations='loc')

        # This might be the only need for the bin locations i.e. loc
        binmean = np.zeros(h['loc'].shape[0])
        binstdv = np.zeros(h['loc'].shape[0])

        ri = h['ri']
        hist = h['histogram']

        for i in np.arange(h['loc'].shape[0]):
            if ri[i + 1] > ri[i]:
                temp_array = flatwsum[ri[ri[i]:ri[i + 1]]]
                binmean[i] = np.mean(temp_array)
                binstdv[i] = np.std(temp_array)

        # Define a fuzzy limit
        limit = binstdv * stdv_mltp
        upper = binmean + limit
        lower = binmean - limit

        print 'upper: ', upper
        print 'lower: ', lower

        #del binmean, binstdv, weight_sum; gc.collect()
        del binmean, binstdv, flatwsum, temp_array
        gc.collect()

        grown_regions = np.zeros(dims, dtype='bool').flatten()

        # Global stats/masking method
        lmin = lower.min()
        umax = upper.max()
        mask = numexpr.evaluate("(weight_sum >= lmin) & (weight_sum <= umax)")
        label_array, num_labels = ndimage.label(mask, structure=s)
        flat_label = label_array.flatten()
        labels = label_array[sindex]
        ulabels = np.unique(labels[labels > 0])

        if ulabels.size > 0:  # only apply if labels are identified
            maxval = np.max(ulabels)  # this will fail if there are no labels
            h = histogram(flat_label, min=0, max=maxval, reverse_indices='ri')

            # assign other variable names, makes it more readable
            # when indexing
            hist = h['histogram']
            ri = h['ri']

            for i in np.arange(ulabels.shape[0]):
                if hist[ulabels[i]] == 0:
                    continue
                grown_regions[ri[ri[ulabels[i]]:ri[ulabels[i] + 1]]] = 1

            cshadow = grown_regions.reshape(dims)

        else:  # if no labels then output no shadow
            cshadow = np.zeros(dims, dtype='byte')

    else:
        cshadow = np.zeros(dims, dtype='byte')
        cshadow[sindex] = 1

    cshadow = majority_filter(array=cshadow, iterations=2)

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

    shadow_percent = (float(cshadow.sum()) / cshadow.size) * 100
    aux_data['%s_cloud_shadow_percent' % (cloud_algorithm, )] = shadow_percent
    logging.info('%s_cloud_shadow_percent: %f',
                 cloud_algorithm, shadow_percent)

    end_time = datetime.datetime.now()
    time = end_time - start_time

    aux_data['%s_cloud_shadow_runtime' % (cloud_algorithm, )] = time

    # Invert array from 1 == "Cloud" to True = "No Cloud"
    return (~cshadow).astype(np.bool)
