"""
ACCA Cloud Masking Algorithm
----------------------------

Follows the ACCA (Automated Cloud Cover Assessment) given by
Irish, R et al, October 2006, Photogrammetric Engineering & Remote
Sensing, Vol. 72, No. 10, October 2006, pp. 1179-1188.

"""
from __future__ import absolute_import, print_function
import datetime
import logging
import gc

import numexpr
import numpy

from scipy import ndimage

NAN = numpy.float32(numpy.NaN)


def water_test(reflectance_stack):
    """
    The WATER_TEST function is used to identify water pixels.

    :param reflectance_stack:
        An nD Numpy array of all bands (ordered).

    :return:
        An ndarray of type bool for water areas where True == water.
    """

    # water = (ndvi < 0.1) or ((nir < 0.04) and (swir1 < 0.05))
    # water = numpy.logical_or(ndvi < -01, numpy.logical_or(array[3,:,:]
    #                             < 0.04, array[4,:,:] < 0.05))

    return numexpr.evaluate("(NDVI < 0.1) & (band5_array < 0.05)",
                            {'NDVI': ndvi(reflectance_stack),
                             'band5_array': reflectance_stack[4]}, locals())


def ndvi(reflectance_stack):
    """
    The NDVI function calculates the Normalised Differenced Vegetation Index.

    :param reflectance_stack:
        An nD Numpy array of all bands (ordered).

    :return:
        A Numpy nD array of NDVI values in the range of -1 to 1.
    """

    return numexpr.evaluate("(a - b) / (a + b)",
                            {'a': reflectance_stack[3],
                             'b': reflectance_stack[2]}, locals())


def ndsi(reflectance_stack):
    """
    Calculates the Normalised Snow Difference Index.

    :param reflectance_stack:
        An nD Numpy array of all bands (ordered).

    :return:
        An Numpy nD array of values in the range of -1 to 1.
    """

    # ndsi = (green - swir1)/(green + swir1)
    return numexpr.evaluate("(a - b) / (a + b)",
                            {'a': reflectance_stack[1],
                             'b': reflectance_stack[4]}, locals())


def filter4(reflectance_stack, thermal_array):
    """
    Calculates band 5/6 composite.

    Clouds are generally cold and highly reflective in band 5.

    :param reflectance_stack:
        An nD Numpy array of all bands (ordered).

    :returns:
        An nD Numpy  array containing the band 5/6 composite.
    """

    return numexpr.evaluate("(1.0 - swir1) * TM6",
                            {'swir1': reflectance_stack[4],
                             'TM6': thermal_array}, locals())


def filter5(reflectance_stack):
    """
    Calculates band 4/3 ratio.

    Eliminates highly reflective vegetation.

    :param reflectance_stack:
        An nD Numpy array of all bands (ordered).

    :return:
        An nD Numpy array (unbound range).
    """

    return numexpr.evaluate("nir / red",
                            {'nir': reflectance_stack[3],
                             'red': reflectance_stack[2]}, locals())


def filter6(reflectance_stack):
    """
    Calculates band 4/2 ratio.

    Eliminates senescing vegetaion.

    :param reflectance_stack:
        An nD Numpy array of all bands (ordered).

    :returns:
        An nD Numpy array (unbound range).
    """

    return numexpr.evaluate("nir / green",
                            {'nir': reflectance_stack[3],
                             'green': reflectance_stack[1]}, locals())


def filter7(reflectance_stack):
    """
    Calculates band 4/5 ratio.

    Eliminates highly reflective soils.

    :param reflectance_stack:
        An nD Numpy array of all bands (ordered).

    :return:
        An nD Numpy array (unbound range).
    """

    return numexpr.evaluate("nir / swir1",
                            {'nir': reflectance_stack[3],
                             'swir1': reflectance_stack[4]}, locals())


def skewness(cloud_thermal_array, mean_temp, stdv_temp, count):
    """
    Calculates the third moment about the mean (aka, skewness).

    :param cloud_thermal_array:
        The thermal band in kelvin masked for cloud.

    :param mean_temp:
        The mean temperature.

    :param stdv_temp:
        The standard deviation of the temperature.

    :param count:
        The count or number of pixels included.

    :return:
        Floating point value.
    """

    cubed_deviates = (cloud_thermal_array - mean_temp) ** 3
    sum_cubed_dv = numpy.sum(cubed_deviates, dtype='float64')
    cubed_stdv = stdv_temp ** 3
    return (sum_cubed_dv / cubed_stdv) / count


def acca_2nd_pass(cloud_mask, ambiguous_array, thermal_array,
                  mean_cloud_temp, pq_const, aux_data=None):
    """
    The second pass of the ACCA algorithm.

    :param cloud_mask:
        The current identifed cloud mask from the 1st pass.

    :param ambiguous_array:
        The pixels labeled as ambiguous from the 1st pass
        (1 = ambiguous).

    :param thermal_array:
        The thermal band of Landsat TM/ETM+ in un-scaled degrees
        Kelvin.

    :param mean_cloud_temp:
        The mean temperature of the currently identifed cloud pixels.

    :param pq_const:
        An instance of PQAConstants applicable to the reflectance stack
        supplied

    :param aux_data:
        A dict into which the function will place interesting intermediate
        results and metrics generated during processing - refer to the 
        code for details

    :return:
        Depending on the result of the second pass, the acca_second_pass
        function can return None (conditions not, therefore
        stick with the cloud identified in the first pass), or a cloud
        mask with with 1 as cloud and 0 as not cloud.
        Note: Any caller-supplied aux_data dict will be updated
    """

    aux_data = aux_data or {}  # initialise aux_data to a dictionary
    logging.info('ACCA Pass Two Engaged')
    aux_data['acca_pass_2'] = 'engaged'

    cloud_stddev = numpy.std(thermal_array[cloud_mask], dtype='float64',
                             ddof=1)
    cloud_count = numpy.sum(cloud_mask, dtype='float64')

    aux_data['acca_pass_2_sdev'] = cloud_stddev

    # Histogram Percentiles for new thermal thresholds
    upper = numpy.percentile(thermal_array[cloud_mask], 97.5)
    lower = numpy.percentile(thermal_array[cloud_mask], 83.5)
    upper_max = numpy.percentile(thermal_array[cloud_mask], 98.75)

    aux_data['acca_pass_2_97_5_percentile'] = upper
    aux_data['acca_pass_2_83_5_percentile'] = lower
    aux_data['acca_pass_2_98_75_percentile'] = upper_max

    # Test for negative skewness
    skew = skewness(thermal_array[cloud_mask], mean_temp=mean_cloud_temp,
                    stdv_temp=cloud_stddev, count=cloud_count)
    logging.debug('skew: %s', skew)

    aux_data['acca_pass_2_skewness'] = skew

    # Calculate threshold shift
    shift = skew * cloud_stddev
    if shift > 1:
        shift = 1

    if skew > 0:
        # change the upper and lower thresholds
        new_upper = upper + shift
        new_lower = lower + shift
        if new_upper > upper_max:
            if new_lower > upper_max:
                new_lower = (upper_max - upper)
        else:
            new_upper = upper_max

        #query  = (thermal_array > new_lower) & (thermal_array <= new_upper)
        query = numexpr.evaluate("((ambiguous_array * thermal_array)"
                                 "> new_lower) & ((ambiguous_array *"
                                 "thermal_array) <= new_upper)")
        query2 = numexpr.evaluate("((ambiguous_array * thermal_array)"
                                  "<= new_lower)")

        # Compute stats for each query/class
        # Max, Mean
        if query.any():
            qmax = thermal_array[query].max()
            qmean = thermal_array[query].mean()
        else:
            qmax = 295
            qmean = 295

        if query2.any():
            qmax2 = thermal_array[query2].max()
            qmean2 = thermal_array[query2].mean()
        else:
            qmax2 = 295
            qmean2 = 295

        aux_data['acca_pass_2_class_1_max'] = qmax
        aux_data['acca_pass_2_class_2_max'] = qmax2

        aux_data['acca_pass_2_class_1_mean'] = qmean
        aux_data['acca_pass_2_class_2_mean'] = qmean2

        # Class percentage of scene
        qpop = (float(query.sum()) / ambiguous_array.size) * 100
        qpop2 = (float(query2.sum()) / ambiguous_array.size) * 100

        aux_data['acca_pass_2_class_1_percent'] = qpop
        aux_data['acca_pass_2_class_2_percent'] = qpop2

        if qpop < pq_const.acca_thermal_effect:
            if qmean < pq_const.acca_cold_cloud_mean:
                # Combine all cloud classes
                return numexpr.evaluate("cloud_mask | query | query2")
            elif qpop2 < pq_const.acca_thermal_effect:
                if qmean2 < pq_const.acca_cold_cloud_mean:
                    # Combine lower threshold clouds and pass 1 clouds
                    return numexpr.evaluate("cloud_mask | query2")
        return None  # Keep first pass cloud

    else:
        query = numexpr.evaluate("((ambiguous_array * thermal_array)"
                                 "> lower) & ((ambiguous_array * "
                                 "thermal_array) <= upper)")
        query2 = numexpr.evaluate("((ambiguous_array * thermal_array)"
                                  "!= 0) & ((ambiguous_array * thermal_array)"
                                  "<= lower)")

        # Compute stats for each query/class
        # Max, Mean
        if query.any():
            qmax = thermal_array[query].max()
            qmean = thermal_array[query].mean()
        else:
            qmax = 295
            qmean = 295

        if query2.any():
            qmax2 = thermal_array[query2].max()
            qmean2 = thermal_array[query2].mean()
        else:
            qmax2 = 295
            qmean2 = 295

        aux_data['acca_pass_2_class_1_max'] = qmax
        aux_data['acca_pass_2_class_2_max'] = qmax2

        aux_data['acca_pass_2_class_1_mean'] = qmean
        aux_data['acca_pass_2_class_2_mean'] = qmean2

        # Class percentage of scene
        qpop = (float(query.sum()) / ambiguous_array.size) * 100
        qpop2 = (float(query2.sum()) / ambiguous_array.size) * 100

        aux_data['acca_pass_2_class_1_percent'] = qpop
        aux_data['acca_pass_2_class_2_percent'] = qpop2

        if qpop < pq_const.acca_thermal_effect:
            if qmean < pq_const.acca_cold_cloud_mean:
                # Combine all cloud classes
                return numexpr.evaluate("cloud_mask | query | query2")
            elif qpop2 < pq_const.acca_thermal_effect:
                if qmean2 < pq_const.acca_cold_cloud_mean:
                    # Combine lower threshold clouds and pass 1 clouds
                    return numexpr.evaluate("cloud_mask | query2")
        return None  # Keep fist pass cloud


def acca(reflectance_stack, thermal_array, potential_cloud_array, pq_const,
         aux_data=None):
    """
    The first pass processing of the ACCA algorithm.

    :param reflectance_stack:
        An nD numpy array containing the reflectance values for each
        band.

    :param thermal_array:
        A 2D Numpy array containing the thermal band.

    :param potential_cloud_array:
        A 2D Numpy array containing values of 1 for valid areas and NAN for
        invalid areas.

    :param pq_const:
        An instance of PQAConstants applicable to the reflectance stack
        supplied

    :param aux_data:
        A dict into which the function will place interesting intermediate
        results and metrics generated during processing - refer to the 
        code for details

    :return:
        An 2D Numpy array cloud mask, with True for non-cloud and False
        for cloud.
        Note: Any caller-supplied aux_data dict will be updated
    """

    aux_data = aux_data or {}  # initialise aux_data to a dictionary
    dims = reflectance_stack.shape

    #===================================================================
    # b1 = image_stack[0,:,:]
    # b2 = image_stack[1,:,:]
    # b3 = image_stack[2,:,:]
    # b4 = image_stack[3,:,:]
    # b5 = image_stack[4,:,:]
    # b6 = image_stack[5,:,:]
    # b7 = image_stack[6,:,:]
    #===================================================================

    # Create the array for Ambiguous Pixels - NAN means not ambiguous
    ambiguous_array = numpy.ones(
        potential_cloud_array.shape, dtype=numpy.float32) * NAN

    # Will add in a water mask, to remove cold water bodies that have been
    # put into the ambigous group. If the water body is high in red
    # reflectance, it will have made it this far.

    water_mask = water_test(reflectance_stack)
    query = numexpr.evaluate("where((potential_cloud_array * b3)"
                             "> thresh_f1, 1, NAN)",
                             {'b3': reflectance_stack[2],
                              'thresh_f1': pq_const.acca_thresh_f1,
                              'NAN': NAN},
                             locals())

    potential_cloud_array *= query

    # Filter 2: NDSI
    ndsi_array = ndsi(reflectance_stack)
    query = numexpr.evaluate("where((potential_cloud_array * ndsi_array)"
                             "< thresh_f2, 1, NAN)",
                             {'thresh_f2': pq_const.acca_thresh_f2,
                              'NAN': NAN}, locals())

    # Find the snow pixels. Sum is used to find the total cloud pixels
    # as valid pixels = 1.  Sum of ones therefore = count
    find = numexpr.evaluate("where(ndsi_array >= thresh_f2, 1, 0)",
                            {'thresh_f2': pq_const.acca_thresh_f2}, locals())
    snow_pixels = find.sum()  # Sum is used as valid pixels = 1
    snow_percent = (float(snow_pixels) / find.size) * 100

    potential_cloud_array *= query

    aux_data['acca_pass_1_snow_percent'] = snow_percent

    # Filter 3; Temp. threshold
    query = numexpr.evaluate("where((potential_cloud_array * thermal_array)"
                             "< thresh_f3, 1, NAN)",
                             {'thresh_f3': pq_const.acca_thresh_f3,
                              'NAN': NAN},
                             locals())
    potential_cloud_array *= query

    # Filter 4; Band 5/6 composite
    _temporary = potential_cloud_array * \
        filter4(reflectance_stack, thermal_array)
    query = numexpr.evaluate("where(_temporary < thresh_f4, 1, NAN)",
                             {'thresh_f4': pq_const.acca_thresh_f4,
                              'NAN': NAN},
                             locals())

    # Get ambiguous pixels
    find = numexpr.evaluate("_temporary >= thresh_f4",
                            {'thresh_f4': pq_const.acca_thresh_f4}, locals())
    ambiguous_array[find] = 1

    potential_cloud_array *= query

    # Will add in a water mask, to remove cold water bodies that have been
    # put into the ambiguous group. If the water body is high in red
    # reflectance, it will have made it this far.

    ambiguous_array[water_mask] = NAN  # All water is unambiguous
    del water_mask
    gc.collect()

    # Filter 5; Band 4/3 ratio (Simple veg ratio)
    _temporary = potential_cloud_array * filter5(reflectance_stack)
    query = numexpr.evaluate("where(_temporary < thresh_f5, 1, NAN)",
                             {'thresh_f5': pq_const.acca_thresh_f5,
                              'NAN': NAN}, locals())

    # Get ambiguous pixels
    find = numexpr.evaluate("_temporary >= thresh_f5",
                            {'thresh_f5': pq_const.acca_thresh_f5}, locals())
    ambiguous_array[find] = 1

    potential_cloud_array *= query

    # Filter 6; Band 4/2 ratio (Dying/senescing veg)
    _temporary = potential_cloud_array * filter6(reflectance_stack)
    query = numexpr.evaluate("where(_temporary < thresh_f6, 1, NAN)",
                             {'thresh_f6': pq_const.acca_thresh_f6,
                              'NAN': NAN}, locals())

    # Tally filter 6 survivors
    f6_surv = numpy.nansum(query)

    # Get ambiguous pixels
    find = numexpr.evaluate("_temporary >= thresh_f6",
                            {'thresh_f6': pq_const.acca_thresh_f6,
                             'NAN': NAN}, locals())
    ambiguous_array[find] = 1

    potential_cloud_array *= query

    # Filter 7; Band 4/5 ratio (Identify highly reflective soils/rocks)
    # The results of this query are clouds at first pass
    _temporary = potential_cloud_array * filter7(reflectance_stack)

    query = numexpr.evaluate("where(_temporary > thresh_f7, 1, NAN)",
                             {'thresh_f7': pq_const.acca_thresh_f7,
                              'NAN': NAN}, locals())

    # Tally filter 7 survivors
    f7_surv = numpy.nansum(query)
    desert_index = float(f7_surv) / f6_surv

    aux_data['acca_pass_1_desert_index'] = desert_index

    # Get ambiguous pixels
    find = numexpr.evaluate("_temporary <= thresh_f7",
                            {'thresh_f7': pq_const.acca_thresh_f7}, locals())
    ambiguous_array[find] = 1

    potential_cloud_array *= query

    # Filter 8; Band 5/6 composite (Separate warm/cold clouds)
    _temporary = potential_cloud_array * \
        filter4(reflectance_stack, thermal_array)
    cold_cloud = numexpr.evaluate("_temporary < thresh_f8",
                                  {'thresh_f8': pq_const.acca_thresh_f8},
                                  locals())
    warm_cloud = numexpr.evaluate("_temporary >= thresh_f8",
                                  {'thresh_f8': pq_const.acca_thresh_f8},
                                  locals())

    cold_cloud_pop = (float(cold_cloud.sum()) / ambiguous_array.size) * 100
    cold_cloud_mean = numpy.mean(thermal_array[cold_cloud], dtype='float64')
    warm_cloud_pop = (float(warm_cloud.sum()) / ambiguous_array.size) * 100
    warm_cloud_mean = numpy.mean(thermal_array[warm_cloud], dtype='float64')

    aux_data['acca_pass_1_cold_cloud_percent'] = cold_cloud_pop
    aux_data['acca_pass_1_cold_cloud_mean'] = cold_cloud_mean
    aux_data['acca_pass_1_warm_cloud_percent'] = warm_cloud_pop
    aux_data['acca_pass_1_warm_cloud_mean'] = warm_cloud_mean

    del query, find, _temporary
    gc.collect()

    # Tests for snow and desert.  If the thresholds aren't breached, Pass two
    # is implemented.

    # REDO of tests for pass two engagement
    if desert_index <= pq_const.acca_desert_index and \
       snow_percent > pq_const.acca_snow_threshold:
        cloud = cold_cloud
        ambiguous_array[warm_cloud] = 1
        logging.debug('cold cloud only: %s', cloud.sum())
    else:
        cloud = cold_cloud | warm_cloud
        logging.debug('combined cloud: %s', cloud.sum())

    if cloud.sum() > 0:
        logging.debug('cold_cloud_pop: %s', cold_cloud_pop)
        logging.debug('desert_index: %s', desert_index)
        logging.debug('Mean temperature: %s', numpy.mean(
            thermal_array[cloud], dtype='float'))
        if ((cold_cloud_pop > pq_const.acca_cold_cloud_pop) and \
            (desert_index > pq_const.acca_desert_index) and \
            (numpy.mean(thermal_array[cloud], dtype='float') <
             pq_const.acca_cold_cloud_mean)):
            # Inititate 2nd Pass Testing
            r_cloud = acca_2nd_pass(cloud_mask=cloud,
                                    ambiguous_array=ambiguous_array,
                                    thermal_array=thermal_array,
                                    mean_cloud_temp=cold_cloud_mean,
                                    pq_const=pq_const, aux_data=aux_data)
            if r_cloud is None:
                return cloud
            return r_cloud

        elif ((desert_index <= pq_const.acca_desert_index) and
              (numpy.mean(thermal_array[cloud], dtype='float') <
               pq_const.acca_cold_cloud_mean)):
            return cold_cloud

        aux_data['acca_desert_index'] = 'failed'
        aux_data['acca_identified_pixels'] = 'all rejected'
        return numpy.zeros((dims[1], dims[2]), dtype='uint8')

    aux_data['acca_identified_pixels'] = 'all rejected'
    return numpy.zeros((dims[1], dims[2]), dtype='uint8')


def majority_filter(array, iterations=1):
    """
    Applies a majority filter to the input array.

    :param array:
        A 2D np array on which to perform majority filtering.

    :param iterations:
        The number of iterations to apply against the input array.
        Default is 1.

    :return:
        A 2D np array of type bool.
    """
    weights_array = [[1, 1, 1], [1, 1, 1], [1, 1, 1]]
    for _ in range(iterations):
        array = ndimage.convolve(array, weights_array)
        array = numexpr.evaluate("array > 4")
    return array


def calc_acca_cloud_mask(blue_dataset, green_dataset, red_dataset,
                         nir_dataset, swir1_dataset, swir2_dataset,
                         kelvin_array, pq_const, contiguity_mask,
                         aux_data=None):
    """
    Identifes the location of clouds.

    Follows the ACCA (Automated Cloud Cover Assessment) given by
    Irish, R et al, October 2006, Photogrammetric Engineering & Remote
    Sensing, Vol. 72, No. 10, October 2006, pp. 1179-1188.

    :param blue_dataset:
        A `NumPy` or `NumPy-like` dataset that allows indexing
        and returns a `NumPy` dataset containing the blue spectral
        data in reflectance units scaled from 0 to 10,000.

    :param green_dataset:
        A `NumPy` or `NumPy-like` dataset that allows indexing
        and returns a `NumPy` dataset containing the green spectral
        data in reflectance units scaled from 0 to 10,000.

    :param red_dataset:
        A `NumPy` or `NumPy-like` dataset that allows indexing
        and returns a `NumPy` dataset containing the red spectral
        data in reflectance units scaled from 0 to 10,000.

    :param nir_dataset:
        A `NumPy` or `NumPy-like` dataset that allows indexing
        and returns a `NumPy` dataset containing the nir spectral
        data in reflectance units scaled from 0 to 10,000.

    :param swri1_dataset:
        A `NumPy` or `NumPy-like` dataset that allows indexing
        and returns a `NumPy` dataset containing the swir1 spectral
        data in reflectance units scaled from 0 to 10,000.

    :param swri2_dataset:
        A `NumPy` or `NumPy-like` dataset that allows indexing
        and returns a `NumPy` dataset containing the swir2 spectral
        data in reflectance units scaled from 0 to 10,000.

    :param kelvin_array:
        A 2D Nump array containing temperature in degrees Kelvin.

    :param contiguity_mask:
        A 2D Numpy array where 0 = Non-contiguous and 1 = Contiguous.

    :param aux_data:
        A dict into which the function will place interesting intermediate
        results and metrics generated during processing - refer to the 
        code for details

    :return:
        A Boolean ndarray with 1 as non-cloud and 0 as cloud.
        Note: Any caller-supplied aux_data dict will be updated
    """

    start_time = datetime.datetime.now()

    aux_data = aux_data or {}  # set aux_data to empty dict if undefined
    dims = (6, kelvin_array.shape[0], kelvin_array.shape[1])

    # Contiguity masking
    null_nan_array = numpy.ones(contiguity_mask.shape, dtype=numpy.float32)
    if contiguity_mask is not None:
        null_nan_array[~contiguity_mask] = NAN

    # reflectance_stack contains surface reflectance in un-scaled units
    reflectance_stack = numpy.zeros(dims, dtype='float32')
    scaling_factor = numpy.float32(0.0001)
    variables = {'scaling_factor': scaling_factor,
                 'null_nan_array': null_nan_array}
    expr = "array * scaling_factor * null_nan_array"
    variables['array'] = blue_dataset
    reflectance_stack[0] = numexpr.evaluate(expr, variables)
    variables['array'] = green_dataset
    reflectance_stack[1] = numexpr.evaluate(expr, variables)
    variables['array'] = red_dataset
    reflectance_stack[2] = numexpr.evaluate(expr, variables)
    variables['array'] = nir_dataset
    reflectance_stack[3] = numexpr.evaluate(expr, variables)
    variables['array'] = swir1_dataset
    reflectance_stack[4] = numexpr.evaluate(expr, variables)
    variables['array'] = swir2_dataset
    reflectance_stack[5] = numexpr.evaluate(expr, variables)

    kelvin_array[~contiguity_mask] = NAN

    cloud = acca(reflectance_stack, kelvin_array, null_nan_array,
                 pq_const, aux_data=aux_data)

    # Apply filtering; gets rid of isolated pixels, and holes.
    if cloud.sum() > 0:
        # Majority filtering
        cloud = majority_filter(cloud, iterations=2)

    # Note this is percent of the array, not just contiguous areas.
    cloud_percent = (float(cloud.sum()) / cloud.size) * 100
    aux_data['final_cloud_layer_percent'] = cloud_percent

    # Calculate cloud percent to return. This is used for input into Fmask.
    # As Fmask calculates percentages of contiguous areas only, the return
    # value from ACCA wll be for contiguous areas if the argument null_mask
    # is set. Otherwise the percent of the entire array is returned.

    cloud_mask = ~cloud.astype(numpy.bool)

    # Upper cloud prob is 22.5
    # Also if low cloud % or desert region, set to original fmask probability
    # Following code is not used (results go out of scope!!) SMR 31/10/14

    # cld_pct = (float(cloud.sum())/contiguity_mask.sum()) * 100
    # if cld_pct > 22.5:
    #     cld_pct = 22.5
    # elif ((cld_pct < 0.03) | (desert_index < 0.5)):
    #     cld_pct = 22.5

    process_time = start_time - datetime.datetime.now()
    aux_data['acca_process_time_secs'] = process_time.total_seconds()
    return cloud_mask
