#! /usr/bin/env python

import datetime, numpy, gc, logging, os, osr, numexpr
from numpy.core.numeric import nan
from osgeo import gdal
from scipy import ndimage
from ULA3.dataset import SceneDataset
from ULA3.utils import dump_array
from ULA3.common.pqa_result import PQAResult
from ULA3.image_processor import ProcessorConfig
from ULA3.image_processor import constants
from ULA3 import DataManager, DataGrid
from ULA3.utils import log_multiline

logger = logging.getLogger('root.' + __name__)

NaN = numpy.float32(nan)

def process(subprocess_list=[], resume=False):
    logger.info('%s.process(%s, %s) called', __name__, subprocess_list, resume)

    CONFIG = ProcessorConfig()
    DATA = DataManager()

    nbar_input_dataset = DATA.get_item(CONFIG.input['nbar']['path'], SceneDataset)
    assert nbar_input_dataset, 'Unable to retrieve SceneDataset object for NBAR input scene dataset'
    logger.debug( 'SceneDataset object for %s retrieved', nbar_input_dataset.pathname)

    nbar_stack = DATA.get_item('nbar_stack', numpy.ndarray)
    assert nbar_stack is not None, 'Unable to retrieve ndarray object for nbar_stack'
    logger.debug( 'ndarray object for nbar_stack retrieved')

    result = DATA.get_item('result.tif', PQAResult)
    assert result, 'Unable to retrieve PQAResult object for result'
    logger.debug( 'PQAResult object for result retrieved')

    pqa_temp_output = DATA.get_item('pqa_temp_output.dat', str)
    assert pqa_temp_output, 'Unable to retrieve string object for pqa_temp_output'
    logger.debug( 'string object for pqa_temp_output retrieved')

    # Get the PQ constants
    pq_const = constants.pqaContants(nbar_input_dataset.sensor)

    #===========================================================================
    # contiguity_mask = DATA.get_item('contiguity_mask', numpy.ndarray)
    # assert contiguity_mask is not None, 'Unable to retrieve ndarray object for contiguity_mask'
    # logger.debug( 'ndarray object for contiguity_mask retrieved')
    #===========================================================================
    #assert CONFIG.pqa_test_index['CONTIGUITY'] in result.test_set, 'Contiguity test not yet run'
    #contiguity_mask = (result.array & (1 << CONFIG.pqa_test_index['CONTIGUITY'])) > 0
    assert pq_const.contiguity in result.test_set, 'Contiguity test not yet run'
    contiguity_mask = (result.array & (1 << pq_const.contiguity)) > 0

    kelvin_grid = DATA.get_item('kelvin.tif', DataGrid)
    assert kelvin_grid is not None, 'Unable to retrieve DataGrid object for kelvin_grid'
    logger.debug( 'DataGrid object for kelvin_grid retrieved')
    kelvin_array = kelvin_grid.array

    def CloudMask(image_stack, kelvin_array, contiguity_mask=None):
        """
        Identifes the location of clouds.

        Follows the ACCA (Automated Cloud Cover Assessment) given by
        Irish, R et al, October 2006, Photogrammetric Engineering & Remote
        Sensing, Vol. 72, No. 10, October 2006, pp. 1179-1188.

        :param image_stack:
            An ordered array of all bands in reflectance units scaled from 0
            to 10,000.

        :param kelvin_array:
            A 2D Nump array containing temperature in degrees Kelvin.

        :param contiguity_mask:
            A 2D Numpy array where 0 = Non-contiguous and 1 = Contiguous.

        :return:
            A Boolean ndarray with 1 as non-cloud and 0 as cloud.
        """
        if len(image_stack) == 0: return None
        assert type(image_stack[0]) == numpy.ndarray, 'Image input is not valid'


    #-------------------------Functions/Filters----------------------------------
        def water_test(reflectance_stack):
            """
            The WATER_TEST function is used to identify water pixels.

            :param reflectance_stack:
                An nD Numpy array of all bands (ordered).

            :return:
                An ndarray of type bool for water areas where True == water.
            """

            # water = (ndvi < 0.1) or ((nir < 0.04) and (swir1 < 0.05))
            #water = numpy.logical_or(ndvi < -01, numpy.logical_or(array[3,:,:]
            #                             < 0.04, array[4,:,:] < 0.05))

            return numexpr.evaluate("(NDVI < 0.1) & (band5_array < 0.05)",
                                     {'NDVI': ndvi(reflectance_stack), 'band5_array': reflectance_stack[4]}, locals())


        def ndvi(reflectance_stack):
            """
            The NDVI function calculates the Normalised Differenced Vegetation Index.

            :param reflectance_stack:
                An nD Numpy array of all bands (ordered).

            :return:
                A Numpy nD array of NDVI values in the range of -1 to 1.
            """

            return numexpr.evaluate("(a - b) / (a + b)", {'a': reflectance_stack[3], 'b': reflectance_stack[2]}, locals())

        # Calculate NSDI; FILTER 2
        def ndsi(reflectance_stack):
            """
            Calculates the Normalised Snow Difference Index.

            :param reflectance_stack:
                An nD Numpy array of all bands (ordered).

            :return:
                An Numpy nD array of values in the range of -1 to 1.
            """

            # ndsi = (green - swir1)/(green + swir1)
            return numexpr.evaluate("(a - b) / (a + b)", {'a': reflectance_stack[1], 'b': reflectance_stack[4]}, locals())

        # Calculate band 5/6 composite; FILTER 4
        def filter4(reflectance_stack):
            """
            Calculates band 5/6 composite.

            Clouds are generally cold and highly reflective in band 5.

            :param reflectance_stack:
                An nD Numpy array of all bands (ordered).

            :returns:
                An nD Numpy  array containing the band 5/6 composite.
            """

            return numexpr.evaluate("(1.0 - swir1) * TM6", {'swir1': reflectance_stack[4], 'TM6': thermal_array}, locals())

        # Calculate band 4/3 ratio; FILTER 5/Simple Ratio
        def filter5(reflectance_stack):
            """
            Calculates band 4/3 ratio.

            Eliminates highly reflective vegetation.

            :param reflectance_stack:
                An nD Numpy array of all bands (ordered).

            :return:
                An nD Numpy array (unbound range).
            """

            return numexpr.evaluate("nir / red", {'nir': reflectance_stack[3], 'red': reflectance_stack[2]}, locals())

        # Calculate band 4/2 ratio; FILTER 6
        def filter6(reflectance_stack):
            """
            Calculates band 4/2 ratio.

            Eliminates senescing vegetaion.

            :param reflectance_stack:
                An nD Numpy array of all bands (ordered).

            :returns:
                An nD Numpy array (unbound range).
            """

            #return nir/green
            #return  numpy.divide(array[3],array[1])
            return numexpr.evaluate("nir / green", {'nir': reflectance_stack[3], 'green': reflectance_stack[1]}, locals())

        # Calculate band 4/5 ratio; FILTER 7
        def filter7(reflectance_stack):
            """
            Calculates band 4/5 ratio.

            Eliminates highly reflective soils.

            :param reflectance_stack:
                An nD Numpy array of all bands (ordered).

            :return:
                An nD Numpy array (unbound range).
            """

            return numexpr.evaluate("nir / swir1", {'nir': reflectance_stack[3], 'swir1': reflectance_stack[4]}, locals())

        # Calculate skewness. Couldn't find a numpy alternative
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

            cubed_deviates = (cloud_thermal_array - mean_temp)**3
            #sum_cubed_dv   = cubed_deviates.sum()
            sum_cubed_dv   = numpy.sum(cubed_deviates, dtype='float64')
            cubed_stdv     = stdv_temp**3
            return (sum_cubed_dv/cubed_stdv)/count

        # ACCA Second Pass
        def acca_2nd_pass(cloud_mask, ambiguous_array, thermal_array, mean_cloud_temp):
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

            :return:
                Depending on the result of the second pass, the acca_second_pass
                function can return None (conditions not, therefore
                stick with the cloud identified in the first pass), or a cloud
                mask with with 1 as cloud and 0 as not cloud.

            :author:
                Josh Sixsmith, joshua.sixsmith@ga.gov.au
            """

            global acca_logfile
            logger.info('ACCA Pass Two Engaged\n')
            acca_logfile.write('Pass Two Engaged\n')

            #cloud_stddev    = thermal_array[cloud_mask].std()
            cloud_stddev    = numpy.std(thermal_array[cloud_mask], dtype='float64', ddof=1)
            #cloud_count     = cloud_mask.sum() # Sum is used as valid pixels = 1
            cloud_count     = numpy.sum(cloud_mask, dtype='float64')

            acca_logfile.write('Standard Deviation: %f\n' %cloud_stddev)

            # Histogram Percentiles for new thermal thresholds
            upper     = numpy.percentile(thermal_array[cloud_mask], 97.5)
            lower     = numpy.percentile(thermal_array[cloud_mask], 83.5)
            upper_max = numpy.percentile(thermal_array[cloud_mask], 98.75)


            acca_logfile.write('97.5 percentile: %f\n' %upper)
            acca_logfile.write('83.5 percentile: %f\n' %lower)
            acca_logfile.write('98.75 percentile: %f\n' %upper_max)

            # Test for negative skewness
            skew = skewness(thermal_array[cloud_mask], mean_temp=mean_cloud_temp, stdv_temp=cloud_stddev,
                            count=cloud_count)
            logger.debug('skew: %s', skew)

            acca_logfile.write('Skewness: %f\n' %skew)

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
                else :
                    new_upper = upper_max


                #query  = (thermal_array > new_lower) & (thermal_array <= new_upper)
                query = numexpr.evaluate("((ambiguous_array * thermal_array) > new_lower) & ((ambiguous_array * thermal_array) <= new_upper)")
                #query2 = ((ambiguous_array * thermal_array) != 0) & ((ambiguous_array * thermal_array) <= new_lower)
                query2 = numexpr.evaluate("((ambiguous_array * thermal_array) != 0) & ((ambiguous_array * thermal_array) <= new_lower)")

                # Compute stats for each query/class
                # Max
                if query.sum() == 0:
                    qmax = 295
                else:
                    qmax  = thermal_array[query].max()
                if query2.sum() == 0:
                    qmax2 = 295
                else:
                    qmax2 = thermal_array[query2].max()

                acca_logfile.write('Class 1 max: %f\n' %qmax)
                acca_logfile.write('Class 2 max: %f\n' %qmax2)


                # Mean
                if query.sum() == 0:
                    qmean = 295
                else:
                    qmean  = thermal_array[query].mean()
                if query2.sum() == 0:
                    qmean2 = 295
                else:
                    qmean2 = thermal_array[query2].mean()

                acca_logfile.write('Class 1 mean: %f\n' %qmean)
                acca_logfile.write('Class 2 mean: %f\n' %qmean2)

                # Class percentage of scene
                qpop  = (float(query.sum())/ambiguous_array.size)*100
                qpop2 = (float(query2.sum())/ambiguous_array.size)*100

                acca_logfile.write('Class 1 percent: %f\n' %qpop)
                acca_logfile.write('Class 2 percent: %f\n' %qpop2)



                if qpop < pq_const.acca_thermal_effect: #CONFIG.pqa_param['acca_thermal_effect']:
                    if qmean < pq_const.acca_coldCloud_mean: #CONFIG.pqa_param['acca_coldcloud_mean']:
                        # Combine all cloud classes
                        return numexpr.evaluate("cloud_mask | query | query2")
                    elif qpop2 < pq_const.acca_thermal_effect: #CONFIG.pqa_param['acca_thermal_effect']:
                        if qmean2 < pq_const.acca_coldCloud_mean: #CONFIG.pqa_param['acca_coldcloud_mean']:
                            # Combine lower threshold clouds and pass 1 clouds
                            return numexpr.evaluate("cloud_mask | query2")
                    else: # Keep first pass cloud
                        return None
                else: # Keep fist pass cloud
                    return None

            else:
                #query  = ((ambiguous_array * thermal_array) > lower) * ((ambiguous_array * thermal_array) <= upper)
                #query2 = ((ambiguous_array * thermal_array) != 0) * ((ambiguous_array * thermal_array) <= lower)
                query = numexpr.evaluate("((ambiguous_array * thermal_array) > lower) & ((ambiguous_array * thermal_array) <= upper)")
                query2 = numexpr.evaluate("((ambiguous_array * thermal_array) != 0) & ((ambiguous_array * thermal_array) <= lower)")


                # Compute stats for each query/class
                # Max
                if query.sum() == 0:
                    qmax = 295
                else:
                    qmax  = thermal_array[query].max()
                qmax2 = thermal_array[query2].max()

                acca_logfile.write('Class 1 max: %f\n' %qmax)
                acca_logfile.write('Class 2 max: %f\n' %qmax2)

                # Mean
                if query.sum() == 0:
                    qmean = 295
                else:
                    qmean  = thermal_array[query].mean()
                qmean2 = thermal_array[query2].mean()

                acca_logfile.write('Class 1 mean: %f\n' %qmean)
                acca_logfile.write('Class 2 mean: %f\n' %qmean2)


                # Class percentage of scene
                qpop  = (float(query.sum())/ambiguous_array.size)*100
                qpop2 = (float(query2.sum())/ambiguous_array.size)*100

                acca_logfile.write('Class 1 percent: %f\n' %qpop)
                acca_logfile.write('Class 2 percent: %f\n' %qpop2)


                if qpop < pq_const.acca_thermal_effect: #CONFIG.pqa_param['acca_thermal_effect']:
                    if qmean < pq_const.acca_coldCloud_mean: #CONFIG.pqa_param['acca_coldcloud_mean']:
                        # Combine all cloud classes
                        return numexpr.evaluate("cloud_mask | query | query2")
                    elif qpop2 < pq_const.acca_thermal_effect: #CONFIG.pqa_param['acca_thermal_effect']:
                        if qmean2 < pq_const.acca_coldCloud_mean: #CONFIG.pqa_param['acca_coldcloud_mean']:
                            # Combine lower threshold clouds and pass 1 clouds
                            return numexpr.evaluate("cloud_mask | query2")
                    else: # Keep first pass cloud
                        return None
                else: # Keep fist pass cloud
                    return None



        def acca(reflectance_stack, thermal_array, potential_cloud_array):
            """
            The first pass processing of the ACCA algorithm.

            :param reflectance_stack:
                An nD numpy array containing the reflectance values for each
                band.

            :param thermal_array:
                A 2D Numpy array containing the thermal band.

            :param potential_cloud_array:
                A 2D Numpy array containing values of 1 for valid areas and NaN for invalid areas.

            :return:
                An 2D Numpy array cloud mask, with True for non-cloud and False
                for cloud.

            :author:
                Josh Sixsmith, joshua.sixsmith@ga.gov.au
            """

            global acca_logfile
            global Desert_Index

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

            # Create the array for Ambiguous Pixels - NaN means not ambiguous
            ambiguous_array = numpy.ones(potential_cloud_array.shape, dtype=numpy.float32) * NaN

            # Keep an un-NaNed copy of the thermal band for later use
            #thermal_array = image_stack[5].copy()
            #thermal_array = image_stack[5]

            # Will add in a water mask, to remove cold water bodies that have been
            # put into the ambigous group. If the water body is high in red
            # reflectance, it will have made it this far.

            water_mask = water_test(reflectance_stack)
#            logger.debug('reflectance_stack[:,1612,2126] %s', reflectance_stack[:,1612,2126])
#            logger.debug('reflectance_stack[:,2252,4294] %s', reflectance_stack[:,2252,4294])
            # Filter 1; brightness threshold (remove dark targets)
            #query = numexpr.evaluate("where((potential_cloud_array * b3) > thresh_f1, 1, NaN)",
            #                         {'b3': reflectance_stack[2],
            #                          'thresh_f1' : CONFIG.pqa_param['acca_thresh_f1'],
            #                          'NaN' : NaN},
            #                         locals())
            query = numexpr.evaluate("where((potential_cloud_array * b3) > thresh_f1, 1, NaN)",
                                     {'b3': reflectance_stack[2],
                                      'thresh_f1' : pq_const.acca_thresh_f1,
                                      'NaN' : NaN},
                                     locals())

            potential_cloud_array *= query

            if CONFIG.debug:
                dump_array(query, os.path.join(CONFIG.work_path, 'acca_filter1.tif'), template_dataset=nbar_input_dataset, dtype=gdal.GDT_Float32)
                dump_array(potential_cloud_array, os.path.join(CONFIG.work_path, 'potential_cloud_array1.tif'), template_dataset=nbar_input_dataset, dtype=gdal.GDT_Float32)

            # Filter 2: NDSI
            ndsi_array = ndsi(reflectance_stack)
            #query   = numexpr.evaluate("where((potential_cloud_array * ndsi_array) < thresh_f2, 1, NaN)",{'thresh_f2' : CONFIG.pqa_param['acca_thresh_f2'], 'NaN' : NaN}, locals())
            query   = numexpr.evaluate("where((potential_cloud_array * ndsi_array) < thresh_f2, 1, NaN)",{'thresh_f2' : pq_const.acca_thresh_f2, 'NaN' : NaN}, locals())

            # Find the snow pixels. Sum is used to find the total cloud pixels
            # as valid pixels = 1.  Sum of ones therefore = count
            find         = numexpr.evaluate("where(ndsi_array >= thresh_f2, 1, 0)",{'thresh_f2' : pq_const.acca_thresh_f2}, locals())
            snow_pixels  = find.sum() # Sum is used as valid pixels = 1
            snow_percent = (float(snow_pixels)/NaN.size) * 100

            potential_cloud_array *= query
            if CONFIG.debug:
                dump_array(query, os.path.join(CONFIG.work_path, 'acca_filter2.tif'), template_dataset=nbar_input_dataset, dtype=gdal.GDT_Float32)
                dump_array(potential_cloud_array, os.path.join(CONFIG.work_path, 'potential_cloud_array2.tif'), template_dataset=nbar_input_dataset, dtype=gdal.GDT_Float32)

            acca_logfile.write('Snow Percent: %f\n' %snow_percent)


            # Filter 3; Temp. threshold
            #query   = numexpr.evaluate("where((potential_cloud_array * thermal_array) < thresh_f3, 1, NaN)",
            #                           {'thresh_f3' : CONFIG.pqa_param['acca_thresh_f3'],
            #                            'NaN' : NaN},
            #                           locals())
            query   = numexpr.evaluate("where((potential_cloud_array * thermal_array) < thresh_f3, 1, NaN)",
                                       {'thresh_f3' : pq_const.acca_thresh_f3,
                                        'NaN' : NaN},
                                       locals())
            potential_cloud_array *= query
            if CONFIG.debug:
                dump_array(query, os.path.join(CONFIG.work_path, 'acca_filter3.tif'), template_dataset=nbar_input_dataset, dtype=gdal.GDT_Float32)
                dump_array(potential_cloud_array, os.path.join(CONFIG.work_path, 'potential_cloud_array3.tif'), template_dataset=nbar_input_dataset, dtype=gdal.GDT_Float32)


            # Filter 4; Band 5/6 composite
            _temporary = potential_cloud_array * filter4(reflectance_stack)
            #query   = numexpr.evaluate("where(_temporary < thresh_f4, 1, NaN)",
            #                           {'thresh_f4' : CONFIG.pqa_param['acca_thresh_f4'],
            #                            'NaN' : NaN},
            #                           locals())
            query   = numexpr.evaluate("where(_temporary < thresh_f4, 1, NaN)",
                                       {'thresh_f4' : pq_const.acca_thresh_f4,
                                        'NaN' : NaN},
                                       locals())

            # Get ambiguous pixels
            #find  = numexpr.evaluate("_temporary >= thresh_f4",{'thresh_f4' : CONFIG.pqa_param['acca_thresh_f4']}, locals())
            find  = numexpr.evaluate("_temporary >= thresh_f4",{'thresh_f4' : pq_const.acca_thresh_f4}, locals())
            ambiguous_array[find] = 1

            potential_cloud_array *= query
            if CONFIG.debug:
                dump_array(query, os.path.join(CONFIG.work_path, 'acca_filter4.tif'), template_dataset=nbar_input_dataset, dtype=gdal.GDT_Float32)
                dump_array(potential_cloud_array, os.path.join(CONFIG.work_path, 'potential_cloud_array4.tif'), template_dataset=nbar_input_dataset, dtype=gdal.GDT_Float32)

            # Will add in a water mask, to remove cold water bodies that have been
            # put into the ambiguous group. If the water body is high in red
            # reflectance, it will have made it this far.

            ambiguous_array[water_mask] = NaN # All water is unambiguous
            del water_mask; gc.collect

            # Filter 5; Band 4/3 ratio (Simple veg ratio)
            _temporary = potential_cloud_array * filter5(reflectance_stack)
            #query   = numexpr.evaluate("where(_temporary < thresh_f5, 1, NaN)",{'thresh_f5' : CONFIG.pqa_param['acca_thresh_f5'], 'NaN' : NaN}, locals())
            query   = numexpr.evaluate("where(_temporary < thresh_f5, 1, NaN)",{'thresh_f5' : pq_const.acca_thresh_f5, 'NaN' : NaN}, locals())

            # Get ambiguous pixels
            #find  = numexpr.evaluate("_temporary >= thresh_f5",{'thresh_f5' : CONFIG.pqa_param['acca_thresh_f5']}, locals())
            find  = numexpr.evaluate("_temporary >= thresh_f5",{'thresh_f5' : pq_const.acca_thresh_f5}, locals())
            ambiguous_array[find] = 1

            potential_cloud_array *= query
            if CONFIG.debug:
                dump_array(query, os.path.join(CONFIG.work_path, 'acca_filter5.tif'), template_dataset=nbar_input_dataset, dtype=gdal.GDT_Float32)
                dump_array(potential_cloud_array, os.path.join(CONFIG.work_path, 'potential_cloud_array5.tif'), template_dataset=nbar_input_dataset, dtype=gdal.GDT_Float32)

            # Filter 6; Band 4/2 ratio (Dying/senescing veg)
            _temporary = potential_cloud_array * filter6(reflectance_stack)
            #query   = numexpr.evaluate("where(_temporary < thresh_f6, 1, NaN)",{'thresh_f6' : CONFIG.pqa_param['acca_thresh_f6'], 'NaN' : NaN}, locals())
            query   = numexpr.evaluate("where(_temporary < thresh_f6, 1, NaN)",{'thresh_f6' : pq_const.acca_thresh_f6, 'NaN' : NaN}, locals())

            # Tally filter 6 survivors
            f6_surv = numpy.nansum(query)

            # Get ambiguous pixels
            #find  = numexpr.evaluate("_temporary >= thresh_f6",{'thresh_f6' : CONFIG.pqa_param['acca_thresh_f6'], 'NaN' : NaN}, locals())
            find  = numexpr.evaluate("_temporary >= thresh_f6",{'thresh_f6' : pq_const.acca_thresh_f6, 'NaN' : NaN}, locals())
            ambiguous_array[find] = 1

            potential_cloud_array *= query
            if CONFIG.debug:
                dump_array(query, os.path.join(CONFIG.work_path, 'acca_filter6.tif'), template_dataset=nbar_input_dataset, dtype=gdal.GDT_Float32)
                dump_array(potential_cloud_array, os.path.join(CONFIG.work_path, 'potential_cloud_array6.tif'), template_dataset=nbar_input_dataset, dtype=gdal.GDT_Float32)

            # Filter 7; Band 4/5 ratio (Identify highly reflective soils/rocks)
            # The results of this query are clouds at first pass
            _temporary = potential_cloud_array * filter7(reflectance_stack)

            #query = numexpr.evaluate("where(_temporary > thresh_f7, 1, NaN)",{'thresh_f7' : CONFIG.pqa_param['acca_thresh_f7'], 'NaN' : NaN}, locals())
            query = numexpr.evaluate("where(_temporary > thresh_f7, 1, NaN)",{'thresh_f7' : pq_const.acca_thresh_f7, 'NaN' : NaN}, locals())

            # Tally filter 7 survivors
            f7_surv      = numpy.nansum(query)
            Desert_Index = float(f7_surv)/f6_surv

            acca_logfile.write('Desert Index: %f\n' %Desert_Index)

            # Get ambiguous pixels
            #find  = numexpr.evaluate("_temporary <= thresh_f7",{'thresh_f7' : CONFIG.pqa_param['acca_thresh_f7']}, locals())
            find  = numexpr.evaluate("_temporary <= thresh_f7",{'thresh_f7' : pq_const.acca_thresh_f7}, locals())
            ambiguous_array[find] = 1

            potential_cloud_array *= query
            if CONFIG.debug:
                dump_array(query, os.path.join(CONFIG.work_path, 'acca_filter7.tif'), template_dataset=nbar_input_dataset, dtype=gdal.GDT_Float32)
                dump_array(potential_cloud_array, os.path.join(CONFIG.work_path, 'potential_cloud_array7.tif'), template_dataset=nbar_input_dataset, dtype=gdal.GDT_Float32)

            # Filter 8; Band 5/6 composite (Separate warm/cold clouds)
            _temporary = potential_cloud_array * filter4(reflectance_stack)
            #cold_cloud    = numexpr.evaluate("_temporary < thresh_f8",{'thresh_f8' : CONFIG.pqa_param['acca_thresh_f8']}, locals())
            cold_cloud    = numexpr.evaluate("_temporary < thresh_f8",{'thresh_f8' : pq_const.acca_thresh_f8}, locals())
            #warm_cloud    = numexpr.evaluate("_temporary >= thresh_f8",{'thresh_f8' : CONFIG.pqa_param['acca_thresh_f8']}, locals())
            warm_cloud    = numexpr.evaluate("_temporary >= thresh_f8",{'thresh_f8' : pq_const.acca_thresh_f8}, locals())

            cold_cloud_pop  = (float(cold_cloud.sum())/ambiguous_array.size) * 100
            cold_cloud_mean =numpy.mean(thermal_array[cold_cloud], dtype='float64')
            warm_cloud_pop  = (float(warm_cloud.sum())/ambiguous_array.size) * 100
            warm_cloud_mean = numpy.mean(thermal_array[warm_cloud], dtype='float64')

            acca_logfile.write('Cold Cloud Percent: %f\n' %cold_cloud_pop)
            acca_logfile.write('Cold Cloud Mean: %f\n' %cold_cloud_mean)
            acca_logfile.write('Warm Cloud Percent: %f\n' %warm_cloud_pop)
            acca_logfile.write('Warm Cloud Mean: %f\n' %warm_cloud_mean)

            del query, find, _temporary; gc.collect()



        #"""
        #             Tests for snow and desert.  If the thresholds aren't
        #             breached, Pass two is implemented.
        #"""

            # REDO of tests for pass two engagement
            #if Desert_Index <= CONFIG.pqa_param['acca_desertindex'] and snow_percent > CONFIG.pqa_param['acca_snow_threshold']:
            if Desert_Index <= pq_const.acca_desertIndex and snow_percent > pq_const.acca_snow_threshold:
                cloud = cold_cloud
                ambiguous_array[warm_cloud] = 1
                logger.debug('cold cloud only: %s', cloud.sum())
            else:
                cloud = cold_cloud | warm_cloud
                logger.debug('combined cloud: %s', cloud.sum())

            if cloud.sum() > 0:
                logger.debug('cold_cloud_pop: %s', cold_cloud_pop)
                logger.debug('Desert_Index: %s', Desert_Index)
                logger.debug('Mean temperature: %s', numpy.mean(thermal_array[cloud], dtype='float'))
                #if ((cold_cloud_pop > CONFIG.pqa_param['acca_coldcloud_pop']) and (Desert_Index > CONFIG.pqa_param['acca_desertindex'])
                #         and (numpy.mean(thermal_array[cloud], dtype='float') < CONFIG.pqa_param['acca_coldcloud_mean'])):
                if ((cold_cloud_pop > pq_const.acca_coldCloud_pop) and (Desert_Index > pq_const.acca_desertIndex)
                         and (numpy.mean(thermal_array[cloud], dtype='float') < pq_const.acca_coldCloud_mean)):
                    # Inititate 2nd Pass Testing
                    r_cloud = acca_2nd_pass(cloud_mask=cloud, ambiguous_array=ambiguous_array, thermal_array=thermal_array,
                                  mean_cloud_temp=cold_cloud_mean)
                    if r_cloud == None:
                        return cloud
                    else:
                        return r_cloud
                #elif ((Desert_Index <= CONFIG.pqa_param['acca_desertindex']) and (numpy.mean(thermal_array[cloud], dtype='float') <
                #                CONFIG.pqa_param['acca_coldcloud_mean'])):
                elif ((Desert_Index <= pq_const.acca_desertIndex) and (numpy.mean(thermal_array[cloud], dtype='float') <
                                pq_const.acca_coldCloud_mean)):
                    return cold_cloud
                else:
                    acca_logfile.write('Desert Index and Temperature Test failed.\n')
                    acca_logfile.write('All identified pixels rejected.\n')
                    return numpy.zeros((dims[1], dims[2]), dtype='uint8')
            else:
                acca_logfile.write('All identified pixels rejected.\n')
                return numpy.zeros((dims[1], dims[2]), dtype='uint8')

        def majority_filter(cloud):
            weights_array = [[1,1,1],[1,1,1],[1,1,1]]
            cloud = ndimage.convolve(cloud, weights_array)
            cloud = numexpr.evaluate("cloud > 4")
            cloud = ndimage.convolve(cloud, weights_array)
            cloud = numexpr.evaluate("cloud > 4")

    #----------------------------Processing Here---------------------------------
        global acca_logfile
        global Desert_Index

        start_time = datetime.datetime.now()
        acca_logfile = open(os.path.join(pqa_temp_output, 'scene01', 'ACCA_LOGFILE.txt'), 'w')

        dim = image_stack.shape

        # Contiguity masking
        null_nan_array = numpy.ones(contiguity_mask.shape, dtype=numpy.float32)
        if contiguity_mask is not None:
            null_nan_array[~contiguity_mask] = NaN
        potential_cloud_array = null_nan_array.copy()
#        if CONFIG.debug:
#            dump_array(null_nan_array, 'null_nan_array.tif', template_dataset=nbar_input_dataset, dtype=gdal.GDT_Float32)

        # reflectance_stack contains surface reflectance in un-scaled units
        scaling_factor = numpy.float32(0.0001)
        reflectance_stack = image_stack.astype(numpy.float32)
#        reflectance_stack *= (scaling_factor * null_nan_array)
        reflectance_stack = numexpr.evaluate("reflectance_stack * scaling_factor * null_nan_array")
        if CONFIG.debug:
            for band_index in range(dim[0]):
                reflectance_stack[band_index].dump(os.path.join(CONFIG.work_path, 'reflectance_stack_' + str(band_index)))

        thermal_array = numexpr.evaluate("kelvin_array * null_nan_array")
        if CONFIG.debug:
            dump_array(thermal_array, os.path.join(CONFIG.work_path, 'thermal_array.tif'), template_dataset=nbar_input_dataset, dtype=gdal.GDT_Float32)

        cloud = acca(reflectance_stack, thermal_array, potential_cloud_array)
        # Apply filtering; gets rid of isolated pixels, and holes.
        if cloud.sum() > 0:
            # Majority filtering
            majority_filter(cloud)

        # Note this is percent of the array, not just contiguous areas.
        Cloud_Percent = (float(cloud.sum())/cloud.size) * 100
        acca_logfile.write('Final Cloud Layer Percent: %f\n' %Cloud_Percent)

        # Calculate cloud percent to return. This is used for input into Fmask.
        # As Fmask calculates percentages of contiguous areas only, the return
        # value from ACCA wll be for contiguous areas if the argument null_mask
        # is set. Otherwise the percent of the entire array is returned.

        cld_pct = (float(cloud.sum())/contiguity_mask.sum()) * 100
        cloud_mask = ~cloud.astype(numpy.bool)

        # Upper cloud prob is 22.5
        # Also if low cloud % or desert region, set to original fmask probability
        if cld_pct > 22.5:
            cld_pct = 22.5
        elif ((cld_pct < 0.03) | (Desert_Index < 0.5)):
            cld_pct = 22.5

        end_time   = datetime.datetime.now()
        time       = end_time - start_time
        acca_logfile.write('ACCA Process Time: %s\n' %time)

        acca_logfile.close()
        return cloud_mask


    mask = CloudMask(nbar_stack, kelvin_array, contiguity_mask)
    
    log_multiline(logger.debug, mask, 'mask', '\t')

    #bit_index = CONFIG.pqa_test_index['ACCA']
    bit_index = pq_const.acca
    result.set_mask(mask, bit_index)
    if CONFIG.debug:
        dump_array(mask,
                   os.path.join(CONFIG.work_path, 'mask_%02d.tif' % bit_index),
                   nbar_input_dataset)
