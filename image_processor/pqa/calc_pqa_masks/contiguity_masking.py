import os, numpy, logging
import numexpr

from osgeo import gdal
from scipy import ndimage

from ULA3.dataset import SceneDataset
from ULA3.utils import dump_array
from ULA3.common.pqa_result import PQAResult
from ULA3.image_processor import ProcessorConfig
from ULA3 import DataManager
from IDL_functions import IDL_Histogram

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

    def linear_percent(array, percent=2):
        """
        Image contrast enhancement.

        A 2D image is ehanced via a specifed percentage (Default 2%).

        :param array:
            A 2D Numpy array of any data type.

        :param perecent:
            A value in the range of 0-100. Default is 2.

        :return:
            A 2D array of the same dimensions as the input array, with values
            scaled by the specified percentage.

        :author:
            Josh Sixsmith, joshua.sixsmith@ga.gov.au
        """

        if len(array.shape) != 2:
            raise Exception('Only 2D arrays are supported.')

        if (percent <= 0) or (percent >= 100):
            raise Exception('Percent must be between 0 and 100')

        low  = (percent/100.)
        high = (1 - (percent/100.))
        nbins = 256.
        imgmin = numpy.min(array).astype('float')
        imgmax = numpy.max(array).astype('float')
        if array.dtype == 'uint8':
            hist, bedge = numpy.histogram(array, bins=nbins, range=(0,255))
            binsize = 1.
            imgmin = 0
            imgmax = 255
        else:
            hist, bedge = numpy.histogram(array, bins=nbins)
            binsize = (imgmax - imgmin)/(nbins - 1)

        cumu = numpy.cumsum(hist, dtype='float')
        n = cumu[-1]

        x1 = numpy.searchsorted(cumu, n * low)
        while cumu[x1] == cumu[x1 + 1]:
            x1 = x1 + 1

        x2 = numpy.searchsorted(cumu, n * high)
        while cumu[x2] == cumu[x2 - 1]:
            x2 = x2 - 1

        minDN = x1 * binsize + imgmin
        maxDN = x2 * binsize + imgmin

        # Scaling in the range 0-255.
        y1 = 0
        y2 = 255
        m = float(y2 - y1)/(maxDN - minDN)
        b = m*(-minDN)
        scl_img = array*m + b
        scl_img[scl_img > 255] = 255
        scl_img[scl_img < 0] = 0
        # Could floor the result before converting to uint8 ?
        scl_img = scl_img.astype('uint8')

        return scl_img


    def Contiguity(image_stack, satellite):
        """
        Determines locations of null values.

        Null values for every band are located in order to create band
        contiguity.

        :param image:
            An nD Numpy array of all bands (ordered).

        :param mask:
            Output array.

        :param slc_off:
            Whether to perform image scaling to flag potential
            saturated/non-contiguous pixels. Should only be applied to Non
            landsat 7 products. (L7 products prior to slc-off could be run).
            Default is False.

        :return:
            A single ndarray determining band/pixel contiguity. 1 for
            contiguous, 0 for non-contiguous.

        :notes:
            Attempts to flag thermal anomolies for Landsat 5TM as well.
        """

        if len(image_stack) == 0:
            return None

        assert type(image_stack[0]) == numpy.ndarray, 'Input is not valid'

        logger.debug('Determining pixel contiguity')
        # Create mask array with True for all pixels which are non-zero in all bands
        #mask = image_stack.all(0)
        mask = numexpr.evaluate('prod(image_stack, 0)') != 0

        # The following is only valid for Landsat 5 images
        if satellite.TAG == 'LS5':
            logger.debug('Finding thermal edge anomalies')
            # Apply thermal edge anomalies
            struct = numpy.ones((7,7), dtype='bool')
            erode  = ndimage.binary_erosion(mask, structure=struct)

            dims = mask.shape
            th_anom = numpy.zeros(dims, dtype='bool').flatten()

            pix_3buff_mask = mask - erode
            pix_3buff_mask[pix_3buff_mask > 0] = 1
            edge = pix_3buff_mask == 1

            low_sat = image_stack[5,:,:] == 1
            low_sat_buff = ndimage.binary_dilation(low_sat, structure=struct)

            s = [[1,1,1],[1,1,1],[1,1,1]]
            low_sat, num_labels = ndimage.label(low_sat_buff, structure=s)

            flat_label = low_sat.flatten()
            labels = low_sat[edge]
            ulabels = numpy.unique(labels[labels > 0])

            # Testing a new method, more code but might be quicker
            #find_lab = numpy.in1d(flat_label, ulabels)
            #th_anom |= find_lab

            # Histogram method, a lot faster
            mx = numpy.max(ulabels)
            h = IDL_Histogram(flat_label, min=0, max=mx, reverse_indices='ri')
            hist = h['histogram']
            ri = h['ri']

            for i in numpy.arange(ulabels.shape[0]):
                if hist[ulabels[i]] == 0:
                    continue
                th_anom[ri[ri[ulabels[i]]:ri[ulabels[i]+1]]] = True

            th_anom = ~(th_anom.reshape(dims))
            mask &= th_anom


        '''
        # Only to apply to Non-Landsat7 products.
        if slc_off == False:

            # The following test only applies to bands [1,2,3,4]. The anomolies are
            # already removed in bands [5,6,7] in the plain contiguity masking section.
            struct = numpy.ones((7,7), dtype='bool')
            erode  = ndimage.binary_erosion(mask, structure=struct) << bitpos

            # Testing just an erode of 3 pixels
            mask = erode

            #pix_3buff_mask = mask - erode
            #pix_3buff_mask[pix_3buff_mask > 0] = 1

            maxBit = 1 << bitpos

            for band in range(4):
                scl_img = linear_percent(l1t_stack[band,:,:], percent=3)
                scl_img *= pix_3buff_mask
                mask *= (scl_img < 255)


                m2 = scl_img == 255
                b = band + 1
                temp_name  = 'Band %d' %b  + '_scaled.tif'
                driver     = gdal.GetDriverByName("GTiff")
                outDataset = driver.Create(temp_name, scl_img.shape[1], scl_img.shape[0], 1, 1)

                outBand = outDataset.GetRasterBand(1)
                outBand.WriteArray(m2)
                #outDataset.SetGeoTransform(geoc)
                #outDataset.SetProjection(proj)
                outDataset = None
                driver = None

        '''

        return mask

    mask = Contiguity(l1t_stack, l1t_input_dataset.satellite)

    bit_index = CONFIG.pqa_test_index['CONTIGUITY']
    result.set_mask(mask, bit_index)
    if CONFIG.debug:
        dump_array(mask,
                   os.path.join(CONFIG.work_path, 'mask_%02d.tif' % bit_index),
                   l1t_input_dataset)
