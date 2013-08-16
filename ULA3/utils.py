import os, logging, unicodedata, subprocess, math, errno, re
import numpy, gdal, gdalconst
import numpy.ma as ma
from pprint import pformat

# note that some of these imports are not used directly in this file,
# but they are imported from this file so must remain.
from meta import print_call
from filter import read_array_int8, read_array_int16, read_array_int32, read_array_float32
from _execute import execute # needs to be imported otherwise we get circular dependency in _gdal_tools
from _gdal_tools import NUMPY_GDAL_TYPE_MAP as DTYPE_MAP
from _gdal_tools import warp
from _gdal_tools import Buffers, ImageShape
from _gdal_tools import default_bounds_getter as get_bounds

logger = logging.getLogger('root.' + __name__)





def log_multiline(log_function, log_text, title=None, prefix=''):
    """Function to log multi-line text
    """
    logger.debug('log_multiline(%s, %s, %s, %s) called', log_function, repr(log_text), repr(title), repr(prefix))

    if type(log_text) == str:
        logger.debug('log_text is type str')
        log_list = log_text.splitlines()
    elif type(log_text) == list and type(log_text[0]) == str:
        logger.debug('log_text is type list with first element of type text')
        log_list = log_text
    else:
        logger.debug('log_text is type ' + type(log_text).__name__)
        log_list = pformat(log_text).splitlines()

    log_function(prefix + '=' * 80)
    if title:
        log_function(prefix + title)
        log_function(prefix + '-' * 80)

    for line in log_list:
        log_function(prefix + line)

    log_function(prefix + '=' * 80)





def unicode_to_ascii(instring):
    """Convert unicode to char string if required and strip any leading/trailing whitespaces
    ToDO: Investigate whether we can just change the encoding of the DOM tree
    """
    result = instring
    if type(result) == unicode:
        result = unicodedata.normalize('NFKD', result).encode('ascii','ignore').strip(""" "'\n\t""")
    return result





@print_call(logger.info)
def generate_md5(directory):
    """
    Use the system command ``md5sum`` to generate :file:`md5sum.txt` for all files in a directory. The output appears
    within the specified directory.

    :param directory:
        The directory to produce the md5 sums for.
    """
    command_string = "rm -f md5sum.txt; for name in `find . -type f | grep -v 'md5sum.txt' | sort`; do md5sum $name 2>/dev/null; done > md5sum.txt"
    logger.info('Invoking: %s', command_string)
    result = execute(command_string=command_string, cwd=directory)

    if result['stdout']:
        log_multiline(logger.info, result['stdout'], command_string)

    if result['returncode']:
        log_multiline(logger.error, result['stderr'], 'stderr from ' + command_string, '\t')
        raise Exception('%s failed', command_string)

    return result;





def find_files(root_dir, filename_pattern='.*', case_insensitive = True):
    """
    List files matching a specified pattern (regex) anywhere under a root directory

    :param root_dir:
        Root directory to search within.

    :param filename_pattern:
        RegEx string for finding files.

    :param case_insensitive:
        Flag specifying whether name matching should be case insensitive.

    :return:
        List containing absolute pathnames of found files or empty list if none found.
    """
    if case_insensitive:
        file_regex = re.compile(filename_pattern, re.IGNORECASE)
    else:
        file_regex = re.compile(filename_pattern)

    filename_list = []

    for root, _dirs, files in os.walk(root_dir):
        for file_path in sorted(files):
            file_path = os.path.abspath(os.path.join(root, file_path))
            m = re.search(file_regex, file_path)
            if m:
                filename_list.append(file_path)

    return filename_list





@print_call(logger.info)
def create_output_image(pfiles, output_image_path, tmp_dir, rasterXSize):

    nrows = ncols = 0
    try:
        logger.info('Calling thumbnail.create_thumbnail(%s, %s, %s, %s, outcols=%d, work_dir=%s, overwrite=%s)',
                             pfiles[0], pfiles[1], pfiles[2],
                             output_image_path,
                             rasterXSize,
                             tmp_dir,
                             False
                             )

        (ncols, nrows) = create_thumbnail(
                             pfiles[0], pfiles[1], pfiles[2],
                             output_image_path,
                             outcols=rasterXSize,
                             work_dir=tmp_dir,
                             overwrite=False
                         )
        logger.info('Image (%s) created (%d x %d)', output_image_path, ncols, nrows)
    except Exception, e:
        logger.error('ERROR: create_thumbnail RAISED EXCEPTION (ignored)')
        logger.error('[ %s ]' % e)

    return (ncols, nrows)





@print_call(logger.info)
def create_thumbnail(red_file, green_file, blue_file, thumbnail_image,
                     outcols=1024, nodata=-999, work_dir=None, overwrite=True):
    """Create JPEG thumbnail image using individual R, G, B images.

    Arguments:
        red_file: red band data file
        green_file: green band data file
        blue_file: blue band data file
        thumbnail_image: Name of thumbnail image file
        outcols: thumbnail image width
        nodata: null/fill data value

    Thumbnail height is adjusted automatically to match the aspect ratio
    of the input images.
    """

    # working files
    file_to = "RGB.vrt"
    warp_to_file = "RGBwarped.vrt"
    outtif = "thumbnail.tif"

    if work_dir:
        file_to = os.path.join(work_dir, file_to)
        warp_to_file = os.path.join(work_dir, warp_to_file)
        outtif = os.path.join(work_dir, outtif)

    # Build the RGB Virtual Raster at full resolution
    subprocess.call(["gdalbuildvrt", "-overwrite", "-separate", file_to, red_file, green_file, blue_file], cwd=work_dir)

    # Determine the pixel scaling to get the correct width thumbnail
    vrt = gdal.Open(file_to)
    intransform = vrt.GetGeoTransform()
    inpixelx = intransform[1]
    #inpixely = intransform[5]
    inrows = vrt.RasterYSize
    incols = vrt.RasterXSize
    #print inrows,incols
    outresx = inpixelx*incols/outcols
    outrows = int(math.ceil((float(inrows)/float(incols))*outcols))
    #print outresx, outcols, outrows

    if (overwrite or not os.path.exists(thumbnail_image)):
        subprocess.call(["gdalwarp", "-of", "VRT", "-tr", str(outresx), str(outresx), "-r", "near", "-overwrite", file_to, warp_to_file], cwd=work_dir)

        # Open VRT file to array
        vrt = gdal.Open(warp_to_file)
        bands = (1,2,3)
        driver = gdal.GetDriverByName ("GTiff")
        outdataset = driver.Create(outtif,outcols,outrows, 3, gdalconst.GDT_Byte)
        #rgb_composite = numpy.zeros((outrows,outcols,3))

        # Loop through bands and apply Scale and Offset
        for bandnum, band in enumerate(bands):
            vrtband = vrt.GetRasterBand(band)
            vrtband_array = vrtband.ReadAsArray()
            nbits=gdal.GetDataTypeSize(vrtband.DataType)
            #print nbits
            dfScaleDstMin,dfScaleDstMax=0.0,255.0

            # Determine scale limits
            #dfScaleSrcMin = dfBandMean - 2.58*(dfBandStdDev)
            #dfScaleSrcMax = dfBandMean + 2.58*(dfBandStdDev)

            if (nbits == 16):
                count = 32767 + int(nodata)
                histogram = vrtband.GetHistogram(-32767, 32767, 65536)
            else:
                count = 0
                histogram = vrtband.GetHistogram()
            total = 0

            cliplower = int(0.01*(sum(histogram)-histogram[count]))
            clipupper = int(0.99*(sum(histogram)-histogram[count]))
            #print sum(histogram)
            #print cliplower,clipupper
            #print histogram[31768]
            while total < cliplower and count < len(histogram)-1:
                count = count+1
                total = int(histogram[count])+total
                dfScaleSrcMin = count
            #print "total",total
            if (nbits == 16):
                count = 32767 + int(nodata)
            else: count = 0
            #print "count for max",count
            total = 0
            while total < clipupper and count < len(histogram)-1:
                count = count+1
                #print count,clipupper,total
                total = int(histogram[count])+total
                dfScaleSrcMax = count

            if (nbits == 16):
                dfScaleSrcMin = dfScaleSrcMin - 32768
                dfScaleSrcMax = dfScaleSrcMax - 32768

            # Determine gain and offset
            dfScale = (dfScaleDstMax - dfScaleDstMin) / (dfScaleSrcMax - dfScaleSrcMin)
            dfOffset = -1 * dfScaleSrcMin * dfScale + dfScaleDstMin

            # Apply gain and offset
            outdataset.GetRasterBand(band).WriteArray((ma.masked_less_equal(vrtband_array, int(nodata))*dfScale)+dfOffset)

        outdataset = None

        # GDAL Create doesn't support JPEG so we need to make a copy of the GeoTIFF
        subprocess.call(["gdal_translate", "-of", "JPEG", outtif, thumbnail_image])

    else:
        logger.warning('File already exists. Skipping creation of %s', thumbnail_image)

    # Clean up work files
    for f in [file_to, warp_to_file, outtif]:
        try:
            os.unlink(f)
        except OSError, e:
            if e.errno != errno.ENOENT:
                raise

    return (outcols, outrows)





@print_call(logger.info)
def write_tif_file(l1t_input_dataset, band_number, dataset_id, temp_output, dtype, sublogger, work_path, debug):
    def write_band_image(image_path, band_data, geotransform, projection):
        if debug:
            sublogger.info('band_data.dtype = %s, DTYPE_MAP[band_data.dtype] = %s',
                        band_data.dtype, DTYPE_MAP[band_data.dtype])
        file_type = 'GTiff'

        nlines, npixels = band_data.shape

        dr = gdal.GetDriverByName(file_type)
        fd = dr.Create(image_path, npixels, nlines, 1, DTYPE_MAP[band_data.dtype])
        fd.SetGeoTransform(geotransform)
        fd.SetProjection(projection)
        band = fd.GetRasterBand(1)
        band.WriteArray(band_data)
        band.SetNoDataValue(-999)
        fd = None

    bin_filename = os.path.join(work_path, 'ref_wbrdf_b%d.bin' % band_number)
    band_file_number = l1t_input_dataset.sensor_band_info(band_number)['NUMBER']
    tif_filename = os.path.join(temp_output, 'scene01', '%s_B%d%s' % (dataset_id, band_file_number, '.tif'))

    band_data = numpy.fromfile(bin_filename, dtype=dtype, count=l1t_input_dataset.RasterYSize * l1t_input_dataset.RasterXSize)
    band_data.resize(l1t_input_dataset.shape)

    write_band_image(tif_filename, band_data, l1t_input_dataset.GetGeoTransform(), l1t_input_dataset.GetProjectionRef())

    sublogger.info('created geotiff product %s from binary file %s', tif_filename, bin_filename)





def load_bin_file(filename, nrow, ncol, dtype):
    """
    Loads data from one of Fuqin's 'bin' files. This is a simple wrapper around a Fortran function that Fuqin uses
    and is intended to read the output files from her programs. No protection against providing the wrong file type
    is provided... so if you get it wrong, you are very likely to seg-fault.

    :param filename:
        The name of the file to load.

    :param nrow:
        The number of rows of data in the file.

    :param ncol:
        The number of columns of data in the file.

    :param dtype:
        The type of data contained in the file.

    :type dtype:
        A numpy data type (e.g. ``numpy.float32``).
    """
    reader = None
    if dtype == numpy.int8:
        reader = read_array_int8
    elif dtype == numpy.int16:
        reader = read_array_int16
    elif dtype == numpy.int32:
        reader = read_array_int32
    elif dtype == numpy.float32:
        reader = read_array_float32

    assert reader, "Cannot read data of type %s" % str(dtype)
    return reader(filename, nrow, ncol)





def as_array(dataset, dtype):
    """
    Get a :py:class:`numpy.ndarray` for the argument ``dataset``.

    :param dataset:
        The 'dataset' to get the array from.

    :param dtype:
        The type to return the array as.

    :type dtype:
        A numpy data type (e.g. ``numpy.float32``).

    The algorithm is:

    * if ``dataset`` is a :py:class:`numpy.ndarray`
        * then if ``dataset.dtype == dtype`` return it
        * else return a copy of type ``dtype``.

    * else, if ``dataset`` is a string, try and open it using :py:func:`gdal.Open` and call :py:meth:`ReadAsArray` on the result.

    * else, call (cross our fingers and) call ``numpy.array(dataset, dtype=dtype)`` on ``dataset``.

    :return:
        A :py:class:`numpy.ndarry` of type ``dtype``.
    """
    if type(dataset) is not numpy.ndarray:
        gd_data = None
        try:
            if type(dataset) == str:
                try:
                    gd_data = gdal.Open(dataset, gdalconst.GA_ReadOnly)
                except:
                    print 'gdal could not open: "%s"' % dataset
            dataset = gd_data.ReadAsArray()
        except:
            return numpy.array(dataset, dtype=dtype)
        finally:
            gd_data = None

    if dataset.dtype != dtype:
        return dataset.astype(dtype)
    else:
        return dataset





def dump_array(
    array,
    output_path,
    template_dataset=None,
    geoc=None,
    proj=None,
    dtype=None,
    no_data_value=None,
    convert_to_byte=True,
    file_format="GTiff"):
    '''
    Save a :py:class:`numpy.ndarray` to a :py:mod:`gdal.Dataset`.

    :param array:
        The data to save.

    :param output_path:
        The path of the dataset to save to.

    :param file_format:
        The string specifying the type of the dataset to save to.

    :return:
        The :py:class:`gdal.Dataset` the data was saved to.
    '''

    geoc = geoc or template_dataset.GetGeoTransform() if template_dataset else None
    proj = proj or template_dataset.GetProjection() if template_dataset else None
    dtype = dtype or DTYPE_MAP.get(array.dtype) or gdal.GDT_Byte

    # Convert float32 to byte to make array visible
    if convert_to_byte and (dtype == gdal.GDT_Float32) and (numpy.nanmax(array) <= 1):
        dtype = gdal.GDT_Byte

    driver = gdal.GetDriverByName(file_format)
    output_dataset = driver.Create(output_path, array.shape[1], array.shape[0], 1, dtype)
    if geoc: output_dataset.SetGeoTransform(geoc)
    if proj: output_dataset.SetProjection(proj)
    if no_data_value: output_dataset.GetRasterBand(1).SetNoDataValue(no_data_value)
    output_dataset.GetRasterBand(1).WriteArray(array)

    output_dataset = None
    driver = None

    logger.debug('%d x %d array written to %s', array.shape[0], array.shape[1], output_path)
