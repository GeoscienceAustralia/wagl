"""
This is the main package developed during the ULA project. Most of the code herein is designed to be useful outside
in its own right, though there are still some dependencies on the specifics of the GA environment that prevent parts
of it from being so. These should be removed in later versions.
"""
import pickle, os, threading, logging, re, numpy
from EOtools import blrb
from osgeo import gdal, gdalconst
from meta import Singleton
from dataset import Dataset
from ULA3.image_processor import ProcessorConfig
from _gdal_tools import NUMPY_GDAL_TYPE_MAP as DTYPE_MAP

logger = logging.getLogger('root.' + __name__)





class DataManager(object):
    """
    A "shelf" for shared data. Each module will be responsible for marshalling its own prerequisites by using the
    :py:meth:`get_item` method on this class. This should be the only way that data can be passed between
    'modules' in the image_processor.

    This class permits us to either run a process end-to-end passing data in memory, or, if required, we
    can cache objects to disk for later retrieval (or debugging).

    :todo:
        This isn't really the appropriate place for this class and it should be moved to
        :py:mod:`ULA3.image_processor`.
    """
    __metaclass__ = Singleton

    def __init__(self, work_path=None):
        self.work_path = work_path or ProcessorConfig().work_path
        self._item_dict = {} # Class dict holding object instances keyed by (item_path, item_type) tuple
        self._data_lock = threading.Lock()

    def get_abs_path(self, item_path):
        """
        Prepend work directory to any relative pathname.
        """
        result = item_path.strip()
        if result[0] != os.pathsep:
            result = os.path.join(self.work_path, result)
        return os.path.abspath(result)

    def get_item(self, item_path, item_type, load_from_disk=True):
        '''
        Function to return object instance or None if not found.

        This first checks whether a reference to the requested object in memory, and if not, will attempt to
        retrieve it from disk.

        :return:
             The requested object if it can be retrieved, ``None`` otherwise.

        '''
        logger.debug('get_item(%s, %s, %s) called', repr(item_path), repr(item_type), repr(load_from_disk))

        item_path = self.get_abs_path(item_path)

        #self._data_lock.acquire()
        with self._data_lock:
            item_instance = self._item_dict.get((item_path, item_type))
            if item_instance is None and load_from_disk: # Try to load item from disk if not in memory
                try:
                    infile = None
                    if issubclass(item_type, Dataset):
                        item_instance = item_type()
                        item_instance.Open(item_path)
                    elif issubclass(item_type, DataGrid):
                        item_instance = item_type(filename=item_path)
                    else: # Default to pickling for unknown types
                        infile = open(item_path, 'r')
                        item_instance = pickle.load(infile)

                    self._item_dict[(item_path, item_type)] = item_instance
                except (Exception), e:
                    logger.warning('Unable to load %s object from file %s: %s', item_type.__name__, item_path, e.message)
                    item_instance = None
                finally:
                    if infile:
                        infile.close()

        logger.debug('%s object %s = %s', item_type.__name__, item_path, item_instance)
        return item_instance

    def set_item(self, item_path, item_instance):
        '''
        Function to register existing object instance.

        '''
        logger.debug('set_item(%s, %s) called', repr(item_path), repr(item_instance))

        item_path = self.get_abs_path(item_path)

        item_type = item_instance.__class__

        with self._data_lock:
            self._item_dict[(item_path, item_type)] = item_instance

        logger.debug('%s object %s set to %s in data_manager', item_type.__name__, item_path, item_instance)
        return item_instance


    def save_item(self, item_path, item_type=None, item_instance=None):
        '''
        Function to save object instance to disk. The object instance can either be provided or specified
        by (item_path, item_type)

        '''
        logger.debug('save_item(%s, %s, %s) called', repr(item_path), repr(item_type), repr(item_instance))

        assert (item_type is None) != (item_instance is None), 'Either item_type or item_instance must be specified (but not both)'

        item_path = self.get_abs_path(item_path)

        if item_instance is not None:
            self.set_item(item_path, item_instance=item_instance) # Register new item
            item_type = item_instance.__class__
        elif item_type is not None:
            item_instance = self._item_dict.get((item_path, item_type))
            assert item_instance is not None, 'No object found to save to disk'

        with self._data_lock:
            outfile = None
            try:
                if issubclass(item_type, Dataset):
                    # Don't do anything with this - Dataset should already exist on disk
                    raise Exception(item_type.__name__ + ' objects are read-only')
                elif issubclass(item_type, DataGrid):
                    item_instance.save(item_path)
                else:
                    outfile = open(item_path, 'w')
                    pickle.dump(item_instance, outfile, pickle.HIGHEST_PROTOCOL)

                logger.info('%s object written to file %s', item_type.__name__, item_path)

                return item_instance
            except (Exception), e:
                logger.warning('Unable to write %s object to file %s: %s', item_type.__name__, item_path, e.message)
                item_instance = None
                if not issubclass(item_type, Dataset): # Do NOT remove Dataset files
                    try:
                        os.remove(item_path) # Clean up failed save
                    except (OSError):
                        pass
            finally:
                if outfile:
                    outfile.close()


    def free_item(self, item_path, item_type, save_first=False):
        '''
        Function to remove object instance.

        '''
        logger.debug('free_item(%s, %s, %s) called', repr(item_path), repr(item_type), repr(save_first))

        item_path = self.get_abs_path(item_path)

        if (item_path, item_type) in self._item_dict:
            if save_first:
                self.save_item(item_path, item_type)
            with self._data_lock:
                del self._item_dict[(item_path, item_type)]
                logger.info('%s object %s deleted from data_manager', item_type.__name__, item_path)
        else:
            assert not save_first, 'Unable to save a non-existent %s object as %s' % (item_type.__name__, item_path)

    def save_all(self):
        """
        Saves all items - useful for debugging or when processing partial workflows.

        """
        logger.debug('save_all() called')
        for item_path, item_type in self._item_dict.keys():
            self.save_item(item_path, item_type)

    def clean(self):
        """
        Free all items. This is useful when memory limits are being reached, and one needs to find a little
        more space.

        """
        for item_path, item_type in self._item_dict.keys():
            try:
                self.free_item(item_path, item_type, True)
            except Exception, e:
                logger.info("problem freeing item: %s" % str(e))





class DataGrid(object):
    """
    Class that models a grid of data and manages serialisation as required.
    """
    def __init__(self, array=None, filename=None, dtype=numpy.float64, shape=None, eval_func=None, depth=1):
        """
        Instantiates a DataGrid object either from an existing array or using a datatype and shape specifier
        """
        logger.debug('DataGrid(%s, %s, %s, %s, %s, %s) called', array, filename, dtype, shape, eval_func, depth)
        assert (array is not None) or filename or (shape and (eval_func or dtype)), 'Insufficient information provided to instantiate DataGrid object'

        self._filename = None

        if (array is not None):
            self.array = array # Use an existing array in memory BY REFERENCE
            logger.debug('grid array set by reference to existing array')
        elif filename:
            self.load(filename) # Load the array from an existing file
            self._filename = filename
            logger.debug('grid array loaded from file ' + filename)
        elif eval_func:
            self.interpolate(depth=depth, shape=shape, eval_func=eval_func, dtype=dtype)
            logger.debug('grid array set by interpolating ' + eval_func.__name__)
        else:
            self.array = numpy.zeros(shape, dtype) # Create new empty array
            logger.debug('grid array set to zeroes')

    def interpolate(self, depth=1, shape=None, eval_func=None, dtype=numpy.float64):
        """
        Construct an interpolated scalar data grid.

        :param depth:
            Bisection recursion depth.
        :type depth:
            int

        :param shape:
            The shape of the grid to be created (see first argument of :py:func:`numpy.zeros`).
        :type shape:
            iterable of ints

        :param interpolator:
            Evaluator function.
        :type interpolator:
            callable

        :param dtype:
            Grid data type.
        :type dtype:
            :py:class:`numpy.dtype`

        Returns:
            Interpolated data grid.

        The eval_func argument should be a function that accepts grid
        indices i and j as arguments and returns a scalar data value.
        """
        logger.debug('DataGrid.interpolate(%s, %s, %s, %s) called', depth, shape, eval_func, dtype)

        assert depth > 0
        assert shape is not None
        assert eval_func is not None

        self.array = numpy.zeros(shape, dtype=dtype)
        blrb.interpolate_grid(depth=depth, origin=(0, 0), shape=shape,
                         eval_func=eval_func, grid=self.array)
        return self.array

    def save_binary(self, filename=None, convert_angles=False, dtype=None):
        """Saves numpy array as a binary dump.
        """
        logger.debug('DataGrid.save_binary(%s, %s, %s) called', filename, convert_angles, repr(dtype))
        filename = filename or self._filename or 'GRID.dat'
        assert filename, 'No filename given to save'

        if not dtype:
            dtype = self.array.dtype
        else:
            dtype = numpy.dtype(dtype)

        if convert_angles: # Convert angles to degrees for output
            out_array = numpy.cast[dtype](self.array)
            numpy.degrees(out_array, out_array)
        else:
            if dtype == self.array.dtype:
                out_array = self.array
            else:
                out_array = numpy.cast[dtype](self.array)

        out_array.tofile(filename)
        logger.debug('DataGrid saved to binary file %s', filename)


    def load_binary(self, filename=None, convert_angles=False, dtype=None):
        """
        """
        logger.debug('DataGrid.load_binary(%s, %s, %s) called', filename, convert_angles, dtype)
        filename = filename or self._filename or 'GRID.dat'
        assert filename, 'No filename given to load'

        if not dtype:
            dtype = numpy.float64
        else:
            dtype = numpy.dtype(dtype)

        self.array = numpy.load(filename, dtype=dtype)
        self._filename = filename
        if convert_angles: # Convert angles to radians after input
            numpy.radians(self.array, self.array)
        logger.debug('DataGrid loaded from binary file %s', filename)


    def save_image(
        self,
        filename=None,
        image_type='GTiff',
        dtype=None,
        projection_ref=None,
        geotransform=None,
        convert_angles=False):
        """Write a numpy array to an image file.

        Arguments:
            filename: image file path
            image_type: image file type (GDAL type specifier)
            dtype: image data type (GDAL data type specifier)
            projection_ref: projection reference
            geotransform: geotransform (6-tuple)
            convert_angles: Should the data be converted to angles.
        """
        logger.debug('DataGrid.save_image(%s, %s, %s, %s, %s, %s) called',
                     filename, image_type, dtype, projection_ref, geotransform, convert_angles)

        filename = filename or self._filename or 'GRID.tif'
        assert filename, 'No filename given to save'

        if not dtype:
            dtype = self.array.dtype
        else:
            dtype = numpy.dtype(dtype)

        logger.debug('numpy dtype = %s', repr(dtype))
        gdal_dtype = DTYPE_MAP[numpy.dtype(dtype)]
        logger.debug('GDAL dtype = %s', repr(gdal_dtype))

        try:
            if convert_angles: # Convert angles to degrees for output
                out_array = numpy.cast[dtype](self.array)
                numpy.degrees(out_array, out_array)
            else:
                if dtype == self.array.dtype:
                    out_array = self.array
                else:
                    out_array = numpy.cast[dtype](self.array)

            gdal_driver = gdal.GetDriverByName(image_type)
            logger.debug('GDAL driver found')
            logger.debug('Calling gdal_driver.Create(%s, %s, %s, %s, %s)', filename, out_array.shape[1], out_array.shape[0], 1, repr(gdal_dtype))
            gdal_dataset = gdal_driver.Create(filename, out_array.shape[1], out_array.shape[0], 1, gdal_dtype)
            logger.debug('GDAL dataset created')
            if projection_ref:
                gdal_dataset.SetProjection(projection_ref)
            if geotransform:
                gdal_dataset.SetGeoTransform(geotransform)
            gdal_dataset.GetRasterBand(1).WriteArray(out_array, 0, 0)
            logger.debug('DataGrid saved to GDAL %s file %s', image_type, filename)
        finally:
            gdal_dataset = None

    def load_image(self, filename, convert_angles=False):
        """Read a numpy array from an image file
        """
        logger.debug('DataGrid.load_image(%s) called', filename)
        filename = filename or self._filename
        assert filename, 'No filename given to open'

        try:
            gdal_dataset = gdal.Open(filename)
            self.array = gdal_dataset.ReadAsArray()
            self._filename = filename
            if convert_angles: # Convert angles to radians
                numpy.radians(self.array, self.array)
            logger.info('DataGrid loaded from GDAL file %s', filename)
        finally:
            gdal_dataset = None

    def load(self, filename=None, convert_angles=False):
        """
        """
        logger.debug('DataGrid.load(%s) called', filename)
        filename = filename or self._filename or 'GRID.dat'
        assert filename, 'No filename given to load'

        s = re.search('\..*$', filename)
        if s and s.group(0).lower() in ['.tif', '.tiff']:
            self.load_image(filename, convert_angles=convert_angles)
        else:
            self.load_binary(filename, convert_angles=convert_angles)

    def save(self, filename=None, convert_angles=False):
        """
        """
        logger.debug('DataGrid.save(%s) called', filename)
        filename = filename or self._filename or 'GRID.dat'
        assert filename, 'No filename given to save'

        s = re.search('\..*$', filename)
        if s and s.group(0).lower() in ['.tif', '.tiff']:
            self.save_image(filename, convert_angles=convert_angles)
        else:
            self.save_binary(filename, convert_angles=convert_angles)
