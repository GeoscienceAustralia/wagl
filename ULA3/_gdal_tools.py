"""
Some tools used to work with GDAL.
"""

import os, sys, logging, numpy
from osgeo import gdal, osr
from _execute import execute
from ULA3.meta import print_call

logger = logging.getLogger('root.' + __name__)

NUMPY_GDAL_TYPE_MAP = {
    # numpy dtype: GDAL dtype
    numpy.dtype('uint8'):   gdal.GDT_Byte,
    numpy.dtype('uint16'):  gdal.GDT_UInt16,
    numpy.dtype('int16'):   gdal.GDT_Int16,
    numpy.dtype('float32'): gdal.GDT_Float32,
    numpy.dtype('float64'): gdal.GDT_Float64,
    numpy.uint8:   gdal.GDT_Byte,
    numpy.uint16:  gdal.GDT_UInt16,
    numpy.int16:   gdal.GDT_Int16,
    numpy.float32: gdal.GDT_Float32,
    numpy.float64: gdal.GDT_Float64}

class Buffers(object):
    """
    Holds some value for each side of an image. This was initially created to hold buffer widths (in pixels) for
    a scene, but does not care about the type of the values passed to the constructer and can hence be used for
    any type.
    """
    def __init__(self, left, right=None, top=None, bottom=None):
        """
        Constructor.

        The arguments are copied directly to members of the same name with the restrictions that:

        - if right is None then top and bottom must also be None, and

        - if right is not None then top and bottom must not be none either.
        """
        if right is None:
            assert top is None and bottom is None, "if right is None then top and bottom must also be None"
            self.left = self.right = self.top = self.bottom = left
        else:
            assert top is not None and bottom is not None, "if right is not None then top and bottom must also not be None"
            self.left = left
            self.right = right
            self.top = top
            self.bottom = bottom

    def __str__(self):
        return "ULA3.utils.Buffers(%s)" % ", ".join(map(str, (self.left, self.right, self.top, self.bottom)))






class ImageShape(object):
    """
    Holds a few values that describe the shape of an image.
    """
    def __init__(self, origin_x, origin_y, dim_x, dim_y, pix_sz_x, pix_sz_y=None, proj=None):
        """
        Constructor

        :param origin_x:
            The X origin of the image (top left).
        :type origin_x:
            float

        :param origin_y:
            The Y origin of the image (top left).
        :type origin_y:
            float

        :param dim_x:
            The width of the image in pixels.
        :type dim_x:
            int

        :param dim_y:
            The height of the image in pixels.
        :type dim_y:
            int

        :param pix_sz_x:
            The width of image pixels (in the appropriate coordinate space).
        :type pix_sz_x:
            float

        :param pix_sz_y:
            The height of image pixels (in the appropriate coordinate space).
        :type pix_sz_y:
            float

        :param proj:
            The projection (WKT).
        :type proj:
            str
        """
        self.__origin = (origin_x, origin_y)
        self.__pix_sz = (pix_sz_x, pix_sz_y or pix_sz_x)
        self.__dim = (dim_y, dim_x)
        self.__spatial_ref = osr.SpatialReference()
        if proj:
            self.__spatial_ref.ImportFromWkt(proj)

    def __str__(self):
        return "ULA3.utils.ImageShape(%s)" % ", ".join(map(str, (self._ImageShape__origin, self._ImageShape__pix_sz, self._ImageShape__dim)))

    @property
    def RasterXSize(self):
        """
        The number of columns in the raster.
        """
        return self._ImageShape__dim[1]

    @property
    def RasterYSize(self):
        """
        The number of rows in the raster.
        """
        return self._ImageShape__dim[0]

    @property
    def RasterCellSize(self):
        """
        The cell size of the raster.

        This asserts that the cell size is the same in both dimensions. If the cell size is different, you must use the
        X and Y getters instead.
        """
        assert self._ImageShape__pix_sz[0] == abs(self._ImageShape__pix_sz[1]), "Cannot use method ImageShape.CellSize when X and Y cell sizes are different."
        return self._ImageShape__pix_sz[0]

    @property
    def RasterXCellSize(self):
        """
        The width of a cell in the raster.
        """
        return self._ImageShape__pix_sz[0]

    @property
    def RasterYCellSize(self):
        """
        The height of a cell in the raster.
        """
        return self._ImageShape__pix_sz[1]

    @property
    def RasterXOrigin(self):
        """
        The X element of the origin of the raster.
        """
        return self._ImageShape__origin[0]

    @property
    def RasterYOrigin(self):
        """
        The Y element of the origin of the raster.
        """
        return self._ImageShape__origin[1]

    @property
    def shape(self):
        """
        Tuple containing the number of columns and rows (in that order).
        """
        return (self.RasterYSize, self.RasterXSize)

    def GetProjection(self, as_proj4=False):
        if as_proj4:
            return self._ImageShape__spatial_ref.ExportToProj4()
        else:
            return self._ImageShape__spatial_ref.ExportToWkt()

    def write_header_slope_file(self, file_name, bounds):
        output = open(file_name, 'w')
        output.write("%i %i\n" % (self.RasterYSize, self.RasterXSize))
        output.write("%i %i\n%i %i\n" % (bounds.left, bounds.right, bounds.top, bounds.bottom))
        output.write("%f\n" % self.RasterCellSize)
        output.write("%f %f\n" % (self.RasterYOrigin, self.RasterXOrigin))
        output.close()







def default_bounds_getter(dataset):
    """
    Get the various bits and pieces needed to determine the 'shape' of an image (origin, cell size and
    number of rows and cols). This function aids in making access to a :py:class:`ImageShape` for an
    object easy and hence making the code that uses it more generic. The algorithm is as follows:

    - If ``dataset`` is an instance of :py:class:`ImageShape`, then it is simply returned.

    - Otherwise, if ``dataset`` has an attribute ``bounds_getter`` it is retrieved and passed ``dataset``
        (see :py:meth:`ULA3.dataset._scene_dataset.SceneDataset.bounds_getter` for an example).

    - Otherwise ``dataset`` must be an instance of :py:class:`gdal.Dataset` (a :py:class:`TypeError` is
        raised otherwise), and the shape is generated with:

        >>> g = dataset.GetGeoTransform()
        >>> return ImageShape(g[0], g[3], dataset.RasterXSize, dataset.RasterYSize, g[1], g[5], dataset.GetProjection())

    :param dataset:
        The dataset to extract the required elements from.

    :return:
        A :py:class:`ImageShape` for ``dataset``.
    """
    if isinstance(dataset, ImageShape):
        return dataset
    elif hasattr(dataset, 'bounds_getter'):
        return dataset.bounds_getter(dataset)
    elif isinstance(dataset, gdal.Dataset):
        g = dataset.GetGeoTransform()
        return ImageShape(g[0], g[3], dataset.RasterXSize, dataset.RasterYSize, g[1], g[5], dataset.GetProjection())
    else:
        raise TypeError("dataset must be instance of ImageShape, gdal.Dataset or ULA3.dataset.Dataset.")





@print_call(logger.info)
def warp(shape_dataset, master_dataset_path, output_filename, buffer_widths, output_format, resammple_method="bilinear", bounds_getter=default_bounds_getter):
    """
    Use the gdalwarp executable to clip (and potentially resample) a region.

    Preconditions on this method are:
        - the directory specified for the output (``output_filename``) exists,
        - that ``master_dataset_path`` exists, and
        - gdalwarp is on the path.

    :param shape_dataset:
        Object to extract the shape of the desired region from. This is done using ``bounds_getter``
        (see :py:func:`default_bounds_getter` for specification of the interface).

    :param master_dataset_path:
        The path to the dataset to clip from. This should be a valid argument to :py:func:`gdal.Open`.
    :type shape_dataset:
        str

    :param output_filename:
        The name of the output file (passed as the output argument to gdalwarp).
    :type shape_dataset:
        str

    :param buffer_widths:
        An object of type :py:class:`ImageShape` (or one that supports the same interface).
    :type buffer_widths:
        :py:class:`Buffers`

    :param output_format:
        The desired format of the clipped dataset. (passed as argument -of to gdalwarp).
    :type output_format:
        str

    :param resample_method:
        The resampling method to be used (passed as argument -r to gdalwarp).
    :type resample_method:
        str

    :param bounds_getter:
        Callable used to extract the bounds from ``shape_dataset``.

    :return:
        The name of the dataset written to disk.
    """
    assert not execute("which gdalwarp")["returncode"], "gdalwarp not available"
    output_dir = os.path.dirname(output_filename)

    assert os.stat(output_dir), "output_dir (%s) must exist." % output_dir
    assert os.stat(master_dataset_path), "master_dataset (%s) must exist" % master_dataset_path

    shape = bounds_getter(shape_dataset)

    xres = shape.RasterXCellSize
    yres = shape.RasterYCellSize
    xmin = shape.RasterXOrigin - xres*buffer_widths.left
    xmax = shape.RasterXOrigin + xres*(shape.RasterXSize + buffer_widths.right)
    ymax = shape.RasterYOrigin - yres*buffer_widths.top # in the cases I've looked at, yres is negative.
    ymin = shape.RasterYOrigin + yres*(shape.RasterYSize + buffer_widths.bottom)

    command_string = 'gdalwarp -overwrite -of %s -t_srs "%s" -r %s -te %f %f %f %f -tr %f %f %s %s' % (
        output_format,
        shape.GetProjection(as_proj4=True),
        resammple_method,
        float(xmin), float(ymin), float(xmax), float(ymax),
        float(xres), float(yres),
        master_dataset_path, output_filename)

    result = execute(command_string)
    if result["returncode"]:
        print "error in executing %s\n\n\tstdout: %s\n\n\tstderr: %s\n" % (command_string, result['stdout'], result['stderr'])

    assert not result["returncode"], "error in executing %s\n\n\tstdout: %s\n\n\tstderr: %s\n" % (command_string, result['stdout'], result['stderr'])
    return output_filename
