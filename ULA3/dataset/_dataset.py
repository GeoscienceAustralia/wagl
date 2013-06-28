"""
Dataset class for generic georeferenced datasets Wraps a single GDAL dataset and provides
additional methods for lat/longs

:author:
    Alex Ip (alex.ip@ga.gov.au)

Includes code from Roger Edberg's original anci_geotiff_loader.py

"""

from osgeo import gdal, gdalconst, osr
import logging

logger = logging.getLogger('root.' + __name__)

class DSException(Exception):
    """
    Custom Exception subclass for GADatasets
    """
    pass

class Dataset(gdal.Dataset):
    """subclass of gdal.Dataset which provides additional methods to access band data
    directly without reference to gdal.Band objects.
    Emulates native C gdal.Dataset behaviour as much as possible
    """
    def __init__(self, pathname=None):
        """Initialise a Dataset instance.

        :param pathname:
            Optional string representing root pathname of scene (i.e parent of scene01)

        """
        self._pathname = None
        self._root_dataset = None
        self._eAccess = None

        if pathname:
            self.Open(pathname)

    def Open(self, pathname, eAccess = gdalconst.GA_ReadOnly):
        """
        Non-GDAL abstract function to open sub-dataset(s)
        Should be overridden in subclasses as required to allow non-GDAL dataset(s) to be read
        Returns self if valid dataset found, otherwise returns None
        Can be called as <instance> = <instance>.Open() to try to simulate gdal.Open behaviour

        :param pathname:
            String representing root pathname of scene (i.e parent of scene01 if present)
        :type pathnaame:
            :py:class:`str`

        :param eAccess:
            File open mode (:py:const:`gdalconst.GA_ReadOnly` or :py:const:`gdalconst.GA_Update`).

        """
        self._root_dataset = gdal.Open(pathname, eAccess)
        self._pathname = pathname
        self._eAccess = eAccess

        return self

    def band_read_as_array(self, nBand=1, *args, **kwargs):
        """
        Non-GDAL abstract function to return specified raster band as array
        Should be overridden in subclasses to allow data in non-GDAL datasets to be read
        Parameter nBand defaults to first available band

        """
        band = self.GetRasterBand(nBand)
        assert band, 'Band ' + str(nBand) + ' does not exist'
        return band.ReadAsArray(*args, **kwargs)

    def band_write_array(self, array, nBand=1, *args, **kwargs):
        """Non-GDAL abstract function to write array to specified raster band
        Should be overridden in subclasses to allow data in non-GDAL datasets to be read
        Parameter nBand defaults to first available band
        """
        band = self.GetRasterBand(nBand)
        assert band, 'Band ' + str(nBand) + ' does not exist'
        return band.WriteArray(array, *args, **kwargs)

    def band_data_type(self, nBand=1):
        """Non-GDAL abstract function to return specified raster band XSize
        Should be overridden in subclasses to allow data in non-GDAL datasets to be read
        Parameter nBand defaults to first available band
        """
        band = self.GetRasterBand(nBand)
        assert band, 'Band ' + str(nBand) + ' does not exist'
        return band.DataType

    def band_x_size(self, nBand=1):
        """Non-GDAL abstract function to return specified raster band XSize
        Should be overridden in subclasses to allow data in non-GDAL datasets to be read
        Parameter nBand defaults to first available band
        """
        band = self.GetRasterBand(nBand)
        assert band, 'Band ' + str(nBand) + ' does not exist'
        return band.XSize

    def band_y_size(self, nBand=1):
        """Non-GDAL abstract function to return specified raster band YSize
        Should be overridden in subclasses to allow data in non-GDAL datasets to be read
        Parameter nBand defaults to first available band
        """
        band = self.GetRasterBand(nBand)
        assert band, 'Band ' + str(nBand) + ' does not exist'
        return band.YSize

    def get_pixel_indices(self, longlat = (0.0, 0.0), nBand=1):
        """Get image pixel coordinates for a given longitude and latitude.

        Arguments:
            longlat: 2-tuple containing (longitude, latitude) in decimal degrees

        Returns:
            Index 2-tuple: (xindex, yindex)
        """
        sref = osr.SpatialReference()
        sref.ImportFromWkt(self.GetProjection())
        sref_lonlat = sref.CloneGeogCS()

        cxform = osr.CoordinateTransformation(sref_lonlat, sref)
        _x, _y, _z = cxform.TransformPoint(longlat[0], longlat[1], 0)

        g = self.GetGeoTransform()
        status, g_inv = gdal.InvGeoTransform(g)

        if status == 0:
            print 'WARNING: InvGeoTransform FAILED for g =', g

        x = int(g_inv[0] + g_inv[1]*_x + g_inv[2]*_y)
        y = int(g_inv[3] + g_inv[4]*_x + g_inv[5]*_y)

        assert x in range(0, self.RasterXSize), 'Pixel X value %d out of range (0 - %d) for longitude %f' % (x, self.RasterXSize, longlat[0])
        assert y in range(0, self.RasterYSize), 'Pixel Y value %d out of range (0 - %d) for latitude %f' % (y, self.RasterYSize, longlat[1])

        return (x, y)

    def get_data_value(self, longlat = (0.0, 0.0), nBand = 1):
        """Extract a data value for a given longitude, latitude and band.

        Arguments:
            longlat: 2-tuple containing (longitude, latitude) in decimal degrees
            band = band number (integer, optional)

        Returns:
            Data value at the specified (lon, lat) coordinates.
        """
        x, y = self.get_pixel_indices(longlat)
        band = self.GetRasterBand(nBand)
        assert band, 'Band %d does not exist' % nBand
        return band.ReadAsArray(x, y, 1, 1)[0, 0]

    @property
    def pathname(self):
        """Non-GDAL property returning pathname of dataset"""
        return self._pathname


    @property
    def pixel_size(self):
        """Gets pixel size.

        Returns:
            pixel size (2-tuple)

        N.B: This is only valid for reflective bands in a scene dataset.
        """
        g = self.GetGeoTransform()
        return (g[1], g[5])

    @property
    def origin(self):
        """Gets scene origin (upper left corner).

        Returns:
            scene origin (2-tuple)
        """
        g = self.GetGeoTransform()
        return (g[0], g[3])

    @property
    def shape(self):
        """Gets scene shape.

        Returns:
            scene shape (2-tuple: nrows, ncols)
        """
        return (self.RasterYSize, self.RasterXSize)


    # Overridden gdal.Dataset properties
    @property
    def RasterCount(self):
        """Overrides gdal.Dataset property"""
        return self._root_dataset.RasterCount

    @property
    def RasterXSize(self):
        """Overrides gdal.Dataset property"""
        return self._root_dataset.RasterXSize

    @property
    def RasterYSize(self):
        """Overrides gdal.Dataset property"""
        return self._root_dataset.RasterYSize

    # Overridden gdal.Dataset methods
    def GetDriver(self):
        """Dummy function to override gdal.Dataset.GetDriver"""
        return None

    def GetRasterBand(self, nBand):
        """Overrides gdal.Dataset method"""
        return self._root_dataset.GetRasterBand(nBand)

    def GetProjection(self):
        """Overrides gdal.Dataset method"""
        return self._root_dataset.GetProjection()

    def GetProjectionRef(self):
        """Overrides gdal.Dataset method"""
        return self._root_dataset.GetProjectionRef()

    def SetProjection(self, *args):
        """Overrides gdal.Dataset method"""
        self._root_dataset.SetProjection(*args)

    def GetGeoTransform(self, *args, **kwargs):
        """Overrides gdal.Dataset method"""
        return self._root_dataset.GetGeoTransform(*args, **kwargs)

    def SetGeoTransform(self, *args):
        """Overrides gdal.Dataset method"""
        self._root_dataset.dataset.SetGeoTransform(*args)

    def BuildOverviews(self, *args, **kwargs):
        """Overrides gdal.Dataset method"""
        return self._root_dataset.BuildOverviews(*args, **kwargs)

    def GetGCPCount(self):
        """Overrides gdal.Dataset method"""
        return self._root_dataset.GetGCPCount()

    def GetGCPProjection(self):
        """Overrides gdal.Dataset method"""
        return self._root_dataset.GetGCPProjection()

    def GetGCPs(self):
        """Overrides gdal.Dataset method"""
        return self._root_dataset.GetGCPs()

    def SetGCPs(self, *args):
        """Overrides gdal.Dataset method"""
        self._root_dataset.SetGCPs(*args);

    def FlushCache(self):
        """Overrides gdal.Dataset method"""
        self._root_dataset.dataset.FlushCache();

#    def AddBand(self, *args, **kwargs)

#    def CreateMaskBand(self, *args)

    def GetFileList(self):
        """Overrides gdal.Dataset method"""
        return self._root_dataset.GetFileList()

    def ReadRaster1(self, *args, **kwargs):
        """Overrides gdal.Dataset method"""
        return self._root_dataset.ReadRaster1(*args, **kwargs)

    def ReadAsArray(self, *args, **kwargs):
        """Overrides gdal.Dataset method"""
        return self._root_dataset.ReadAsArray(*args, **kwargs)

    def WriteRaster(self, *args, **kwargs):
        """Overrides gdal.Dataset method"""
        return self._root_dataset.WriteRaster(*args, **kwargs)

    def ReadRaster(self, *args, **kwargs):
        """Overrides gdal.Dataset method"""
        return self._root_dataset.ReadRaster(*args, **kwargs)

    def GetSubDatasets(self):
        """Overrides gdal.Dataset method"""
        return list(self._root_dataset)

#    def BeginAsyncReader(self, *args, **kwargs):
#        """Overrides gdal.Dataset method"""
#        return self._root_dataset.BeginAsyncReader(*args, **kwargs)

#    def EndAsyncReader(self, *args, **kwargs):
#        """Overrides gdal.Dataset method"""
#        return self._root_dataset.EndAsyncReader(*args, **kwargs)

    # Overridden gdal.MajorObject methods
    def GetDescription(self):
        """Overrides gdal.MajorObject method"""
        return self._pathname

    def GetMetadata(self, domain = ''):
        """Overrides gdal.MajorObject method"""
        return self._root_dataset.GetMetadata(domain)

    def GetMetadataItem(self, name, domain = ''):
        """Overrides gdal.MajorObject method"""
        return self._root_dataset.GetMetadataItem(name, domain)

    def GetMetadata_Dict(self, domain = ''):
        """Overrides gdal.MajorObject method"""
        return self._root_dataset.GetMetadata_Dict(domain)

    def GetMetadata_List(self, domain = ''):
        """Overrides gdal.MajorObject method"""
        return self._root_dataset.GetMetadata_List(domain)

    def SetDescription(self, newDesc):
        """Overrides gdal.MajorObject method"""
        return None # Todo

    def SetMetadata(self, metadata, domain):
        """Overrides gdal.MajorObject method"""
        return None # Todo

    def SetMetadataItem(self, *args):
        """Overrides gdal.MajorObject method"""
        return None # Todo
