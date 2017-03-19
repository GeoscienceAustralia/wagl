"""
Gridded Data
"""
from __future__ import absolute_import, print_function, unicode_literals
import math
from math import radians
import gdal
import h5py
import rasterio as rio
import osr
import affine
from affine import Affine
from gaip.calculate_angles import setup_spheroid
from gaip.vincenty import vinc_dist


# Landsat tranforms have very small determinants
# the following setting is required, and there is a
# bug in rasterio.set_epsilon() ver < 1.0.5
affine.EPSILON = 1e-9
affine.EPSILON2 = 1e-18


# WGS84
CRS = "EPSG:4326"


class GriddedGeoBox(object):

    """
    Represents a north up rectangular region on the Earth's surface which
    has been divided into equal size retangular pixels for the purpose of 
    data processing.

    The position and extent of the GriddedGeoBox on the Earth's
    surface are defined with respect to a specific Coordinate Reference
    System (CRS) and comprise:

      - the Origin (x, y) of the North West corner of the box,
      - the pixel size (sizeX, sizeY) expressed in CRS units
      - the shape of the GGB (yPixels, xPixels) where each value is
        a positive integer

    Pixels within a GriddedGeoBox may be reference by their (x,y) grid
    position relative to the origin pixel in the upper (northwest) corner.
    x values increase to the east, y value increase to the south.

    All pixels are grid-aligned. Origin and corner will always be 
    located on grid corners
    """

    @staticmethod
    def from_dataset(dataset):
        """
        Return the GriddedGeoBox that encloses the full extent of 
        the supplied Rasterio or GDAL dataset.

        :param dataset:
            An open rasterio or GDAL Dataset.

        :return:
            A GriddedGeoBox that encloses the full extent of
            the supplied dataset.
        """
        if isinstance(dataset, gdal.Dataset):
            return GriddedGeoBox.from_gdal_dataset(dataset)
        elif isinstance(dataset, rio.io.DatasetReader):
            return GriddedGeoBox.from_rio_dataset(dataset)
        elif isinstance(dataset, h5py.Dataset):
            return GriddedGeoBox.from_h5_dataset(dataset)
        else:
            raise ValueError("GriddedGeoBox.from_dataset() expects"
                             " GDAL or rasterio dataset, not %s" %
                             type(dataset))

    @staticmethod
    def from_rio_dataset(dataset):
        """
        Return the GriddedGeoBox that encloses the full extent of 
        the supplied Rasterio dataset.

        :param dataset:
            An open rasterio Dataset.

        :return:
            A GriddedGeoBox that encloses the full extent of
            the supplied dataset.
        """
        bbshape = dataset.shape
        origin = (dataset.affine[2], dataset.affine[5])
        pixelsize = dataset.res
        crsString = bytes(dataset.crs.wkt)

        return GriddedGeoBox(bbshape, origin, pixelsize, crsString)

    @staticmethod
    def from_gdal_dataset(dataset):
        """
        Return the GriddedGeoBox that encloses the full extent of 
        the supplied GDAL dataset.

        :param dataset:
            An open GDAL Dataset.

        :return:
            A GriddedGeoBox that encloses the full extent of
            the supplied dataset.
        """
        bbshape = (dataset.RasterYSize, dataset.RasterXSize)
        transform = dataset.GetGeoTransform()
        origin = (transform[0], transform[3])
        pixelsize = (abs(transform[1]), abs(transform[5]))
        crsString = str(dataset.GetProjection())

        return GriddedGeoBox(bbshape, origin, pixelsize, crsString)

    @staticmethod
    def from_h5_dataset(dataset):
        """
        Return the GriddedGeoBox that encloses the full extent of 
        the supplied Rasterio dataset.

        :param dataset:
            An open rasterio Dataset.

        :return:
            A GriddedGeoBox that encloses the full extent of
            the supplied dataset.
        """
        bbshape = dataset.shape
        transform = dataset.attrs['geotransform']
        origin = (transform[0], transform[3])
        crs = dataset.attrs['crs_wkt']
        pixelsize = (abs(transform[1]), abs(transform[5]))

        return GriddedGeoBox(bbshape, origin, pixelsize, crs)

    @staticmethod
    def from_corners(origin, corner, pixelsize=(0.00025, 0.00025),
                     crs='EPSG:4326'):
        """
        Return a GriddedGeoBox defined by the the two supplied
        corners.

        :param origin:
            A tuple (in CRS coordinates) representing the positon of
            the NW corner of the GriddedGeoBox.

        :param corner:
            A tuple (in CRS coordinates) representing the positon of
            the NW corner of the GriddedGeoBox.

        :param pixelsize:
             Pixel (xSize, ySize) in CRS units.

        :param crs:
            The SpatialReferenceSystem in which both origin and
            pixelsize are expressed (supports various text formats).
        """
        a = Affine(pixelsize[0], 0, origin[0], 0, -pixelsize[1], origin[1])
        shapeXY = tuple([int(math.ceil(v)) for v in ~a * corner])
        shapeYX = (shapeXY[1], shapeXY[0])
        return GriddedGeoBox(shapeYX, origin, pixelsize, crs)

    def equals(self, geobox):
        """
        Compare this GriddedGeoBox with the supplied geobox. Return
        true if they are actually or functionally the same object
        """

        if self == geobox:
            return True

        if self.origin != geobox.origin:
            return False

        if self.shape != geobox.shape:
            return False

        if self.pixelsize != geobox.pixelsize:
            return False

        if self.crs.ExportToWkt() != geobox.crs.ExportToWkt():
            return False

        return True


    def __init__(self, shape=(1, 1), origin=(0.0, 0.0),
                 pixelsize=(0.00025, 0.00025), crs='EPSG:4326'):
        """
        Create a new GriddedGeoBox.

        :param shape:
            (ySize, xSize) 2-tuple defining the shape of the GGB.

            * Use get_shape_xy() to get (xSize, ySize).

        :param origin:
            (xPos, yPos) 2-tuple difining the upper left GGB corner.

        :param pixelsize:
            Pixel (xSize, ySize) in CRS units.

        :param crs:
            The SpatialReferenceSystem in which both origin and
            pixelsize are expressed (supports various text formats).
        """
        self.pixelsize = pixelsize
        self.shape = tuple([int(v) for v in shape])
        self.origin = origin
        if isinstance(crs, osr.SpatialReference):
            self.crs = crs
        else:
            self.crs = osr.SpatialReference()
            if self.crs == self.crs.SetFromUserInput(crs):
                raise ValueError("Invalid crs: %s" % (crs, ))
        self.affine = Affine(self.pixelsize[0], 0, self.origin[0], 0,
                             -self.pixelsize[1], self.origin[1])
        self.corner = self.affine * self.get_shape_xy()

    def get_shape_xy(self):
        """Get the shape as a tuple (x,y)."""
        return (self.shape[1], self.shape[0])

    def get_shape_yx(self):
        """Get the shape as a tuple (y,x)."""
        return self.shape

    def transform_point(self, transformation, point):
        """Transform the point."""
        (x, y, _) = transformation.TransformPoint(point[0], point[1])
        return (x, y)

    def copy(self, crs='EPSG:4326'):
        """
        Create a copy of this GriddedGeoBox transformed to the supplied
        Coordinate Reference System. The new GGB will have idential shape
        to the old and will be grid aligned to the new CRS. Pixel size 
        may change to accommodate the new CRS.
        """
        newCrs = osr.SpatialReference()
        newCrs.SetFromUserInput(crs)
        old2New = osr.CoordinateTransformation(self.crs, newCrs)
        newOrigin = self.transform_point(old2New, self.origin)
        newCorner = self.transform_point(old2New, self.corner)
        newPixelSize = tuple([
            abs((newOrigin[0] - newCorner[0]) / self.get_shape_xy()[0]),
            abs((newOrigin[1] - newCorner[1]) / self.get_shape_xy()[1])
        ])

        return GriddedGeoBox(self.shape, newOrigin, newPixelSize, crs=crs)

    def __str__(self):
        return 'GriddedGeoBox(origin=%s,shape=%s,pixelsize=%s,crs: %s)' % \
            (self.origin, self.shape, str(self.pixelsize),
             self.crs.ExportToProj4())

    def window(self, enclosedGGB):
        """
        Return the window, (expressed in the grid coordinates of this
        GriddedGeoBox) of the supplied GriddedGeoBox. Self must be a
        GriddedGeoBox which fully encloses the enclosed GGB.

        :param enclosedGGB:
            A GGB which is a subset of this GGB and is fully enclosed
            by it.

        :return:
            A pair of range tuples defining the enclosed retangular subset
            ((row_start, row_stop), (col_start, col_stop)).
        """
        # transform to map enclosedGGB coords to self.crs coordinates
        enclosed2self = osr.CoordinateTransformation(enclosedGGB.crs,
                                                     self.crs)

        pos = enclosed2self.TransformPoint(enclosedGGB.origin[0],
                                           enclosedGGB.origin[1])
        originT = (pos[0], pos[1])

        areaCorner = enclosedGGB.corner
        pos = enclosed2self.TransformPoint(areaCorner[0], areaCorner[1])
        areaCornerT = (pos[0], pos[1])

        areaOriginXY = [int(math.floor(v)) for v in ~self.affine * originT]
        areaCornerXY = [int(math.ceil(v)) for v in ~self.affine * areaCornerT]

        return ((areaOriginXY[1], areaCornerXY[1]),
                (areaOriginXY[0], areaCornerXY[0]))

    def x_size(self):
        """The x-axis size."""
        return self.shape[1]

    def y_size(self):
        """The y-axis size."""
        return self.shape[0]

    def convert_coordinates(self, xy, to_map=True, centre=False):
        """
        Given a tuple containing an (x, y) co-ordinate pair, convert
        the co-ordinate pair to either image/array co-ordinates or
        real world (map) co-ordinates.

        :param xy:
            A tuple containing an (x, y) co-ordinate pair. The pair
            can be either image/array co-ordinates or map
            co-ordinates. If image co-ordinates are input, then set
            to_map=True. If map co-ordinates are input, then set
            to_map=False.

        :param to_map:
            A boolean indicating if the conversion should be image to
            map or map to image. Default is True (image to map).

        :param centre:
            A boolean indicating if the returned co-ordinate pair
            should be offset by 0.5 indicating the centre of a pixel.
            Default is False.

        :return:
            A tuple containing an (x, y) co-ordinate pair.
            The returned type will be int if to_map=False and float
            if to_map=True (Default).
        """
        if to_map:
            if centre:
                xy = tuple(v + 0.5 for v in xy)
            x, y = xy * self.affine
        else:
            inv = ~self.affine
            x, y = [int(v) for v in inv * xy]

        return (x, y)

    def transform_coordinates(self, xy, to_crs):
        """
        Transform a tuple co-ordinate pair (x, y) from one CRS to
        another.

        :param xy:
            A tuple containing an (x, y) co-ordinate pair of real
            world co-ordinates.

        :param to_crs:
            An instance of a defined osr.SpatialReference object.

        :return:
            A tuple (x, y) floating point co-ordinate pair.
        """
        if not isinstance(to_crs, osr.SpatialReference):
            err = 'Err: to_crs is not an instance of osr.SpatialReference: {}'
            err = err.format(type(to_crs))
            raise TypeError(err)

        # Define the transform we are transforming to
        transform = osr.CoordinateTransformation(self.crs, to_crs)

        x, y = self.transform_point(transform, xy)

        return (x, y)
    
    def get_pixelsize_metres(self, xy=None):
        """
        Compute the size (in metres) of the pixel at the specified xy position.

        :param xy:
            A tuple containing an (x, y) grid co-ordinates of a pixel. Defaults \
            to the central pixel in the grid.

        :return:
            A tuple (x_size, y_size) gives the size of the pixel in metres
        """
        if xy is None:
            xy = (self.shape[1]/2, self.shape[0]/2)

        (x, y) = xy

        spheroid, _ = setup_spheroid(self.crs.ExportToWkt())

        (lon1, lat1) = self.affine * (x, y+0.5)
        (lon2, lat2) = self.affine * (x+1, y+0.5)
        x_size, _az_to, _az_from = vinc_dist(spheroid[1], spheroid[0],
                                             radians(lat1), radians(lon1),
                                             radians(lat2), radians(lon2))

        (lon1, lat1) = self.affine * (x+0.5, y)
        (lon2, lat2) = self.affine * (x+0.5, y+1)
        y_size, _az_to, _az_from = vinc_dist(spheroid[1], spheroid[0],
                                             radians(lat1), radians(lon1),
                                             radians(lat2), radians(lon2))

        return (x_size, y_size)

    def get_all_pixelsize_metres(self):
        """
        Compute the size (in metres) of each pixel in this GriddedGeoBox. \
        Only one longitude column is returned from which as all other pixel \
        sizes can be derived.

        :return:
            An array of tuples (x_size, y_size) each gives the size of one pixel \
            in metres beginning with the pixel at the NW corner and extending to the \
            one at SW corner.
        """
        result = []
        for y_val in range(0, self.shape[1]):
            result.append(self.get_pixelsize_metres(xy=(0, y_val)))
        return result

    @property
    def ul(self):
        """
        Return the upper left corner co-ordinate in the units
        defined by the geobox's co-ordinate reference frame.
        """
        return self.origin

    @property
    def ur(self):
        """
        Return the upper right corner co-ordinate in the units
        defined by the geobox's co-ordinate reference frame.
        """
        ur = self.convert_coordinates((self.shape[1], 0))
        return ur

    @property
    def lr(self):
        """
        Return the lower right corner co-ordinate in the units
        defined by the geobox's co-ordinate reference frame.
        """
        return self.corner

    @property
    def ll(self):
        """
        Return the lower left corner co-ordinate in the units
        defined by the geobox's co-ordinate reference frame.
        """
        ll = self.convert_coordinates((0, self.shape[0]))
        return ll

    @property
    def centre(self):
        """
        Return the centre co-ordinate in the units defined by the
        geobox's co-ordinate reference frame.
        """
        x = (self.ul[0] + self.ur[0] + self.lr[0] + self.ll[0]) / 4.0
        y = (self.ul[1] + self.ur[1] + self.lr[1] + self.ll[1]) / 4.0
        return (x, y)

    @property
    def ul_lonlat(self):
        """
        Return the upper left corner co-ordinate in geographical
        longitude and latitude degrees based on the WGS84 datum.
        """
        sr = osr.SpatialReference()
        sr.SetFromUserInput(CRS)
        ul = self.transform_coordinates(self.origin, sr)
        return ul

    @property
    def ur_lonlat(self):
        """
        Return the upper right corner co-ordinate in geographical
        longitude and latitude degrees based on the WGS84 datum.
        """
        sr = osr.SpatialReference()
        sr.SetFromUserInput(CRS)
        ur = self.transform_coordinates(self.ur, sr)
        return ur

    @property
    def lr_lonlat(self):
        """
        Return the lower right corner co-ordinate in geographical
        longitude and latitude degrees based on the WGS84 datum.
        """
        sr = osr.SpatialReference()
        sr.SetFromUserInput(CRS)
        lr = self.transform_coordinates(self.corner, sr)
        return lr

    @property
    def ll_lonlat(self):
        """
        Return the lower left corner co-ordinate in geographical
        longitude and latitude degrees based on the WGS84 datum.
        """
        sr = osr.SpatialReference()
        sr.SetFromUserInput(CRS)
        ll = self.transform_coordinates(self.ll, sr)
        return ll

    @property
    def centre_lonlat(self):
        """
        Return the centre co-ordinate in geographical longitude
        and latitude degrees based on the WGS84 datum.
        """
        sr = osr.SpatialReference()
        sr.SetFromUserInput(CRS)
        centre = self.transform_coordinates(self.centre, sr)
        return centre

