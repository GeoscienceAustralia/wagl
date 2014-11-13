"""
Gridded Data
"""
import gdal
import rasterio as rio
import math
import osr

from affine import Affine

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

    :author: Steven Ring, Oct 2014
    """

    @staticmethod
    def from_dataset(dataset):
        """
        Return the GriddedGeoBox that encloses the full extent of 
        the supplied Rasterio or GDAL dataset

        Arguments:
            dataset: an open rasterio or GDAL Dataset

        Returns:
            A GriddedGeoBox that encloses the full extent of
            the supplied dataset
        """
        if isinstance(dataset, gdal.Dataset):
            return GriddedGeoBox.from_gdal_dataset(dataset)
        elif isinstance(dataset, rio._base.DatasetReader):
            return GriddedGeoBox.from_rio_dataset(dataset)
        else:
            raise ValueError("GriddedGeoBox.from_dataset() expects"\
                " GDAL or rasterio dataset, not %s" % (type(dataset))) 

    @staticmethod
    def from_rio_dataset(dataset):
        """
        Return the GriddedGeoBox that encloses the full extent of 
        the supplied Rasterio dataset

        Arguments:
            dataset: an open rasterio Dataset

        Returns:
            A GriddedGeoBox that encloses the full extent of
            the supplied dataset
        """
        bbshape = dataset.shape
        origin = (dataset.affine[2], dataset.affine[5])
        pixelsize = dataset.res
        crsName = str(dataset.crs[u'init'])

        return GriddedGeoBox(bbshape, origin, pixelsize, crsName)

    @staticmethod
    def from_gdal_dataset(dataset):
        """
        Return the GriddedGeoBox that encloses the full extent of 
        the supplied GDAL dataset

        Arguments:
            dataset: an open GDAL Dataset

        Returns:
            A GriddedGeoBox that encloses the full extent of
            the supplied dataset
        """
        bbshape = (dataset.RasterYSize, dataset.RasterXSize)
        transform = dataset.GetGeoTransform()
        origin = (transform[0], transform[3])
        pixelsize = (abs(transform[1]), abs(transform[5]))
        crsString = str(dataset.GetProjection())

	return GriddedGeoBox(bbshape, origin, pixelsize, crsString)

    @staticmethod
    def from_corners(origin, corner,  pixelsize=(0.00025,0.00025), crs='EPSG:4326'):
        """
        Return a GriddedGeoBox defined by the the two supplied corners

        Arguments:
            origin: a tuple (in CRS coordinates) representing the positon
                    of the NW corner of the GriddedGeoBox 
            corner: a tuple (in CRS coordinates) representing the positon
                    of the NW corner of the GriddedGeoBox
            pixelsize: pixel (xSize, ySize) in CRS units
            crs: the SpatialReferenceSystem in which both origin and pixelsize
                are expressed (supports various text formats)
        """
        a = Affine(pixelsize[0], 0, origin[0], 0, -pixelsize[1], origin[1])
        shapeXY = tuple([int(math.ceil(v)) for v in ~a * corner])
        shapeYX = (shapeXY[1], shapeXY[0])
        return GriddedGeoBox(shapeYX, origin, pixelsize, crs)

    def __init__(self, shape=(1, 1), origin=(0.0, 0.0), pixelsize=(0.00025,0.00025),
                 crs='EPSG:4326'):
        """
        Create a new GriddedGeoBox

        Arguments:
            shape: (ySize, xSize) 2-tuple defining the shape of the GGB
                   Note: Use getShapeXY() to get (xSize, ySize)
            origin: (xPos, yPos) 2-tuple difining the upper left GGB corner
            pixelsize: pixel (xSize, ySize) in CRS units
            crs: the SpatialReferenceSystem in which both origin and
                 pixelsize are expressed (supports various text formats)
        """
        self.pixelsize = pixelsize
        self.shape = tuple([int(v) for v in shape])
        # the origin must be aligned with the grid
        self.origin = tuple([
            int(origin[0] / self.pixelsize[0]) * self.pixelsize[0],
            int(origin[1] / self.pixelsize[1]) * self.pixelsize[1]
            ])
        self.crs = osr.SpatialReference()
        if self.crs == self.crs.SetFromUserInput(crs):
            raise ValueError("Invalid crs: %s" % (crs, ))
        self.affine = Affine(self.pixelsize[0], 0, self.origin[0], 0,
                             -self.pixelsize[1], self.origin[1])
        self.corner = self.affine * self.getShapeXY()

    def getShapeXY(self):
        return (self.shape[1], self.shape[0])

    def getShapeYX(self):
        return self.shape

    def copy(self, crs='EPSG:4326'):
        """
        Create a copy of this GriddedGeoBox transformed to the supplied
        Coordinate Reference System. The new GGB will have idential shape
        to the old and will be grid aligned to the new CRS. Pixel size 
        may change to accommodate the new CRS
        """
        def transformPoint(transformation, point):
            (x, y, z) = transformation.TransformPoint(point[0], point[1])
            return (x, y)

        newCrs = osr.SpatialReference()
        newCrs.SetFromUserInput(crs)
        old2New = osr.CoordinateTransformation(self.crs, newCrs)
        newOrigin = transformPoint(old2New, self.origin)
        newCorner = transformPoint(old2New, self.corner)
        newPixelSize = tuple([
            abs((self.origin[0]-newCorner[0])/self.shape[0]),
            abs((self.origin[1]-newCorner[1])/self.shape[0])
            ])

        return GriddedGeoBox(self.shape, newOrigin, newPixelSize, crs=crs)

    def __str__(self):
        return 'GriddedGeoBox(origin=%s,shape=%s,pixelsize=%s,crs: %s)' % (
                    self.origin, self.shape, str(self.pixelsize),
                    self.crs.ExportToProj4())

    def window(self, enclosedGGB):
        """
        Return the window, (expressed in the grid coordinates of this
        GriddedGeoBox) of the supplied GriddedGeoBox. Self must be a
        GriddedGeoBox which fully encloses the enclosed GGB.

        Arguments:
            enclosedGGB: a GGB which is a subset of this GGB and is fully
                        enclosed by it

        Returns:
            A pair of range tuples defining the enclosed retangular subset
            ((row_start, row_stop), (col_start, col_stop))
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
        return self.shape[1]

    def y_size(self):
        return self.shape[0]
