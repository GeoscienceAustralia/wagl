#!/usr/bin/env python

def write_img(array, name='', format='ENVI', projection=None, geotransform=None):
    """
    Just a small function to write float32 2D arrays.
    """

    dims = array.shape
    dtype = 6
    bands = 1

    drv = gdal.GetDriverByName(format)
    outds = drv.Create(name, dims[1], dims[0], bands, dtype)

    if prj:
        outds.SetProjection(prj)

    if geot:
        outds.SetGeoTransform(geot)

    band = outds.GetRasterBand(1)
    band.WriteArray(array)
    band.FlushCache()

    outds = None

def read_img(fname):
    """
    
    """

    ds = gdal.Open(fname)

    img = ds.ReadAsArray()

    return img

def find_file(dir, file):
    """
    
    """
    fname = os.path.join(dir, file)
    if os.path.isfile(fname):
        return fname
    else:
        print "Error! file not found: %s" %fname
        sys.exit(2)

