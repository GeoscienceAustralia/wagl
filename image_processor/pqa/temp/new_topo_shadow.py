#! /usr/bin/env python

import sys, os, numpy, gc, numexpr, ogr, osr, pdb, datetime
from scipy import ndimage
from osgeo import gdal

#import cast_shadow
#import topo2
from frtn_topo_shad import topo_shadow
#import frtn_topo_shad_v2

create_shapefile = False

def TopographicShadow(image, img_geoT, img_prj, DEM, metafile, lat, sensor_azi_angle, sensor_view_angle, solar_azi, solar_zen, bitpos=14):
    '''Creates a 2D array of topographic shadow.

       Executes "gdaldem hillshade" from the command line.

       Args:
           image: either a multiband or single band ndarray.
           img_geoT: The geo-transformation co-ordinates of the image.
           img_prj: The projection information of the image.
           DEM: A string file path to the location of the Digital Elevation
                Model that is to be used for surface topography.
           metafile: Either a full file path string name of the metadata file,
                     or a dictionary containing the relevant parameters
                     (Sun azimuth and elevation).
           bitpos: The bit position for specifying the resulting mask (Default
                   value is 14).

       Returns:
           An ndarray with 1 for no shadow and 0 for shadow specified by the
           bit position.

    '''

    def img2map(geoTransform, pixel):
        '''Converts a pixel (image) co-ordinate into a map co-ordinate.

        '''

        mapx = pixel[1] * geoTransform[1] + geoTransform[0]
        mapy = geoTransform[3] - (pixel[0] * (numpy.abs(geoTransform[5])))
        return (mapx,mapy)

    def map2img(geoTransform, location):
        '''Converts a map co-ordinate into a pixel (image) co-ordinate.

        '''

        imgx = int(numpy.round((location[0] - geoTransform[0])/geoTransform[1]))
        imgy = int(numpy.round((geoTransform[3] - location[1])/numpy.abs(geoTransform[5])))
        return (imgy,imgx)

    def map2img_array(mapx, mapy, geoTransform):
        '''Converts the x and y map locations in image co-ordinates.

           Rather than operating over single co-ordinates, this will operate
           over an array.

           Args:
               new_mapx: The x projected shadow location.
               new_mapy: The y projected shadow location.
               geoTransform: Is the Image co-ordinate information (upper left
                  coords, offset and pixel sizes)

          Returns:
               A tuple containing the indices of the locations.
        '''

        dict = { 'a': numpy.float32(geoTransform[0]), 'b': numpy.float32(geoTransform[1]) }
        imgx = numpy.round(numexpr.evaluate("(mapx - a)/b", dict, locals())).astype('int32')

        dict = { 'a': numpy.float32(geoTransform[3]), 'b': numpy.float32(abs(geoTransform[5])) }
        imgy = numpy.round(numexpr.evaluate("(a - mapy)/b", dict, locals())).astype('int32')

        return (imgy, imgx)


    # Returns the required line from a list of strings
    def linefinder(array, string = ""):
        '''Searches a list for the specified string.

           Args:
               array: A list containing searchable strings.
               string: User input containing the string to search.

           Returns:
               The line containing the found sting.
        '''

        for line in array:
            if string in str(line):
                return line

    # Reads the metadata file in order to extract the needed parameters
    def read_metafile(metafile):
        '''Opens the metadata file and extracs relevant parameters.

           Args:
               metafile: A full string path name to the metadata file.

           Returns:
               Dictionary containing the parameters.
        '''

        f         = open(metafile, 'r')
        met_array = f.readlines()
        f.close()

        sfind   = linefinder(met_array, 'SUN_AZIMUTH')
        s_azi   = float(sfind.split()[2])
        sfind   = linefinder(met_array, 'SUN_ELEVATION')
        s_elev  = float(sfind.split()[2])



        params = {
                     'Sun_Azimuth'      : s_azi,
                     'Sun_Elevation'    : s_elev,
                 }

        return params

    def slope_aspect(array, pxsize, pysize):
        '''Calculates the slope and aspect of an array.

           Args:
               array: A 2D numpy array, generally a DEM.
               pix_size: The size of a pixel, eg 25 for 25 metres.

           Returns:
               Two seperate 2D numpy arrays; the first containing the slope
               and the second containing the aspect.
        '''
        dzdx = numpy.zeros(array.shape, dtype='float32')
        dzdy = numpy.zeros(array.shape, dtype='float32')
        slp  = numpy.zeros(array.shape, dtype='float32')
        asp  = numpy.zeros(array.shape, dtype='float32')
        ndimage.sobel(array, axis=1, output=dzdx)
        ndimage.sobel(array, axis=0, output=dzdy)
        dzdx /= (8.*pxsize)
        dzdy /= (8.*pysize)
        numpy.arctan(numpy.hypot(dzdx,dzdy), out=slp)
        numpy.arctan2(dzdy, -dzdx, out=asp)
        return slp, asp


    def cal_pole(zenith, azimuth, slope, aspect):

        '''
          The zenith argument is not necessarily the zenith. The first call to this
          function will be solar_zenith, solar_azimuth, slope & aspect. The second
          call to this function is the sensor_view_angle, sensor_azimuth, slope &
          aspect.
        '''

        eps= 0.000001
        pi = numpy.pi
        d2r = pi/180
        r2d = 180/pi


        pdiff = azimuth - aspect

        # convert arrays to radians, do the calcs then convert back to degrees
        # (if needed) for the queries
        zenith *= d2r
        azimuth *= d2r
        slope *= d2r
        aspect *= d2r


        pdiff *= d2r

        costhp = numexpr.evaluate("cos(zenith) * cos(slope) + sin(zenith) * sin(slope) * cos(pdiff)")

        thp = numexpr.evaluate("arccos(costhp)") # convert to degrees later !!!


        slope *= r2d
        zenith *= r2d
        azimuth *= r2d
        aspect *= r2d

        thp   *= r2d

        query = numpy.abs(slope) <= eps
        thp[query] = zenith[query]

        query = costhp >= (1.0 - eps)
        thp[query] = 0.0

        thp *= d2r

        return thp
#--------------Processing Here-------------------------------

    st = datetime.datetime.now()
    d_box = []
    d_co_ords = []

    img_ref = osr.SpatialReference()
    dem_ref = osr.SpatialReference()
    img_ref.ImportFromWkt(img_prj)

    dims = image.shape
    if len(dims) >2:
        ncols = dims[2]
        nrows = dims[1]
        dims  = (nrows,ncols)

    # Fuqin states that in order to account for changes in elevation and

    # solar zenith angles for a scene; a buffer of 250 pixels should be
    # sufficient for Australian conditions. The 250 pixel buffer is based on a
    # scene resolution of 0.00025 degrees (approx 25m). This equates to 6.25km.
    # The scene's bounding co-ordinates will be added/subtracted accordingly.

    # Two options: get the bounding map co-ordinates then add/subtract
    # 6250m, 'OR' add/subtract 250 pixels. The latter will work for both
    # projected and geographic reference frames whereas the former will only
    # work for projected frames but it can account for images that don't have
    # a 25m resolution.
    # Will implement the latter method (add/subtract 250 pixels).

    d_box.append(img2map(geoTransform=img_geoT, pixel=(0-250,0-250))) # UL
    d_box.append(img2map(geoTransform=img_geoT, pixel=(0-250,dims[1]+250))) #UR
    d_box.append(img2map(geoTransform=img_geoT, pixel=(dims[0]+250,dims[1]+250))) #LR
    d_box.append(img2map(geoTransform=img_geoT, pixel=(dims[0]+250,0-250))) # LL

    print 'd_box: ', d_box
    sys.stdout.flush()

    for corner in d_box:
        d_co_ords.append(corner[0])
        d_co_ords.append(corner[1])

    print 'd_co_ords: ', d_co_ords
    print 'img_geoT: ', img_geoT
    sys.stdout.flush()

    if os.path.isfile(DEM):
        dem_obj = gdal.Open(DEM, gdal.gdalconst.GA_ReadOnly)
        assert dem_obj
        dem_prj  = dem_obj.GetProjection()
        dem_geoT = dem_obj.GetGeoTransform()
        dem_ref.ImportFromWkt(dem_prj)
        dem_cols = dem_obj.RasterXSize
        dem_rows = dem_obj.RasterYSize
        R = dem_ref.GetSemiMajor()
    else:
        raise Exception('DEM needs to be a string pathname to a valid file')

    # Retrieve the image bounding co-ords and create a vector geometry set
    if type(d_co_ords[0]) == int:
        d_wkt = 'MULTIPOINT(%d %d, %d %d, %d %d, %d %d)' %(d_co_ords[0],d_co_ords[1],d_co_ords[2],d_co_ords[3],d_co_ords[4],d_co_ords[5],d_co_ords[6],d_co_ords[7])
    else:
        d_wkt = 'MULTIPOINT(%f %f, %f %f, %f %f, %f %f)' %(d_co_ords[0],d_co_ords[1],d_co_ords[2],d_co_ords[3],d_co_ords[4],d_co_ords[5],d_co_ords[6],d_co_ords[7])

    if create_shapefile:
        if type(d_co_ords[0]) == int:
            dp_wkt = 'POLYGON((%d %d, %d %d, %d %d, %d %d))' %(d_co_ords[0],d_co_ords[1],d_co_ords[2],d_co_ords[3],d_co_ords[4],d_co_ords[5],d_co_ords[6],d_co_ords[7])
        else:
            dp_wkt = 'POLYGON((%f %f, %f %f, %f %f, %f %f))' %(d_co_ords[0],d_co_ords[1],d_co_ords[2],d_co_ords[3],d_co_ords[4],d_co_ords[5],d_co_ords[6],d_co_ords[7])

        d_out_name = '/short/v10/tmp/temp_stuff/topotest/Extended_DEM_coverage.shp'
        pobj = ogr.GetDriverByName('ESRI Shapefile').CreateDataSource(d_out_name)
        layer = pobj.CreateLayer('Extended_DEM_Cover',geom_type=ogr.wkbPolygon, srs=img_ref)
        feature = ogr.Feature(layer.GetLayerDefn())
        polygon = ogr.CreateGeometryFromWkt(dp_wkt)
        feature.SetGeometry(polygon)
        layer.CreateFeature(feature)
        pobj = None
        layer = None


    # Create the vector geometry set and transform the co-ords to match
    # the DEM file
    d_box_geom = ogr.CreateGeometryFromWkt(d_wkt)
    d_tform    = osr.CoordinateTransformation(img_ref, dem_ref)
    d_box_geom.Transform(d_tform)

    d_new_box = []
    for p in range(d_box_geom.GetGeometryCount()):
        d_point = d_box_geom.GetGeometryRef(p)
        d_new_box.append(d_point.GetPoint_2D())

    print 'd_new_box: ', d_new_box
    sys.stdout.flush()

    d_x = []
    d_y = []
    for c in d_new_box:
        d_x.append(c[0])
        d_y.append(c[1])

    d_xmin = numpy.min(d_x)
    d_xmax = numpy.max(d_x)
    d_ymin = numpy.min(d_y)
    d_ymax = numpy.max(d_y)

    # Retrieve the image co_ords of the DEM
    d_UL = map2img(geoTransform=dem_geoT, location=d_new_box[0])
    d_UR = map2img(geoTransform=dem_geoT, location=d_new_box[1])
    d_LR = map2img(geoTransform=dem_geoT, location=d_new_box[2])
    d_LL = map2img(geoTransform=dem_geoT, location=d_new_box[3])

    # Compute the offsets in order to read only the portion of the DEM that
    # covers the extents of the image file.
    d_ix    = numpy.array([d_UL[1],d_UR[1],d_LR[1],d_LL[1]])
    d_iy    = numpy.array([d_UL[0],d_UR[0],d_LR[0],d_LL[0]])
    d_ixmin = int(numpy.min(d_ix))
    d_ixmax = int(numpy.max(d_ix))
    d_iymin = int(numpy.min(d_iy))
    d_iymax = int(numpy.max(d_iy))
    d_xoff  = d_ixmin
    d_yoff  = d_iymin
    d_xsize = d_ixmax - d_ixmin
    d_ysize = d_iymax - d_iymin



    # TODO if extents go outside the DEM, get another DEM.
    # Test if the number of columns and rows to read exceeds the DEM extents.
    if (d_xoff > dem_cols) | (d_yoff > dem_rows):
        return 'Topo Shadow not performed; Image has no assocciated DEM.'
    elif d_xoff + d_xsize > dem_cols:
        #if xoff + xsize > dem_cols:
            #xsize = dem_cols - xoff
            #return 'Topo Shadow not performed; Image has no assocciated DEM.'
        return 'Topo Shadow not performed; Image has no assocciated DEM.'
    elif d_yoff + d_ysize > dem_rows:
        #if yoff + ysize > dem_rows:
            #ysize = dem_rows - yoff
            #return 'Topo Shadow not performed; Image has no assocciated DEM.'
        return 'Topo Shadow not performed; Image has no assocciated DEM.'
    else:

        # Read the data when creating the gdal image memory object
        dem_subs_geoT = (d_xmin,dem_geoT[1],0.0,d_ymax,0.0,dem_geoT[5])

        # As we're now dependent upon the sensor view angle, and other grid
        # angle information stuff, we have to project the DEM to the landsat
        # space, ie projected to metres.
        # In order to properly project the larger DEM extent to a space that
        # will match (pixel wise) too a smaller landsat extent, we need to
        # create a dummy array, with a dummy geoTransform variable.

        # The d_box variable contains the buffered extents
        buff_geoT = (d_box[0][0], img_geoT[1], img_geoT[2], d_box[0][1], img_geoT[4], img_geoT[5])

        et = datetime.datetime.now()

        print 'time to retrieve the dem: ', et - st
        print 'check for 6.25km buffer!'
        print 'buff_geoT: ', buff_geoT
        print 'img_geoT: ', img_geoT
        sys.stdout.flush()

        st = datetime.datetime.now()

        memdriver = gdal.GetDriverByName("MEM")
        memdem = memdriver.Create("", d_xsize, d_ysize, 1, gdal.GDT_Float32)
        memdem.SetGeoTransform(dem_subs_geoT)
        memdem.SetProjection(dem_prj)
        outband = memdem.GetRasterBand(1)
        outband.SetNoDataValue(outband.GetNoDataValue())
        outband.WriteArray(dem_obj.ReadAsArray(d_xoff, d_yoff, d_xsize, d_ysize))

        NoDataValue = outband.GetNoDataValue()

        # at this stage the buffer is 250 pixels on either side of the array.
        # This could change!!!
        outds = memdriver.Create("",  dims[1]+500, dims[0]+500, 1, gdal.GDT_Float32)
        outds.SetGeoTransform(buff_geoT)
        outds.SetProjection(img_prj)

        proj = gdal.ReprojectImage(memdem, outds, None, None, gdal.GRA_Bilinear)
        memdem = None
        #del memdem, dem_arr; gc.collect()
        del memdem; gc.collect()
        prjDEM = outds.ReadAsArray()
        outds = None
        del outds; gc.collect()

        et = datetime.datetime.now()
        print 'time to re-project the dem: ', et - st
        sys.stdout.flush()

        l_box = []
        l_box.append(img2map(geoTransform=img_geoT, pixel=(0,0))) # UL
        l_box.append(img2map(geoTransform=img_geoT, pixel=(0,dims[1]))) #UR
        l_box.append(img2map(geoTransform=img_geoT, pixel=(dims[0],dims[1]))) #LR
        l_box.append(img2map(geoTransform=img_geoT, pixel=(dims[0],0))) # LL

        # Calculate offsets between the DEM and the landsat scene
        sub_UL = map2img(geoTransform=buff_geoT, location=l_box[0])
        sub_UR = map2img(geoTransform=buff_geoT, location=l_box[1])
        sub_LR = map2img(geoTransform=buff_geoT, location=l_box[2])
        sub_LL = map2img(geoTransform=buff_geoT, location=l_box[3])

        x_off = sub_UL[1]
        y_off = sub_UL[0]
        x_end = sub_LR[1]
        y_end = sub_LR[0]

        # In order to index the DEM for the landsat extents
        # +1 as numpy upper extent is 'up to but not including end index'
        #prjDEM_L_subs = prjDEM[y_off:y_end+1, x_off:x_end+1]
        prjDEM_L_subs = prjDEM[y_off:y_end, x_off:x_end]

        '''
        # Write out the projected DEM
        driver = gdal.GetDriverByName("ENVI")
        outds  = driver.Create('projected_dem', dims[1], dims[0],1, gdal.GDT_Float32)
        outds.SetGeoTransform(img_geoT)
        outds.SetProjection(img_prj)
        outband = outds.GetRasterBand(1)
        outband.WriteArray(prjDEM_L_subs)
        outds = None
        '''

        # Lat/lon grids will now be generated from the same code as that which
        # is in the NBAR algorithm.


        # If the landsat dataset is not in a geographic space, ie in lat/lon
        # then the pixel size is contant and retrieved from the geoTransform
        # variable. The cast shadow algorithm is expecting an array though.
        # Solution is to just create an array of constant pixel sizes.
        # There may be a case of a sensor having non-uniform x and y pixel
        # sizes. Solution, check if x_pix_size == y_pix_size (ABS needed).
        # If they're the same just create one array and feed it in twice to the
        # cast shadow function, else create two arrays, which would be the
        # default for lat/lon based projection (i.e. geographics).

        # Calculate arrays of pixel size in metres at respective lat/lon
        # If lat/lon grids are needed as radians, then convert from the start
        # and leave as is. FIND OUT !!!
        # Probably implement this as a function

        pi = numpy.float32(numpy.pi)
        d2r = pi/180
        r2d = 180/pi

        print 'img_ref.IsGeographic: ', img_ref.IsGeographic()

        if (img_ref.IsGeographic() == 1):

            lat *= d2r
            SemiMajor = dem_ref.GetSemiMajor()
            SemiMinor = dem_ref.GetSemiMinor()
            InvFlattening = dem_ref.GetInvFlattening()

            print 'InvFlattening: ', InvFlattening

            bb = SemiMajor * (1 - 1/InvFlattening)
            print 'bb: ', bb
            print 'SemiMinor: ', SemiMinor
            sys.stdout.flush()

            cc = numexpr.evaluate("SemiMajor * cos(lat)")
            dd = numexpr.evaluate("SemiMinor * sin(lat)")
            rr = numexpr.evaluate("sqrt((SemiMajor**2 * cc**2 + SemiMinor**2 * dd**2) / (cc**2 + dd**2))")
            ddy = (numpy.abs(img_geoT[5])) * d2r
            ddx = numexpr.evaluate("arccos(sin(lat)**2 + cos(lat)**2 * cos(ddy))")
            dy = numexpr.evaluate("rr * ddy")
            dx = numexpr.evaluate("rr * ddx")


            del cc, dd, rr, lat, rad_lat_grid; gc.collect()
            sys.stdout.flush()

        else:
            dy = numpy.abs(img_geoT[5])
            dx = img_geoT[1]
            del lat; gc.collect()


        mask = numpy.ones(dims, dtype='int8')

        # Self shadow algorithm !!!!!!!!

        '''
        This routine needs the solar azimuth, view azimuth, solar zenith, view
        zenith angles.  Will need to take the code from NBAR to duplicate these
        rasters. As these rasters are calculated for the extents of the landsat
        raster; The DEM will need to be projected to the landsat grid. However
        as a buffer is needed for the cast shadow calculation, will need to
        create a dummy image that contains the buffer. The projection info will
        be extracted from the landsat file, but the geotransform info will be
        created. This dummy image will be used as a base from which the DEM
        will be projected to. This should ensure that the DEM will match up
        pixel by pixel with the landsat file.

        Potentially, the DEM will have No-Data values, eg -32767.00 as is with
        the current DEM. Specifying the no-data value will allow GDAL to ignore
        it in the projection process. The value will be returned as a zero,
        even when specifying the no-data value during the creation of the
        projected object (as 13/08/2012). This maybe a bug within GDAL.

        The algorithms created by Fuqin don't really account for negative
        values. So for large negative values, the algorithm could be searching          outside the buffer extent. Mostly only a problem for negative null
        values.
        '''


        # Fortran input order (not including header and DEM)
        # inputs required: solar_zenith, sensor_view_angle, view_azimuth, solar_azimuth
        # fortran names: solar, view, solar_azimuth, view_azimuth

        # Python input (different order to above)
        #          solar_zenith, solar_azimuth, sensor_view_angle, view_azimuth
        # named inputs: solar_zen, solar_azi, sensor_view_angle, sensor_view_azi

        print 'prjDEM.shape: ', prjDEM.shape
        print 'y_off, y_end ', y_off, y_end
        print 'x_off, x_end ', x_off, x_end
        print 'sub_UL: ', sub_UL
        print 'sub_UR: ', sub_UR
        print 'sub_LR: ', sub_LR
        print 'sub_LL: ', sub_LL

        dem_1pix_buf = prjDEM[y_off-1:y_end+1, x_off-1:x_end+1]
        print 'dem_1pix_buf.shape: ', dem_1pix_buf.shape
        sys.stdout.flush()

        st = datetime.datetime.now()

        # Check if the landsat image is in geographics
        if (img_ref.IsGeographic() == 1):
            # need to modify the pixel size array (add a 1 pixel border)
            # We could leave the numpy to do a reflect at the border, for the
            # DEM, but some may consider this an inaccurate result. The pixel
            # sizes shouldn't change much (if at all) per lattitude for a
            # single pixel at the edge, so we'll reflect it at the border to
            # match the DEM extent. Pixel sizes will differ over 1 pix if the
            # change in lattitude is great, i.e. big pixels.
            new_pxsize = numpy.zeros(dem_1pix_buf.shape, dtype='float32')
            new_pysize = numpy.zeros(dem_1pix_buf.shape, dtype='float32')
            new_pxsize[1:-1,1:-1] = dx
            new_pysize[1:-1,1:-1] = dy
            new_pxsize[0,1:-1] = dx[0,:] # 1st row
            new_pysize[0,1:-1] = dy[0,:]
            new_pxsize[1:-1,0] = dx[:,0] # 1st col
            new_pysize[1:-1,0] = dy[:,0]
            new_pxsize[1:-1,-1] = dx[:,-1] # last col
            new_pysize[1:-1,-1] = dy[:,-1]
            new_pxsize[-1,1:-1] = dx[-1,:] # last row
            new_pysize[-1,1:-1] = dy[-1,:]
            new_pxsize[0,0] = dx[0,0] # UL
            new_pysize[0,0] = dy[0,0]
            new_pxsize[-1,-1] = dx[0,-1] # UR
            new_pysize[0,-1] = dy[0,-1]
            new_pxsize[-1,-1] = dx[-1,-1] # LR
            new_pysize[-1,-1] = dy[-1,-1]
            new_pxsize[-1,0] = dx[-1,0] # LL
            new_pysize[-1,0] = dy[-1,0]

            slope, aspect = slope_aspect(dem_1pix_buf, pxsize=new_pxsize, pysize=new_pysize)
            del new_pxsize; gc.collect()

        else:
            # The dataset is projected, and the dx and dy are single values
            slope, aspect = slope_aspect(dem_1pix_buf, pxsize=dx, pysize=dy)


        # get landsat_dims
        # subset the slope and aspect rasters based on the landsat extents
        print 'dims: ', dims
        print 'slope.shape: ', slope.shape
        print 'aspect.shape: ', aspect.shape
        sys.stdout.flush()
        slope = slope[1:-1,1:-1]
        aspect = aspect[1:-1,1:-1]
        print 'remove 1 pix buffer'
        print 'slope.shape: ', slope.shape
        print 'aspect.shape: ', aspect.shape
        sys.stdout.flush()
        slope *= r2d
        aspect *= r2d

        et = datetime.datetime.now()
        print 'slope & aspect calc time: ', et - st

        st = datetime.datetime.now()
        print 'calculating incident angles'
        sys.stdout.flush()

        incident_t = cal_pole(solar_zen, solar_azi, slope, aspect)

        et = datetime.datetime.now()
        print 'incident angle calc time: ', et - st

        st = datetime.datetime.now()
        print 'calculating exiting angles'
        sys.stdout.flush()

        exiting_t = cal_pole(sensor_view_angle, sensor_azi_angle, slope, aspect)

        et = datetime.datetime.now()
        print 'exiting angle calc time: ', et - st

        st = datetime.datetime.now()
        print 'finding self shadow'
        sys.stdout.flush()


        query = numpy.cos(incident_t) <= 0.0
        mask[query] = 0

        query = numpy.cos(exiting_t) <= 0.0
        mask[query] = 0

        et = datetime.datetime.now()
        print 'find self shadow calc time: ', et - st

        print 'num self shadow pixels: ', numpy.sum(mask[mask == 0])
        query = mask < 1
        sum = numpy.sum(mask[query])
        print 'num self shadow pixels: ', sum

        print 'writing out the self shadow mask'
        sys.stdout.flush()

        #pdb.set_trace()

        # Write out the self shadow mask
        driver = gdal.GetDriverByName("ENVI")
        outds  = driver.Create('self_shadow', mask.shape[1], mask.shape[0],1, gdal.GDT_Byte)
        outds.SetGeoTransform(img_geoT)
        outds.SetProjection(img_prj)
        outband = outds.GetRasterBand(1)
        outband.WriteArray(mask)
        outds = None

        del slope, aspect, sensor_azi_angle, sensor_view_angle; gc.collect()
        del incident_t; gc.collect()
        del dem_1pix_buf; gc.collect()


        # Cast shadow algorithm!!!!!!
        print 'starting cast shadow'
        sys.stdout.flush()


        # zenith + 3 degrees
        solar_zen += 3
        solar_zen = numpy.radians(solar_zen)


        # max height in DEM
        zmax = numpy.max(prjDEM)

        k_setting = 1500
        htol = 1.0

        # Check for angles > 360 and < 0
        query = (solar_azi > 360.0) | (solar_azi < 0.0)
        if (query.sum() > 0):
           return 'Error: Azimuth must be in 0 to 360 degrees!'

        # Accounting for different azimuthal cases
        # Case 1
        az_case = numpy.zeros(dims, dtype='int8')
        solar_azi_rad = numpy.zeros(dims, dtype='float32')
        query = (solar_azi >= 0.0) & (solar_azi <= 90.0)
        solar_azi_rad[query] = pi/2 - numpy.radians(solar_azi[query])
        az_case[query] = 1
        # Case 2
        query = (solar_azi >= 270.0) & (solar_azi <= 360.0)
        solar_azi_rad[query] = numpy.radians(solar_azi[query]) - (3*pi/2)
        az_case[query] = 2
        # Case 3
        query = (solar_azi >= 180.0) & (solar_azi <= 270.0)
        solar_azi_rad[query] = (3*pi/2) - numpy.radians(solar_azi[query])
        az_case[query] = 3
        # Case 4
        query = (solar_azi >= 90.0) & (solar_azi <= 180.0)
        solar_azi_rad[query] = numpy.radians(solar_azi[query]) - pi/2
        az_case[query] = 4

        sinphc = numexpr.evaluate("sin(solar_azi_rad)")
        cosphc = numexpr.evaluate("cos(solar_azi_rad)")


        # d is the max plane distance from the pixel that will occlude a pixel
        # will create an array of distances
        d = numexpr.evaluate("(zmax - prjDEM_L_subs)/tan(pi/2 - solar_zen)")

        # Two cases for d0. Will step in the larger of the x,y direction.
        # d0 is the basic step in the d direction that increments
        # 0.5 pixel in the larger of the two directions

        d0 = numpy.zeros(dims, dtype='float32')
        sinphc_x = sinphc * dx
        cosphc_y = cosphc * dy
        query = (cosphc_y > sinphc_x)
        del sinphc_x, cosphc_y; gc.collect()

        # Check whether to run with arrays of pixels sizes
        if (img_ref.IsGeographic() == 1):
            # Case 1
            d0[query] = (0.5*dx[query])/cosphc[query]
            # Case 2
            d0[~(query)] = (0.5*dy[~(query)])/sinphc[~(query)]
        else:
            # Case 1
            d0[query] = (0.5*dx)/cosphc[query]
            # Case 2
            d0[~(query)] = (0.5*dy)/sinphc[~(query)]


        # k_max is the number of steps to get out to d
        # n_inc(k) and m_inc(k) are the sample and line increments
        # h_offset(k) on exit is the altitude above the start of the
        # sun beam as you move towards the sun

        #k_max = ((d/d0 + 0.5) + 1).astype('int16')
        k_max = numexpr.evaluate("(d / d0 + 0.5) + 1)")
        k_max = k_max.astype('int16')

        # Test for steps > k_setting (set at 1500)
        query = k_max > k_setting
        if (query.sum() > 0):
            return 'Error: One or more pixels has a maximum distance larger than k_setting!'

        # the increments to search for occlusions have sign depending
        # on the azimuth case
        # The array size of the increments can be up to 1500, (hard coded in
        # original code).

        del query; gc.collect()

        '''
        Have re-written the fortran module.
        The following is a list of what is needed.
        DEM, dem_cols, dem_rows, sub_col, sub_row, mask, dy, dx, zenith,
        lat, lon, k_max, az_case, d0, h_offset, pi, zmax, l_xoff, l_yoff,
        cosphc, sinphc.

        Arrays will need to be converted to a fortran contiguous array.
        DEM, mask, dy, dx, zenith, lat, lon, k_max, az_case, d0, cosphc, sinphc.
        DEM, az_case, mask, k_max are integer arrays.
        cosphc, sinphc, zenith, dx, dy, lat, lon, d0 are float32.
        pi, zmax are float32.
        dem_cols, dem_rows, sub_col, sub_row, l_xoff, l_yoff are integers.
        '''

        print 'Converting to fortran contiguous arrays'
        print 'zmax: ', zmax
        print 'type(zmax): ', type(zmax)
        #print 'dem_dims: ', dem_dims
        sys.stdout.flush()

        x_off = sub_UL[1]
        y_off = sub_UL[0]

        # Convert the az_case to 1 or 2. 1 = case 1 or 4. 2 = case 2 or 3.
        # Just helps set up the loops for the fortran code
        az_case[az_case == 4] = 1
        az_case[az_case == 3] = 2

        st = datetime.datetime.now()


        # Check for geographic.
        # The only difference is geographic will have x & y pixel sizes
        # as arrays, whereas the projected version will have single values.
        if (img_ref.IsGeographic() == 1):

            # Convert arrays to fortran contiguous arrays
            prjDEM = numpy.asfortranarray(prjDEM)
            mask = numpy.asfortranarray(mask)
            dy = numpy.asfortranarray(dy)
            dx = numpy.asfortranarray(dx)
            solar_zen = numpy.asfortranarray(solar_zen)
            az_case = numpy.asfortranarray(az_case)
            cosphc = numpy.asfortranarray(cosphc)
            sinphc = numpy.asfortranarray(sinphc)
            k_max = numpy.asfortranarray(k_max)
            d0 = numpy.asfortranarray(d0)

            et = datetime.datetime.now()
            print 'time to convert to fortran contiguous: ', et - st
            sys.stdout.flush()
            st = datetime.datetime.now()

            topo_shadow.cast_shadow_geo(prjDEM, mask, dy, dx, solar_zen, k_max, az_case, d0, pi, zmax, x_off, y_off, cosphc, sinphc)

        else:
            # Convert arrays to fortran contiguous arrays
            prjDEM = numpy.asfortranarray(prjDEM)
            mask = numpy.asfortranarray(mask)
            solar_zen = numpy.asfortranarray(solar_zen)
            az_case = numpy.asfortranarray(az_case)
            cosphc = numpy.asfortranarray(cosphc)
            sinphc = numpy.asfortranarray(sinphc)
            k_max = numpy.asfortranarray(k_max)
            d0 = numpy.asfortranarray(d0)

            print 'prjDEM.dtype: ', prjDEM.dtype
            print 'mask.dtype: ', mask.dtype
            print 'solar_zen.dtype: ', solar_zen.dtype
            print 'az_case.dtype: ', az_case.dtype
            print 'cosphc.dtype: ', cosphc.dtype
            print 'sinphc.dtype: ', sinphc.dtype
            print 'k_max.dtype: ', k_max.dtype
            print 'd0.dtype: ', d0.dtype

            et = datetime.datetime.now()
            print 'time to convert to fortran contiguous: ', et - st
            sys.stdout.flush()
            st = datetime.datetime.now()

            topo_shadow.cast_shadow_prj(prjDEM, mask, dy, dx, solar_zen, k_max, az_case, d0, pi, zmax, x_off, y_off, cosphc, sinphc)

            dem_x = prjDEM.shape[1]
            dem_y = prjDEM.shape[0]
            msk_x = mask.shape[1]
            msk_y = mask.shape[0]

            #frtn_topo_shad_v2.cast_shadow_prj(prjDEM, mask, dy, dx, solar_zen, k_max, az_case, d0, pi, zmax, x_off, y_off, cosphc, sinphc, dem_x, dem_y, msk_x, msk_y)
            #frtn_topo_shad_v2.cast_shadow_prj(prjDEM, mask, dy, dx, solar_zen, k_max, az_case, d0, pi, zmax, x_off, y_off, cosphc, sinphc)



        #pdb.set_trace()

        et = datetime.datetime.now()
        print 'fortran cast shadow calc time: ', et - st


        mask = numpy.array(mask, order='C')

        #return topo_shad
    return mask


# The following is the original Fortran code for the cast shadow algorithm.
# Filename: shade_main_landsat_pixel.f
'''
      program shade_main
c
c     Program to calculate cast shadow for a standard Landsat scene
c     the program was originally written by DLB Jupp in Oct. 2010
c     for a small sub_matrix and was modified by Fuqin Li in Oct.
c     2010 so that the program can be used for large landsat scene.
c
c     Basically, a sub-matrix A is embedded in a larger DEM image
c     and the borders must be large enough to find the shaded pixels.
c     If we assume the solar azimuth and zenith angles change very
c     little within the sub-matrix A, then the Landsat scene can be
c     divided into several sub_matrix.
c     For Australian region, with 0.00025 degree resolution, the
c     sub-marix A is set to 500x500
c
c     we also need to set extra DEM lines/columns to run the Landsat
c     scene. This will change with elevation difference within the
c     scene and solar zenith angle. For Australian region and Landsat
c     scene with 0.00025 degree resolution, the maximum extra lines
c     are set to 250 pixels/lines for each direction. This figure
c     shold be sufficient for everywhere and anytime in Australia.
c     thus the DEM image will be larger than landsat image for
c     500 lines x 500 columns
c
c     In this program, these satellite related parameters are
c     in the input file read below. For other satellite, they can
c     be set differently.
c
c     Current program operates in all 4 Azimuth cases
c     (the four quadrants)
c
c     command line shade_main_landsat <input.dat> <dem.dat> <solar_angle.bin>
c                                     <sazi_angle.bin> <castshadow.img>

      parameter (k_setting=1500)

	integer ns,nl,nchf,nrow,ncol
c     NOTE: n_inc and m_inc are Floating Point arrays
	real n_inc(k_setting),m_inc(k_setting),hx,hy
	integer Aoff_x1, Aoff_x2,Aoff_y1,Aoff_y2
        integer nsA,nlA,nlA_ori,nsA_ori
	real h_offset(k_setting),zmax,zmin
        real solar(k_setting,k_setting*10)
        real sazi(k_setting,k_setting*10)
	real phi_sun,sun_zen,htol
	real a(k_setting,k_setting),dem(k_setting,k_setting*10)
        double precision alon(k_setting*10),alat(k_setting)
        logical exists
        double precision dres
        double precision alat1,alon1
	integer*2 mask(k_setting,k_setting)
        integer*2 mask_all(k_setting,k_setting*10)
	character*128 fname
c
      real pi,r2d,d2r,eps
      common/base/pi,r2d,d2r,eps

      if (IARGC() .ne. 5) then
      write(*,*) 'Error: Required parameters not specified properly!'
      stop 11
      endif
c
c     set basic constants
c
        pi=4.0*atan(1.0)
	r2d=180.0/pi
	d2r=pi/180.0
	eps=1.0e-7
c      open input file. the input file include some satellite
c      related parameters
c
       call GETARG(1, fname)
       open (1,file=fname,status='old')
c--------------------------------------------
c      read lines and columns for DEM image
       read(1,*)nl,ns
c--------------------------------------------
c      read lines and columns for Landsat image
       read(1,*) nrow,ncol
c------------------------------------------------
c      read Aoff_x1 and Aoff_x2 where Aoff_x1 is the last pixel
c      number before the Landsat image and Aoff_x2 is the first pixel
c      number after teh Landsat image
       read(1,*)Aoff_x1,Aoff_x2
c-------------------------------------------------
c      read Aoff_y1 and Aoff_y2 where Aoff_y1 is the last line
c      number before the Landsat image and Aoff_y2 is the first line
c      number after teh Landsat image
       read(1,*)Aoff_y1,Aoff_y2
c------------------------------------------------------
c      read sub-matrix lines and columns
       read(1,*)nlA_ori, nsA_ori
c----------------------------------------------
c      upple left latitude and lonitude for Landsat (in degree)
       read(1,*)alat1,alon1
c-----------------------------------------
c      read spatial resolution (in degree)
       read(1,*)dres
c
c     write a message to the console
c
      write(0,*) 'Program Shade Main'
	write(0,*) ''
c
c     set the tolerance for occlusion in metres
c     (usually >0 and less than 10.0m)
c
	htol=1.0

c----------------------------------------------------------------------

c
c    open the DEM file and read in the data (assumed real)

25      call GETARG(2, fname)
      call stripbz(fname,nchf)
      open(7,file=fname(1:nchf),status='old',
     .  access='direct',form='unformatted',recl=4*ns)
c---------------------------------------------------------
c      open solar zenith angle file
       call GETARG(3, fname)
       call stripbz(fname,nchf)
       open(8,file=fname(1:nchf),status='old',
     .  access='direct',form='unformatted',recl=4*ncol)
c--------------------------------------------------------
c      open solar azimuth angle file
       call GETARG(4, fname)
       call stripbz(fname,nchf)
       open(9,file=fname(1:nchf),status='old',
     .  access='direct',form='unformatted',recl=4*ncol)
c--------------------------------------------------------
c     Now write out the Mask file.

      call GETARG(5, fname)
      call stripbz(fname,nchf)
      open(unit=55,file=fname(1:nchf),access='direct',
     . form='unformatted',recl=2*ncol)
c--------------------------------------------------------------

c      kky for line and kkx for column
c      kky and kkx are the sub_marix number
       kky=int(nrow/nlA_ori)
       kkx=int(ncol/nsA_ori)
c
c      calculate longitude for each pixel of teh line
        do j=1,ncol
        alon(j)=alon1+(j-1)*dres
        enddo

c
        do k=1, kky
c      calculate sub_marix DEM dimensions
        nlA=nlA_ori
        mmax_sub=nlA+Aoff_y1+Aoff_y2
c
c       read DEM data
	do i=1,mmax_sub
	  read(7,rec=(k-1)*nlA_ori+i)(dem(i,j),j=1,ns)
	enddo

        zmax=maxval(dem(1:mmax_sub,1:ns))
        zmin=minval(dem(1:mmax_sub,1:ns))
c
c       read solar zenith and azimuth angle
        do i=1,nlA
        read(8,rec=(k-1)*nlA_ori+i)(solar(i,j),j=1,ncol)
        read(9,rec=(k-1)*nlA_ori+i)(sazi(i,j),j=1,ncol)
        enddo

        print*,zmax,zmin,solar(500,8800),sazi(500,8800)
c
c       calculate latitude for each line
        do i=1,nlA
        alat(i)=alat1-((k-1)*nlA_ori+i-1)*dres
        enddo
        ii=nlA/2
        call pixelsize(alat(ii),dres,hx,hy)
c
c       divide seveal sub_matrix according to columns
c
        do l=1,kkx
        nsA=nsA_ori
        nmax_sub=nsA+Aoff_x1+Aoff_x2


       jj=(l-1)*nsA_ori+nsA/2
       write(98,*) k,l,ii,jj
       write(98,*)alat(ii),alon(jj)

       phi_sun=sazi(ii,jj)
c     NOTE zenith + 3 degrees
       sun_zen=solar(ii,jj)+3
       write(98,*)hx,hy,phi_sun,sun_zen

       do i=1,mmax_sub
       do j=1,nmax_sub
       a(i,j)=dem(i,(l-1)*nsA_ori+j)
       enddo
       enddo
       if (k .eq.7 .and. l .eq. 15) then
       print*,k,l,ii,jj,a(1,1),phi_sun,sun_zen,hx,hy
       print*,nsA,nlA,aoff_x1,aoff_y1
       endif

c
      call get_proj_shadows(hx,hy,nmax_sub,mmax_sub,
     . htol,phi_sun,sun_zen,zmax,zmin,
     .  a,mask,h_offset,n_inc,m_inc,Aoff_x1,Aoff_y1,nsA,nlA,
     .  k_setting)
       do i=1,nlA
       do j=1,nsA
       mask_all(i,(l-1)*nsA_ori+j)=mask(i,j)
       enddo
       enddo
c
       enddo

       if (ncol .gt. kkx*nsA_ori) then

       nsA=ncol-kkx*nsA_ori
       nmax_sub=nsA+Aoff_x1+Aoff_x2
       jj=kkx*nsA_ori+nsA/2

       phi_sun=sazi(ii,jj)
c     NOTE zenith + 3 degrees
       sun_zen=solar(ii,jj)+3

       do i=1,mmax_sub
       do j=1,nmax_sub
       a(i,j)=dem(i,kkx*nsA_ori+j)
       enddo
       enddo
       if (k .eq.7) then
       print*,k,a(1,1),phi_sun,sun_zen,hx,hy
       print*,nsA,nlA,Aoff_x1,Aoff_y1
       endif


c
      call get_proj_shadows(hx,hy,nmax_sub,mmax_sub,
     . htol,phi_sun,sun_zen,zmax,zmin,
     .  a,mask,h_offset,n_inc,m_inc,Aoff_x1,Aoff_y1,nsA,nlA,
     .  k_setting)
       do i=1,nlA
       do j=1,nsA
       mask_all(i,kkx*nsA_ori+j)=mask(i,j)
       enddo
       enddo
       endif

c----------------------------------------------------------
       do i=1,nlA
         write(unit=55,rec=(k-1)*nlA_ori+i)(mask_all(i,j),j=1,ncol)
         write(99,*)(mask_all(i,j),j=1,ncol)
        enddo
       enddo

       if (nrow .gt. kky*nlA_ori) then
        nlA=nrow-kky*nlA_ori
        mmax_sub=nlA+Aoff_y1+Aoff_y2
c
c       read DEM data
        do i=1,mmax_sub
          read(7,rec=kky*nlA_ori+i)(dem(i,j),j=1,ns)
        enddo

        zmax=maxval(dem(1:mmax_sub,1:ns))
        zmin=minval(dem(1:mmax_sub,1:ns))
        write(98,*) 'zmin, zmax in image=',zmin,zmax
c
c       read solar zenith and azimuth angle
        do i=1,nlA
        read(8,rec=kky*nlA_ori+i)(solar(i,j),j=1,ncol)
        read(9,rec=kky*nlA_ori+i)(sazi(i,j),j=1,ncol)
        enddo
        print*,zmax,zmin,solar(1,8800),sazi(1,8800)
c
c       calculate latitude and longitude for sub_matrix
        do i=1,nlA
        alat(i)=alat1-(kky*nlA_ori+i-1)*dres
        enddo

        ii=nlA/2

        call pixelsize(alat(ii),dres,hx,hy)

c
c       divide seveal sub_matrix according to columns
c
        do l=1,kkx
        nsA=nsA_ori
        nmax_sub=nsA+Aoff_x1+Aoff_x2

       jj=(l-1)*nsA_ori+nsA/2

       phi_sun=sazi(ii,jj)
c     NOTE zenith + 3 degrees
       sun_zen=solar(ii,jj)+3
       if (l.eq. 15) then
       print*,kky+1,l,phi_sun,sun_zen,hx,hy
       endif

       do i=1,mmax_sub
       do j=1,nmax_sub
       a(i,j)=dem(i,(l-1)*nsA_ori+j)
       enddo
       enddo
       if (l.eq. 2) then
       print*,kky+1,l,a(1,1),phi_sun, sun_zen
       endif
c
      call get_proj_shadows(hx,hy,nmax_sub,mmax_sub,
     . htol,phi_sun,sun_zen,zmax,zmin,
     .  a,mask,h_offset,n_inc,m_inc,Aoff_x1,Aoff_y1,nsA,nlA,
     .  k_setting)
       do i=1,nlA
       do j=1,nsA
       mask_all(i,(l-1)*nsA_ori+j)=mask(i,j)
       enddo
       enddo
c
       enddo

c----------------------------------------------------------
       if (ncol .gt. kkx*nsA_ori) then

       nsA=ncol-kkx*nsA_ori
       nmax_sub=nsA+Aoff_x1+Aoff_x2
       jj=kkx*nsA_ori+nsA/2

       phi_sun=sazi(ii,jj)
c     NOTE zenith + 3 degrees
       sun_zen=solar(ii,jj)+3
       print*,kky+1,kkx+1,phi_sun,sun_zen,hx,hy

       do i=1,mmax_sub
       do j=1,nmax_sub
       a(i,j)=dem(i,kkx*nsA_ori+j)
       enddo
       enddo
c
      call get_proj_shadows(hx,hy,nmax_sub,mmax_sub,
     . htol,phi_sun,sun_zen,zmax,zmin,
     .  a,mask,h_offset,n_inc,m_inc,Aoff_x1,Aoff_y1,nsA,nlA,
     .  k_setting)
       do i=1,nlA
       do j=1,nsA
       mask_all(i,kkx*nsA_ori+j)=mask(i,j)
       enddo
       enddo
       endif

c----------------------------------------------------------
       do i=1,nlA
         write(unit=55,rec=kky*nlA_ori+i)(mask_all(i,j),j=1,ncol)
        enddo
       endif
99      continue

	stop
	end
c-------------------------------------------------------------

      SUBROUTINE STRIPBz( CB, NCH)
C
c     microBRIAN routine with some modifications
c     stripbz strips off leading and trailing blanks
c     from a character string and appends a null after
c     the last non-blank character
c
c     the search is made to the length of the string
c     or the first null
c     the string is shifted to remove the leading blanks
c     nch is the number of non blank characters before the null
c
      CHARACTER*(*) CB
      INTEGER ICH, NCH
      NCH = 1
C
      DO 10 ICH = 1, LEN( CB)
         if (cb(ich:ich).eq.char(0)) then
             go to 15
         else IF (CB(ICH:ICH) .EQ. ' ' ) THEN
            IF (NCH .GT. 1 .AND. CB(NCH-1:NCH-1) .NE. ' ') THEN
               CB(NCH:NCH) = ' '
               NCH = NCH + 1
            ENDIF
         ELSE
            CB(NCH:NCH) = CB(ICH:ICH)
            NCH = NCH + 1
         ENDIF
   10 CONTINUE
15    NCH = NCH - 1
      IF (NCH .GT. 0 .AND. CB(NCH:NCH) .EQ. ' ') NCH = NCH - 1
      IF (NCH .le. 0) then
          nch=0
      endif
      cb(nch+1:nch+1)=char(0)
      RETURN
      END
c
c-----------------------------------------------------------------
       subroutine get_proj_shadows(hx,hy,ns,nl,htol,
     . phi_sun,sun_zen,zmax,zmin,a,mask,h_offset,n_inc,m_inc,
     .  aoff_x,aoff_y,nsA,nlA,k_setting)
c
        integer k_max,n_add,m_add,set_border,err,az_case
c     NOTE: n_inc and m_inc are Floating Point arrays
        real n_inc(k_setting),m_inc(k_setting)
        integer nchf,ncho,nmax_sub,mmax_sub,aoff_x,aoff_y,nsA,nlA
        integer t_aoff_x,t_aoff_y,t_nsA,t_nlA
        integer tval(8)
        real h_offset(k_setting)
        real phi_sun,zmax,zmin,sun_zen,hx,hy,htol
        real d,d0,a(k_setting,k_setting)
        logical status
        integer*2 mask(k_setting,k_setting),rmax,rmin

         real pi,r2d,d2r,eps
         common/base/pi,r2d,d2r,eps
c
c     calculate the border info for the sun position
c     In Australia and in the south in particular case=1
c     for Landsat since Landsat overpass is around 10:am
c     local time.
c     That is the sun azimuth is between East and North
c
      status=set_border(phi_sun,zmax,zmin,sun_zen,hx,hy,az_case,
     .        d,d0,k_max,h_offset,n_inc,m_inc,n_add,m_add,k_setting,
     .        k_setting,err)
c
      if (.not.status) then
          write(98,*) 'Error in set_border. Err=',err
          goto 99
      else
          write(98,*) 'Set_Border ran successfully'
        endif
c
c     define the maximum sized subset that can be processed
c
      if (az_case.eq.1) then
          t_aoff_x=0
          t_aoff_y=m_add
          t_nsA=ns-n_add
          t_nlA=nl-m_add
        else if (az_case.eq.2) then
          t_aoff_x=n_add
          t_aoff_y=m_add
          t_nsA=ns-n_add
          t_nlA=nl-m_add
        else if (az_case.eq.3) then
          t_aoff_x=n_add
          t_aoff_y=0
          t_nsA=ns-n_add
          t_nlA=nl-m_add
        else if (az_case.eq.4) then
          t_aoff_x=0
          t_aoff_y=0
          t_nsA=ns-n_add
          t_nlA=nl-m_add
        endif
c
c     Set the sub-matrix A where shade is to be found
c     aoff_x is the offset in samples
c     aoff_y is the offset in lines
c       pos in line in image for A(i,j) is aoff_x+j
c       pos in lines in image for A(i,j) is aoff_y+i
c
c     check the submatrix A is valid in various ways
c
c     first that the individual components are valid
c
      if (((aoff_x.lt.0) .or. (aoff_x.ge.ns)) .or.
     .    ((aoff_y.lt.0) .or. (aoff_y.ge.nl)) .or.
     .    ((nsA.lt.1) .or. (nsA.gt.ns)) .or.
     .    ((nlA.lt.1) .or. (nlA.gt.nl))) then
          write(0,*) 'Parameters defining A are invalid'
          write(0,*) 'Check parameters!'
          goto 99
        endif
c
c     Check A is embedded in the whole image
c
        if((aoff_x+nsA.gt.ns) .or. (aoff_y+nlA.gt.nl)) then
          write(0,*) 'Matrix A not embedded in image'
          write(0,*) 'Check parameters!'
          goto 99
        endif
c
c     check the sub-image A plus the border area
c     needed to test A still fits inside the main image
c     with the buffer available
c
c     NOTE: treat the four cases in pairs
c
      if(az_case.eq.1.or.az_case.eq.2) then
          if (aoff_y.lt.m_add) then
            write(0,*) 'matrix A does not have sufficient y buffer'
            write(0,*) 'check numbers'
            goto 99
          endif
        else if (az_case.eq.3.or.az_case.eq.4) then
          if (aoff_y+nlA+m_add.gt.nl) then
            write(0,*) 'matrix A does not have sufficient y buffer'
            write(0,*) 'check numbers'
            goto 99
          endif
        endif
      if(az_case.eq.2.or.az_case.eq.3) then
          if (aoff_x.lt.n_add) then
            write(0,*) 'matrix A does not have sufficient x buffer'
            write(0,*) 'check numbers'
            goto 99
          endif
        else if (az_case.eq.1.or.az_case.eq.4) then
          if (aoff_x+nsA+n_add.gt.ns) then
            write(0,*) 'matrix A does not have sufficient x buffer'
            write(0,*) 'check numbers'
            goto 99
          endif
        endif
c
c     now set up the mask image to record shade pixels in
c     A NOTE: mask has the dimensions of A and not the
c     DEM - set as a 1-D array so that it can be indexed
c     in the subroutine without worrying about set bounds
c
c     Mask is a long integer (integer*4) due to recl issue
c     First set to 1 so zero will represent deep shadow
c
      do i=1,nlA
          do j=1,nsA
            mask(i,j)=1
          enddo
        enddo

        rmax=maxval(mask(1:nlA,1:nsA))
        rmin=minval(mask(1:nlA,1:nsA))
c
c     proj_terrain does the job of checking for occlusion
c     along the vector from a pixel to the sun
c     if any terrain obstructs the path the pixel is set to
c     zero in Mask
c

      call proj_terrain(ns,nl,nsA,nlA,a,mask,
     . Aoff_x,Aoff_y,k_max,
     . n_inc,m_inc,h_offset,zmax,htol,k_setting)
c
c     Description of proj_terrain
c
c      call proj_terrain(n_max,m_max,n,m,z,mask,n_off,m_off,k_max,
c     . n_inc,m_inc,h_offset,zmax,htol,k_setting)
c
c     subroutine to construct mask of shade pixels
c
c     z(m_max,n_max) is the main array of heights
c     A(m,n) is the (sub-)matrix of target heights in:
c       z(m_off+1,n_off+1) to z(m_off+m,n_off+n)
c     mask(m,n) is the output mask
c       on input assumed to be 1 where valid data exist 0 else
c     k_max is the number of lags in the projection
c     n_inc(k_max) real increments column number for projection
c     m_inc(k_max)  real increments row number for the projection
c     h_offset(k_max) is the height of the projection hor each lag
c     zmax is the maximum altitude in the whole array
c     htol is a tolerance (m) for the test for a hit (RMS error in z)
c
99    continue
      return
c
      end
c
c---------------------------------------------------------------------------
      logical function set_border(phi_sun,zmax,zmin,sun_zen,hx,hy,
     .        az_case,d,d0,k_max,h_offset,n_inc,m_inc,n_add,m_add,
     .        k_setting,add_max,err)
c
c     set_border defines the buffer that is needed to find the
c     occluding terrain. This will be larger for low sun elevations
c
      integer k_max,k,n_add,m_add,k_setting,add_max,err,az_case
      real n_inc(k_setting),m_inc(k_setting)
        real h_offset(k_setting)
        real phi_sun,zmax,zmin,sun_zen,hx,hy,phc
        real sinphc,cosphc,d,d0
c
      real pi,r2d,d2r,eps
      common/base/pi,r2d,d2r,eps
c
      set_border=.true.
      err=0
c
c     the case is dependent on the sun azimuth
c     this defines which borders need to be available
c
c     NOTE: azimuth must be in degrees between 0 and 360 and clockwise from N
c
      if (phi_sun.ge.0.0 .and. phi_sun.le.90.0) then
          az_case=1
          phc=pi/2.0-phi_sun*d2r
        else if (phi_sun.ge.270.0 .and. phi_sun.le.360.0) then
          az_case=2
          phc=phi_sun*d2r-3.0*pi/2.0
        else if (phi_sun.ge.180.0 .and. phi_sun.le.270.0) then
          az_case=3
          phc=3.0*pi/2.0-phi_sun*d2r
        else if (phi_sun.ge.90.0 .and. phi_sun.le.180.0) then
          az_case=4
          phc=phi_sun*d2r-pi/2.0
        else
          write(0,*) 'error - azimuth case not possible'
          write(0,*) 'phi_sun=',phi_sun,' must be in 0 to 360 deg'
          err=3
          goto 99
        endif
c
c     calculate the border info
c     and the increments for the projection line
c
      sinphc=sin(phc)
        cosphc=cos(phc)
c
c     d is the plane distance from the pixel
c     with zmin where if the pixel at the
c     distance in the sun direction is at Zmax
c     then it can just occlude the first pixel
c
        d=(zmax-zmin)/tan(pi/2.0-sun_zen*d2r)
        if(cosphc*hy.gt.sinphc*hx) then
          d0=0.5*hx/cosphc
        else
          d0=0.5*hy/sinphc
        endif
c
c     d0 is the basic step in the d direction that increments
c     0.5 pixel in the larger of the two directions
c
c     k_max is the number of steps to get out to d
c     n_inc(k) and m_inc(k) are the sample and line increments
c     h_offset(k) on exit is the altitude above the start of the
c       sun beam as you move towards the sun
c
        k_max=ifix(d/d0+0.5)+1
        if (k_max.gt.k_setting) then
          write(0,*) 'k_max=',k_max
          write(0,*) 'Error... k_max gt k_setting'
          err=1
          goto 99
        endif
        do 10 k=1,k_max
          h_offset(k)=float(k)*d0
c
c     the increments to search for occlusions have sign depending
c     on the azimuth case
c
          if (az_case.eq.1.or.az_case.eq.4) then
            n_inc(k)=h_offset(k)*cosphc/hx
          else
            n_inc(k)=-h_offset(k)*cosphc/hx
          endif
          if (az_case.eq.3.or.az_case.eq.4) then
            m_inc(k)=h_offset(k)*sinphc/hy
          else
            m_inc(k)=-h_offset(k)*sinphc/hy
          endif
          h_offset(k)=h_offset(k)*tan(pi/2.0-sun_zen*d2r)
10    continue
c
c     n_add and m_add are the sizes of the pixel buffer
c     in sample and line dimensions
c
c     the two buffers will be on sides of the target image defined
c     by the azimuth case
c
      n_add=ifix(d*cosphc/hx+1.5)
        m_add=ifix(d*sinphc/hy+1.5)

        if ((n_add.gt.add_max .or. m_add.gt.add_max)
     .    .or. (n_add.lt.0.or.m_add.lt.0)) then
          write(0,*) 'n_add,m_add=',n_add,m_add
          write(0,*) 'Error... add outside add_max ranges'
          err=2
          goto 99
        endif
       write(98,*)az_case,d,d0,k_max,n_add,m_add
      return
99    set_border=.false.
      return
c
      end
c
c---------------------------------------------------------------------
      subroutine proj_terrain(n_max,m_max,n,m,z,mask,n_off,m_off,k_max,
     . n_inc,m_inc,h_offset,zmax,htol,k_setting)
c
c     subroutine to construct mask of shade pixels
c
c     z(m_max,n_max) is the main array of heights
c     A(m,n) is the (sub-)matrix of target heights in:
c       z(m_off+1,n_off+1) to z(m_off+m,n_off+n)
c     mask(m,n) is the output mask
c       on input assumed to be 1 where valid data exist 0 else
c     k_max is the number of lags in the projection
c     n_inc(k_max) real increments column number for projection
c     m_inc(k_max) real increments row number for the projection
c     h_offset(k_max) is the height of the projection hor each lag
c     zmax is the maximum altitude in the whole array
c     htol is a tolerance (m) for the test for a hit (RMS error in z)
c
      integer n_max, m_max, n, m, n_off, m_off, k_max
c     NOTE: n_inc and m_inc are Floating Point arrays
        real m_inc(k_max), n_inc(k_max)
        real h_offset(k_max)
        real z(k_setting,k_setting)
        integer*2 mask(k_setting,k_setting)
       real zmax, t, tt, test, htol, xd, yd, zpos
       integer i, j, ii, jj, k, ipos, jpos
        integer*4 itot
      real pi,r2d,d2r,eps
      common/base/pi,r2d,d2r,eps
c
c     loop over object matrix A and project out into buffer
c
c
c     The main 2 loops are over the pixels in the submatrix A
c
c     For given A(i,j) the search for occluding terrain is done along a "line"
c     in the sun direction. The search can stop when the terrain would
c     have to be higher than the maximum value to occlude the current test
c     pixel (i,j) in the sub-matrix A
c
      write(98,*) 'in proj_terrain'
        write(98,*) 'htol=',htol
        write(98,*) 'input zmax',zmax
        write(98,*) 'zmax from data=',maxval(z(1:m_max,1:n_max))
        write(98,*) 'zmin from data=',minval(z(1:m_max,1:n_max))
        write(98,*) 'mask_max from data=',maxval(mask(1:m,1:n))
        write(98,*) 'mask_min from data=',minval(mask(1:m,1:n))
      itot=0
        write(98,*) 'z(250,250)=',z(250,250)
      do 100 i=1,m
          ii=m_off+i
          do 110 j=1,n
              jj=n_off+j
              t=z(ii,jj)
              do 120 k=1,k_max
                tt=t+h_offset(k)
               if(tt .le. zmax+htol) then
                 ipos=ifix(float(ii)+m_inc(k))
                 jpos=ifix(float(jj)+n_inc(k))
                 yd=float(ii)+m_inc(k)-float(ipos)
                 xd=float(jj)+n_inc(k)-float(jpos)
                 zpos=z(ipos,jpos)*(1.0-yd)*(1.0-xd)+
     .               z(ipos,jpos+1)*xd*(1.0-yd)+
     .               z(ipos+1,jpos+1)*xd*yd+
     .               z(ipos+1,jpos)*(1.0-xd)*yd
                 test=tt-zpos

                  if (test.le.htol) then
                            mask(i,j)=0
                    itot=itot+1
                    goto 125
                  endif
                else
                  go to 125
                endif
120         continue
125       continue
110     continue
100   continue
c      write(0,*) 'itot=',itot
c
      return
        end
c------------------------------------------------
c
      subroutine pixelsize(rlat,dres,dx,dy)

c     subroutine is used to calculate pixel size (in meters) at latitude and
c     longitude projection
      double precision aa,bb,cc,dd,ff,rlat,rlon,pi
      double precision pia,rr,dres,ddx,ddy
      real dx,dy
c     set projection parameters. here WGS84 is used
c     semi-major axis
      aa=6.3781370d6
c     flattening
      ff=2.98257223563d2
c     semi-minor axis
      bb=aa*(1.-1/ff)
      pi=4.0*atan(1.0)
      pia=pi/180.0
      cc=aa*cos(rlat*pia)
      dd=bb*sin(rlat*pia)
      rr=sqrt((aa**2*cc**2+bb**2*dd**2)/(cc**2+dd**2))
      ddy=dres*pia
      ddx=acos(sin(rlat*pia)**2+cos(rlat*pia)**2*cos(dres*pia))
      dy=rr*ddy
      dx=rr*ddx
      return
      end
'''

# The following is the original Fortran code for the self shadow algorithm.
# Filename: slope_pixelsize_newpole.f
'''
      program slope
c     this program is used to calculate slope and aspect angles
c     using Sobel filter and then calculate incident and
c     exiting angles as well as their azimuth angles.
c     note: the row and column of DEM data must be larger
c     than the image (extra each line and column for the four sides.
c     it is needed for sobel filter.
c     command line  slope <header> <dem> <solar> <view>
c                   <soalr_azimuth> <view_azimuth>
c                   <slope> <aspect> <incident>
c                   <inci_azi> <exiting> <exit_azi>
c                  <rela_angle> <mask>
c
      character fname*80
      real dem(20000,20000)
      double precision alat(20000),alon(20000)
      real theta(20000),phit(20000),dem_ori(20000)
      real solar(20000),sazi(20000)
      real view(20000),azi(20000)
      real it(20000),et(20000)
      real azi_it(20000),azi_et(20000),rela(20000)
	real offset
      double precision dx,dy
      double precision UTMN1,UTMN2,UTMN3
      double precision UTME1,UTME2,UTME3
      double precision pi,pia,pib,dres,rlat1,rlon1
      double precision p,q
	integer iargc,i,j
      integer ierr,Aoff_x1,Aoff_x2
      integer Aoff_y1,Aoff_y2,nrow,ncol
      integer*2 mask(20000)
      if (IARGC() .ne. 14) then
      write(*,*) 'Error: Required parameters not specified properly!'
      stop 11
      endif
c--------------------------------------------------------------------
c      open header file
      call GETARG(1, fname)
      open(1,file=fname,status='old')
c      read header file
c     read row and column (Landsat image)
      read(1,*)nrow,ncol
c     read extra DEM pixel number for each of x directions
      read(1,*)Aoff_x1,Aoff_x2
c     read extra line number for each of y directions
      read(1,*)Aoff_y1,Aoff_y2
c     read spatial resolution
      read(1,*)dres
c     read upperleft lat and lon
      read(1,*)rlat1,rlon1
c---------------------------------------------------------------
c      open input file
c      open DEM file. DEM file is the same as teh one used for
c      cast shadow.
      call GETARG(2, fname)
      open (2,file=fname,status='old',access='direct',
     ! form='unformatted',recl=4*(ncol+Aoff_x1+Aoff_x2))
c------------------------------------------------------------
c      open other input files (solar, view)
       do i=3,4
       call GETARG(i, fname)
      open (i,file=fname,status='old',access='direct',
     ! form='unformatted',recl=4*ncol)
       enddo
c-------------------------------------------------------------
c      open other input files (solar_azi,view_azi)
       do i=5,6
        call GETARG(i, fname)
      open (i+2,file=fname,status='old',access='direct',
     ! form='unformatted',recl=4*ncol)
       enddo
c---------------------------------------------------------------
c      open output files (slope,apsect,incident,inci_azi,
c      exiting,exit_azi,rela_angle
       do i=7,13
        call GETARG(i, fname)
       open (i+6,file=fname,access='direct',
     ! form='unformatted',recl=4*ncol)
       enddo
c      open shadow mask
        call GETARG(14, fname)
       open (20,file=fname,access='direct',
     ! form='unformatted',recl=2*ncol)

c---------------------------------------------------------------------
c      upper left coordinator for satellite image.
c      here we assume that Landsat and DEM have some spatial resolution
c      and coordinator. The information can be obtained from satellite
c      header file. DEM header should add an extra pixel.
       ierr=0
       pi=4.0d0*atan(1.0d0)
       pia=pi/180.0d0
       pib=180.0d0/pi
c----------------------------------------------------------------
c      calculate latitude and longitude for each pixel (DEM data)
c      for first pixel of DEM data
       alat(1)=rlat1+dres
       alon(1)=rlon1-dres
       do i=2,nrow+2
        alat(i)=rlat1-dble(i-2)*dres
        enddo
       do j=2,ncol+2
       alon(j)=rlon1+dble(j-2)*dres
       enddo
c----------------------------------------------------------------------
c      read DEM data

25       do i=1,nrow+2
       read(2,rec=Aoff_y1-1+i)(dem_ori(j),j=1,ncol+Aoff_x1+Aoff_x2)
       do j=1,ncol+2
       dem(i,j)=dem_ori(Aoff_x1-1+j)
        enddo
       if (i .eq. 1) then
       print*,dem(1,1)
       endif
       enddo
c-----------------------------------------------------------------------
c       start to calculate angle
       do i=2, nrow+1
       read(3,rec=i-1)(solar(j),j=2,ncol+1)
       read(4,rec=i-1)(view(j),j=2,ncol+1)
       read(7,rec=i-1)(sazi(j),j=2,ncol+1)
       read(8,rec=i-1)(azi(j),j=2,ncol+1)
c      calculate pixel size in meters
       call pixelsize(alat(i),dres,dx,dy)
c-------------------------------------------------------------------
       do j=2,ncol+1
       mask(j)=1
       if (sazi(j) .le. -180.0) sazi(j)=sazi(j)+360.0
       if (sazi(j) .gt. 180.0) sazi(j)=sazi(j)-360.0
       if (azi(j) .le. -180.0) azi(j)=azi(j)+360.0
       if (azi(j) .gt. 180.0) azi(j)=azi(j)-360.0
       p=(dble(dem(i-1,j+1))-dble(dem(i-1,j-1))+
     .    2.0d0*(dble(dem(i,j+1))-dble(dem(i,j-1)))+
     .    dble(dem(i+1,j+1))-dble(dem(i+1,j-1)))/(8.0d0*dx)
       q=(dble(dem(i-1,j-1))-dble(dem(i+1,j-1))+
     .    2.0d0*(dble(dem(i-1,j))-dble(dem(i+1,j)))+
     .    dble(dem(i-1,j+1))-dble(dem(i+1,j+1)))/(8.0d0*dy)
c----------------------------------------------------------------
       theta(j)=sngl(atan(sqrt(p**2+q**2))*pib)
       phit(j)=sngl(atan2(-p,-q)*pib)
       if (phit(j) .le. -180.0) phit(j)=phit(j)+360.0
       if (phit(j) .gt. 180.0) phit(j)=phit(j)-360.0
       call cal_pole(solar(j),sazi(j),theta(j),phit(j),
     # it(j),azi_it(j),ierr,offset)
       call cal_pole(view(j),azi(j),theta(j),phit(j),
     # et(j),azi_et(j),ierr,offset)
       rela(j)=azi_it(j)-azi_et(j)
       if (rela(j) .le. -180.0) rela(j)=rela(j)+360.0
       if (rela(j) .gt. 180.0) rela(j)=rela(j)-360.0
       if (cos(it(j)*pia) .le. 0.0) then
        mask(j)=0
c        it(j)=90
        endif
        if (cos(et(j)*pia) .le.0.0 ) then
        mask(j)=0
c        et(j)=90
        endif
c----------------------------------------------------------------
       if (i.eq.16 .and. j.eq.5269) then
       print*,alat(i),alon(j)
       print*,dx,dy,solar(j),view(j),sazi(j),azi(j)
       print*,dem(i-1,j-1),dem(i-1,j),dem(i-1,j+1)
       print*,dem(i,j-1),dem(i,j),dem(i,j+1)
       print*,dem(i+1,j-1),dem(i+1,j),dem(i+1,j+1)
       print*,p,q,sqrt(p*p+q*q)
       print*,theta(j),phit(j),it(j),et(j),
     # azi_it(j),azi_et(j)
       endif
c----------------------------------------------------------------

c       print*,p,q,atan(tantheta(j))*180/pi,phit(j)
        enddo
        write(13,rec=i-1)(theta(j),j=2,ncol+1)
        write(14,rec=i-1)(phit(j),j=2,ncol+1)
        write(15,rec=i-1)(it(j),j=2,ncol+1)
        write(16,rec=i-1)(azi_it(j),j=2,ncol+1)
        write(17,rec=i-1)(et(j),j=2,ncol+1)
        write(18,rec=i-1)(azi_et(j),j=2,ncol+1)
        write(19,rec=i-1)(rela(j),j=2,ncol+1)
        write(20,rec=i-1)(mask(j),j=2,ncol+1)
        enddo
       stop
       end

      subroutine pixelsize(rlat,dres,dx,dy)

c     subroutine is used to calculate pixel size (in meters) at latitude and
c     longitude projection
      double precision aa,bb,cc,dd,ff,rlat,rlon,pi
      double precision pia,rr,dres,dx,dy,ddx,ddy
c     set projection parameters. here WGS84 is used
c     semi-major axis
      aa=6.3781370d6
c     flattening
      ff=2.98257223563d2
c     semi-minor axis
      bb=aa*(1.0d0-1.0d0/ff)
      pi=4.0*atan(1.0d0)
      pia=pi/180.0
      cc=aa*cos(rlat*pia)
      dd=bb*sin(rlat*pia)
      rr=sqrt((aa**2*cc**2+bb**2*dd**2)/(cc**2+dd**2))
      ddy=dres*pia
      ddx=acos(sin(rlat*pia)**2+cos(rlat*pia)**2*cos(dres*pia))
      dy=rr*ddy
      dx=rr*ddx
      return
      end
c
      subroutine cal_pole(theta,phi,theta_p,phi_p,thp,php,ierr,offset)
c
        real theta,phi,theta_p,phi_p,thp,php,offset
        real pdiff,costhp,sinphp,cosphp,tnum,tden
        real eps,pi,d2r
        integer ierr
c
        ierr=0
        offset=0.0
c
        eps=1.0e-6
        pi=4.0*atan(1.0)
        d2r=pi/180.0
c
        if (abs(theta_p).le.eps) then
          thp=theta
	  php=phi
          return
	endif
c
        offset=atan(tan(pi-d2r*phi_p)*cos(d2r*theta_p))
        pdiff=d2r*(phi-phi_p)
c
        costhp=cos(d2r*theta)*cos(d2r*theta_p)+sin(d2r*theta)*
     .  sin(d2r*theta_p)*cos(pdiff)
        if (costhp.ge.1.0-eps) then
          thp=0.0
          php=0.0-offset/d2r
          return
        else
          thp=acos(costhp)/d2r
          sinphp=sin(d2r*theta)*sin(pdiff)/sin(d2r*thp)
          cosphp=(cos(d2r*theta)*sin(d2r*theta_p)-sin(d2r*theta)*
     .  cos(d2r*theta_p)*cos(pdiff))/sin(d2r*thp)
c
          if (abs(sinphp).le.eps) then
            if (cosphp.gt.eps) then
              php=0.0-offset/d2r
	    else
              php=180.0-offset/d2r
	    endif
	    return
	  else if (abs(cosphp).le.eps) then
	    if (sinphp.gt.eps) then
	      php=90.0-offset/d2r
	    else
             php=-90.0-offset/d2r
            endif
            return
          endif
        endif
c
        tnum=sin(d2r*theta)*sin(pdiff)
        tden=(cos(d2r*theta)*sin(d2r*theta_p)-sin(d2r*theta)*
     .	cos(d2r*theta_p)*cos(pdiff))
        php=(atan2(tnum,tden)-offset)/d2r
        if (php.gt.180.0) php=php-360.0
        if (php.lt.-180.0) php=php+360.0
c
        return
        end
'''
