import numexpr, fnmatch, os, errno, numpy, gc, re, logging
from osgeo import gdal, gdalconst
from ULA3.dataset import SceneDataset
from ULA3.meta import print_call
from ULA3.fc.utils import datatype, create_dir, unmix

from endmembers import EndMember, EndMemberFactory

"""
Fractional cover calculations. The only function that is likely to be of interest to users
is :py:func:`fractional_cover`. The other modules under fc contain utility functions used internally.
Some of these are likely to be moved to other modules, renamed or removed completely in the future and
hence should not be used directly.
"""

logger = logging.getLogger('root.' + __name__)

class FractionalCoverResult(object):
    """
    Holds the results from a fracional cover analysis.
    """
    def __init__(self, green, dead, bare, unmix_err):
        self.green = green #! Josh - please document this.
        self.dead = dead #! Josh - please document this.
        self.bare =bare #! Josh - please document this.
        self.unmix_err = unmix_err #! Josh - please document this.





@print_call(logger.info)
def fractional_cover(nbar_data_path, asfloat32=False, fc_data_path=None, single_tif=False):
    '''
    The process flow function for running fractional cover.

    :param nbar_data_path: The input directory containing the reflectance images. This must be a string
        that returns points to a :py:class:`ULA3.dataset.SceneDataset`.

    :param asfloat32: Should the results be returned (and possibly written to disk) as :py:class:`numpy.float32`
        (``True``) or :py:class:`numpy.int16` (``False``). Note that the results are scaled differently depending on
        the return type. See details below.

    :param fc_data_path: The output directory for fractional cover. If this is ``None`` the the results
        are not written to disk.

    :return:
        An instance of :py:class:`FractionalCoverResult`.

    :return type:
        * Int16: Fractions are expressed as a percent with a scale factor of 1000.
        * Float32: Fractions are expressed as a percent with a scale factor of 0.1.

    :warning:
        Developers - some IDEs will highlight unused variables... be careful here, as the variables
        in this function that appear to be unused are referenced in expressions used by :py:mod:`numexpr`... so
        please don't delete them!
    '''

    # PRODUCTION TWEAK
    # Create output filename that complies with GA dataset spec.
    # TODO We'll need to modify to support writing one FC file per band.

    # TODO We shouldn't really be trusting the file name - best to get from metadata.
    # N.B: Product TIF name is taken to the the same as the FC dataset directory
    _name = os.path.basename(nbar_data_path)
    _product_id = os.path.basename(fc_data_path)
    _output_fname = _product_id + '.tif'
    res_outdir = os.path.join(fc_data_path, 'scene01')
    outname = os.path.join(res_outdir, _output_fname)

    nbar_dataset = SceneDataset(nbar_data_path)

#    b_match = [2,3,4,5,7] # Arbitrary list for LS5 & LS7 NBAR - this will break for LS8

    # b_match is a list of all reflective bands excluding visible blue
    b_match = [band_no for band_no in nbar_dataset.bands('REFLECTIVE')
               if nbar_dataset.sensor_band_info(band_no)['WAVELENGTH'][1] > 0.52]

    iprj = nbar_dataset.GetProjection()
    igeot = nbar_dataset.GetGeoTransform()
    band  = nbar_dataset.GetRasterBand(1)
    datype = datatype(band.DataType)
    r_data = numpy.zeros((len(b_match), nbar_dataset.RasterYSize, nbar_dataset.RasterXSize), dtype=datype)

    for i in range(len(b_match)):
        r_data[i,:,:] = nbar_dataset.band_read_as_array(b_match[i])

    # Get the dimensions, only interested in the x and y dims
    dims = r_data.shape
    if len(dims) >2:
        ncols = dims[2]
        nrows = dims[1]
        dims  = (nrows,ncols)

    # Find the minimum value. This is tailored to GA's NBAR. Won't work if data is
    # not -999, or if a valid value is < -999
    #min_ = r_data.min() #GA NBAR data has a fill of -999
    data_ignore = band.GetNoDataValue()
    if not(data_ignore):
        # set a default of -999
        data_ignore = -999

    #null_val = numpy.int16(32766) # This was used for MODIS data
    null_val = numpy.int16(0)

    # Setting the no_data values to zero. It was just easier to do this when
    # feeding into the pixel unmixing algorithm.
    image = numexpr.evaluate("(r_data == data_ignore) * null_val + (r_data != data_ignore) * r_data")
    del r_data; gc.collect()

    # 2013_01_08_version produces green, dead1, dead2, bare, unmixing error
    frac = unmix(image)
    wh = numexpr.evaluate("image == null_val")

    # Need to change the null data values back to the original -999 null value
    wh_any = numpy.any(wh, axis=0)
    del wh; gc.collect()

    # scale factors
    sf2 = numpy.float32(0.01)
    sf3 = 10000

    # new model
    green = frac[0,:,:]
    dead1 = frac[1,:,:]
    dead2 = frac[2,:,:]
    bare  = frac[3,:,:]
    err   = frac[4,:,:] # unmixing error

    # Creating the output image
    driver = gdal.GetDriverByName("GTiff")

    create_dir(fc_data_path)
    create_dir(res_outdir)

    # **** Ideally would like to set band names ****
    #bnames = {'Band_1' : 'Green', 'Band_2' : 'Non Green', 'Band_3' : 'Bare', 'Band_4' : 'Unmixing Error'}
    #bnames = {'Band_1' : 'Green', 'Band_2' : 'Dead', 'Band_3' : 'Bare', 'Band_4' : 'Unmixing Error'}

    outds = None
    if asfloat32:
        green[wh_any] = data_ignore
        dead = numexpr.evaluate("(dead1 + dead2)").astype('float32')
        dead[wh_any] = data_ignore
        bare[wh_any] = data_ignore
        unmix_err = numexpr.evaluate("err * sf2").astype('float32')
        unmix_err[wh_any] = data_ignore

        # output as float32
        gdal_datatype = gdalconst.GDT_Float32
    else:
        green = numexpr.evaluate("green * sf3").astype('int16')
        green[wh_any] = data_ignore
        dead = numexpr.evaluate("(dead1 + dead2) * sf3").astype('int16')
        dead[wh_any] = data_ignore
        bare = numexpr.evaluate("bare * sf3").astype('int16')
        bare[wh_any] = data_ignore
        unmix_err = numexpr.evaluate("err * sf2 * sf3").astype('int16')
        unmix_err[wh_any] = data_ignore

        # output as int16
        gdal_datatype = gdalconst.GDT_Int16

    if fc_data_path:
        _fc_bands = [('PV', green), ('NPV', dead), ('BS', bare), ('UE', unmix_err)]

        if single_tif:
            outds = driver.Create(outname,  dims[1], dims[0], len(_fc_bands), gdal_datatype)
            outds.SetGeoTransform(igeot)
            outds.SetProjection(iprj)
            #outds.SetMetadata(bnames) # unfortunately doesn't actually set band names

            for band_index in range(len(_fc_bands)):
                outband = outds.GetRasterBand(band_index + 1)
                outband.SetNoDataValue(data_ignore)
                outband.WriteArray(_fc_bands[band_index][1])

        else: # Multiple single-band TIF files required
            for band_index in range(len(_fc_bands)):
                outds = driver.Create(re.sub('\.tif$', '_' + _fc_bands[band_index][0] + '.tif', outname),
                                      dims[1], dims[0], 1, gdal_datatype)
                outds.SetGeoTransform(igeot)
                outds.SetProjection(iprj)

                outband = outds.GetRasterBand(1)
                outband.SetNoDataValue(data_ignore)
                outband.WriteArray(_fc_bands[band_index][1])

        outds = None


    return FractionalCoverResult(green, dead, bare, unmix_err)

