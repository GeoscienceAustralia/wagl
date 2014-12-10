"""
Utilties used in :ref:`tc-algorithm-label`. Most of the functions are flexible as to what they take as input
datasets - they can be strings that are valid :py:class:`gdal.Dataset`s, open :py:class:`gdal.Dataset`s,
or :py:class:`numpy.ndarray`s.
"""

import os, logging
import numpy
import gdal, gdalconst
from ULA3.meta import print_call
from ULA3.utils import warp, get_bounds, DTYPE_MAP, as_array, dump_array
from gaip import filter
from gaip import shade_main_landsat_pixel
from gaip import slope_pixelsize_newpole
from gaip import terrain_correction
# from _brdf_terrain_newdiff_all_LS8 import terrain_correction_ls8

logger = logging.getLogger('root.' + __name__)





class FortranError(Exception):
    """
    Base class for errors thrown from the Fortran code used in this module.
    """
    def __init__(self, function_name, code, msg):
        self.function_name = function_name #! The name of the Fortran function called.
        self.code = code #! The error code produced by the Fortran code.
        self.msg = msg or "Unknown error" #! The Message corresponding to ``code``.

    def __str__(self):
        """
        Return a string representation of this Error.
        """
        return "Error in Fotran code %s (code %i): %s" % (self.function_name, self.code, self.msg)





@print_call(logger.info)
def clip_dsm(shape_dataset, dsm_filename, output_filename, buffer_widths, output_format):
    """
    Clip a region out of a DSM. This function calls :py:func:`ULA3.utils.warp` (which imports it from
        :py:func:`ULA3._gdal_tools.warp`). See the latter method for more details.

    :param shape_dataset:
        Dataset used to one which to base the width, height, cell size and projection of the clipped region
        on (note that buffers, as specified by ``buffer_widths``, are added to the width and height).

    :param dsm_filename:
        Dataset to clip from. This must be a string suitable for a call to gdal.Open.

    :param output_filename:
        The name of the file to write the clipped data to.

    :param buffer_widths:
        An object with members left, right, top and bottom that describe the buffer (in pixels)
        to add to the corresponding edges of the clipped region.

    :param output_format:
        A GDAL format string specifying the format of the output.

    :return:
        The clipped data (as a :py:class:`numpy.ndarray`).
    """
    output_dir = os.path.dirname(output_filename)
    assert output_dir, "output_filname must be a fully qualified file name."
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    warped_filename = warp(shape_dataset, dsm_filename, output_filename, buffer_widths, output_format)
    warped_dataset = gdal.Open(warped_filename, gdalconst.GA_ReadOnly)
    warped_data = as_array(warped_dataset.ReadAsArray(), dtype=numpy.float32)
    warped_dataset = None # close the dataset.

    return warped_data





@print_call(logger.info)
def filter_dsm(clipped_dem):
    """
    Apply Fuqin's (3x3 Gausian) smoothing filter to a clipped DEM. This code is an interface to the fortran code
    :file:`filter.f90` (which is compiled to a Python module using F2py).

    :param clipped_dem:
        Data to be smoothed (converted to an array via a call to :py:meth:`ULA3.utils.as_array`).

    :return:
        The smoothed data (as a :py:class:`numpy.ndarray` of data type :py:class:`numpy.float32`).
    """
    return filter(as_array(clipped_dem, dtype=numpy.float32))





class SlopeResultSet(object):
    """
    Holds the results of a call to :py:func:`run_slope`.
    """
    def __init__(self, mask_self, slope, aspect, incident, exiting, azi_incident, azi_exiting, rela_slope):
        """
        All arguments are :py:class:`numpy.ndarray`s. These correspond to the arguments I have no idea what these
        actually represent, so someone who does should document this (Fuqin will know).
        """
        self.mask_self = mask_self
        self.slope = slope
        self.aspect = aspect
        self.incident = incident
        self.exiting = exiting
        self.azi_incident = azi_incident
        self.azi_exiting = azi_exiting
        self.rela_slope = rela_slope

    def dump_arrays(self, output_path, l1t_input_dataset, file_type='ENVI', file_extension='.img'):
        def dumper(data, name):
            dump_array(
                array=data,
                output_path=os.path.join(output_path, name + file_extension),
                template_dataset=l1t_input_dataset,
                no_data_value=-999,
                convert_to_byte=False,
                file_format=file_type)

        dumper(self.mask_self, 'mask_self')
        dumper(self.slope, 'slope')
        dumper(self.aspect, 'aspect')
        dumper(self.incident, 'incident')
        dumper(self.exiting, 'exiting')
        dumper(self.azi_incident, 'azi_incident')
        dumper(self.azi_exiting, 'azi_exiting')
        dumper(self.rela_slope, 'rela_slope')




class SlopeError(FortranError):
    """
    Class that deals with errors from :py:func:`run_slope`.
    """
    def __init__(self, code):
        super(SlopeError, self).__init__("slope_pixelsize_newpole", code, SlopeError.get_error_message(code))

    @staticmethod
    def get_error_message(code):
        """
        Gennerate an error message for a specific code (It is OK that this have non-returning control paths).
        """
        if code in 10:
            return "X dimensions of scene and DEM not correct."
        if code == 11:
            return "Y dimensions of scene and DEM not correct."

#TODO: ensure that we are reading from the bottom, not the top (see second arg).
@print_call(logger.info)
def run_slope(
    shape_dataset,
    dem_data,
    solar_zenith_data,
    view_zenith_data,
    solar_azimuth_data,
    view_azimuth_data,
    pix_buf,
    is_utm,
    spheroid,
    output_type = "ENVI",
    slope_dataset = None,
    aspect_dataset = None,
    incident_dataset = None,
    azi_incident_dataset = None,
    exiting_dataset = None,
    azi_exiting_dataset = None,
    rela_slope_dataset = None,
    mask_self_dataset = None):
    """
    Calculate the slope and angles for a region. This code is an interface to the fortran code slope_pixel_newpole.f90
    written by Fuqin (which was modified to work with F2py).

    The following was taken from the top of the Fotran program: "slope_pixelsize_newpole.f90:
    This program is used to calculate slope and aspect angles
    using Sobel filter and then calculate incident and
    exiting angles as well as their azimuth angles.
    note: the row and column of DEM data must be larger
    than the image (extra each line and column for the four sides.
    it is needed for sobel filter.

    :param shape_dataset:
        Object used to specify the 'shape' of the region (origin, cell size and dimensions).

    :param dem_data:
        The DEM data for the region. This must have the same dimensions as ``shape_dataset`` plus buffers as specified
        by ``pix_buf`` (see :py:func:`clip_dsm`). This is converted to a :py:class:`gdal.ndarray` of type
        :py:class:`gdal.float32` using the function :py:func:`ULA3.utils.as_array`.

    :param solar_zenith_data:
        The solar zenith angle data for the region.

    :param view_zenith_data:
        The satellite zenith angle data for the region.

    :param solar_azimuth_data:
        The solar azimuth angle data for the region.

    :param view_azimuth_data:
        The satellite azimuth angle data for the region.

    :param pix_buf
        An object with members top, bottom, left and right giving the size of the buffer (in pixels) which
        have been added to the corresponding sides of ``dem_data``.

    :param is_utm:
        Boolean specifying whether the data is in UTM coordinates. This is passed to the Fortran code which
        uses it to determine whether the edge length of the pixels (if it is not true, then the resolution
        is assumed to be in degrees, and the pixel size is calculated internally).

    :param spheroid:
        A 4 element floating point array containing the Earth
        spheroidal paramaters.
        Index 0 contains the spheroid Major Axis.
        Index 1 contains the spheroid Inverse Flattening.
        Index 2 contains the spheroid Squared Eccentricity.
        Index 3 contains the Earth rotational angular velocity in
        radians/second.

    :param output_type:
        The output types for any datasets written to disk (see the following arguments).

    :param slope_dataset:
        Defines the output type of the slope results (see details below).

    :param aspect_dataset:
        Defines the output type of the aspect angle data to (see details below).

    :param incident_dataset:
        Defines the output type of the incident angle data (see details below).

    :param exiting_dataset:
        Defines the output of the exiting (zenith?) angle data (see details below).

    :param azi_exiting_dataset:
        The name of the file to write the exiting azimuth angle data to (see details below).

    :param rela_slope_dataset:
        The name of the file to write the relative slope data to (see details below).

    :param _dataset: The name of the file to write the 'self mask' data to (see details below).

    The parameters ``solar_zenith_data, view_zenith_data, solar_azimuth_data`` and ``view_azimuth_data``
    must have the same dimensions as ``shape_dataset``. They argument are converted to
    :py:class:`gdal.ndarray`s of type :py:class:`gdal.float32` using the function :py:func:`ULA3.utils.as_array`.

    :todo: This documentation would benefit from having someone who understands the
        concepts describe the outputs more meaningfully.
    """
    bounds = get_bounds(shape_dataset)

    nrow = bounds.RasterYSize + 2
    ncol = bounds.RasterXSize + 2
    dres = bounds.RasterCellSize # assumes that the x and y res are the same.
    rlat = bounds.RasterYOrigin

    # x & y pixel resolution (This should handle cases of non-square pixels.)
    dresx = abs(bounds.RasterXCellSize)
    dresy = abs(bounds.RasterYCellSize)

    dem_dat = as_array(dem_data, dtype=numpy.float32)[(pix_buf.top-1):-(pix_buf.bottom-1),(pix_buf.left-1):-(pix_buf.right-1)]
    assert dem_dat.shape == (nrow, ncol), "dem_data not of correct shape " + str((nrow, ncol)) + " != " + str(dem_dat.shape)

    # This will be ignored if is_utm == True
    alat = numpy.array([rlat-i*dresy for i in range(-1, nrow-1)], dtype=numpy.float64) # yes, I did mean float64.

    mask, theta, phit, it, et, azi_it, azi_et, rela, ierr = slope_pixelsize_newpole(
        dresx, dresy, spheroid, alat, is_utm,
        dem_dat,
        as_array(solar_zenith_data, dtype=numpy.float32),
        as_array(view_zenith_data, dtype=numpy.float32),
        as_array(solar_azimuth_data, dtype=numpy.float32),
        as_array(view_azimuth_data, dtype=numpy.float32))

    if ierr:
        raise SlopeError(ierr)

    driver = gdal.GetDriverByName(output_type)

    def make_dataset(dataset, data):
        if type(dataset) == str:
            dataset = driver.Create(dataset, ncol, nrow, 1, DTYPE_MAP[numpy.float32])
            dataset.GetRasterBand(1).WriteArray(data)
            dataset = None
        return data

    return SlopeResultSet(
        mask_self = make_dataset(mask_self_dataset, mask),
        slope = make_dataset(slope_dataset, theta),
        aspect = make_dataset(aspect_dataset, phit),
        incident = make_dataset(incident_dataset, it),
        exiting = make_dataset(exiting_dataset, et),
        azi_incident = make_dataset(azi_incident_dataset, azi_it),
        azi_exiting = make_dataset(azi_exiting_dataset, azi_et),
        rela_slope = make_dataset(rela_slope_dataset, rela))





class CastShadowError(FortranError):
    """
    Class that deals with errors from :py:func:`run_castshadow`.
    """
    def __init__(self, code):
        super(CastShadowError, self).__init__("shade_main_landsat_pixel", code, CastShadowError.get_error_message(code))

    @staticmethod
    def get_error_message(code):
        """
        Generate an error message for a specific code. It is OK for this have non-returning control paths,
            as this will results in ``None``, which is handled in the super class.
        """
        def tmpt(d, n):
            return "attempt to access invalid %s of %s" % (d, n)

        if code == 20:
            return tmpt('x', 'dem')
        if code == 21:
            return tmpt('x', 'dem_data')
        if code == 22:
            return tmpt('x', 'solar and sazi')
        if code == 23:
            return tmpt('x', 'solar_data')
        if code == 24:
            return tmpt('x', 'a')
        if code == 25:
            return tmpt('y', 'dem_data')
        if code == 26:
            return tmpt('y', 'a')
        if code == 27:
            return tmpt('x', 'mask_all')
        if code == 28:
            return tmpt('y', 'mask_all')
        if code == 29:
            return tmpt('x', 'mask')
        if code == 30:
            return tmpt('y', 'mask')
        if code == 31:
            return tmpt('X', 'dem and a')
        if code == 32:
            return tmpt('y', 'a')
        if code == 33:
            return tmpt('y', 'dem')
        if code == 34:
            return tmpt('x', 'mask_all')
        if code == 35:
            return tmpt('x', 'mask')
        if code == 36:
            return tmpt('y', 'mask_all')
        if code == 37:
            return tmpt('y', 'mask')
        if code == 38:
            return tmpt('x', 'dem')
        if code == 39:
            return tmpt('x', 'dem_data')
        if code == 40:
            return tmpt('x', 'solar')
        if code == 41:
            return tmpt('x', 'solar_data')
        if code == 42:
            return tmpt('x', 'a and dem')
        if code == 43:
            return tmpt('y', 'a')
        if code == 44:
            return tmpt('y', 'dem')
        if code == 45:
            return tmpt('x', 'mask_all')
        if code == 46:
            return tmpt('x', 'mask')
        if code == 47:
            return tmpt('y', 'mask_alll')
        if code == 48:
            return tmpt('y', 'mask')
        if code == 49:
            return tmpt('x', 'a and dem')
        if code == 50:
            return tmpt('y', 'a')
        if code == 51:
            return tmpt('y', 'dem')
        if code == 52:
            return tmpt('x', 'mask_all')
        if code == 53:
            return tmpt('x', 'mask')
        if code == 54:
            return tmpt('y', 'mask_all')
        if code == 55:
            return tmpt('y', 'mask')
        if code == 61:
            return "azimuth case not possible - phi_sun must be in 0 to 360 deg"
        if code == 62:
            return "k_max gt k_setting"
        if code == 63:
            return "add outside add_max ranges"
        if code == 71:
            return "Parameters defining A are invalid"
        if code == 72:
            return "Matrix A not embedded in image"
        if code == 73:
            return "matrix A does not have sufficient y buffer"
        if code == 74:
            return "matrix A does not have sufficient x buffer"





@print_call(logger.info)
def run_castshadow(
    shape_dataset,
    dem_data,
    solar_angle_data,
    sazi_angle_data,
    pix_buf,
    block_height,
    block_width,
    is_utm,
    spheroid):
    """
    This code is an interface to the fortran code shade_main_landsat_pixel.f90 written by Fuqin
    (and modified to work with F2py).

    The following was taken from the top of the Fotran program: "shade_main_landsat_pixel.f90":

    Creates a shadow mask for a standard Landsat scene
    the program was originally written by DLB Jupp in Oct. 2010
    for a small sub_matrix and was modified by Fuqin Li in Oct.
    2010 so that the program can be used for large landsat scene.

    Basically, a sub-matrix A is embedded in a larger DEM image
    and the borders must be large enough to find the shaded pixels.
    If we assume the solar azimuth and zenith angles change very
    little within the sub-matrix A, then the Landsat scene can be
    divided into several sub_matrix.
    For Australian region, with 0 .00025 degree resolution, the
    sub-marix A is set to 500x500

    we also need to set extra DEM lines/columns to run the Landsat
    scene (see parameter ``pix_buf``. This will change with elevation difference within the
    scene and solar zenith angle. For Australian region and Landsat
    scene with 0.00025 degree resolution, the maximum extra lines
    are set to 250 pixels/lines for each direction. This figure
    shold be sufficient for everywhere and anytime in Australia.
    thus the DEM image will be larger than landsat image for
    500 lines x 500 columns

    :param shape_dataset:
        Dataset that defines the shape of the region.
    :type shape_dataset:
        anything that can be passed to :py:class:`ULA3._gdal_tools.default_bounds_getter`.

    :param dem_data:
        A DEM of the region. This must have the same dimensions as ``shape_dataset`` plus a buffer of
        widths specified by ``pix_buf``.
    :type dem_data:
        Anything that can be passed to :py:func:`ULA3.utils.as_array` along with a dtype of :py:class:`numpy.float32`.

    :param solar_angle_data:
        Array of solar zenith angles (in degrees).
    :type solar_angle_data:
        Anything that can be passed to :py:func:`ULA3.utils.as_array` along with a dtype of :py:class:`numpy.float32`.

    :param sazi_angle_data:
        Array of solar azimuth angles (in degrees).
    :type sazi_angle_data:
        Anything that can be passed to :py:func:`ULA3.utils.as_array` along with a dtype of :py:class:`numpy.float32`.

    :param pix_buf:
        Object describing the buffers around ``dem_data``.
    :type pix_buf:
        :py:class:`ULA3._gdal_tools.ImageShape`.

    :param block_height:
        The height of the sub-array to be embedded (see notes above).
    :type block_height:
        int

    :param block_width:
        The width of the sub-array to be embedded (see notes above).
    :type block_width:
        int

    :param spheroid:
        A 4 element floating point array containing the Earth
        spheroidal paramaters.
        Index 0 contains the spheroid Major Axis.
        Index 1 contains the spheroid Inverse Flattening.
        Index 2 contains the spheroid Squared Eccentricity.
        Index 3 contains the Earth rotational angular velocity in
        radians/second.

    :warning:
        The parameters ``solar_angle_data`` and ``sazi_angle_data`` require inputs that are in degrees.
        This is different to most other functions in ULA3. This is the case because this is just a thin wrapper
        around Fuqin's Fortran code, which expects arrays in degrees.

    :warning:
        The Fortran code cannot be compiled with ``-O3`` as it produces incorrect results if it is.

    :todo:
        Perhaps the functions ``solar_angle_data`` and ``sazi_angle_data`` should accept that arrays are
        in radians rather than degrees for consistency with the rest of ULA3.
    """
    # save the inputs (this has been useful for debugging - they get read again from
    # ULA3.tests.RunCastShadowTestCase.test3).
    #numpy.save(file="/short/v10/tmp/ula3_tests/nbar/work/cs_dem_data.npy", arr=dem_data)
    #numpy.save(file="/short/v10/tmp/ula3_tests/nbar/work/cs_solar_data.npy", arr=solar_angle_data)
    #numpy.save(file="/short/v10/tmp/ula3_tests/nbar/work/cs_sazi_data.npy", arr=sazi_angle_data)

    bounds = get_bounds(shape_dataset)

    # x & y pixel resolution (This should handle cases of non-square pixels.)
    dresx = abs(bounds.RasterXCellSize)
    dresy = abs(bounds.RasterYCellSize)

    #print "bounds = %s" % str(bounds)
    #print "dem_data.shape = %s" % str(dem_data.shape)
    #print "solar_angle_data.shape = %s" % str(solar_angle_data.shape)
    #print "sazi_angle_data.shape = %s" % str(sazi_angle_data.shape)
    #print "as_array(dem_data, dtype=numpy.float32).shape = %s" % str(as_array(dem_data, dtype=numpy.float32).shape)

    ierr, mask_all = shade_main_landsat_pixel(
        as_array(dem_data, dtype=numpy.float32),
        as_array(solar_angle_data, dtype=numpy.float32),
        as_array(sazi_angle_data, dtype=numpy.float32),
        dresx,
        dresy,
        spheroid,
        bounds.RasterYOrigin,
        bounds.RasterXOrigin,
        pix_buf.left,
        pix_buf.right,
        pix_buf.top,
        pix_buf.bottom,
        block_height,
        block_width,
        is_utm)

    if ierr:
        raise CastShadowError(ierr)

    return mask_all







"""
@print_call(logger.info)
def run_brdfterrain(
    rori, # threshold for terrain correction
    brdf0, brdf1, brdf2, # BRDF parameters
    bias, slope_ca, esun, dd, # satellite calibration coefficients
    ref_adj, # average reflectance for terrain correction
    #line,
    istart,
    #imid,
    iend,
    #ii,
    dn_1, # raw image
    mask_self, # mask
    mask_castsun, # self shadow mask
    mask_castview, # cast shadow mask
    solar_angle, # solar zenith angle
    sazi_angle, # solar azimuth angle
    view_angle, # view angle (for flat surface)
    rela_angle, # relative azimuth angle (for flat surface)
    slope_angle, # slop angle
    aspect_angle, # aspect angle
    it_angle, # incident angle (for inclined surface)
    et_angle, # exiting angle (for inclined surface)
    rela_slope, # relative angle (for inclined surface)
    a_mod, # MODTRAN output (a)
    b_mod, # MODTRAN output (b)
    s_mod, # MODTRAN output (s)
    fs, # MODTRAN output (fs)
    fv, # MODTRAN output (fv)
    ts, # MODTRAN output (ts)
    edir_h, # MODTRAN output (direct irradiance)
    edif_h # MODTRAN output (diffuse irradiance)
    ):
    """"""
    BRDF correction including terrain correction. This code is an interface to the fortran code brdf_terrain_newdiff.f90
    (which is compiled to a Python module using F2py). The parameters have the same names as those used in that code...
    so please see Fuqin for information on what they mean!

    :param rori:
        (type: float) Threshold for terrain correction.
    :type rori:
        float

    :param brdf0:
        (type: float) BRDF parameter.
    :type brdf0:
        float

    :param brdf1:
        (type: float) BRDF parameter.
    :type brdf1:
        float

    :param brdf2:
        (type: float) BRDF parameter.
    :type brdf2:
        float

    :param bias:
        (type: float) Satellite calibration coefficient.
    :type bias:
        float

    :param slope_ca:
        (type: float) Satellite calibration coefficient.
    :type slope_cs:
        float

    :param esun:
        (type: float) Satellite calibration coefficient.
    :param esun:
        float

    :param dd:
        (type: float) Satellite calibration coefficients.
    :type dd:
        float

    :param ref_adj:
        (type: float) Average reflectance for terrain correction.
    :type ref_adj:
        float

    :param istart:
        ???
    :type istart:
        One dimensional :py:class:`numpy.ndarray` which can be cast to type :py:const:`numpy.int4` with length equal to the
        number of rows in ``dn_1``.

    :param iend:
        ???
    :type iend:
        One dimensional :py:class:`numpy.ndarray` which can be cast to type :py:const:`numpy.int4` with length equal to the
        number of rows in ``dn_1``.

    :param dn_1:
        Raw image data.
    :type dn_1:
        Two dimensional :py:class:`numpy.ndarray` which can be cast to type :py:const:`numpy.int8`. The dimensions are
        unspecified and are used to determine the dimensions of ``istart``, ``iend`` the remaining (all following)
        array arguments.

    :param mask_self:
        Mask of pixels where the incident angle is greater than 90 degrees. These pixels are excluded as there is no
        illumination of the scene at these locations.
    :type mask_self:
        Array with the same dimensions as ``dn_1`` which can be cast to type :py:const:`numpy.int16`.

    :param mask_castsun:
        Mask of pixels which are shaded by other objects. These pixels are excluded as there is no illumination of the
        scene at these locations.
    :type mask_castsun:
        Array with the same dimensions as ``dn_1`` which can be cast to type :py:const:`numpy.int16`.

    :param mask_castview:
        Mask of pixels which are not visible to the satelite.
    :type mask_castview:
        Array with the same dimensions as ``dn_1`` which can be cast to type :py:const:`numpy.int16`.

    :param solar_angle:
        The solar zenith angle.
    :type solar_angle:
        Array with the same dimensions as ``dn_1`` which can be cast to type :py:const:`numpy.float32`.

    :param sazi_angle:
        solar azimuth angle.
    :type sazi_angle:
        Array with the same dimensions as ``dn_1`` which can be cast to type :py:const:`numpy.float32`.

    :param view_angle:
        view angle (for flat surface).
    :type view_angle:
        Array with the same dimensions as ``dn_1`` which can be cast to type :py:const:`numpy.float32`.

    :param rela_angle:
        relative azimuth angle (for flat surface).
    :type rela_angle:
        Array with the same dimensions as ``dn_1`` which can be cast to type :py:const:`numpy.float32`.

    :param slope_angle:
        slope angle.
    :type slope_angle:
        Array with the same dimensions as ``dn_1`` which can be cast to type :py:const:`numpy.float32`.

    :param aspect_angle:
        aspect angle.
    :type aspect_angle:
        Array with the same dimensions as ``dn_1`` which can be cast to type :py:const:`numpy.float32`.

    :param it_angle:
        incident angle (for inclined surface).
    :type it_angle:
        Array with the same dimensions as ``dn_1`` which can be cast to type :py:const:`numpy.float32`.

    :param et_angle:
        exiting angle (for inclined surface).
    :type et_angle:
        Array with the same dimensions as ``dn_1`` which can be cast to type :py:const:`numpy.float32`.

    :param rela_slope:
        relative angle (for inclined surface).
    :type rela_slope:
        Array with the same dimensions as ``dn_1`` which can be cast to type :py:const:`numpy.float32`.

    :param a_mod:
        MODTRAN output (a).
    :type a_mod:
        Array with the same dimensions as ``dn_1`` which can be cast to type :py:const:`numpy.float32`.

    :param b_mod:
        MODTRAN output (b).
    :type b_mod:
        Array with the same dimensions as ``dn_1`` which can be cast to type :py:const:`numpy.float32`.

    :param s_mod:
        MODTRAN output (s).
    :type s_mod:
        Array with the same dimensions as ``dn_1`` which can be cast to type :py:const:`numpy.float32`.

    :param fs:
        MODTRAN output (fs).
    :type fs:
        Array with the same dimensions as ``dn_1`` which can be cast to type :py:const:`numpy.float32`.

    :param fv:
        MODTRAN output (fv).
    :type fv:
        Array with the same dimensions as ``dn_1`` which can be cast to type :py:const:`numpy.float32`.

    :param ts:
        MODTRAN output (ts).
    :type ts:
        Array with the same dimensions as ``dn_1`` which can be cast to type :py:const:`numpy.float32`.

    :param edir_h:
        MODTRAN output (direct irradiance).
    :type edir_h:
        Array with the same dimensions as ``dn_1`` which can be cast to type :py:const:`numpy.float32`.

    :param edif_h:
        MODTRAN output (diffuse irradiance).
    :type edif_h:
        Array with the same dimensions as ``dn_1`` which can be cast to type :py:const:`numpy.float32`.

    Parameters: ``mask_self``, ``mask_castsun``, ``mask_castview`` can be generated using :py:func:`run_castshadow`.

    Parameters ``mask_self``, ``slope_angle``, ``aspect_angle``, ``it_angle``, ``et_angle``, and ``rela_slope``
    can be generated using the function :py:func:`run_slope`.

    All parameters after ``ref_adj`` are passed through :py:func:`ULA3.utils.as_array` with the appropriate argument
    type and hence, can have types that will work as arguments to that function; inparticular, they can be
    :py:class:`gdal.Dataset`s or paths to files that can be opened using :py:func:`gdal.Open`.

    :return: A tuple of three :py:class:`numpy.ndarray`s:

        - (index 0) Atmospheric corrected lambertial reflectance,

        - (index 1) Atmospheric and brdf corrected reflectance, and

        - (index 2) Atmospheric and brdf and terrain corrected reflectance

    :todo:
        This documentation should be reviewed by someone whom understands the process more thoroughly, and better
        descriptions of the arguments provided.

    """"""
    return terrain_correction(
        rori,
        brdf0, brdf1, brdf2,
        bias, slope_ca, esun, dd,
        ref_adj,
        as_array(istart, dtype=numpy.int32),
        as_array(iend, dtype=numpy.int32),
        as_array(dn_1, dtype=numpy.int8),
        as_array(mask_self, dtype=numpy.int16),
        as_array(mask_castsun, dtype=numpy.int16),
        as_array(mask_castview, dtype=numpy.int16),
        as_array(solar_angle, dtype=numpy.float32),
        as_array(sazi_angle, dtype=numpy.float32),
        as_array(view_angle, dtype=numpy.float32),
        as_array(rela_angle, dtype=numpy.float32),
        as_array(slope_angle, dtype=numpy.float32),
        as_array(aspect_angle, dtype=numpy.float32),
        as_array(it_angle, dtype=numpy.float32),
        as_array(et_angle, dtype=numpy.float32),
        as_array(rela_slope, dtype=numpy.float32),
        as_array(a_mod, dtype=numpy.float32),
        as_array(b_mod, dtype=numpy.float32),
        as_array(s_mod, dtype=numpy.float32),
        as_array(fs, dtype=numpy.float32),
        as_array(fv, dtype=numpy.float32),
        as_array(ts, dtype=numpy.float32),
        as_array(edir_h, dtype=numpy.float32),
        as_array(edif_h, dtype=numpy.float32))
"""

@print_call(logger.info)
def run_brdfterrain(
    rori, # threshold for terrain correction
    brdf0, brdf1, brdf2, # BRDF parameters
    bias, slope_ca, esun, dd, # satellite calibration coefficients
    ref_adj, # average reflectance for terrain correction
    dn_1, # raw image
    mask_self, # mask
    mask_castsun, # self shadow mask
    mask_castview, # cast shadow mask
    solar_angle, # solar zenith angle
    sazi_angle, # solar azimuth angle
    view_angle, # view angle (for flat surface)
    rela_angle, # relative azimuth angle (for flat surface)
    slope_angle, # slop angle
    aspect_angle, # aspect angle
    it_angle, # incident angle (for inclined surface)
    et_angle, # exiting angle (for inclined surface)
    rela_slope, # relative angle (for inclined surface)
    a_mod, # MODTRAN output (a)
    b_mod, # MODTRAN output (b)
    s_mod, # MODTRAN output (s)
    fs, # MODTRAN output (fs)
    fv, # MODTRAN output (fv)
    ts, # MODTRAN output (ts)
    edir_h, # MODTRAN output (direct irradiance)
    edif_h # MODTRAN output (diffuse irradiance)
    ):
    """
    BRDF correction including terrain correction. This code is an interface to the fortran code brdf_terrain_newdiff_LS8.f90
    (which is compiled to a Python module using F2py). The parameters have the same names as those used in that code...
    so please see Fuqin for information on what they mean!

    :param rori:
        (type: float) Threshold for terrain correction.
    :type rori:
        float

    :param brdf0:
        (type: float) BRDF parameter.
    :type brdf0:
        float

    :param brdf1:
        (type: float) BRDF parameter.
    :type brdf1:
        float

    :param brdf2:
        (type: float) BRDF parameter.
    :type brdf2:
        float

    :param bias:
        (type: float) Satellite calibration coefficient.
    :type bias:
        float

    :param slope_ca:
        (type: float) Satellite calibration coefficient.
    :type slope_cs:
        float

    :param esun:
        (type: float) Satellite calibration coefficient.
    :param esun:
        float

    :param dd:
        (type: float) Satellite calibration coefficients.
    :type dd:
        float

    :param ref_adj:
        (type: float) Average reflectance for terrain correction.
    :type ref_adj:
        float

    :param istart:
        ???
    :type istart:
        One dimensional :py:class:`numpy.ndarray` which can be cast to type :py:const:`numpy.int4` with length equal to the
        number of rows in ``dn_1``.

    :param iend:
        ???
    :type iend:
        One dimensional :py:class:`numpy.ndarray` which can be cast to type :py:const:`numpy.int4` with length equal to the
        number of rows in ``dn_1``.

    :param dn_1:
        Raw image data.
    :type dn_1:
        Two dimensional :py:class:`numpy.ndarray` which can be cast to type :py:const:`numpy.int8`. The dimensions are
        unspecified and are used to determine the dimensions of ``istart``, ``iend`` the remaining (all following)
        array arguments.

    :param mask_self:
        Mask of pixels where the incident angle is greater than 90 degrees. These pixels are excluded as there is no
        illumination of the scene at these locations.
    :type mask_self:
        Array with the same dimensions as ``dn_1`` which can be cast to type :py:const:`numpy.int16`.

    :param mask_castsun:
        Mask of pixels which are shaded by other objects. These pixels are excluded as there is no illumination of the
        scene at these locations.
    :type mask_castsun:
        Array with the same dimensions as ``dn_1`` which can be cast to type :py:const:`numpy.int16`.

    :param mask_castview:
        Mask of pixels which are not visible to the satelite.
    :type mask_castview:
        Array with the same dimensions as ``dn_1`` which can be cast to type :py:const:`numpy.int16`.

    :param solar_angle:
        The solar zenith angle.
    :type solar_angle:
        Array with the same dimensions as ``dn_1`` which can be cast to type :py:const:`numpy.float32`.

    :param sazi_angle:
        solar azimuth angle.
    :type sazi_angle:
        Array with the same dimensions as ``dn_1`` which can be cast to type :py:const:`numpy.float32`.

    :param view_angle:
        view angle (for flat surface).
    :type view_angle:
        Array with the same dimensions as ``dn_1`` which can be cast to type :py:const:`numpy.float32`.

    :param rela_angle:
        relative azimuth angle (for flat surface).
    :type rela_angle:
        Array with the same dimensions as ``dn_1`` which can be cast to type :py:const:`numpy.float32`.

    :param slope_angle:
        slope angle.
    :type slope_angle:
        Array with the same dimensions as ``dn_1`` which can be cast to type :py:const:`numpy.float32`.

    :param aspect_angle:
        aspect angle.
    :type aspect_angle:
        Array with the same dimensions as ``dn_1`` which can be cast to type :py:const:`numpy.float32`.

    :param it_angle:
        incident angle (for inclined surface).
    :type it_angle:
        Array with the same dimensions as ``dn_1`` which can be cast to type :py:const:`numpy.float32`.

    :param et_angle:
        exiting angle (for inclined surface).
    :type et_angle:
        Array with the same dimensions as ``dn_1`` which can be cast to type :py:const:`numpy.float32`.

    :param rela_slope:
        relative angle (for inclined surface).
    :type rela_slope:
        Array with the same dimensions as ``dn_1`` which can be cast to type :py:const:`numpy.float32`.

    :param a_mod:
        MODTRAN output (a).
    :type a_mod:
        Array with the same dimensions as ``dn_1`` which can be cast to type :py:const:`numpy.float32`.

    :param b_mod:
        MODTRAN output (b).
    :type b_mod:
        Array with the same dimensions as ``dn_1`` which can be cast to type :py:const:`numpy.float32`.

    :param s_mod:
        MODTRAN output (s).
    :type s_mod:
        Array with the same dimensions as ``dn_1`` which can be cast to type :py:const:`numpy.float32`.

    :param fs:
        MODTRAN output (fs).
    :type fs:
        Array with the same dimensions as ``dn_1`` which can be cast to type :py:const:`numpy.float32`.

    :param fv:
        MODTRAN output (fv).
    :type fv:
        Array with the same dimensions as ``dn_1`` which can be cast to type :py:const:`numpy.float32`.

    :param ts:
        MODTRAN output (ts).
    :type ts:
        Array with the same dimensions as ``dn_1`` which can be cast to type :py:const:`numpy.float32`.

    :param edir_h:
        MODTRAN output (direct irradiance).
    :type edir_h:
        Array with the same dimensions as ``dn_1`` which can be cast to type :py:const:`numpy.float32`.

    :param edif_h:
        MODTRAN output (diffuse irradiance).
    :type edif_h:
        Array with the same dimensions as ``dn_1`` which can be cast to type :py:const:`numpy.float32`.

    Parameters: ``mask_self``, ``mask_castsun``, ``mask_castview`` can be generated using :py:func:`run_castshadow`.

    Parameters ``mask_self``, ``slope_angle``, ``aspect_angle``, ``it_angle``, ``et_angle``, and ``rela_slope``
    can be generated using the function :py:func:`run_slope`.

    All parameters after ``ref_adj`` are passed through :py:func:`ULA3.utils.as_array` with the appropriate argument
    type and hence, can have types that will work as arguments to that function; inparticular, they can be
    :py:class:`gdal.Dataset`s or paths to files that can be opened using :py:func:`gdal.Open`.

    :return: A tuple of three :py:class:`numpy.ndarray`s:

        - (index 0) Atmospheric corrected lambertial reflectance,

        - (index 1) Atmospheric and brdf corrected reflectance, and

        - (index 2) Atmospheric and brdf and terrain corrected reflectance

    :todo:
        This documentation should be reviewed by someone whom understands the process more thoroughly, and better
        descriptions of the arguments provided.

    """
    return terrain_correction(
        rori,
        brdf0, brdf1, brdf2,
        bias, slope_ca, esun, dd,
        ref_adj,
        as_array(dn_1, dtype=numpy.int16),
        as_array(mask_self, dtype=numpy.int16),
        as_array(mask_castsun, dtype=numpy.int16),
        as_array(mask_castview, dtype=numpy.int16),
        as_array(solar_angle, dtype=numpy.float32),
        as_array(sazi_angle, dtype=numpy.float32),
        as_array(view_angle, dtype=numpy.float32),
        as_array(rela_angle, dtype=numpy.float32),
        as_array(slope_angle, dtype=numpy.float32),
        as_array(aspect_angle, dtype=numpy.float32),
        as_array(it_angle, dtype=numpy.float32),
        as_array(et_angle, dtype=numpy.float32),
        as_array(rela_slope, dtype=numpy.float32),
        as_array(a_mod, dtype=numpy.float32),
        as_array(b_mod, dtype=numpy.float32),
        as_array(s_mod, dtype=numpy.float32),
        as_array(fs, dtype=numpy.float32),
        as_array(fv, dtype=numpy.float32),
        as_array(ts, dtype=numpy.float32),
        as_array(edir_h, dtype=numpy.float32),
        as_array(edif_h, dtype=numpy.float32))

