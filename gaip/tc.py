import os
import numpy
from scipy import ndimage
from gaip import shade_main_landsat_pixel
from gaip import slope_pixelsize_newpole
from gaip import terrain_correction
from gaip import as_array
from gaip import write_img


def filter_dsm(array):
    """
    Applies a gaussian filter to array.

    :param array:
        A 2D NumPy array.

    :return:
        A 2D NumPy array.
    """
    # Define the kernel
    kernel = [0.009511, 0.078501, 0.009511, 0.078501, 0.647954, 0.078501,
        0.009511, 0.078501, 0.009511]
    kernel = numpy.array(kernel).reshape((3,3))

    filtered = ndimage.convolve(array, kernel)
    return filtered


class Buffers(object):
    """
    Holds some value for each side of an image. This was initially
    created to hold buffer widths (in pixels) for a scene, but does
    not care about the type of the values passed to the constructer
    and can hence be used for any type.
    """
    def __init__(self, left, right=None, top=None, bottom=None):
        """
        Constructor.

        The arguments are copied directly to members of the same name
        with the restrictions that:

        - if right is None then top and bottom must also be None, and

        - if right is not None then top and bottom must not be none
          either.
        """
        if right is None:
            msg = "if right is None then top and bottom must also be None"
            assert top is None and bottom is None, msg
            self.left = self.right = self.top = self.bottom = left
        else:
            msg = ("if right is not None then top and bottom must "
                   "also not be None")
            assert top is not None and bottom is not None, msg
            self.left = left
            self.right = right
            self.top = top
            self.bottom = bottom

    def __str__(self):
        msg = "Buffers({left}, {right}, {top}, {bottom})"
        msg = msg.format(left=self.left, right=self.right, top=self.top,
            bottom=self.bottom)
        return msg


def write_new_brdf_file(file_name, *args):
    with open(file_name, 'w') as src:
        out_string = "{0}\n{1} {2} {3}\n{4} {5} {6} {7}\n{8}\n"
        out_string = out_string.format(*args)
        src.write(out_string)


def write_header_slope_file(file_name, bounds, geobox):
    with open(file_name, 'w') as output:
        # get dimensions, resolution and pixel origin
        rows, cols = geobox.shape
        res = geobox.res
        origin = geobox.origin

        # Now output the details
        output.write("{nr} {nc}\n".format(nr=rows, nc=cols))
        output.write("{0} {1}\n{2} {3}\n".format(bounds.left, bounds.right,
            bounds.top, bounds.bottom))
        output.write("{resy} {resx}\n".format(resx=res[0], resy=res[1]))
        output.write("{yorigin} {xorigin}\n".format(xorigin=origin[0],
            yorigin=origin[1]))


class FortranError(Exception):
    """
    Base class for errors thrown from the Fortran code used in this module.
    """
    def __init__(self, function_name, code, msg):
        # The name of the Fortran function called.
        self.function_name = function_name
        self.code = code # The error code produced by the Fortran code.
        # The Message corresponding to ``code``.
        self.msg = msg or "Unknown error"

    def __str__(self):
        """
        Return a string representation of this Error.
        """
        err = "Error in Fotran code {0} (code {1}): {2}"
        err = err.format(self.function_name, self.code, self.msg)
        return err


class SlopeResultSet(object):
    """
    Holds the results of a call to :py:func:`run_slope`.
    """
    def __init__(self, mask_self, slope, aspect, incident, exiting,
        azi_incident, azi_exiting, rela_slope):
        """
        All arguments are :py:class:`numpy.ndarray`s. These
        correspond to the arguments I have no idea what these
        actually represent, so someone who does should document this
        (Fuqin will know).
        """
        self.mask_self = mask_self
        self.slope = slope
        self.aspect = aspect
        self.incident = incident
        self.exiting = exiting
        self.azi_incident = azi_incident
        self.azi_exiting = azi_exiting
        self.rela_slope = rela_slope

    def write_arrays(self, output_path, geobox, file_type='ENVI',
        file_extension='.img'):

        # Filenames
        fname_mask_self = os.path.join(output_path, 'mask_self' +
            file_extension)
        fname_slope = os.path.join(output_path, 'slope' + file_extension)
        fname_aspect = os.path.join(output_path, 'aspect' + file_extension)
        fname_incident = os.path.join(output_path, 'incident' +
            file_extension)
        fname_exiting = os.path.join(output_path, 'exiting' + file_extension)
        fname_azimuth_incident = os.path.join(output_path, 'azi_incident' +
            file_extension)
        fname_azimuth_exiting = os.path.join(output_path, 'azi_exiting' +
            file_extension)
        fname_relative_slope = os.path.join(output_path, 'rela_slope' +
            file_extension)

        # Write
        write_img(self.mask_self, fname_mask_self, format=file_type,
            geobox=geobox, nodata=-999)
        write_img(self.slope, fname_slope, format=file_type, geobox=geobox,
            nodata=-999)
        write_img(self.aspect, fname_aspect, format=file_type, geobox=geobox,
            nodata=-999)
        write_img(self.incident, fname_incident, format=file_type,
            geobox=geobox, nodata=-999)
        write_img(self.exiting, fname_exiting, format=file_type,
            geobox=geobox, nodata=-999)
        write_img(self.azi_incident, fname_azimuth_incident, format=file_type,
            geobox=geobox, nodata=-999)
        write_img(self.azi_exiting, fname_azimuth_exiting, format=file_type,
            geobox=geobox, nodata=-999)
        write_img(self.rela_slope, fname_relative_slope, format=file_type,
            geobox=geobox, nodata=-999)


class SlopeError(FortranError):
    """
    Class that deals with errors from :py:func:`run_slope`.
    """
    def __init__(self, code):
        super(SlopeError, self).__init__("slope_pixelsize_newpole", code,
            SlopeError.get_error_message(code))

    @staticmethod
    def get_error_message(code):
        """
        Gennerate an error message for a specific code (It is OK
        that this have non-returning control paths).
        """
        if code in 10:
            return "X dimensions of scene and DEM not correct."
        if code == 11:
            return "Y dimensions of scene and DEM not correct."


def run_slope(
    acquisition,
    DEM,
    solar_zenith,
    satellite_zenith,
    solar_azimuth,
    satellite_azimuth,
    buffer,
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
    Calculate the slope and angles for a region. This code is an
    interface to the fortran code slope_pixel_newpole.f90 written by
    Fuqin (which was modified to work with F2py).

    The following was taken from the top of the Fotran program:
    "slope_pixelsize_newpole.f90:
    This program is used to calculate slope and aspect angles
    using Sobel filter and then calculate incident and
    exiting angles as well as their azimuth angles.
    note: the row and column of DEM data must be larger
    than the image (extra each line and column for the four sides.
    it is needed for sobel filter.

    :param acquisition:
        An instance of an acquisition object.

    :param DEM:
        A DEM of the region. This must have the same dimensions as
        zenith_angle plus a buffer of widths specified by buffer.

    :param solar_zenith:
        The solar zenith angle data for the region.

    :param satellite_zenith:
        The satellite zenith angle data for the region.

    :param solar_azimuth:
        The solar azimuth angle data for the region.

    :param satellite_azimuth:
        The satellite azimuth angle data for the region.

    :param buffer:
        An object with members top, bottom, left and right giving the
        size of the buffer (in pixels) which have been added to the
        corresponding sides of DEM.

    :param is_utm:
        Boolean specifying whether the data is in UTM coordinates.
        This is passed to the Fortran code which uses it to determine
        whether the edge length of the pixels (if it is not true,
        then the resolution is assumed to be in degrees, and the
        pixel size is calculated internally).

    :param spheroid:
        A 4 element floating point array containing the Earth
        spheroidal paramaters.
        Index 0 contains the spheroid Major Axis.
        Index 1 contains the spheroid Inverse Flattening.
        Index 2 contains the spheroid Squared Eccentricity.
        Index 3 contains the Earth rotational angular velocity in
        radians/second.

    :return:
        A SlopeResultSet Class with the following NumPy 2D arrays:
        mask_self
        slope
        aspect
        incident
        exiting
        azi_incident
        azi_exiting
        rela_slope

    :notes:
    The parameters ``solar_zenith, satellite_zenith, solar_azimuth``
    and ``satellite_azimuth_data`` must have the same dimensions.
    """
    # Perform datatype checks
    if DEM.dtype.name != 'float32':
        msg = 'DEM datatype must be float32! Datatype: {dtype}'
        msg = msg.format(dtype=DEM.dtype.name)
        raise TypeError(msg)

    if solar_zenith.dtype.name != 'float32':
        msg = 'Solar zenith datatype must be float32! Datatype: {dtype}'
        msg = msg.format(dtype=solar_zenith.dtype.name)
        raise TypeError(msg)

    if satellite_zenith.dtype.name != 'float32':
        msg = 'Satellite zenith datatype must be float32! Datatype: {dtype}'
        msg = msg.format(dtype=satellite_zenith.dtype.name)

    if solar_azimuth.dtype.name != 'float32':
        msg = 'Solar azimuth datatype must be float32! Datatype: {dtype}'
        msg = msg.format(dtype=solar_azimuth.dtype.name)
        raise TypeError(msg)

    if satellite_azimuth.dtype.name != 'float32':
        msg = 'Satellite azimuth datatype must be float32! Datatype: {dtype}'
        msg = msg.format(dtype=satellite_azimuth.dtype.name)
        raise TypeError(msg)

    # Get the x and y pixel sizes
    geobox = acquisition.gridded_geo_box()
    x_origin, y_origin = geobox.origin
    x_res, y_res = geobox.pixelsize
    dresx = x_res + 2
    dresy = y_res + 2

    # Get acquisition dimensions and add 1 pixel top, bottom, left & right
    cols, rows = geobox.getShapeXY()
    ncol = cols + 2
    nrow = rows + 2

    dem_dat = DEM[(buffer.top-1):-(buffer.bottom-1),(buffer.left-1):
        -(buffer.right-1)]

    # Check that the dimensions match
    if dem_dat.shape != (nrow, ncol):
        msg = ('DEM index not of correct shape ({row}, {col}) '
              '!= ({drow}, {dcol})')
        msg = msg.format(row=nrow, col=ncol, drow=dem_dat.shape[0],
            dcol=dem_dat.shape[1])
        raise IndexError(msg)

    # This will be ignored if is_utm == True
    alat = numpy.array([y_origin-i*dresy for i in range(-1, nrow-1)],
        dtype=numpy.float64) # yes, I did mean float64.

    (mask, theta, phit, it, et, azi_it,
     azi_et, rela, ierr) = slope_pixelsize_newpole(
        dresx, dresy, spheroid, alat, is_utm,
        dem_dat,
        solar_zenith,
        satellite_zenith,
        solar_azimuth,
        satellite_azimuth)

    if ierr:
        raise SlopeError(ierr)

    slope_results_set = SlopeResultSet(mask_self=mask, slope=theta,
        aspect=phit, incident=it, exiting=et, azi_incident=azi_it,
        azi_exiting=azi_et, rela_slope=rela)

    return slope_results_set


class CastShadowError(FortranError):
    """
    Class that deals with errors from :py:func:`run_castshadow`.
    """
    def __init__(self, code):
        super(CastShadowError, self).__init__("shade_main_landsat_pixel",
            code, CastShadowError.get_error_message(code))

    @staticmethod
    def get_error_message(code):
        """
        Generate an error message for a specific code. It is OK for this have
        non-returning control paths, as this will results in ``None``, which
        is handled in the super class.
        """
        def tmpt(d, n):
            err = "attempt to access invalid {0} of {1}".format(d, n)
            return err

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





def run_castshadow(
    acquisition,
    DEM,
    zenith_angle,
    azimuth_angle,
    buffer,
    block_height,
    block_width,
    spheroid):
    """
    This code is an interface to the fortran code
    shade_main_landsat_pixel.f90 written by Fuqin (and modified to
    work with F2py).

    The following was taken from the top of the Fotran program:
    "shade_main_landsat_pixel.f90":

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
    scene (see parameter ``pix_buf``. This will change with elevation
    difference within the scene and solar zenith angle. For
    Australian region and Landsat scene with 0.00025 degree
    resolution, the maximum extra lines are set to 250 pixels/lines
    for each direction. This figure shold be sufficient for everywhere
    and anytime in Australia. Thus the DEM image will be larger than
    landsat image for 500 lines x 500 columns

    :param acquisition:
        An instance of an acquisition object.
    :type acquisition:
        Class, Acquisition

    :param DEM:
        A DEM of the region. This must have the same dimensions as
        zenith_angle plus a buffer of widths specified by buffer.
    :type DEM:
        A 2D NumPy float32 array.

    :param zenith_angle:
        Array of zenith angles (in degrees). Must be of the same
        dimensions as azimuth_angle.
    :type zenith_angle:
        A 2D NumPy float32 array.

    :param azimuth_angle:
        Array of azimuth angles (in degrees). Must be of the same
        dimensions as zenith_angle.
    :type azimuth_angle:
        A 2D NumPy float32 array.

    :param buffer:
        Object describing the pixel buffers around the azimuth_angle
        and the zenith_angle arrays.
    :type buffer:
        Class, Buffers with properties left, right, top & bottom.

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

    :return:
        A 2D NumPy array containing the shadow mask.

    :warning:
        The Fortran code cannot be compiled with ``-O3`` as it
        produces incorrect results if it is.
    """
    # Get the x and y pixel sizes
    geobox = acquisition.gridded_geo_box()
    x_res, y_res = geobox.pixelsize
    x_origin, y_origin = geobox.origin

    # Are we in UTM or geographics?
    is_utm = not geobox.crs.IsGeographic()

    # Perform datatype checks
    if DEM.dtype.name != 'float32':
        msg = 'DEM datatype must be float32! Datatype: {dtype}'
        msg = msg.format(dtype=DEM.dtype.name)
        raise TypeError(msg)

    if zenith_angle.dtype.name != 'float32':
        msg = 'Zenith angle datatype must be float32! Datatype: {dtype}'
        msg = msg.format(dtype=zenith_angle.dtype.name)
        raise TypeError(msg)

    if azimuth_angle.dtype.name != 'float32':
        msg = 'Azimuth angle datatype must be float32! Datatype: {dtype}'
        msg = msg.format(dtype=azimuth_angle.dtype.name)
        raise TypeError(msg)

    ierr, mask = shade_main_landsat_pixel(DEM, zenith_angle, azimuth_angle,
        x_res, y_res, spheroid, y_origin, x_origin, buffer.left, buffer.right,
        buffer.top, buffer.bottom, block_height, block_width, is_utm)

    if ierr:
        raise CastShadowError(ierr)

    return mask








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
    BRDF correction including terrain correction. This code is an
    interface to the fortran code brdf_terrain_newdiff_LS8.f90
    (which is compiled to a Python module using F2py).
    The parameters have the same names as those used in that code...
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

    :param dn_1:
        Raw image data.
    :type dn_1:
        Two dimensional :py:class:`numpy.ndarray` of type
        :py:const:`numpy.int16`.

    :param mask_self:
        Mask of pixels where the incident angle is greater than 90
        degrees. These pixels are excluded as there is no illumination
        of the scene at these locations.
    :type mask_self:
        Array with the same dimensions as ``dn_1`` of type
        :py:const:`numpy.int16`.
        If not of the required datatype, then the input will be cast
        to the required datatype.

    :param mask_castsun:
        Mask of pixels which are shaded by other objects. These pixels
        are excluded as there is no illumination of the scene at these
        locations.
    :type mask_castsun:
        Array with the same dimensions as ``dn_1`` of type
        :py:const:`numpy.int16`.
        If not of the required datatype, then the input will be cast
        to the required datatype.

    :param mask_castview:
        Mask of pixels which are not visible to the satelite.
    :type mask_castview:
        Array with the same dimensions as ``dn_1`` of type
        :py:const:`numpy.int16`.
        If not of the required datatype, then the input will be cast
        to the required datatype.

    :param solar_angle:
        The solar zenith angle.
    :type solar_angle:
        Array with the same dimensions as ``dn_1``of type
        :py:const:`numpy.float32`.
        If not of the required datatype, then the input will be cast
        to the required datatype.

    :param sazi_angle:
        solar azimuth angle.
    :type sazi_angle:
        Array with the same dimensions as ``dn_1`` of type
        :py:const:`numpy.float32`.
        If not of the required datatype, then the input will be cast
        to the required datatype.

    :param view_angle:
        view angle (for flat surface).
    :type view_angle:
        Array with the same dimensions as ``dn_1`` of type
        :py:const:`numpy.float32`.
        If not of the required datatype, then the input will be cast
        to the required datatype.

    :param rela_angle:
        relative azimuth angle (for flat surface).
    :type rela_angle:
        Array with the same dimensions as ``dn_1`` of type
        :py:const:`numpy.float32`.
        If not of the required datatype, then the input will be cast
        to the required datatype.

    :param slope_angle:
        slope angle.
    :type slope_angle:
        Array with the same dimensions as ``dn_1`` of type
        :py:const:`numpy.float32`.
        If not of the required datatype, then the input will be cast
        to the required datatype.

    :param aspect_angle:
        aspect angle.
    :type aspect_angle:
        Array with the same dimensions as ``dn_1`` of type
        :py:const:`numpy.float32`.
        If not of the required datatype, then the input will be cast
        to the required datatype.

    :param it_angle:
        incident angle (for inclined surface).
    :type it_angle:
        Array with the same dimensions as ``dn_1`` of type
        :py:const:`numpy.float32`.
        If not of the required datatype, then the input will be cast
        to the required datatype.

    :param et_angle:
        exiting angle (for inclined surface).
    :type et_angle:
        Array with the same dimensions as ``dn_1`` of type
        :py:const:`numpy.float32`.
        If not of the required datatype, then the input will be cast
        to the required datatype.

    :param rela_slope:
        relative angle (for inclined surface).
    :type rela_slope:
        Array with the same dimensions as ``dn_1`` of type
        :py:const:`numpy.float32`.
        If not of the required datatype, then the input will be cast
        to the required datatype.

    :param a_mod:
        MODTRAN output (a).
    :type a_mod:
        Array with the same dimensions as ``dn_1`` of type
        :py:const:`numpy.float32`.
        If not of the required datatype, then the input will be cast
        to the required datatype.

    :param b_mod:
        MODTRAN output (b).
    :type b_mod:
        Array with the same dimensions as ``dn_1`` of type
        :py:const:`numpy.float32`.
        If not of the required datatype, then the input will be cast
        to the required datatype.

    :param s_mod:
        MODTRAN output (s).
    :type s_mod:
        Array with the same dimensions as ``dn_1`` of type
        :py:const:`numpy.float32`.
        If not of the required datatype, then the input will be cast
        to the required datatype.

    :param fs:
        MODTRAN output (fs).
    :type fs:
        Array with the same dimensions as ``dn_1`` of type
        :py:const:`numpy.float32`.
        If not of the required datatype, then the input will be cast
        to the required datatype.

    :param fv:
        MODTRAN output (fv).
    :type fv:
        Array with the same dimensions as ``dn_1`` of type
        :py:const:`numpy.float32`.
        If not of the required datatype, then the input will be cast
        to the required datatype.

    :param ts:
        MODTRAN output (ts).
    :type ts:
        Array with the same dimensions as ``dn_1`` of type
        :py:const:`numpy.float32`.
        If not of the required datatype, then the input will be cast
        to the required datatype.

    :param edir_h:
        MODTRAN output (direct irradiance).
    :type edir_h:
        Array with the same dimensions as ``dn_1`` of type
        :py:const:`numpy.float32`.
        If not of the required datatype, then the input will be cast
        to the required datatype.

    :param edif_h:
        MODTRAN output (diffuse irradiance).
    :type edif_h:
        Array with the same dimensions as ``dn_1`` of type
        :py:const:`numpy.float32`.
        If not of the required datatype, then the input will be cast
        to the required datatype.

    Parameters: ``mask_self``, ``mask_castsun``, ``mask_castview``
        can be generated using :py:func:`run_castshadow`.

    Parameters ``mask_self``, ``slope_angle``, ``aspect_angle``,
         ``it_angle``, ``et_angle``, and ``rela_slope`` can be
         generated using the function :py:func:`run_slope`.

    :notes:
        All parameters after ``ref_adj`` will be transposed to FORTRAN
        ordering. The same parameters will also be cast to the
        required datatype if necessary.

    :return:
        A tuple of three :py:class:`numpy.ndarray`s:

        - (index 0) Atmospheric corrected lambertial reflectance,

        - (index 1) Atmospheric and brdf corrected reflectance, and

        - (index 2) Atmospheric and brdf and terrain corrected
          reflectance
    """
    # TODO check for correct array datatypes
    # We want to limit any potential copies
    return terrain_correction(
        rori,
        brdf0, brdf1, brdf2,
        bias, slope_ca, esun, dd,
        ref_adj,
        as_array(dn_1, dtype=numpy.int16, transpose=True),
        as_array(mask_self, dtype=numpy.int16, transpose=True),
        as_array(mask_castsun, dtype=numpy.int16, transpose=True),
        as_array(mask_castview, dtype=numpy.int16, transpose=True),
        as_array(solar_angle, dtype=numpy.float32, transpose=True),
        as_array(sazi_angle, dtype=numpy.float32, transpose=True),
        as_array(view_angle, dtype=numpy.float32, transpose=True),
        as_array(rela_angle, dtype=numpy.float32, transpose=True),
        as_array(slope_angle, dtype=numpy.float32, transpose=True),
        as_array(aspect_angle, dtype=numpy.float32, transpose=True),
        as_array(it_angle, dtype=numpy.float32, transpose=True),
        as_array(et_angle, dtype=numpy.float32, transpose=True),
        as_array(rela_slope, dtype=numpy.float32, transpose=True),
        as_array(a_mod, dtype=numpy.float32, transpose=True),
        as_array(b_mod, dtype=numpy.float32, transpose=True),
        as_array(s_mod, dtype=numpy.float32, transpose=True),
        as_array(fs, dtype=numpy.float32, transpose=True),
        as_array(fv, dtype=numpy.float32, transpose=True),
        as_array(ts, dtype=numpy.float32, transpose=True),
        as_array(edir_h, dtype=numpy.float32, transpose=True),
        as_array(edif_h, dtype=numpy.float32, transpose=True))
