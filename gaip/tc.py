"""
Terrain Correction
------------------
"""
import numpy
from scipy import ndimage
from gaip import cast_shadow_main
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
    kernel = numpy.array(kernel).reshape((3, 3))

    filtered = ndimage.convolve(array, kernel)
    return filtered


class FortranError(Exception):

    """
    Base class for errors thrown from the Fortran code used in this module.
    """

    def __init__(self, function_name, code, msg):
        self.function_name = function_name
        self.code = code
        self.msg = msg or "Unknown error"

    def __str__(self):
        """
        Return a string representation of this Error.
        """
        err = "Error in Fotran code {0} (code {1}): {2}"
        err = err.format(self.function_name, self.code, self.msg)
        return err


class CastShadowError(FortranError):

    """
    Class that deals with errors from :py:func:`run_castshadow`.
    """

    def __init__(self, code):
        super(CastShadowError,
              self).__init__("cast_shadow_main",
                             code,
                             CastShadowError.get_error_message(code))

    @staticmethod
    def get_error_message(code):
        """
        Generate an error message for a specific code. It is OK for this have
        non-returning control paths, as this will results in ``None``, which
        is handled in the super class.
        """
        def tmpt(d, n):
            """Generate message."""
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
            return "matrix A does not have sufficient y margin"
        if code == 74:
            return "matrix A does not have sufficient x margin"


def run_castshadow(acquisition, dem, zenith_angle, azimuth_angle, margin,
                   block_height, block_width, spheroid):
    """
    This code is an interface to the fortran code
    cast_shadow_main.f90 written by Fuqin (and modified to
    work with F2py).

    The following was taken from the top of the Fotran program:
    "cast_shadow_main.f90":

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

    :param dem:
        A DEM of the region. This must have the same dimensions as
        zenith_angle plus a margin of widths specified by margin.
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

    :param margin:
        Object describing the pixel margins around the azimuth_angle
        and the zenith_angle arrays.
    :type margin:
        Class, ImageMargins with properties left, right, top & bottom.

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
    if dem.dtype.name != 'float32':
        msg = 'DEM datatype must be float32! Datatype: {dtype}'
        msg = msg.format(dtype=dem.dtype.name)
        raise TypeError(msg)

    if zenith_angle.dtype.name != 'float32':
        msg = 'Zenith angle datatype must be float32! Datatype: {dtype}'
        msg = msg.format(dtype=zenith_angle.dtype.name)
        raise TypeError(msg)

    if azimuth_angle.dtype.name != 'float32':
        msg = 'Azimuth angle datatype must be float32! Datatype: {dtype}'
        msg = msg.format(dtype=azimuth_angle.dtype.name)
        raise TypeError(msg)

    ierr, mask = cast_shadow_main(dem, zenith_angle, azimuth_angle, x_res,
                                  y_res, spheroid, y_origin, x_origin,
                                  margin.left, margin.right, margin.top,
                                  margin.bottom, block_height, block_width,
                                  is_utm)

    if ierr:
        raise CastShadowError(ierr)

    return mask
