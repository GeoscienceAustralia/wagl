from gaip.__bilinear_interpolate import bilinear_interpolation as fortran
from scipy.interpolate import Rbf
import numpy as np
import math


def fortran_bilinear_interpolate(cols, rows, locations, samples,
                                 row_start, row_end, row_centre):
    """
    Original NBAR interpolation scheme.
    Sheared 4-cell bilinear, implemented in fortran.
    """
    assert len(samples) == 3*3
    assert len(locations) == len(samples)

    s1 = samples[[0,1,3,4]]
    s2 = samples[[1,2,4,5]]
    s3 = samples[[3,4,6,7]]
    s4 = samples[[4,5,7,8]]

    output = np.empty((rows, cols))

    fortran(rows, cols, locations, s1, s2, s3, s4,
            row_start-1, row_end-1, row_centre-1, output.T)

    return output


def rbf_interpolate(cols, rows, locations, samples, _, _, _):
    """scipy linear radial basis function interpolation"""

    rbf = Rbf(locations[:,0], locations[:,1], samples, function='linear')

    return rbf(*np.mgrid[:rows, :cols])


def sheared_bilinear_interpolate(cols, rows, locations, samples,
                                 row_start, row_end, row_centre):
    """
    Generalisation of the original NBAR interpolation scheme
    """
    raise NotImplementedError


interpolate = fortran_bilinear_interpolate