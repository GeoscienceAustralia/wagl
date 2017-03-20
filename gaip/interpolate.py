from scipy.interpolate import Rbf
import numpy as np
import math


def fortran_bilinear_interpolate(cols, rows, locations, samples,
                                 row_start, row_end, row_centre):
    """
    Original NBAR interpolation scheme.
    Sheared 4-cell bilinear, implemented in fortran.
    """
    from gaip.__bilinear_interpolate import bilinear_interpolation as fortran

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


def rbf_interpolate(cols, rows, locations, samples, *_):
    """scipy linear radial basis function interpolation"""

    rbf = Rbf(locations[:,0], locations[:,1], samples, function='linear')

    return rbf(*np.mgrid[:rows, :cols])


def sheared_bilinear_interpolate(cols, rows, locations, samples,
                                 row_start, row_end, row_centre):
    """
    Generalisation of the original NBAR interpolation scheme
    """
    #raise NotImplementedError

    n = len(samples)
    grid_size = int(math.sqrt(n)) - 1

    print(grid_size,n)

    assert (grid_size+1)**2 == n and not (grid_size & 1) and not grid_size % 1
    # Assume count of samples is 9 or 25, 49, 81.. (Grid size is 2, 4, 6, ..)

    vertex_shape = (grid_size+1,)*2
    locations = locations.reshape(vertex_shape+(2,))
    samples = samples.reshape(vertex_shape)

    once = lambda f: f() # only for code structure
    @once
    def lines():
        """Place parcel boundaries alongside track (by 1D linear interpolation)"""
        L = np.empty((grid_size+1, rows), dtype=np.uint64)

        middle_vertex = grid_size//2

        L[0] = row_start
        L[middle_vertex] = row_centre
        L[-1] = row_end

        for i in range(1, middle_vertex):
            L[i] = row_start + (row_centre - row_start) * (i/middle_vertex)
            L[i+middle_vertex] = row_centre + (row_end - row_centre) * (i/middle_vertex)

        return L

    x, y = np.ogrid[:rows, :cols]

    @once
    def parcellation():
        """Generate parcellation map"""
        # Would use e.g. a matplotlib poly drawing routine,
        # except if curved sides are required.

        zones = np.full((rows, cols), 0, dtype=np.int8)

        # first axis
        for line in lines[1:-1]:
            zones += (x >= line)
            # note, needn't be concerned with how edges are zoned
            # because they will be masked out downstream

        # second axis
        for line in locations[1:-1,0,1]: # y component of left-edge samples
            zones += (y >= line) * grid_size

        return zones

    def shear(i, j, both_sides=False):
        """Warp to straighten edges of trapezoid"""
        left = lines[j]
        width = lines[j+1] - left

        xx = x - left

        return xx if not both_sides else xx / width

    def patch(i, j, x=x, shear=False):
        """bilinear cell"""
        vertices = locations[i:i+1, j:j+1].reshape(4, 2)
        values = samples[i:i+1, j:j+1].reshape(4)

        if shear: # then re-map coordinates
            x = shear(i, j)
            vertices = vertices.copy()
            vertices[:,0] = x[vertices[:,0]]

        matrix = np.ones((4,4))
        matrix[:, 1:3] = vertices
        matrix[:, 3] = vertices[:,0] * vertices[:,1]

        a = np.linalg.solve(matrix, values) # determine coefficients

        return a[0] + a[1]*x + a[2]*y + a[3]*x*y # broadcast as raster

    result = np.full((rows, cols), np.nan, dtype=np.float)
    for i in range(grid_size):
        for j in range(grid_size):
            subset = parcellation == i*grid_size + j
            result[subset] = patch(i, j)[subset]

    return result


interpolate = fortran_bilinear_interpolate

