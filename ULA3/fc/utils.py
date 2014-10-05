import os, errno, sys, numpy, numexpr
import unmiximage, endmembers
import logging

"""
Utility functions used in fractional cover. These should not be used outside of this package as they
may be moved, renamed or deleted in the future. If you need to add more functions that only for internal
use, this is the place to put them.
"""





def datatype(dtype):
    """
    Maps a gdal datatype to a numpy datatype.

    :param dtype:
        The gdal datatype to get the numpy datatype for.

    :return:
        The numpy datatype corresponding to ``gdal_data_type``.
    """
    instr = str(dtype)
    return {
        '1' : 'uint8',
        '2' : 'uint16',
        '3' : 'int16',
        '4' : 'uint32',
        '5' : 'int32',
        '6' : 'float32',
        '7' : 'float64',
        '8' : 'complex64',
        '9' : 'complex64',
        '10': 'complex64',
        '11': 'complex128',
        }.get(instr, 'float64')





def create_dir(path):
    """
    Create a directory.
    """
    try:
        os.makedirs(path)
    except OSError, e:
        if e.errno != errno.EEXIST:
            raise





def unmix(landsatReflectance):
    # NNLS Unmixing v1.0
    # Scarth 20090810 14:06:35 CEST
    # This implements a constrained unmixing process to recover the fraction images from
    # a synthetic reflectance generated from a large number of interactive
    # terms produced from the original and log-transformed landsat bands
    # It uses pre-defined endmembers that are defined below
    # Note that this relies on the more recent builds of Scipy that have the nnls.f code wrapped

    # The tiling routine that generates a list of tiles to process
    def get_tile2(array, xtile=100,ytile=100):
        dims = array.shape
        ncols = dims[1]
        nrows = dims[0]
        l = []
        if len(dims) >2:
            ncols = dims[2]
            nrows = dims[1]
            dims  = (nrows,ncols)
        xstart = numpy.arange(0,ncols,xtile)
        ystart = numpy.arange(0,nrows,ytile)
        for ystep in ystart:
            if ystep + ytile < nrows:
                yend = ystep + ytile
            else:
                yend = nrows
            for xstep in xstart:
                if xstep + xtile < ncols:
                    xend = xstep + xtile
                else:
                    xend = ncols
                l.append((ystep,yend,xstep,xend))
        return l


    # Define the weight of the sum to one constraint
    # This value determined how well the resulting fractions will sum to 100%
    # I typically determine this by running the unmixing against field data for a number of values, picking the best one
    sumToOneWeight = endmembers.sum_weight('2013_01_08')

    # 2009v gives green, dead, bare1 and bare2
    # 2012v gives green, dead, bare1 and bare2
    # 2013v gives green, dead1, dead2 and bare fractions
    # Note the last row is the sum to one constraint value
    endmembers_array = endmembers.endmember_version('2013_01_08')

    dims = landsatReflectance.shape
    if len(dims) >2:
        ncols = dims[2]
        nrows = dims[1]
        dims  = (nrows,ncols)

    # Output can be float32, but the input and calculations have to be
    # performed in float64, otherwise results will be zero.

    # For use with the fortran function
    fractions = numpy.zeros((5,dims[0],dims[1]), dtype='float32')

    # *******************************************
    # Use a tiling procedure, as the array that is created for input into the
    # unmixing algorithm is [55,y,x]. Not sure how QDERM evaluates their
    # components, but a [55,8000,8000] float64 is too large to be held in
    # memory all at once.
    # [58,y,x] array is calculated for 2013_01_08 version
    # *******************************************

    tiles = get_tile2(array=landsatReflectance)

    # Calculate the [55,y,x] array to feed into the unmixing algorithm
    # Loops over the number of tiles. Potentially could be run in parallel as
    # each tile is independent, and no tiles overlap.
    # [58,y,x] array is calculated for 2013_01_08 version
    for tile in tiles:
        ystart = tile[0]
        yend   = tile[1]
        xstart = tile[2]
        xend   = tile[3]

        subset = landsatReflectance[:,ystart:yend,xstart:xend]
        subset = numexpr.evaluate("(1.0 + subset) * 0.0001")

        band2 = subset[0]
        band3 = subset[1]
        band4 = subset[2]
        band5 = subset[3]
        band7 = subset[4]

        b_logs = numexpr.evaluate("log(subset)")
        logb2 = b_logs[0]
        logb3 = b_logs[1]
        logb4 = b_logs[2]
        logb5 = b_logs[3]
        logb7 = b_logs[4]

        b2b3  = numexpr.evaluate("band2 * band3")
        b2b4  = numexpr.evaluate("band2 * band4")
        b2b5  = numexpr.evaluate("band2 * band5")
        b2b7  = numexpr.evaluate("band2 * band7")
        b2lb2 = numexpr.evaluate("band2 * logb2")
        b2lb3 = numexpr.evaluate("band2 * logb3")
        b2lb4 = numexpr.evaluate("band2 * logb4")
        b2lb5 = numexpr.evaluate("band2 * logb5")
        b2lb7 = numexpr.evaluate("band2 * logb7")

        b3b4  = numexpr.evaluate("band3 * band4")
        b3b5  = numexpr.evaluate("band3 * band5")
        b3b7  = numexpr.evaluate("band3 * band7")
        b3lb2 = numexpr.evaluate("band3 * logb2")
        b3lb3 = numexpr.evaluate("band3 * logb3")
        b3lb4 = numexpr.evaluate("band3 * logb4")
        b3lb5 = numexpr.evaluate("band3 * logb5")
        b3lb7 = numexpr.evaluate("band3 * logb7")

        b4b5  = numexpr.evaluate("band4 * band5")
        b4b7  = numexpr.evaluate("band4 * band7")
        b4lb2 = numexpr.evaluate("band4 * logb2")
        b4lb3 = numexpr.evaluate("band4 * logb3")
        b4lb4 = numexpr.evaluate("band4 * logb4")
        b4lb5 = numexpr.evaluate("band4 * logb5")
        b4lb7 = numexpr.evaluate("band4 * logb7")

        b5b7  = numexpr.evaluate("band5 * band7")
        b5lb2 = numexpr.evaluate("band5 * logb2")
        b5lb3 = numexpr.evaluate("band5 * logb3")
        b5lb4 = numexpr.evaluate("band5 * logb4")
        b5lb5 = numexpr.evaluate("band5 * logb5")
        b5lb7 = numexpr.evaluate("band5 * logb7")

        b7lb2 = numexpr.evaluate("band7 * logb2")
        b7lb3 = numexpr.evaluate("band7 * logb3")
        b7lb4 = numexpr.evaluate("band7 * logb4")
        b7lb5 = numexpr.evaluate("band7 * logb5")
        b7lb7 = numexpr.evaluate("band7 * logb7")

        lb2lb3 = numexpr.evaluate("logb2 * logb3")
        lb2lb4 = numexpr.evaluate("logb2 * logb4")
        lb2lb5 = numexpr.evaluate("logb2 * logb5")
        lb2lb7 = numexpr.evaluate("logb2 * logb7")

        lb3lb4 = numexpr.evaluate("logb3 * logb4")
        lb3lb5 = numexpr.evaluate("logb3 * logb5")
        lb3lb7 = numexpr.evaluate("logb3 * logb7")

        lb4lb5 = numexpr.evaluate("logb4 * logb5")
        lb4lb7 = numexpr.evaluate("logb4 * logb7")

        lb5lb7 = numexpr.evaluate("logb5 * logb7")

        band_ratio1 = numexpr.evaluate("(band4-band3) / (band4+band3)")
        band_ratio2 = numexpr.evaluate("(band4-band5) / (band4+band3)")
        band_ratio3 = numexpr.evaluate("(band5-band3) / (band5+band3)")

        # The 2009_08_10 and 2012_12_07 versions use a different interactive
        # terms array compared to the 2013_01_08 version
        # 2013_01_08 uses 59 endmebers
        # 2009_08_10 uses 56 endmebers
        # 2012_12_07 uses 56 endmebers
        # TODO write an interface that can retrieve the correct
        # interactiveTerms array according to the specified version.

        interactiveTerms=numpy.array([
                                      b2b3, b2b4, b2b5, b2b7, b2lb2, b2lb3, b2lb4,
                                      b2lb5, b2lb7, b3b4, b3b5, b3b7, b3lb2, b3lb3,
                                      b3lb4, b3lb5, b3lb7, b4b5, b4b7, b4lb2, b4lb3,
                                      b4lb4, b4lb5, b4lb7, b5b7, b5lb2, b5lb3, b5lb4,
                                      b5lb5, b5lb7, b7lb2, b7lb3, b7lb4, b7lb5, b7lb7,
                                      lb2lb3, lb2lb4, lb2lb5, lb2lb7, lb3lb4, lb3lb5,
                                      lb3lb7, lb4lb5, lb4lb7, lb5lb7, band2, band3,
                                      band4 ,band5, band7, logb2, logb3, logb4, logb5,
                                      logb7, band_ratio1, band_ratio2, band_ratio3
                                     ])

        # Now add the sum to one constraint to the interactive terms
        # First make a zero array of the right shape
        weightedSpectra = numpy.zeros((interactiveTerms.shape[0]+1,)+interactiveTerms.shape[1:])
        # Insert the interactive terms
        weightedSpectra[:-1, ...] = interactiveTerms
        # Last element is special weighting
        weightedSpectra[-1] = sumToOneWeight

        inNullValDN = 0.0001
        outUnmixNullVal = 0

        logging.info("Calling unmiximage() for tile")

        fractions[:,ystart:yend,xstart:xend] = unmiximage.unmiximage(weightedSpectra, endmembers_array, inNullValDN, outUnmixNullVal) # The fortan method

    # 2013v gives green, dead1, dead2 and bare fractions
    # the last band should be the unmixing error
    logging.info("returning fractions")
    return fractions






# def linefinder(array, string = ""):
#     '''
#     Searches a list for the specified string.
#
#     :param array: A list containing searchable strings.
#
#     :param string: User input containing the string to search.
#
#     :return: The line containing the found sting.
#     '''
#
#     for line in array:
#         if string in str(line):
#             return line





# def locate(pattern, root):
#     '''
#     Finds files that match the given pattern. This will not search sub-directories.
#
#     :param pattern: A string containing the pattern to search, eg '*.csv'
#
#     :param root: The path directory to search
#
#     :return: A list of file-path name strings of files that match the given pattern.
#     '''
#
#     matches = []
#     for fl in os.listdir(root):
#         if fnmatch.fnmatch(fl, pattern):
#             matches.append(os.path.join(root, fl))
#
#     return matches

