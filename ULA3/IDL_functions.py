import numpy
import datetime
import _idl_histogram

def IDL_Histogram(data, binsize=None, max=None, min=None, nbins=None, omax=None, omin=None, reverse_indices=None, locations=None, input=None, NaN=False):
    """
    Replicates the histogram function avaiable within IDL (Interactive Data Language, EXELISvis).

    :param data:
        A 1-Dimensional array to calculate the histogram for.

    :param binsize:
        (Optional) The binsize (Default is 1) to be used for creating the histogram.

    :param max:
        (Optional) The maximum value to be used in creating the histogram.
        If not specified the array will be searched for max.

    :param min:
        (Optional) The minimum value to be used in creating the histogram.
        If not specified the array will be searched for min.

    :param nbins:
        (Optional) The number of bins to be used for creating the histogram.
        If set binsize is calculated as ((max - min) / (nbins - 1)), and the
        max value will be adjusted to (nbins*binsize + min).

    :param omax:
        (Optional) A string name used to refer to the dictionary key
        that will contain the maximum value used in generating the histogram.

    :param omin:
        (Optional) A string name used to refer to the dictionary key
        that will contain the minimum value used in generating the histogram.

    :param reverse_indices:
        (Optional) A string name used to refer to the
        dictionary key that will contain the reverse indices of the histogram.

    :param locations:
        (Optional) A string name used to refer to the dictionary
        key that will contain the starting locations of each bin.

    :param input:
        (Optional) Used to specify an input array that will be added to the
        result of the histogram. Useful for tiling mechanisms that only handle
        portions of an array at a time. The input array must be 1-Dimensional
        and contain at least as many elements as are required to construct
        the histogram.

    :param NaN:
        If set to True (Default is False) then NaN values will be
        accounted for and treated as missing data.

    :return:
        A dictionary containing the histogram and other optional components.
        The dictionary key name for the histogram is 'histogram'.

    Example:

        >>> h = IDL_Histogram(data, min=0, max=max, omin='omin', omax='omax', reverse_indices='ri')
        >>> hist = h['histogram']
        >>> ri = h['ri']
        >>> loc = loc['ri']
        >>> data_at_ith_bin_indices = data[ri[ri[i]:ri[i+1]]]

    :author:
        Josh Sixsmith, joshua.sixsmith@ga.gov.au

    :history:
       * 04/02/2013: Created
       *  05/04/2013: Added NaN keyword
       *  05/06/2013: Now checks for max value of 256 and datatype of 'uint8'
       *  12/06/2013: Added input keyword

    """
    def hist_int(data, n, min, max, binsize, nbins, max_bin, ri):
        # increase the size by one. When specifying a min and max, it shouldn't
        # be included in the histogram. Stuff not to be included gets dumped
        # into the 1st position then removed prior to returning to the user.

        nbins_ = nbins + 1
        hist = numpy.zeros(nbins_, dtype='uint32')

        _idl_histogram.idl_histogram.histogram_int(data, hist, n, nbins_, min, max, max_bin, binsize)

        if ri:
            return hist
        else:
            return hist[1:]

    def hist_long(data, n, min, max, binsize, nbins, max_bin, ri):
        # increase the size by one. When specifying a min and max, it shouldn't
        # be included in the histogram. Stuff not to be included gets dumped
        # into the 1st position then removed prior to returning to the user.

        nbins_ = nbins + 1
        hist = numpy.zeros(nbins_, dtype='uint32')

        _idl_histogram.idl_histogram.histogram_long(data, hist, n, nbins_, min, max, max_bin, binsize)

        if ri:
            return hist
        else:
            return hist[1:]

    def hist_dlong(data, n, min, max, binsize, nbins, max_bin, ri):
        # increase the size by one. When specifying a min and max, it shouldn't
        # be included in the histogram. Stuff not to be included gets dumped
        # into the 1st position then removed prior to returning to the user.

        nbins_ = nbins + 1
        hist = numpy.zeros(nbins_, dtype='uint32')

        _idl_histogram.idl_histogram.histogram_dlong(data, hist, n, nbins_, min, max, max_bin, binsize)

        if ri:
            return hist
        else:
            return hist[1:]

    def hist_float(data, n, min, max, binsize, nbins, max_bin, ri):
        # increase the size by one. When specifying a min and max, it shouldn't
        # be included in the histogram. Stuff not to be included gets dumped
        # into the 1st position then removed prior to returning to the user.

        nbins_ = nbins + 1
        hist = numpy.zeros(nbins_, dtype='uint32')

        _idl_histogram.idl_histogram.histogram_float(data, hist, n, nbins_, min, max, max_bin, binsize)

        if ri:
            return hist
        else:
            return hist[1:]

    def hist_dfloat(data, n, min, max, binsize, nbins, max_bin, ri):
        # increase the size by one. When specifying a min and max, it shouldn't
        # be included in the histogram. Stuff not to be included gets dumped
        # into the 1st position then removed prior to returning to the user.

        nbins_ = nbins + 1
        hist = numpy.zeros(nbins_, dtype='uint32')

        _idl_histogram.idl_histogram.histogram_dfloat(data, hist, n, nbins_, min, max, max_bin, binsize)

        if ri:
            return hist
        else:
            return hist[1:]

    def ri_int(data, hist, nbins, n, ri_sz, min, max, max_bin, binsize):
        # increase the size by one. When specifying a min and max, it shouldn't
        # be included in the reverse indices. Stuff not to be included gets
        # dumped into the 1st position then removed prior to returning to the
        # user.

        nbins_ = nbins + 1
        ri = numpy.zeros(ri_sz, dtype='uint32')

        _idl_histogram.idl_histogram.reverse_indices_int(data, hist, ri, nbins_, n, ri_sz, min, max, max_bin, binsize)

        return (hist[1:], ri[1:])

    def ri_long(data, hist, nbins, n, ri_sz, min, max, max_bin, binsize):
        # increase the size by one. When specifying a min and max, it shouldn't
        # be included in the reverse indices. Stuff not to be included gets
        # dumped into the 1st position then removed prior to returning to the
        # user.

        nbins_ = nbins + 1
        ri = numpy.zeros(ri_sz, dtype='uint32')

        _idl_histogram.idl_histogram.reverse_indices_long(data, hist, ri, nbins_, n, ri_sz, min, max, max_bin, binsize)

        return (hist[1:], ri[1:])

    def ri_dlong(data, hist, nbins, n, ri_sz, min, max, max_bin, binsize):
        # increase the size by one. When specifying a min and max, it shouldn't
        # be included in the reverse indices. Stuff not to be included gets
        # dumped into the 1st position then removed prior to returning to the
        # user.

        nbins_ = nbins + 1
        ri = numpy.zeros(ri_sz, dtype='uint32')

        _idl_histogram.idl_histogram.reverse_indices_dlong(data, hist, ri, nbins_, n, ri_sz, min, max, max_bin, binsize)

        return (hist[1:], ri[1:])

    def ri_float(data, hist, nbins, n, ri_sz, min, max, max_bin, binsize):
        # increase the size by one. When specifying a min and max, it shouldn't
        # be included in the reverse indices. Stuff not to be included gets
        # dumped into the 1st position then removed prior to returning to the
        # user.

        nbins_ = nbins + 1
        ri = numpy.zeros(ri_sz, dtype='uint32')

        _idl_histogram.idl_histogram.reverse_indices_float(data, hist, ri, nbins_, n, ri_sz, min, max, max_bin, binsize)

        return (hist[1:], ri[1:])

    def ri_dfloat(data, hist, nbins, n, ri_sz, min, max, max_bin, binsize):
        # increase the size by one. When specifying a min and max, it shouldn't
        # be included in the reverse indices. Stuff not to be included gets
        # dumped into the 1st position then removed prior to returning to the
        # user.

        nbins_ = nbins + 1
        ri = numpy.zeros(ri_sz, dtype='uint32')

        _idl_histogram.idl_histogram.reverse_indices_dfloat(data, hist, ri, nbins_, n, ri_sz, min, max, max_bin, binsize)

        return (hist[1:], ri[1:])


    def datatype(val):
        instr = str(val)
        return {
            'int8' : '1',
            'uint8' : '1',
            'int16' : '2',
            'uint16' : '12',
            'int32' : '3',
            'uint32' : '13',
            'int64' : '13',
            'uint64' : '15',
            'int' : '13',
            'float32' : '4',
            'float64' : '5',
            }.get(instr, 'Error')

    def data_convert(val, b):
        instr = str(val)
        return {
            'int8' : numpy.int8(b),
            'uint8' : numpy.uint8(b),
            'int16' : numpy.int16(b),
            'uint16' : numpy.uint16(b),
            'int32' : numpy.int32(b),
            'uint32' : numpy.uint32(b),
            'int64' : numpy.int64(b),
            'uint64' : numpy.uint64(b),
            'int' : numpy.int64(b),
            'float32' : numpy.float32(b),
            'float64' : numpy.float64(b),
            }.get(instr, 'Error')


    dtype = datatype(data.dtype.name)
    if (dtype == 'Error'):
       raise Exception('Error. Incompatable Data Type. Compatable Data Types Include: int8, uint8, int16, uint16, int32, uint32, int64, uint64, float32, float64')

    if len(data.shape) != 1:
        raise Exception('Error. Array must be 1 dimensional. Use .flatten()')

    if ((max != None) & (binsize != None) & (nbins != None)):
        raise Exception('Error. Conflicting Keywords. Max cannot be set when both binsize and nbins are set.')

    if ((input != None) & (reverse_indices != None)):
        raise Exception('Error. Conflicting Keywords. Both input and reverse_indices cannot be set at the same time.')

    if (max == None):
        if NaN:
            max = numpy.nanmax(data)
        else:
            max = numpy.max(data)

    if (min == None):
        if NaN:
            min = numpy.nanmin(data)
        else:
            min = numpy.min(data)

    min = data_convert(data.dtype.name, min)
    max = data_convert(data.dtype.name, max)

    if (binsize == None) & (nbins == None):
        #print 'binsize=None & nbins=None'
        binsize = 1
        nbins = (max - min) + 1
    elif (binsize == None):
        #print 'binsize=None & nbins= ', nbins
        binsize = (max - min) / (nbins - 1)
        max = nbins * binsize + min
    elif (binsize != None) & (nbins == None):
        #print 'binsize= ', binsize, 'nbins=None'
        nbins = numpy.floor((max - min) / binsize) + 1
    else:
        #print 'binsize= ', binsize, 'nbins= ', nbins
        max = nbins * binsize + min

    binsize = data_convert(data.dtype.name, binsize)
    min = data_convert(data.dtype.name, min)

    # If nbins is set to 256 and the array datatype is uint8, then the max
    # value will be adjusted to 256, however due to datatype conversions, the
    # max value of 256 will change to 0
    # This fix conforms with IDL.
    if ((max == 256) & (data.dtype.name == 'uint8')):
        max = 255
    max = data_convert(data.dtype.name, max)

    #probably also need to pass in a max binvalue into the fortran code
    # the max bin value is non-inclusive, but also check that the data
    #values are <= the max value
    # eg max value = 1.0, but max bin = 1.08, therefore a value of 1.04
    # will not be included
    max_bin = nbins * binsize + min

    if (binsize == 0):
        raise Exception("Error. Binsize = 0, histogram can't be computed.")

    # Probably unessessary to include the max and max_bin equality warning
    #if (max == max_bin):
    #    print "\n!!!!!Warning!!!!! \nMax is equal to the last bin's right edge, maximum value will not be included in the histogram."

    if (input != None):
        # Check that input is 1-Dimensional
        if (len(input.shape) != 1):
            raise Exception("Error. Input array must be 1 dimensional. Use array.flatten()")
        # Check that input is at least nbins in length
        if (input.shape[0] < nbins):
            print 'Number of elements of input: ', input.shape[0]
            print 'Minimum number of elemets required: ', nbins
            raise Exception("Error. Input array does not have enough elements.")

    n = numpy.size(data)

    # Some unsigned data types will be promoted as Fortran doesn't handle
    # unsigned data types.
    get_hist = {
                 'int8' : hist_int,
                 'uint8' : hist_int,
                 'int16' : hist_int,
                 'uint16' : hist_long,
                 'int32' : hist_long,
                 'uint32' : hist_dlong,
                 'int64' : hist_dlong,
                 'uint64' : hist_dlong,
                 'int' : hist_dlong,
                 'float32' : hist_float,
                 'float64' : hist_dfloat,
                }

    ri = False

    if (type(reverse_indices) == str):
        ri = True
        #print 'Creating histogram'
        #st = datetime.datetime.now()
        hist = get_hist[data.dtype.name](data, n, min, max, binsize, nbins, max_bin, ri)
        #et = datetime.datetime.now()
        #print 'histogram time taken: ', et - st
        cum_sum = numpy.sum(hist[1:])
        #print 'histogram', hist[1:]
        #st = datetime.datetime.now()
        ri_sz = nbins + cum_sum + 1+1

        get_ri = {
                   'int8' : ri_int,
                   'uint8' : ri_int,
                   'int16' : ri_int,
                   'uint16' : ri_int,
                   'int32' : ri_long,
                   'uint32' : ri_long,
                   'int64' : ri_dlong,
                   'uint64' : ri_dlong,
                   'int' : ri_dlong,
                   'float32' : ri_float,
                   'float64' : ri_dfloat,
                  }
        #print 'get reverse indices'
        hri = get_ri[data.dtype.name](data, hist, nbins, n, ri_sz, min, max, max_bin, binsize)

        results = {'histogram': hri[0]}
        results[reverse_indices] = hri[1]
        #et = datetime.datetime.now()
        #print 'reverse indices time taken: ', et - st
    else:
        #print 'Creating histogram'
        #st = datetime.datetime.now()
        hist = get_hist[data.dtype.name](data, n, min, max, binsize, nbins, max_bin, ri)
        #et = datetime.datetime.now()
        if (input != None):
            # Now to add the input array to the histogram.
            # The result will take the shape of the larger of the two arrays.
            if (input.shape[0] == hist.shape[0]):
               hist += input
               results = {'histogram': hist}
            else:
               temp = numpy.zeros(input.shape, dtype='uint32')
               temp[0:hist.shape[0]] = hist
               temp += input
               results = {'histogram': temp}
        else:
            results = {'histogram': hist}
        #print 'histogram time taken: ', et - st

    if (type(omax) == str):
        results[omax] = max

    if (type(omin) == str):
        results[omin] = min

    if (type(locations) == str):
        loc = numpy.zeros(nbins, dtype=data.dtype.name)
        for i in numpy.arange(nbins):
            loc[i] = min + i * binsize

        results[locations] = loc

    return results

