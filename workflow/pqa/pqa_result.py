import logging
import numpy


class PQAResult(object):
    '''
    Represents the PQA result
    '''

    def __init__(self, shape, dtype=numpy.uint16):
        '''
        Constructor
        '''
        assert shape is not None 

        self.test_set = set()
        self.array = numpy.zeros(shape, dtype=dtype)
        self.bitcount = self.array.itemsize * 8

    def set_mask(self, mask, bit_index, unset_bits=False):
        '''Takes a boolean mask array and sets the bit in the result array.
        '''
        assert mask.shape == self.array.shape, \
            "Mask shape %s does not match result array %s" % (mask.shape, self.array.shape)
        assert mask.dtype == numpy.bool, 'Mask must be of type bool'
        assert 0 <= bit_index < self.bitcount, 'Invalid bit index'
        assert bit_index not in self.test_set, 'Bit %d already set' % bit_index
        self.test_set.add(bit_index)

        logging.debug('Setting result for bit %d', bit_index)
        numpy.bitwise_or(self.array, (mask << bit_index), self.array) # Set any 1 bits
        if unset_bits:
            numpy.bitwise_and(self.array, ~(~mask << bit_index), self.array) # Clear any 0 bits

    def get_mask(self, bit_index):
        """
        Return boolean mask for specified bit index
        """
        assert 0 <= bit_index < self.bitcount, 'Invalid bit index'
        assert bit_index in self.test_set, 'Test %d not run' % bit_index
        return (self.array & (1 << bit_index)) > 0

    @property
    def test_list(self):
        '''Returns a sorted list of all bit indices which have been set
        '''
        return sorted(self.test_set)

    @property
    def test_string(self):
        '''Returns a string showing all bit indices which have been set
        '''
        bit_list = ['0'] * self.bitcount
        for test_index in self.test_set:
            bit_list[test_index] = '1' # Show bits as big-endian
            # bit_list[15 - test_index] = '1' # Show bits as little-endian

        return ''.join(bit_list)
