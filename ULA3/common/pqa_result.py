import threading, logging, numpy
from ULA3 import DataGrid

logger = logging.getLogger('root.' + __name__)

class PQAResult(DataGrid):
    '''
    DataGrid descendant with locking and tracking of PQA tests run
    '''

    def __init__(self, image_array=None, filename=None, dtype='uint16'):
        '''
        Constructor
        '''
        assert image_array is not None or filename is not None # N.B: Filename will have priority over image_array

        # Lock object to prevent multiple instances of brdf_dict
        self.result_lock = threading.Lock()

        self.bitcount = 0
        self.test_set = set()

        logger.debug('Acquiring PQAResult lock')
        self.result_lock.acquire()
        try:
            if filename: # Restore result array from existing file
                DataGrid.__init__(self, filename=filename)
                self.bitcount = 8 * self.array.itemsize

                # This is a bit nasty - need to reconstitute self.test_set from array
                # Operation will fail if all pixels have failed any one test, but should only happen when debugging
                for bit_index in range(self.bitcount):
                    if (self.array & (1 << bit_index)).any():
                        self.test_set.add(bit_index)
            else: # Create empty result array
                DataGrid.__init__(self, array=numpy.zeros_like(image_array[0]).astype(numpy.dtype(dtype)))
                self.bitcount = 8 * self.array.itemsize
        finally:
            logger.debug('Releasing PQAResult lock')
            self.result_lock.release()

    def get_lock(self):
        logger.debug('Acquiring PQAResult lock')
        self.result_lock.acquire()
        return self.array

    def unlock(self):
        logger.debug('Releasing PQAResult lock')
        self.result_lock.release()

    def set_mask(self, mask, bit_index, unset_bits=False):
        '''Takes a boolean mask array and sets the bit in the result array.
        '''
        assert mask.shape == self.array.shape, 'Mismatched shape between mask and result arrays'
        assert mask.dtype == numpy.bool, 'Mask must be of type bool'
        assert 0 <= bit_index < self.bitcount, 'Invalid bit index'
        assert bit_index not in self.test_set, 'Bit %d already set' % bit_index
        self.test_set.add(bit_index)

        self.result_lock.acquire()
        logger.debug('Acquired PQAResult lock for bit %d', bit_index)
        try:
            logger.debug('Setting result for bit %d', bit_index)
#            self.array += (mask << bit_index) # Can't do addition more than once
            numpy.bitwise_or(self.array, (mask << bit_index), self.array) # Set any 1 bits
            if unset_bits:
                numpy.bitwise_and(self.array, ~(~mask << bit_index), self.array) # Clear any 0 bits
        finally:
            logger.debug('Releasing PQAResult lock for bit %d', bit_index)
            self.result_lock.release()

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
