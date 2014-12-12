import logging
import logging.config
import os

class MemuseFilter(logging.Filter):
    def filter(self, record):
        """ This function overrides logging.Filter, adds memuse as a field
        """ 
        record.memuse = self.str_mem()
        return True

    # Following code from http://stackoverflow.com/a/938800/819110: 
    _proc_status = '/proc/%d/status' % os.getpid()
    _scale = {'kB': 1024.0, 'mB': 1024.0*1024.0,
              'KB': 1024.0, 'MB': 1024.0*1024.0}

    def _VmB(self,VmKey):
        """Private.
        """
        # get pseudo file  /proc/<pid>/status
        try:
            t = open(self._proc_status)
            v = t.read()
            t.close()
        except:
            return 0.0  # non-Linux?
        # get VmKey line e.g. 'VmRSS:  9999  kB\n ...'
        i = v.index(VmKey)
        v = v[i:].split(None, 3)  # whitespace
        if len(v) < 3:
            return 0.0  # invalid format?
        # convert Vm value to bytes
        return float(v[1]) * self._scale[v[2]]

    def memory(self,since=0.0):
        """Return memory usage in bytes.
        """
        return self._VmB('VmSize:') - since

    def resident(self,since=0.0):
        """Resident memory in bytes.
        """
        return self._VmB('VmRSS:') - since

    def swapsize(self,since=0.0):
        """Return swap size in bytes.
        """
        return self._VmB('VmSwap:') - since

    def byte_to_mb(self,byte):
        """return size in MB (being lazy)
        """
        return byte/(1024*1024)

    def str_mem(self):
        """Return a string with the total memuse and swap size in MB
        """
        return "MemTotal:%.0fM,Resident:%.0fM,Swap:%.0fM"%(\
            self.byte_to_mb(self.memory()),self.byte_to_mb(self.resident()),self.byte_to_mb(self.swapsize()) )
