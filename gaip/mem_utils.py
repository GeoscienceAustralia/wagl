import os

# Following code from http://stackoverflow.com/a/938800/819110: 
_scale = {'kB': 1024.0, 'mB': 1024.0*1024.0,
          'KB': 1024.0, 'MB': 1024.0*1024.0}

def _VmB(VmKey):
    # get pseudo file  /proc/<pid>/status
    with open('/proc/%d/status' % os.getpid()) as t:
        v = t.read()
    
    # get VmKey line e.g. 'VmRSS:  9999  kB\n ...'
    i = v.index(VmKey)
    v = v[i:].split(None, 3)  # whitespace
    if len(v) < 3:
        return 0.0  # invalid format?
    # convert Vm value to bytes
    return float(v[1]) * _scale[v[2]]

def memory(since=0.0):
    """Return memory usage in bytes.
    """
    return _VmB('VmSize:') - since

def resident(since=0.0):
    """Resident memory in bytes.
    """
    return _VmB('VmRSS:') - since

def resident_mb(since=0.0):
    """Resident memory in bytes.
    """
    return _VmB('VmRSS:') / _scale['MB'] - since

def swapsize(since=0.0):
    """Return swap size in bytes.
    """
    return _VmB('VmSwap:') - since

def byte_to_mb(byte):
    """return size in MB (being lazy)
    """
    return byte/(1024*1024)
