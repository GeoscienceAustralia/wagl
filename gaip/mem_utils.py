#!/bin/env python

"""
Memory tracking functions
Following code from http://stackoverflow.com/a/938800/819110:a
"""

import os

SCALE = {'kB': 1024.0, 'mB': 1024.0*1024.0, \
    'KB': 1024.0, 'MB': 1024.0*1024.0}
PROC_MEMINFO = "/proc/meminfo"

def get_proc_value(proc_key, fname=None):
    """
    Extract the value from /proc/<pid>/status that
    matches the supplied proc_key
    """
    if fname is None:
        fname = '/proc/%d/status' % os.getpid()
    with open(fname) as t:
        v = t.read()

    # get pooc_key line e.g. 'VmRSS:  9999  kB\n ...'
    i = v.index(proc_key)
    v = v[i:].split(None, 3)  # whitespace
    if len(v) < 3:
        return 0.0  # invalid format?
    # convert Vm value to bytes
    return float(v[1]) * SCALE[v[2]]

def get_mem_total_mb(since=0.0):
    """
    retrun total physical memory on this host in MB
    """
    return get_proc_value('MemTotal:', fname=PROC_MEMINFO) / SCALE['MB'] - since

def get_vm_bytes(since=0.0):
    """
    Return memory usage in bytes.
    """
    return int(get_proc_value('VmSize:') - since)


def get_rss_bytes(since=0.0):
    """
    Resident memory in bytes.
    """
    return int(get_proc_value('VmRSS:') - since)


def get_rss_mbytes(since=0.0):
    """
    Resident memory in bytes.
    """
    return get_proc_value('VmRSS:') / SCALE['MB'] - since


def get_swap_bytes(since=0.0):
    """
    Return swap size in bytes.
    """
    return int(get_proc_value('VmSwap:') - since)


if __name__=='__main__':
    print get_rss_mbytes()
    print get_mem_total_mb()
