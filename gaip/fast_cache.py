import os
import subprocess

def fast_read(origin_path, ssd_env_var='PBS_JOBFS', cache_scope=os.getpid()):

    """
    Provide fast, read-only access to the specified a file or directory
    (origin_path) by first copying data to solid state drive defined
    by the supplied enviroment variable PBS_JOBFS

    Copy contents of origin_path to temporary directory on SSD
    and return the path to the ssd_path. Copy is only performed once
    and subsequent calls simply return the ssd_path.

    If no SSD device is available (PBS_JOBSFS not defined) then
    no copy is performed and ssd_path == origin_path

    The cached copy of the data is shared within the specified scope which 
    defaults to the current process. Use:

        path = fast_read(some_path, cache_scope=os.environ['PBS_JOBID'])

    for a cache which operates across the current PBS job.
    """

    # origin must exist, if not then NO-OP and let the caller handle it
    if not os.path.exists(origin_path):
        return origin_path

    # PBS_JOBFS must be defined, otherwise NO-OP as we have no SSD
    if ssd_env_var not in os.environ:
        return origin_path

    # compute the path to the root of the cache on SSD
    cache_root = os.path.join(os.environ[ssd_env_var], 'fast_read_cache', \
        str(cache_scope))
    if not os.path.exists(cache_root):
        os.makedirs(cache_root)

    # calcuate the ssd_path
    ssd_path = "%s/%s" % (cache_root, origin_path)
    ssd_path = os.path.abspath(ssd_path)

    # return it if already in cache
    if os.path.exists(ssd_path):
        return ssd_path

    # not there so copy it into the cache

    # first make the target directory
    target_dir = os.path.dirname(ssd_path)

    if not os.path.exists(target_dir):
        os.makedirs(target_dir)

    # now copy the data to the cache
    recurse = ''
    if os.path.isdir(origin_path):
        recurse = '-r'

    cmd = "cp %s %s %s" % (recurse, origin_path, target_dir)

    # execute the copy in a child process
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    for line in p.stdout.readlines():
        pass
    retval = p.wait()

    # If we can't copy for ANY reason then use original file/dir in-situ 
    if retval != 0:
        return origin_path

    # got this far, then return the path to the SSD cache copy

    return ssd_path
