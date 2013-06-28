#!/usr/bin/env python

import os, pprint, argparse, errno, re, time, fnmatch, datetime
import ConfigParser
import subprocess
import logging.config

ppr = pprint.pprint

# Absolute path for job runner and its log.

#ROOT = '/projects/v10/dev/rae/scale-test/jobs'
ROOT = '/home/547/axi547/ULA3'
EXE = os.path.normpath(os.path.join(ROOT, __file__))

# Set up logging

logging.config.fileConfig(os.path.join(ROOT, 'logging.config'))
ROOT_LOGGER_NAME = os.path.split(__file__)[1].strip('.py')
log = logging.getLogger(ROOT_LOGGER_NAME)

# Globals

JOB_NAME = '_ula3'
CONFIG_DICT = {}
STATUS_FILE_OK = 'STATUS.OK'
STATUS_FILE_ERROR = 'STATUS.ERROR'


# Template for dummy job to check template substitution and submission logic.

TEST_SCRIPT_TEMPLATE = """
#!/bin/bash
#PBS -P v10
#PBS -q normal
#PBS -l walltime=00:05:00
#PBS -wd

module load python/2.7.1
module load geotiff/1.3.0
module load hdf4/4.2.6_2012
module load proj/4.7.0
module load gdal/1.9.0

echo "-------- TEST_SCRIPT_TEMPLATE --------"
echo "__dataset_path       = %(__dataset_path)s"
echo "__product_root       = %(__product_root)s"
echo "__tmp_root           = %(__tmp_root)s"
echo "__status_file_ok     = %(__status_file_ok)s"
echo "__status_file_error  = %(__status_file_error)s"
echo "__link_job           = %(__link_job)s"

/bin/true
test $? -eq 0 && touch %(__status_file_ok)s || touch %(__status_file_error)s

sleep 120
%(__link_job)s
"""


# Template for NBAR production job.
# Set walltime=90min to prevent job links from breaking when things get slow.
# If the system is really and truly hosed, the job will time out anyway...

NBAR_SCRIPT_TEMPLATE = """
#!/bin/bash
#PBS -P v10
#PBS -q normal
#@#PBS -l walltime=00:45:00,vmem=8778MB,ncpus=1,jobfs=16GB
#PBS -l walltime=00:60:00,vmem=8778MB,ncpus=1
#PBS -wd
#@#PBS -m e
#@#PBS -M %(__notification_email)s

module load python/2.7.1
module load geotiff/1.3.0
module load hdf4/4.2.6_2012
module load proj/4.7.0
module load gdal/1.9.0

export ULA_NBAR_ROOT=%(__nbar_root)s
export PATH=$ULA_NBAR_ROOT:$ULA_NBAR_ROOT/bin:${PATH}
PYTHONPATH=$ULA_NBAR_ROOT:$ULA_NBAR_ROOT/common:$PYTHONPATH
PYTHONPATH=/projects/v10/GASoftware/pyephem-3.7.5.1/lib/python2.7/site-packages:$PYTHONPATH

# Only process dataset if path is not an empty string
if [ "%(__dataset_path)s" != "" ]
then
  # Use explicit work directory
  %(__nbar_exe)s --input=%(__dataset_path)s --output-root=%(__product_root)s --work=%(__tmp_root)s/$PBS_JOBID --repackage

  # Use $PBS_JOBFS work directory
  #@%(__nbar_exe)s --input=%(__dataset_path)s --output-root=%(__product_root)s --work=$PBS_JOBFS --repackage

  # Put status marker file in the log directory
  test $? -eq 0 && touch %(__status_file_ok)s || touch %(__status_file_error)s
fi

# Submit the next job
sleep 2
%(__link_job)s
"""



def __subdirs(path):

    try:
        return sorted(os.listdir(path))
    except OSError, e:
        if e.errno != errno.ENOENT:
            raise
    return []


def products():

    try:
        return __subdirs(CONFIG_DICT['product_root'])
    except KeyError:
        return []


def logged():

    try:
        return __subdirs(CONFIG_DICT['log_root'])
    except KeyError:
        return []


def unprocessed():

    try:
        paths = CONFIG_DICT['input_paths']
    except KeyError:
        return []

    names = [os.path.split(x)[1] for x in paths]
    # TODO assert inputs actually exist
    return [paths[names.index(n)] for n in
            [x for x in names if x not in logged()]]


def next():

    try:
        return unprocessed().pop(0)
    except IndexError:
        return None


def create_dir(path):

    log = logging.getLogger(ROOT_LOGGER_NAME + '.create_dir')

    try:
        os.makedirs(path)
        log.debug('path=%s' % path)
        return True
    except OSError, e:
        if e.errno == errno.EEXIST:
            return False
        else:
            log.error(str(e))
            raise


def get_inputs(root, fpattern='*_L1?.xml'):
    def filter(filelist, include_pattern='.*', exclude_pattern='.*_SYS_.*'):
        return [file for file in filelist if re.search(include_pattern, file) and not re.search(exclude_pattern, file)]

    log = logging.getLogger(ROOT_LOGGER_NAME + '.get_inputs')
    log.debug('walking %s (pattern=%s)...' % (root, fpattern))
    return filter([re.sub('/scene01$', '', r) for r,d,f in os.walk(root) if fnmatch.filter(f, fpattern)])


def submit_job(dataset_path=None, config_path=None, link=False, test=False):

    log = logging.getLogger(ROOT_LOGGER_NAME + '.submit_job')
    log.debug('%s link=%s config_path=%s'
              % (dataset_path, link, config_path))

    if config_path is None:
        log.error('config_path is undefined')
        return None

    if dataset_path is None:
        if link: # If no more datasets to process, then move onto next config file (if defined)
            next_config = CONFIG_DICT['next_config']
            if next_config:
                config_path = next_config # Move onto next config file
                dataset_path = '' # Set path to empty string to bypass processing stage

    if dataset_path is None:
        log.error('dataset_path is undefined')
        return None

    # Ensure product and log root directories exist.

    create_dir(CONFIG_DICT['product_root'])
    create_dir(CONFIG_DICT['log_root'])
    create_dir(CONFIG_DICT['tmp_root'])

    # Create the job log directory.

    if dataset_path: # Skip this for empty string
        dataset_name = os.path.split(dataset_path)[1]
        log_dir = os.path.join(CONFIG_DICT['log_root'], dataset_name)

        # Recursively fetch next scene if another job has already started this one
        if create_dir(log_dir) is False:
            return submit_job(next(), config_path, link, test)
    else:
        dataset_name = ''
        log_dir = CONFIG_DICT['log_root']

    # Define locals required to populate the job template.
    # TODO Debugging errors in template substitution is ugly!!
    # TODO Add a sanity check: template keys vs locals

    if link:
        # DEBUG ...allow --link to spawn one job only
        #__link_job = '%s --jobs=1 --config=%s' % (EXE, config_path)

        __link_job = '%s --jobs=1 --config=%s --link' % (EXE, config_path)
    else:
        __link_job = '# link=False'

    __dataset_path = dataset_path
    __product_root = CONFIG_DICT['product_root']
    __tmp_root = CONFIG_DICT['tmp_root']
    __nbar_exe = CONFIG_DICT['nbar_exe']
    __notification_email = CONFIG_DICT['notification_email']
    __nbar_root = CONFIG_DICT['nbar_root']
    __status_file_ok = STATUS_FILE_OK
    __status_file_error = STATUS_FILE_ERROR

    # Write the job script into the log directory.

    scr_file = os.path.join(log_dir, JOB_NAME)

    if test:
        scr_text = TEST_SCRIPT_TEMPLATE.lstrip('\n') % locals()
    else:
        scr_text = NBAR_SCRIPT_TEMPLATE.lstrip('\n') % locals()

    with open(scr_file, 'w') as _file:
        _file.write(scr_text)
        log.debug('created job script: %s' % scr_file)

    # Check the state of the input dataset (online, offline, ...)
    if dataset_path: # Skip this for empty string
        p = subprocess.Popen('dmls -l %s' % dataset_path, shell=True,
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        pout, perr = p.communicate()
        log.debug('input dataset status\n%s' % (pout + perr))

    # Submit the job.

    p = subprocess.Popen('qsub %s' % JOB_NAME, shell=True, cwd=log_dir,
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    pout, perr = p.communicate()
    job_id = None
    if p.returncode == 0 and re.search('^(\d+).(\w+)', pout):
        job_id = pout.strip('\n')
        log.info('%s %s (link=%s)' % (job_id, dataset_path, link))
    else:
        log.error('SUBMIT FAILED\n\n%s' % '\n'.join([pout, perr, '\n']))

    return job_id


def report(log_root):

    log = logging.getLogger(ROOT_LOGGER_NAME + '.report')

    if not os.path.exists(log_root):
        log.warning('log root not found: %s' % log_root)
        return

    # Generate a report for the run.

    job_status = { 'pass': [], 'fail': [], 'unknown': [] }
    job_process = { 'process': [], 'package': [] }
    process_etimes = []
    package_etimes = []

    # Count passed, failed, unknown jobs and gather elapsed times
    # for jobs that completed successfully.

    for log in sorted(os.listdir(log_root)):
        d = os.path.join(log_root, log)
        if os.path.isdir(d):
            files = sorted(os.listdir(d))

            state = 'unknown'
            if STATUS_FILE_OK in files:
                state = 'pass'
            elif STATUS_FILE_ERROR in files:
                state = 'fail'

            job_status[state].append(d)

            if state == 'pass':
                out_file = os.path.join(d, fnmatch.filter(files, '%s.o*' % JOB_NAME)[0])

                with open(out_file, 'r') as _file:
                    out_text = _file.read()
                if re.search('process completed successfully', out_text):
                    job_process['process'].append(d)
                    m = re.search('Elapsed time:(\s+)(\d+):(\d+):(\d+)', out_text)
                    if m:
                        _etime = datetime.timedelta(hours=int(m.group(2)),
                                                       minutes=int(m.group(3)),
                                                       seconds=int(m.group(4)))
                        process_etimes.append(_etime)
                if re.search('Repackaging completed successfully', out_text):
                    job_process['package'].append(d)
                    m = re.search('Elapsed time:(\s+)(\d+):(\d+):(\d+)', out_text)
                    if m:
                        _etime = datetime.timedelta(hours=int(m.group(2)),
                                                       minutes=int(m.group(3)),
                                                       seconds=int(m.group(4)))
                        package_etimes.append(_etime)

    # Calculate mean and max elapsed times.

    if process_etimes:
        process_time_min = min(process_etimes)
        process_time_max = max(process_etimes)

        process_seconds_sum = 0.0
        for t in process_etimes:
            process_seconds_sum += t.total_seconds()
        process_seconds_mean = process_seconds_sum / len(process_etimes)
    else:
        process_time_min = datetime.timedelta(seconds=0)
        process_time_max = datetime.timedelta(seconds=0)
        process_seconds_mean = 0

    if package_etimes:
        package_time_min = min(package_etimes)
        package_time_max = max(package_etimes)

        package_seconds_sum = 0.0
        for t in package_etimes:
            package_seconds_sum += t.total_seconds()
        package_seconds_mean = package_seconds_sum / len(package_etimes)
    else:
        package_time_min = datetime.timedelta(seconds=0)
        package_time_max = datetime.timedelta(seconds=0)
        package_seconds_mean = 0

    # Print report.

    _stats_txt = (
        '\nSUMMARY [%s]\n'
        '\tpass    %d\t(%d processed, %d packaged)\n'
        '\tfail    %d\n'
        '\tunknown %d\n'
        '\tmin / mean / max elapsed processing times: %.1f / %.1f / %.1f sec\t(%s / %s / %s)\n'
        '\tmin / mean / max elapsed packaging times:  %.1f / %.1f / %.1f sec\t(%s / %s / %s)\n'
    ) % (
        log_root,
        len(job_status['pass']), len(job_process['process']), len(job_process['package']),
        len(job_status['fail']),
        len(job_status['unknown']),
        process_time_min.total_seconds(), process_seconds_mean, process_time_max.total_seconds(),
        process_time_min, datetime.timedelta(seconds=process_seconds_mean), process_time_max,
        package_time_min.total_seconds(), package_seconds_mean, package_time_max.total_seconds(),
        package_time_min, datetime.timedelta(seconds=package_seconds_mean), package_time_max
    )

    print _stats_txt


def list_incomplete_logs(log_root):

    ilog_list = []
    for log in sorted(os.listdir(log_root)):
        logdir = os.path.join(log_root, log)
        logfiles = sorted(os.listdir(logdir))
        if (STATUS_FILE_OK not in logfiles and
            STATUS_FILE_ERROR not in logfiles and
            JOB_NAME in logfiles):
            for __r, __d, __f in os.walk(logdir, topdown=False):
                ilog_list.append(__r)

    print 'INCOMPLETE LOGS'
    for x in ilog_list:
        print x


def foo():

    # DEBUG for testing logging

    log = logging.getLogger(ROOT_LOGGER_NAME + '.foo')
    log.info('fooey!')


def parse_config(config_path):

    log = logging.getLogger(ROOT_LOGGER_NAME + '.parse_config')
    log.debug('config file: %s' % config_path)

    c = ConfigParser.ConfigParser(allow_no_value=True)
    c.optionxform = str
    c.readfp(open(config_path))

    _ipaths = []

    try:

        # Read [INPUTS] section if present.
        _ipaths = [x[0] for x in c.items('INPUTS')]
        log.debug('found [INPUTS] in config file')

    except ConfigParser.NoSectionError:

        # No [INPUTS], so gather input list...
        _iroot = c.get('CONFIGURATION', 'input_root')
        _ipaths = get_inputs(_iroot, fpattern='*_L1?.xml')

        # Add [INPUTS] section to config instance.
        c.add_section('INPUTS')
        for i in _ipaths:
            c.set('INPUTS', i, '')

        # Write modified config file.
        with open(config_path, 'w') as _file:
            c.write(_file)

        log.debug('created [INPUTS] in config file %s' % config_path)

    return {
        'product_root': c.get('CONFIGURATION', 'product_root'),
        'log_root':     c.get('CONFIGURATION', 'log_root'),
        'tmp_root':     c.get('CONFIGURATION', 'tmp_root'),
        'nbar_exe':     c.get('CONFIGURATION', 'nbar_exe'),
        'nbar_root':     c.get('CONFIGURATION', 'nbar_root'),
        'notification_email':     c.get('CONFIGURATION', 'notification_email'),
        'next_config':     c.get('CONFIGURATION', 'next_config'),
        'input_paths':  _ipaths,
        'input_names':  [os.path.split(x)[1] for x in _ipaths],
        'config_path':  config_path
    }


def status():

    print
    print 'LOGGED (%s)' % CONFIG_DICT['log_root']
    ppr(logged())
    print
    print 'PRODUCTS (%s)' % CONFIG_DICT['product_root']
    ppr(products())
    print
    print 'UNPROCESSED'
    ppr(unprocessed())
    print
    print 'NEXT'
    print next()
    print


def parse_args():

    argp = argparse.ArgumentParser(__file__)
    # TODO require --config??
    argp.add_argument('-c', '--config', default=None, required=True,
                      help='Specify configuration file')
    argp.add_argument('-s', '--status', action='store_true', default=False,
                      help='Display NBAR processing status')
    argp.add_argument('-j', '--jobs', type=int, default=0,
                      help='Specify number of jobs to submit')
    argp.add_argument('-l', '--link', action='store_true', default=False,
                      help='Enable linked jobs')
    argp.add_argument('-t', '--test', action='store_true', default=False,
                      help='Submit dummy jobs for testing')
    argp.add_argument('-r', '--report', action='store_true', default=False, help='Summarise completed jobs')
    argp.add_argument('--incomplete', action='store_true', default=False,
                      help='List log directories that are incomplete (active job or no result)')
    #argp.add_argument('--foo', action='store_true', default=False,
    #                  help='dummy flag (logging tests)')
    return argp.parse_args()


def main():
    global CONFIG_DICT

    log = logging.getLogger(ROOT_LOGGER_NAME + '.main')

    args = parse_args()
    CONFIG_DICT = parse_config(os.path.join(ROOT, args.config))

    if args.status:
        status()
    elif args.report:
        report(CONFIG_DICT['log_root'])
    elif args.incomplete:
        list_incomplete_logs(CONFIG_DICT['log_root'])
    else:
        # Submit PBS job(s)
        log.debug('jobs=%d link=%s' % (args.jobs, args.link))
        for j in xrange(args.jobs):
            flag = submit_job(dataset_path=next(), config_path=args.config,
                              link=args.link, test=args.test)
            if flag is None:
                log.info('ALL DATASETS PROCESSED')
                break
            time.sleep(2)



if __name__ == '__main__':
    main()
