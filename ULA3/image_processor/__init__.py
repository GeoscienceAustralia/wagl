'''
Created on Jun 14, 2012

@author: Alex Ip (alex.ip@ga.gov.au)

Implements functions used by workflow management framework
'''
#TODO: Implement generic subprocess list completion function to remove duplicated code
import os, sys, logging, re, threading, socket, traceback, argparse, ConfigParser
from datetime import datetime, date, time
from glob import glob
from os import path
from ULA3.meta import Singleton
from ULA3.utils import log_multiline, execute

# Set top level standard output
console_handler = logging.StreamHandler(sys.stdout)
console_handler.setLevel(logging.WARNING)
console_formatter = logging.Formatter('%(message)s')
console_handler.setFormatter(console_formatter)

# Set standard error output
#error_handler = logging.StreamHandler(sys.stderr)
#error_handler.setLevel(logging.WARNING)
#error_formatter = logging.Formatter('%(levelname)s: %(message)s')
#error_handler.setFormatter(error_formatter)

root_logger = logging.getLogger('root')
if not root_logger.level:
    root_logger.setLevel(logging.INFO) # Default logging level for all modules
    root_logger.addHandler(console_handler)
#    root_logger.addHandler(error_handler)

# N.B: Local logger must be set manually. Overridden by conf file value (if any)
logger = logging.getLogger('root.processor_config')
if not logger.level:
    logger.setLevel(logging.INFO)

thread_exception = None

def get_submodule_list(module_directory):
    """Function to collate a list of all (potentially) valid Python modules in current module directory
    Should be called by directory module's __init__.py

    Parameter:
        module_directory: Directory to search for sub-modules
    """
    submodule_list = [path.splitext(path.basename(_lib))[0] for _lib in glob(path.join(module_directory, '[a-z]*.py'))]
    # Append all valid package directories
    for pathname in [path.dirname(pathname) for pathname in glob(path.join(module_directory, '[a-z]*/__init__.py'))]:
        if path.isdir(pathname) and path.exists(path.join(pathname, '__init__.py')):
            submodule_list.append(path.basename(pathname))

    return sorted(submodule_list)


def flatten_process_list(process_list_item):
    """ Function to flatten a nested list
    """
    result = []
    if type(process_list_item) == str:
        return [process_list_item]
    elif type(process_list_item) == list:
        for item in process_list_item:
            result += flatten_process_list(item)
        return result


def get_process_list_from_text(list_text, level=0):
    """ Recursive function to create parse a text file and create a nested list
    with [[process1, process2...processn]] indicating that parallel processing is possible.

    Parameters:
        list_text: text to process
        level: integer indicating the level of recursion. DO NOT SET THIS WHEN CALLING
    """
    def replace_placeholder(_list, placeholder, sublist):
        """Helper function to replace a text placeholder with the previously extracted sub-list
        """
        logger.debug('  replace_placeholder(%s, %s, %s) called', repr(_list), repr(placeholder), repr(sublist))
        for index in range(len(_list)):
#            logger.debug('    _list[index] = %s', repr(_list[index]))
            if type(_list[index]) == list:
                replace_placeholder(_list[index], placeholder, sublist)
            elif _list[index] == placeholder:
                _list[index] = sublist
                break

    logger.debug('get_process_list_from_text(%s, %s) called', repr(list_text), repr(level))
    if list_text:
        _list_text = re.sub('#.*\n|\s|\"|\'', '', list_text) # Remove comments, whitespaces and quotes
        s = re.search('(.*\[)(.*?)(\].*)', _list_text) # Find innermost bracketed list (i.e. with no sublists)
        if s: # Bracketed list found
            if s.group(1) == '[' and s.group(3) == ']': # Flat bracketed list - don't process again
                sublist = s.group(2).split(',')
            else: # Nested list found - process again
                placeholder = '<sublist%s>' % (level)
                sublist = get_process_list_from_text(re.sub('\[$', '', s.group(1)) +
                                                     placeholder +
                                                     re.sub('^\]', '', s.group(3)), level + 1)
                replace_placeholder(sublist, placeholder, get_process_list_from_text(s.group(2)))
        else: # Flat unbracketed list
            sublist = _list_text.split(',')
    else: # No text - empty list
        sublist = []

    logger.debug('sublist = %s', repr(sublist))
    return sublist


def get_process_list_from_file(file_name):
    """Function to read a list of tasks from a given file
    """
    if file_name:
        infile = open(file_name, 'r')
        try:
            return get_process_list_from_text(infile.read())
        finally:
            infile.close()
    else:
        return []

def get_subprocess_list(module_name):
    """
    Function to collate a list of all tasks which the calling module should run by default if no single task is specified
    Should be called by directory module's __init__.py
    """
    CONFIG = ProcessorConfig()
    #return get_process_list_from_file(path.join(path.dirname(path.dirname(__file__)), module_name.replace('.', path.sep), '__init__.txt'))
    return get_process_list_from_file(path.join(CONFIG.code_root, module_name.replace('.', path.sep), '__init__.txt'))

def process(module_name, submodule_dict, subprocess_list=[], resume=False):
    """Function called by directory modules __init__.py to propagate calls down through the workflow tree.
    Note that individual python modules must implement an identical method interface in order to work within the
    workflow management framework.

    Parameters:
        module_name: name of calling module
        submodule_dict: dict containing submodules keyed by name
        subprocess_list: list of subpprocesses to be run in parallel

        resume: Boolean flag to indicate whether execution should continue from the specified
        sub-process. N.B: Only valid when ONE subprocess is specified.
    """

    CONFIG = ProcessorConfig()
    def check_thread_exception():
        """"Check for exception raised by previous thread and raise it if found.
        Note that any other threads already underway will be allowed to finish normally.
        """
        global thread_exception
        logger.debug('thread_exception: %s', thread_exception)
        # Check for exception raised by previous thread and raise it if found
        if thread_exception:
            logger.error('Thread error: ' + thread_exception.message)
            raise thread_exception # Raise the exception in the main thread

    def sequential_process(subprocess_list, resume):
        global thread_exception
        logger.debug('sequential_process(%s, %s) called', subprocess_list, resume)
        for subprocess in subprocess_list:
            check_thread_exception()
            if type(subprocess) == str:
                submodule_list = subprocess.split('.')

                submodule = submodule_list.pop(0)
                assert submodule_dict.get(submodule), 'Sub-module %s does not exist in directory %s' % (submodule, module_name)

                if submodule_list:
                    _subsubprocess_list = ['.'.join(submodule_list)]
                else:
                    _subsubprocess_list = []

                submodule_dict[submodule].process(_subsubprocess_list, resume) # Call individual submodule
                result = submodule # Return last submodule called

            elif type(subprocess) == list: # Process sub-list as required
                process(module_name, submodule_dict, subprocess, resume)
                result = None

        return result

    def parallel_process(subprocess_list, resume):
        '''Function to (eventually) launch separate threads to parallel-process all tasks in the
        supplied subprocess_list

        Parameters:
            subprocess_list: list of subprocesses to be run in parallel

            resume: Boolean flag to indicate whether execution should continue from the specified
            sub-process. N.B: Only valid when ONE subprocess is specified.
        '''
        logger.debug('parallel_process(%s, %s) called', subprocess_list, resume)
        global thread_exception

        def thread_execute(process, subprocess_list, resume):
            """Helper function to capture exception within the thread and set a global
            variable to be checked in the main thread
            N.B: THIS FUNCTION RUNS WITHIN THE SPAWNED THREAD
            """
            global thread_exception
            try:
                process(subprocess_list, resume)
            except Exception, e:
                thread_exception = e
                log_multiline(logger.error, traceback.format_exc(), 'Error in thread: ' + e.message, '\t')
                raise thread_exception # Re-raise the exception within the thread
            finally:
                logger.debug('Thread finished')

        thread_list = []
        for subprocess in subprocess_list:
            check_thread_exception()
            if type(subprocess) == str:
                submodule_list = subprocess.split('.')

                submodule = submodule_list.pop(0)
                assert submodule_dict.get(submodule), 'Sub-module %s does not exist in directory %s' % (submodule, module_name)

                if submodule_list:
                    _subsubprocess_list = ['.'.join(submodule_list)]
                else:
                    _subsubprocess_list = []

#                submodule_dict[submodule].process(_subsubprocess_list, resume) # Call individual submodule
                process_thread = threading.Thread(
                    target=thread_execute,
                    args=(submodule_dict[submodule].process, _subsubprocess_list, resume),
#                    name=submodule + '.' + _subsubprocess_list
                    )
                thread_list.append(process_thread)
                process_thread.setDaemon(False)
                process_thread.start()
                logger.debug('Started thread for %s.process(%s, %s)', submodule, _subsubprocess_list, resume)

                result = submodule # Return last submodule called

            elif type(subprocess) == list:
#                process(module_name, submodule_dict, subprocess, resume) # Process sub-list as required
                process_thread = threading.Thread(
                    target=thread_execute,
                    args=(process, module_name, submodule_dict, subprocess, resume), # Process sub-list as required
#                    name=submodule + '.' + _subsubprocess_list
                    )
                thread_list.append(process_thread)
                process_thread.setDaemon(False)
                process_thread.start()
                logger.debug('Started thread for process(%s, %s, %s, %s)', module_name, submodule_dict, subprocess, resume)
                result = None

        # Wait for all threads to finish
        for process_thread in thread_list:
            check_thread_exception()
            process_thread.join()

        check_thread_exception()
        logger.debug('All threads finished in %s.', subprocess_list)

        return result


    logger.debug('%s.process(%s, %s, %s, %s) called', __name__, module_name, submodule_dict, subprocess_list, resume)

    # Change string argument to list if necessary
    if type(subprocess_list) == str:
        subprocess_list = subprocess_list.split(',')

    #TODO: Make this check more thorough
    assert not resume or len(subprocess_list) <= 1, 'resume only valid when single or no subprocess specified'

    _subprocess_list = subprocess_list or get_subprocess_list(module_name) # Use full process list if none provided
    logger.debug('_subprocess_list = %s', _subprocess_list)

    submodule = None
    if len(_subprocess_list) == 1 and type(_subprocess_list[0]) == list: # [[process1, process2...processn]] indicates parallel processing is possible
        while len(_subprocess_list) == 1 and type(_subprocess_list[0]): # Use innermost nested list
            _subprocess_list = _subprocess_list[0]
        if CONFIG.sequential: # Process sequentially
            submodule = sequential_process(_subprocess_list, resume)
        else: # Process in parallel
            submodule = parallel_process(_subprocess_list, resume)
    else: # Process sequentially
        submodule = sequential_process(_subprocess_list, resume)

    # Resume
    # TODO: Need to respect parallel processing instead of just flattening list of remaining tasks
    if resume and submodule:
        # Process all downstream tasks if resume is required
        _subprocess_list = flatten_process_list(get_subprocess_list(module_name))
        logger.debug('submodule = %s', submodule)
        logger.debug('Original _subprocess_list = %s', _subprocess_list)
#        while _subprocess_list and re.sub('\..*', '', _subprocess_list.pop(0)) != submodule:
        while _subprocess_list and _subprocess_list.pop(0) != submodule:
            pass

        if _subprocess_list:
            logger.info('Resuming execution from %s', _subprocess_list[0])
            sequential_process([_subprocess_list[0]], True)





class ProcessorConfig(object):
    """
    ProcessorConfig class to provide access to global NBAR configuration options held in a
    single conf file. Options are accessed as class attributes ensuring that the same values
    are read by all class instances. Note that options are looked up within the __setattr__
    function the first time they are accessed and then they are stored as "real" class
    attributes for subsequent accesses.

    Care must be taken to ensure that this class is imported with exactly the same module
    name (i.e. "from python.processor_config import ProcessorConfig") to ensure that the class is added
    to the same namespace.

    Do NOT attempt to access any attributes of this class until after the config file has
    been opened, otherwise an assertion error will be raised. Note that only one config file
    can be opened - repeated open_config calls will raise an assertion error.

    """

    __metaclass__ = Singleton

    _TYPEDICT = {'bool': bool, 'boolean': bool, 'int': int, 'integer': int,
                  'float': float, 'datetime': datetime, 'date': date, 'time': time,
                  'list': list, 'dict': dict
                  }

    _FORMATDICT = {datetime: ['%Y-%m-%d %H:%M:%S', '%Y-%m-%d %H:%M:%S.%f'],
                   date: ['%Y-%m-%d'],
                   time: ['%H:%M:%S', '%H:%M:%S.%f']
                   }

    PROCESS_LEVELS = ['l1t', 'nbar', 'pqa', 'fc']

    def __init__(self, config_file=None):
        """
        Constructor with optional parameter for config file path

        """

        self._config_file = None # Name of configuration file opened
        self._configParser = None
        self._section_name = None

        self.log_file = None

        self.station_code = {} # Dict containing all ground station IDs and codes

        self.pqa_test_consts = [] # List containing all PQA test constant names defined in conf file
        self.pqa_test_index = {} # Dict containing all PQA test bit indices keyed by constant names
        self.pqa_test_description = {} # Dict containing all PQA test descriptions keyed by constant names
        self.pqa_param = {} # Dict containing all PQA parameters and corresponding values

        # Dict of non-string configuration parameters
        self._option_types = {}

        # Allow for multiple input paths
        self.input_levels = None
        self.input = {} # Dict keyed by level containing input root, path & name

        self.work_path = None
        self.code_root = None

        self._args = None

        self.svn_revision = 'Unknown'


        def parse_args():
            """
            Parse arguments to the NBAR command line driver.

            :return:
                argparse namespace object

            """
            logger.debug('  Calling parse_args()')

            _arg_parser = argparse.ArgumentParser('process')

            # N.B: modtran_root is a direct overrides of config entries
            # and its variable name must be prefixed with "_" to allow lookup in conf file
            _arg_parser.add_argument('-C', '--config', dest='config_file',
                default=os.path.join(os.path.dirname(__file__), 'processor.conf'),
                help='Image processor configuration file')

            _arg_parser.add_argument('-1', '--l1t', dest='l1t_path', default=None,
                help='L1T dataset directory')
            _arg_parser.add_argument('-n', '--nbar', dest='nbar_path', default=None,
                help='NBAR dataset directory')
            _arg_parser.add_argument('-q', '--pqa', dest='pqa_path', default=None,
                help='PQA dataset directory')
            _arg_parser.add_argument('-f', '--fc', dest='fc_path', default=None,
                help='FC dataset directory')
            #_arg_parser.add_argument('-t', '--tc', dest='tc_path', default=None,
            #    help='TC dataset directory')

            # Optional parameters used with relative or auto-calculated dataset paths
            _arg_parser.add_argument('--l1t-root', dest='l1t_root', default=None,
                help='L1T root directory')
            _arg_parser.add_argument('--nbar-root', dest='nbar_root', default=None,
                help='NBAR root directory')
            _arg_parser.add_argument('--pqa-root', dest='pqa_root', default=None,
                help='PQA root directory')
            _arg_parser.add_argument('--fc-root', dest='fc_root', default=None,
                help='FC root directory')
            #_arg_parser.add_argument('--tc-root', dest='tc_root', default=None,
            #    help='TC root directory')

            _arg_parser.add_argument('-w', '--work', dest='work_path', default=None,
                help='Work directory')
            _arg_parser.add_argument('-l', '--log', dest='log_file', default=None,
                help='Log File')
            _arg_parser.add_argument('--modtran-root', dest='_modtran_root', default=None,
                help='MODTRAN root directory')
            _arg_parser.add_argument('-P', '--process_level', dest='_process_level', default='nbar',
                help='Process level to run')
            _arg_parser.add_argument('--subprocess', dest='subprocess', default=None,
                help='Individual subrocess(es) to run')
            _arg_parser.add_argument('--subprocess-file', dest='subprocess_file', default=None,
                help='Subprocess(es) to run defined in a text file')
            _arg_parser.add_argument('--resume', dest='resume', default=False,
                action='store_const', const=True,
                help='Resume execution from specified process')
            _arg_parser.add_argument('--debug', dest='debug', default=False,
                action='store_const', const=True,
                help='Create debug files in working directory')
            _arg_parser.add_argument('--sequential', dest='sequential', default=False,
                action='store_const', const=True,
                help='Process all tasks sequentially (i.e. no parallel processing)')
            _arg_parser.add_argument('--repackage', dest='repackage', default=False,
                action='store_const', const=True,
                help='Repackage existing output only (i.e. do not re-create output)')

            # Command line options for production NBAR
            _arg_parser.add_argument('-c', '--constraintid', dest='constraint_id', default=None,
                help='constraint ID')
            _arg_parser.add_argument('-L', '--li-source-description', dest='li_source_description', default=None,
                help='Lineage source description')
            _arg_parser.add_argument('-p', '--purpose', dest='purpose', default=None,
                help='Purpose')
            _arg_parser.add_argument('-d', '--dataset-id', dest='dataset_id', default=None,
                help='Full dataset ID')

            return _arg_parser.parse_args()

        logger.debug('Calling ProcessorConfig(%s)', repr(config_file))

        # Set code_root to parent directory of common
        self.code_root = os.path.normpath(os.path.join(os.path.dirname(__file__), "../.."))

        if not self._configParser:
            self._args = parse_args()

            # Default conf file is processor.conf - show absolute pathname in error messages
            config_file = os.path.abspath(config_file or
                                          self._args.config_file or
                                          os.path.join(os.path.dirname(__file__), 'processor.conf'))

            if not self._configParser:
                self.open_config(config_file)

    def open_config(self, config_file):
        """Function to open configuration file.
        Needs to be called before attempting to access data and should only ever be called once.
        Optional parameter config_file for configuration file path - defaults to processor.conf if not provided"""


        ## Nested function not exploiting the closure isn't the best idiom.
        ## Need to do a review of the functionality that bleongs on this class
        ## suitabliy split this class up into a class and set of functions for
        ## handling config.
        ##
        ## Michael O

        def _read_option_types():
            """Function to create a dict containing (option_name: option_type) pairs.
            This function will ignore any invalid types
            """
            logger.debug('  Calling _read_option_types()')

            assert self._configParser.has_section("option_types"), "Section [option_types] does not exist in config_file"

            for option in self._configParser.options('option_types'):
                option_type = ProcessorConfig._TYPEDICT.get(self._configParser.get('option_types', option))
                if type(option_type) == type:
                    self._option_types[option] = option_type

        def _try_host_section(hostname):
            """Function to check for config file section for host name or alias.
            Sets and returns ProcessorConfig._section_name if matching section found.
            """
            logger.debug('  Calling _try_host_section(%s)', repr(hostname))
            if hostname:
                if self._configParser.has_section(hostname):
                    # Try raw hostname
                    self._section_name = hostname
                else:
                    # Try stripping trailing digits from hostname
                    host_alias = re.match('(.*?)(\d*)(\.|$)', hostname).group(1)
                    if self._configParser.has_section(host_alias):
                        self._section_name = host_alias

            return self._section_name

        def _set_logging():
            """Function to set logging for individual modules
            N.B: Logger names are of the format "nbar.<module_name>"
            """
            logger.debug('  Calling _set_logging()')

            if not self._configParser.has_section("logging"): # Section [logging] does not exist in config_file
                return

#            if self.__getattr__('debug'):
#                console_handler.setLevel(logging.DEBUG) # Show debug messages on console

            for module_name in self._configParser.options('logging'):
                log_level_string = self._configParser.get('logging', module_name).upper()
                log_level = getattr(logging, log_level_string) # Convert from string to numeric
                assert isinstance(log_level, int), 'Invalid logging level: ' + log_level_string

                module_logger = logging.getLogger('root.' + module_name)
                module_logger.setLevel(log_level)
                logger.debug('    Logger %s set to log level %s', module_logger.name, log_level_string)

        def read_ground_stations():
            assert self._configParser.has_section('ground_stations'), 'No ground_stations section in conf file'

            for _ground_station in self._configParser.options('ground_stations'):
                self.station_code[_ground_station.upper()] = self._configParser.get('ground_stations', _ground_station)
            log_multiline(logger.debug, self.station_code, 'Ground stations & codes: ', '\t')

        def read_pqa_tests():
            ''' Reads PQA test constants with corresponding bit indices and descriptions into dicts.
            Saturation indices are read into a single sub-dict
            This is necessary because of the effective re-use of bit 5 for both saturation_60 & saturation_61
            '''
            assert self._configParser.has_section('pqa_tests'), 'No pqa_tests section in conf file'

            saturation_dict = {}
            self.pqa_test_index['SATURATION'] = saturation_dict
            for _pqa_test_const in [_pqa_test_const.upper() for _pqa_test_const in self._configParser.options('pqa_tests')]:
                # Get integer index and string description
                _index_desc = self._configParser.get('pqa_tests', _pqa_test_const).split(',')
                _bit_index = int(_index_desc[0])
                _description = _index_desc[1]

                m = re.match('SATURATION_(\d{2})', _pqa_test_const)
                if m: # Saturation for individual bands
                    saturation_dict[int(m.group(1))] = _bit_index

                self.pqa_test_index[_pqa_test_const] = _bit_index
                self.pqa_test_consts.append(_pqa_test_const)
                self.pqa_test_description[_pqa_test_const] = _description

            log_multiline(logger.debug, self.pqa_test_index, 'PQA Tests & bit indices: ', '\t')
            log_multiline(logger.debug, self.pqa_test_consts, 'PQA Test constants: ', '\t')
            log_multiline(logger.debug, self.pqa_test_description, 'PQA Test descriptions: ', '\t')

        def read_pqa_params():
            assert self._configParser.has_section('pqa_params'), 'No pqa_params section in conf file'

            for _pqa_param in self._configParser.options('pqa_params'):
                self.pqa_param[_pqa_param.lower()] = self._configParser.getfloat('pqa_params', _pqa_param)
            log_multiline(logger.debug, self.pqa_param, 'PQA Params and Values: ', '\t')

        def get_svn_revision():
            """
            Returns a string representing the current SVN revision number or None if unknown
            """
            s = re.search('\nRevision: (.*)\n',
                          execute(command_string='svn info ' + os.path.abspath(os.path.dirname(os.path.dirname(__file__))))['stdout']
                          )
            if s and s.group(1):
                return s.group(1)
            else:
                return None

        def get_git_version():
            """
            Returns a string representing the current Git version or None if unknown
            """
            result = None
            s = execute(cwd=os.path.dirname(os.path.dirname(os.path.dirname(__file__))), # __file__ should be <repository>/ULA3/image_processor/__init.py
                        command_string='git show --pretty=oneline ')['stdout']

            if s and len(s) > 40:
                m = re.match('^(\w{40})\s+.*', s)
                if m:
                    result = m.group(1)

            return result

        logger.debug('ProcessorConfig.open_config(%s)', repr(config_file))

        assert os.path.exists(config_file), config_file + " does not exist"

        logger.debug('  Opening conf file %s', repr(config_file))
        self._configParser = ConfigParser.SafeConfigParser(allow_no_value=True)
        self._configParser.read(config_file)

        _read_option_types()

        # Read PBS_O_HOST environment variable if current hostname not found
        if not _try_host_section(socket.gethostname()):
            if not _try_host_section(os.getenv('PBS_O_HOST', None)):
                # For jobs at NCI the hostname will be a node name, and not
                # 'vayu', 'dcc' or 'xe'. Try to work out the hostname from
                # $PBS_JOBID: dcc => '123.dccpbs', vayu => '123.vu-pbs'.
                host = os.getenv('PBS_JOBID', None) # Try deriving cluster name from job ID
                if host: # If this is a PBS job
                    host = re.match('(.*?).([a-z\-]+)pbs$', host).group(2)
                    host = ('vayu' if host == 'vu-' else host)
                    _try_host_section(host)
                else:
                    _try_host_section('default_config') # Fall back to universal Gdefault config section


        assert self._section_name, 'No configuration section found for hostname or alias'
        logger.debug('  Host section name = %s', self._section_name)
#        print("ProcessorConfig successfully opened config file " + config_file)
        self._config_file = config_file

        _set_logging() # Apply desired logging level to all specified modules


        print '__name__ == %s' % __name__
        if __name__ == '__main__':
            # Define which datasets are required as input for each processing level. It might be good to refactor this
            # into level-specific modules rather than having input requirements hard-coded here. (AI)
            self.process_level = self._args._process_level.lower() or 'nbar'
            if self.process_level == 'nbar' or self.process_level == 'tc':
                self.input_levels=['l1t']
            elif self.process_level == 'pqa':
                self.input_levels=['l1t','nbar']
            elif self.process_level == 'fc':
                self.input_levels=['nbar']
            else:
                assert False, 'Invalid process ' + self.process + ' specified. Must be nbar, pqa, fc or tc.'
    
            self.modtran_root = self._args._modtran_root or self.__getattr__('modtran_root')
    
            # Set ProcessorConfig.BIN_DIR from conf file or next to 'common' dir
            # Simon K: Why is this in upper case (seems to imply constant)?
            try:
                self.__getattr__('BIN_DIR')
            except(NameError):
                self.BIN_DIR = os.path.join(self.code_root, 'bin')
    
            # Set up inputs
            for _level in self.input_levels:
                _input_path = vars(self._args).get(_level + '_path')
                assert _input_path, 'No ' + _level + ' input path provided'
    
                _input_root = vars(self._args).get(_level + '_root') or self.__getattr__(_level + '_data_root')
                assert _input_root, 'No ' + _level + ' data root provided'
    
                _input_name = os.path.basename(_input_path)
    
                if not os.path.exists(_input_path) and _input_path[0] != os.pathsep: # Try relative path
                    _input_path = os.path.join(_input_root, _input_path)
                _input_path = os.path.abspath(_input_path)
    
                assert os.path.exists(_input_path), _input_path + ' does not exist'
                assert os.path.isdir(_input_path), _input_path + ' is not a directory'
    
                if not self.work_path:
                    # Determine work directory path
                    if not self._args.work_path:
                        _work_path = os.path.join(self.__getattr__('work_root'), _input_name)
                    else:
                        _work_path = self._args.work_path
                        if not os.path.exists(os.path.dirname(_work_path)) and self._args.work_path[0] != os.pathsep:
                            _work_path = os.path.join(self.__getattr__('work_root'), _work_path)
    
                    _work_path = os.path.abspath(_work_path)
                    assert os.path.exists(os.path.dirname(_work_path)), 'Work path ' + _work_path + ' does not exist'
    
                    if not os.path.exists(_work_path):
                        os.mkdir(_work_path)
                        logger.info('created work directory %s', _work_path)
    
                    self.work_path = _work_path
                    logger.debug('work_path = %s', self.work_path)
    
                self.input[_level] = {}
                logger.debug('data root for %s = %s', _level, _input_root)
                self.input[_level]['root'] = _input_root
                logger.debug('input path for %s = %s', _level, _input_path)
                self.input[_level]['path'] = _input_path
                logger.debug('input name for %s = %s', _level, _input_name)
                self.input[_level]['name'] = _input_name
    
                # NBAR & TC both use NBAR paths for output - this could be removed if TC replaces NBAR
                if self.process_level == 'tc':
                    self.output_level = 'nbar'
                else:
                    self.output_level = self.process_level
    
            # Output path values will be used in downstream modules after metadata has been read
            # N.B: output_path may be left as None to be filled in downstream
            self.output_path = vars(self._args).get(self.output_level + '_path')
            #TODO: Fix this so that we don't need to define tc_root in the conf file
            self.output_root = vars(self._args).get(self.output_level + '_root') or self.__getattr__(self.output_level + '_data_root')
    
            # Pre-pend root if nonexistent relative directory specified for output
            if self.output_path:
                if not os.path.exists(os.path.dirname(self.output_path)) and self.output_path[0] != os.pathsep:
                    self.output_path = os.path.join(self.output_root, self.output_path)

        if self.work_path:
            # Set up log file - defaults to process.log in work directory
            self.log_file = (self._args.log_file or
                                        os.path.join(self.work_path, 'process.log'))

            file_handler = logging.FileHandler(self.log_file)
            file_handler.setLevel(logging.DEBUG)
            file_formatter = logging.Formatter('%(levelname)s|%(name)s|%(asctime)s|%(message)s')
            file_handler.setFormatter(file_formatter)
            root_logger.addHandler(file_handler)

        self.svn_revision = get_svn_revision() or 'Unknown'
        self.git_version = get_git_version() or 'Unknown'

        read_ground_stations()

        read_pqa_tests()
        read_pqa_params()

    def __getattr__(self, attribute_name):
        """return value of specified attribute from cache or from previously
        read config file"""

        def parse_list(list_string, element_type=str):
            return [element_type(element.strip())
                    for element in re.search('\s*(?<=\[)(.*)(?=\])\s*', list_string).group(0).split(',')]

        def parse_dict(list_string, key_type=str, value_type=str):
            result = dict()
            for key, value in [pair.strip().split(':') for pair in re.search('\s*(?<=\{)(.*)(?=\})\s*', list_string).group(0).split(',')]:
                result[key_type(key.strip())] = value_type(value.strip())
            return result

        logger.debug('Calling ProcessorConfig.__getattr__(%s)', repr(attribute_name))

        assert self._configParser is not None, "No configuration file opened"
        assert self._section_name, "No host section defined"

        try: # Check command line arguments first
            attribute_value = self._args.__dict__[attribute_name.lower()]
        except (KeyError): # Then check .conf file
            logger.debug('%s not found in command line arguments', repr(attribute_name))
            attribute_type = self._option_types.get(attribute_name.lower(), str)

            if attribute_type == bool:
                get_function = self._configParser.getboolean
            elif attribute_type == int:
                get_function = self._configParser.getint
            elif attribute_type == float:
                get_function = self._configParser.getfloat
            else:
                get_function = self._configParser.get

            try:
                # Check host-specific section first
                attribute_value = get_function(self._section_name, attribute_name)
                logger.debug('  %s found in conf file as %s', attribute_name, repr(attribute_value))
            except ConfigParser.NoOptionError:
                # Check global_config section if not found in host-specific section
                try:
                    attribute_value = get_function('global_config', attribute_name)
                except ConfigParser.NoOptionError:
                    # Raise NameError instead of ConfigParser.NoOptionError if option not found
                    raise NameError("option '" + attribute_name + "' is not defined in " + self._config_file)

            if attribute_type == list:
                # Attempt to read in a list of integers, floats or strings
                for element_type in [int, float, str]:
                    try:
                        # Convert value from string to list of specified element type
                        attribute_value = parse_list(attribute_value, element_type)
                        break
                    except ValueError:
                        continue

            if attribute_type == dict:
                # Convert value from string to dict
                attribute_value = parse_dict(attribute_value, int, int)
            elif attribute_type in ProcessorConfig._FORMATDICT:
                # Post-process string for date/time conversion
                formatted_value = None
                for format_string in ProcessorConfig._FORMATDICT[attribute_type]:
                    try:
                        logger.debug('Trying format %s for %s object %s',
                                     format_string,
                                     attribute_type.__name__,
                                     repr(format_string)
                                     )
                        if attribute_type == datetime:
                            formatted_value = datetime.strptime(attribute_value, format_string)
                        elif attribute_type == date:
                            formatted_value = datetime.strptime(attribute_value, format_string).date()
                        elif attribute_type == time:
                            formatted_value = datetime.strptime(attribute_value, format_string).time()
                        break
                    except(ValueError):
                        pass

                assert formatted_value, 'Invalid ' + attribute_type.__name__ + ' format : ' + attribute_value
                attribute_value = formatted_value

        # Cache retrieved value as new class attribute
        self.__setattr__(attribute_name, attribute_value)
        logger.debug('  Cached configuration value: %s = %s', attribute_name, attribute_value)

        return attribute_value

    @property
    def config_file(self):
        """Get the current config_file."""
        return self._config_file
