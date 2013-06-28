'''
__init__.py - Created on Jun 14, 2012

@author: Alex Ip (alex.ip@ga.gov.au)

This is a generic __init__.py module which should exist in all workflow management subdirectories.
Its purpose is to recursively import all sub-modules and also to implement a generic
process() method which cascades calls to the process() method implemented in each sub-module.

Note that this module requires an associated __init__.txt file which contains a list defining the tasks
in the default workflow for this module. A double-nested list (e.g. [[task1, task2,...,task<n>]] indicates
that the enclosed tasks are independent and able to be processed in parallel. Note that each task may be
a sub-list for either sequential or parallel processing.
'''

import ULA3.image_processor as process_manager
import logging

logger = logging.getLogger('root.' + __name__)

_submodule_dict = {}

def __import_submodules(submodule_list):
    """Function to import all modules in a list into the current module's namespace
    and populate _submodule_dict with the imported modules keyed by their name.
    Note that this function must be implemented in each directory module's __init__.py
    module in order to import the sub-modules into the module's namespace.

    Parameter:
        submodule_list: A list containing the names of all (potentially) valid Python
        submodules to be imported.
    """
    submodule_dict = {}
    for submodule_name in submodule_list:
        try:
            submodule_dict[submodule_name] = __import__('%s.%s' % (__name__, submodule_name), fromlist = [__name__])
            logger.debug('module %s.%s imported', __name__, submodule_name)
        except (Exception), e:
            logger.error('Unable to import %s.%s: %s', __name__, submodule_name, e.message)
            raise e
    return submodule_dict

def process(subprocess_list=[], resume=False):
    """Generic method to cascade calls to the process() method in sub-modules.
    Calls process_manager.process() function to avoid duplicating code in all __init__.py modules

    If the subprocess_list is not given or empty, then the default workflow for this module
    will be read from the associated __init__.txt file by the process_manager.process() function.

    Parameters:
        subprocess_list: A list representing a path through the module namespace to a module
        which implements a process() method.

        resume: Boolean flag to indicate whether execution should continue from the specified
        sub-process. N.B: Only valid when ONE subprocess is specified.
    """
    logger.debug('%s.process(%s, %s) called', __name__, subprocess_list, resume)

    process_manager.process(__name__, _submodule_dict, subprocess_list, resume)


# Import all valid Python sub-modules in this module's directory
_submodule_dict = __import_submodules(process_manager.get_submodule_list(__path__[0]))
