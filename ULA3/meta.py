"""
Provides utilities for logging and meta programming to ULA3.
"""

import threading, os
from functools import wraps
from inspect import getcallargs

RUNNING_SPHINX = os.environ.has_key('RUNNING_SPHINX_BUILD')

class Singleton(type):
    """
    Metaclass for Singletons.

    We could also keep the singletons in a dictionary in this class with keys of type
    class. I prefer, however, to keep them in the actual class.

    """
    def __new__(cls, name, bases, dct):
        dct['_instance_for_singleton_ssfusousoifusos'] = None
        dct['_lock_for_singleton_ssfusousoifusos'] = threading.Lock()
        return super(Singleton, cls).__new__(cls, name, bases, dct)

    def __call__(cls, *args, **kwargs):
        if not cls._instance_for_singleton_ssfusousoifusos:
            with cls._lock_for_singleton_ssfusousoifusos:
                if not cls._instance_for_singleton_ssfusousoifusos:
                    cls._instance_for_singleton_ssfusousoifusos = super(Singleton, cls).__call__(*args, **kwargs)
        return cls._instance_for_singleton_ssfusousoifusos





def create_arg_string(func, *args, **kwargs):
    """
    Constructs a string of the arguments passed to a function on a given invocation.

    :param func: The function for which the string is to be constructed.

    :param args: The positional arguments passed in the call to ``func``.

    :param kwargs: The keyword arguments passed in the call to ``func``.

    """
    return ',\n\t'.join(['='.join([str(y) for y in item]) for item in getcallargs(func, *args, **kwargs).iteritems()])





def print_call(logger):
    """
    Decorator which prints the call to a function, including all the arguments passed.

    :param func: The function to be decorated.

    :param logger: Callable which will be passed the string representation of the function call. Then no
        logging is performed (the decorated is simply returned.

    """

    def wrap(func):
        if RUNNING_SPHINX or not logger:
            return func

        @wraps(func)
        def wrapper(*args, **kwargs):
            arg_str = create_arg_string(func, *args, **kwargs)
            logger('\ncalling %s.%s(\n\t%s\n)\n' % (func.__module__, func.__name__, arg_str))
            res = func(*args, **kwargs)
            logger("SUCCESS!")
            return res

        return wrapper

    return wrap
