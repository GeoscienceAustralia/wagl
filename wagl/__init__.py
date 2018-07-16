from __future__ import print_function
from __future__ import absolute_import
from ._version import get_versions

from logging import NullHandler

__version__ = get_versions()['version']
del get_versions
