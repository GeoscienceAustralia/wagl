from __future__ import print_function
from __future__ import absolute_import

try:
    from importlib import metadata as _md
except ImportError:
    # Running on pre-3.8 Python; use importlib-metadata package
    import importlib_metadata as _md

try:
    __version__ = _md.version(__name__)
except _md.PackageNotFoundError:
    __version__ = "Not Installed"
