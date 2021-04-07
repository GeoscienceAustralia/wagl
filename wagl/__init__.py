from __future__ import print_function
from __future__ import absolute_import

try:
    from importlib import metadata
except ImportError:
    # Running on pre-3.8 Python; use importlib-metadata package
    import importlib_metadata as metadata

try:
    __version__ = metadata.version(__name__)
except metadata.PackageNotFoundError:
    __version__ = "Not Installed"
