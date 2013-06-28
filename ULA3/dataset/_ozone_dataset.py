#!/usr/bin/env python
"""
OzoneDataset class for ancillary data
Subclassed from Dataset

Author: Alex Ip (alex.ip@ga.gov.au)
Includes code from Roger Edberg's original anci_geotiff_loader.py
"""
import os, logging

from . import Dataset

logger = logging.getLogger('root.' + __name__)

class OzoneDataset(Dataset):
    """
    Dataset for handle Ozone ancillary data
    """
    def __init__(self, ozone_data_root, ancillary_datetime):
        """
        Constructor.

        :param ozone_data_root:
            Root directory for Ozone ancillary data.
        :type ozone_data_root:
            str

        :param ancillary_datetime:
            ???
        :type ancillary_datetime:
            :py:class:`datetime.datetime`
        """
        self.ancillary_datetime = ancillary_datetime
        pathname = os.path.join(ozone_data_root, "%s.tif" % ancillary_datetime.strftime('%b').lower())
        assert os.path.exists(pathname), 'Ozone ancillary data file ' + pathname + ' does not exist'
        super(OzoneDataset, self).__init__(pathname)
