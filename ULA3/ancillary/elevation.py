'''
Utilities for getting elevation data for a scene.

:note:
    The other place that elevation data is provided is in :py:func:`ULA3.tc.clip_dsm`, which is used (as
    its location suggests) in the terrain correction algorithm :ref:`tc-algorithm-label`.
'''

import logging
from ULA3.dataset import ElevationDataset
from ULA3.meta import print_call

logger = logging.getLogger('root.' + __name__)

@print_call(logger.info)
def get_elevation_data(lon_lat, dem_dir):
    """
    Get elevation data for a scene.

    :param lon_lat:
        The latitude, longitude of the scene center.
    :type lon_lat:
        float (2-tuple)

    :dem_dir:
        The directory in which the DEM can be found.
    :type dem_dir:
        str
    """
    elevation_dataset = ElevationDataset(dem_dir)
    elevation_value = elevation_dataset.get_data_value(lon_lat)

    elevation_data = {
        'data_source': 'Elevation',
        'data_file': elevation_dataset.pathname,
        'value': elevation_value
        }

    return elevation_data

