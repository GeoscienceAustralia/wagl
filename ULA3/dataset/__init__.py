"""
Various types of dataset used within ULA3 and the image_processor. These are documented in their
respective submodules, but should be imported from here, not directly from those submodules.

:todo:
    These should be moved to the modules in which they are used. that is:

    - ElevationDataset -> ULA3.ancillary.elevation
    - OzoneDataset -> ULA3.ancillary.ozone
    - WaterVapourDataset -> ULA3.ancillary.water
    - SceneDataset -> ULA3
    - Dataset -> ULA3
"""

from _dataset import Dataset, DSException
from _elevation_dataset import ElevationDataset
from _ozone_dataset import OzoneDataset
from _scene_dataset import SceneDataset
from _water_vapour_dataset import WaterVapourDataset
