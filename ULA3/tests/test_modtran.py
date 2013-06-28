import shutil, os
import unittest as ut
from osgeo import gdal, gdalconst
from ULA3.dataset import SceneDataset
from ULA3.modtran import clip_dsm

class ClipDSMTestCase(ut.TestCase):
    """
    Test case for testing the clipping of the DSM for Terrain Correction.

    This test uses the 'real' (i.e. an actual scene and DSM) for the test. The DSM is a file of approximately
    70GB and the scene itself is also a relitively large dataset. Hence this test is specific to a particular
    system configuration... the NCI as it was on 05/04/2013 and hence can only be run on a system which has
    the same filesystem structure.

    To change this, you should only have to change the setup method, which is where the file system specifics
    are encapsulated.

    Hopefully, in the near future, we will have a test harness that can be more portable, and testing will be
    possible elsewhere.

    """
    def setUp(self): #needed to ensure that tearDown is called.
        self.l1t_path = "/g/data/v10/L1/2009-02/LS5_TM_OTH_P51_GALPGS01-002_091_085_20090225"
        self.dsm_path = "/short/v10/sok547/tc/dsm/dsm1sv1_0_Clean.img"
        self.output_dir = "/short/v10/sok547/tc/dsm_clip_test"
        self.original_clip_example_path = "/short/v10/sok547/tc/091_085_20090225_dem.img"
        self.output_format = "ENVI"

    def tearDown(self):
        shutil.rmtree(self.output_dir, True)

    def test_clip1(self):
        self.output_path = os.path.join(self.output_dir, "clipped.img")
        self.buffer_size = 250
        self.scene_data = SceneDataset(self.l1t_path)
        return clip_dsm(self.scene_data, self.dsm_path, self.output_path, self.buffer_size, self.output_format)

    def test_clipCompare(self):
        self.output_path = os.path.join(self.output_dir, "clipped_compare.img")
        self.buffer_size = 0
        self.scene_data = gdal.Open(self.original_clip_example_path, gdalconst.GA_ReadOnly)
        return clip_dsm(self.scene_data, self.dsm_path, self.output_path, self.buffer_size, self.output_format)





def suite():
    tests = ['test_clipCompare', 'test_clip1']
    return ut.TestSuite(map(ClipDSMTestCase, tests))
