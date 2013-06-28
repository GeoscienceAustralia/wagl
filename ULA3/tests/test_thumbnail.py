import unittest, os
from ULA3.utils import create_thumbnail

class ThumbnailTestCase(unittest.TestCase):
    def runTest(self):
        #outdir = '/short/v10/tmp/NBAR-DEBUG'
        outdir = '/short/v10/tmp/NBAR-DEBUG2'
        f_fmt = 'LS7_ETM_NBAR_P54_GANBAR04_091_081_20090609_B%d.tif'

        def fpath(b):
            return os.path.join(outdir, f_fmt % (b*10))

        create_thumbnail(fpath(7), fpath(4), fpath(1), 'TEST.jpg')