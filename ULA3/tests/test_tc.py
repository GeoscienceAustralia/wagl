"""
Tests for the Terrain Correction module.

Since we are lacking a portable test harness, these tests are configured to run on the NCI, but should
also work on DCC.


"""

import os, gc
import numpy as np
import unittest as ut
from osgeo import gdal, gdalconst
from ULA3.utils import Buffers, ImageShape
from ULA3.dataset import SceneDataset
from ULA3.tc import clip_dsm, run_castshadow, run_slope, run_brdfterrain
from ULA3.utils import load_bin_file, read_array_int8, as_array, dump_array
from ULA3.filter import filter_float as filter


SYS_BASE_DIR = '/g/data/v10'

TEST_OUTPUT_BASE_DIR = '/short/v10'

# the location of the test L1T dataset to use
L1T_PATH = os.path.join(SYS_BASE_DIR, 'ULA3-TESTDATA', 'L1', 'UTM', 'LS5_TM_OTH_P51_GALPGS01-002_091_085_20090225')

# the path to the dsm data (a national grid)
DSM_PATH = os.path.join(SYS_BASE_DIR, "eoancillarydata", "elevation", "tc_aus_3sec", "dsm1sv1_0_Clean.img")

# where the data initially provided by Fuqin lives
TEST_DATA_DIR = os.path.join(SYS_BASE_DIR, 'ULA3-TESTDATA', 'TC-FUQINS-REGION')

# the base directory under which the outputs and intermediate files will be written
TEST_OUTPUT_BASE = os.path.join(TEST_OUTPUT_BASE_DIR, os.getenv("USER"), "tc")

# the number of rows in the test region (i.e. the region in ``TEST_DATA_DIR``)
TEST_REGION_NROW = 8800

# the number of columns in the test region (i.e. the region in ``TEST_DATA_DIR``)
TEST_REGION_NCOL = 8800

# the prefix used for files in the test region (i.e. the region in ``TEST_DATA_DIR``)
# this reflects just how Fuqin named her original files.
REGION_DATE_PREFIX = "091_085_20090225"

# the path that the test terrain corrected image will be written to.
TEST_TC_OUTPUT_PATH = os.path.join(TEST_OUTPUT_BASE, "ref_terrain_%i.img")

# should we be clipping the DSM for the tests using Fuqin's original data? If True, then
# TEST_DSM_CLIP_OUTPUT_PATH must exist, which is most easily achieved by running
# ClipDSMTestCase.test_fuqin_test_data_smooth, with RUN_WITH_NEW_DSM = True.
RUN_WITH_NEW_DSM = True

# the path that the clipped (and smoothed) DSM will be written to
TEST_DSM_CLIP_OUTPUT_PATH = os.path.join(TEST_OUTPUT_BASE, "clipped_compare_to_orig_smoothed.img")





def compare_images(a, b, name):
    """
    Compare two arrays for any differences.
    """
    dif = np.max(np.abs(a-b))
    print "max differences in %s: %f" % (name, dif)
    return dif




def run_slope_fuq():
    is_utm = False
    pb = Buffers(left=251, right=249, top=249, bottom=251)
    shape = ImageShape(
        origin_x=146.750125, origin_y=-34.916875,
        pix_sz_x=0.00025, pix_sz_y=0.00025,
        dim_y=TEST_REGION_NROW, dim_x=TEST_REGION_NCOL)

    dem_path = TEST_DSM_CLIP_OUTPUT_PATH if RUN_WITH_NEW_DSM else os.path.join(TEST_DATA_DIR, REGION_DATE_PREFIX + "_dem.img")
    sol_zenith_angle_path = os.path.join(TEST_DATA_DIR, REGION_DATE_PREFIX + "_solar_angle.bin")
    sol_azimuth_angle_path = os.path.join(TEST_DATA_DIR, REGION_DATE_PREFIX + "_sazi_angle.bin")
    view_zenith_angle_path = os.path.join(TEST_DATA_DIR, REGION_DATE_PREFIX + "_view_angle.bin")
    view_azimuth_angle_path = os.path.join(TEST_DATA_DIR, REGION_DATE_PREFIX + "_azi_angle.bin")

    dem_data = as_array(dem_path, dtype=np.float32)
    solar_zenith_data = load_bin_file(sol_zenith_angle_path, TEST_REGION_NROW, TEST_REGION_NCOL, dtype=np.float32)
    solar_azimuth_data = load_bin_file(sol_azimuth_angle_path, TEST_REGION_NROW, TEST_REGION_NCOL, dtype=np.float32)
    view_zenith_data = load_bin_file(view_zenith_angle_path, TEST_REGION_NROW, TEST_REGION_NCOL, dtype=np.float32)
    view_azimuth_data = load_bin_file(view_azimuth_angle_path, TEST_REGION_NROW, TEST_REGION_NCOL, dtype=np.float32)

    return run_slope(
        shape,
        dem_data,
        solar_zenith_data,
        view_zenith_data,
        solar_azimuth_data,
        view_azimuth_data,
        pb,
        is_utm)





class TCTestCase(ut.TestCase):
    """
    A :py:class:`unittest.TestCase` used as a base class for all other test cases defined herein.
    its function is to create the base directories for the test case, and (if ``tearDown`` is
    uncommented) remove them afterwards.
    """
    def baseSetUp(self, output_leaf_dir):
        self.output_dir = os.path.join(TEST_OUTPUT_BASE, output_leaf_dir)
        self.output_format = "ENVI"
        try:
            os.makedirs(self.output_dir)
        except Exception, e:
            if not os.path.exists(self.output_dir):
                # all sorts of other things might be wrong here... but those things will manifest later.
                raise e

    #def tearDown(self):
    #    shutil.rmtree(self.output_dir, True)



class ClipDSMTestCase(TCTestCase):
    """
    Test case for testing the clipping of the DSM for Terrain Correction.

    This test uses the 'real' (i.e. an actual scene and DSM) for the test. The DSM is a file of approximately
    70GB and the scene itself is also a relitively large dataset. Hence this test is specific to a particular
    system configuration... the NCI as it was on 05/04/2013 and hence can only be run on a system which has
    the same filesystem structure.

    To change this, you should only have to change the setup method, which is where the file system specifics
    are encapsulated.
    """

    def setUp(self): #needed to ensure that tearDown is called.
        super(ClipDSMTestCase, self).baseSetUp("dsm_clip_test")


    def test_clip_scene(self):
        """
        Test of clipping a scene for an L1T image.
        """
        output_path = os.path.join(self.output_dir, "clipped_to_l1t.img")
        buffer_size = Buffers(250)
        scene_data = SceneDataset(L1T_PATH)
        clipped_data = clip_dsm(scene_data, DSM_PATH, output_path, buffer_size, self.output_format)
        expected_shape = (
            scene_data.shape[0] + buffer_size.top + buffer_size.bottom,
            scene_data.shape[1] + buffer_size.left + buffer_size.right)
        self.assertTupleEqual(clipped_data.shape, expected_shape)


    def __do_clip(self, output_name):
        original_clip_example_path = os.path.join(TEST_DATA_DIR, REGION_DATE_PREFIX + "_dem.img")
        output_path = os.path.join(self.output_dir, output_name)
        buffer_size = Buffers(0)
        scene_data = gdal.Open(original_clip_example_path, gdalconst.GA_ReadOnly)
        return (scene_data, clip_dsm(scene_data, DSM_PATH, output_path, buffer_size, self.output_format))


    def test_fuqin_test_data_unsmooth(self):
        """
        Test clipping of Fuqin's region.
        """
        scene_data, clipped_data = self._ClipDSMTestCase__do_clip("clipped_compare_to_orig_unsmoothed.img")
        self.assertTupleEqual(clipped_data.shape, (scene_data.RasterYSize, scene_data.RasterXSize))
        #TODO: assert that the values are close enough to what is expected.


    def test_fuqin_test_data_smooth(self):
        """
        Test clipping and smoothing of Fuqin's region.
        """
        scene_data, clipped_data = self._ClipDSMTestCase__do_clip("clipped_compare_to_orig_unsmoothed.img")
        smoothed = filter(clipped_data)

        if RUN_WITH_NEW_DSM:
            dump_array(
                array=smoothed,
                output_path=TEST_DSM_CLIP_OUTPUT_PATH,
                template_dataset=scene_data,
                file_format="ENVI")

        #TODO: assert that the values are close enough to what is expected... once that has been specified.





class RunCastShadowTestCase(TCTestCase):
    def setUp(self):
        super(RunCastShadowTestCase, self).baseSetUp("cast_shadow_test")


    def __do_random(self, nrow, ncol, buffer_width, block_height, block_width):
        is_utm = False

        shape = ImageShape(origin_x=0., origin_y=0., pix_sz_x=0.00025, pix_sz_y=0.00025, dim_y=nrow, dim_x=ncol)

        pb = Buffers(buffer_width)
        dem_data = np.random.ranf((nrow+pb.top+pb.bottom, ncol+pb.left+pb.right)).astype(np.float32)
        solar_zenith_data = np.random.ranf((nrow, ncol)).astype(np.float32)
        solar_azimuth_data = np.random.ranf((nrow, ncol)).astype(np.float32)

        mask_all = run_castshadow(
            shape,
            dem_data,
            solar_zenith_data,
            solar_azimuth_data,
            pb,
            block_height,
            block_width,
            is_utm)

        self.assertTupleEqual(shape.shape, mask_all.shape, "dimensions of input and output data do not match.")


    #@ut.skip("too slow for general use")
    def test_small_random(self):
        """
        Test run_castshadow on some random data.
        """
        self._RunCastShadowTestCase__do_random(50, 75, 11, 6, 7)


    def test_scene_size_random(self):
        """
        Runs run_castshadow on some random data approximately the size of a scene.
        """
        l1t = SceneDataset(L1T_PATH)
        self._RunCastShadowTestCase__do_random(l1t.shape[0], l1t.shape[1], 250, 500, 500)


    def test_fuqin_data(self):
        """
        Runs run_castshadow on some intermediate output from an NBAR run
        """
        block_height = block_width = 500
        is_utm = False
        pb = Buffers(left=251, right=249, top=249, bottom=251)
        l1t = ImageShape(
            origin_x=146.750125, origin_y=-34.916875,
            pix_sz_x=0.00025, pix_sz_y=0.00025,
            dim_y=TEST_REGION_NROW, dim_x=TEST_REGION_NCOL)

        dem_path = RUN_WITH_NEW_DSM if TEST_DSM_CLIP_OUTPUT_PATH else os.path.join(TEST_DATA_DIR, REGION_DATE_PREFIX + "_dem.img")
        sol_zenith_angle_path = os.path.join(TEST_DATA_DIR, REGION_DATE_PREFIX + "_solar_angle.bin")
        sol_azimuth_angle_path = os.path.join(TEST_DATA_DIR, REGION_DATE_PREFIX + "_sazi_angle.bin")

        dem_data = as_array(dem_path, dtype=np.float32)
        sol_zenith_angle_data = load_bin_file(sol_zenith_angle_path, TEST_REGION_NROW, TEST_REGION_NCOL, dtype=np.float32)
        sol_azimuth_angle_data = load_bin_file(sol_azimuth_angle_path, TEST_REGION_NROW, TEST_REGION_NCOL, dtype=np.float32)

        mask_all = run_castshadow(
            l1t,
            dem_data,
            sol_zenith_angle_data,
            sol_azimuth_angle_data,
            pb,
            block_height,
            block_width,
            is_utm)

        fuq_mask_all = load_bin_file(
            os.path.join(TEST_DATA_DIR, '%s_castshadow_s.bin' % REGION_DATE_PREFIX),
            TEST_REGION_NROW, TEST_REGION_NCOL,
            np.int16)

        compare_images(mask_all, fuq_mask_all, "cast shadow (solar)")





class RunSlopeAngleTestCase(TCTestCase):
    def setUp(self):
        super(RunSlopeAngleTestCase, self).baseSetUp("dsm_slope_angle_test")


    def test_random_with_files(self):
        """
        Test run_slop on some random data.
        """
        is_utm = True
        nrow = 50
        ncol = 75
        buffer_width = 11
        shape = ImageShape(origin_x=0., origin_y=0., dim_x=ncol, dim_y=nrow, pix_sz_x=0.00025)
        pb = Buffers(buffer_width)
        dem_data = np.random.ranf((nrow+pb.top+pb.bottom, ncol+pb.left+pb.right)).astype(np.float32)
        solar_zenith_data = np.random.ranf((nrow, ncol)).astype(np.float32)
        view_zenith_data = np.random.ranf((nrow, ncol)).astype(np.float32)
        solar_azimuth_data = np.random.ranf((nrow, ncol)).astype(np.float32)
        view_azimuth_data = np.random.ranf((nrow, ncol)).astype(np.float32)

        def filename(prefix):
            return os.path.join(self.output_dir, prefix + "_random.img")

        run_slope(
            shape,
            dem_data,
            solar_zenith_data,
            view_zenith_data,
            solar_azimuth_data,
            view_azimuth_data,
            pb,
            is_utm,
            #outputs
            self.output_format,
            filename('slope'),
            filename('aspect'),
            filename('incident'),
            filename('inci_azi'),
            filename('exiting'),
            filename('exit_azi'),
            filename('rela_angle'),
            filename('mask'))


    def test_random_without_files(self):
        """
        Test run_slop on some random data.
        """
        is_utm = True
        nrow = 50
        ncol = 75
        buffer_width = 11
        shape = ImageShape(origin_x=0., origin_y=0., dim_x=ncol, dim_y=nrow, pix_sz_x=0.00025)
        pb = Buffers(buffer_width)
        dem_data = np.random.ranf((nrow+pb.top+pb.bottom, ncol+pb.left+pb.right)).astype(np.float32)
        solar_zenith_data = np.random.ranf((nrow, ncol)).astype(np.float32)
        view_zenith_data = np.random.ranf((nrow, ncol)).astype(np.float32)
        solar_azimuth_data = np.random.ranf((nrow, ncol)).astype(np.float32)
        view_azimuth_data = np.random.ranf((nrow, ncol)).astype(np.float32)

        run_slope(
            shape,
            dem_data,
            solar_zenith_data,
            view_zenith_data,
            solar_azimuth_data,
            view_azimuth_data,
            pb,
            is_utm)


    def test_fuqin_data(self):
        slope_results = run_slope_fuq()

        def load_target_file(postfix, dtype=np.float32):
            return load_bin_file(
                os.path.join(TEST_DATA_DIR, '%s_%s.bin' % (REGION_DATE_PREFIX, postfix)),
                TEST_REGION_NROW, TEST_REGION_NCOL,
                dtype)

        compare_images(slope_results.slope, load_target_file('slope'), 'slope')
        compare_images(slope_results.aspect, load_target_file('aspect'), 'aspect')
        compare_images(slope_results.incident, load_target_file('incident'), 'incident')
        compare_images(slope_results.azi_incident, load_target_file('azi_incident'), 'azi_incident')
        compare_images(slope_results.exiting, load_target_file('exiting'), 'exiting')
        compare_images(slope_results.azi_exiting, load_target_file('azi_exiting'), 'azi_exiting')
        compare_images(slope_results.rela_slope, load_target_file('rela_slope'), 'rela_slope')
        compare_images(slope_results.mask_self, load_target_file('mask_self', np.int16), 'mask_self')





class RunTerrainCorrectionTestCase(TCTestCase):
    def setUp(self):
        super(RunTerrainCorrectionTestCase, self).baseSetUp("tc_test")


    def test_fuqin_data(self):
        start_end_file = os.path.join(TEST_DATA_DIR, 'startend_' + REGION_DATE_PREFIX + '_ortho.txt')
        rori = 0.52
        brdf0, brdf1, brdf2 = 0.058812, 0.023160, 0.010383
        bias, slope_ca, esun, dd = -1.52000, 0.762824, 1957.0, 0.98989
        ref_adj = 0.0365
        istart = np.zeros(TEST_REGION_NROW, dtype=np.int32)
        iend = np.zeros(TEST_REGION_NROW, dtype=np.int32)

        start_end = open(start_end_file)
        line_no = 0
        for line in start_end:
            bits = line.split()
            istart[line_no] = int(bits[1])
            iend[line_no] = int(bits[3])
            line_no = line_no + 1
        start_end.close()

        def load_file(rpath, dtype=np.float32, band=''):
            return load_bin_file(
                os.path.join(TEST_DATA_DIR, '%s_%s%s.bin' % (REGION_DATE_PREFIX, rpath, str(band))),
                TEST_REGION_NROW, TEST_REGION_NCOL,
                dtype)

        def load_target_file(fpref, band):
            return load_bin_file(
                os.path.join(TEST_DATA_DIR, '%s_%s_b%i.bin' % (fpref, REGION_DATE_PREFIX, band)),
                TEST_REGION_NROW, TEST_REGION_NCOL,
                np.int16)

        if RUN_WITH_NEW_DSM:
            slope_results = run_slope_fuq()
            def load_slope_file(rpath, dtype=np.float32, band=''):
                return slope_results.__getattribute__(rpath)
        else:
            def load_slope_file(rpath, dtype=np.float32, band=''):
                return load_file(rpath, dtype, band)

        gc.collect()

        for band in (1,2,3,4):
            dn_1_path = os.path.join(TEST_DATA_DIR, 'band1.dat')
            dn_1 = read_array_int8(dn_1_path, TEST_REGION_NROW, TEST_REGION_NCOL)

            ref_lm, ref_brdf, ref_terrain = run_brdfterrain(
                rori,
                brdf0, brdf1, brdf2,
                bias, slope_ca, esun, dd,
                ref_adj,
                istart, iend,
                dn_1,
                load_slope_file('mask_self', np.int16),
                load_file('castshadow_s', np.int16),
                load_file('castshadow_v', np.int16),
                load_file('solar_angle'),
                load_file('sazi_angle'),
                load_file('view_angle'),
                load_file('rela_angle'),
                load_slope_file('slope'),
                load_slope_file('aspect'),
                load_slope_file('incident'),
                load_slope_file('exiting'),
                load_slope_file('rela_slope'),
                load_file('a_band', band=band),
                load_file('b_band', band=band),
                load_file('s_band', band=band),
                load_file('fs_band', band=band),
                load_file('fv_band', band=band),
                load_file('ts_band', band=band),
                load_file('dir_band', band=band),
                load_file('dif_band', band=band))

            self.assertTupleEqual(dn_1.shape, ref_lm.shape, "ref_lm not of expected shape")
            self.assertTupleEqual(dn_1.shape, ref_brdf.shape, "ref_brdf not of expected shape")
            self.assertTupleEqual(dn_1.shape, ref_terrain.shape, "ref_terrain not of expected shape")

            # load Fuqin's results for compareson
            fuq_lm = load_target_file('ref_nobrdf', band)
            fuq_brdf = load_target_file('ref_wbrdf', band)
            fuq_terrain = load_target_file('ref_wterrain', band)

            compare_images(ref_lm, fuq_lm, "lambertial reflectance")
            compare_images(ref_brdf, fuq_brdf, "brdf")
            compare_images(ref_terrain, fuq_terrain, "terrain correted")

            dump_array(array=ref_terrain, output_path=TEST_TC_OUTPUT_PATH % band, file_format="ENVI")

            ref_lm = ref_brdf = ref_terrain = fuq_lm = fuq_brdf = fuq_terrain = None

            gc.collect()
