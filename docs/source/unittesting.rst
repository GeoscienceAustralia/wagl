
=================
Running unittests
=================

Purpose
-------
The purpose behind the unittesting environment is to ensure consistancy of the product outputs
between different versions of the code.

Semi-automated tests have been written for the NBAR-TC algorithm, for not just the final product
outputs, but also for many of the intermediate files that get generated during the computational
phase.

NBAR-TC
-------
The following tests are available for the NBAR-TC algorithm:

* `test_calculate_angles.py` Which tests the main intermediate files that get produced during the angular computational phase (Option to compute).
* `test_cast_shadow_satellite.py` Test the cast shadow array generated from the satellite look angles (Option to compute).
* `test_cast_shadow_sun.py` Test the cast shadow array generated from the solar look angles (Option to compute).
* `test_smoothe_dsm.py` Test the sobel smoothing result of the DEM subset (Option to compute).
* `test_exiting_angles.py` Test the exiting and azimuthal exiting angle arrays (Option to compute).
* `test_incident_angles.py` Test the incident and azimuthal incident angle arrays (Option to compute).
* `test_relative_slope.py` Test the relative azimithal angle on the slope surface array (Option to compute).
* `test_slope_aspect.py` Test the slope and aspect arrays (Option to compute).
* `test_self_shadow_mask.py` Test the self shadow mask array (Option to compute).
* `test_nbar_tc.py` Test the lambertian, NBAR-BRDF, and terrain corrected NBAR-BRDF output products, for each band.

For those tests with an option to compute, the required outputs will be computed using pre-existing inputs from the reference directory before
any tests are evaluated against the reference data.


Fractional Cover
----------------
The following test was written for the products that are output as part of the agdc fractional cover code base:

* `test_fc_products.py` Tests each fractional cover component contained within the output product.


Running each script
-------------------
Each script has in-built help available. To see the available help execute:
`python {script} --help` from the command line. For example:

* `python test_nbar_tc.py --help` which outputs:


usage: test_nbar_tc.py [-h] --reference_dir REFERENCE_DIR --test_dir TEST_DIR
                       [--int_precision INT_PRECISION]

Test and compare the Lambertian, BRDF and TC outputs.

optional arguments:
  -h, --help            show this help message and exit
  --reference_dir REFERENCE_DIR
                        A directory path of a nbar reference output.
  --test_dir TEST_DIR   A directory path to the test output.
  --int_precision INT_PRECISION
                        The integer precision used for array comparison


* `python test_calculate_angles.py --help` outputs:

usage: test_calculate_angles.py [-h] --L1T_dir L1T_DIR --nbar_work_dir
                                NBAR_WORK_DIR --outdir OUTDIR
                                [--dec_precision DEC_PRECISION]
                                [--int_precision INT_PRECISION] [--compute]

Calculates satellite and solar angles.

optional arguments:
  -h, --help            show this help message and exit
  --L1T_dir L1T_DIR     A directory path of a L1T scene.
  --nbar_work_dir NBAR_WORK_DIR
                        A directory path to the associated NBAR working
                        directory.
  --outdir OUTDIR       A directory path that will contain the output files.
  --dec_precision DEC_PRECISION
                        The decimal precision used for array comparison
  --int_precision INT_PRECISION
                        The integer precision used for array comparison
  --compute             If set then the solar and sateliite angles will be
                        computed as will the CENTRELINE text file before
                        running the unittests.

* `python test_cast_shadow_satellite.py --help` outputs:

usage: test_cast_shadow_satellite.py [-h] --L1T_dir L1T_DIR --nbar_work_dir
                                     NBAR_WORK_DIR --outdir OUTDIR
                                     [--dec_precision DEC_PRECISION]
                                     [--int_precision INT_PRECISION]
                                     [--compute] [--buffer BUFFER]
                                     [--block_x BLOCK_X] [--block_y BLOCK_Y]

Perform unittesting for the satellite view shadow mask and optionaly
calculates the satellite view shadow mask.

optional arguments:
  -h, --help            show this help message and exit
  --L1T_dir L1T_DIR     A directory path of a L1T scene.
  --nbar_work_dir NBAR_WORK_DIR
                        A directory path to the associated NBAR working
                        directory.
  --outdir OUTDIR       A directory path that will contain the output files.
  --dec_precision DEC_PRECISION
                        The decimal precision used for array comparison
  --int_precision INT_PRECISION
                        The integer precision used for array comparison
  --compute             If set then the satellite view shadow mask will be
                        computed before running the unittests.
  --buffer BUFFER       The buffer in pixels to be used in calculating the
                        satellite view shadow.
  --block_x BLOCK_X     The x block size in pixels (Twice the buffer).
  --block_y BLOCK_Y     The y block size in pixels (Twice the buffer)..

* `python test_cast_shadow_sun.py --help` outputs:

usage: test_cast_shadow_sun.py [-h] --L1T_dir L1T_DIR --nbar_work_dir
                               NBAR_WORK_DIR --outdir OUTDIR
                               [--dec_precision DEC_PRECISION]
                               [--int_precision INT_PRECISION] [--compute]
                               [--buffer BUFFER] [--block_x BLOCK_X]
                               [--block_y BLOCK_Y]

Perform unittesting for the cast shadow sun mask and optionaly calculates the
cast shadow sun mask.

optional arguments:
  -h, --help            show this help message and exit
  --L1T_dir L1T_DIR     A directory path of a L1T scene.
  --nbar_work_dir NBAR_WORK_DIR
                        A directory path to the associated NBAR working
                        directory.
  --outdir OUTDIR       A directory path that will contain the output files.
  --dec_precision DEC_PRECISION
                        The decimal precision used for array comparison
  --int_precision INT_PRECISION
                        The integer precision used for array comparison
  --compute             If set then the self shadow array will be computed
                        before running the unittests.
  --buffer BUFFER       The buffer in pixels to be used in calculating the
                        cast shadow sun mask.
  --block_x BLOCK_X     The x block size in pixels (Twice the buffer).
  --block_y BLOCK_Y     The y block size in pixels (Twice the buffer)..

The unittests for `test_smoothe_dsm.py`, `test_exiting_angles.py`, `test_incident_angles.py`,
`test_relative_slope.py`, `test_slope_aspect.py`, and `test_self_shadow_mask.py` all use the same commandline arguments:

optional arguments:
  -h, --help            show this help message and exit
  --L1T_dir L1T_DIR     A directory path of a L1T scene.
  --nbar_work_dir NBAR_WORK_DIR
                        A directory path to the associated NBAR working
                        directory.
  --outdir OUTDIR       A directory path that will contain the output files.
  --dec_precision DEC_PRECISION
                        The decimal precision used for array comparison
  --int_precision INT_PRECISION
                        The integer precision used for array comparison
  --compute             If set then the self shadow array will be computed
                        before running the unittests.

* `python test_fc_products.py --help` outputs:

usage: test_fc_products.py [-h] --reference_filename REFERENCE_FILENAME
                           --test_filename TEST_FILENAME
                           [--int_precision INT_PRECISION]

Test and compare the fractional cover outputs.

optional arguments:
  -h, --help            show this help message and exit
  --reference_filename REFERENCE_FILENAME
                        The filename for the referenc dataset.
  --test_filename TEST_FILENAME
                        The filename for the test dataset.
  --int_precision INT_PRECISION
                        The integer precision used for array comparison

Each script more or less follows the same principle. i.e. provide the reference directory and the test directory.
The fractional cover is slightly different, the user is expected to point to an actual file rather than the directory.
