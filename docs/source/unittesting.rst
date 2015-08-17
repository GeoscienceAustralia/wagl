
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

* `test_calculate_angles.py` Which tests the main intermediate files that get produced during the angular computational phase.
* `test_cast_shadow_satellite.py` Test the cast shadow array generated from the satellite look angles.
* `test_cast_shadow_sun.py` Test the cast shadow array generated from the solar look angles.
* `test_smoothe_dsm.py` Test the sobel smoothing result of the DEM subset.
* `test_nbar_tc.py` Test the lambertian, NBAR-BRDF, and terrain corrected NBAR-BRDF output products, for each band.


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


* `python test_nbar_tc.py --help` outputs:

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
