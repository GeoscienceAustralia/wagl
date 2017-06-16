gaip Directory Structure
========================

* gaip
        * Contains Python modules and FORTRAN 90 code built as Python modules using F2Py.
        * f90_sources
                * Contains the FORTRAN 90 code to be built as Python modules using F2Py.
        * spectral_response
                * Conatains the spectral responses for various satellites and sensors.
        * tests
                * data
                        Contains small data files required for some unittests.
                * Contains unittesting modules.
        * scripts
                * Contains the source code for command line scripts and utilities.

* docs
        * Documentation build scripts.
        * source
                * Contains documentations source trees.

* configs
        * Sample configuration files for overriding luigi parameters, and logging.

* utils
        * Commandline scripts and utilities
