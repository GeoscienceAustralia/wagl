wagl Directory Structure
========================

* wagl
        * Contains Python modules and FORTRAN 90 code built as Python modules using F2Py.
        * acquisition
                * Contains the acquisition base class and information regarding the satellite/sensors
                  supported by wagl.
        * f90_sources
                * Contains the FORTRAN 90 code to be built as Python modules using F2Py.
        * spectral_response
                * Conatains the spectral responses for various satellites and sensors.
        * scripts
                * Contains the source code for command line scripts and utilities.

* docs
        * Documentation build scripts.
        * source
                * Contains documentations source trees.

* configs
        * Sample configuration files for overriding luigi parameters, and logging.

* tests
        * data
                Contains small data files required for some unittests.
        * Contains unittesting modules.

* utils
        * Commandline scripts and utilities
