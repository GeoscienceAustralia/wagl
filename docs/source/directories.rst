Geoscience Australia Image Processor Directory Structure
========================

* bin
	* Contains both Fortran and C executables.

* src
        * Contains C code.

* fortran
        * Contains Fortran 77 code.
        * mock
                * Contains mocking routines for the Fortran executables.

* gaip
        * Contains Python modules and Fortran 90 code built as Python modules using F2Py.
        * tests
                * Contains unittesting modules.

* workflow
        * Contains the Luigi workflow and configuration script.

* docs
        * Documentation build scripts.
        * source
                * Contains documentations source trees.

* test
        * Sample shell scripts for job execution
