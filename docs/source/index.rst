.. ULA3 documentation master file, created by
   sphinx-quickstart on Mon Apr 15 11:54:15 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to ULA3's documentation!
================================

Contents:

.. toctree::
   :maxdepth: 2

   invoking
   dependencies
   directories
   history
   algorithms
   tc_implementation
   testing

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Bulding ULA3
============

ULA3 uses several Fortran programs, C programs and modules built using `F2py <http://www.scipy.org/F2py>`_. The build process is performed using make, though if make is not available, it is pretty clear what needs to be done from the make files (perhaps we should provide a batch file for windows?). Once the required bits and pieces are installed (see :doc:`dependencies </dependencies>`), building and generating documentation should be as easy as typing ``make`` in the top level directory.

ULA3 testing
============

Integration testing can be performed on VAYU using the scripts found under :file:`test` in the top level folder of the projects code. These scripts are designed to be submitted as jobs via the command :command:`qsub` (they will not run directly as the memory requirements are to large).

Unit tests written using :py:mod:`unittest` for some modules can be found in :py:mod:`ULA3.tests`. Some of these tests also depend on data on VAYU and will not run on other systems. Please see the :py:mod:`unittest` documentation on how to run these.
