ULA3 dependencies
=================

On VAYU, one can load all the required dependencies (for building and running) with the following commands (but not all of these may be necessary):

:command:`MODULEPATH=/g/data/v10/opt/modules/modulefiles:$MODULEPATH`

:command:`module load intel-fc`

:command:`module load intel-cc`

:command:`module load intel-mkl/10.3.0`

:command:`module load geotiff/1.3.0`

:command:`module load python/2.7.2`

:command:`module load hdf4/4.2.6_2012`

:command:`module load gdal/1.9.0_HDF5`

:command:`module load proj/4.7.0`

:command:`module load pyephem/3.7.5.1`

:command:`module load numexpr/2.0.1`

:command:`module load geotiff/1.3.0`

:command:`module load sphinx/1.1.3`


Python packages
---------------

* `NumPy <http://www.numpy.org/>`_,
* `SciPy <http://www.scipy.org/>`_,
* `GDAL <https://pypi.python.org/pypi/GDAL/>`_,
* `numexpr <https://code.google.com/p/numexpr/>`_,
* hdf4,
* proj,
* `PyEphem <http://rhodesmill.org/pyephem/>`_, and
* osr.


Documentation
-------------

* `Sphinx <http://sphinx-doc.org/>`_ (on VAYU, you can load this with :command:`module load sphinx`), and
* `LaTeX <http://www.latex-project.org/>`_ (if you want to build PDF documentation).


Building
--------

* `F2py <http://www.scipy.org/F2py>`_ (which should be included with `SciPy <http://www.scipy.org/>`_), and
* `GCC <http://gcc.gnu.org/>`_ or some other suitable suite of (C and Fortran) compilers.
