gaip dependencies
=================

At NCI, one can load all the required dependencies (for building and running) with the following commands (but not all of these may be necessary):

:command:`MODULEPATH=/projects/u46/opt/modules/modulefiles:$MODULEPATH`

:command:`module load intel-fc`

:command:`module load intel-cc`

:command:`module load intel-mkl/12.1.9.293`

:command:`module load python/2.7.6`

:command:`module load hdf4`

:command:`module load gdal/1.9.2`

:command:`module load proj/4.7.0`

:command:`module load luigi-mpi`

:command:`module load pyephem/3.7.5.1`

:command:`module load numexpr/2.2`

:command:`module load rasterio/0.15`

:command:`module load sphinx/1.3b2`

:command:`module load scikit-image/0.8.2`

:command:`module load IDL_functions`

:command:`module load EOtools`


Python packages
---------------

* `NumPy <http://www.numpy.org/>`_,
* `SciPy <http://www.scipy.org/>`_,
* `GDAL <https://pypi.python.org/pypi/GDAL/>`_,
* `Luigi <https://github.com/spotify/luigi/>`_,
* `numexpr <https://code.google.com/p/numexpr/>`_,
* `hdf4 <http://www.hdfgroup.org/products/hdf4/>`_,
* `proj <http://trac.osgeo.org/proj/>`_,
* `PyEphem <http://rhodesmill.org/pyephem/>`_,
* `scikit-image <http://scikit-image.org/>`_,
* `rasterio <https://github.com/mapbox/rasterio/>`_,
* `IDL_functions <https://github.com/sixy6e/IDL_functions/>`_,
* `EOtools <https://github.com/GeoscienceAustralia/EO_tools/>`_.


FORTRAN packages
----------------

* `MODTRAN <http://www.ontar.com/software/productdetails.aspx?item=modtran/>`_.


Ancillary
---------

* Aerosol; `AATSR <http://www.leos.le.ac.uk/aatsr/howto/index.html>`_,
* Ozone; In the form <month>.tif eg jan.tif,
* Water Vapour; In the form pr_wtr.eatm.<year>.tif eg pr_wtr.eatm.2009.tif. Each band in the image represents an 8 day period,
* Earth to Sun Distance in AU; A text file with a Julian day look-up,
* Solar Irradiance; Satellite-Sensor named text files of the form solar_irrad_<satellite-sensor>.txt eg solar_irrad_landsat8.txt,
* Satellite filter function; Satellite-spectral cover named text file of the form <satellite>_<spectral>.flt eg landsat8_vsir.flt,
* TLE (Two Line Element set) `TLE <http://en.wikipedia.org/wiki/Two-line_element_set>`_. Directory structure in the form /Satellite/TLE/Satellite_YEAR eg LANDSAT7/TLE/LS7_YEAR,
* BRDF; CSIRO Mosaics in the form 2005.02.02/MCD43A1.2005.033.aust.005.b01.500m_<wavelength>nm_brdf_par_f<factor>.hdf.gz eg 2005.02.02/MCD43A1.2005.033.aust.005.b01.500m_0620_0670nm_brdf_par_fiso.hdf.gz,
* Pre-MODIS BRDF can be obtained by producing decadal averages for a given 8 day peroid. The directory structure is a Julian day eg 001, 009, 017 for the first 2 8 day periods,
* DSM; A National scale Digital Surface Model, or at least one covering the area of interest,


Documentation
-------------

* `Sphinx <http://sphinx-doc.org/>`_ (at NCI, you can load this with :command:`module load sphinx`), and
* `LaTeX <http://www.latex-project.org/>`_ (if you want to build PDF documentation).


Building
--------

* `F2py <http://www.scipy.org/F2py>`_ (which should be included with `SciPy <http://www.scipy.org/>`_), and
* `GCC <http://gcc.gnu.org/>`_ or some other suitable suite of (C and Fortran) compilers.
