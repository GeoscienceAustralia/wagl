wagl dependencies
=================


Software requirements
----------------------

* `NumPy <http://www.numpy.org/>`_,
* `SciPy <http://www.scipy.org/>`_,
* `F2py <http://www.scipy.org/F2py>`_ (which should be included with `SciPy <http://www.scipy.org/>`_)
* `GDAL <https://pypi.python.org/pypi/GDAL/>`_,
* `Luigi <https://github.com/spotify/luigi/>`_,
* `numexpr <https://github.com/pydata/numexpr>`_,
* `hdf5 <https://support.hdfgroup.org/HDF5/>`_,
* `hdf4 <http://www.hdfgroup.org/products/hdf4/>`_,
* `proj <http://trac.osgeo.org/proj/>`_,
* `PyEphem <http://rhodesmill.org/pyephem/>`_,
* `scikit-image <http://scikit-image.org/>`_,
* `rasterio <https://github.com/mapbox/rasterio/>`_,
* `h5py <https://github.com/h5py/h5py>`_,
* `idl-functions <https://github.com/sixy6e/idl-functions>`_,
* `bitshuffle <https://github.com/kiyo-masui/bitshuffle>`_,
* `mafisc <https://wr.informatik.uni-hamburg.de/research/projects/icomex/mafisc>`_,
* `pytables <https://github.com/PyTables/PyTables>`_,
* `pandas <https://github.com/pandas-dev/pandas>`_,
* `geopandas <https://github.com/geopandas/geopandas>`_
* `attr <https://github.com/python-attrs/attrs>`_


Radiative Transfer packages
---------------------------

* `MODTRAN <http://modtran.spectral.com/>`_.


Ancillary
---------

* Aerosol; `AATSR <http://www.leos.le.ac.uk/aatsr/howto/index.html>`_,
* Ozone; In the form <month>.tif eg jan.tif,
* Water Vapour; In the form pr_wtr.eatm.<year>.tif eg pr_wtr.eatm.2009.tif. Each band in the image represents an 8 day period,
* Earth to Sun Distance in AU; A text file with a Julian day look-up,
* Solar Irradiance; Satellite-Sensor named text files of the form solar_irrad_<satellite-sensor>.txt eg solar_irrad_landsat8.txt,
* TLE (Two Line Element set) `TLE <http://en.wikipedia.org/wiki/Two-line_element_set>`_. Directory structure in the form /Satellite/TLE/Satellite_YEAR eg LANDSAT7/TLE/LS7_YEAR,
* BRDF; CSIRO Mosaics in the form 2005.02.02/MCD43A1.2005.033.aust.005.b01.500m_<wavelength>nm_brdf_par_f<factor>.hdf.gz eg 2005.02.02/MCD43A1.2005.033.aust.005.b01.500m_0620_0670nm_brdf_par_fiso.hdf.gz,
* Pre-MODIS BRDF can be obtained by producing decadal averages for a given 8 day peroid. The directory structure is a Julian day eg 001, 009, 017 for the first 2 8 day periods,
* DSM; A National scale Digital Surface Model, or at least one covering the area of interest,


Documentation
-------------

* `Sphinx <http://sphinx-doc.org/>`_
* `LaTeX <http://www.latex-project.org/>`_ (if you want to build PDF documentation).


Building
--------

* python setup.py install --prefix=<prefix>
