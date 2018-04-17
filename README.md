# wagl
------


wagl is a Python package for producing standarised imagery in the form of:

* Nadir Bi-directional Reflectance Distribution Function Adjusted Reflectance (NBAR)
* NBART; NBAR with Terrain Illumination correction
* Surface Brightness Temperature
* Pixel Quality (per pixel metadata)

The luigi task workflow for producing NBAR for a Landsat 5TM scene is given below.

![](docs/source/diagrams/luigi-task-visualiser-reduced.png)

## Supported Satellites and Sensors
-----------------------------------
* Landsat 5 TM
* Landsat 7 ETM
* Landsat 8 OLI
* Landsat 8 TIRS
* Sentinel-2a

## Requirements
---------------
* [luigi](https://github.com/spotify/luigi)
* [numpy](https://github.com/numpy/numpy)
* [scipy](https://github.com/scipy/scipy)
* [numexpr](https://github.com/pydata/numexpr)
* [pyephem](http://rhodesmill.org/pyephem/)
* [proj](https://github.com/OSGeo/proj.4)
* [h5py](https://github.com/h5py/h5py)
* [tables](https://github.com/PyTables/PyTables)
* [pandas](https://github.com/pandas-dev/pandas)
* [scikit-image](https://github.com/scikit-image/scikit-image)
* [GDAL](https://github.com/OSGeo/gdal)
* [rasterio](https://github.com/mapbox/rasterio)
* [fiona](https://github.com/Toblerity/Fiona)
* [shapely](https://github.com/Toblerity/Shapely)
* [geopandas](https://github.com/geopandas/geopandas)
* [pyyaml](https://github.com/yaml/pyyaml)
* [attr](https://github.com/python-attrs/attrs)

## Installation
---------------

### wagl Package
The wagl pacakage can be installed via:

`$ python setup.py install --prefix=<prefix>`

### Additional HDF5 compression filters (optional)
Additional compression filters can be used via HDF5's
[dynamically loaded filters](https://support.hdfgroup.org/HDF5/doc/Advanced/DynamicallyLoadedFilters/HDF5DynamicallyLoadedFilters.pdf).
Essentially the filter needs to be compiled against the HDF5 library, and
installed into HDF5's plugin path, or a path of your choosing, and set the
HDF5_PLUGIN_PATH environment variable. The filters are then automatically
accessible by HDF5 via the [integer code](https://support.hdfgroup.org/services/contributions.html)
assigned to the filter.

#### Mafisc compression filter
Mafisc combines both a bitshuffling filter and lzma compression filter in order
to get the best compression possible at the cost of lower compression speeds.
To install the `mafisc` compression filter, follow these [instructions](https://wr.informatik.uni-hamburg.de/research/projects/icomex/mafisc).

#### Bitshuffle
The [bitshuffle filter](https://github.com/kiyo-masui/bitshuffle) can be installed
from source, or conda via the supplied [conda recipe](https://github.com/kiyo-masui/bitshuffle/tree/master/conda-recipe).
It utilises a bitshuffling filter on top of either a lz4 or lzf compression filter.

## Basic command line useage
--------------------------
Using the [local scheduler](http://luigi.readthedocs.io/en/stable/command_line.html):

    $ luigi --module wagl.multifile_workflow ARD --model NBAR --level1-list scenes.txt --outdir /some/path --local-scheduler --workers 4

Using the [central scheduler](http://luigi.readthedocs.io/en/stable/central_scheduler.html):

    $ luigid --background --pidfile <PATH_TO_PIDFILE> --logdir <PATH_TO_LOGDIR> --state-path <PATH_TO_STATEFILE>

    $ luigi --module wagl.multifile_workflow ARD --level1-list scenes.txt --model STANDARD --outdir /some/path --workers 4

    $ luigi --module wagl.multifile_workflow ARD --level1-list scenes.txt --model NBAR --outdir /some/path --workers 4

    $ luigi --module wagl.multifile_workflow ARD --level1-list scenes.txt --model SBT --outdir /some/path --workers 4
