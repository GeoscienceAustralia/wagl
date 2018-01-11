HDF5
====

As of wagl v5.0.0, the storage backend for all dataset outputs uses HDF5, whether the datasets be imagery, tables, or scalars.
Currently datasets can be stored according to the *IMAGE* Class, or *TABLE* Class.

* **IMAGE** Class `image specification <https://support.hdfgroup.org/HDF5/doc/ADGuide/ImageSpec.html>`_
* **TABLE** Class `table specification <https://support.hdfgroup.org/HDF5/doc/HL/H5TB_Spec.html>`_

Moving to HDF5
--------------

Some of the main reasons behind moving to HDF5 are:

* Reduce the number of files being written to disk. Yes imagery could be stored multi-band, but band ordering would need to be explicity defined, and we have the added benefit of including metadata such as Dataset descriptions, thereby allowing an easier knowledge transfer for future developers. For some variations on processing a scene, the number of files output to disk was reduced by approximately 22% when using the *wagl.multifile_workflow*. When using the *wagl.singlefile_workflow* the number of files output is reduced to **1** (the MODTRAN working directory is automatically cleaned up after the results are saved into the HDF5 file).
* Consequently, a reduction in the number of files, also reduces the number of parameters being passed back and forth between `luigi <https://github.com/spotify/luigi>`_ *Tasks*.
* Consolodate the numerous text files that could be comma, tab, space, and variable space delimited variations, into a consistant *TABLE* structure that allows compression.
* *TABLE* and *IMAGE* datasets can be stored and accessed together from the same file.
* HDF5 has an `ExternalLink <http://docs.h5py.org/en/latest/high/group.html#group-extlinks>`_ feature, which for the case of the multitude of interpolated factors for each spectral band, all the outputs can be combined and accessed from the same file, thereby simplifying the workflow and functional architecture of wagl.
* Various `compression <https://support.hdfgroup.org/services/contributions.html>`_ filters are available, ranging from extremely high I/O speed, to slower I/O but very high compression ratios.
  High compression ratios allow for better archiving, and with conjuction of being able to store *TABLES*, *IMAGES*, *SCALARS*, variable resolutions, all within the same file, can aid in keeping the number of inodes down on HPC systems.
  This presents an archiving solution without the need to create tarballs, whilst still maintaining the capacity to operate with the archives. Essentially creating a kind of operational archive.
* Datasets can be stored in a hierarchical fashion, for instance, creating a hierarchical layout based on resolution, or the same dataset name stored in different locations based on a given attribute.
  For example, the ancillary data required for surface brightness temperature, gathers data at a minimum of 25 points across a scene or granule.
  The point id is used as the group label to differentiate between the same ancillary data gathered at different locations. i.e.:

    * /POINT-0/TEMPERATURE
    * /POINT-1/TEMPERATURE
    * /POINT-2/TEMPERATURE

* Metadata; wagl can now store a lot more metadata such as longitude and latitude information with each ancillary point location, as opposed to a plain text label.
  Parameter settings used for a given algorithm such as for the satellite and solar angles calculation can be stored alongside the results, potentially making it easier for validation, and archive comparisons to be undertaken. Dataset descriptions have been useful for new people working with the code base.
* Utilise a consistant library for data I/O, rather than a dozen or so different libraries. This helps to simplify the wagl data model.
* It simplified the workflow **A LOT**. By writing multiple datasets within the same file, the parameter passing bewteen Task's and functions, was reduced in some cases from a dozen parameters, to a single parameter. Some Tasks would act as a helper task whose sole purpose was to manage a bunch of individual Tasks, and combine the results into a single file for easy access by downstream Tasks which then didn't have to open several dozen files. This is achieved by HDF5's ExternalLink feature, which allows the linking of Datasets contained within other files to be readable from a single file that acts as a Table Of Contents (TOC). It is similar to a UNIX symbolic link, or Windows shortcut.
* A simpler structure for testing and evaluation; eg compare the same scene but different versions of wagl. And have the results stored directly alongside the inputs. That way, it is easier to track exactly what datasets were used to determine the comparison, and have it immediately in a form suitable for archiving and immediate access without having to decompress a tarball containing hundreds or thousands of scenes.

Singlefile workflow
-------------------

The singlefile workflow is useful in situations where you don't want several hundred files being output per scene, which can clog the filesystem, or if submitting a large list of scenes for processing and rather than clogging the scheduler with hundreds of tasks per scene, there'll be a single task per scene.
This approach is useful for large scale production, and the results are stored in a fashion that serves as both an archive and operational use. One could call it an *operational archive*.
Routine testing and comparisons between different versions of *wagl* (assuming that the base naming structure is the same), is also more easily evaluated.
However, the *wagl.singlefile_workflow* isn't able to pickup from where the workflow left off, as is the case in the *wagl.multifile_workflow*, instead it starts off from scratch. The retry count is set to 1 for *wagl.singlefile_workflow*, rather than the default of 3. This allows luigi to attempt to reprocesss the scene one more time, before flagging it as an error.

The contents for a Landsat 8 scene going through the nbar model and (3, 3) vertices for the radiative transfer is list in :ref:`appendix_a`.

Dataset Names
-------------

Dataset names for each output are as follows:

* ancillary.h5
    * /ANCILLARY/AEROSOL
    * /ANCILLARY/COORDINATOR
    * /ANCILLARY/ELEVATION
    * /ANCILLARY/OZONE
    * /ANCILLARY/WATER-VAPOUR
    * /ANCILLARY/BRDF-{band_name}-GEO *(combinations based on band_name)*
    * /ANCILLARY/BRDF-{band_name}-ISO *(combinations based on band_name)*
    * /ANCILLARY/BRDF-{band_name}-VOL *(combinations based on band_name)*
* atmospheric-inputs.h5
    * /ATMOSPHERIC-INPUTS/POINT-0/ALBEDO-0/TP5_DATA *(combinations based on point label, albedo label)*
    * /ATMOSPHERIC-INPUTS/POINT-1/ALBEDO-1/TP5_DATA *(combinations based on point label, albedo label)*
    * /ATMOSPHERIC-INPUTS/POINT-2/ALBEDO-T/TP5_DATA *(combinations based on point label, albedo label)*
    * /ATMOSPHERIC-INPUTS/POINT-3/ALBEDO-TH/TP5_DATA *(combinations based on point label, albedo label)*
* atmospheric-results.h5
    * /ATMOSPHERIC-RESULTS/POINT-0/ALBEDO-0/FLUX *(combinations based on point label, albedo label)*
    * /ATMOSPHERIC-RESULTS/POINT-0/ALBEDO-0/ALTITUDES *(combinations based on point label, albedo label)*
    * /ATMOSPHERIC-RESULTS/POINT-0/ALBEDO-0/CHANNEL *(combinations based on point label, albedo label)*
    * /ATMOSPHERIC-RESULTS/POINT-0/ALBEDO-0/SOLAR-IRRADIANCE *(combinations based on point label, albedo label)*
    * /ATMOSPHERIC-RESULTS/POINT-0/ALBEDO-TH/UPWARD-RADIATION-CHANNEL **SBT Only** *(combinations based on point label)*
    * /ATMOSPHERIC-RESULTS/POINT-0/ALBEDO-TH/DOWNWARD-RADIATION-CHANNEL **SBT Only** *(combinations based on point label)*
* components.h5
    * /ATMOSPHERIC-COMPONENTS/NBAR-COMPONENTS
    * /ATMOSPHERIC-COMPONENTS/SBT-COMPONENTS
* interpolated-components.h5
    * /INTERPOLATED-ATMOSPHERIC-COMPONENTS/A/{band_name} *(combinations based on the band_name)*
    * /INTERPOLATED-ATMOSPHERIC-COMPONENTS/B/{band_name} *(combinations based on the band_name)*
    * /INTERPOLATED-ATMOSPHERIC-COMPONENTS/DIF/{band_name} *(combinations based on the band_name)*
    * /INTERPOLATED-ATMOSPHERIC-COMPONENTS/DIR/{band_name} *(combinations based on the band_name)*
    * /INTERPOLATED-ATMOSPHERIC-COMPONENTS/FS/{band_name} *(combinations based on the band_name)*
    * /INTERPOLATED-ATMOSPHERIC-COMPONENTS/FV/{band_name} *(combinations based on the band_name)*
    * /INTERPOLATED-ATMOSPHERIC-COMPONENTS/S/{band_name} *(combinations based on the band_name)*
    * /INTERPOLATED-ATMOSPHERIC-COMPONENTS/TS/{band_name} *(combinations based on the band_name)*
* dsm-extract.h5
    * /ELEVATION/DSM
    * /ELEVATION/DSM-SMOOTHED
* exiting-angles.h5
    * /EXITING-ANGLES/AZIMUTHAL-EXITING
    * /EXITING-ANGLES/EXITING
* incident-angles.h5
    * /INCIDENT-ANGLES/AZIMUTHAL-INCIDENT
    * /INCIDENT-ANGLES/INCIDENT
* longitude-latitude.h5
    * /LONGITUDE-LATITUDE/LONGITUDE
    * /LONGITUDE-LATITUDE/LATITUDE
* relative-slope.h5
   * /RELATIVE-SLOPE/RELATIVE-SLOPE
* satellite-solar.h5
    * /SATELLITE-SOLAR/BOXLINE
    * /SATELLITE-SOLAR/CENTRELINE
    * /SATELLITE-SOLAR/PARAMETERS/ORBITAL-ELEMENTS
    * /SATELLITE-SOLAR/PARAMETERS/SATELLITE-MODEL
    * /SATELLITE-SOLAR/PARAMETERS/SATELLITE-TRACK
    * /SATELLITE-SOLAR/PARAMETERS/SPHEROID
    * /SATELLITE-SOLAR/RELATIVE-AZIMUTH
    * /SATELLITE-SOLAR/SATELLITE-AZIMUTH
    * /SATELLITE-SOLAR/SATELLITE-VIEW
    * /SATELLITE-SOLAR/SOLAR-AZIMUTH
    * /SATELLITE-SOLAR/SOLAR-ZENITH
* shadow-masks.h5
    * /SHADOW-MASKS/CAST-SHADOW-SATELLITE
    * /SHADOW-MASKS/CAST-SHADOW-SUN
    * /SHADOW-MASKS/COMBINED-SHADOW
    * /SHADOW-MASKS/SELF-SHADOW
* slope-aspect.h5
    * /SLOPE-ASPECT/ASPECT
    * /SLOPE-ASPECT/SLOPE
* standardised-products.h5
    * /METADATA/NBAR-METADATA
    * /METADATA/PQ-METADATA
    * /METADATA/SBT-METADATA
    * /STANDARDISED-PRODUCTS/REFLECTANCE/LAMBERTIAN/{band_name} *(combinations based on the band_name)*
    * /STANDARDISED-PRODUCTS/REFLECTANCE/NBAR/{band_name} *(combinations based on the band_name)*
    * /STANDARDISED-PRODUCTS/REFLECTANCE/NBART/{band_name} *(combinations based on the band_name)*
    * /STANDARDISED-PRODUCTS/PIXEL-QUALITY/PIXEL-QUALITY
    * /STANDARDISED-PRODUCTS/THERMAL/SURFACE-BRIGHTNESS-TEMPERATURE/{band_name} *(combinations based on the band_name)*

Geospatial Information
----------------------

Geospatial information for *IMAGE* Class datasets can be stored in various different ways. For wagl, we attach 2 attributes specifically related to geospatial context:

* transform (GDAL like GeoTransform; 6 element array)
* crs_wkt (CRS stored as a variable length string using the Well Known Text specification

This approach is very simple, and similar to lots of other mainstream formats such as `ENVI <https://www.harrisgeospatial.com/docs/ENVIHeaderFiles.html>`_,
`KEA <http://kealib.org/>`_. The geospatial information can automatically be interpreted using *wagl.geobox.GriddedGeoBox*.

Tables
------

Tabulated data created by wagl is stored in HDF5 using the compound datatype, and read back into memory as either a custom *NumPy* datatype, or directly into a *pandas.DataFrame*.
Datatypes are mapped between HDF5 and NumPy as best as possible. Additional attached attributes inlcuded by wagl can aid in the transitional mapping.
`PyTables <http://www.pytables.org/>`_ could've been used to store the tables, as well as the imagery, however `h5py <http://www.h5py.org/>`_ provides a simpler api, as well as optional mpi driver mode for when the case arises (HDF5 must be compiled with the MPI switch turned on).

An example table is the *coordinator* table used to define the point locations at which to run the atmospheric calculations.

+-----------+--------------+------------+------------+---------+--------+
| row_index | column_index | latitude   | longitude  | map_y   | map_x  |
|           |              |            |            |         |        |
+===========+==============+============+============+=========+========+
|    0      | 1395         | -33.636477 | 147.233989 | 6278125 | 521700 |
+-----------+--------------+------------+------------+---------+--------+
|    0      | 4299         | -33.632518 | 148.016761 | 6278125 | 594300 |
+-----------+--------------+------------+------------+---------+--------+
|    0      | 9729         | -33.611835 | 149.479600 | 6278125 | 730050 |
+-----------+--------------+------------+------------+---------+--------+
| 4299      |  339         | -34.605977 | 146.948739 | 6170650 | 495300 |
+-----------+--------------+------------+------------+---------+--------+
| 4299      | 4299         | -34.601653 | 148.028427 | 6170650 | 594300 |
+-----------+--------------+------------+------------+---------+--------+
| 4299      | 9395         | -34.582043 | 149.417061 | 6170650 | 721700 |
+-----------+--------------+------------+------------+---------+--------+
| 8598      |    0         | -35.575035 | 146.854595 | 6063175 | 486825 |
+-----------+--------------+------------+------------+---------+--------+
| 8598      | 4299         | -35.570630 | 148.040664 | 6063175 | 594300 |
+-----------+--------------+------------+------------+---------+--------+
| 8598      | 8337         | -35.555872 | 149.154192 | 6063175 | 695250 |
+-----------+--------------+------------+------------+---------+--------+

An example of how to read the coordinator table into a *pandas.DataFrame*:

       .. code-block:: python

          >>> from wagl.hdf5 import read_h5_table
          >>> import h5py
          >>> fid = h5py.File('coordinator.h5', 'r')
          >>> df = read_h5_table(fid, 'nbar-coordinator')


Attributes (metadata)
---------------------

All datasets created by *wagl* have attributes attached to them. Each dataset class type eg *SCALAR*, *TABLE*, *IMAGE*, has its own unique attribute set, as well as some common attribute labels.
The attributes can be printed to screen using the *wagl_ls --filename my-file.h5 --verbose* utility script, or the *wagl.hdf5.h5ls* function and setting the *verbose=True* parameter. Additionally one can also use HDF5's h5ls command line utility which *wagl's* version is fashioned afer.
The attributes can also be extracted and written to disk using the `yaml <https://en.wikipedia.org/wiki/YAML>`_ format, using the *wagl_convert* utility script, which converts Images to GeoTiff, Tables to csv, and Scalars to yaml.


TABLE attributes
~~~~~~~~~~~~~~~~

TABLE datasets will be written following the HDF5 `table specification <https://support.hdfgroup.org/HDF5/doc/HL/H5TB_Spec.html>`_
with just the base amount of information such as:

* CLASS
* VERSION
* TITLE
* field/column names

Most table datasets sourced from a NumPy `structured array <https://docs.scipy.org/doc/numpy/user/basics.rec.html>`_
will be of this simpler form, and might have an additional attribute such as *description*.

If the source of the table was a pandas `DataFrame <https://pandas.pydata.org/pandas-docs/stable/generated/pandas.DataFrame.html>`_,
then additional attributes will be attached such as:

* datatype mappings between HDF5 and pandas
* number of row the table contains
* column(s) to be used as the index

An example of the attributes attached to a table dataset whose source is a Numpy structured array is given below as a *yaml* document which is what would be yielded if using the command line utility *wagl_convert*:

.. code-block:: yaml

   CLASS: TABLE
   description: Contains the array, latitude and longitude coordinates of the satellite
       track path.
   FIELD_0_NAME: row_index
   FIELD_1_NAME: col_index
   FIELD_2_NAME: n_pixels
   FIELD_3_NAME: latitude
   FIELD_4_NAME: longitude
   TITLE: Centreline
   VERSION: '0.2'
   array_coordinate_offset: 0

An example of the attributes attached to a table dataset whose source is a pandas Dataframe, once again as a yaml document, is given below:

.. code-block:: yaml

   Albedo: '0'
   CLASS: TABLE
   description: Accumulated solar irradiation for point 0 and albedo 0.
   FIELD_0_NAME: index
   FIELD_1_NAME: diffuse
   FIELD_2_NAME: direct
   FIELD_3_NAME: direct_top
   Point: 0
   TITLE: Table
   VERSION: '0.2'
   diffuse_dtype: float64
   direct_dtype: float64
   direct_top_dtype: float64
   index_dtype: object
   index_names:
   - index
   lonlat:
   - 125.79006336856054
   - -33.65767449909174
   metadata: '`Pandas.DataFrame` converted to HDF5 compound datatype.'
   nrows: 8
   python_type: '`Pandas.DataFrame`'


IMAGE attributes
~~~~~~~~~~~~~~~~

IMAGE datasets will be written following the HDF5 `image specification <https://support.hdfgroup.org/HDF5/doc/ADGuide/ImageSpec.html>`_
with just the base amount of information such as:

* CLASS
* IMAGE_VERSION
* DISPLAY_ORIGIN

Images written as a whole at once using the *wagl.hdf5.write_h5_image* routine will attach *IMAGE_MINMAXRANGE* as an additional attribute.
All images with geospatial context, which within *wagl* should be all images, wiil attach the following two additional attributes:

* transform (GDAL like GeoTransform; 6 element array)
* crs_wkt (CRS stored as a variable length string using the Well Known Text specification

As mentioned previously, is a simple method similar to other geospatial formats for storing the corner tie point of the array with a real
world coordinate, along with the coordinate reference system. Both items are easily parsed to GDAL or rasterio for interpretation.
A *description* attribute is generally attached to every Image dataset as a means of easier understanding for anyone working with the code and wondering what a given image is representing.

An example of the yaml document, as extracted using *wagl_convert*, for an IMAGE dataset written tile by tile is given as follows:

.. code-block:: yaml

   CLASS: IMAGE
   DISPLAY_ORIGIN: UL
   description: Contains the solar azimuth angle in degrees.
   IMAGE_VERSION: '1.2'
   crs_wkt: PROJCS["GDA94 / MGA zone 52",GEOGCS["GDA94",DATUM["Geocentric_Datum_of_Australia_1994",SPHEROID["GRS
       1980",6378137,298.257222101,AUTHORITY["EPSG","7019"]],TOWGS84[0,0,0,0,0,0,0],AUTHORITY["EPSG","6283"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4283"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",129],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],AXIS["Northing",NORTH],AUTHORITY["EPSG","28352"]]
   geotransform:
   - 202325.0
   - 25.0
   - 0.0
   - 6271175.0
   - 0.0
   - -25.0
   no_data_value: -999

An example of the yaml document, as extracted using *wagl_convert*, for an IMAGE dataset written using the *wagl.hdf5.write_h5_image* routine is given as follows:

.. code-block:: yaml

   CLASS: IMAGE
   DISPLAY_ORIGIN: UL
   description: Contains the interpolated result of factor a for band 6 from sensor Landsat-8.
   IMAGE_MINMAXRANGE:
   - -999.0
   - 32.484375
   IMAGE_VERSION: '1.2'
   crs_wkt: PROJCS["GDA94 / MGA zone 52",GEOGCS["GDA94",DATUM["Geocentric_Datum_of_Australia_1994",SPHEROID["GRS
       1980",6378137,298.257222101,AUTHORITY["EPSG","7019"]],TOWGS84[0,0,0,0,0,0,0],AUTHORITY["EPSG","6283"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4283"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",129],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],AXIS["Northing",NORTH],AUTHORITY["EPSG","28352"]]
   geotransform:
   - 202325.0
   - 25.0
   - 0.0
   - 6271175.0
   - 0.0
   - -25.0
   interpolation_method: linear
   no_data_value: -999


Compression filters
-------------------

By default, wagl (via) h5py, will provide access to sever filters:

* lzf
* gzip
* shuffle

The default filters used by wagl is the shuffle filter and lzf compression filter.
The lzf filter is geared around speed, whilst still having modest compression.
The shuffle filter is designed to reorganise the data so that similar bytes are
closer together, thus potentially gaining better compression ratios.

Additional HDF5 compression filters (optional)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Additional compression filters can be used via HDF5's
`dynamically loaded filters <https://support.hdfgroup.org/HDF5/doc/Advanced/DynamicallyLoadedFilters/HDF5DynamicallyLoadedFilters.pdf>`_.
Essentially the filter needs to be compiled against the HDF5 library, and
installed into HDF5's plugin path, or a path of your choosing, and set the
HDF5_PLUGIN_PATH environment variable. The filters are then automatically
accessible by HDF5 via the `integer code  <https://support.hdfgroup.org/services/contributions.html>`_
assigned to the filter.

Mafisc compression filter
~~~~~~~~~~~~~~~~~~~~~~~~~

Mafisc combines both a bitshuffling filter and lzma compression filter in order
to get the best compression possible at the cost of lower compression speeds.
To install the *mafisc* compression filter, follow these `instructions <https://wr.informatik.uni-hamburg.de/research/projects/icomex/mafisc>`_.

Bitshuffle
~~~~~~~~~~

The `bitshuffle filter <https://github.com/kiyo-masui/bitshuffle>`_ can be installed
from source, or conda via the supplied `conda recipe <https://github.com/kiyo-masui/bitshuffle/tree/master/conda-recipe>`_.
It utilises a bitshuffling filter on top of either a lz4 or lzf compression filter.
