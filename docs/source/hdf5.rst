HDF5
====

As of gaip v5.0.0, the storage backend for all dataset outputs uses HDF5, whether the datasets be imagery, tables, or scalars.
Currently datasets can be stored according to the *IMAGE* Class, or *TABLE* Class.

* **IMAGE** Class `specification <https://support.hdfgroup.org/HDF5/doc/ADGuide/ImageSpec.html>`_
* **TABLE** Class `specification <https://support.hdfgroup.org/HDF5/doc/HL/H5TB_Spec.html>`_

Moving to HDF5
--------------

Some of the main reasons behind moving to HDF5 are:

* Reduce the number of files being written to disk. Yes imagery could be stored multi-band, but band ordering would need to be explicity defined.
* Consequently, a reduction in the number of files, also reduces the number of parameters being passed back and forth between `luigi <https://github.com/spotify/luigi>`_ *Tasks*.
* Consolodate the numerous text files that could be comma, tab, space, and variable space delimited variations, into a consistant *TABLE* structure that allows compression.
* *TABLE* and *IMAGE* datasets can be stored and accessed together from the same file.
* HDF5 has an `ExternalLink <http://docs.h5py.org/en/latest/high/group.html#group-extlinks>`_ feature, which for the case of the multitude of interpolated factors for each spectral band, all the outputs can be combined and accessed from the same file, thereby simplifying the workflow and functional architecture of gaip.
* Various `compression <https://support.hdfgroup.org/services/contributions.html>`_ filters are available, ranging from extremely high I/O speed, to slower I/O but very high compression ratios.
  High compression ratios allow for better archiving, and with conjuction of being able to store *TABLES*, *IMAGES*, *SCALARS*, variable resolutions, all within the same file, can aid in keeping the number of inodes down on HPC systems.
  This presents an archiving solution without the need to create tarballs, whilst still maintaining the capacity to operate with the archives. Essentially creating a kind of operational archive.
* Datasets can be stored in a hierarchical fashion, for instance, creating a hierarchical layout based on resolution, or the same dataset name stored in different locations based on a given attribute.
  For example, the ancillary data required for surface brightness temperature, gathers data at a minimum of 25 points across a scene or granule.
  The point id is used as the group label to differentiate between the same ancillary data gathered at different locations. i.e.:
    * /point-0/temperature
    * /point-1/temperature
    * /point-2/temperature
* Metadata; gaip can now store a lot more metadata such as longitude and latitude information with each ancillary point location, as opposed to a plain text label.
  Parameter settings used for a given algorithm such as for the satellite and solar angles calculation can be stored alongside the results, potentially making it easier for validation, and archive comparisons to be undertaken.
* Utilise a consistant library for data I/O, rather than a dozen or so different libraries. This helps to simplify the gaip data model, and have fewer library dependencies.
* It simplified the workflow **A LOT**

Dataset Names
-------------

Dataset names for each output are as follows:

* ancillary.h5
    * /aerosol
    * /coordinator
    * /elevation
    * /ozone
    * /water-vapour
    * /BRDF-Band-{band_id}-geo (combinations based on band_id)
    * /BRDF-Band-{band_id}-iso
    * /BRDF-Band-{band_id}-vol
    * /brdf-image-datasets/Band_{band_id)_0459_0479nm_geo (combinations based on band_id and lower/upper brdf wavelength)
    * /brdf-image-datasets/Band_{band_id)_0459_0479nm_iso
    * /brdf-image-datasets/Band_{band_id)_0459_0479nm_vol
* atmospheric-inputs.h5
    * /modtran-inputs/point-0/albedo-0/tp5_data (combinations based on point labe;, albedo label)
    * /modtran-inputs/point-1/albedo-1/tp5_data
    * /modtran-inputs/point-2/albedo-t/tp5_data
    * /modtran-inputs/point-3/albedo-th/tp5_data
* accumulated-solar-irradiance.h5
    * /point-0/albedo-0/channel (combinations based on point label, albedo label)
    * /point-0/albedo-0/solar-irradiance (combinations based on point label, albedo label)
* coefficients.h5
    * /coefficients
* bilinearly-interpolated-data.h5
    * /a-band-{band_id} (combinations based on the band_id)
    * /b-band-{band_id} (combinations based on the band_id)
    * /dif-band-{band_id} (combinations based on the band_id)
    * /dir-band-{band_id} (combinations based on the band_id)
    * /fs-band-{band_id} (combinations based on the band_id)
    * /fv-band-{band_id} (combinations based on the band_id)
    * /s-band-{band_id} (combinations based on the band_id)
    * /ts-band-{band_id} (combinations based on the band_id)
* dsm-extract.h5
    * /dsm
    * /dsm-smoothed
* exiting-angles.h5
    * /azimuthal-exiting
    * /exiting
* incident-angles.h5
    * /azimuthal-incident
    * /incident
* longitude-latitude.h5
    * /longitude
    * /latitude
* relative-slope.h5
   * /relative-slope
* satellite-solar.h5
    * /boxline
    * /centreline
    * /parameters/orbital-elements
    * /parameters/satellite-model
    * /parameters/satellite-track
    * /parameters/spheroid
    * /relative-azimuth
    * /satellite-azimuth
    * /satellite-view
    * /solar-azimuth
    * /solar-zenith
* shadow-masks.h5
    * /cast-shadow-satellite
    * /cast-shadow-sun
    * /combined-shadow
    * /self-shadow
* slope-aspect.h5
    * /aspect
    * /slope
* standardised-data.h5
    * /brdf-reflectance-band-{band_id} (combinations based on the band_id)
    * /lambertian-reflectance-band-{band_id} (combinations based on the band_id)
    * /terrain-reflectance-band-{band_id} (combinations based on the band_id)
    * /surface-brightness-temperature-band-{band_id} (combinations based on the band_id)

Geospatial Information
----------------------

Geospatial information for *IMAGE* Class datasets can be stored in various different ways. For gaip, we attach 2 attributes specifically related to geospatial context:

* transform (GDAL like GeoTransform; 6 element array)
* crs_wkt (CRS stored as a variable length string using the Well Known Text specification

This approach is very simple, and similar to lots of other mainstream formats such as `ENVI <https://www.harrisgeospatial.com/docs/ENVIHeaderFiles.html>`_,
`KEA <http://kealib.org/>`_. The geospatial information can automatically be interpreted using *gaip.geobox.GriddedGeoBox*.

Tables
------

Tabulated data created by gaip is stored in HDF5 using the compound datatype, and read back into memory as either a custom *NumPy* datatype, or directly into a *pandas.DataFrame*.
Datatypes are mapped between HDF5 and NumPy as best as possible. Additional attached attributes inlcuded by gaip can aid in the transitional mapping.
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

          >>> from gaip.hdf5 import read_table
          >>> import h5py
          >>> fid = h5py.File('coordinator.h5', 'r')
          >>> df = read_table(fid, 'coordinator')
