HDF5
====

As of gaip v5.0.0, the storage backend for all dataset outputs uses HDF5, whether the datasets be imagery, tables, or scalars.
Currently datasets can be stored according to the *IMAGE* Class, or *TABLE* Class.

* **IMAGE** Class `specification <https://support.hdfgroup.org/HDF5/doc/ADGuide/ImageSpec.html>`_.
* **TABLE** Class `specification <https://support.hdfgroup.org/HDF5/doc/HL/H5TB_Spec.html>`_

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

Dataset names for each output are as follows:
