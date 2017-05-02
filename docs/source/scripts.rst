Command line utility scripts
============================

There are several command line scripts for either extracting datasets or unittesting.

    * gaip_convert:  *An unpacking/converting utility that converts HDF5 Tables to CSV, HDF5 images to GeoTiff, and metadata to yaml files.*
    * test_calculate_angles: *Compares and evaluates each dataset contained within a* **satellite-solar.h5** *file against the same datasets contained within another file.*
    * test_dsm: *Compare and evaluate two* **dsm-extract.h5** *files.*
    * test_exiting_angles: *Compare and evaluate two* **exiting-angles.h5** *files.*
    * test_incident_angles: *Compare and evaluate two* **incident-angles.h5** *files.*
    * test_relative_slope *Compare and evaluate two* **relative-slope.h5** *files.*
    * test_shadow_masks: *Compare and evaluate two* **shadow-masks** *files.*
    * test_slope_aspect: *Compare and evaluate two* **slope-aspect.h5** *files.*

**gaip_convert**

Unpacks/converts SCALAR, TABLE & IMAGE HDF5 datasets to yaml, csv, and GeoTiff file respectively, along with an attributes to a yaml file. Any hierarchial structure will be replicated on disk as directories. Individual Groups or Datasets can be specified using the *--pathname* argument. Thereby only Dataset's contained under that Group, or the selected Dataset will be converted.

Example of use:

.. code-block:: bash

   $ gaip_convert --filename satellite-solar.h5 --outdir /some/output/directory

   $ gaip_convert --filename satellite-solar.h5 --outdir /some/output/directory --pathname /solar-zenith

**gaip_ls**

Lists the contents of a HDF5 file, printing the full pathname of each Group and Dataset from the root level of the HDF5 file, as well as the Class, eg *Group, Dataset, IMAGE Dataset, TABLE Dataset*. Optionally, attributes for the Group or Dataset can be output as well by specifying the *--verbose* argument.

Example of use:

.. code-block:: bash

   $ gaip_ls --filename satellite-solar.h5

   $ gaip_ls --filename satellite-solar.h5 --verbose

**test_calculate_angles**

Compares and evaluates the following datasets:

    * /boxline
    * /centreline
    * /relative-azimuth
    * /satellite-azimuth
    * /satellite-view
    * /solar-azimuth
    * /solar-zenith

.. code-block:: bash

   $ test_calculate_angles --reference_fname /reference/satellite-solar.h5 --test_fname /test/satellite-solar.h5

**test_dsm**

Compares and evaluates the following datasets:

    * /dsm
    * /dsm-smoothed

.. code-block:: bash

   $ test_dsm --reference_fname /reference/dsm-extract.h5 --test_fname /test/dsm-extract.h5

**test_exiting_angles**

Compares and evaluates the following datasets:

    * /azimuthal-exiting
    * /exiting

.. code-block:: bash

   $ test_exiting_angles --reference_fname /reference/exiting-angles.h5 --test_fname /test/exiting-angles.h5

**test_incident_angles**

Compares and evaluates the following datasets:

    * /azimuthal-incident
    * /incident

.. code-block:: bash

   $ test_incident_angles --reference_fname /reference/incident-angles.h5 --test_fname /test/incident-angles.h5

**test_relative_slope**

Compares and evaluates the following datasets:

   * /relative-slope

.. code-block:: bash

   $ test_relative_slope --reference_fname /reference/relative-slope.h5 --test_fname /test/relative-slope.h5

**test_shadow_masks**

Compares and evaluates the following datasets:

    * /cast-shadow-satellite
    * /cast-shadow-sun
    * /combined-shadow
    * /self-shadow

.. code-block:: bash

   $ test_shadow_masks --reference_fname /reference/shadow-masks.h5 --test_fname /test/shadow-masks.h5

**test_slope_aspect**

Compares and evaluates the following datasets:

    * /aspect
    * /slope

.. code-block:: bash

   $ test_slope_aspect --reference_fname /reference/slope-aspect.h5 --test_fname /test/slope-aspect.h5
