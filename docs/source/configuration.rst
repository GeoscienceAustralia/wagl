Configuration
=============

While the ARD (Analysis Ready Data) workflow can be executed at the command line, a config file named *luigi.cfg* can be used to override command line arguments, or even to supply 
parameters to Tasks that are higher up in the dependency tree. They can also be used to supply default parameters rather than specifying in the code directly.

An example *luigi.cfg* can be found `here <http://github.com/GeoscienceAustralia/wagl/blob/develop/configs/luigi.cfg>`_.

The easiest way to get the workflow to find the config file, is to place it in the working directory that the luigi workflow will be executed from.

Luigi configurable options
--------------------------

For configurable options for luigi see the `luigi docs <http://luigi.readthedocs.io/en/stable/configuration.html>`_, but here are a few that *wagl* makes use of:

Logging
~~~~~~~

.. code-block:: cfg

   [core]
   logging_conf_file = /path/to/logging.cfg

Task history
~~~~~~~~~~~~

The luigi history and task status can also be stored in a database, thus providing another source of status information that in conjuction with the logs could be useful to any number of users.
To output the task history to an sqlite database, modify your luigi.cfg file with the following contents:

.. code-block:: cfg

   [scheduler]
   record_task_history = True
   
   [task_history]
   db_connection = sqlite:///luigi-task-history.db

And luigi will output the task history to an sqlite database named *luigi-task-history.db* relative to the directory that the scheduler was launched from.


Workflow configurable options
-----------------------------

Most of the options aren't required as they're either internal defaults if you're fine with them, or they're inhertied parameters from dependent tasks.
If you're a user that would like to intersect any part of the dependency tree, then you'll have more control over the parameter inputs, and you'll be required to provide some but not all parameter inputs. A large portion of the Luigi tasks have default values for the parameters, which can be overidden by specifying values via the *luigi.cfg* file.

.. code-block:: cfg

   [DEFAULT]
   # Acts as a way of specifying a default parameter that will be created for each task
   # I'd only recommend it's use for tasks that actually have these items as parameters

   # The compression filter to use (internally the code defaults to use *lzf*)
   compression = lzf

   [AncillaryData]
   # These parameters determine the directory locations required for ancillary retrieval

   # A dict defining a user input value or a pathname to ancillary sources
   # internally the code defaults to {"user": 0.5}
   aerosol = {"user": 0.5}
   aerosol = {"pathname": "/path/to/data"}

   # File path name to the directory containing the MODIS BRDF data
   brdf_path = 

   # File path name to the directory containing the Pre-MODIS BRDF data
   brdf_premodis_path = 

   # File path name to the directory containing the ozone data
   ozone_path = 

   # A dict defining a user input value or a pathname to ancillary sources
   # internally the code defaults to {"user": 1.5}
   water_vapour = {"user": 1.5}
   water_vapour = {"pathname": "/path/to/data"}

   # File path name to the directory containing the world 1 degree DEM data
   dem_path = 

   # File path name to the directory containing the dewpoint data
   dewpoint_path = 

   # File path name to the directory containing the 2m (surface) temperature data
   temp_2m_path = 

   # File path name to the directory containing the surface pressure data
   surface_pressure_path = 

   # File path name to the directory containing the atmospheric layers geopotential data
   geopotential_path = 

   # File path name to the directory containing the atmospheric layers temperature data
   temperature_path = 

   # File path name to the directory containing the atmospheric layers relative humidity data
   relative_humidity_path = 

   # File path name to a file containing the invariant geopotential height data
   invariant_height_fname = 

   [CalculateLonLatGrids]
   # The compression filter to use (internally the code defaults to use *lzf*)
   compression = lzf

   [CalculateSatelliteAndSolarGrids]
   # File path name to the directory containing the Two-line-element data
   tle_path = 

   # The compression filter to use (internally the code defaults to use *lzf*)
   compression = lzf

   [WriteTp5]
   This controls the tp5 file creation required for input into MODTRAN.

   # A name indicating the base directory to output the result to
   # internally defaults to _atmospherics
   base_dir = _atmospherics

   # The compression filter to use (internally the code defaults to use *lzf*)
   compression = lzf

   # The number of vertices required for evaluating the radiative transfer over
   vertices = (5, 5)

   # The model run to use; *standard*, *nbar*, or *sbt*
   model = standard

   [AtmosphericsCase]
   # This controls the running of MODTRAN
   # most of the parameters are inherited from the *WriteTp5* task

   # A name indicating the base directory to output the result to
   # internally defaults to _atmospherics
   base_dir = _atmospherics

   # The compression filter to use (internally the code defaults to use *lzf*)
   compression = lzf

   # The number of vertices required for evaluating the radiative transfer over
   vertices = (5, 5)

   # The point id to be run
   point = 

   # A *list* containing the albedo factor to be run
   albedos = 

   # A file path name to the MODTRAN executable
   exe = 

   [Atmospherics]
   # This controls the submition of *AtmosphericsCase* taks, and most of the
   # parameters are inherited from the *WriteTp5* task

   # A name indicating the base directory to output the result to
   # internally defaults to _atmospherics
   base_dir = _atmospherics

   # The compression filter to use (internally the code defaults to use *lzf*)
   compression = lzf

   # The number of vertices required for evaluating the radiative transfer over
   # internally defaults to (5, 5)
   vertices = (5, 5)

   # The model run to use; *standard*, *nbar*, or *sbt*
   # internally defaults to standard
   model = standard

   # A *boolean* to indicate whether MODTRAN evaluations for a single point should
   # be issued as separate tasks, or combined together in a single process
   # internally defaults to False
   separate = false

   [CalculateCoefficients]
   # Same options as the *Atmospherics* task.

   [InterpolateCoefficient]
   # A name indicating the base directory to output the results to
   # internally defaults to _interpolation
   base_dir = _interpolation

   # The compression filter to use (internally the code defaults to use *lzf*)
   compression = lzf

   # The number of vertices required for evaluating the radiative transfer over
   # internally defaults to (5, 5)
   vertices = (5, 5)

   # The model run to use; *standard*, *nbar*, or *sbt*
   # internally defaults to standard
   model = standard

   # The factor id to run
   factor = 

   # The band number to run
   band_name = 

   # The interpolation method to use;
   # *bilinear*, *fbilinear*, *shear*, *shearb*, or *rbf*
   # internally defaults to shear
   method = shear

   [InterpolateCoefficients]
   # The number of vertices required for evaluating the radiative transfer over
   vertices = (5, 5)

   # The model run to use; *standard*, *nbar*, or *sbt*
   model = standard

   # The compression filter to use (internally the code defaults to use *lzf*)
   compression = lzf

   # The interpolation method to use;
   # *bilinear*, *fbilinear*, *shear*, *shearb*, or *rbf*
   method = shear

   [DEMExctraction]
   # The compression filter to use (internally the code defaults to use *lzf*)
   compression = lzf
   # The distance in units by which to buffer an image's extents by
   # (internally defaults to 8000)
   buffer_distance = 8000

   [SlopeAndAspect]
   # The compression filter to use (internally the code defaults to use *lzf*)
   compression = lzf

   [IncidentAngles]
   # The compression filter to use (internally the code defaults to use *lzf*)
   compression = lzf

   [ExitingAngles]
   # The compression filter to use (internally the code defaults to use *lzf*)
   compression = lzf

   [RelativeAzimuthSlope]
   # The compression filter to use (internally the code defaults to use *lzf*)
   compression = lzf

   [SelfShadow]
   # A name indicating the base directory to output the results to
   # internally defaults to _shadow
   base_dir = _shadow

   # The compression filter to use (internally the code defaults to use *lzf*)
   compression = lzf

   [CalculateCastShadowSun]
   # A name indicating the base directory to output the results to
   # internally defaults to _shadow
   base_dir = _shadow

   # The compression filter to use (internally the code defaults to use *lzf*)
   compression = lzf

   [CalculateCastShadowSatellite]
   # A name indicating the base directory to output the results to
   # internally defaults to _shadow
   base_dir = _shadow

   # The compression filter to use (internally the code defaults to use *lzf*)
   compression = lzf

   [CalculateShadowMasks]
   # The compression filter to use (internally the code defaults to use *lzf*)
   compression = lzf

   [SurfaceReflectance]
   # A floating point value for surface reflectance adjustment (Fuqin to document)
   # internally defaults to 0.52
   rori = 0.52

   # A name indicating the base directory to output the results to
   # internally defaults to _standardised
   base_dir = _standardised

   [SurfaceTemperature]
   # A name indicating the base directory to output the results to
   # internally defaults to _standardised
   base_dir = _standardised

   [Standard]
   # A boolean indicating whether or not to run the pixel quality workflow
   # default is false
   pixel_quality = false

   [LinkwaglOutputs]
   # The path to the level-1 dataset
   level1

   # The model run to use; *standard*, *nbar*, or *sbt*
   model = standard

   # The number of vertices required for evaluating the radiative transfer over
   vertices = (5, 5)

   # A boolean indicating whether or not to run the pixel quality workflow
   # default is false
   pixel_quality = false

   # The interpolation method to use;
   # *bilinear*, *fbilinear*, *shear*, *shearb*, or *rbf*
   method = shear

   [ARD]
   # The model run to use; *standard*, *nbar*, or *sbt*
   model = standard

   # The number of vertices required for evaluating the radiative transfer over
   vertices = (5, 5)

   # A boolean indicating whether or not to run the pixel quality workflow
   # default is false
   pixel_quality = false

   # The interpolation method to use;
   # *bilinear*, *fbilinear*, *shear*, *shearb*, or *rbf*
   method = shear
