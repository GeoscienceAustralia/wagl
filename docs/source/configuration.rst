Configuration
=============

While the ARD (Analysis Ready Data) workflow can be executed at the command line, a config file named *luigi.cfg* can be used to override command line arguments, or even to supply 
parameters to Tasks that are higher up in the dependency tree. They can also be used to supply default parameters rather than specifying in the code directly.

An example *luigi.cfg* can be found at config/luigi.cfg.

The easiest way to get the workflow to find the config file, is to place it in the working directory that the luigi workflow was executed from.


Configurable options
--------------------

Most of the options aren't required as they're either internal defaults if you're fine with them, or they're inhertied parameters from dependent tasks.
If you're a user that would like to intersect any part of the dependency tree, then you'll have more control over the parameter inputs, and you'll be required to provide some.


[DEFAULT]
---------

Acts as a way of specifying a default parameter that will be created for each task. I'd only recommend it's use for tasks that actually have these as parameters.

compression
  The compression filter to use (internally the code defaults to use *lzf*.

y_tile
  The y-axis tile size to be used for processing.


[AncillaryData]
---------------

These parameters determine the directory locations required for ancillary retrieval.

aerosol_fname
  File path name to a file containing the aerosol data.

brdf_path
  File path name to the directory containing the MODIS BRDF data.

brdf_premodis_path
  File path name to the directory containing the Pre-MODIS BRDF data.

ozone_path
  File path name to the directory containing the ozone data.

water_vapour_path
  File path name to the directory containing the water vapour data.

dem_path
  File path name to the directory containing the world 1 degree DEM data.

dewpoint_path
  File path name to the directory containing the dewpoint data.

temp_2m_path
  File path name to the directory containing the 2m (surface) temperature data.

surface_pressure_path
  File path name to the directory containing the surface pressure data.

geopotential_path
   File path name to the directory containing the atmospheric layers geopotential data.

temperature_path
  File path name to the directory containing the atmospheric layers temperature data.

relative_humidity_path
  File path name to the directory containing the atmospheric layers relative humidity data.

invariant_height_fname
  File path name to a file containing the invariant geopotential height data.


[CalculateLonLatGrids]
----------------------

compression
  The compression filter to use (internally the code defaults to use *lzf*.


[CalculateSatelliteAndSolarGrids]
---------------------------------

tle_path
  File path name to the directory containing the Two-line-element data.

compression
  The compression filter to use (internally the code defaults to use *lzf*.


[WriteTp5]
----------

This controls the tp5 file creation required for input into MODTRAN.

base_dir
  A name indicating the base directory to output the result to.

compression
  The compression filter to use (internally the code defaults to use *lzf*.

vertices
  The number of vertices required for evaluating the radiative transfer over.

model
  The model run to use; *standard*, *nbar*, or *sbt*.


[AtmosphericsCase]
------------------

This controls the running of MODTRAN, and most of the parameters are inherited
from the *WriteTp5* task.

base_dir
  A name indicating the base directory to output the result to.

compression
  The compression filter to use (internally the code defaults to use *lzf*.

vertices
  The number of vertices required for evaluating the radiative transfer over.

point
  The point id to be run.

albedos
  A *list* containing the albedo factor to be run.

exe
  A file path name to the MODTRAN executable.


[Atmospherics]
--------------

This controls the submition of *AtmosphericsCase* taks, and most of the parameters are inherited
from the *WriteTp5* task.

base_dir
  A name indicating the base directory to output the result to.

compression
  The compression filter to use (internally the code defaults to use *lzf*.

vertices
  The number of vertices required for evaluating the radiative transfer over.

model
  The model run to use; *standard*, *nbar*, or *sbt*.

combined
  A *boolean* to indicate whether MODTRAN evaluations for a single point should
  be combined together in a single process.


[CalculateCoefficients]
-----------------------

Same options as the *Atmospherics* task.


[BilinearInterpolationBand]
---------------------------

base_dir
  A name indicating the base directory to output the results to.
  Internally defaults to _bilinear.

compression
  The compression filter to use (internally the code defaults to use *lzf*.

vertices
  The number of vertices required for evaluating the radiative transfer over.

model
  The model run to use; *standard*, *nbar*, or *sbt*.

factor
  The factor id to run.

band_num
  The band number to run.

method
  THe interpolation method to use; *linear*, *shear* or *rbf*. The default is *shear*.


[BilinearInterpolation]
-----------------------

vertices
  The number of vertices required for evaluating the radiative transfer over.

model
  The model run to use; *standard*, *nbar*, or *sbt*.

compression
  The compression filter to use (internally the code defaults to use *lzf*.


[DEMExctraction]
----------------

compression
  The compression filter to use (internally the code defaults to use *lzf*.


[SlopeAndAspect]
----------------

compression
  The compression filter to use (internally the code defaults to use *lzf*.

y_tile
  The y-axis tile size to be used for processing.


[IncidentAngles]
----------------

compression
  The compression filter to use (internally the code defaults to use *lzf*.

y_tile
  The y-axis tile size to be used for processing.


[ExitingAngles]
---------------

compression
  The compression filter to use (internally the code defaults to use *lzf*.

y_tile
  The y-axis tile size to be used for processing.


[RelativeAzimuthSlope]
----------------------

compression
  The compression filter to use (internally the code defaults to use *lzf*.

y_tile
  The y-axis tile size to be used for processing.


[SelfShadow]
------------

base_dir
  A name indicating the base directory to output the results to.
  Internally defaults to _shadow.

compression
  The compression filter to use (internally the code defaults to use *lzf*.

y_tile
  The y-axis tile size to be used for processing.


[CalculateCastShadowSun]
------------------------

base_dir
  A name indicating the base directory to output the results to.
  Internally defaults to _shadow.

compression
  The compression filter to use (internally the code defaults to use *lzf*.

y_tile
  The y-axis tile size to be used for processing.


[CalculateCastShadowSatellite]
------------------------------

base_dir
  A name indicating the base directory to output the results to.
  Internally defaults to _shadow.

compression
  The compression filter to use (internally the code defaults to use *lzf*.

y_tile
  The y-axis tile size to be used for processing.


[CalculateShadowMasks]
----------------------

compression
  The compression filter to use (internally the code defaults to use *lzf*.

y_tile
  The y-axis tile size to be used for processing.


[SurfaceReflectance]
--------------------
rori
  A floating point value. Internally defaults to 0.52.

base_dir
  A name indicating the base directory to output the results to.
  Internally defaults to _standardised.


[SurfaceTemperature]
--------------------

base_dir
  A name indicating the base directory to output the results to.
  Internally defaults to _standardised.


[Standard]
----------

pixel_quality
  A boolean indicating whether or not to run the pixel quality workflow.


[ARD]
-----

work_extension
  The extension to append to the working directory.

model
  The model run to use; *standard*, *nbar*, or *sbt*.

vertices
  The number of vertices required for evaluating the radiative transfer over.

pixel_quality
  A boolean indicating whether or not to run the pixel quality workflow.
