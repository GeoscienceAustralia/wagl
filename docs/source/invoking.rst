Invoking the Geoscience Australia Image Processor
=================================================

The gaip image processor can be called directly from the command line.
It uses `luigi's command line tool <http://luigi.readthedocs.io/en/stable/command_line.html>`_ to execute the workflow for producing Analysis Ready Data (ARD).

Three options are available for the gaip model workflow:

* **standard**; Executes both the *NBAR* (Surface Reflectance) and the *SBT* (Surface Brightness Temperature) workflows.
* **nbar**; Executes only the *NBAR* (Surface Reflectance) workflow
* **sbt**; Executes only the *SBT* (Surface Brightness Temperature) workflow

To run the entire standard ARD workflow using luigi's local scheduler, set the --model parameter to *standard*:

.. code-block:: bash

   $ luigi --module gaip.standard_workflow ARD --level1-csv scenes.txt --model standard --output-directory /some/path --workers 4

To run the entire *nbar* ARD workflow using luigi's local scheduler, set the --model parameter to *nbar*:

.. code-block:: bash

   $ luigi --module gaip.standard_workflow ARD --level1-csv scenes.txt --model nbar --output-directory /some/path --workers 4

To run the entire *sbt* ARD workflow using luigi's local scheduler, set the --model parameter to *sbt*:

.. code-block:: bash

   $ luigi --module gaip.standard_workflow ARD --level1-csv scenes.txt --model sbt --output-directory /some/path --workers 4 --local-scheduler

To run using luigi's `central scheduler <http://luigi.readthedocs.io/en/stable/central_scheduler.html>`_:

.. code-block:: bash

   $ luigid --background --pidfile <PATH_TO_PIDFILE> --logdir <PATH_TO_LOGDIR> --state-path <PATH_TO_STATEFILE>

   $ luigi --module gaip.standard_workflow ARD --level1-csv scenes.txt --model standard --output-directory /some/path --workers 4

To include the pixel quality workflow as part of the main workflow, you need to set the pixel-quality switch as such:

.. code-block:: bash

   $ luigi --module gaip.standard_workflow ARD --level1-csv scenes.txt --model standard --pixel-quality --output-directory /some/path --workers 4

Luigi will then execute, and manage, the entire ARD (Analysis Ready Data) workflow for every Level-1 scene listed in *scenes.txt*.

If however, you want to run just a specific part of the workflow, say for example *CalculateCoefficients*, then you would need to
specify the following arguments:

--level1         /path/to/level1/scene
--work-root      path/to/working/directory
--granule        granule id name; Default is None; and can be ignored for Landsat
--vertices       the number of points atmospherical calculations will run across the scene; Default is 25
--base-dir       the base directory to contain the atmospheric calculations; Default is _atmospherics
--model          the ARD (Analysis Ready Data) workflow to produce *standard*, *nbar* or *sbt*
--compression    the compression filter used when writing the outputs to disk; Default is lzf
--pixel-quality  if set then the workflow will execute the pixel quality workflow, if a scene can be processed through the pixel quality workflow

An example of running the *CalculateCoefficients* Task using the local scehduler is:

.. code-block:: bash

   $ luigi --module gaip.standard_workflow CalculateCoefficients \
     --level1 /path/to/LS5_TM_OTH_P51_GALPGS01-007_111_068_20000707 \
     --work-root /my/work/LS5_TM_OTH_P51_GALPGS01-007_111_068_20000707.gaip-work --workers 4 --local-scheduler
   
The Tasks callable from the command line are:

* **ARD** (Issues full NBAR and/or SBT workflows for each level-1 in a list)
* **Standard** (Issues SurfaceReflectance and SurfaceTemerature Tasks for each band in a level-1 scene)
* **SurfaceReflectance** (Calculates terrain corrected surface reflectance for a given band in a level-1 scene)
* **CalculateShadowMasks** (Issues *CalculateCastShadowSun*, *CalculateCastShadowSatellite*, and *SelfShadow* Tasks for a level-1 scene)
* **CalculateCastShadowSatellite** (Executes the cast shadow algorithm for the satellite direction, for a level-1 scene)
* **CalculateCastShadowSun** (Executes the cast shadow algorithm for the solar direction, for a level-1 scene)
* **SelfShadow** (Executes the self shadoe algorithm for a level-1 scene)
* **RelativeAzimuthSlope** (Calculates the relative azimuth on the sloping surface for a level-1 scene)
* **ExitingAngles** (Calculates the exiting angles for a level-1 scene)
* **IncidentAngles** (Calculates the incident angles for a level-1 scene)
* **SlopeAndAspect** (Calculates the slope and aspect for a level-1 scene)
* **DEMExtraction** (Extracts the DEM for a level-1 scene)
* **BilinearInterpolation** (Issues *BilinearInterpolationBand* Tasks for each band, for each factor for a level-1 scene)
* **BilinearInterpolationBand** (Executes the bilinear interpolation for a given band for a given factor)
* **CalculateCoefficients** (Calculates the atmospheric coefficients derived from running a radiative transfer algorithm such as `MODTRAN <http://modtran.spectral.com/>`_)
* **Atmospherics** (Issues AtmosphericsCase Tasks, for each point/vertex for each albedo)
* **AtmosphericsCase** (Executes `MODTRAN <http://modtran.spectral.com/>`_ for a given point location and albedo factor)
* **WriteTp5** (Creates the Tape5 files for each point location and albedo factor required by `MODTRAN <http://modtran.spectral.com/>`_)
* **CalculateSatelliteAndSolarGrids** (Calculates the satellite and solar angles for a given level-1 scene)
* **CalculateLonLatGrids** (Calculates the longitude  and latitude grids for a given level-1 scene)
* **AncillaryData** (Retrieves the ancillary data for a given level-1 scene)

The added bonus is that luigi will take care of all prior dependencies required to run the chosen Task. To execute the same Task again, simply remove the output file,
and luigi will re-run the task without re-running any of the prior dependencies, unless those outputs are removed as well.

Help on executing a Task can be retrieved, for example:

.. code-block:: bash

   $ luigi --module gaip.standard_workflow CalculateCoefficients --help

   $ luigi --module gaip.standard_workflow CalculateCoefficients --help-all

The number of workers to assign to the Task tree *--workers* tells luigi how many Tasks to run in parallel (for those tasks that don't depend on each other).
While not making the best use of luigi (for such a quick and simple workflow), it does aid in quick research and development for a single scene to 100's of scenes,
using this simple workflow.

For even larger numbers of scenes, say several thousand or tens of thousands to be exectued as a single workflow, then an alternate luigi workflow can be implemented
such as the PBS task flow. In this example, luigi issues and monitors PBS jobs, each job kicking off an MPI scheduler.
