Invoking wagl
=============

The wagl workflow can be called directly from the command line.
It uses `luigi's command line tool <http://luigi.readthedocs.io/en/stable/command_line.html>`_ to execute the workflow for producing Analysis Ready Data (ARD).

Three options are available for the wagl workflow:

* **STANDARD**; Executes both the *NBAR* (Surface Reflectance) and the *SBT* (Surface Brightness Temperature) workflows.
* **NBAR**; Executes only the *NBAR* (Surface Reflectance) workflow
* **SBT**; Executes only the *SBT* (Surface Brightness Temperature) workflow


Local scheduler
---------------

To run the entire standard ARD workflow using luigi's local scheduler, set the --workflow parameter to *STANDARD*:

.. code-block:: bash

   $ luigi --module wagl.multifile_workflow ARD --level1-list datasets.txt --workflow STANDARD --outdir /some/path --workers 4

To run the entire *NBAR* ARD workflow using luigi's local scheduler, set the --workflow parameter to *NBAR*:

.. code-block:: bash

   $ luigi --module wagl.multifile_workflow ARD --level1-list datasets.txt --workflow NBAR --outdir /some/path --workers 4

To run the entire *SBT* ARD workflow using luigi's local scheduler, set the --workflow parameter to *SBT*:

.. code-block:: bash

   $ luigi --module wagl.multifile_workflow ARD --level1-list datasets.txt --workflow SBT --outdir /some/path --workers 4 --local-scheduler


Central scheduler
-----------------

To run using luigi's `central scheduler <http://luigi.readthedocs.io/en/stable/central_scheduler.html>`_:

.. code-block:: bash

   $ luigid --background --pidfile <PATH_TO_PIDFILE> --logdir <PATH_TO_LOGDIR> --state-path <PATH_TO_STATEFILE>

   $ luigi --module wagl.multifile_workflow ARD --level1-list datasets.txt --workflow STANDARD --outdir /some/path --workers 4

To include the pixel quality workflow as part of the main workflow, you need to set the pixel-quality switch as such:

.. code-block:: bash

   $ luigi --module wagl.multifile_workflow ARD --level1-list datasets.txt --workflow STANDARD --pixel-quality --outdir /some/path --workers 4

Luigi will then execute, and manage, the entire ARD (Analysis Ready Data) workflow for every Level-1 dataset listed in *datasets.txt*.


Intersecting the workfow at a given task
----------------------------------------

If however, you want to run just a specific part of the workflow, say for example *CalculateCoefficients*, then you would need to
specify the following arguments:

--level1         /path/to/level1/dataset
--work-root      path/to/working/directory
--granule        granule id name; Default is None; and can be ignored for Landsat
--vertices       the number of points atmospherical calculations will run across the dataset; Default is 25
--base-dir       the base directory to contain the atmospheric calculations; Default is _atmospherics
--workflow          the ARD (Analysis Ready Data) workflow to produce *STANDARD*, *NBAR* or *SBT*
--compression    the compression filter used when writing the outputs to disk; Default is lzf
--separate       if set, then issue separate atmospheric calculation tasks for each albedo effect. The default is to combine multiple albedo calculations into a single task, with the advantage of having less files on disk.

An example of running the *CalculateCoefficients* Task using the local scehduler is:

.. code-block:: bash

   $ luigi --module wagl.multifile_workflow CalculateCoefficients \
     --level1 /path/to/LS5_TM_OTH_P51_GALPGS01-007_111_068_20000707 \
     --work-root /my/work/LS5_TM_OTH_P51_GALPGS01-007_111_068_20000707.wagl-work --workers 4 --local-scheduler
   

Tasks available via the command line
------------------------------------

The Tasks callable from the command line are:

* **ARD** (Issues full NBAR and/or SBT workflows for each level-1 in a list)
* **LinkwaglOutputs** (Issues DataStandardisation Tasks for each level-1 dataset and links the results into a single file)
* **DataStandardisation** (Issues SurfaceReflectance and SurfaceTemerature Tasks for each band in a level-1 dataset)
* **SurfaceReflectance** (Calculates terrain corrected surface reflectance for a given band in a level-1 dataset)
* **CalculateShadowMasks** (Issues *CalculateCastShadowSun*, *CalculateCastShadowSatellite*, and *SelfShadow* Tasks for a level-1 dataset)
* **CalculateCastShadowSatellite** (Executes the cast shadow algorithm for the satellite direction, for a level-1 dataset)
* **CalculateCastShadowSun** (Executes the cast shadow algorithm for the solar direction, for a level-1 dataset)
* **SelfShadow** (Executes the self shadoe algorithm for a level-1 dataset)
* **RelativeAzimuthSlope** (Calculates the relative azimuth on the sloping surface for a level-1 dataset)
* **ExitingAngles** (Calculates the exiting angles for a level-1 dataset)
* **IncidentAngles** (Calculates the incident angles for a level-1 dataset)
* **SlopeAndAspect** (Calculates the slope and aspect for a level-1 dataset)
* **DEMExtraction** (Extracts the DEM for a level-1 dataset)
* **InterpolateCoefficients** (Issues *InterpolateCoefficient* Tasks for each band, for each factor for a level-1 dataset)
* **InterpolateCoefficient** (Executes interpolation for a given band for a given factor)
* **CalculateCoefficients** (Calculates the atmospheric coefficients derived from running a radiative transfer algorithm such as `MODTRAN <http://modtran.spectral.com/>`_)
* **Atmospherics** (Issues AtmosphericsCase Tasks, for each point/vertex for each albedo)
* **AtmosphericsCase** (Executes `MODTRAN <http://modtran.spectral.com/>`_ for a given point location and albedo factor)
* **WriteTp5** (Creates the Tape5 files for each point location and albedo factor required by `MODTRAN <http://modtran.spectral.com/>`_)
* **CalculateSatelliteAndSolarGrids** (Calculates the satellite and solar angles for a given level-1 dataset)
* **CalculateLonLatGrids** (Calculates the longitude  and latitude grids for a given level-1 dataset)
* **AncillaryData** (Retrieves the ancillary data for a given level-1 dataset)
* **CallTask** (Call any task and submit a list of level-1 datasets to that task)

The added bonus is that luigi will take care of all prior dependencies required to run the chosen Task. To execute the same Task again, simply remove the output file,
and luigi will re-run the task without re-running any of the prior dependencies, unless those outputs are removed as well.


Help an a specific Task
-----------------------

Help on executing a Task can be retrieved, for example:

.. code-block:: bash

   $ luigi --module wagl.multifile_workflow CalculateCoefficients --help

   $ luigi --module wagl.multifile_workflow CalculateCoefficients --help-all

The number of workers to assign to the Task tree *--workers* tells luigi how many Tasks to run in parallel (for those tasks that don't depend on each other).
While not making the best use of luigi (for such a quick and simple workflow), it does aid in quick research and development for a single dataset to 100's of datasets,
using this simple workflow.

For even larger numbers of datasets, say several thousand or tens of thousands to be exectued as a single workflow, then an alternate luigi workflow can be implemented
such as the PBS task flow. In this example, luigi issues and monitors PBS jobs, each job kicking off an MPI scheduler.

PBS submission
--------------

For users on a system that utilises a `PBS <https://en.wikipedia.org/wiki/Portable_Batch_System>`_ scheduler, wagl provides a command line tool *wagl_pbs* for automatic job submission into a PBS queue. The tool can partition the list of datasets into roughly equally sized chunks, based on the number of nodes requested. For example, a list containing 600 datasets, and a job requesting 10 nodes, will partition the list into 10 blocks each containing 60 scenes that a given node will process. Two flavours of jobs can be submitted to the PBS queue in this way:

1. Individual single node jobs; i.e. A single node represents a single submitted job.

  * Advantages:

    * If a node finishes its block of datasets earlier, the whole job doesn't have to wait for the other nodes to finish, therefore higher CPU utilisation can be sustained for the jobs duration.

  * Disadvantages:

    * More jobs to monitor.
    * Queue limits can be quickly reached.
    * Single node jobs tend to stay in the PBS queue for longer than multi-node jobs.
    * Have to wait for all submitted jobs to finish, which is dependent on how well the PBS queue can allocate the resources.

2. A single batch job is submitted to the queue, and each requested node executes a job using PBSDSH.

  * Advantages:

    * A single job to monitor.
    * PBS tends to allocate large single job resources quite well.

  * Disadvantages:

    * Whilst the blocks of datasets allocated to each node are roughly equal, the time taken to process a dataset is not. Some scenes may not have the required ancillary and will be skipped or fail (filtering the list of scenes prior to job submission can help with this), partial scenes can also process quicker. This means that while 1 or more of the nodes in the enitire job request have finished, the whole job has to wait until other nodes have finished their jobs. This can result in lower CPU utilisation over the jobs duration.

The arguments for *wagl_pbs* are:

--level1-list        The input level1 dataset list.
--vertices           Number of vertices to evaluate the radiative transfer at. JSON styled string is required, eg '(3, 3)'.
--workflow              The type of ARD workflow to invoke, eg STANDARD, NBAR, SBT.
--method             The interpolation method to invoke, eg BILINEAR, SHEAR, RBF.
--pixel-quality      Whether to run the pixel quality workflow, if applicable, or not.
--buffer-distance    The distance in units by which to buffer an image's extents by.
--outdir             The base output directory.
--logdir             The base logging and scripts output directory.
--env                Environment script to source.
--nodes              The number of nodes to request.
--project            Project code to run under.
--queue              The type of queue to submit the job into, eg normal, express.
--hours              Job walltime in hours.
--email              Notification email address.
--local-scheduler    Use a local scheduler instead of a central scheduler.
--dsh                Run using PBS Distributed Shell.
--task               A luigi task defined within the wagl.multifile_workflow; eg *CalculateCoefficients*
--test               Test job execution (Don't submit the job to the PBS queue).

An example of submitting individual jobs to the PBS queue using the following specifications:

  * Run using the *NBAR* workflow.
  * The *BILINEAR* interpolation function.
  * Specify a 3x3 point grid location to calculate the radiative transfer at.
  * 10 nodes.
  * Use the nx200 project allocation code identifier.
  * Submit to the express queue.
  * Maximum job runtime of 2 hours.

.. code-block:: bash

   $ wagl_pbs --level1-list /path/to/level1-datasets.txt --vertices '(3, 3)' --workflow NBAR --method BILINEAR --outdir /path/to/the/output/directory --logdir /path/to/the/logs/directory --env /path/to/the/environment/script --nodes 10 --project nx200 --queue express --hours 2 --email your.name@something.com

The same job resources, but use PBSDSH instead of individual jobs being submitted to the PBS queue.

.. code-block:: bash

   $ wagl_pbs --level1-list /path/to/level1-datasets.txt --vertices '(3, 3)' --workflow NBAR --method BILINEAR --outdir /path/to/the/output/directory --logdir /path/to/the/logs/directory --env /path/to/the/environment/script --nodes 10 --project v10 --queue express --hours 2 --email your.name@something.com --dsh

Each call to *wagl_pbs* will generate a new batch id, and each node will be assigned a job id. In this way each node will have its logs and output data contained in its own directory structure.  For example:

.. code-block:: bash

  $ /base/logs/directory/batchid-b6cbadbe98/jobid-074cb6/
  $ /base/logs/directory/batchid-b6cbadbe98/jobid-113f33/
  $ /base/logs/directory/batchid-b6cbadbe98/jobid-5b00d6/
  $ /base/output/directory/batchid-b6cbadbe98/jobid-074cb6/
  $ /base/output/directory/batchid-b6cbadbe98/jobid-113f33/
  $ /base/output/directory/batchid-b6cbadbe98/jobid-5b00d6/


Intersecting the wagl workflow, and have it execute across a list of datasets
-----------------------------------------------------------------------------

The *--task* command line option for *wagl_pbs* allows the user to have specific control of the workflow, whilst still retaining the capability of running it in bulk over a list of datasets.
The example below only executes the workflow up to the end of CalculateCoefficients, and only for a single dataset. This is because most of the luigi tasks defined in wagl.multifie_workflow are for a given dataset's group and granules.

.. code-block:: bash

   $ luigi --module wagl.multifile_workflow CalculateCoefficients \
     --level1 /path/to/LS5_TM_OTH_P51_GALPGS01-007_111_068_20000707 \
     --work-root /my/work/LS5_TM_OTH_P51_GALPGS01-007_111_068_20000707.wagl-work --workers 4 --local-scheduler
   
The bulk submission workflow entrypoint is defined in the luigi Task named *ARD*, which initialise the entire wagl.multifile_workflow tree. In order to submit a list of datasets but only execute a partial workflow such as *CalculateCoefficients*, then a generic luigi task class named *CallTask* has been defined for this very purpose.

The example below will run the *CalculateCoefficients* for each input dataset:

.. code-block:: bash

   $ luigi --module wagl.multifile_workflow CallTask --level1-list /path/to/level1-datasets.txt --outdir /path/to/the/output/directory --task CalculateCoefficients

The example below is using the *wagl_pbs* command line utility:

.. code-block:: bash

   $ wagl_pbs --level1-list /path/to/level1-datasets.txt --outdir /path/to/the/output/directory --logdir /path/to/the/logs/directory --env /path/to/the/environment/script --nodes 10 --project v10 --queue express --hours 2 --email your.name@something.com --dsh --task CalculateCoefficients

You might notice that no arguments such as *--workflow*, *--vertices* or *--method* are present. This is because in order for the CallTask to be generic, it's easier to let any parameters that need parsing, and specify them using the *luigi.cfg* file and have luigi do all the work of parsing additional parameters.

An example configuration for executing the *CalculateCoefficients* task and its dependencies, for a list of datasets is given by:

.. code-block:: cfg

   [CalculateCoefficients]
   vertices = (15, 15)
   workflow = NBAR

This will parse in a 15x15 point grid at which to evaluate the radiative transfer, and only for the *NBAR* workflow.

The *CallTask* luigi task will work for any task in the *wagl.multifile_workflow* if the first 3 arguments of a task are:
[level1 (file pathname), work_root (directory pathname), granule]

or for tasks that contain a datasets *group* parameter, the first 4 arguments of a task should be:
[level1 (file pathname), work_root (directory pathname), granule, group]


Singlefile workflow
-------------------

There are two module workflows available from the command line, each sharing much the same command line arguments:

* multifile_workflow
* singlefile_workflow

The above examples present the multifile workflow which is suitable for testing, experimentation, or simply rapidly producing a small number of outputs.
The singlefile workflow presents a case more suitable for mass routine production, which as the name of the module suggests, outputs a single file.
This makes it less demanding on the filesystem, eg more bands, more resolutions, and more points, equal more files, and easier for any scheduler to track, and easier to distribute the single file to other people.
It could be thought of as an operational archive, that doesn't need to be untarred, or uncompressed, as the file can be accessed quite easily via `h5py <http://www.h5py.org/>`_ without decompressing the entire file.
It also makes it easier for automatic testing and evaluation to occur between different version of the same dataset, to not just test and compare the final outputs of the algorithm, but also for all the intermediate images, tables, constants that are calculated.
