
====================================================
Standard Operating Procedure (AGDC Fractional Cover)
====================================================

AGDC Fractional Cover
---------------------

Purpose
-------
The fractional cover program (agdc_fc.py_) takes NBAR images from the datacube  and produces a Fractional Cover product as output. 

.. _agdc_fc.py: https://github.com/GeoscienceAustralia/ga-neo-landsat-processor/blob/develop/workflow/agdc_fc.py

This document describes how to run the Fractional Cover job via PBS script to process production data. 

The Fractional Cover PBS job uses Luigi to run multiple instances of the fractional cover program in parallel. By adjusting the system resources used by the PBS job, hundreds (or thousands) of NBAR input files may be processed in a short period of time.

* `End Members`_ data are maintained in code (changes are applied through revision control)
* output scenes (Fraction Cover scenes) are written to the FC Output Directory
* log files are written to the logs directory

.. _End Members: https://github.com/GeoscienceAustralia/ga-neo-landsat-processor/blob/develop/gaip/endmembers.py

Obtaining the code
------------------
The production code can be found in the GitHub repository. Obtain a copy of the production code into a convenient directory by following this procedure:

.. code-block:: bash

 cd <base directory of your choice>
 git clone git@github.com:GeoscienceAustralia/ga-neo-landsat-processor.git
 cd ga-neo-landsat-processor
 git checkout develop

Now change directory into the GA image processing directory (gaip)

.. code-block:: bash

 cd ga-neo-landsat-processor/workflow

Key scripts
-----------
The `fc.cfg`_ config file is used to configure the processing job. It is used to
query the Australian Geoscience Datacube, and dynamically generate the PBS
submission script and execute it, thereby submitting it to the queue.
Within the config file, the section titled `[agdc]` is used to configure the
agdc query. Items such as `min_date`, `max_date` and `satellites` are input
directly into the agdc query. The option `chunks_per_node` is used to configure
how many fractional cover files are to be produced on a given node. Depending
on the size of the query, this option can increase or decrease the processing
time.
The section titled `[pbs]` is used to configure the PBS submission script.
Here the user can provide items such as the `queue` type, `project` code,
`walltime` and any required modules.

Procedure
---------
To run the Fractional Cover job, follow these steps:

1. Specify agdc query
2. Execute the query -> `python query_agdc.py`
4. Review results
5. Handling errors

Specify inputs and outputs
--------------------------
The fc.cfg configuration file can be edited to control how the agdc fractional cover workflow will operate.
The fc.configuration file covers both the agdc and the scene level workflows.
The section titled `[work]` is used by both the agdc and the scene level workflows. Currently there are only
two parameters: `x_tile_size` & `y_tile_size`. These tell the program the spatial tile size to use for I/O and
processing, when the entire image is too large to fit into memory all at once.
The section titled `[scene]` caters for the scene level workflow. Here the final output filename format can
be defined, as well as the range of wavelengths for determining which bands of an acquisition are to be used
for deriving the fractional cover product.
The section titled `[agdc]` caters for the agdc workflow. Here the user specifies the output directory, and the path
for individual log files. Database query options include `min_date`, `max_date`, `satellites`, `cell_xmin`, `cell_ymin`, 
`cell_xmax`, `cell_ymax`. The result of the database query will be saved to a file defined by the `query_result` parameter.
The section titled `[core]` points to the logging configuration file that determines how luigi will log the processing.
`gaip` comes with a default log file named logging.cfg. Here the logging level such as `ERROR`, `WARNING`, `DEBUG` can be defined.

Set job resources
-----------------
The current version of the Fraction Cover program requires the following resources to process one input scene:

+------------------------+---------------------------+
| Resource               | Quantity                  |
+========================+===========================+
| CPUs                   + 1                         |
+------------------------+---------------------------+
| Wallclock time         | 27 minutes                 |
+------------------------+---------------------------+
| Memory                 | 4 GBytes                  |
+------------------------+---------------------------+
| Job file system        | 1 MB                      |
| (solid state disk)     | (for log files)           |
+------------------------+---------------------------+
| Luigi Workers          | 1                         |
+------------------------+---------------------------+

Luigi parallel processing
-------------------------
Luigi employs multiple CPUs to run many instances of the fractional cover program at the same time, within the 
context of a single PBS job. This is quite different from the previous way of doing FC processing.
Operations staff are required to adjust the PBS job resource parameters by editing the submit_FC.sh 
script so that the input workload can be processed efficiently and in a reasonable timeframe.

Scaling up
----------
Where there are many input scenes to processed additional resources need to be allocated to the PBS job to 
allow processing to complete in a reasonable (wallclock) timeframe. The following table provides a guide to the 
resources that should be allocated.

+---------+----+----+-----+-------+--------+---------+---------+
| Scenes  | 1  | 10 | 100 | 1,000 | 10,000 | 100,000 | 200,000 |
+=========+====+====+=====+=======+========+=========+=========+
| CPUs    | 1  | 16 | 32  |  128  |  1024  |   4096  |  4096   |
+---------+----+----+-----+-------+--------+---------+---------+
| Nodes   | 1  | 1  |  2  |   16  |  64    |   256   |   256   |
+---------+----+----+-----+-------+--------+---------+---------+
| Wall    | 30 | 30 | 60  | 8 hrs | 10 hrs | 24 hrs  | 24 hrs  |
| clock   |    |    |     |       |        |         |         |
+---------+----+----+-----+-------+--------+---------+---------+
| Memory  | 8  | 32 | 64  |  512  |  2048  |  8192   |  8192   |
| (GB)    |    |    |     |       |        |         |         |
+---------+----+----+-----+-------+--------+---------+---------+
| Job FS  | 1  | 1  | 1   |  1    |   100  |   500   |  1GB    |
| (GB)    |    |    |     |       |        |         |         |
+---------+----+----+-----+-------+--------+---------+---------+
| Luigi   |    |    |     |       |        |         |         |
| Worker  |  1 |  8 |  8  |   8   |    8   |     8   |    8    |
| / node  |    |    |     |       |        |         |         |
+---------+----+----+-----+-------+--------+---------+---------+

Key constraints to note:

1. A maximum of 8 Luigi workers per node is allowed (8 workers X 4GB per worker = 32GB = max memory available per node)
2. For any production workload, NCPUS (number of CPUs) should always be a multiple of 16 (so that whole Nodes will be allocated to the PBS job)

Specify scale of job
--------------------
Edit the following two lines in the submit_FC.sh script file

.. code-block:: bash

 WALLCLOCK='01:30:00'    # <---- Change this depending on workload (see SOP)
 NODES=2                 # <---- Change this depending on workload (see SOP)

using the information above as a guide to the number of CPUs and wallclock time required to process the current workload.

Submit and monitor job
----------------------
Once the job script submit_FC.sh has been edited and the correct entries inserted, run the script so that the PBS job will be submitted:

.. code-block:: bash

 ./submit_agdc_FC.sh

Check that the job is queued and, after some short delay is executing

.. code-block:: bash

 nqstat | grep run_agdc_fc

Review Results
--------------
Reviewing the results involves:

1. Checking output files
2. Reviewing exit code of PBS job
3. Check PBS standard error file
4. Checking Luigi Worker Logs

Checking output files
---------------------
Check that the expected number of pixel quality files have been written to the output directory.

Checking log files
------------------
Review the files in the log directory. An example is shown below.

.. code-block:: bash

 run_FC_raijin4_4596.stderr  run_FC_r82_7646.log   run_FC_r83_29470.log
 run_FC_raijin4_4596.stdout  run_FC_r82_7648.log   run_FC_r83_29472.log
 run_FC_r82_7642.log         run_FC_r83_29466.log
 run_FC_r82_7644.log         run_FC_r83_29468.log

Three types of files are present,  job STDOUT, job STDERR (recognised by the familiar file suffix). The remaining files (with the .log suffix) are Luigi Worker log files.
Reviewing exit code of PBS job
The job STDOUT file should be inspected to ensure that the Exit Status: 0 message is present as shown below. Any other status should be investigated.


.. code-block:: bash

 ============================================================================
               Resource Usage on 2014-12-24 11:16:38.991116:
 JobId:  8538551.r-man2
 Project: v10
 Exit Status: 0 (Linux Signal 0)
 Service Units: 4.34
 NCPUs Requested: 32                             NCPUs Used: 32
                                                 CPU Time Used: 00:23:20
 Memory Requested: 65536mb                       Memory Used: 25536mb
                                                 Vmem Used: 33966mb
 Walltime requested: 01:30:00                    Walltime Used: 00:08:08
 jobfs request: 200mb                            jobfs used: 2mb
 =============================================================================

Check PBS standard error and output files
-----------------------------------------
Both the job STDERR file and the STDOUT file in the logs directory should be checked for errors and warnings. They should be free of errors and can be checked using:

.. code-block:: bash

 cd <log directory>
 grep ERROR *.std*
 grep WARN *.std*

Look carefully at these files particularly if the job terminated with a non-zero exit status (see previous section)

Check Luigi Worker Logs
-----------------------
Each Luigi Work (up to 8 per Node) will produce a log file recording all events that the worker has encountered. A set of typical work log files looks like:


.. code-block:: bash

 run_fc_r2393_2767.log   run_fc_r2942_11499.log  run_fc_r2944_31469.log
 run_fc_r2393_2769.log   run_fc_r2942_11501.log  run_fc_r2944_31471.log
 run_fc_r2393_2771.log   run_fc_r2942_11503.log  run_fc_r2945_27573.log
 run_fc_r2393_2773.log   run_fc_r2942_11505.log  run_fc_r2945_27575.log
 run_fc_r2393_2775.log   run_fc_r2942_11507.log  run_fc_r2945_27577.log
 run_fc_r2393_2777.log   run_fc_r2944_31457.log  run_fc_r2945_27579.log
 run_fc_r2393_2779.log   run_fc_r2944_31459.log  run_fc_r2945_27581.log
 run_fc_r2393_2781.log   run_fc_r2944_31461.log  run_fc_r2945_27583.log
 run_fc_r2942_11493.log  run_fc_r2944_31463.log  run_fc_r2945_27585.log
 run_fc_r2942_11495.log  run_fc_r2944_31465.log  run_fc_r2945_27587.log
 run_fc_r2942_11497.log  run_fc_r2944_31467.log

Each log file includes the host name of the Node on which the job ran (e.g. “r2393”)
as will as the process ID of the worker on that host (e.g. “2777”)

Check for error messages in these file by:

.. code-block:: bash

 cd <log directory>
 grep ERROR *.log
 grep WARN *.log

Investigate any errors found by this process.

Handling errors
---------------
It is impossible to predict the various types of error that may occur during a processing run. Evaluate each error and decide on the appropriate actions to fix the error.

As a general rule, Fractional Cover jobs are completely re-runnable. So once errors have been fixed (and offending data files have been fixed or deleted), simply re-submit the Fractional Cover job and allow it to re-run.

When a Fractional Cover job is re-run, Luigi ensures that steps that previously completed without error will not be re-run. This property allows a strategy of “run, fix and rerun” to be employed until the workload has been fully processed.



Appendix A - Scene input file formats
-------------------------------------

Scene input data (NBAR) used by the Fractional Cover job are stored in directories, one scene per directory. The directory names subscribe to the following convention demonstrated here by example.
 
 Directory name: ``LS5_TM_NBAR_P54_GANBAR01-002_092_086_20090115``

The name is broken into fields using the underscore “_” character as a field delimiter. The following table describes the fields:


+---------------------------+--------------------+------------------------------------------------+
| Field                     | Example            |  Comment                                       |
+===========================+====================+================================================+
| Satellite                 | LS5                |                                                |
+---------------------------+--------------------+------------------------------------------------+
| Sensor                    | TM                 |                                                |
+---------------------------+--------------------+------------------------------------------------+
| Product                   | NBAR               |                                                |
+---------------------------+--------------------+------------------------------------------------+
| Product ID                | P54                |                                                |
+---------------------------+--------------------+------------------------------------------------+
| Product code and version  | GANBAR01           |                                                |
+---------------------------+--------------------+------------------------------------------------+
| Station ID                | 002                |                                                |
+---------------------------+--------------------+------------------------------------------------+
| Path                      | 092                |                                                |
+---------------------------+--------------------+------------------------------------------------+
| Row                       | 086                |                                                |
+---------------------------+--------------------+------------------------------------------------+
| Acquisition Date          | 20090205           |                                                |
+---------------------------+--------------------+------------------------------------------------+
