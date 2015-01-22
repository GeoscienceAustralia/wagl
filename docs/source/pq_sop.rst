
===========================================
Standard Operating Procedure (Pixel Quality)
===========================================

Pixel Quality
-------------

Purpose
-------
A pixel quality program (pq.py_) takes Landsat L1T and NBAR scene pairs and produces a Pixel Quality scene mask as output. 

.. _pq.py: https://github.com/GeoscienceAustralia/ga-neo-landsat-processor/blob/develop/workflow/pq.py

This document describes how to run the Pixel Quality job via PBS script to process production data. 

The Pixel Quality PBS job uses Luigi to run multiple instances of the pixel quality program in parallel. By adjusting the system resources used by the PBS job, hundreds (or thousands) of L1T input files may be processed in a short period of time.

Schematic
---------
The diagram below shows the basic operation of the Pixel Quality Job

.. image:: https://raw.githubusercontent.com/smr547/ga-neo-landsat-processor/sops/docs/source/diagrams/pq.png?token=AEIJ1yV9RFDWooW_7muuFban7sdT91G9ks5UtXyJwA%3D%3D

Key elements are:

* shell script submit_PQ.sh is used to submit the run_PQ.pbs script to the PBS job queue
* run_PQ.pbs reads scene data from the L1T (product code “OTH”) and * NBAR input directories (product code “NBAR”)
* data are also read from the ancillary data directory
* output scenes (Pixel Quality scenes) are written to the PQ Output Directory
* log files are written to the logs directory

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

PBS Scripts
-----------
Two scripts are used to run Pixel Quality processing

* submit_PQ.sh_ – a convenience script to bootstrap the PBS job
* run_PQ.pbs_ – main PBS job scipt

.. _submit_PQ.sh: https://github.com/GeoscienceAustralia/ga-neo-landsat-processor/blob/develop/workflow/submit_PQ.sh
.. _run_PQ.pbs: https://github.com/GeoscienceAustralia/ga-neo-landsat-processor/blob/develop/workflow/run_PQ.pbs

You can view these scripts by clicking on the links above.

Procedure
---------
To run the Pixel Quality job, follow these steps:

1. Specify inputs and outputs
2. Set job resources
3. Submit and monitor job
4. Review results
5. Handling errors

Specify inputs and outputs
--------------------------
Edit the submit_PQ.sh script file and specify the following paths

+---------------+----------------------------------------------------------------------------+
| Variable      | Description                                                                |
+===============+============================================================================+
| L1T_PATH      | Path to directory containing L1T files. Files in this directory must       |
|               | comply with the file naming standards described in Appendix A              |
+---------------+----------------------------------------------------------------------------+
| NBAR_PATH     | Path to directory containing NBAR files. Files in this directory must      |
|               | comply with the file naming standards described in Appendix A              |
+---------------+----------------------------------------------------------------------------+
| LAND_SEA_PATH | Path to directory containing the LAND/SEA mask Geotiff files. Refer to     |
|               | Appendix B                                                                 |
+---------------+----------------------------------------------------------------------------+
| OUTPUT_PATH   | Path to the directory where output files are to be placed. The PQA files   |
|               | output by this job will have a filename compliant with the standard        |
|               | described in Appendix A. The product code will be “PQA”. The remaining     |
|               | components of the filename will be those of the L1T input filename.        |
+---------------+----------------------------------------------------------------------------+
| LOG_PATH      | Path to the directory where log files will be placed.                      |
+---------------+----------------------------------------------------------------------------+

Set job resources
-----------------
The current version of the Pixel Quality program requires the following resources to process one input scene:

+------------------------+---------------------------+
| Resource               | Quantity                  |
+========================+===========================+
| CPUs                   + 1                         |
+------------------------+---------------------------+
| Wallclock time         | 4 minutes                 |
+------------------------+---------------------------+
| Memory                 | 8 GBytes                  |
+------------------------+---------------------------+
| Job file system        | 1 MB                      |
| (solid state disk)     | (for log files)           |
+------------------------+---------------------------+
| Luigi Workers          | 1                         |
+------------------------+---------------------------+

Luigi parallel processing
-------------------------
Luigi employs multiple CPUs to run many instances of the pixel quality program at the same time, within the 
context of a single PBS job. This is quite different from the previous way of doing PQ processing.
Operations staff are required to adjust the PBS job resource parameters by editing the submit_PQ.sh 
script so that the input workload can be processed efficiently and in a reasonable timeframe.

Scaling up
----------
Where there are many input scenes to processed additional resources need to be allocated to the PBS job to 
allow processing to complete in a reasonable (wallclock) time. The following table provides a guide to the 
resources that should be allocated.

+---------+----+----+-----+-------+--------+---------+---------+---------+
| Scenes  | 1  | 10 | 100 | 1,000 | 10,000 | 100,000 | 500,000 | 500,000 |
+=========+====+====+=====+=======+========+=========+=========+=========+
| CPUs    | 1  | 16 | 32  |  128  |   512  |   1024  |  3072   |   4096  |
+---------+----+----+-----+-------+--------+---------+---------+---------+
| Nodes   | 1  | 1  |  2  |    8  |   32   |    64   |   192   |   256   |
+---------+----+----+-----+-------+--------+---------+---------+---------+
| Wall    | 5  | 5  | 50  |  125  |  312   |  1562   |  2604   |  1953   |
| clock   |    |    |     | 2 hrs |  5 hrs | 26 hrs  |  43 hrs | 33 hrs  |
+---------+----+----+-----+-------+--------+---------+---------+---------+
| Memory  | 8  | 32 | 64  |  256  |  1024  |  2048   |  6144   |  8192   |
| (GB)    |    |    |     |       |        |         |         |         |
+---------+----+----+-----+-------+--------+---------+---------+---------+
| Job FS  | 1  | 1  | 1   |  1    |    10  |   100   |  500    |  500    |
| (GB)    |    |    |     |       |        |         |         |         |
+---------+----+----+-----+-------+--------+---------+---------+---------+
| Luigi   |    |    |     |       |        |         |         |         |
| Worker  |  1 |  4 |  4  |   4   |    4   |     4   |    4    |    4    |
| / node  |    |    |     |       |        |         |         |         | 
+---------+----+----+-----+-------+--------+---------+---------+---------+

Key constraints to note:

1. A maximum of 4 Luigi workers per node is allowed (4 workers X 8GB per worker = 32GB = max memory available per node)
2. For any production workload, NCPUS (number of CPUs) should always be a multiple of 16 (so that whole Nodes will be allocated to the PBS job)

Specify scale of job
--------------------
Edit the following two lines in the submit_PQ.sh script file

.. code-block:: bash

 WALLCLOCK='01:30:00'    # <---- Change this depending on workload (see SOP)
 NODES=2                 # <---- Change this depending on workload (see SOP)

using the information above as a guide to the number of CPUs and wallclock time required to process the current workload.

Submit and monitor job
----------------------
Once the job script submit_PQ.sh has been edited and the correct entries inserted, run the script so that the PBS job will be submitted:

.. code-block:: bash

 ./submit_PQ.sh

Check that the job is queued and, after some short delay is executing

.. code-block:: bash

 nqstat | grep run_pq

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

 run_PQ_raijin4_4596.stderr  run_pq_r82_7646.log   run_pq_r83_29470.log
 run_PQ_raijin4_4596.stdout  run_pq_r82_7648.log   run_pq_r83_29472.log
 run_pq_r82_7642.log         run_pq_r83_29466.log
 run_pq_r82_7644.log         run_pq_r83_29468.log

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
Each Luigi Work (up to 4 per Node) will produce a log file recording all events that the worker has encountered. A set of typical work log files looks like:


.. code-block:: bash

 run_pq_r82_15108.log  run_pq_r83_11591.log  run_pq_r85_25905.log 
 run_pq_r82_15110.log  run_pq_r83_11593.log  run_pq_r85_25907.log 
 run_pq_r82_15112.log  run_pq_r84_3376.log   run_pq_r85_25909.log 
 run_pq_r82_15114.log  run_pq_r84_3378.log   run_pq_r85_25911.log 
 run_pq_r83_11587.log  run_pq_r84_3380.log 
 run_pq_r83_11589.log  run_pq_r84_3382.log

Each log file includes the host name of the Node on which the job ran (e.g. “r82”) as will as the process ID of the worker on that host (e.g. “15108”)

Check for error messages in these file by:

.. code-block:: bash

 cd <log directory>
 grep ERROR *.log
 grep WARN *.log

Investigate any errors found by this process.

Handling errors
---------------
It is impossible to predict the various types of error that may occur during PQ processing. Evaluate each error and decide on the appropriate actions to fix the error.

As a general rule, Pixel Quality jobs are completely re-runnable. So once errors have been fixed (and offending data files have been fixed or deleted), simply re-submit the Pixel Quality job and allow it to re-run.

When a Pixel Quality job is re-run, Luigi ensures that steps that previously completed without error will not be re-run. This property allows a strategy of “run, fix and rerun” to be employed until the workload has been fully processed.



Appendix A - Scene input file formats
-------------------------------------

Scene data (both L1T and NBAR) used by the Pixel Quality job are stored in directories, one scene per directory. The directory names subscribe to the following convention demonstrated here by example.
 
 Directory name: ``LS5_TM_NBAR_P54_GANBAR01-002_092_086_20090115``

The name is broken into fields using the underscore “_” character as a field delimiter. The following table describes the fields:


+---------------------------+--------------------+------------------------------------------------+
| Field                     | Example            |  Comment                                       |
+===========================+====================+================================================+
| Satellite                 | LS5                |                                                |
+---------------------------+--------------------+------------------------------------------------+
| Sensor                    | TM                 |                                                |
+---------------------------+--------------------+------------------------------------------------+
| Product                   | NBAR               |  "OTH" for L1T scenes                          |
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




Appendix B - Land/Sea data files
--------------------------------

Land sea raster files are currently stored in ``/g/data1/v10/eoancillarydata/Land_Sea_Rasters``

and have a filename format like ``WORLDzone57.tif``, where, in this case, 57 is the UTM zone.





