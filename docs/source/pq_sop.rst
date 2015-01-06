
============================
Standard Operating Procedure
============================

Pixel Quality
-------------

Purpose
-------
A pixel quality program (run_pq.py) takes Landsat L1T and NBAR scene pairs and produces a Pixel Quality scene mask as output. 

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

* submit_PQ.sh – a convenience script to bootstrap the PBS job
* run_PQ.pbs – main PBS job scipt
* `Python <http://www.python.org/>`_

You can view these scripts by clicking on the links above.

Procedure
---------
To run the Pixel Quality job, follow these steps:

# Specify inputs and outputs
# Set job resources
# Submit and monitor job
# Review results
# Handling errors

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





