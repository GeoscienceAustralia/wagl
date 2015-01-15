Invoking the Geoscience Australia Image Processor
=================================================

The Geoscience Australia Image Processor is invoked as ...TODO 

.. raw:: html

	<pre>
	usage: process [-h] [-C CONFIG_FILE] [-i INPUT_PATH] [-o OUTPUT_PATH]
               	[-r OUTPUT_ROOT] [-w WORK_PATH] [-l LOG_FILE]
               	[--nbar-root _NBAR_ROOT] [--modtran-root _MODTRAN_ROOT]
               	[--process PROCESS] [--process-file PROCESS_FILE] [--resume]
               	[--debug] [--sequential] [--repackage] [-c CONSTRAINT_ID]
               	[-L LI_SOURCE_DESCRIPTION] [-p PURPOSE] [-d DATASET_ID]
	
	optional arguments:
  	-h, --help           show this help message and exit
  	-C CONFIG_FILE, --config CONFIG_FILE
                        NBAR configuration file
  	-i INPUT_PATH, --input INPUT_PATH
                        Input directory
  	-o OUTPUT_PATH, --output OUTPUT_PATH
                        Output directory
  	-r OUTPUT_ROOT, --output-root OUTPUT_ROOT
                        Output root directory
  	-w WORK_PATH, --work WORK_PATH
                        Work directory
  	-l LOG_FILE, --log LOG_FILE
                        Log File
  	--modtran-root _MODTRAN_ROOT
                        MODTRAN root directory
  	--process PROCESS     Process(es) to run
  	--process-file PROCESS_FILE
                        Process(es) to run defined in a text file
  	--resume              Resume execution from specified process
  	--debug               Create debug files in working directory
  	--sequential          Process all tasks sequentially (i.e. no parallel
                        processing)
  	--repackage           Repackage existing output only (i.e. do not re-create
                        output)
  	-c CONSTRAINT_ID, --constraintid CONSTRAINT_ID
                        constraint ID
  	-L LI_SOURCE_DESCRIPTION, --li-source-description LI_SOURCE_DESCRIPTION
                        Lineage source description
  	-p PURPOSE, --purpose PURPOSE
                        Purpose
  	-d DATASET_ID, --dataset-id DATASET_ID
                        Full dataset ID
	</pre>


To invoke the Pixel Quality workflow see :doc:`pq_sop </pq_sop>`.
