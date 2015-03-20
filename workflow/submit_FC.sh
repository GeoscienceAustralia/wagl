#!/bin/bash
#
# Start a FC run by submitting the run_FC.pbs script to the PBS queue
# This script initialises all the appropriate logging and job parameters
#
# Edit this file to get your run parameters (Inputs, outputs and resources)
#
# usage: ./submit_FC.sh 
#
#================================================================================
umask 007

# INPUTS (change any of these)
# ======

NBAR_PATH=/g/data1/v10/projects/Luigi_work_flow_test/NBAR

# OUTPUTS (change any of these)
# =======

OUTPUT_ROOT=/short/v10/$USER/FC_luigi_new_workflow_testing
OUTPUT_PATH=$OUTPUT_ROOT/output
LOG_PATH=$OUTPUT_ROOT/logs

# RESOURCES
# =========

PROJECT='v10'           # <---- Change this as required
QUEUE='normal'          # <---- Change this as required
WALLCLOCK='01:30:00'    # <---- Change this depending on workload (see SOP)
NODES=2                 # <---- Change this depending on workload (see SOP)

let MEM_GB=32*$NODES    # don't touch
let NCPUS=16*$NODES     # don't touch
let JOBFS_MB=100*$NODES # don't touch

# JOB SPECS
# =========

MAIL_LIST="steven.ring@ga.gov.au,smr@southsky.com.au"

# Create ouput directory

rm -r ${OUTPUT_ROOT}
mkdir -p ${OUTPUT_PATH}
mkdir -p ${LOG_PATH}


#================================================================================
HOSTNAME=`hostname`
PID=$$
LSPEC="walltime=${WALLCLOCK},ncpus=${NCPUS},mem=${MEM_GB}GB,jobfs=${JOBFS_MB}MB"
echo $LSPEC

# Submit the job

qsub -v NBAR_PATH=$NBAR_PATH,OUTPUT_PATH=$OUTPUT_PATH,LOG_PATH=$LOG_PATH \
 -M $MAIL_LIST \
 -e $LOG_PATH/run_FC_${HOSTNAME}_${PID}.stderr \
 -o $LOG_PATH/run_FC_${HOSTNAME}_${PID}.stdout \
 -l $LSPEC \
 -q $QUEUE \
 run_FC.pbs
