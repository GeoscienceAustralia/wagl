#!/bin/bash
#PBS -P v10
#PBS -q express
#PBS -l walltime=00:50:00,mem=32GB,ncpus=16
#PBS -l wd
#PBS -me
#PBS -M joshua.sixsmith@ga.gov.au

module load python/2.7.6

module use /projects/u46/opt/modules/modulefiles
module load gdal
module load gaip/test_aspect
module load agdc-api


mpirun -n $NCPUS python agdc_fc.py
