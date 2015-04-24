#!/bin/bash
#PBS -P v10
#PBS -q normal
#PBS -l walltime=00:50:00,mem=64GB,ncpus=32
#PBS -l wd
#PBS -me
#PBS -M joshua.sixsmith@ga.gov.au

module load python/2.7.6

module use /projects/u46/opt/modules/modulefiles
module load gdal/1.10.1
module load gaip
module load agdc-api


mpirun -n 32 python dc_fc.py
