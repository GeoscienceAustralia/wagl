#!/bin/bash
#PBS -P v10
#PBS -q normal
#PBS -l walltime=00:30:00,ncpus=1,vmem=10700MB
#PBS -wd
#PBS -me
#PBS -M Simon.Knapp@ga.gov.au

module load geotiff/1.3.0
module load python/2.7.2
module load hdf4/4.2.6_2012
#module load gdal/1.9.0_HDF5
module load gdal/1.9.0
module load proj

module use /short/v10/mitch-sys/opt/modules --append
module load py-dev-tools

export IMAGEPROCESSOR_ROOT=$(readlink -f ${0%/*}/..)

PYTHONPATH=$PYTHONPATH:$IMAGEPROCESSOR_ROOT:/short/v10/nbar/pyephem-3.7.5.1/lib/python2.7/site-packages

# NOTE: if you are going to test Fuqin's code with a new DSM, you should run the first of these.

python -m unittest -v test_tc.ClipDSMTestCase.test_fuqin_test_data_smooth
#python -m unittest -v test_tc.RunCastShadowTestCase.test_fuqin_data
#python -m unittest -v test_tc.RunSlopeAngleTestCase.test_fuqin_data
python -m unittest -v test_tc.RunTerrainCorrectionTestCase.test_fuqin_data

