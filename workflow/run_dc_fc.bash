#!/bin/bash
#PBS -P v10
#PBS -q express
#PBS -l walltime=00:50:00,mem=64GB,ncpus=32
#PBS -l wd
#PBS -me
#PBS -M joshua.sixsmith@ga.gov.au

#ls

module load python/2.7.6

#python /home/547/jps547/git_repositories/my_code/Python/testing/test_parallel3.py

module use /projects/u46/opt/modules/modulefiles
#module load gdal/1.10.1
#module load ga_geotiff/1.3.0
module load gaip/test
module load agdc-api

#/projects/v10/geo_assess/tool/image-gverify_v0.11 -b /g/data/v10/eoancillarydata/GCP/GLS2000_GCP_SCENE/wrs2/101/078/p101r078_7dt20000911_z53_50.tif -m /g/data1/v10/testing_ground/jps547/gqa7/output/LS8_OLITIRS_GQA_P51_GALPGS01-015_101_078_20141012/test_warp.tif  -w /g/data1/v10/testing_ground/jps547/gqa7/output/LS8_OLITIRS_GQA_P51_GALPGS01-015_101_078_20141012 -l /g/data1/v10/testing_ground/jps547/gqa7/output/LS8_OLITIRS_GQA_P51_GALPGS01-015_101_078_20141012 -o /g/data1/v10/testing_ground/jps547/gqa7/output/LS8_OLITIRS_GQA_P51_GALPGS01-015_101_078_20141012 -p 5 -n 4 -nv 0 -c -.75 -r BI -cs 27 -g 22 -sws 28

#/projects/v10/geo_assess/tool/image-gverify_v0.11 -b /g/data/v10/eoancillarydata/GCP/GLS2000_GCP_SCENE/wrs2/105/080/p105r080_7dt20000721_z52_50.tif -m /g/data/v10/testing_ground/jps547/gqa7/output/LS5_TM_GQA_P51_GALPGS01-002_105_080_20100215/tmp/test_warp.tif -w /g/data/v10/testing_ground/jps547/gqa7/output/LS5_TM_GQA_P51_GALPGS01-002_105_080_20100215/tmp/ -l /g/data/v10/testing_ground/jps547/gqa7/output/LS5_TM_GQA_P51_GALPGS01-002_105_080_20100215/tmp/ -o /g/data/v10/testing_ground/jps547/gqa7/output/LS5_TM_GQA_P51_GALPGS01-002_105_080_20100215/tmp/ -p 5 -n 4 -nv 0 -c -.75 -r BI -cs 27 -g 22 -sws 28

######################/projects/v10/geo_assess/tool/image-gverify_v0.11 -b /g/data/v10/eoancillarydata/GCP/GLS2000_GCP_SCENE/wrs2/100/081/p100r081_7dt20001006_z53_50.tif -m /g/data1/v10/testing_ground/jps547/gqa7/output/LS5_TM_GQA_P51_GALPGS01-002_100_081_20100228/tmp/test_warp.tif -w /g/data1/v10/testing_ground/jps547/gqa7/output/LS5_TM_GQA_P51_GALPGS01-002_100_081_20100228/tmp/ -l /g/data1/v10/testing_ground/jps547/gqa7/output/LS5_TM_GQA_P51_GALPGS01-002_100_081_20100228/tmp -o /g/data1/v10/testing_ground/jps547/gqa7/output/LS5_TM_GQA_P51_GALPGS01-002_100_081_20100228/tmp -p 5 -n 4 -nv 0 -c -.75 -r BI -cs 27 -g 22 -sws 28

mpirun -n 32 python /home/547/jps547/git_repositories/sixy6e/ga-neo-landsat-processor/workflow/dc_fc.py
