#PBS -P v10
#PBS -q normal
#PBS -l walltime=15:00:00,ncpus=1,mem=6GB,jobfs=1GB
#PBS -l wd
#PBS -me

module use /projects/u46/opt/modules/modulefiles
module load gcc/5.2.0 core

python compare_files.py --reference_dir /g/data/v10/testing_ground/jps547/test-mod-base2/output/ --test_dir /g/data/v10/testing_ground/jps547/test-mod-input3/output --files image-files.txt --scenes scenes.txt
