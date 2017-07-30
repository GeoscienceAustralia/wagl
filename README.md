# GAIP

Geoscience Australia Image Processor

## Instructions (quickstart)

Note this is an older branch of gaip, which is preserved as an interim solution for maintaining a downstream product (WOfS).

These instructioned are tailored to NCI/Raijin system.

First, open terminal (ssh raijin).

Obtain and prepare gaip software. E.g., git checkout and make, and configure as a loadable module. Prepare the environment by loading modules containing dependencies:

``` 
module use /g/data/v10/private/modules/modulefiles
module load gcc/5.2.0 core
```

If there are multiple versions of gaip installed, ensure that the correct scripts will be invoked. For example, work in a directory containing `pq_script_generator.py` and `run_PQ.pbs`, checking that the latter pbs script contains lines which load *the wofs pq variant of* the gaip module (on the worker nodes). Avoid actually loading the gaip module in the interactive terminal to ensure no ambiguity between script versions.

`cd /g/data/v10/testing_ground/4.2.5-pq-wofs`

The necessary inputs are the level1 and NBAR scenes (both are necessary). There is also one ancillary product, the land-sea rasters, and its configuration is hard-coded.

The `python pq_script_generator.py` will produce a bash script to kick off the job, and will immediately source (execute) that script unless the `--test` argument is supplied (e.g., permitting modification of PBS settings for just one job). 

The script only handles one platform, year, and month at a time. It can easily be scripted to queue e.g. an entire year:

`for p in ls7 ls8; do for m in {1..12}; do echo $p $m; done; done`

`$ for i in {1..7}; do python pq_script_generator.py --l1t_path /g/data/v10/reprocess/ --nbar_path /short/v10/scenes/nbar-scenes-tmp --out_path /g/data/v10/testing_ground/4.2.5-pq-wofs/outputs --year 2016 --month $i --sensor ls7 ; done`

`for p in ls7 ls8 ; do for i in {1..12}; do python pq_script_generator.py --l1t_path /g/data/v10/reprocess/ --nbar_path /short/v10/scenes/nbar-scenes-tmp --out_path /g/data/v10/testing_ground/4.2.5-pq-wofs/outputs --year 2015 --month $i --sensor $p ; done ; done`

Any job will result in log files and, for each scene, a .completed file (small), a packaged pqa output (e.g. roughly 8 gig/month), and a folder of intermediate outputs (e.g. >100gig/month; presumably uncompressed rasters). The latter should be cleaned up for sake of storage space. The .log files are worth checking for signal that any issues were detected. It is also easy to confirm the count of successful completions.

```
/g/data/v10/testing_ground/4.2.5-pq-wofs
$ find outputs/*/PQ/*/*/logs -maxdepth 1 -mindepth 1 -name *.log | xargs grep 'progress looks :)' -L | tee /dev/tty | wc --lines
$ ls outputs/*/PQ/*/*/output/*.completed | wc --lines
$ ls outputs/*/PQ/*/*/output/pqa/* -d | wc --lines
$ ls outputs/*/PQ/*/*/output/*.completed | sed s:\.completed$:: | xargs rm -rf
$ ls outputs/*/PQ/*/*/output/LS* -d | grep -v .completed
```

Certain naming conventions are expected by gaip. If the input folders contain exceptions, then a staging area should be used as input and be populated with a cleaned set of symbolic links. For example, applying a rule of only processing definitive (not predictive) datasets:

`for file in /g/data/v10/reprocess_interim/ls7/level1/2016/12/L*OTH*; do ln -s $file $(basename $file); done && rm *PREDICT*`

```
cd /g/data/v10/reprocess
ls ls5/level1/*/* -d --color=never | xargs -I {} mkdir -p /g/data/v10/testing_ground/4.2.5-pq-wofs/staging/{}
find ls5/level1 -maxdepth 3 -mindepth 3 -iname LS\*OTH\* | xargs -I {} ln -s $PWD{} /g/data/v10/testing_ground/4.2.5-pq-wofs/staging/{}
cd /g/data/v10/testing_ground/4.2.5-pq-wofs/staging
rm ls5/level1/*/*/*PREDICT*
for p in ls7 ls8 ; do mkdir -p $p ; ln -s /g/data/v10/reprocess/$p/level1 $p/level1 ; done
```

Ideally, if the input data are spread over multiple locations, they could all be unified together under a single consistent facade, the staging area. Consistently also facilitates greater automation.

```
for y in {1999..2009} ; do  for i in {1..12}; do python pq_script_generator.py --l1t_path /g/data/v10/testing_ground/4.2.5-pq-wofs/staging/ --nbar_path /short/v10/scenes/nbar-scenes-tmp --out_path /g/data/v10/testing_ground/4.2.5-pq-wofs/outputs --year $y --month $i --sensor ls7 ; done ; done

```

May be an issue if too many jobs are submitted simultaneously. This can be resolved.

```

for y in {1986..2011} ; do  for i in {1..12}; do python pq_script_generator.py --l1t_path /g/data/v10/testing_ground/4.2.5-pq-wofs/staging/ --nbar_path /short/v10/scenes/nbar-scenes-tmp --out_path /g/data/v10/testing_ground/4.2.5-pq-wofs/outputs --year $y --month $i --sensor ls5 ; while [ $( qstat | grep $USER | wc --lines ) -gt 100 ] ; do sleep 30 ; done ; done ; done

```

