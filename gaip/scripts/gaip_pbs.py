#!/usr/bin/env python

"""
PBS submission scripts.
"""

import os
from os.path import join as pjoin, dirname, exists
import subprocess
import uuid
import argparse

from gaip.tiling import scatter


PBS_TEMPLATE = ("""#!/bin/bash
#PBS -P {project}
#PBS -q {queue}
#PBS -l walltime={hours}:00:00,mem={memory}GB,ncpus={ncpus}
#PBS -l wd
#PBS -me
#PBS -M {email}

source {env}

{daemon}

luigi --module gaip.{workflow} ARD --model {model} --level1-list {scene_list} --outdir {outdir} --workers 16{scheduler} --vertices '{vertices}' --method {method}{pq}
""")

DSH_TEMPLATE = ("""#!/bin/bash
#PBS -P {project}
#PBS -q {queue}
#PBS -l walltime={hours}:00:00,mem={memory}GB,ncpus={ncpus}
#PBS -l wd
#PBS -me
#PBS -M {email}

FILES=({files})

DAEMONS=({daemons})

OUTDIRS=({outdirs})

for i in "${{!FILES[@]}}"; do
  X=$(($i+1))
  pbsdsh -n $((16 *$X)) -- bash -l -c "source {env}; ${{DAEMONS[$i]}}; luigi \\
    --module gaip.{workflow} ARD \\
    --model {model} \\
    --level1-list ${{FILES[$i]}} \\
    --outdir ${{OUTDIRS[$i]}} \\
    --workers 16 \\
    --vertices '{vertices}' \\
    --method {method}{pq}" &
done;
wait
""")

FMT1 = 'level1-scenes-{jobid}.txt'
FMT2 = '{model}-ard-{jobid}.bash'
DAEMON_FMT = 'luigid --background --logdir {}'


def _submit_dsh(scattered, vertices, model, method, pq, batchid, batch_logdir,
                batch_outdir, project, queue, memory, ncpus, hours, email, env,
                test, workflow):
    """Submit a single PBSDSH formatted job."""
    files = []
    daemons = []
    outdirs = []
    jobids = []

    # setup each block of scenes for processing
    for block in scattered:
        jobid = uuid.uuid4().hex[0:6]
        jobids.append(jobid)
        jobdir = pjoin(batch_logdir, 'jobid-{}'.format(jobid))
        job_outdir = pjoin(batch_outdir, 'jobid-{}'.format(jobid))

        if not exists(jobdir):
            os.makedirs(jobdir)

        if not exists(job_outdir):
            os.makedirs(job_outdir)

        # write the block of scenes to process
        out_fname = pjoin(jobdir, FMT1.format(jobid=jobid))
        with open(out_fname, 'w') as src:
            src.writelines(block)

        files.append(out_fname)
        daemons.append(DAEMON_FMT.format(jobdir))
        outdirs.append(job_outdir)

    files = ['"{}"\n'.format(f) for f in files]
    daemons = ['"{}"\n'.format(f) for f in daemons]
    outdirs = ['"{}"\n'.format(f) for f in outdirs]

    pbs = DSH_TEMPLATE.format(project=project, queue=queue, hours=hours,
                              memory=memory, ncpus=ncpus,
                              email=email, files=''.join(files), env=env,
                              daemons=''.join(daemons), model=model, pq=pq,
                              outdirs=''.join(outdirs), vertices=vertices,
                              method=method, workflow=workflow)

    out_fname = pjoin(batch_logdir, FMT2.format(model=model, jobid=batchid))
    with open(out_fname, 'w') as src:
        src.write(pbs)

    print("Job ids:\n{}".format(jobids))
    if test:
        print("qsub {}".format(out_fname))
    else:
        os.chdir(dirname(out_fname))
        subprocess.call(['qsub', out_fname])


def _submit_multiple(scattered, vertices, model, method, pq, batchid,
                     batch_logdir, batch_outdir, project, queue, memory, ncpus,
                     hours, email, local_scheduler, env, test, workflow):
    """Submit multiple PBS formatted jobs."""
    print("Executing Batch: {}".format(batchid))
    # setup and submit each block of scenes for processing
    for block in scattered:
        jobid = uuid.uuid4().hex[0:6]
        jobdir = pjoin(batch_logdir, 'jobid-{}'.format(jobid))
        job_outdir = pjoin(batch_outdir, 'jobid-{}'.format(jobid))

        if not exists(jobdir):
            os.makedirs(jobdir)

        if not exists(job_outdir):
            os.makedirs(job_outdir)

        # local or central scheduler
        if local_scheduler:
            daemon = ''
            scheduler = ' --local-scheduler'
        else:
            daemon = DAEMON_FMT.format(jobdir)
            scheduler = ''

        out_fname = pjoin(jobdir, FMT1.format(jobid=jobid))
        with open(out_fname, 'w') as src:
            src.writelines(block)

        pbs = PBS_TEMPLATE.format(project=project, queue=queue,
                                  hours=hours, memory=memory, ncpus=ncpus,
                                  email=email, env=env, daemon=daemon,
                                  model=model, pq=pq, scene_list=out_fname,
                                  outdir=job_outdir, scheduler=scheduler,
                                  vertices=vertices, method=method,
                                  workflow=workflow)

        out_fname = pjoin(jobdir, FMT2.format(model=model, jobid=jobid))
        with open(out_fname, 'w') as src:
            src.write(pbs)

        if test:
            print("Mocking... Submitting Job: {} ...Mocking".format(jobid))
            print("qsub {}".format(out_fname))
        else:
            os.chdir(dirname(out_fname))
            print("Submitting Job: {}".format(jobid))
            subprocess.call(['qsub', out_fname])


def run(level1, vertices='(5, 5)', model='standard', method='linear',
        pixel_quality=False, outdir=None, logdir=None, env=None, nodes=10,
        project=None, queue='normal', hours=48,
        email='your.name@something.com', local_scheduler=False, dsh=False,
        test=False, singlefile=False):
    """Base level program."""
    with open(level1, 'r') as src:
        scenes = src.readlines()

    # scattered = scatter(filter_scenes(scenes), nodes)
    scattered = scatter(scenes, nodes)

    batchid = uuid.uuid4().hex[0:10]
    batch_logdir = pjoin(logdir, 'batchid-{}'.format(batchid))
    batch_outdir = pjoin(outdir, 'batchid-{}'.format(batchid))

    # compute resources
    memory = 32 * nodes
    ncpus = 16 * nodes

    pq = ' --pixel-quality' if pixel_quality else ''

    workflow = 'singlefile_workflow' if singlefile else 'multifile_workflow'

    if test:
        print("Mocking... Submitting Batch: {} ...Mocking".format(batchid))
    else:
        print("Submitting Batch: {}".format(batchid))

    if dsh:
        _submit_dsh(scattered, vertices, model, method, pq, batchid,
                    batch_logdir, batch_outdir, project, queue, memory, ncpus,
                    hours, email, env, test, workflow)
    else:
        _submit_multiple(scattered, vertices, model, method, pq, batchid,
                         batch_logdir, batch_outdir, project, queue, memory,
                         ncpus, hours, email, local_scheduler, env, test,
                         workflow)


def _parser():
    """ Argument parser. """
    description = ("Equally partition a list of scenes in n nodes and submit "
                   "into the PBS queue. Optionally submit as multiple "
                   "jobs into the PBS queue, or as a single job "
                   "and executed using PBSDSH.")
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--level1-list", help="The input level1 scene list.",
                        required=True)
    parser.add_argument("--vertices", default="(5, 5)", type=str,
                        help=("Number of vertices to evaluate the radiative "
                              "transfer at. JSON styled string is required."))
    parser.add_argument("--model", default="standard",
                        help=("The type of ARD workflow to invoke, "
                              "eg standard, nbar, sbt."))
    parser.add_argument("--method", default="shear",
                        help=("The interpolation method to invoke, "
                              "eg linear, shear, rbf."))
    parser.add_argument("--pixel-quality", action='store_true',
                        help=("Whether to run the pixel quality workflow, "
                              "if applicable, or not."))
    parser.add_argument("--outdir", help="The base output directory.",
                        required=True)
    parser.add_argument("--logdir", required=True,
                        help="The base logging and scripts output directory.")
    parser.add_argument("--env", help="Environment script to source.",
                        required=True)
    parser.add_argument("--nodes", type=int,
                        help="The number of nodes to request.", default=10)
    parser.add_argument("--project", help="Project code to run under.",
                        required=True)
    parser.add_argument("--queue", default='normal',
                        help=("Queue to submit the job into, "
                              "eg normal, express."))
    parser.add_argument("--hours", help="Job walltime in hours.", default=48)
    parser.add_argument("--email", help="Notification email address.",
                        default="your.name@something.com")
    parser.add_argument("--local-scheduler", action='store_true',
                        help=("Use a local scheduler instead of a central "
                              "scheduler."))
    parser.add_argument("--dsh", help="Run using PBS Distributed Shell.",
                        action='store_true')
    parser.add_argument("--singlefile", action='store_true',
                        help=("Run gaip using the single-file workflow. "
                              "Default is to output multiple files."))
    parser.add_argument("--test", action='store_true',
                        help=("Test job execution (Don't submit the job to "
                              "the PBS queue)."))
    return parser


def main():
    """ Main execution. """
    parser = _parser()
    args = parser.parse_args()
    run(args.level1_list, args.vertices, args.model, args.method,
        args.pixel_quality, args.outdir, args.logdir, args.env, args.nodes,
        args.project, args.queue, args.hours, args.email, args.local_scheduler,
        args.dsh, args.test, args.singlefile)


if __name__ == '__main__':
    main()
