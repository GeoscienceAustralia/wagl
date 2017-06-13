#!/usr/bin/env python

import os
from os.path import join as pjoin, dirname, exists
import subprocess
import uuid
import argparse

# from gaip.acquisition import acquisitions
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

luigi --module gaip.standard_workflow ARD --model {model} --level1-list {scene_list} --outdir {outdir} --workers 16{scheduler} --vertices '{vertices}' --method {method}
""")

DSH_TEMPLATE = ("""#!/bin/bash
#PBS -P {project}
#PBS -q {queue}
#PBS -l walltime={hours}:00:00,mem={memory}GB,ncpus={ncpus}
#PBS -l wd
#PBS -me
#PBS -M {email}

{files}

{daemons}

{outdirs}

for i in "${{!FILES[@]}}"; do
  X=$(($i+1))
  pbsdsh -n $((16 *$X)) -- bash -l -c "source {env}; ${{DAEMONS[$i]}}; luigi \\
    --module gaip.standard_workflow ARD \\
    --model {model} \\
    --level1-list ${{FILES[$i]}} \\
    --outdir ${{OUTDIRS[$i]}} \\
    --workers 16 \\
    --vertices '{vertices}' \\
    --method {method}" &
done;
wait
""")


def run(level1, vertices='(5, 5)', model='standard', method='linear',
        outdir=None, logdir=None, env=None, nodes=10, project=None,
        queue='normal', hours=48, email='your.name@something.com',
        local_scheduler=False, dsh=False, test=False):
    """Base level program."""
    with open(level1, 'r') as src:
        scenes = src.readlines()

    # scattered = scatter(filter_scenes(scenes), nodes)
    scattered = scatter(scenes, nodes)

    fmt1 = 'level1-scenes-{jobid}.txt'
    fmt2 = '{model}-ard-{jobid}.bash'

    batchid = uuid.uuid4().hex[0:10]
    batch_logdir = pjoin(logdir, 'batchid-{}'.format(batchid))
    batch_outdir = pjoin(outdir, 'batchid-{}'.format(batchid))

    # compute resources
    memory = 32 * nodes
    ncpus = 16 * nodes

    if dsh:
        files = []
        daemons = []
        outdirs = []
        daemon_fmt = 'luigid --background --logdir {}'

        # setup each block of scenes for processing
        for block in scattered:
            jobid = uuid.uuid4().hex[0:6]
            jobdir = pjoin(batch_logdir, 'jobid-{}'.format(jobid))
            job_outdir = pjoin(batch_outdir, 'jobid-{}'.format(jobid))

            if not exists(jobdir):
                os.makedirs(jobdir)

            if not exists(job_outdir):
                os.makedirs(job_outdir)

            # write the block of scenes to process
            out_fname = pjoin(jobdir, fmt1.format(jobid=jobid))
            with open(out_fname, 'w') as src:
                src.writelines(block)

            files.append(out_fname)
            daemons.append(daemon_fmt.format(jobdir))
            outdirs.append(job_outdir)

        # overwrite the contents of the first and last items
        # for an ugly styled list
        files[0] = 'FILES=("{}'.format(files[0])
        daemons[0] = 'DAEMONS=("{}'.format(daemons[0])
        outdirs[0] = 'OUTDIRS=("{}'.format(outdirs[0])
        files[-1] = '{}")'.format(files[-1])
        daemons[-1] = '{}")'.format(daemons[-1])
        outdirs[-1] = '{}")'.format(outdirs[-1])

        files = ['"{}"\n'.format(f) for f in files]
        daemons = ['"{}"\n'.format(f) for f in daemons]
        outdirs = ['"{}"\n'.format(f) for f in outdirs]

        files = ''.join(files)
        daemons = ''.join(daemons)
        outdirs = ''.join(outdirs)

        pbs = DSH_TEMPLATE.format(project=project, queue=queue, hours=hours,
                                  memory=memory, ncpus=ncpus,
                                  email=email, files=files, env=env,
                                  daemons=daemons, model=model,
                                  outdirs=outdirs, vertices=vertices,
                                  method=method)

        out_fname = pjoin(batch_logdir,
                          fmt2.format(model=model, jobid=batchid))
        with open(out_fname, 'w') as src:
            src.write(pbs)

        if test:
            print("Testing execution")
            print("qsub {}".format(out_fname))
        else:
            os.chdir(dirname(out_fname))
            subprocess.call(['qsub', out_fname])
    else:
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
                daemon = 'luigid --background --logdir {}'.format(jobdir)
                scheduler = ''


            pbs = PBS_TEMPLATE.format(project=project, queue=queue,
                                      hours=hours, memory=memory, ncpus=ncpus,
                                      email=email, env=env, daemon=daemon,
                                      model=model, scene_list=level1,
                                      outdir=job_outdir, scheduler=scheduler,
                                      vertices=vertices, method=method)

            out_fname = pjoin(jobdir, fmt1.format(jobid=jobid))
            with open(out_fname, 'w') as src:
                src.writelines(block)

            out_fname = pjoin(jobdir, fmt2.format(model=model, jobid=jobid))
            with open(out_fname, 'w') as src:
                src.write(pbs)

        if test:
            print("Testing execution")
            print("qsub {}".format(out_fname))
        else:
            os.chdir(dirname(out_fname))
            subprocess.call(['qsub', out_fname])


def _parser():
    """ Argument parser. """
    description = "qsub nbar jobs into n nodes."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--level1-list", help="The input level1 scene list.",
                        required=True)
    parser.add_argument("--vertices", default="(5, 5)", type=str,
                        help=("Number of vertices to evaluate the radiative "
                              "transfer at. JSON styled string is required."))
    parser.add_argument("--model", default="standard",
                        help="The type of ARD workflow to invoke.")
    parser.add_argument("--method", default="linear",
                        help="The interpolation method to invoke.")
    parser.add_argument("--outdir", help="The output directory.",
                        required=True)
    parser.add_argument("--logdir", required=True,
                        help="The base logging and scripts output directory.")
    parser.add_argument("--env", help="Environment script to source.",
                        required=True)
    parser.add_argument("--nodes", type=int,
                        help="The number of nodes to request.", default=10)
    parser.add_argument("--project", help="Project code to run under.",
                        required=True)
    parser.add_argument("--queue", help="Queue to submit the job into.",
                        default='normal')
    parser.add_argument("--hours", help="Job walltime in hours.", default=48)
    parser.add_argument("--email", help="Notification email address.",
                        default="your.name@something.com")
    parser.add_argument("--local-scheduler", help="Use a local scheduler.",
                        action='store_true')
    parser.add_argument("--dsh", help="Run using PBS Distributed Shell.",
                        action='store_true')
    parser.add_argument("--test", help="Test job execution.",
                        action='store_true')
    return parser


def main():
    """ Main execution. """
    parser = _parser()
    args = parser.parse_args()
    run(args.level1_list, args.vertices, args.model, args.method, args.outdir,
        args.logdir, args.env, args.nodes, args.project, args.queue,
        args.hours, args.email, args.local_scheduler, args.dsh, args.test)


if __name__ == '__main__':
    main()
