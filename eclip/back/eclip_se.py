#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
A pipeline for processing eCLIP data to identify genomic locations of RNA binding proteins (RBPs).
"""

import argparse
import glob
import os
import gzip
import sys
import time
import subprocess
import json

import yaml
import datetime
import tempfile
import shutil
import itertools
import math

import pysam
import cmder
from seqflow import Flow, task, logger
from Bio.SeqIO.QualityIO import FastqGeneralIterator

from reporter import report


def file_path(p):
    if not os.path.isfile(p):
        logger.error(f'File "{p}" may not be a file or does not exist.')
        sys.exit(1)
    return p


def dir_path(p):
    if not os.path.isdir(p):
        logger.error(f'Path "{p}" may not be a directory or does not exist.')
        sys.exit(1)
    return p


parser = argparse.ArgumentParser(description=__doc__, prog='se_fastq_to_bam')
parser.add_argument('--fastq', required=True, help='Path to a UMI extracted FASTQ file.', type=file_path)
parser.add_argument('--bam', required=True,  help='Path to the output BAM file (must ends with .bam).')
parser.add_argument('--adapters_fasta', help="Path to the fasta file contains adapters and their sequences (for "
                                             "single-end dataset only.", required=True,  type=file_path)
parser.add_argument('--genome', help="Path to STAR reference genome index directory.", type=dir_path)
parser.add_argument('--repeat', help="Path to STAR repeat elements index directory.", type=dir_path)
parser.add_argument('--cpus', type=int, help='Maximum number of CPU cores can be used for your job.', default=16)
parser.add_argument('--dry_run', action='store_true',
                    help='Print out steps and files involved in each step without actually running the pipeline.')


args = parser.parse_args()
fastq, bam = args.fastq, args.bam
if not bam.endswith('.bam'):
    logger.error(f'Output BAM file "{bam}" does not end with .bam extension.')
    sys.exit(1)
name = args.name if args.name else bam.replace('.bam', '')
outdir = os.path.dirname(bam) or os.getcwd()
if not os.path.isdir(outdir):
    logger.error(f'Cane not set "{outdir}" as output directory.')
    sys.exit(1)
os.chdir(options.outdir)

ips, inputs, names = options.ip_fastqs, options.input_fastqs, options.labels
if len(ips) == len(inputs) == len(names):
    pass
else:
    logger.error('Number of items in ip_fastqs, input_fastqs, and names are not equal.')
    sys.exit(1)


@task(inputs=fastq, outputs=f'{name}.trim.fastq.gz')
def cut_adapt(fastq, output):
    cmd = ['eclip_cut_adapt', '-i', fastq, '-o', output,
           '-m', f'{name}.trim.metrics,{name}.trim.trim.metrics', '-a', args.adapters_fasta]
    cmder.run(cmd, msg=f'Cutting adapters for single-end read {fastq} ...')


@task(inputs=cut_adapt, outputs=f'{name}.map.to.repbase.mate1')
def map_to_repbase(fastq, mate):
    cmd = ['eclip_star_map', fastq, 'none', args.repeat, 100, f'{name}.map.to.repbase', args.cpus, 'repbase']
    cmder.run(cmd, stdout=sys.stdout, stderr=sys.stderr)


@task(inputs=map_to_repbase, outputs=f'{name}.map.to.genome.bam')
def map_to_reference_genome(mate1, bam):
    cmd = ['eclip_star_map', fastq, 'none', args.repeat, 1, f'{name}.map.to.genome', args.cpus, 'genome']
    cmder.run(cmd, stdout=sys.stdout, stderr=sys.stderr)


def schedule():
    sbatch = """#!/usr/bin/env bash

#SBATCH -n {cores}                        # Number of cores (-n)
#SBATCH -N 1                        # Ensure that all cores are on one Node (-N)
#SBATCH -t {runtime}                  # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH --mem={memory}G                   # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --job-name={job}          # Short name for the job
"""
    sbatch_email = """
#SBATCH --mail-user={email}
#SBATCH --mail-type=ALL
"""
    pbs = """ #!/usr/bin/env bash

#PBS -l nodes=1:ppn={cores}
#PBS -l walltime={runtime}
#PBS -l vmem={memory}gb
#PBS -j oe
#PBS -N {jobname}
"""
    pbs_email = """
#PBS -M {email}
#PBS -m abe
"""
    code = r"""
export TMPDIR={project}/tmp
export TEMP={project}/tmp
export TMP={project}/tmp

{program} \
    --verbose {debug}\
    --outdir {outdir} \
    --genome {genome} \
    --repeat {repeat} \
    --gtf {gtf} \
    --strand_direction {strand_direction} \
    --strand_specific {strand_specific} \
    --cores {cores} \
    {MANIFEST}
"""
    if options.scheduler.upper() in ('PBS', 'QSUB'):
        runtime, directive, exe, mail = f'{options.time}:00:00', pbs, 'qsub', pbs_email
        project = '/project/vannostrand'
    elif options.scheduler.upper() in ('SLURM', 'SBATCH'):
        days, hours = divmod(options.time, 24)
        runtime, directive, exe, mail = f'{days}-{hours:02}:00', sbatch, 'sbatch', sbatch_email
        project = '/storage/vannostrand'
    else:
        raise ValueError(f'Unsupported scheduler: {options.scheduler}, see help for supported schedulers.')
    root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    setattr(options, 'runtime', runtime)
    setattr(options, 'project', project)
    if options.debug:
        setattr(options, 'debug', '--debug ')
        setattr(options, 'program', os.path.join(root, 'rnaseq', parser.prog))
    else:
        setattr(options, 'debug', '')
        setattr(options, 'program', os.path.join(root, parser.prog))
    text = [directive, mail, code] if options.email else [directive, code]
    text = ''.join(text).format(**vars(options))

    submitter = os.path.join(options.outdir, 'submit.sh')
    with open(submitter, 'w') as o:
        o.write(text)

    print(f'Job submit script was saved to:\n    {submitter}')
    if options.hold_submit:
        print(f'Job {options.job} was not submitted yet, submit it after carefully review the submit script using:')
        print(f'    {exe} {submitter}')
    else:
        subprocess.run([exe, submitter])
        print(f'Job {options.job} was successfully submitted with the following settings:')
        data = {'Job name:': options.job, 'Output directory:': options.outdir,
                'Number of cores:': options.cores, 'Job memory:': options.memory,
                'Job runtime:': f'{runtime} (D-HH:MM)'}
        for k, v in data.items():
            print(f'{k:>20} {v}')


FOOTER = """
+======================================================================================================================+
|                                                                                                                      |
|                                                 MISSION ACCOMPLISHED                                                 |
|{hh_mm_ss}|
|                                                                                                                      |
+======================================================================================================================+
"""


@logger.catch()
def main():
    if options.scheduler:
        schedule()
    else:
        keys = ('MANIFEST', 'outdir', 'genome', 'repeat', 'processes')
        d = vars(options).copy()
        d['processes'] = options.cores
        setting = '\n'.join([f'{k.title():>22}: {v}' for k, v in d.items() if k in keys])
        logger.debug(HEADER)
        logger.trace(f'\nRunning eclip using the following settings:\n\n{setting}\n')
        flow = Flow('eCLIP', description=__doc__.strip())
        flow.run(dry=options.dry_run, processes=options.cores, verbose=options.verbose)
        logger.debug('')
        run_time = str(datetime.timedelta(seconds=int(time.perf_counter() - START_TIME)))
        logger.trace(FOOTER.format(hh_mm_ss=f'time consumed: {run_time}'.upper().center(118)))


if __name__ == '__main__':
    main()
