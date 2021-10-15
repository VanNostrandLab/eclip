#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
A pipeline for processing SE eCLIP fastq files and map them to both repeat elements and reference genome.
"""

import argparse
import os
import sys
import subprocess

import pysam
import cmder
from seqflow import Flow, task, logger


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


def basename(path, exts=None):
    path = os.path.basename(path)
    exts = exts if exts else ('.fq', 'fq.gz', '.fastq', '.fastq.gz')
    for ext in exts:
        if path.endswith(ext):
            path = path.rsplit(ext, maxsplit=1)[0]
            break
    return path
    

parser = argparse.ArgumentParser(description=__doc__, prog='eclip_se_fastq_to_bam')
parser.add_argument('--fastq', required=True, nargs='+', help='Path to FASTQ file(s).', type=file_path)
parser.add_argument('--adapters_fasta', nargs='+', required=True, type=file_path,
                    help="Path to the fasta file contains adapters and their sequences (for single-end dataset only.")
parser.add_argument('--name', nargs='+', help='Path to a UMI extracted FASTQ file.', type=file_path)
parser.add_argument('--outdir', help='Path to the output BAM file (must ends with .bam).')
parser.add_argument('--genome', required=True, help="Path to STAR reference genome index directory.", type=dir_path)
parser.add_argument('--genome_label', help="A short label for reference genome, e.g., hg19, mm10, default: genome.",
                    default='genome')
parser.add_argument('--repeat', required=True, help="Path to STAR repeat elements index directory.", type=dir_path)
parser.add_argument('--repeat_label', help="A short label for repeat elements, default: repbase.", default='repbase')
parser.add_argument('--cpus', type=int, help='Maximum number of CPU cores can be used for your job.', default=16)
parser.add_argument('--dry_run', action='store_true',
                    help='Print out steps and files involved in each step without actually running the pipeline.')
parser.add_argument('--verbose', action='store_true',
                    help='Print out processing details of each task.')
parser.add_argument('--job', help="Name of your job", default='eCLIP')
parser.add_argument('--email', help='Email address for notifying you the start, end, and abort of you job.')
parser.add_argument('--time', type=int, help='Time (in integer hours) for running your job.', default=36)
parser.add_argument('--memory', type=int, help='Amount of memory (in GB) for all CPU cores.', default=32)
parser.add_argument('--scheduler', choices=('pbs', 'qsub', 'slurm', 'sbatch'),
                    help='Name (case insensitive) of the scheduler on your cluster.')
parser.add_argument('--hold_submit', action='store_true',
                    help='Generate the submit script but hold it without submitting to the job scheduler.')


args = parser.parse_args()
outdir = args.outdir or os.getcwd()
dir_path(outdir)
os.chdir(outdir)

fastq, adapters_fasta, name = args.fastq, args.adapters_fasta, args.name
rtag, gtag = args.repeat_label, args.genome_label
if len(fastq) != len(adapters_fasta):
    logger.error('Number of items for fastq and adapters_fasta are not equal.')
    sys.exit(1)
if name:
    if not len(args.fastq) != len(args.name):
        logger.error('Number of items for fastq and name are not equal.')
    sys.exit(1)
else:
    name = [basename(n) for n in fastq]
name = [os.path.join(outdir, n) for n in name]
fastq_to_name = {fq: n for fq, n in zip(fastq, name)}
fastq_to_adapters = {n: adapter for n, adapter in zip(name, adapters_fasta)}


@task(inputs=fastq, outputs=lambda i: f'{fastq_to_name[i]}.umi.fastq.gz',
      cmd=['eclip_umi_extract', 'input', '-o', 'output'], cpus=args.cpus)
def extract_umi(fastq, output):
    pass


@task(inputs=extract_umi, outputs=lambda i: f'{i.rsplit(".umi.fastq.gz", maxsplit=1)[0]}.trim.fastq.gz')
def cut_adapt(fastq, output):
    key = fastq.rsplit(".umi.fastq.gz", maxsplit=1)[0]
    cmd = ['eclip_cut_adapt', fastq, '-o', output, '-a', fastq_to_adapters[key], '-c', args.cpus]
    cmder.run(cmd, stdout=sys.stdout, stderr=sys.stderr, log_cmd=False)


@task(inputs=cut_adapt, outputs=lambda i: f'{i.rsplit(".trim.fastq.gz", maxsplit=1)[0]}.{rtag}.mate1.gz')
def map_to_repbase(fastq, mate):
    bam = f'{mate.rsplit(".mate1.gz", maxsplit=1)[0]}.bam'
    cmd = ['star_repbase_map', fastq, '-x', args.repeat, '-c', args.cpus, '-o', bam]
    cmder.run(cmd, stdout=sys.stdout, stderr=sys.stderr, log_cmd=False)


@task(inputs=map_to_repbase, outputs=lambda i: f'{i.rsplit(f".{rtag}.mate1.gz", maxsplit=1)[0]}.{gtag}.bam')
def map_to_reference_genome(mate1, bam):
    cmd = ['star_genome_map', mate1, '-x', args.genome, '-c', args.cpus, '-o', bam]
    cmder.run(cmd, stdout=sys.stdout, stderr=sys.stderr, log_cmd=False)


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
    if args.scheduler.upper() in ('PBS', 'QSUB'):
        runtime, directive, exe, mail = f'{args.time}:00:00', pbs, 'qsub', pbs_email
        project = '/project/vannostrand'
    elif args.scheduler.upper() in ('SLURM', 'SBATCH'):
        days, hours = divmod(args.time, 24)
        runtime, directive, exe, mail = f'{days}-{hours:02}:00', sbatch, 'sbatch', sbatch_email
        project = '/storage/vannostrand'
    else:
        raise ValueError(f'Unsupported scheduler: {args.scheduler}, see help for supported schedulers.')
    root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    setattr(args, 'runtime', runtime)
    setattr(args, 'project', project)
    if args.debug:
        setattr(args, 'debug', '--debug ')
        setattr(args, 'program', os.path.join(root, 'rnaseq', parser.prog))
    else:
        setattr(args, 'debug', '')
        setattr(args, 'program', os.path.join(root, parser.prog))
    text = [directive, mail, code] if args.email else [directive, code]
    text = ''.join(text).format(**vars(args))

    submitter = os.path.join(args.outdir, 'submit.sh')
    with open(submitter, 'w') as o:
        o.write(text)

    print(f'Job submit script was saved to:\n    {submitter}')
    if args.hold_submit:
        print(f'Job {args.job} was not submitted yet, submit it after carefully review the submit script using:')
        print(f'    {exe} {submitter}')
    else:
        subprocess.run([exe, submitter])
        print(f'Job {args.job} was successfully submitted with the following settings:')
        data = {'Job name:': args.job, 'Output directory:': args.outdir,
                'Number of cores:': args.cores, 'Job memory:': args.memory,
                'Job runtime:': f'{runtime} (D-HH:MM)'}
        for k, v in data.items():
            print(f'{k:>20} {v}')


@logger.catch()
def main():
    if args.scheduler:
        schedule()
    else:
        flow = Flow('eCLIP', description=__doc__.strip())
        flow.run(dry_run=args.dry_run, cpus=args.cpus, verbose=args.verbose)


if __name__ == '__main__':
    main()
