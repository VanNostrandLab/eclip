#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
A pipeline for processing eCLIP data to identify genomic locations of RNA binding proteins (RBPs).
"""

import collections
import glob
import os
import time
import subprocess
import json
import yaml
import datetime
import sys
import tempfile
import shutil
import math
from collections import defaultdict

import pysam as pysam
import scipy.stats as stats
import pandas as pd
import numpy as np
import ruffus
from loguru import logger

ignored_args = ['target_tasks', 'jobs', 'use_threads', 'just_print', 'log_file', 'verbose', 'forced_tasks']
parser = ruffus.cmdline.get_argparse(description=__doc__, prog='rnaseq',
                                     ignored_args=ignored_args)
parser.add_argument('MANIFEST', type=str,
                    help='Path to the manifest file that specifies paths for RNA-Seq data.')
parser.add_argument('--outdir', type=str,
                    help="Path to the output directory. Default to the current work "
                         "directory and if the specified path does not exist, it will "
                         "try to create it first.")
parser.add_argument('--genome', type=str,
                    help="Path to STAR reference genome index directory.")
parser.add_argument('--repeat', type=str,
                    help="Path to STAR repeat elements index directory.")
parser.add_argument('--track', type=str,
                    help="Name for the UCSC Genome Browser track, default: RNA-Seq",
                    default='RNA-Seq')
parser.add_argument('--track_label', type=str,
                    help="Label for the UCSC Genome Browser track, default: RNA-Seq",
                    default='RNA-Seq')
parser.add_argument('--track_genome', type=str,
                    help="Genome name for the UCSC Genome Browser track, default: RNA-Seq",
                    default='hg19')
parser.add_argument('--l2fc', type=int,
                    help="Only consider peaks at or above this log2 fold change cutoff.",
                    default=3)
parser.add_argument('--l10p', type=int,
                    help="Only consider peaks at or above this log10p value cutoff.",
                    default=3)
parser.add_argument('--enrichment_filter', type=int,
                    help="Pre-filter peaks that are enriched over input (default: 0).",
                    default=3)
parser.add_argument('--job', type=str,
                    help="Name of your job, default: RNA-Seq",
                    default='RNA-Seq')
parser.add_argument('--email', type=str,
                    help='Email address for notifying you the start, end, and abort of you job.')
parser.add_argument('--scheduler', type=str,
                    help='Name of the scheduler on your cluster, '
                         'e.g., PBS (or QSUB) or SBATCH (or SLURM), case insensitive.')
parser.add_argument('--time', type=int,
                    help='Time (in integer hours) for running your job, default: 24.',
                    default=24)
parser.add_argument('--memory', type=int,
                    help='Amount of memory (in GB) for all cores needed for your job, default: 32.',
                    default=32)
parser.add_argument('--cores', type=int,
                    help='Number of CPU cores can be used for your job, default: 8.',
                    default=8)
parser.add_argument('--verbose', action='store_true',
                    help='Print out detailed processing messages.')
parser.add_argument('--verbose_path', type=int, default=2,
                    help='Whether input and output paths are abbreviated.')
parser.add_argument('--debug', action='store_true',
                    help='Invoke debug mode.')
parser.add_argument('--dry_run', action='store_true',
                    help='Print out steps and files involved in each step without actually '
                         'running the pipeline.')
parser.add_argument('--hold_submit', action='store_true',
                    help='Generate the submit script but hold it without submitting to the job scheduler. '
                         'Useful when you want to further review the submit script to make sure everything '
                         'looks good and ready to submit.')
parser.add_argument('--target_tasks', metavar='TASKS', action="append", type=str, default=[],
                    help='Target task(s) of pipeline..')
parser.add_argument('--forced_tasks', metavar='TASKS', action="append", type=str, default=[],
                    help='Task(s) which will be included even if they are up to date.')

START_TIME = time.perf_counter()


def parse_and_sanitize_options():
    options = parser.parse_args()
    try:
        with open(options.MANIFEST) as stream:
            if options.MANIFEST.endswith('.json'):
                manifest = json.load(stream)
            elif options.MANIFEST.endswith('.yaml'):
                manifest = yaml.safe_load(stream)
            else:
                raise ValueError(f'Unknown manifest file format, manifest filename needs to end with .yaml or .json.')
    except IOError:
        raise ValueError(f'Manifest {options.MANIFEST} may not exist or not be file.')
    setattr(options, 'manifest', manifest)
    setattr(options, 'outdir', options.outdir or manifest.get('outdir', '') or os.getcwd())
    if not os.path.isdir(options.outdir):
        try:
            os.mkdir(options.outdir)
        except OSError:
            raise OSError(f'Outdir {options.outdir} does not exist and cannot be created!')
    setattr(options, 'genome', options.genome or manifest.get('genome', ''))
    if not os.path.isdir(options.genome):
        raise ValueError(f'Reference genome index {options.genome} is not a directory or does not exist.')
    
    setattr(options, 'repeat', options.repeat or manifest.get('repeat', ''))
    if not os.path.isdir(options.repeat):
        raise ValueError(f'Repeat element index {options.repeat} is not a directory or does not exist.')
    
    barcodes_fasta = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'barcodes.fasta')
    setattr(options, 'barcodes_fasta', manifest.get('barcodes_fasta', '') or barcodes_fasta)
    if not os.path.isfile(options.barcodes_fasta):
        raise ValueError(f'Barcodes fasta {options.barcodes_fasta} is not a file or does not exist.')
    setattr(options, 'randomer_length', manifest.get('randomer_length', '5'))
    setattr(options, 'blacklist_bed', manifest.get('blacklist_bed', ''))
    
    setattr(options, 'dataset', manifest.get('dataset', 'eclip'))
    setattr(options, 'species', manifest.get('species', 'hg19'))
    setattr(options, 'track', options.dataset)
    setattr(options, 'track_label', options.dataset)
    setattr(options, 'track_genome', options.species)
    
    setattr(options, 'verbose', 3 if options.dry_run else 0)
    setattr(options, 'jobs', options.cores)
    setattr(options, 'multiprocess', options.cores)
    setattr(options, 'use_threads', False)
    setattr(options, 'just_print', options.dry_run)
    setattr(options, 'verbose_abbreviated_path', options.verbose_path)
    setattr(options, 'exceptions_terminate_immediately', True)
    setattr(options, 'one_second_per_job', True)
    setattr(options, 'log_exceptions', True)
    return options


options = parse_and_sanitize_options()
os.chdir(options.outdir)
ECLIP, HUB, QC = 'eclip', 'hub', 'qc'


class Read:
    def __init__(self, record, read_type):
        self.type = read_type.upper()
        self.name = record.get('name', '')
        if not self.name:
            raise KeyError(f'No name was assigned for {self.type}\nread {record}.')
        self.fastq1 = record.get('fastq', '') or record.get('fastq1', '')
        if self.fastq1:
            if not os.path.isfile(self.fastq1):
                raise ValueError(f'FASTQ1 {self.fastq1} may not be a file or does not exist.')
        else:
            raise KeyError(f'No fastq or fastq1 was assigned for {self.type} read\n{record}.')
        self.fastq2 = record.get('fastq2', '')
        if self.fastq2 and not os.path.isfile(self.fastq2):
            raise ValueError(f'FASTQ2 {self.fastq2} may not be a file or does not exist.')
        self.barcodes = record.get('barcodes', ['', ''])
        self.adapters = record.get('adapters', '')
        
        self.key = f'{options.dataset}.{self.name}'
        self.paired = True if self.fastq2 else False
        self.link1 = f'{ECLIP}/{self.key}.r1.fastq.gz'
        self.link2 = f'{ECLIP}/{self.key}.r2.fastq.gz' if self.fastq2 else ''
        self.bam = f'{ECLIP}/{self.key}.bam'
        self.bed = f'{ECLIP}/{self.key}.peak.clusters.bed' if read_type == 'ip' else ''
        self.full_bed = f'{ECLIP}/{self.key}.peak.clusters.full.bed' if read_type == 'ip' else ''
        self.normalized_bed = f'{ECLIP}/{self.key}.peak.clusters.normalized.bed' if read_type == 'ip' else ''
        self.compressed_bed = f'{ECLIP}/{self.key}.peak.clusters.compressed.bed' if read_type == 'ip' else ''


class Sample:
    def __init__(self, ip_read, input_read):
        self.ip_read = Read(ip_read, 'ip')
        self.input_read = Read(input_read, 'input')


def size(file, formatting=True):
    if os.path.isfile(file):
        s = os.path.getsize(file)
        if formatting:
            if s < 1000:
                s = f'[{s:.2f}B]'
            elif s < 1000 * 1000:
                s = f'[{s / 1000:.2f}KB]'
            elif s < 1000 * 1000 * 1000:
                s = f'[{s / 1000 / 1000:.2f}MB]'
            else:
                s = f'[{s / 1000 / 1000 / 1000:.2f}GB]'
    else:
        s = '[?KB]' if formatting else 0
    return s


def parse_and_validate_samples():
    def estimate_max_processes():
        max_size = max(sizes) / (1000 * 1000 * 1000) * 5
        n = int(options.memory / max_size)
        if n == 0:
            n = 1
        elif n > options.cores:
            n = options.cores
        return n
    
    samples, reads, sizes, types = [], {}, [], []
    for sample in options.manifest.get('samples', []):
        ip_read, input_read = sample.get('ip_read', {}), sample.get('input_read', {})
        if not ip_read:
            raise KeyError(f'No key ip_read was found in sample\n{sample}.')
        if not input_read:
            raise KeyError('No key input_read was found in sample\n{sample}.')
        ip_read, input_read = Read(ip_read, 'ip'), Read(input_read, 'input')
        if ip_read.key in samples:
            raise ValueError(f'Duplicated names found for ip_read\n{ip_read}.')
        if input_read.key in samples:
            raise ValueError(f'Duplicated names found for input_read\n{input_read}.')
        samples.append([ip_read, input_read])
        reads[ip_read.key], reads[input_read.key] = ip_read, input_read
        sizes.extend([size(ip_read.fastq1, formatting=False), size(ip_read.fastq2, formatting=False),
                      size(input_read.fastq1, formatting=False), size(input_read.fastq2, formatting=False)])
        types.extend([ip_read.paired, input_read.paired])
    if not samples:
        raise KeyError('Manifest file misses key samples or does not contain any samples.')
    assert len(set(types)) == 1, ValueError('Reads mixed with single-end and paired-end fastq file specifications.')
    return samples, reads, estimate_max_processes(), 'paired' if types[0] else 'single'


HASH = 100000
SAMPLES, READS, PROCESSES, TYPE = parse_and_validate_samples()
VERSION = 1.0
HEADER = fr"""
                           ____   _       ___   ____
                    ___   / ___| | |     |_ _| |  _ \
                   / _ \ | |     | |      | |  | |_) |
                  |  __/ | |___  | |___   | |  |  __/
                   \___|  \____| |_____| |___| |_|
                               
                               VERSION {VERSION}
"""

logger.remove()
out = sys.stderr if options.debug else sys.stdout
logger.add(out, format="<light-green>[{time:YYYY-MM-DD HH:mm:ss}]</light-green> <level>{message}</level>",
           filter=lambda record: record["level"].name == "TRACE",
           level="TRACE")
logger.add(out, format="<level>{message}</level>", filter=lambda record: record["level"].name == "DEBUG")
logger.add(out, format="<light-green>[{time:HH:mm:ss}]</light-green> <level>{message}</level>", level="INFO")


@logger.catch(onerror=lambda _: sys.exit(1))
def cmding(cmd, **kwargs):
    """ Run cmd or raise exception if run fails. """
    def format_cmd(command):
        if isinstance(command, str):
            exe = command.split()[0]
        else:
            command = [str(c) for c in command]
            exe = command[0]
            command = ' '.join([f'\\\n  {c}' if c.startswith('-') or '<' in c or '>' in c else c for c in command])
            command = command.splitlines()
            commands = []
            for i, c in enumerate(command):
                if i == 0:
                    commands.append(c)
                else:
                    if len(c) <= 80:
                        commands.append(c)
                    else:
                        items = c.strip().replace(' \\', '').split()
                        commands.append(f'  {items[0]} {items[1]} \\')
                        for item in items[2:]:
                            commands.append(' ' * (len(items[0]) + 3) + item + ' \\')
            command = '\n'.join(commands)
            if command.endswith(' \\'):
                command = command[:-2]
        return exe, command
    
    def parse_profile():
        try:
            with open(profile_output) as f:
                t, m = f.read().strip().split()
                t = t.split(".")[0]
                try:
                    hh, mm, ss = t.split(':')
                except ValueError:
                    hh, (mm, ss) = 0, t.split(':')
                t = f'{int(hh):02d}:{int(mm):02d}:{int(ss):02d}'
                m = float(m)
                if m < 1000:
                    m = f'{m:.2f}KB'
                elif m < 1000 * 1000:
                    m = f'{m / 1000:.2f}MB'
                else:
                    m = f'{m / 1000 / 1000:.2f}GB'
                s = f'{t} {m}'
        except FileNotFoundError:
            s = '00:00:00 0.00KB'
        return s
    
    message, start_time = kwargs.pop('message', ''), time.perf_counter()
    program, cmd = format_cmd(cmd)
    if message:
        logger.info(message)
    logger.debug(cmd)
    cwd = kwargs.pop('cwd', options.outdir)
    profile_output = tempfile.mktemp(suffix='.txt', prefix='.profile.', dir=cwd)
    try:
        cmd = f'time -f "%E %M" -o {profile_output} {cmd}'
        process = subprocess.Popen(cmd, universal_newlines=True, shell=True, cwd=cwd,
                                   stdout=subprocess.PIPE, stderr=subprocess.PIPE, **kwargs)
        stdout, stderr = process.communicate()
        if process.returncode:
            logger.error(f'Failed to run {program} (exit code {process.returncode}):\n{stderr or stdout}')
            sys.exit(process.returncode)
        if message:
            logger.info(message.replace(' ...', f' complete [{parse_profile()}].'))
    finally:
        os.unlink(profile_output)


@ruffus.jobs_limit(1)
@ruffus.follows(ruffus.mkdir(ECLIP))
@ruffus.originate([read.link1 for read in READS.values()])
def soft_link(link1):
    """ Create soft links for original fastq files. """
    
    def make_link(path, link):
        if path:
            if path == os.path.abspath(link):
                message = "No symbolic link was made for {path}! You are directly working on the original file!"
                logger.warning(message)
            else:
                if not os.path.exists(link):
                    message = f'Soft link fastq: {os.path.basename(path)} ...'
                    cmding(f'ln -s {path} {link}', message=message)
                    
    read = READS[os.path.basename(link1.replace('.r1.fastq.gz', ''))]
    link1, link2 = read.link1, read.link2
    make_link(read.fastq1, link1)
    if link2:
        make_link(read.fastq2, link2)


def cleanup_demux():
    fastq, rm = set(), [options.dataset]
    for key, read in READS.items():
        if read.link2:
            fastq.add(f'{ECLIP}/{read.key}.{read.barcodes[0]}.r1.fq.gz')
            fastq.add(f'{ECLIP}/{read.key}.{read.barcodes[0]}.r2.fq.gz')
            fastq.add(f'{ECLIP}/{read.key}.{read.barcodes[1]}.r1.fq.gz')
            fastq.add(f'{ECLIP}/{read.key}.{read.barcodes[1]}.r2.fq.gz')
            rm.append(read.name)
            rm.extend(read.barcodes)
        else:
            fastq.add(f'{ECLIP}/{read.key}.umi.r1.fq.gz')
    rm = [os.path.join(ECLIP, r) for r in rm] + [x for x in glob.iglob(f'{ECLIP}/*.fq.gz') if x not in fastq]
    _ = [os.unlink(x) for x in rm if os.path.isfile(x)]


@ruffus.jobs_limit(PROCESSES)
@ruffus.posttask(cleanup_demux)
@ruffus.transform(soft_link, ruffus.suffix('.fastq.gz'), '.prepare.done')
def prepare_reads(link, output):
    """Extract UMIs for single-end reads or demultiplex paired-end reads."""
    read = READS[os.path.basename(link.replace('.r1.fastq.gz', ''))]
    fastq1, fastq2 = read.link1, read.link2
    if fastq2:
        barcodes, cwd = read.barcodes, os.path.join(options.outdir, 'eclip')
        message = f'Demultiplexing paired-end reads {fastq1} {size(fastq1)} and\n{" " * 43}{fastq2} {size(fastq2)} ...'
        cmd = ['demux',
               '--fastq_1', os.path.basename(fastq1),
               '--fastq_2', os.path.basename(fastq2),
               '--newname', read.name,
               '--expectedbarcodeida', barcodes[0],
               '--expectedbarcodeidb', barcodes[1],
               '--dataset', options.dataset,
               '--barcodesfile', options.barcodes_fasta,
               '--length', options.randomer_length,
               '--metrics', f'{options.dataset}.{read.name}.demux.metrics']
        cmding(cmd, message=message, cwd=cwd)
    else:
        message = f'Extract UMIs for single-end read {fastq1} {size(fastq1)} ...'
        cmd = ['umi_tools', 'extract',
               '--random-seed', 1,
               '--stdin', fastq1,
               '--bc-pattern', 'NNNNNNNNNN',
               '--log', fastq1.replace('.fastq.gz', '.extract.metrics'),
               '--stdout', fastq1.replace('.r1.fastq.gz', '.umi.r1.fq.gz')]
        cmding(cmd, message=message)
    
    with open(output, 'w') as o:
        o.write('')


def prepare_adapter_cut():
    fastq = set()
    for key, read in READS.items():
        name, barcodes = read.name, read.barcodes
        if TYPE == 'paired':
            fastq.add(f'{ECLIP}/{key}.{barcodes[0]}.r1.clean.fastq')
            fastq.add(f'{ECLIP}/{key}.{barcodes[1]}.r1.clean.fastq')
        else:
            fastq.add(f'{ECLIP}/{key}.umi.r1.clean.fastq')
    return fastq


@ruffus.jobs_limit(1)
@ruffus.follows(prepare_reads)
@ruffus.originate(prepare_adapter_cut())
def cut_adapt(fastq):
    def parse_adapters(flag, fasta):
        adapters = []
        with open(fasta) as f:
            for line in f:
                if line.strip() and not line.startswith('>'):
                    adapters.extend([flag, line.strip()])
        return adapters
    
    def parse_overlap(txt):
        with open(txt) as f:
            return f.read().strip()
    
    def get_adapters(read):
        adapter, adapters = read.adapters, []
        if adapter:
            adapters1 = parse_adapters('-a', adapter)
            adapters2 = adapters1
            overlap1, overlap2 = '1', '5'
        else:
            cmd = ['parse_barcodes.sh', options.randomer_length, options.barcodes_fasta] + read.barcodes
            folder = tempfile.mkdtemp(dir=os.path.join(options.outdir, 'eclip'))
            cmding(cmd, message='Parsing barcodes and finding adapters ...', cwd=folder)
            adapters = parse_adapters('-a', os.path.join(folder, 'a_adapters.fasta'))
            adapters += parse_adapters('-A', os.path.join(folder, 'A_adapters.fasta'))
            adapters += parse_adapters('-g', os.path.join(folder, 'g_adapters.fasta'))
            adapters1 = adapters
            adapters2 = parse_adapters('-A', os.path.join(folder, 'A_adapters.fasta'))
            overlap1 = parse_overlap(os.path.join(folder, 'trim_first_overlap_length.txt'))
            overlap2 = parse_overlap(os.path.join(folder, 'trim_again_overlap_length.txt'))
            shutil.rmtree(folder)
        return adapters1, adapters2, overlap1, overlap2
    
    def get_ios(out1):
        if '.umi.' in out1:
            input11, input12 = out1.replace('.clean.', '.').replace('.fastq', '.fq.gz'), ''
            ios1 = ['-o', out1.replace('.clean.', '.trim.').replace('.fastq', '.fastq.gz'),
                    input11, '>', out1.replace('.r1.clean.fastq', '.trim.metrics')]
            msg1 = f"Cutting adapters for single read {out1.replace('.clean.', '.')} {{s1}} (first round) ..."
            input21, input22 = out1.replace('.clean.', '.trim.').replace('.fastq', '.fastq.gz'), ''
            ios2 = ['-o', out1, input21, '>', out1.replace('.r1.clean.fastq', '.clean.metrics')]
            msg2 = f"Cutting adapters for single read {out1.replace('.clean.', '.trim.')} {{s1}} (second round) ..."
        else:
            input11 = out1.replace('.clean.', '.').replace('.fastq', '.fq.gz')
            input12 = out1.replace('.clean.', '.').replace('.fastq', '.fq.gz').replace('.r1.', '.r2.')
            ios1 = ['-o', out1.replace('.clean.', '.trim.').replace('.fastq', '.fastq.gz'),
                    '-p', out1.replace('.clean.', '.trim.').replace('.fastq', '.fastq.gz').replace('.r1.', '.r2.'),
                    input11, input12, '>', out1.replace('.r1.clean.fastq', '.trim.metrics')]
            msg1 = (f"Cutting adapters for paired reads {out1.replace('.clean.', '.')} {{s1}} and\n{' ' * 45}"
                    f"{out1.replace('.clean.', '.').replace('.r1.', '.r2.')} {{s2}} (first round) ...")
            input21 = out1.replace('.clean.', '.trim.').replace('.fastq', '.fastq.gz')
            input22 = out1.replace('.clean.', '.trim.').replace('.fastq', '.fastq.gz').replace('.r1.', '.r2.')
            ios2 = ['-o', out1, '-p', out1.replace('.r1.', '.r2.'),
                    input21, input22, '>', out1.replace('.r1.clean.fastq', '.clean.metrics')]
            msg2 = (f"Cutting adapters for paired reads {out1.replace('.clean.', '.trim.')} {{s1}} and\n{' ' * 45}"
                    f"{out1.replace('.clean.', '.trim.').replace('.r1.', '.r2.')} {{s2}} (second round) ...")
        return ios1, ios2, msg1, msg2, input11, input12, input21, input22
    
    def trim_adapters(adapters, overlap, ios, message):
        cmd = ['cutadapt', '-O', overlap, '--times', '2', '-e', '0.0', '-j', options.cores, '-m', '18',
               '--quality-cutoff', '6', '--match-read-wildcards'] + adapters + ios
        cmding(cmd, message=message)
    
    link = '.'.join(fastq.split('.')[:-4])
    read = READS[os.path.basename(link)]
    adapters1, adapters2, overlap1, overlap2 = get_adapters(read)
    ios1, ios2, msg1, msg2, input11, input12, input21, input22 = get_ios(fastq)
    trim_adapters(adapters1, overlap1, ios1, msg1.format(s1=size(input11), s2=size(input12)))
    trim_adapters(adapters2, overlap2, ios2, msg2.format(s1=size(input21), s2=size(input22)))
    cmding(f'rm {ECLIP}/*.trim.fastq.gz')


@ruffus.jobs_limit(PROCESSES)
@ruffus.transform(cut_adapt, ruffus.suffix('.clean.fastq'), '.clean.sort.fastq')
def sort_fastq(fastq, output):
    cmd = ['fastq-sort', '--id', fastq, '>', output]
    cmding(cmd, message=f'Sort fastq {fastq} {size(fastq)} ...')
    
    fastq2 = fastq.replace('.r1.clean.fastq', '.r2.clean.fastq')
    if os.path.isfile(fastq2):
        output2 = output.replace('.r1.clean.sort.fastq', '.r2.clean.sort.fastq')
        cmd = ['fastq-sort', '--id', fastq2, '>', output2]
        cmding(cmd, message=f'Sort fastq {fastq2} {size(fastq2)} ...')


@ruffus.jobs_limit(1)
@ruffus.follows(ruffus.mkdir(f'{ECLIP}/repeat.elements.map'))
@ruffus.transform(sort_fastq,
                  ruffus.formatter(r'.+/(?P<BASENAME>.*).r1.clean.sort.fastq'),
                  ECLIP + '/repeat.elements.map/{BASENAME[0]}/Unmapped.out.mate1')
def map_to_repeat_elements(fastq, mate1):
    fastq1, fastq2 = fastq, fastq.replace('.r1.clean.sort.fastq', '.r2.clean.sort.fastq')
    prefix = os.path.dirname(mate1)
    if not os.path.isdir(prefix):
        os.mkdir(prefix)
    cmd = ['STAR',
           '--runMode', 'alignReads',
           '--runThreadN', options.cores,
           '--alignEndsType', 'EndToEnd',
           '--genomeDir', options.repeat,
           '--genomeLoad', 'NoSharedMemory',
           '--outBAMcompression', 10,
           '--outFileNamePrefix', f"{prefix}/",
           '--outFilterMultimapNmax', 30,
           '--outFilterMultimapScoreRange', 1,
           '--outFilterScoreMin', 10,
           '--outFilterType', 'BySJout',
           '--outReadsUnmapped', 'Fastx',
           '--outSAMattrRGline', 'ID:foo',
           '--outSAMattributes', 'All',
           '--outSAMmode', 'Full',
           '--outSAMtype', 'BAM', 'Unsorted',
           '--outSAMunmapped', 'Within',
           '--outStd', 'Log',
           '--readFilesIn', fastq1]
    if os.path.exists(fastq2):
        cmd.append(fastq2)
        message = (f'Map paired reads {fastq1} {size(fastq1)} and\n{28 * " "}'
                   f'{fastq2} {size(fastq2)} to repeat elements ...')
    else:
        message = f'Map single read {fastq1} {size(fastq1)} to repeat elements ...'
    cmding(cmd, message=message)


@ruffus.jobs_limit(PROCESSES)
@ruffus.transform(map_to_repeat_elements, ruffus.suffix('.out.mate1'), '.out.mate1.sort.fastq')
def sort_mate(mate1, output):
    cmd = ['fastq-sort', '--id', mate1, '>', output]
    cmding(cmd, message=f'Sort mate1 {mate1} {size(mate1)} ...')
    
    mate2 = mate1.replace('.out.mate1', '.out.mate2')
    if os.path.isfile(mate2):
        output2 = output.replace('.out.mate1.sort.fastq', '.out.mate2.sort.fastq')
        cmd = ['fastq-sort', '--id', mate2, '>', output2]
        cmding(cmd, message=f'Sort mate2 {mate2} {size(mate2)} ...')


@ruffus.jobs_limit(1)
@ruffus.follows(ruffus.mkdir(ECLIP + '/reference.genome.map'))
@ruffus.transform(sort_mate,
                  ruffus.formatter(r'.+/Unmapped.out.mate1.sort.fastq'),
                  ECLIP + '/reference.genome.map/{subdir[0][0]}/Aligned.out.bam')
def map_to_reference_genome(mate1, bam):
    # '--outSAMunmapped' flag needs to be set to 'Within', otherwise barcode_collapse.py for duplication removal will 
    # throw out name not match error.
    prefix = os.path.dirname(bam)
    if not os.path.isdir(prefix):
        os.mkdir(prefix)
    mate2 = mate1.replace('Unmapped.out.mate1.sort.fastq', 'Unmapped.out.mate2.sort.fastq')
    cmd = ['STAR',
           '--runMode', 'alignReads',
           '--runThreadN', options.cores,
           '--alignEndsType', 'EndToEnd',
           '--genomeDir', options.genome,
           '--genomeLoad', 'NoSharedMemory',
           '--outBAMcompression', 10,
           '--outFileNamePrefix', f"{prefix}/",
           '--outFilterMultimapNmax', 1,
           '--outFilterMultimapScoreRange', 1,
           '--outFilterScoreMin', 10,
           '--outFilterType', 'BySJout',
           '--outReadsUnmapped', 'Fastx',
           '--outSAMattrRGline', 'ID:foo',
           '--outSAMattributes', 'All',
           '--outSAMmode', 'Full',
           '--outSAMtype', 'BAM', 'Unsorted',
           '--outSAMunmapped', 'Within',
           '--outStd', 'Log',
           '--readFilesIn', mate1]
    if os.path.exists(mate2):
        cmd.append(mate2)
        message = f'Map paired mates {mate1} {size(mate1)} and\n{28 * " "}{mate2} {size(mate2)} to reference genome ...'
    else:
        message = f'Map single mate {mate1} {size(mate1)} to reference genome ...'
    cmding(cmd, message=message)


def name_sort_bam(bam, out):
    cmding(f'samtools sort -n -@ {options.cores} -m {min(4, int(options.memory / PROCESSES))}G -o {out} {bam}',
           message=f'Name sorting {bam} {size(bam)} ...')


def position_sort_bam(bam, out):
    cmding(f'samtools sort -@ {options.cores} -m {min(4, int(options.memory / PROCESSES))}G -o {out} {bam}',
           message=f'Sorting {bam} {size(bam)} ...')
    
    
def index_sorted_bam(bam):
    cmding(f'samtools index -@ {options.cores} {bam}', message=f'Indexing {bam} {size(bam)} ...')


def merge_paired_bam(bam, out):
    if not os.path.exists(out):
        key = out.replace(".merge.bam", "")
        barcodes = READS[os.path.basename(key)].barcodes
        if barcodes[0] == 'NIL':
            cmding(f'mv {bam} {out}')
        else:
            b1, b2 = barcodes
            if b1 in bam:
                b1, b2 = bam, bam.replace(b1, b2)
            else:
                b1, b2 = bam.replace(b2, b1), bam
            cmding(f'samtools merge -@ {options.cores} {out} {b1} {b2}',
                   message=f'Merging {b1} {size(b1)} and {b2} {size(b2)} ...')


@ruffus.jobs_limit(1)
@ruffus.transform(map_to_reference_genome, ruffus.formatter(r'.+/Aligned.out.bam'), ECLIP + '/{subdir[0][0]}.sort.bam')
def prepare_bam(bam, out):
    if TYPE == 'single':
        name_sort = out.replace('.sort.bam', '.name.sort.bam')
        name_sort_bam(bam, name_sort)
        position_sort_bam(name_sort, out)
        index_sorted_bam(out)
        cmding(f'rm {name_sort}')
    else:
        name_sort_bam(bam, out)


@ruffus.jobs_limit(PROCESSES)
@ruffus.transform(prepare_bam, ruffus.suffix('.sort.bam'), '.sort.dedup.bam')
def dedup_bam(bam, out):
    """Collapse barcodes of paired-end bam or umi_tools dedup single-end bam."""
    if TYPE == 'single':
        cmd = ['umi_tools', 'dedup', '--random-seed', 1, '--stdin', bam, '--method', 'unique',
               '--output-stats', out.replace(".dedup.bam", ".dedup.metrics"), '--stdout', out]
        message = f'Deduplicating {bam} {size(bam)} by umi_tools dedup ...'
    else:
        cmd = f'barcode_collapse.py -o {out} -m {out.replace(".bam", ".collapse.metrics")} -b {bam}'
        message = f'Deduplicating {bam} {size(bam)} by collapsing barcodes ...'
    cmding(cmd, message=message)


@ruffus.jobs_limit(1)
@ruffus.transform(dedup_bam, ruffus.suffix('.sort.dedup.bam'), '.sort.dedup.sort.bam')
def sort_bam(bam, out):
    position_sort_bam(bam, out)


@ruffus.jobs_limit(1)
@ruffus.transform(sort_bam, ruffus.regex(r"^(.+)\.(.+)\.sort.dedup.sort.bam$"), r'\1.merge.bam')
def merge_bam(bam, out):
    if TYPE == 'single':
        cmding(f'cp {bam} {out}')
    else:
        merge_paired_bam(bam, out)


@ruffus.jobs_limit(1)
@ruffus.transform(merge_bam, ruffus.suffix('.merge.bam'), '.bam')
def index_bam(bam, out):
    if TYPE == 'paired':
        cmding(f'samtools view -f 128 -@ {options.cores} -b -o {out} {bam}',
               message=f'Extracting r2 reads from {bam} {size(bam)} ...')
    else:
        cmding(f'cp {bam} {out}')
    if not os.path.exists(f'{bam}.bai'):
        index_sorted_bam(out)
        

@ruffus.jobs_limit(PROCESSES)
@ruffus.follows(ruffus.mkdir('hub'))
@ruffus.follows(index_bam)
@ruffus.transform([read.bam for read in READS.values()],
                  ruffus.formatter(r'.+/(?P<BASENAME>.*).bam'), HUB + '/{BASENAME[0]}.plus.bw')
def make_bigwig_files(bam, bigwig):
    def bam_to_bigwig(bam, scale, strand, bw):
        bg, bg_sort = f'{bw[:-3]}.sort.bg', f'{bw[:-3]}.bg'
        cmd = f'genomeCoverageBed -ibam {bam} -bg -scale {scale} -strand {strand} -du -split > {bg}'
        cmding(cmd, message=f'Calculating genome coverage for {bam} {size(bam)} ...')
        cmd = f'bedSort {bg} {bg_sort}'
        cmding(cmd, message=f'Sorting {bg} {size(bg)} ...')
        cmd = f'bedGraphToBigWig {bg_sort} {options.genome}/chrNameLength.txt {bw}'
        cmding(cmd, message=f'Converting {bg_sort} {size(bg_sort)} to {bw} ...')
        cmding(f'rm {bg}')
    
    if not os.path.exists(bigwig):
        message, start_time = f'Make BigWig files for {bam} ...', time.perf_counter()
        logger.info(message)
        pos_bw, neg_bw = bigwig, bigwig.replace('.plus.bw', '.minus.bw')
        with pysam.AlignmentFile(bam, 'rb') as sam:
            total_reads = sam.mapped
        total_reads = total_reads if TYPE == 'single' else total_reads / 2
        if total_reads:
            scale = 1000000.0 / total_reads
        else:
            logger.error(f'No reads was found in BAM {bam}')
            ruffus.touch_file(bigwig)
            return
        if TYPE == 'single':
            bam_to_bigwig(bam, scale, '+', pos_bw)
            bam_to_bigwig(bam, -1 * scale, '-', neg_bw)
        else:
            bam_to_bigwig(bam, -1 * scale, '-', pos_bw)
            bam_to_bigwig(bam, scale, '+', neg_bw)
        run_time = int(time.perf_counter() - start_time)
        message = message.replace(' ...', f' completed in [{str(datetime.timedelta(seconds=run_time))}].')
        logger.info(message)


@ruffus.merge(make_bigwig_files, HUB + '/hub.txt')
def make_hub_files(inputs, output):
    message, start_time = 'Make hub track file ...', time.perf_counter()
    logger.info(message)
    header = f"""hub {options.track.replace(' ', '_')}
shortLabel {options.track_label}
longLabel {options.track_label}
useOneFile on
email {options.email if options.email else 'fei.yuan@bcm.edu'}

genome {options.track_genome}

track {options.track.replace(' ', '_')}
shortLabel {options.track_label}
longLabel {options.track_label}
type bigWig
superTrack on
"""
    block = """
track {basename}
shortLabel {basename}
longLabel {basename}
type bigWig
visibility full
alwaysZero on
autoScale on
aggregate transparentOverlay
showSubtrackColorOnUi on
parent {track}
container multiWig

    track {name1}
    bigDataUrl {plus}
    shortLabel {basename} Plus strand
    longLabel {basename} Plus strand
    type bigWig
    color 0,100,0
    parent {basename}

    track {name2}
    bigDataUrl {minus}
    shortLabel {basename} Minus strand
    longLabel {basename} Minus strand
    type bigWig
    color 0,100,0
    parent {basename}
    """
    
    track = options.track.replace(' ', '_')
    with open(output, 'w') as o:
        o.write(header)
        for bw in inputs:
            plus = os.path.basename(bw)
            name1 = plus.replace('.bw', '').replace('.', '_')
            name2 = name1.replace('plus', 'minus')
            basename = plus.replace('.plus.bw', '')
            minus = plus.replace('.plus.bw', '.minus.bw')
            o.write(block.format(track=track, name1=name1, name2=name2, basename=basename, plus=plus, minus=minus))
    run_time = int(time.perf_counter() - start_time)
    message = message.replace(' ...', f' completed in [{str(datetime.timedelta(seconds=run_time))}].')
    logger.info(message)


@ruffus.jobs_limit(1)
@ruffus.follows(make_hub_files)
@ruffus.transform([read.bam for read in READS.values() if read.type == 'IP'],
                  ruffus.suffix('.bam'), '.peak.clusters.bed')
def clipper(bam, bed):
    cmd = f'clipper --species {options.species} --processors {options.cores} --bam {bam} --outfile {bed}'
    cmding(cmd, message=f'Calling peaks from {bam} {size(bam)} using clipper ...')


@ruffus.jobs_limit(1)
@ruffus.follows(clipper)
@ruffus.transform([read.bam for read in READS.values() if read.type == 'IP'],
                  ruffus.suffix('.bam'), '.crosslink.sites.bed')
def pureclip(bam, bed):
    ip_bam, input_bam = [[ip_read.bam, input_read.bam] for (ip_read, input_read) in SAMPLES
                         if ip_read.bam == bam][0]
    cmd = ['pureclip', '-i', ip_bam, '-bai', f'{ip_bam}.bai', '-g', f'{options.genome}/genome.fa',
           '-iv', "'chr1;chr2;chr3'", '-nt', options.cores, '-ibam', input_bam, '-ibai', f'{input_bam}.bai',
           '-o', bed, '-or', bed.replace('.crosslink.sites.bed', '.binding.regions.bed'),
           '>', bed.replace('.crosslink.sites.bed', '.pureclip.log')]
    cmding(cmd, message=f'Calling peaks from {bam} {size(bam)} using pureCLIP ...')


def mapped_read_count(bam):
    with pysam.AlignmentFile(bam, 'rb') as sam:
        return sam.mapped
    
    
@ruffus.jobs_limit(PROCESSES)
@ruffus.follows(pureclip)
@ruffus.transform(clipper, ruffus.suffix('.peak.clusters.bed'), '.peak.clusters.normalized.bed')
def overlap_peaks(peak, norm_peak):
    ip_bam, input_bam = [[ip_read.bam, input_read.bam] for (ip_read, input_read) in SAMPLES
                         if ip_read.bed == peak][0]
    cmd = ['overlap_peak.pl', ip_bam, input_bam, peak,
           mapped_read_count(ip_bam), mapped_read_count(input_bam), norm_peak]
    cmding(cmd, message=f'Normalizing peaks in {peak} {size(peak)} ...')


@ruffus.transform(overlap_peaks, ruffus.suffix('.normalized.bed'), '.normalized.compressed.bed')
def compress_peaks(bed, output):
    bed = bed.replace('.bed', '.full.bed')
    cmd = ['compress_peak.pl', bed, output]
    cmding(cmd, message=f'Compressing peaks in {bed} {size(bed)} ...')


@ruffus.transform(compress_peaks, ruffus.suffix('.normalized.compressed.bed'), '.entropy.bed')
def calculate_entropy(bed, output):
    message, start_time = f'Calculating entropy for {bed} ...', time.perf_counter()
    logger.info(message)
    peak_bed = bed.replace('.normalized.compressed.bed', '.bed')
    ip_bam, input_bam = [[ip_read.bam, input_read.bam] for (ip_read, input_read) in SAMPLES
                         if ip_read.bed == peak_bed][0]
    ip_mapped_read_count, input_mapped_read_count = mapped_read_count(ip_bam), mapped_read_count(input_bam)
    bed = bed.replace('.bed', '.full.bed')
    columns = ['chrom', 'start', 'end', 'peak', 'ip_read_number', 'input_read_number',
               'p', 'v', 'method', 'status', 'l10p', 'l2fc']
    df = pd.read_csv(bed, sep='\t', header=None, names=columns)
    df['pi'] = df['ip_read_number'] / ip_mapped_read_count
    df['qi'] = df['input_read_number'] / input_mapped_read_count
    df['entropy'] = df.apply(lambda row: 0 if row.pi <= row.qi else row.pi * math.log2(row.pi / row.qi), axis=1)
    df['excess_reads'] = df['pi'] - df['qi']
    entropy = output.replace('.entropy.bed', '.entropy.tsv')
    df.to_csv(entropy, index=False, columns=columns + ['entropy'], sep='\t', header=False)
    excess_read = output.replace('.entropy.tsv', 'excess.reads.entropy.tsv')
    df.to_csv(excess_read, index=False, columns=columns + ['excess_reads'], sep='\t', header=False)
    df['strand'] = df.peak.str.split(':', expand=True)[2]
    df['l2fc'] = df['l2fc'].map('{:.15f}'.format)
    df['entropy'] = df['entropy'].map('{:.10f}'.format)
    # For IDR 2.0.2, columns 'excess_reads', 'pi', and 'qi' need to be excluded for .entropy.bed
    # For IDR 2.0.3, columns 'excess_reads', 'pi', and 'qi' need to be retained for .entropy.bed
    columns = ['chrom', 'start', 'end', 'l2fc', 'entropy', 'strand', 'excess_reads', 'pi', 'qi']
    df.to_csv(output, index=False, columns=columns, sep='\t', header=False)
    run_time = int(time.perf_counter() - start_time)
    message = message.replace(' ...', f' completed in [{str(datetime.timedelta(seconds=run_time))}].')
    logger.info(message)
    
    
@ruffus.merge(calculate_entropy, ECLIP + '/01v02.idr.out.bed')
def idr_peaks(inputs, output):
    (peak1, peak2), idr_out = inputs, output.replace('.bed', '')
    cmd = ['idr', '--sample', peak1, peak2, '--input-file-type', 'bed', '--rank', 5,
           '--peak-merge-method', 'max', '--plot', '-o', idr_out]
    cmding(cmd, message=f'Running IDR to rank peaks in {peak1} and\n{" " * 40}{peak2} ...')
    cmd = ['parse_idr_peaks.pl', idr_out, peak1, peak2, output]
    cmding(cmd, message=f'Parsing IDR peaks in {idr_out} ...')


def overlap_idr_peaks(idr_peak_bed, bed):
    for ip_read, input_read in SAMPLES:
        cmd = ['overlap_peak.pl', ip_read.bam, input_read.bam, idr_peak_bed,
               mapped_read_count(ip_read.bam), mapped_read_count(input_read.bam), bed]
        cmding(cmd, message=f'Normalizing IDR peaks in {idr_peak_bed} {size(idr_peak_bed)} ...')


def reproducible_peaks():
    pass

    
def sort_compressed_peaks():
    pass


def remove_peaks_on_blacklist():
    pass


def narrow_peaks():
    pass


def make_big_bed():
    pass


def rescue_ratio():
    pass


def consistency_ratio():
    pass


def cleanup():
    logger.info('Deleting soft links ...')
    cmding(f'rm {ECLIP}/*.fastq.gz')
    logger.info(f'Deleting {"UMI extract" if TYPE == "single" else "barcodes demultiplex"} files ...')
    cmding(f'rm {ECLIP}/*.fq.gz {ECLIP}/*.prepare.done')
    logger.info('Deleting cut adapter files ...')
    cmding(f'rm {ECLIP}/*.clean.fastq')
    logger.info('Compressing sorted clean fastq files ...')
    for fastq in glob.iglob(f'{ECLIP}/*.clean.sort.fastq'):
        cmding(f'pigz -p 8 {fastq}')
    logger.info('Compressing sorted mate files ...')
    for mate in glob.iglob(f'{ECLIP}/repeat.elements.map/*/Unmapped.out.mate1'):
        sample = os.path.basename(os.path.dirname(mate))
        cmding(f'pigz -p 8 -c {mate} > {ECLIP}/{sample}.mate1.sort.fastq.gz')
    for mate in glob.iglob(f'{ECLIP}/repeat.elements.map/*/Unmapped.out.mate2'):
        sample = os.path.basename(os.path.dirname(mate))
        cmding(f'pigz -p 8 -c {mate} > {ECLIP}/{sample}.mate2.sort.fastq.gz')
    logger.info('Deleting map directories ...')
    cmding(f'rm -rf {ECLIP}/repeat.elements.map {ECLIP}/reference.genome.map')
    logger.info('Deleting BAM dedup temporary files ...')
    cmding(f'rm {ECLIP}/*.dedup.*')
    logger.info('Deleting temporary BAM files ...')
    cmding(f'rm {ECLIP}/*.sort.bam')
    if TYPE == 'paired':
        cmding(f'rm {ECLIP}/*.merge.bam')


# @ruffus.posttask(cleanup)
# @ruffus.follows(ruffus.mkdir(QC))
# @ruffus.merge([DESeq2, rMATS, sort_index_bam], 'qc/rnaseq.qc.html')
def qc(inputs, output):
    logger.info('Moving cutadapt metrics to qc ...')
    cmding(f'mv {ECLIP}/*.clean.metrics {ECLIP}/*.trim.metrics {QC}/')
    logger.info('Moving STAR logs to qc ...')
    logs = glob.iglob(f'{ECLIP}/*/*/Log.final.out')
    for log in logs:
        _, group, basename, _ = log.split('/')
        cmding(f'mv {log} {QC}/{basename}.{group}.log')

    config = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'multiqc.config.yaml')
    cmd = ['multiqc', '--config', config, '--filename', output, '--force', 'qc']
    cmding(cmd, message='Running MultiQC to generating QC summary ...')


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
        d['processes'] = f'{PROCESSES} / {options.cores}'
        setting = '\n'.join([f'{k.title():>22}: {v}' for k, v in d.items() if k in keys])
        logger.debug(HEADER)
        logger.trace(f'\nRunning eclip using the following settings:\n\n{setting}\n')
        ruffus.cmdline.run(options)
        logger.debug('')
        run_time = str(datetime.timedelta(seconds=int(time.perf_counter() - START_TIME)))
        logger.trace(FOOTER.format(hh_mm_ss=f'time consumed: {run_time}'.upper().center(118)))


if __name__ == '__main__':
    main()
