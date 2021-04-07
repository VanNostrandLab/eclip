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
import datetime
import sys
import tempfile
import shutil

import pysam as pysam
import ruffus
import pandas as pd
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
options = parser.parse_args()
try:
    with open(options.MANIFEST) as f:
        manifest = json.load(f)
except IOError:
    raise ValueError(f'Manifest {options.MANIFEST} may not exist or not be file.')
setattr(options, 'outdir', options.outdir or manifest.get('outdir', '') or os.getcwd())
if not os.path.isdir(options.outdir):
    try:
        os.mkdir(options.outdir)
    except OSError:
        raise OSError(f'Failed to create outdir {options.outdir}.')
os.chdir(options.outdir)

setattr(options, 'genome', options.genome or manifest.get('genome', ''))
if not os.path.isdir(options.genome):
    raise ValueError(f'Reference genome index {options.genome} is not a directory or does not exist.')

setattr(options, 'repeat', options.repeat or manifest.get('repeat', ''))
if not os.path.isdir(options.repeat):
    raise ValueError(f'Repeat element index {options.repeat} is not a directory or does not exist.')

setattr(options, 'dataset', manifest.get('dataset', 'eclip'))
setattr(options, 'species', manifest.get('species', 'hg19'))
setattr(options, 'track', options.dataset)
setattr(options, 'track_label', options.dataset)
setattr(options, 'track_genome', options.species)

blacklist_bed = manifest.get('blacklist_bed', '')
if blacklist_bed and not os.path.isfile(blacklist_bed):
    raise ValueError(f'Blacklist bed {blacklist_bed} is not a file or does not exist.')
else:
    setattr(options, 'blacklist_bed', blacklist_bed)

barcodes_fasta = manifest.get('barcodes_fasta', '')
if barcodes_fasta:
    if os.path.isfile(barcodes_fasta):
        setattr(options, 'barcodes_fasta', barcodes_fasta)
    else:
        raise ValueError(f'Barcodes fasta {barcodes_fasta} is not a file or does not exist.')
else:
    setattr(options, 'barcodes_fasta', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'barcodes.fasta'))
setattr(options, 'randomer_length', manifest.get('randomer_length', '5'))


def validate_sample(i, sample):
    ip_dict, input_dict = {}, {}
    ip_read, input_read = sample.get('ip_read', {}), sample.get('input_read', {})
    if not ip_read:
        raise ValueError(f'Sample {i} does not have ip_read.')
    if not input_read:
        raise ValueError(f'Sample {i} does not have input_read.')
    assert 'name' in ip_read and 'name' in input_read, ValueError('Missing name in ip_read and/or input_read.')
    ip_dict['name'], input_dict['name'] = ip_read['name'], input_read['name']
    ip_dict['type'], input_dict['type'] = 'ip', 'input'
    ip_dict['link'] = f'eclip/{options.dataset}.{ip_read["name"]}'
    input_dict['link'] = f'eclip/{options.dataset}.{input_read["name"]}'
    ip_dict['link1'] = f'eclip/{options.dataset}.{ip_read["name"]}.r1.fastq.gz'
    input_dict['link1'] = f'eclip/{options.dataset}.{input_read["name"]}.r1.fastq.gz'
    
    if 'fastq' in ip_read:
        read_type = 'single'
        if not os.path.isfile(ip_read['fastq']):
            raise ValueError(f'The ip_read fastq {ip_read["fastq"]} may not be a file or does not exist!')
        assert 'fastq' in input_read, ValueError('Single-End sample input_read missing fastq.')
        if not os.path.isfile(input_read['fastq']):
            raise ValueError(f'The input_read fastq {input_read["fastq"]} may not be a file or does not exist!')
        assert 'adapters' in ip_read and 'adapters' in input_read, ValueError('Single-End sample missing adapters.')
        if not os.path.isfile(ip_read['adapters']):
            raise ValueError(f'Adapters {ip_read["adapters"]} may not be a file or does not exist!')
        
        ip_dict['fastq1'], input_dict['fastq1'] = ip_read['fastq'], input_read['fastq']
        ip_dict['fastq2'], input_dict['fastq2'] = '', ''
        ip_dict['link2'], input_dict['link2'] = '', ''
        ip_dict['adapters'], input_dict['adapters'] = ip_read['adapters'], input_read['adapters']
        ip_dict['barcodes'], input_dict['barcodes'] = ['', ''], ['', '']
    else:
        read_type = 'paired'
        assert 'fastq1' in ip_read and 'fastq2' in ip_read, ValueError('Paired-End sample ip_read missing '
                                                                       'fastq1 and/or fastq2.')
        if not os.path.isfile(ip_read['fastq1']):
            raise ValueError(f'The ip_read fastq1 {ip_read["fastq1"]} may not be a file or does not exist!')
        if not os.path.isfile(ip_read['fastq2']):
            raise ValueError(f'The ip_read fastq2 {ip_read["fastq2"]} may not be a file or does not exist!')
        
        assert 'fastq1' in input_read and 'fastq2' in input_read, ValueError('Paired-End sample input_read missing '
                                                                             'fastq1 and/or fastq2.')
        if not os.path.isfile(input_read['fastq1']):
            raise ValueError(f'The input_read fastq1 {input_read["fastq1"]} may not be a file or does not exist!')
        if not os.path.isfile(input_read['fastq2']):
            raise ValueError(f'The input_read fastq2 {input_read["fastq2"]} may not be a file or does not exist!')
        
        assert 'barcodes' in ip_read, ValueError('Paired-End sample in_read missing barcodes.')
        assert 'barcodes' in input_read, ValueError('Paired-End sample input_read missing barcodes.')
        ip_dict['link2'] = f'eclip/{options.dataset}.{ip_read["name"]}.r2.fastq.gz'
        input_dict['link2'] = f'eclip/{options.dataset}.{input_read["name"]}.r2.fastq.gz'
        ip_dict['fastq1'], input_dict['fastq1'] = ip_read['fastq1'], input_read['fastq1']
        ip_dict['fastq2'], input_dict['fastq2'] = ip_read['fastq2'], input_read['fastq2']
        ip_dict['adapters'], input_dict['adapters'] = '', ''
        ip_dict['barcodes'], input_dict['barcodes'] = ip_read['barcodes'], input_read['barcodes']
    return ip_dict, input_dict, read_type


SAMPLES = manifest.get('samples', [])
if not SAMPLES:
    raise KeyError('Manifest file does not contain any samples.')
FASTQS, SIZES, TYPE = {}, [], set()
for i, sample in enumerate(SAMPLES, 1):
    ips, inputs, read_type = validate_sample(i, sample)
    TYPE.add(read_type)
    ip_link, input_link = ips['link'], inputs['link']
    if ip_link in FASTQS:
        raise ValueError('Duplicated names found for ip_read.')
    if input_link in FASTQS:
        raise ValueError('Duplicated names found for input_read.')
    FASTQS[ip_link], FASTQS[input_link] = ips, inputs
    SIZES.extend([os.path.getsize(ips['fastq1']), os.path.getsize(inputs['fastq1'])])
    if ips['fastq2']:
        SIZES.extend([os.path.getsize(ips['fastq2']), os.path.getsize(inputs['fastq2'])])
        
assert len(TYPE) == 1, ValueError('Manifest mixed with single-end and paired-end fastq file specification.')
TYPE = list(TYPE)[0]

VERBOSE = options.verbose
setattr(options, 'verbose', 3 if options.dry_run else 0)
setattr(options, 'jobs', options.cores)
setattr(options, 'multiprocess', options.cores)
setattr(options, 'use_threads', False)
setattr(options, 'just_print', options.dry_run)
setattr(options, 'verbose_abbreviated_path', options.verbose_path)
setattr(options, 'exceptions_terminate_immediately', True)
setattr(options, 'one_second_per_job', True)
setattr(options, 'log_exceptions', True)
setattr(options, 'logger', ruffus.black_hole_logger)

VERSION = 1.0
logo = fr"""
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
            for c in command:
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
            os.unlink(profile_output)
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
    cmd = f'time -f "%E %M" -o {profile_output} {cmd}'
    process = subprocess.Popen(cmd, universal_newlines=True, shell=True, cwd=cwd,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE, **kwargs)
    stdout, stderr = process.communicate()
    if process.returncode:
        logger.error(f'Failed to run {program} (exit code {process.returncode}):\n{stderr or stdout}')
        sys.exit(process.returncode)
    if message:
        logger.info(message.replace(' ...', f' complete [{parse_profile()}].'))


def estimate_process():
    """Estimate number of processes based on the maximum size of fastq file."""
    
    size = max(SIZES) / (1000 * 1000 * 1000) * 10
    n = int(options.memory / size)
    if n == 0:
        n = 1
    elif n > options.cores:
        n = options.cores
    return n


PROCESS = estimate_process()


@ruffus.jobs_limit(1)
@ruffus.follows(ruffus.mkdir('eclip'))
@ruffus.originate([v['link1'] for v in FASTQS.values()])
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
                    
    read = FASTQS[link1.replace('.r1.fastq.gz', '')]
    link1, link2 = read['link1'], read['link2']
    make_link(read['fastq1'], link1)
    if link2:
        make_link(read['fastq2'], link2)


def cleanup_demux():
    fastq, rm = set(), [options.dataset]
    for link, read in FASTQS.items():
        name, barcodes = read['name'], read['barcodes']
        if read['link2']:
            fastq.add(f'eclip/{options.dataset}.{name}.{barcodes[0]}.r1.fq.gz')
            fastq.add(f'eclip/{options.dataset}.{name}.{barcodes[0]}.r2.fq.gz')
            fastq.add(f'eclip/{options.dataset}.{name}.{barcodes[1]}.r1.fq.gz')
            fastq.add(f'eclip/{options.dataset}.{name}.{barcodes[1]}.r2.fq.gz')
            rm.extend([name] + barcodes)
        else:
            fastq.add(f'eclip/{options.dataset}.{name}.umi.r1.fq.gz')
    rm = [os.path.join('eclip', r) for r in rm] + [x for x in glob.iglob('eclip/*.fq.gz') if x not in fastq]
    _ = [os.unlink(x) for x in rm if os.path.isfile(x)]


@ruffus.jobs_limit(PROCESS)
@ruffus.posttask(cleanup_demux)
@ruffus.transform(soft_link, ruffus.suffix('.fastq.gz'), '.prepare.done')
def prepare_reads(link, output):
    """Extract UMIs for single-end reads or demultiplex paired-end reads."""
    read = FASTQS[link.replace('.r1.fastq.gz', '')]
    fastq1, fastq2 = read['link1'], read['link2']
    if fastq2:
        barcodes, cwd = read['barcodes'], os.path.join(options.outdir, 'eclip')
        message = f'Demultiplexing paired-end reads {fastq1} {size(fastq1)} and\n{" " * 43}{fastq2} {size(fastq2)} ...'
        cmd = ['demux',
               '--fastq_1', os.path.basename(fastq1),
               '--fastq_2', os.path.basename(fastq2),
               '--newname', read['name'],
               '--expectedbarcodeida', barcodes[0],
               '--expectedbarcodeidb', barcodes[1],
               '--dataset', options.dataset,
               '--barcodesfile', options.barcodes_fasta,
               '--length', options.randomer_length,
               '--metrics', f'{options.dataset}.{read["name"]}.demux.metrics']
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
    for link, read in FASTQS.items():
        name, barcodes = read['name'], read['barcodes']
        if TYPE == 'paired':
            fastq.add(f'eclip/{options.dataset}.{name}.{barcodes[0]}.r1.clean.fastq')
            fastq.add(f'eclip/{options.dataset}.{name}.{barcodes[1]}.r1.clean.fastq')
        else:
            fastq.add(f'eclip/{options.dataset}.{name}.umi.r1.clean.fastq')
    return fastq


def size(file):
    if os.path.isfile(file):
        s = os.path.getsize(file)
        if s < 1000:
            s = f'[{s:.2f}B]'
        elif s < 1000 * 1000:
            s = f'[{s / 1000:.2f}KB]'
        elif s < 1000 * 1000 * 1000:
            s = f'[{s / 1000 / 1000:.2f}MB]'
        else:
            s = f'[{s / 1000 / 1000 / 1000:.2f}GB]'
    else:
        s = '[?KB]'
    return s


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
        adapter, adapters = read['adapters'], []
        if adapter:
            adapters1 = parse_adapters('-a', adapter)
            adapters2 = adapters1
            overlap1, overlap2 = '1', '5'
        else:
            env = os.environ.copy()
            env['PATH'] = f'{os.path.dirname(os.path.abspath(__file__))}:{env["PATH"]}'
            cmd = ['parse_barcodes.sh', options.randomer_length, options.barcodes_fasta] + read['barcodes']
            folder = tempfile.mkdtemp(dir=os.path.join(options.outdir, 'eclip'))
            cmding(cmd, message='Parsing barcodes and finding adapters ...', cwd=folder, env=env)
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
            ios1 = ['-o', out1.replace('.clean.', '.trim.').replace('.fastq', '.fastq.gz'),
                    out1.replace('.clean.', '.').replace('.fastq', '.fq.gz'),
                    '>', out1.replace('.r1.clean.fastq', '.trim.metrics')]
            msg1 = (f"Cutting adapters for single read {out1.replace('.clean.', '.')} "
                    f"{size(out1.replace('.clean.', '.'))} (first round) ...")
            ios2 = ['-o', out1,
                    out1.replace('.clean.', '.trim.').replace('.fastq', '.fastq.gz'),
                    '>', out1.replace('.r1.clean.fastq', '.clean.metrics')]
            msg2 = (f"Cutting adapters for single read {out1.replace('.clean.', '.trim.')} "
                    f"{size({out1.replace('.clean.', '.trim.')})} (second round) ...")
        else:
            ios1 = ['-o', out1.replace('.clean.', '.trim.').replace('.fastq', '.fastq.gz'),
                    '-p', out1.replace('.clean.', '.trim.').replace('.fastq', '.fastq.gz').replace('.r1.', '.r2.'),
                    out1.replace('.clean.', '.').replace('.fastq', '.fq.gz'),
                    out1.replace('.clean.', '.').replace('.fastq', '.fq.gz').replace('.r1.', '.r2.'),
                    '>', out1.replace('.r1.clean.fastq', '.trim.metrics')]
            msg1 = (f"Cutting adapters for paired reads {out1.replace('.clean.', '.')} and\n{' ' * 45}"
                    f"{out1.replace('.clean.', '.').replace('.r1.', '.r2.')} "
                    f"{size({out1.replace('.clean.', '.').replace('.r1.', '.r2.')})} (first round) ...")
            ios2 = ['-o', out1, '-p', out1.replace('.r1.', '.r2.'),
                    out1.replace('.clean.', '.trim.').replace('.fastq', '.fastq.gz'),
                    out1.replace('.clean.', '.trim.').replace('.fastq', '.fastq.gz').replace('.r1.', '.r2.'),
                    '>', out1.replace('.r1.clean.fastq', '.clean.metrics')]
            msg2 = (f"Cutting adapters for paired reads {out1.replace('.clean.', '.trim.')} "
                    f"{size(out1.replace('.clean.', '.trim.'))} and\n{' ' * 45}"
                    f"{out1.replace('.clean.', '.trim.').replace('.r1.', '.r2.')} "
                    f"{size(out1.replace('.clean.', '.trim.').replace('.r1.', '.r2.'))} (second round) ...")
        return ios1, ios2, msg1, msg2
    
    def trim_adapters(adapters, overlap, ios, message):
        cmd = ['cutadapt', '-O', overlap, '--times', '2', '-e', '0.0', '-j', options.cores, '-m', '18',
               '--quality-cutoff', '6', '--match-read-wildcards'] + adapters + ios
        cmding(cmd, message=message)
    
    link = '.'.join(fastq.split('.')[:-4])
    read = FASTQS[link]
    adapters1, adapters2, overlap1, overlap2 = get_adapters(read)
    ios1, ios2, msg1, msg2 = get_ios(fastq)
    trim_adapters(adapters1, overlap1, ios1, msg1)
    trim_adapters(adapters2, overlap2, ios2, msg2)
    cmding('rm eclip/*.trim.fastq.gz')


@ruffus.jobs_limit(PROCESS)
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
@ruffus.follows(ruffus.mkdir('eclip/repeat.elements.map'))
@ruffus.transform(sort_fastq,
                  ruffus.formatter(r'.+/(?P<BASENAME>.*).r1.clean.sort.fastq'),
                  'eclip/repeat.elements.map/{BASENAME[0]}/Unmapped.out.mate1')
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


@ruffus.jobs_limit(PROCESS)
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
@ruffus.follows(ruffus.mkdir('eclip/reference.genome.map'))
@ruffus.transform(sort_mate,
                  ruffus.formatter(r'.+/Unmapped.out.mate1.sort.fastq'),
                  'eclip/reference.genome.map/{subdir[0][0]}/Aligned.out.bam')
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
    cmding(f'samtools sort -n -@ {options.cores} -m {min(4, int(options.memory / PROCESS))}G -o {out} {bam}',
           message=f'Name sorting {bam} {size(bam)} ...')


def position_sort_bam(bam, out):
    cmding(f'samtools sort -@ {options.cores} -m {min(4, int(options.memory / PROCESS))}G -o {out} {bam}',
           message=f'Sorting {bam} {size(bam)} ...')
    
    
def index_sorted_bam(bam):
    cmding(f'samtools index -@ {options.cores} {bam}', message=f'Indexing {bam} {size(bam)} ...')


def merge_paired_bam(bam, out):
    if not os.path.exists(out):
        link = out.replace(".merge.bam", "")
        barcodes = FASTQS[link]['barcodes']
        if barcodes[0] == 'NIL':
            cmding(f'mv {bam} {out}')
        else:
            b1, b2 = barcodes
            if b1 in bam:
                b1, b2 = bam, bam.replace(b1, b2)
            else:
                b1, b2 = bam.replace(b2, b1), bam
            message = f'Merging {b1} {size(b1)} and {b2} {size(b2)} ...'
            cmding(f'samtools merge -@ {options.cores} -m {min(4, int(options.memory / PROCESS))}G {out} {b1} {b2}',
                   message=message)


@ruffus.jobs_limit(1)
@ruffus.transform(map_to_reference_genome, ruffus.formatter(r'.+/Aligned.out.bam'), 'eclip/{subdir[0][0]}.sort.bam')
def prepare_bam(bam, out):
    if TYPE == 'single':
        name_sort = out.replace('.sort.bam', '.name.sort.bam')
        name_sort_bam(bam, name_sort)
        position_sort_bam(name_sort, out)
        index_sorted_bam(out)
        cmding(f'rm {name_sort}')
    else:
        name_sort_bam(bam, out)


@ruffus.jobs_limit(PROCESS)
@ruffus.transform(prepare_bam, ruffus.suffix('.sort.bam'), '.sort.dedup.bam')
def dedup_bam(bam, out):
    """Collapse barcodes of paired-end bam or umi_tools dedup single-end bam."""
    # env = os.environ.copy()
    if TYPE == 'single':
        cmd = ['umi_tools', 'dedup', '--random-seed', 1, '--stdin', bam, '--method', 'unique',
               '--output-stats', out.replace(".dedup.bam", ".dedup.metrics"), '--stdout', out]
        message = f'Deduplicating {bam} {size(bam)} by umi_tools dedup ...'
    else:
        # env['PATH'] = f'{os.path.dirname(os.path.abspath(__file__))}:{env["PATH"]}'
        cmd = f'barcode_collapse.py -o {out} -m {out.replace(".bam", ".collapse.metrics")} -b {bam}'
        message = f'Deduplicating {bam} {size(bam)} by collapsing barcodes ...'
    # cmding(cmd, message=message, env=env)
    cmding(cmd, message=message)


@ruffus.transform(dedup_bam, ruffus.suffix('.sort.dedup.bam'), '.sort.dedup.sort.bam')
def sort_bam(bam, out):
    position_sort_bam(bam, out)


@ruffus.jobs_limit(1)
@ruffus.transform(sort_bam, ruffus.regex(r"^(.+)\.(.+)\.sort.dedup.sort.bam$"), r'\1.merge.bam')
def merge_bam(bam, out):
    if TYPE == 'single':
        cmding(f'mv {bam} {out}')
    else:
        merge_paired_bam(bam, out)


@ruffus.jobs_limit(1)
@ruffus.transform(merge_bam, ruffus.suffix('.merge.bam'), '.bam')
def index_bam(bam, out):
    if TYPE == 'paired':
        cmding(f'samtools view -f 128 -@ {options.cores} -b -o {out} {bam}',
               message='Extracting r2 reads from {bam} {size(bam)} ...')
    else:
        cmding(f'mv {bam} {out}')
    if not os.path.exists(f'{bam}.bai'):
        index_sorted_bam(out)
        

@ruffus.jobs_limit(PROCESS)
@ruffus.follows(ruffus.mkdir('hub'))
@ruffus.follows(index_bam)
@ruffus.transform([f'{link}.bam' for link in FASTQS],
                  ruffus.formatter(r'.+/(?P<BASENAME>.*).bam'), 'hub/{BASENAME[0]}.plus.bw')
def make_bigwig_files(bam, bigwig):
    def bam_to_bigwig(bam, scale, strand, bw):
        bg, bg_sort = f'{bw[:-3]}.bg', f'{bw[:-3]}.sort.bg'
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
        link = bigwig.replace('.plus.bw', '').replace('hub/', 'eclip/')
        read_type = 'paired' if FASTQS[link]['link2'] else 'single'
        pos_bw, neg_bw = bigwig, bigwig.replace('.plus.bw', '.minus.bw')
        with pysam.AlignmentFile(bam, 'rb') as sam:
            total_reads = sam.mapped
        total_reads = total_reads if read_type == 'single' else total_reads / 2
        if total_reads:
            scale = 1000000.0 / total_reads
        else:
            logger.error(f'No reads was found in BAM {bam}')
            ruffus.touch_file(bigwig)
            return
        if read_type == 'single':
            bam_to_bigwig(bam, scale, '+', pos_bw)
            bam_to_bigwig(bam, -1 * scale, '-', neg_bw)
        else:
            bam_to_bigwig(bam, -1 * scale, '-', pos_bw)
            bam_to_bigwig(bam, scale, '+', neg_bw)
        run_time = int(time.perf_counter() - start_time)
        message = message.replace(' ...', f' completed in [{str(datetime.timedelta(seconds=run_time))}].')
        logger.info(message)


@ruffus.merge(make_bigwig_files, 'hub/hub.txt')
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


@ruffus.jobs_limit(PROCESS)
@ruffus.follows(make_hub_files)
@ruffus.transform([f'{k}.bam' for k, v in FASTQS.items() if v['type'] == 'ip'],
                  ruffus.suffix('.bam'), '.peak.clusters.bed')
def clipper(bam, bed):
    cmd = f'clipper --species {options.species} --bam {bam} --outfile {bed}'
    cmding(cmd, message=f'Calling peaks from {bam} {size(bam)} ...')


@ruffus.jobs_limit(PROCESS)
@ruffus.transform(clipper, ruffus.suffix('.peak.clusters.bed'), '.peak.clusters.normalized.bed')
def overlap_peaks(peak, norm_peak):
    def peak_to_bams(peak):
        for sample in SAMPLES:
            ip_name = [s['name'] for s in sample if 'ip_read' in s][0]
            input_name = [s['name'] for s in sample if 'input_read' in s][0]
            ip_bam = f'eclip/{options.dataset}.{ip_name}.bam'
            input_bam = f'eclip/{options.dataset}.{input_name}.bam'
            if peak == f'eclip/{options.dataset}.{ip_name}.peak.clusters.bed':
                return ip_bam, input_bam
        else:
            raise ValueError(f'Failed to find corresponding bam files for {peak}.')
        
    def mapped_read_number(bam):
        cmd = f'samtools view -c -F 4 -o {bam}.mapped.read.number {bam}'
        cmding(cmd, message=f'Counting mapped read number of {bam} {size(bam)} ...')
        
    ip_bam, input_bam = peak_to_bams(peak)
    ip_count, input_count = mapped_read_number(ip_bam), mapped_read_number(input_bam)
    cmd = ['overlap_peak.pl', ip_bam, input_bam, peak, ip_count, input_count, norm_peak]
    cmding(cmd, message=f'Normalizing peaks in {peak} {size(peak)} ...')
    cmding(f'rm {ip_count} {input_count}')
    

def calculate_entropy():
    pass


def compress_peaks():
    pass


# def cleanup():
#     logger.info('Deleting soft links ...')
#     cmding('rm fastq_to_bam/*.fastq.gz')
#     logger.info('Deleting fastq files ...')
#     cmding('rm fastq_to_bam/*.clean.fastq')
#     logger.info('Compressing sorted fastq files ...')
#     for fastq in glob.iglob('fastq_to_bam/*.clean.sort.fastq'):
#         cmding(f'pigz -p 8 {fastq}')
#     logger.info('Compressing sorted mate files ...')
#     for mate in glob.iglob('fastq_to_bam/repeat.elements.map/*/Unmapped.out.mate1'):
#         sample = os.path.basename(os.path.dirname(mate))
#         cmding(f'pigz -p 8 -c {mate} > fastq_to_bam/{sample}.mate1.sort.fastq.gz')
#     for mate in glob.iglob('fastq_to_bam/repeat.elements.map/*/Unmapped.out.mate2'):
#         sample = os.path.basename(os.path.dirname(mate))
#         cmding(f'pigz -p 8 -c {mate} > fastq_to_bam/{sample}.mate2.sort.fastq.gz')
#     logger.info('Deleting map directories ...')
#     cmding('rm -rf fastq_to_bam/repeat.elements.map fastq_to_bam/reference.genome.map')
#     logger.info('Deleting rMATS temporary directory ...')
#     cmding('rm alternative_splicing/rmats.summary.txt')
#     cmding('rm -rf alternative_splicing/tmp')
#
#
# @ruffus.posttask(cleanup)
# @ruffus.follows(ruffus.mkdir('qc'))
# @ruffus.merge([DESeq2, rMATS, sort_index_bam], 'qc/rnaseq.qc.html')
# def qc(inputs, output):
#     logger.info('Moving cutadapt metrics to qc ...')
#     for metric in glob.iglob('fastq_to_bam/*.cutadapt.metrics'):
#         cmding(f'mv {metric} qc/')
#     logger.info('Moving STAR logs to qc ...')
#     logs = glob.iglob('fastq_to_bam/*/*/Log.final.out')
#     for log in logs:
#         _, group, basename, _ = log.split('/')
#         cmding(f'mv {log} qc/{basename}.{group}.log')
#     logger.info('Moving feature count summary to qc ...')
#     cmding(f'mv differential_expression/feature.count.summary qc/')
#
#     config = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'multiqc.config.yaml')
#     cmd = ['multiqc', '--config', config, '--filename', output, '--force', 'qc']
#     cmding(cmd, message='Running MultiQC to generating QC summary ...')
#
#
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


@logger.catch()
def main():
    if options.scheduler:
        schedule()
    else:
        keys = ('MANIFEST', 'outdir', 'genome', 'repeat', 'processes')
        d = vars(options).copy()
        d['processes'] = f'{PROCESS} / {options.cores}'
        setting = '\n'.join([f'{k.title():>22}: {v}' for k, v in d.items() if k in keys])
        logger.debug(logo)
        logger.trace(f'Running eclip using the following settings:\n\n{setting}\n')
        try:
            ruffus.cmdline.run(options)
        except OSError:
            pass
        finally:
            logger.debug('')
            run_time = int(time.perf_counter() - START_TIME)
            logger.info(f'Mission Accomplished!')
            logger.trace(f'Time consumed: {str(datetime.timedelta(seconds=run_time))}.')


if __name__ == '__main__':
    main()
