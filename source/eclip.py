#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
A pipeline for processing eCLIP data to identify genomic locations of RNA binding proteins (RBPs).
"""

import argparse
import os
import glob
import gzip
import sys
import time
import subprocess
import json
from collections import defaultdict

import yaml
import datetime
import tempfile
import shutil
import itertools
import math

import pysam
import cmder
import pandas as pd
from seqflow import Flow, task, logger
from Bio.SeqIO.QualityIO import FastqGeneralIterator


class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass


epilog = """
Note:
For sample process argument, they can be set either by command line argument using --argument_name or manifest
entry using argument_name. When set in both places,the command line argument will always take precedence.
Argument outdir will use current work directory if neither the command line argument --outdir nor the manifest outdir
entry was set.
"""
parser = argparse.ArgumentParser(description=__doc__, prog='eclip', epilog=epilog,
                                 formatter_class=CustomFormatter)
sp = parser.add_argument_group('Sample Process')
sp.add_argument('MANIFEST', type=str, help='Path to the manifest file that specifies paths for eCLIP data.')
sp.add_argument('--dataset', type=str, help='Name of the eCLIP dataset.')
sp.add_argument('--outdir', type=str, help="Path to the output directory.")
sp.add_argument('--genome', type=str, help="Path to STAR reference genome index directory.")
sp.add_argument('--repeat', type=str, help="Path to STAR repeat elements index directory.")
sp.add_argument('--species', type=str, help="Species name (short name code) the dataset associated with.")
sp.add_argument('--barcodes_fasta', type=str, help="Path to the fasta file contains barcodes and their sequences.")
sp.add_argument('--blacklist_bed', type=str, help="Path to the bed file contains blacklist.")
sp.add_argument('--randomer_length', type=int, help="Length (int) of the randomer.")
sp.add_argument('--allow_mismatch', type=int, help="Allowed mismatch among barcodes for demultiplexing.")
sp.add_argument('--track', type=str, help="Name for the UCSC Genome Browser track, default: eCLIP", default='eCLIP')
sp.add_argument('--track_label', type=str, help="Label for the track.", default='eCLIP')
sp.add_argument('--track_genome', type=str, help="Genome name (a short name code) for the track.", default='hg19')
sp.add_argument('--l2fc', type=int, help="Only consider peaks at or above this log2 fold change cutoff.", default=3)
sp.add_argument('--l10p', type=int, help="Only consider peaks at or above this log10p value cutoff.", default=3)
sp.add_argument('--enrichment_filter', type=int, help="Pre-filter peaks that are enriched over input.", default=0)

js = parser.add_argument_group('Job Submit')
js.add_argument('--job', type=str, help="Name of your job", default='eCLIP')
js.add_argument('--email', type=str, help='Email address for notifying you the start, end, and abort of you job.')
js.add_argument('--time', type=int, help='Time (in integer hours) for running your job.', default=36)
js.add_argument('--memory', type=int, help='Amount of memory (in GB) for all CPU cores.', default=32)
js.add_argument('--cores', type=int, help='Maximum number of CPU cores can be used for your job.', default=16)
js.add_argument('--scheduler', type=str, choices=('pbs', 'qsub', 'slurm', 'sbatch'),
                help='Name (case insensitive) of the scheduler on your cluster.')
js.add_argument('--hold_submit', action='store_true',
                help='Generate the submit script but hold it without submitting to the job scheduler.')

op = parser.add_argument_group('Options')
op.add_argument('--verbose', action='store_true', help='Print out detailed processing messages.')
op.add_argument('--debug', action='store_true', help='Invoke debug mode.')
op.add_argument('--dry_run', action='store_true',
                help='Print out steps and files involved in each step without actually running the pipeline.')
op.add_argument('--task_list', action="store_true",
                help='Print out task list of the pipeline without actually running the pipeline.')
op.add_argument('--target_tasks', metavar='TASKS', action="append", type=str, default=[],
                help='Target task(s) of pipeline.')
op.add_argument('--forced_tasks', metavar='TASKS', action="append", type=str, default=[],
                help='Task(s) which will be included and force re-run even if they are already up to date.')

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
    setattr(options, 'outdir', options.outdir or manifest.get('outdir', os.getcwd()))
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
    setattr(options, 'barcodes_fasta', options.barcodes_fasta or manifest.get('barcodes_fasta', barcodes_fasta))
    if not os.path.isfile(options.barcodes_fasta):
        raise ValueError(f'Barcodes fasta {options.barcodes_fasta} is not a file or does not exist.')

    setattr(options, 'randomer_length', options.randomer_length or int(manifest.get('randomer_length', 10)))
    setattr(options, 'allow_mismatch', options.allow_mismatch or int(manifest.get('allow_mismatch', 1)))

    blacklist_bed = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'hg19.blacklist.bed')
    setattr(options, 'blacklist_bed', options.blacklist_bed or manifest.get('blacklist_bed', blacklist_bed))
    if not os.path.isfile(options.blacklist_bed):
        raise ValueError(f'Blacklist bed {options.blacklist_bed} is not a file or does not exist.')

    setattr(options, 'dataset', options.dataset or manifest.get('dataset', 'eclip'))
    setattr(options, 'species', options.species or manifest.get('species', 'hg19'))
    setattr(options, 'track', options.track or options.dataset or 'eCLIP')
    setattr(options, 'track_label', options.track or options.dataset or 'eCLIP')
    setattr(options, 'track_genome', options.species or manifest.get('species', 'hg19'))

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
        # if self.fastq1:
        #     if not os.path.isfile(self.fastq1):
        #         raise ValueError(f'{self.fastq1} may not be a file or does not exist.')
        # else:
        #     raise KeyError(f'No fastq or fastq1 was assigned for {self.type} read\n{record}.')
        self.fastq2 = record.get('fastq2', '')
        # if self.fastq2 and not os.path.isfile(self.fastq2):
        #     raise ValueError(f'{self.fastq2} may not be a file or does not exist.')
        self.barcodes = record.get('barcodes', ['', ''])
        self.adapters = record.get('adapters', '')

        self.key = self.name
        self.paired = True if self.fastq2 else False
        self.link1 = f'{ECLIP}/{self.key}.r1.fastq.gz'
        self.link2 = f'{ECLIP}/{self.key}.r2.fastq.gz' if self.fastq2 else ''
        self.bam = f'{ECLIP}/{self.key}.bam'
        self.bed = f'{ECLIP}/{self.key}.peak.clusters.bed' if read_type == 'ip' else ''
        self.full_bed = f'{ECLIP}/{self.key}.peak.clusters.full.bed' if read_type == 'ip' else ''
        self.normalized_bed = f'{ECLIP}/{self.key}.peak.clusters.normalized.bed' if read_type == 'ip' else ''
        self.normalized_full_bed = f'{ECLIP}/{self.key}.peak.clusters.normalized.full.bed' if read_type == 'ip' else ''
        self.compressed_bed = f'{ECLIP}/{self.key}.peak.clusters.compressed.bed' if read_type == 'ip' else ''
        self.compressed_full_bed = f'{ECLIP}/{self.key}.peak.clusters.compressed.full.bed' if read_type == 'ip' else ''
        self.entropy_bed = f'{ECLIP}/{self.key}.peak.clusters.entropy.bed' if read_type == 'ip' else ''
        self.entropy_full_bed = f'{ECLIP}/{self.key}.peak.clusters.entropy.full.bed' if read_type == 'ip' else ''


class Sample:
    def __init__(self, ip_read, input_read):
        self.ip_read = ip_read
        self.input_read = input_read
        self.key = self.ip_read.key
        self.bed = self.ip_read.bed


def size(file, formatting=True):
    if os.path.isfile(file):
        size_bytes = os.path.getsize(file)
        if formatting:
            if size_bytes == 0:
                return '0B'
            size_name = ("B", "KB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB")
            i = int(math.floor(math.log(size_bytes, 1024)))
            size_bytes = size_bytes / math.pow(1024, i)
            size_bytes = f'[{size_bytes:.2f}{size_name[i]}]'
            # if s < 1024:
            #     s = f'[{s:.2f}B]'
            # elif s < 1024 * 1024:
            #     s = f'[{s / 1024:.2f}KB]'
            # elif s < 1024 * 1024 * 1024:
            #     s = f'[{s / 1024 / 1024:.2f}MB]'
            # else:
            #     s = f'[{s / 1024 / 1024 / 1024:.2f}GB]'
    else:
        size_bytes = '[?KB]' if formatting else 0
    return size_bytes


def parse_and_validate_samples():
    def estimate_max_processes():
        # max_size = max(sizes) / (1024 * 1024 * 1024) * 2
        # n = int(options.memory / max_size)
        # if n == 0:
        #     n = 1
        # elif n > options.cores:
        #     n = options.cores
        return options.cores

    samples, reads, sizes, types, names = [], {}, [], [], set()
    for i, sample in enumerate(options.manifest.get('samples', []), 1):
        ip_read, input_read = sample.get('ip_read', {}), sample.get('input_read', {})
        if not ip_read:
            raise KeyError(f'No key ip_read was found in sample\n{sample}.')
        if not input_read:
            raise KeyError('No key input_read was found in sample\n{sample}.')
        ip_read, input_read = Read(ip_read, 'ip'), Read(input_read, 'input')
        if ip_read.key in names:
            raise ValueError(f'Duplicated name {ip_read.key} found for ip_read\n{ip_read}.')
        else:
            names.add(ip_read.key)
        if input_read.key in names:
            raise ValueError(f'Duplicated name {input_read.key} found for input_read\n{input_read}.')
        else:
            names.add(input_read.key)
        samples.append(Sample(ip_read, input_read))
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


@task(inputs=[read.fastq1 for read in READS.values()], outputs=[read.link1 for read in READS.values()],
      mkdir_before_run=['eclip'])
def soft_link(fastq1, link1):
    """ Create soft links for original fastq files. """

    def make_link(path, link):
        if path == os.path.abspath(link):
            message = "No symbolic link was made for {path}! You are directly working on the original file!"
            logger.warning(message)
        else:
            cmder.run(f'ln -s {path} {link}', msg=f'Soft link fastq: {os.path.basename(path)} ...')

    read = READS[os.path.basename(link1.replace('.r1.fastq.gz', ''))]
    link1, link2 = read.link1, read.link2
    make_link(read.fastq1, link1)
    if link2:
        make_link(read.fastq2, link2)
    return link1


def umi_extract_or_barcode_demux(link):
    key = os.path.basename(link).replace('.r1.fastq.gz', '')
    read = READS[key]
    barcode = read.barcodes[0] if read.paired else 'umi'
    return f'{ECLIP}/{key}.{barcode}.r1.fastq.gz'


def demux(fastq1, fastq2, basename, barcodes):
    """Demultiplex paired-end reads."""
    def hamming(key, barcode, seq, allow_mismatch):
        mismatch = len(barcode) - sum(x == y or x == 'N' or y == 'N' for x, y in zip(barcode, seq))
        return (key, len(barcode), mismatch) if mismatch <= allow_mismatch else None

    logger.info(f'Demultiplexing {fastq1} and {fastq2} with barcodes {" and ".join(barcodes)} ...')
    barcodes_dict = {'A01': 'AAGCAAT',
                     'A03': 'ATGACCNNNNT',
                     'A04': 'CAGCTTNNNNT',
                     'B06': 'GGCTTGT',
                     'C01': 'ACAAGTT',
                     'D8f': 'TGGTCCT',
                     'F05': 'GGATACNNNNT',
                     'G07': 'TCCTGTNNNNT',
                     'X1A': 'NNNNNCCTATAT',
                     'X1B': 'NNNNNTGCTATT',
                     'X2A': 'NNNNNTATACTT',
                     'X2B': 'NNNNNATCTTCT'}
    allow_mismatch, randomer_length = options.allow_mismatch, options.randomer_length
    print(allow_mismatch, randomer_length)
    max_barcode_length = max(len(barcode) for barcode in barcodes_dict.values())
    writers = {barcode: (gzip.open(f'{basename}.{barcode}.r1.fastq.gz', 'wt'),
                         gzip.open(f'{basename}.{barcode}.r2.fastq.gz', 'wt'))
               for barcode in set(barcodes)}
    with gzip.open(fastq1, 'rt') as f1, gzip.open(fastq2, 'rt') as f2:
        for i, (read1, read2) in enumerate(zip(FastqGeneralIterator(f1), FastqGeneralIterator(f2))):
            (name1, seq1, quality1), (name2, seq2, quality2) = read1, read2
            n1, n2 = name1.split()[0], name2.split()[0]
            assert n1 == n2, ValueError(f'Paired-End reads have mismatch names: {name1} != {name2}')

            matches = (hamming(key, barcode, seq1[:max_barcode_length], allow_mismatch)
                       for key, barcode in barcodes_dict.items())
            matches = [match for match in matches if match]
            if matches:
                barcode, barcode_length, _ = sorted(matches, key=lambda x: x[2])[0]
                r1 = f'@{seq2[:randomer_length]}:{name1}\n{seq1[barcode_length:]}\n+\n{quality1[barcode_length:]}\n'
            else:
                barcode = 'NIL'
                r1 = f'@{seq2[:randomer_length]}:{name1}\n{seq1}\n+\n{quality1}\n'
            r2 = f'@{seq2[:randomer_length]}:{name2}\n{seq2[randomer_length:]}\n+\n{quality2[randomer_length:]}\n'

            if barcode in writers:
                writer1, writer2 = writers[barcode]
                writer1.write(r1)
                writer2.write(r2)
            if i > 50000:
                break
    _ = [[v[0].close(), v[1].close()] for v in writers.values()]
    logger.info(f'Demultiplexing {fastq1} and {fastq2} with barcodes {" and ".join(barcodes)} complete.')


@task(inputs=soft_link, outputs=lambda i: umi_extract_or_barcode_demux(i), processes=options.cores)
def prepare_reads(link, output):
    """Extract UMIs for single-end reads or demultiplex paired-end reads."""
    read = READS[os.path.basename(link.replace('.r1.fastq.gz', ''))]
    fastq1, fastq2 = read.link1, read.link2
    if fastq2:
        demux(fastq1, fastq2, fastq1.replace('.r1.fastq.gz', ''), read.barcodes)
    else:
        message = f'Extract UMIs for single-end read {fastq1} {size(fastq1)} ...'
        cmd = ['umi_tools', 'extract',
               '--random-seed', 1,
               '--stdin', fastq1,
               '--bc-pattern', 'NNNNNNNNNN',
               '--log', fastq1.replace('.fastq.gz', '.extract.metrics'),
               '--stdout', fastq1.replace('.r1.fastq.gz', '.umi.r1.fastq.gz')]
        cmder.run(cmd, msg=message, pmt=True)


def umi_extract_or_barcode_demux_outputs():
    outputs = set()
    for r in READS.values():
        barcodes = r.barcodes if r.paired else ['umi', 'umi']
        for barcode in barcodes:
            outputs.add(f'{ECLIP}/{r.key}.{barcode}.r1.fastq.gz')
    return list(outputs)


@task(inputs=umi_extract_or_barcode_demux_outputs(), outputs=lambda i: i.replace('.r1.fastq.gz', '.trim.r1.fastq'),
      parent=prepare_reads)
def cut_adapt(r1, fastq):
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
            folder = tempfile.mkdtemp(dir=ECLIP)
            cmder.run(cmd, msg='Parsing barcodes and finding adapters ...', cwd=folder)
            adapters = parse_adapters('-g', os.path.join(folder, 'g_adapters.fasta'))
            adapters += parse_adapters('-A', os.path.join(folder, 'A_adapters.fasta'))
            adapters += parse_adapters('-a', os.path.join(folder, 'a_adapters.fasta'))
            adapters1 = adapters
            adapters2 = parse_adapters('-A', os.path.join(folder, 'A_adapters.fasta'))
            overlap1 = parse_overlap(os.path.join(folder, 'trim_first_overlap_length.txt'))
            overlap2 = parse_overlap(os.path.join(folder, 'trim_again_overlap_length.txt'))
            shutil.rmtree(folder)
        return adapters1, adapters2, overlap1, overlap2

    def get_ios(input1, input2, out1, out2):
        if '.umi.' in out1:
            tmp = out1.replace('.trim.', '.clean.')
            ios1 = ['-o', tmp, input1, '>', input1.replace('.r1.fastq', '.trim.first.metrics')]
            msg1 = f"Cutting adapters for single read {input1} {size(input1)} (first round) ..."

            ios2 = ['-o', out1, tmp, '>', input1.replace('.r1.fastq', '.trim.second.metrics')]
            msg2 = f"Cutting adapters for single read {tmp} {size(tmp)} (second round) ..."
        else:
            tmp1, tmp2 = out1.replace('.trim.', '.clean.'), out2.replace('.trim.', '.clean.')
            tmp_metrics = input1.replace('.r1.fastq.gz', '.trim.first.metrics')
            ios1 = ['-o', tmp1, '-p', tmp2, input1, input2, '>', tmp_metrics]
            msg1 = (f"Cutting adapters for paired reads {input1} {size(input1)} and\n{' ' * 45}"
                    f"{input2} {size(input2)} (first round) ...")

            metrics = input1.replace('.r1.fastq.gz', '.trim.second.metrics')
            ios2 = ['-o', out1, '-p', out2, tmp1, tmp2, '>', metrics]
            msg2 = (f"Cutting adapters for paired reads {tmp1} and\n{' ' * 45}"
                    f"{tmp2} (second round) ...")
        return ios1, ios2, msg1, msg2

    def trim_adapters(adapters, overlap, ios, message):
        cmd = ['cutadapt', '-O', overlap, '-j', options.cores, '--match-read-wildcards', '--times', 1,
               '-e', 0.1, '--quality-cutoff', 6, '-m', 18] + adapters + ios
        cmder.run(cmd, msg=message, pmt=True)

    key = os.path.basename(r1).rsplit('.', maxsplit=4)[0]
    read = READS[key]
    adapters1, adapters2, overlap1, overlap2 = get_adapters(read)
    ios1, ios2, msg1, msg2 = get_ios(r1, r1.replace('.r1.', '.r2.') if read.paired else '',
                                     fastq, fastq.replace('.r1.', '.r2.') if read.paired else '')
    trim_adapters(adapters1, overlap1, ios1, msg1)
    trim_adapters(adapters2, overlap2, ios2, msg2)
    return fastq


@task(inputs=cut_adapt, outputs=lambda i: i.replace('.trim.r1.', '.trim.sort.r1.'), processes=options.cores)
def sort_fastq(fastq, output):
    tmp = tempfile.mkdtemp(suffix='_sort', prefix='fastq_', dir=ECLIP)
    cmd = f'fastq-sort --id --temporary-directory {tmp} -S 2G {fastq} > {output}'
    try:
        cmder.run(cmd, msg=f'Sorting {fastq} {size(fastq)} ...', fmt_cmd=False, pmt=True)
    finally:
        shutil.rmtree(tmp)

    fastq2 = fastq.replace('.r1.', '.r2.')
    if os.path.isfile(fastq2):
        output2 = output.replace('.r1.', '.r2.')
        tmp = tempfile.mkdtemp(suffix='_sort', prefix='fastq_', dir=ECLIP)
        cmd = f'fastq-sort --id --temporary-directory {tmp} -S 2G {fastq2} >  {output2}'
        try:
            cmder.run(cmd, msg=f'Sort fastq {fastq2} {size(fastq2)} ...', fmt_cmd=False)
        finally:
            shutil.rmtree(tmp)
    return output


@task(inputs=[i.replace('.r1.', '.trim.sort.r1.').replace('.gz', '') for i in umi_extract_or_barcode_demux_outputs()],
      parent=sort_fastq, mkdir_before_run=[f'{ECLIP}/repeat.elements.map'],
      outputs=lambda i: (f'{ECLIP}/repeat.elements.map/{os.path.basename(i).replace(".trim.sort.r1.fastq", "")}'
                         f'/Unmapped.out.mate1'))
def map_to_repeat_elements(fastq, mate1):
    fastq1, fastq2 = fastq, fastq.replace('.r1.', '.r2.')
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
    cmder.run(cmd, msg=message, pmt=True)
    return mate1


@task(map_to_repeat_elements, lambda i: i.replace('.out.mate1', '.out.sort.mate1'), processes=options.cores)
def sort_mate(mate1, output):
    tmp = tempfile.mkdtemp(suffix='_sort', prefix='fastq_', dir=os.path.dirname(mate1))
    cmd = f'fastq-sort --id --temporary-directory {tmp} -S 2G {mate1} > {output}'
    cmder.run(cmd, msg=f'Sort mate1 {mate1} {size(mate1)} ...', fmt_cmd=False)
    cmder.run(f'rm -r {tmp}')

    mate2 = mate1.replace('.mate1', '.mate2')
    if os.path.isfile(mate2):
        output2 = output.replace('.mate1', '.mate2')
        tmp = tempfile.mkdtemp(suffix='_sort', prefix='fastq_', dir=os.path.dirname(mate1))
        cmd = f'fastq-sort --id --temporary-directory {tmp} -S 2G {mate2} > {output2}'
        cmder.run(cmd, msg=f'Sort mate2 {mate2} {size(mate2)} ...', fmt_cmd=False)
        cmder.run(f'rm -r {tmp}')
    return output


@task(inputs=sort_mate, mkdir_before_run=[f'{ECLIP}/reference.genome.map'],
      outputs=lambda i: i.replace('repeat.elements.map', 'reference.genome.map')
      .replace('Unmapped.out.sort.mate1', 'Aligned.out.bam'))
def map_to_reference_genome(mate1, bam):
    # '--outSAMunmapped' flag needs to be set to 'Within', otherwise barcode_collapse.py for duplication removal will
    # throw out name not match error.
    prefix = os.path.dirname(bam)
    if not os.path.isdir(prefix):
        os.mkdir(prefix)
    mate2 = mate1.replace('.mate1', '.mate2')
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
    cmder.run(cmd, msg=message, pmt=True)
    return bam


def name_sort_bam(bam, out):
    cmder.run(f'samtools sort -n -@ {options.cores} -m 2G -o {out} {bam}',
              msg=f'Name sorting {bam} {size(bam)} ...')


def position_sort_bam(bam, out):
    cmder.run(f'samtools sort -@ {options.cores} -m 2G -o {out} {bam}',
              msg=f'Sorting {bam} {size(bam)} ...')


def index_sorted_bam(bam):
    cmder.run(f'samtools index -@ {options.cores} {bam}', msg=f'Indexing {bam} {size(bam)} ...')


def merge_paired_bam(bam, out):
    if not os.path.exists(out):
        key = out.replace(".merge.bam", "")
        barcodes = READS[os.path.basename(key)].barcodes
        if barcodes[0] == 'NIL':
            cmder.run(f'cp {bam} {out}')
        else:
            b1, b2 = barcodes
            if b1 in bam:
                b1, b2 = bam, bam.replace(b1, b2)
            else:
                b1, b2 = bam.replace(b2, b1), bam
            cmder.run(f'samtools merge -@ {options.cores} {out} {b1} {b2}',
                      msg=f'Merging {b1} {size(b1)} and {b2} {size(b2)} ...', pmt=True)


@task(inputs=map_to_reference_genome,
      outputs=lambda i: i.replace('/reference.genome.map', '').replace('/Aligned.out', '.sort'))
def prepare_bam(bam, out):
    if TYPE == 'single':
        name_sort = out.replace('.sort.bam', '.name.sort.bam')
        name_sort_bam(bam, name_sort)
        position_sort_bam(name_sort, out)
        index_sorted_bam(out)
        cmder.run(f'rm {name_sort}')
    else:
        name_sort_bam(bam, out)


def collapse_barcode(bam, out):
    logger.info(f'Deduplicating {bam} {size(bam)} by collapsing barcodes ...')
    verbosity = pysam.set_verbosity(0)
    with pysam.AlignmentFile(bam, 'rb') as b1, pysam.AlignmentFile(bam, 'rb') as b2:
        results = {}
        for read1, read2 in zip(itertools.islice(b1, 0, None, 2), itertools.islice(b2, 1, None, 2)):
            if read1.query_name != read2.query_name:
                raise ValueError(f'Read names do not match: {read1.query_name} != {read2.query_name}.')
            if read1.is_unmapped or read2.is_unmapped or read1.reference_name != read2.reference_name:
                continue
            if not read1.is_read1:
                read1, read2 = read2, read1
            randomer = read1.query_name.split(':')[0]
            start = read1.positions[-1] if read1.is_reverse else read1.pos
            stop = read2.positions[-1] if read2.is_reverse else read2.pos
            strand = '-' if read1.is_reverse else '+'
            location = (read1.reference_name, start, stop, strand, randomer)
            if location in results:
                continue
            results[location] = (read1, read2)
        with pysam.AlignmentFile(out, 'wb', template=b1) as o:
            for (read1, read2) in results.values():
                o.write(read1)
                o.write(read2)
        logger.info(f'Deduplicating {bam} {size(bam)} by collapsing barcodes complete.')
    pysam.set_verbosity(verbosity)


@task(inputs=prepare_bam, outputs=lambda i: i.replace('.sort.bam', '.sort.dedup.bam'), processes=options.cores)
def dedup_bam(bam, out):
    """Collapse barcodes of paired-end bam or umi_tools dedup single-end bam."""
    if TYPE == 'single':
        cmd = ['umi_tools', 'dedup', '--random-seed', 1, '--stdin', bam, '--method', 'unique', '--stdout', out]
        message = f'Deduplicating {bam} {size(bam)} by umi_tools dedup ...'
        cmder.run(cmd, msg=message, pmt=True)
    else:
        collapse_barcode(bam, out)


@task(inputs=dedup_bam, outputs=lambda i: i.replace('.sort.dedup.bam', '.sort.dedup.sort.bam'))
def sort_bam(bam, out):
    position_sort_bam(bam, out)


@task(inputs=[f'{ECLIP}/{k}.{v.barcodes[0]}.sort.dedup.sort.bam' for k, v in READS.items()],
      outputs=lambda i: f"{i.rsplit('.', maxsplit=5)[0]}.merge.bam", parent=sort_bam)
def merge_bam(bam, out):
    if TYPE == 'single':
        cmder.run(f'cp {bam} {out}')
    else:
        merge_paired_bam(bam, out)


@task(inputs=merge_bam, outputs=lambda i: i.replace('.merge.bam', '.bam'))
def index_bam(bam, out):
    if TYPE == 'paired':
        cmder.run(f'samtools view -f 128 -@ {options.cores} -b -o {out} {bam}',
                  msg=f'Creating bam {bam} {size(bam)} with r2 reads only ...')
    else:
        cmder.run(f'cp {bam} {out}')
    if not os.path.exists(f'{bam}.bai'):
        index_sorted_bam(out)


@task(inputs=index_bam, processes=options.cores, mkdir_before_run=[HUB],
      outputs=lambda i: f'{HUB}/{os.path.basename(i).replace(".bam", ".plus.bw")}')
def make_bigwig_files(bam, bigwig):
    def bam_to_bigwig(bam, scale, strand, bw):
        bg, bg_sort = bw.replace('.bw', '.bg'), bw.replace('.bw', '.sort.bg')
        cmd = f'genomeCoverageBed -ibam {bam} -bg -scale {scale} -strand {strand} -du -split > {bg}'
        cmder.run(cmd)
        cmd = f'bedSort {bg} {bg_sort}'
        cmder.run(cmd)
        cmd = f'bedGraphToBigWig {bg_sort} {options.genome}/chrNameLength.txt {bw}'
        cmder.run(cmd)
        cmder.run(f'rm {bg}')

    message, start_time = f'Make BigWig files for {bam} ...', time.perf_counter()
    logger.info(message)
    pos_bw, neg_bw = bigwig, bigwig.replace('.plus.bw', '.minus.bw')
    with pysam.AlignmentFile(bam, 'rb') as sam:
        total_reads = sam.mapped
    total_reads = total_reads / 2 if TYPE == 'paired' else total_reads
    try:
        scale = 1000000.0 / total_reads
    except ZeroDivisionError:
        logger.error(f'No reads was found in BAM {bam}, empty BigWig file was created.')
        with open(bigwig, 'w') as o:
            o.write('')
        return bigwig
    if TYPE == 'single':
        bam_to_bigwig(bam, scale, '+', pos_bw)
        bam_to_bigwig(bam, -1 * scale, '-', neg_bw)
    else:
        bam_to_bigwig(bam, -1 * scale, '-', pos_bw)
        bam_to_bigwig(bam, scale, '+', neg_bw)
    run_time = int(time.perf_counter() - start_time)
    message = message.replace(' ...', f' completed in [{str(datetime.timedelta(seconds=run_time))}].')
    logger.info(message)
    return bigwig


@task(inputs=make_bigwig_files, outputs='hub/hub.txt', kind='merge')
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
        for key in READS:
            plus = f'{key}.plus.bw'
            name1 = plus.replace('.bw', '').replace('.', '_')
            name2 = name1.replace('plus', 'minus')
            minus = f'{key}.minus.bw'
            o.write(block.format(track=track, name1=name1, name2=name2, basename=key, plus=plus, minus=minus))
    run_time = int(time.perf_counter() - start_time)
    message = message.replace(' ...', f' completed in [{str(datetime.timedelta(seconds=run_time))}].')
    logger.info(message)


def clipper_peaks(bam, bed=''):
    bed = bed if bed else bam.replace('.bam', '.peak.clusters.bed')
    if os.path.isfile(bed):
        logger.info(f'Clipper bed {bed} already exists.')
    else:
        cmd = f'clipper --species {options.species} --processors {options.cores} --bam {bam} --outfile {bed}'
        cmder.run(cmd, msg=f'Calling peaks from {bam} {size(bam)} using clipper ...', pmt=True)
    return bed


@task(inputs=[read.bam for read in READS.values() if read.type == 'IP'],
      outputs=lambda i: i.replace('.bam', '.peak.clusters.bed'), parent=index_bam)
def clipper(bam, bed):
    return clipper_peaks(bam, bed)


@task(inputs=[read.bam for read in READS.values() if read.type == 'IP'],
      outputs=lambda i: i.replace('.bam', '.crosslink.sites.bed'), parent=index_bam)
def pureclip(bam, bed):
    ip_bam, input_bam = [[sample.ip_read.bam, sample.input_read.bam] for sample in SAMPLES
                         if sample.ip_read.bam == bam][0]
    cmd = ['pureclip', '-i', ip_bam, '-bai', f'{ip_bam}.bai', '-g', f'{options.genome}/genome.fa',
           '-iv', "'chr1;chr2;chr3'", '-nt', options.cores, '-ibam', input_bam, '-ibai', f'{input_bam}.bai',
           '-o', bed, '-or', bed.replace('.crosslink.sites.bed', '.binding.regions.bed'),
           '>', bed.replace('.crosslink.sites.bed', '.pureclip.log')]
    cmder.run(cmd, msg=f'Calling peaks from {bam} {size(bam)} using pureCLIP ...', pmt=True)


def count_mapped_reads(bam):
    count = int(cmder.run(f'samtools view -c -F 0x4 {bam}', msg='')[0])
    logger.info(f'Found {count:,} mapped reads in {bam}.')
    return count


def calculate_entropy(bed, output, ip_read_count, input_read_count):
    logger.info(f'Calculating entropy for {bed} ...')
    columns = ['chrom', 'start', 'end', 'peak', 'ip_read_number', 'input_read_number',
               'p', 'v', 'method', 'status', 'l10p', 'l2fc',
               'ensg_overlap', 'feature_type', 'feature_ensg', 'gene', 'region']
    df = pd.read_csv(bed, sep='\t', header=None, names=columns)
    df = df[df.l2fc > 0]
    df['pi'] = df['ip_read_number'] / ip_read_count
    df['qi'] = df['input_read_number'] / input_read_count
    if df.empty:
        logger.error(f'No valid peaks found in {bed} (l2fc > 0 failed).')
        sys.exit(1)
    df['entropy'] = df.apply(lambda row: 0 if row.pi <= row.qi else row.pi * math.log2(row.pi / row.qi), axis=1)
    df['excess_reads'] = df['pi'] - df['qi']
    entropy = output.replace('.entropy.bed', '.entropy.full.bed')
    df.to_csv(entropy, index=False, columns=columns + ['entropy'], sep='\t', header=False)
    excess_read = output.replace('.bed', 'excess.reads.tsv')
    df.to_csv(excess_read, index=False, columns=columns + ['excess_reads'], sep='\t')
    df['strand'] = df.peak.str.split(':', expand=True)[2]
    df['l2fc'] = df['l2fc'].map('{:.15f}'.format)
    df['entropy'] = df['entropy'].map('{:.10f}'.format)
    # For IDR 2.0.2, columns 'excess_reads', 'pi', and 'qi' need to be excluded for .entropy.bed
    # For IDR 2.0.3, columns 'excess_reads', 'pi', and 'qi' need to be retained for .entropy.bed
    columns = ['chrom', 'start', 'end', 'l2fc', 'entropy', 'strand', 'excess_reads', 'pi', 'qi']
    df.to_csv(output, index=False, columns=columns, sep='\t', header=False)
    logger.info(f'Calculating entropy for {bed} complete.')
    return output


def peak(ip_bams, input_bams, peak_beds, reproducible_bed):
    entropy_beds, folder = {}, os.path.dirname(reproducible_bed)
    for ip_bam, input_bam, peak_bed in zip(ip_bams, input_bams, peak_beds):
        logger.info(f'Processing {peak_bed} ...')
        name = os.path.basename(peak_bed.replace('.peak.clusters.bed', ''))
        normalized_bed = peak_bed.replace('.peak.clusters.bed', '.peak.clusters.normalized.bed')
        ip_read_count, input_read_count = count_mapped_reads(ip_bam), count_mapped_reads(input_bam)
        if os.path.isfile(normalized_bed):
            logger.info(f'Normalized peaks {normalized_bed} already exist.')
        else:
            cmd = ['overlap_peak.pl', ip_bam, input_bam, peak_bed, ip_read_count, input_read_count, normalized_bed]
            cmder.run(cmd, msg=f'Normalizing peaks in {peak_bed} ...', pmt=True)

        compressed_bed = normalized_bed.replace('.normalized.bed', '.normalized.compressed.bed')
        if os.path.isfile(compressed_bed):
            logger.info(f'Compressed peaks {compressed_bed} already exist.')
        else:
            cmd = ['compress_peak.pl', normalized_bed.replace('.bed', '.full.bed'), compressed_bed]
            cmder.run(cmd, msg=f'Compressing peaks in {normalized_bed} ...', pmt=True)

        annotated_bed = compressed_bed.replace('.compressed.bed', '.compressed.annotated.full.bed')
        if os.path.isfile(annotated_bed):
            logger.info(f'Annotated peaks {annotated_bed} already exist.')
        else:
            species = 'hg19' if options.species == 'hg19chr19' else options.species
            cmd = ['annotate_peak.pl', compressed_bed.replace('.bed', '.full.bed'), annotated_bed, species, 'full']
            cmder.run(cmd, msg=f'Annotating peaks in {compressed_bed} ...', pmt=True)

        entropy_bed = annotated_bed.replace('.annotated.full.bed', '.annotated.entropy.bed')
        annotated_bed = compressed_bed.replace('.compressed.bed', '.compressed.annotated.full.bed')
        if os.path.isfile(entropy_bed):
            logger.info(f'Entropy bed {annotated_bed} already exist.')
        else:
            calculate_entropy(annotated_bed, entropy_bed, ip_read_count, input_read_count)
        entropy_beds[name] = entropy_bed

    for key1, key2 in itertools.combinations(entropy_beds.keys(), 2):
        idr_out, idr_bed = f'{folder}/{key1}.vs.{key2}.idr.out', f'{folder}/{key1}.vs.{key2}.idr.out.bed'
        entropy_bed1, entropy_bed2 = entropy_beds[key1], entropy_beds[key2]
        if os.path.isfile(idr_out):
            logger.info(f'IDR out {idr_out} already exist.')
        else:
            cmd = ['idr', '--sample', entropy_bed1, entropy_bed2, '--input-file-type', 'bed', '--rank', '5',
                   '--peak-merge-method', 'max', '--plot', '-o', idr_out]
            cmder.run(cmd, msg=f'Running IDR to rank peaks in {entropy_bed1} and\n{" " * 40}{entropy_bed2} ...',
                      pmt=True)
        if os.path.isfile(idr_bed):
            logger.info(f'IDR bed {idr_bed} already exist.')
        else:
            cmd = ['parse_idr_peaks.pl', idr_out,
                   entropy_bed1.replace('.bed', '.full.bed'), entropy_bed2.replace('.bed', '.full.bed'), idr_bed]
            cmder.run(cmd, msg=f'Parsing IDR peaks in {idr_out} ...', pmt=True)

    idr_bed = f'{folder}/{".vs.".join(entropy_beds.keys())}.idr.out.bed'
    if len(entropy_beds) == 2:
        script = 'reproducible_peaks.pl'
    elif len(entropy_beds) == 3:
        script = 'reproducible_peaks_3.pl'
        if os.path.isfile(idr_bed):
            logger.info(f'IDR bed {idr_bed} already exist.')
        else:
            bed1, bed2, bed3 = [f'{folder}/{key1}.vs.{key2}.idr.out.bed'
                                for key1, key2 in itertools.combinations(entropy_beds.keys(), 2)]
            cmder.run(f'bedtools intersect -a {bed1} -b {bed2} {bed3} > {idr_bed}', msg='Intersecting IDR beds ...')
    else:
        raise ValueError('Method for handling more than 3 replicates has not been implemented yet.')

    if os.path.isfile(reproducible_bed):
        logger.info(f'Reproducible peaks {reproducible_bed} already exist.')
    else:
        custom_bed = reproducible_bed.replace('.peaks.bed', '.peaks.custom.bed')
        idr_normalized_full_beds, entropy_full_beds, reproducible_full_beds = [], [], []
        for ip_bam, input_bam, peak_bed in zip(ip_bams, input_bams, peak_beds):
            name = os.path.basename(peak_bed.replace('.peak.clusters.bed', ''))
            idr_normalized_bed = f'{folder}/{name}.idr.normalized.bed'
            if os.path.isfile(idr_normalized_bed):
                logger.info(f'IDR normalized bed {idr_normalized_bed} already exist.')
            else:
                cmd = ['overlap_peak.pl', ip_bam, input_bam, idr_bed,
                       count_mapped_reads(ip_bam), count_mapped_reads(input_bam), idr_normalized_bed]
                cmder.run(cmd, msg=f'Normalizing IDR peaks for sample {name} ...', pmt=True)
            idr_normalized_full_beds.append(idr_normalized_bed.replace('.bed', '.full.bed'))
            entropy_full_beds.append(f'{folder}/{name}.peak.clusters.normalized.compressed.annotated.entropy.full.bed')
            reproducible_full_beds.append(f'{folder}/{name}.reproducible.peaks.full.bed')
    
        cmd = [script, ] + idr_normalized_full_beds + reproducible_full_beds
        cmd += [reproducible_bed, custom_bed] + entropy_full_beds
        cmd += [f'{folder}/{".vs.".join(entropy_beds.keys())}.idr.out']
        cmder.run(cmd, msg='Identifying reproducible peaks ...', pmt=True)
    return reproducible_bed


@task(inputs=[], outputs=f'{ECLIP}/{".vs.".join([s.key for s in SAMPLES])}.reproducible.peaks.bed', parent=clipper)
def reproducible_peaks(inputs, outputs):
    ip_bams, input_bams, peak_beds = [], [], []
    for sample in SAMPLES:
        ip_bams.append(sample.ip_read.bam)
        input_bams.append(sample.input_read.bam)
        peak_beds.append(sample.bed)
    peak(ip_bams, input_bams, peak_beds,  outputs)


def split_bam(bam, bam1, bam2):
    if os.path.isfile(bam1) and os.path.isfile(bam2):
        logger.info(f'BAMs {bam1} and {bam2} already exist.')
    else:
        half_lines = int(count_mapped_reads(bam) / 2) + 1
        cmd = f'samtools view {bam} | shuf | split -d -l {half_lines} - {bam}'
        cmder.run(cmd, msg=f'Shuffling and splitting {bam} ...')
        tmp_bam1, tmp_bam2 = bam1.replace('.bam', '.tmp.bam'), bam2.replace('.bam', '.tmp.bam')
        cmd = f'samtools view -H {bam} | cat - {bam}00 | samtools view -bS - > {tmp_bam1}'
        
        cmder.run(cmd, msg=f'Creating headers for {bam1} ...')
        cmder.run(f'samtools sort -@ {options.cores} -m 2G -o {bam1} {tmp_bam1}')
        cmd = f'samtools view -H {bam} | cat - {bam}01 | samtools view -bS - > {tmp_bam2}'
        
        cmder.run(cmd, msg=f'Creating headers for {bam2} ...')
        cmder.run(f'samtools sort -@ {options.cores} -m 2G -o {bam2} {tmp_bam2}')
        cmder.run(f'rm {bam}00 {bam}01 {tmp_bam1} {tmp_bam2}')
    return bam1, bam2


def count_lines(file):
    lines = int(cmder.run(f'wc -l {file}')[0].split()[0])
    logger.info(f'Found {lines:,} lines in {file}.')
    return lines


@task(inputs=[], outputs=f'{ECLIP}/rescue.ratio.txt', parent=clipper, mkdir_before_run=['rescue'])
def rescue_ratio(inputs, outputs):
    def prepare_pseudo_bam(bam1, bam2, basename):
        pseudo_bam = f'{basename}.bam'
        tmp_pseudo_bam = pseudo_bam.replace('.bam', '.tmp.bam')
        cmd = f'samtools merge {tmp_pseudo_bam} {bam1} {bam2}'
        cmder.run(cmd, msg=f'Merging {bam1} and {bam2} ...')

        cmder.run(f'samtools sort -@ {options.cores} -m 2G -o {pseudo_bam} {tmp_pseudo_bam}')
        cmder.run(f'rm {tmp_pseudo_bam}')

        bam1, bam2 = split_bam(pseudo_bam, f'{basename}.pseudo.01.bam', f'{basename}.pseudo.02.bam')
        return bam1, bam2

    pseudo_ip_bams, pseudo_input_bams, pseudo_peak_beds = [], [], []
    for i, (sample1, sample2) in enumerate(itertools.combinations(SAMPLES, 2), start=1):
        pseudo_ip_bam = prepare_pseudo_bam(sample1.ip_read.bam, sample2.ip_read.bam,
                                           f'rescue/{sample1.ip_read.name}.{sample2.ip_read.name}')
        pseudo_ip_bams.extend(pseudo_ip_bam)

        pseudo_input_bam = prepare_pseudo_bam(sample1.input_read.bam, sample2.input_read.bam,
                                              f'rescue/{sample1.input_read.name}.{sample2.input_read.name}')
        pseudo_input_bams.extend(pseudo_input_bam)
        
        pseudo_peak_beds.extend([clipper_peaks(bam) for bam in pseudo_ip_bam])

    key = ".vs.".join([sample.ip_read.name for sample in SAMPLES])
    pseudo_reproducible_bed = f'rescue/{key}.reproducible.peaks.bed'
    peak(pseudo_ip_bams, pseudo_input_bams, pseudo_peak_beds, pseudo_reproducible_bed)
    pseudo_count = count_lines(pseudo_reproducible_bed)
    
    count = count_lines(f'{ECLIP}/{key}.reproducible.peaks.bed')
    try:
        ratio = max(count, pseudo_count) / min(count, pseudo_count)
    except ZeroDivisionError:
        ratio = 0
        logger.error(f'No peaks found in reproducible peaks or pseudo reproducible peaks, return ratio 0.')
    with open(outputs, 'w') as o1, open(outputs.replace(f'{ECLIP}/', 'rescue/'), 'w') as o2:
        o1.write(f'{ratio}\n')
        o2.write(f'{ratio}\n')


@task(inputs=[], outputs=f'{ECLIP}/consistency.ratio.txt', parent=clipper, mkdir_before_run=['consistency'])
def consistency_ratio(inputs, outputs):
    counts = []
    for sample in SAMPLES:
        split_ip_bams = split_bam(sample.ip_read.bam,
                                  f'consistency/{sample.ip_read.name}.split.01.bam',
                                  f'consistency/{sample.ip_read.name}.split.02.bam')
        split_input_bams = split_bam(sample.input_read.bam,
                                     f'consistency/{sample.input_read.name}.split.01.bam',
                                     f'consistency/{sample.input_read.name}.split.02.bam')
        split_peak_beds = [clipper_peaks(split_ip_bams[0]), clipper_peaks(split_ip_bams[1])]

        reproducible_bed = f'consistency/{sample.key}.split.reproducible.peaks.bed'
        peak(split_ip_bams, split_input_bams, split_peak_beds, reproducible_bed)
        counts.append(count_lines(reproducible_bed))

    try:
        ratio = counts[0] / counts[1]
    except ZeroDivisionError:
        ratio = 0
        logger.error(f'No peaks found in one of the split reproducible peaks, return ratio 0.')
    with open(outputs, 'w') as o1, open(outputs.replace(f'{ECLIP}/', 'consistency/'), 'w') as o2:
        o1.write(f'{ratio}\n')
        o2.write(f'{ratio}\n')


def prepare_fastqs():
    fastqs = []
    for key, read in READS.items():
        fastq1, fastq2 = read.link1, read.link2
        fastqs.append(fastq1)
        if fastq2:
            fastqs.append(fastq2)
            barcodes = read.barcodes
            for barcode in barcodes:
                fastqs.extend([f'{ECLIP}/{key}.{barcode}.r1.fastq.gz', f'{ECLIP}/{key}.{barcode}.r2.fastq.gz'])
    return fastqs


@task(inputs=prepare_fastqs(), outputs=lambda i: i.replace(f'{ECLIP}/', f'{QC}/').replace('.fastq.gz', '.fastq.qc.txt'),
      parent=consistency_ratio, mkdir_before_run=['qc'], processes=options.cores)
def falco(fastq, txt):
    tmp = tempfile.mkdtemp(suffix='_qc', prefix='falco_', dir=QC)
    cmd = f'falco --outdir {tmp} --skip-html {fastq}'
    try:
        cmder.run(cmd, msg=f'Checking reads in {fastq} {size(fastq)} using falco ...', fmt_cmd=False)
        cmder.run(f'mv {tmp}/fastqc_data.txt {txt}')
    finally:
        shutil.rmtree(tmp)


@task(inputs=[], outputs='summary.html', kind='create', parent=falco)
def summary(inputs, outputs):
    def parse_falco_metrics(txt):
        sample, reads, duplicate = '', 0, 0.0
        with open(txt) as f:
            for line in f:
                if line.startswith('Filename'):
                    sample = line.strip().split('\t')[1].replace('.fastq.gz', '')
                elif line.startswith('Total Sequences'):
                    reads = int(line.strip().split('\t')[1])
                elif line.startswith('#Total Deduplicated Percentage'):
                    duplicate = float(line.strip().split('\t')[1])
                    break
        return sample, reads, duplicate

    def parse_cutadapt_metrics(metrics):
        reads, reads_too_short = 0, 0
        with open(metrics) as f:
            for line in f:
                if ('Reads written (passing filters)' in line) or ('Pairs written (passing filters)' in line):
                    reads = int(line.strip().split()[-2].replace(',', ''))
                elif ('Reads that were too short' in line) or ('Pairs that were too short' in line):
                    reads_too_short = int(line.strip().split()[-2].replace(',', ''))
        return os.path.basename(metrics).split('.trim.')[0], reads, reads_too_short

    def parse_star_log(log):
        counts = [os.path.basename(os.path.dirname(log))]
        with open(log) as f:
            for line in f:
                if 'Number of input reads' in line:
                    reads = int(line.strip().split('\t')[-1])
                elif 'Uniquely mapped reads number' in line:
                    counts.append(int(line.strip().split('\t')[-1]))
                elif 'Number of reads mapped to multiple loci' in line:
                    counts.append(int(line.strip().split('\t')[-1]))
                elif 'Number of reads mapped to too many loci' in line:
                    counts.append(int(line.strip().split('\t')[-1]))
                elif '% of reads unmapped: too many mismatches' in line:
                    counts.append(int(float(line.strip().split('\t')[-1].replace('%', '')) * reads / 100))
                elif '% of reads unmapped: too short' in line:
                    counts.append(int(float(line.strip().split('\t')[-1].replace('%', '')) * reads / 100))
                elif '% of reads unmapped: other' in line:
                    counts.append(int(float(line.strip().split('\t')[-1].replace('%', '')) * reads / 100))
        return counts

    def get_dedup_metrics(bam1, bam2):
        n, _ = cmder.run(f'samtools view -c -F 0x4 {bam1}', msg='')
        n2, _ = cmder.run(f'samtools view -c -F 0x4 {bam2}', msg='')
        return os.path.basename(bam1).replace('.sort.bam', ''), int(n) - int(n2), int(n2)

    def get_usable_reads(bam):
        return os.path.basename(bam1).replace('.bam', ''), int(cmder.run(f'samtools view -c -F 0x4 {bam}', msg='')[0])

    raw_reads, demux_reads, cut1, cut2 = defaultdict(list), defaultdict(list), defaultdict(list), defaultdict(list)
    repeat_map, genome_map = defaultdict(list), defaultdict(list)
    dedup_reads, usable_reads = defaultdict(list), defaultdict(list)
    for key, read in READS.items():
        sample, unique_reads, duplicate_reads = parse_falco_metrics(f'{QC}/{key}.r1.fastq.qc.txt')
        raw_reads['Sample'].append(sample)
        raw_reads['Unique Reads'].append(unique_reads)
        raw_reads['Duplicate Reads'].append(duplicate_reads)
        barcodes = read.barcodes if read.link2 else ['umi']
        if read.link2:
            sample, unique_reads, duplicate_reads = parse_falco_metrics(f'{QC}/{key}.r2.fastq.qc.txt')
            raw_reads['Sample'].append(sample)
            raw_reads['Unique Reads'].append(unique_reads)
            raw_reads['Duplicate Reads'].append(duplicate_reads)
            for barcode in set(read.barcodes):
                sample, unique_reads, duplicate_reads = parse_falco_metrics(f'{QC}/{key}.{barcode}.r1.fastq.qc.txt')
                demux_reads['Sample'].append(sample)
                demux_reads['Unique Reads'].append(unique_reads)
                demux_reads['Duplicate Reads'].append(duplicate_reads)
                sample, unique_reads, duplicate_reads = parse_falco_metrics(f'{QC}/{key}.{barcode}.r2.fastq.qc.txt')
                demux_reads['Sample'].append(sample)
                demux_reads['Unique Reads'].append(unique_reads)
                demux_reads['Duplicate Reads'].append(duplicate_reads)

        for barcode in set(barcodes):
            cut1_metrics = f'{ECLIP}/{key}.{barcode}.trim.first.metrics'
            sample, passed_reads, failed_reads = parse_cutadapt_metrics(cut1_metrics)
            cut1['Sample'].append(sample)
            cut1['Passed Reads'].append(passed_reads)
            cut1['Failed Reads'].append(failed_reads)
            # shutil.move(cut1_metrics, cut1_metrics.replace(f'{ECLIP}/', f'{QC}/'))

            cut2_metrics = f'{ECLIP}/{key}.{barcode}.trim.second.metrics'
            sample, passed_reads, failed_reads = parse_cutadapt_metrics(cut2_metrics)
            cut2['Sample'].append(sample)
            cut2['Passed Reads'].append(passed_reads)
            cut2['Failed Reads'].append(failed_reads)
            # shutil.move(cut2_metrics, cut2_metrics.replace(f'{ECLIP}/', f'{QC}/'))

            names = ['Sample', 'Uniquely mapped', 'Mapped to multiple loci', 'Mapped to too many loci',
                     'Unmapped: too many mismatches', 'Unmapped: too short', 'Unmapped: other']
            repeat_map_log = f'{ECLIP}/repeat.elements.map/{key}.{barcode}/Log.final.out'
            counts = parse_star_log(repeat_map_log)
            for i, name in enumerate(names):
                repeat_map[name].append(counts[i])
            # shutil.move(repeat_map_log, f'{QC}/{key}.{barcode}.repeat.element.map.log')

            genome_map_log = f'{ECLIP}/reference.genome.map/{key}.{barcode}/Log.final.out'
            counts = parse_star_log(genome_map_log)
            for i, name in enumerate(names):
                genome_map[name].append(counts[i])
            # shutil.move(genome_map_log, f'{QC}/{key}.{barcode}.reference.genome.map.log')

            bam1, bam2 = f'{ECLIP}/{key}.{barcode}.sort.bam', f'{ECLIP}/{key}.{barcode}.sort.dedup.sort.bam'
            sample, unique_reads, duplicate_reads = get_dedup_metrics(bam1, bam2)
            dedup_reads['Sample'].append(sample)
            dedup_reads['Unique Reads'].append(unique_reads)
            dedup_reads['Duplicate Reads'].append(duplicate_reads)

            counts = get_usable_reads(f'{ECLIP}/{key}.bam')
            usable_reads['Sample'] = counts[0]
            usable_reads['Usable Reads'] = counts[1]

    counts = {'raw_reads': raw_reads, 'demux_reads': demux_reads, 'cut1': cut1, 'cut2': cut2,
              'repeat_map': repeat_map, 'genome_map': genome_map,
              'dedup_reads': dedup_reads, 'usable_reads': usable_reads}
    with open(f'{QC}/summary.json', 'w') as o:
        json.dump(counts, o, indent=4)


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
        
        
def dataset_basic():
    data = [['Data Type', 'Sample', 'Description', 'Replicate', 'Library Type', 'Barcodes']]
    library_type = 'PE' if TYPE == 'paired' else 'SE'
    for i, sample in enumerate(SAMPLES, 1):
        data.append(['eCLIP', sample.ip_read.key, 'IP', i, library_type, ', '.join(sample.ip_read.barcodes)])
        data.append(['eCLIP', sample.input_read.key, 'INPUT', i, library_type, ', '.join(sample.input_read.barcodes)])
    return data


# from bokeh.io import output_file, show
# from bokeh.models import ColumnDataSource, FactorRange, NumeralTickFormatter, Legend
# from bokeh.plotting import figure
# from bokeh.models.widgets import Tabs, Panel
#
# import datapane as dp


def bar_plot(data, data_type, tooltips, colors):
    samples, categories = data['Samples'], [key for key in data.keys if key != 'Samples']
    fig = figure(y_range=FactorRange(factors=samples[::-1]), plot_height=250, toolbar_location=None, tooltips=tooltips)
    fig.add_layout(Legend(), 'right')
    fig.hbar_stack(categories, y='Samples', height=0.8, color=colors, legend_label=categories,
                   source=ColumnDataSource(data))
    fig.x_range.start = 0
    fig.x_range.range_padding = 0.1
    formatter = NumeralTickFormatter(format="0 a") if data_type == 'counts' else NumeralTickFormatter(format="0,0")
    fig.xaxis[0].formatter = formatter
    fig.xaxis.axis_label = 'Number of reads' if data_type == 'counts' else 'Percent of reads (%)'
    fig.xaxis.axis_line_color = None
    fig.y_range.range_padding = 0.1
    fig.yaxis.axis_line_color = None
    fig.ygrid.grid_line_color = None
    fig.legend.border_line_color = None
    fig.axis.minor_tick_line_color = None
    fig.axis.major_tick_line_color = None
    fig.outline_line_color = None
    return fig


def count_percent_plot(counts, tooltips=None, colors=("#718dbf", "#e84d60")):
    if tooltips:
        count_tooltips, percent_tooltips = tooltips
    else:
        count_tooltips = [("Sample", "@samples"), ("Unique Reads", "@{Unique Reads}{0.00 a}"),
                          ("Duplicate Reads", "@{Duplicate Reads}{0.00 a}")]
        percent_tooltips = [("Sample", "@samples"), ("Unique Reads (%)", "@{Unique Reads}{0.00}"),
                            ("Duplicate Reads (%)", "@{Duplicate Reads}{0.00}")]
    count_figure = bar_plot(counts, 'counts', count_tooltips, colors)
    df = pd.DataFrame(counts)
    df = df.set_index('Samples')
    df = df.div(df.sum(axis=1), axis=0)
    df.reset_index()
    counts = {column: df[column].tolist() for column in df.columns}
    percent_figure = bar_plot(counts, 'percent', percent_tooltips, colors)

    count_panel = Panel(child=count_figure, title='Count')
    percent_panel = Panel(child=percent_figure, title='Percentage')

    tabs = Tabs(tabs=[count_panel, percent_panel])
    
    return tabs
    

# def plot_raw_reads(counts):
#     samples, categories, uniques, duplicates = [], ['Unique Reads', 'Duplicate Reads'], [], []
#     unique_percents, duplicate_percents = [], []
#     for sample, count, percent in counts:
#         samples.append(sample)
#         duplicate = int(count * percent / 100)
#         duplicates.append(duplicate)
#         uniques.append(count - duplicate)
#         unique_percents.append(100 - percent)
#         duplicate_percents.append(percent)
#     counts = {'samples': samples, 'Unique Reads': uniques, 'Duplicate Reads': duplicates}
#     percents = {'samples': samples, 'Unique Reads': unique_percents, 'Duplicate Reads': duplicate_percents}
#
#     # output_file("/Users/fei/Downloads/work/software/tmp/stacked_split.html")
#
#     tooltips = [("Sample", "@samples"), ("Unique Reads", "@{Unique Reads}{0.00 a}"),
#                 ("Duplicate Reads", "@{Duplicate Reads}{0.00 a}")]
#     count_fig = figure(y_range=FactorRange(factors=samples[::-1]), plot_height=250,
#                        toolbar_location=None, tooltips=tooltips)
#
#     count_fig.add_layout(Legend(), 'right')
#     count_fig.hbar_stack(categories, y='samples', height=0.8, color=["#718dbf", "#e84d60"], legend_label=categories,
#                          source=ColumnDataSource(counts))
#
#     count_fig.x_range.start = 0
#     count_fig.x_range.range_padding = 0.1
#     count_fig.xaxis[0].formatter = NumeralTickFormatter(format="0 a")
#     count_fig.xaxis.axis_label = 'Number of reads'
#     count_fig.xaxis.axis_line_color = None
#     count_fig.y_range.range_padding = 0.1
#     count_fig.yaxis.axis_line_color = None
#     count_fig.ygrid.grid_line_color = None
#     count_fig.legend.border_line_color = None
#     count_fig.axis.minor_tick_line_color = None
#     count_fig.axis.major_tick_line_color = None
#     count_fig.outline_line_color = None
#
#     tooltips = [("Sample", "@samples"), ("Unique Reads (%)", "@{Unique Reads}{0.00}"),
#                 ("Duplicate Reads (%)", "@{Duplicate Reads}{0.00}")]
#     percent_fig = figure(y_range=FactorRange(factors=samples[::-1]), plot_height=250, toolbar_location=None,
#                          tooltips=tooltips)
#     percent_fig.add_layout(Legend(), 'right')
#     percent_fig.hbar_stack(categories, y='samples', height=0.8, color=["#718dbf", "#e84d60"],
#                            legend_label=categories, source=ColumnDataSource(percents))
#     percent_fig.x_range.start = 0
#     percent_fig.x_range.range_padding = 0.1
#     percent_fig.xaxis[0].formatter = NumeralTickFormatter(format="0,0")
#     percent_fig.xaxis.axis_label = 'Percent of reads (%)'
#     percent_fig.xaxis.axis_line_color = None
#     percent_fig.y_range.range_padding = 0.1
#     percent_fig.yaxis.axis_line_color = None
#     percent_fig.ygrid.grid_line_color = None
#     percent_fig.legend.border_line_color = None
#     percent_fig.axis.minor_tick_line_color = None
#     percent_fig.axis.major_tick_line_color = None
#     percent_fig.outline_line_color = None
#
#     count_panel = Panel(child=count_fig, title='Count')
#     percent_panel = Panel(child=percent_fig, title='Percentage')
#
#     tabs = Tabs(tabs=[count_panel, percent_panel])
#
#     # dp.Report("## Vaccination Report",
#     #           dp.Plot(tabs, caption="Vaccinations by manufacturer over time"),
#     #           ).save(path='/Users/fei/Downloads/work/software/tmp/report.html', open=True)
#
#     # show(tabs)
#     return tabs


if __name__ == '__main__':
    pass
    collapse_barcode('/storage/vannostrand/analysis/20210528_Deniz/eclip/K002000195.80891.X1A.r1.fq.genome-mappedSo.bam',
                     '/storage/vannostrand/software/eclip/test/fastd/dedup.bam')