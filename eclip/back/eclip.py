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


def fastq_file(s):
    s = s if ',' in s else f'{s},'
    if not len(s) == 2:
        logger.error(f'Invalid FASTQ file(s) specified, {len(s)} files were given while only accepts 1 or 2 files.')
        sys.exit(1)
    return [file_path(p) if p else p for p in s.split(',')]


parser = argparse.ArgumentParser(description=__doc__, prog='eclip')
parser.add_argument('--input_fastqs', nargs='+', required=True, type=fastq_file, help='Path(s) to INPUT FASTQ file(s).')
parser.add_argument('--ip_fastqs', nargs='+', required=True, type=fastq_file, help='Path(s) to IP FASTQ file(s).')
parser.add_argument('--names', nargs='+', required=True, help='Shortnames for each sample, e.g., rep1, rep2.')
parser.add_argument('--barcodes_pattern', help="Pattern of barcodes for umi-tools using to extract UMIs (for "
                                               "single-end dataset only, default: NNNNNNNNNN", default='NNNNNNNNNN')
parser.add_argument('--adapters_fasta', help="Path to the fasta file contains adapters and their sequences (for "
                                             "single-end dataset only.", type=file_path)
parser.add_argument('--barcodes_fasta', help="Path to the fasta file contains barcodes and their sequences (for "
                                             "paired-end dataset only.", type=file_path)
parser.add_argument('--randomer_length', type=int, help="Length (int) of the randomer, default: 10", default=10)
parser.add_argument('--allow_mismatch', type=int,
                    help="Allowed mismatch among barcodes for demultiplexing, default: 1", default=1)
parser.add_argument('--dataset', help='Name of the eCLIP dataset, default: eCLIP.', default='eCLIP')
parser.add_argument('--outdir', help="Path to the output directory, an pre-existing directory or current directory.")
parser.add_argument('--genome', help="Path to STAR reference genome index directory.", type=dir_path)
parser.add_argument('--repeat', help="Path to STAR repeat elements index directory.", type=dir_path)
parser.add_argument('--species', help="Species name (short name code) the dataset associated with, e.g., hg19, mm10.")
parser.add_argument('--blacklist_bed', help="Path to the bed file contains blacklist.")
parser.add_argument('--track', help="Name for the UCSC Genome Browser track, default: eCLIP", default='eCLIP')
parser.add_argument('--track_genome', help="Genome name (a short name code) for the track.", default='hg19')
parser.add_argument('--l2fc', type=int, help="Only consider peaks at or above this log2 fold change cutoff.", default=3)
parser.add_argument('--l10p', type=int, help="Only consider peaks at or above this log10 p value cutoff.", default=3)
parser.add_argument('--enrichment_filter', type=int, help="Pre-filter peaks that are enriched over input.", default=0)

parser.add_argument('--job', help="Name of your job", default='eCLIP')
parser.add_argument('--email', help='Email address for notifying you the start, end, and abort of you job.')
parser.add_argument('--time', type=int, help='Time (in integer hours) for running your job.', default=36)
parser.add_argument('--memory', type=int, help='Amount of memory (in GB) for all CPU cores.', default=32)
parser.add_argument('--cpus', type=int, help='Maximum number of CPU cores can be used for your job.', default=16)
parser.add_argument('--scheduler', choices=('pbs', 'qsub', 'slurm', 'sbatch'),
                    help='Name (case insensitive) of the scheduler on your cluster.')
parser.add_argument('--hold_submit', action='store_true',
                    help='Generate the submit script but hold it without submitting to the job scheduler.')

parser.add_argument('--debug', action='store_true', help='Invoke debug mode (for development use only).')
parser.add_argument('--dry_run', action='store_true',
                    help='Print out steps and files involved in each step without actually running the pipeline.')

START_TIME = time.perf_counter()


options = parser.parse_args()
setattr(options, 'outdir', options.outdir or os.getcwd())
dir_path(options.outdir)
os.chdir(options.outdir)

ips, inputs, names = options.ip_fastqs, options.input_fastqs, options.labels
if len(ips) == len(names):
    if len(ips) == len(inputs):
        input_type = 'single-input'
    else:
        input_type = 'multiple-inputs'
        if len(inputs) != 1:
            logger.error('Wrong number of input_fastqs were provided.')
            sys.exit(1)
else:
    logger.error('Number of items in ip_fastqs and names are not equal.')
    sys.exit(1)


class Read:
    def __init__(self, fastq1, fastq2, read_name, read_type):
        self.fastq1 = fastq1
        self.fastq2 = fastq2
        self.read_name = read_name
        self.read_type = read_type
        self.key = read_name if read_name else fastq1.replace('.fastq.gz', '').replace('.fq.gz', '')
        self.paired = True if self.fastq2 else False
        self.link1 = f'{self.key}.r1.fastq.gz'
        self.link2 = f'{self.key}.r2.fastq.gz' if self.fastq2 else ''


class Sample:
    def __init__(self, ip_read, input_read):
        self.key = self.ip_read.key
        self.ip_read = ip_read
        self.input_read = input_read
        self.bed = self.ip_read.bed


READS, SAMPLES, NEED_TO_REMOVE = {}, {}, []
inputs = inputs if len(inputs) == len(ips) else inputs * len(ips)
for (ip_fastq1, ip_fastq2), (input_fastq1, input_fastq2), name in zip(ips, inputs, names):
    ip_read = READ(ip_fastq1, ip_fastq2, name, 'IP')
    input_read = READ(input_fastq1, input_fastq2, name, 'INPUT')
    READS[ip_fastq1] = ip_read
    READS[input_fastq1] = input_read
    SAMPLES[name] = (ip_read, input_read)

HASH = 100000
VERSION = 1.0
HEADER = fr"""
                           ____   _       ___   ____
                    ___   / ___| | |     |_ _| |  _ \
                   / _ \ | |     | |      | |  | |_) |
                  |  __/ | |___  | |___   | |  |  __/
                   \___|  \____| |_____| |___| |_|

                               VERSION {VERSION}
"""


@task(inputs=[read.fastq1 for read in READS.values()], outputs=lambda i: READS[i].link1)
def soft_link(fastq1, link1):
    """ Create soft links for original fastq files. """

    def make_link(path, link):
        if path == os.path.abspath(link):
            logger.warning("No symbolic link was made for {path}! You are directly working on the original file!")
        else:
            os.symlink(path, link)

    make_link(fastq1, link1)
    NEED_TO_REMOVE.append(link1)
    read = READS[fastq1]
    if read2.link2:
        make_link(read.fastq2, link2)
        NEED_TO_REMOVE.append(link2)
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
    max_barcode_length = max(len(barcode) for barcode in barcodes_dict.values())
    writers = {}
    for barcode in set(barcodes):
        file1, file2 = f'{basename}.{barcode}.r1.fastq.gz', f'{basename}.{barcode}.r2.fastq.gz'
        writers[barcode] = (gzip.open(file1, 'wt'), gzip.open(file2, 'wt'))
        NEED_TO_REMOVE.extend([file1, file2])

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
        NEED_TO_REMOVE.append(fastq1.replace('.r1.fastq.gz', '.umi.r1.fastq.gz'))


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
            NEED_TO_REMOVE.extend([tmp, out1])
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
            NEED_TO_REMOVE.extend([tmp1, tmp2, out1, out2])
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
        cmder.run(cmd, msg=f'Sorting {fastq} {size(fastq)} ...', pmt=True)
        NEED_TO_ZIP.append(output)
    finally:
        shutil.rmtree(tmp)

    fastq2 = fastq.replace('.r1.', '.r2.')
    if os.path.isfile(fastq2):
        output2 = output.replace('.r1.', '.r2.')
        tmp = tempfile.mkdtemp(suffix='_sort', prefix='fastq_', dir=ECLIP)
        cmd = f'fastq-sort --id --temporary-directory {tmp} -S 2G {fastq2} >  {output2}'
        try:
            cmder.run(cmd, msg=f'Sort fastq {fastq2} {size(fastq2)} ...')
            NEED_TO_ZIP.append(output2)
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
    cmder.run(cmd, msg=f'Sort mate1 {mate1} {size(mate1)} ...')
    cmder.run(f'rm -r {tmp}')

    mate2 = mate1.replace('.mate1', '.mate2')
    if os.path.isfile(mate2):
        output2 = output.replace('.mate1', '.mate2')
        tmp = tempfile.mkdtemp(suffix='_sort', prefix='fastq_', dir=os.path.dirname(mate1))
        cmd = f'fastq-sort --id --temporary-directory {tmp} -S 2G {mate2} > {output2}'
        cmder.run(cmd, msg=f'Sort mate2 {mate2} {size(mate2)} ...')
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
        cmder.run(f'rm {bg} {bg_sort}')

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
    color 100,0,0
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
    # '-iv', "'chr1;chr2;chr3'", Genomic chromosomes to learn HMM parameters
    cmd = ['pureclip', '-i', ip_bam, '-bai', f'{ip_bam}.bai', '-g', f'{options.genome}/genome.fa',
           '-nt', options.cores, '-ibam', input_bam, '-ibai', f'{input_bam}.bai',
           '-o', bed, '-or', bed.replace('.crosslink.sites.bed', '.binding.regions.bed'),
           '>', bed.replace('.crosslink.sites.bed', '.pureclip.log')]
    cmder.run(cmd, msg=f'Calling peaks from {bam} {size(bam)} using pureCLIP ...', pmt=True)


def peak(ip_bams, input_bams, peak_beds, reproducible_bed, outdir, cwd=''):
    cmd = ['peak', '--ip_bams', ' '.join(ip_bams),
           '--input_bam', ' '.join(input_bams),
           '--peak_beds', ' '.join(peak_beds),
           '--read_type', 'PE' if TYPE == 'paired' else 'SE',
           '--species', 'hg19' if options.species in ('hg19', 'hg19chr19') else options.species,
           '--outdir', outdir, '--cores', options.cores]
    cwd = cwd if cwd else os.path.dirname(reproducible_bed)
    cmder.run(cmd, cwd=cwd, stdout=sys.stdout, stderr=sys.stderr)
    return reproducible_bed


@task(inputs=[], outputs=f'{ECLIP}/{".vs.".join([s.key for s in SAMPLES])}.reproducible.peaks.bed', parent=clipper)
def reproducible_peaks(inputs, outputs):
    ip_bams, input_bams, peak_beds = [], [], []
    for sample in SAMPLES:
        ip_bams.append(os.path.basename(sample.ip_read.bam))
        input_bams.append(os.path.basename(sample.input_read.bam))
        peak_beds.append(os.path.basename(sample.bed))
    peak(ip_bams, input_bams, peak_beds, outputs, os.path.abspath(ECLIP))


def split_bam(bam, bam1, bam2):
    def count_mapped_reads(bam):
        count = int(cmder.run(f'samtools view -c -F 0x4 {bam}', msg='').stdout.read())
        logger.info(f'Found {count:,} mapped reads in {bam}.')
        return count

    if os.path.isfile(bam1) and os.path.isfile(bam2):
        logger.info(f'BAMs {bam1} and {bam2} already exist.')
    else:
        half_lines = int(count_mapped_reads(bam) / 2) + 1
        cmd = f'samtools view {bam} | shuf | split -d -l {half_lines} - {bam}'
        cmder.run(cmd, msg=f'Shuffling and splitting {bam} ...')
        tmp_bam1, tmp_bam2 = bam1.replace('.bam', '.tmp.bam'), bam2.replace('.bam', '.tmp.bam')
        cmd = f'samtools view -H {bam} | cat - {bam}00 | samtools view -bS - > {tmp_bam1}'
        
        cmder.run(cmd, msg=f'Creating headers for {bam1} ...')
        cmder.run(f'samtools sort -@ {options.cores} -o {bam1} {tmp_bam1}')
        cmd = f'samtools view -H {bam} | cat - {bam}01 | samtools view -bS - > {tmp_bam2}'
        
        cmder.run(cmd, msg=f'Creating headers for {bam2} ...')
        cmder.run(f'samtools sort -@ {options.cores} -o {bam2} {tmp_bam2}')
        cmder.run(f'rm {bam}00 {bam}01 {tmp_bam1} {tmp_bam2}')
    return bam1, bam2


def count_lines(file):
    lines = int(cmder.run(f'wc -l {file}').stdout.read().split()[0])
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

    key = ".".join([sample.ip_read.name for sample in SAMPLES])
    pseudo_reproducible_bed = f'rescue/{key}.pseudo.01.vs.{key}.pseudo.02.reproducible.peaks.bed'
    peak(pseudo_ip_bams, pseudo_input_bams, pseudo_peak_beds, pseudo_reproducible_bed, 'rescue', cwd=options.outdir)
    pseudo_count = count_lines(pseudo_reproducible_bed)

    key = ".vs.".join([sample.ip_read.name for sample in SAMPLES])
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
    for i, sample in enumerate(SAMPLES, start=1):
        split_ip_bams = split_bam(sample.ip_read.bam,
                                  f'consistency/{sample.ip_read.name}.split.01.bam',
                                  f'consistency/{sample.ip_read.name}.split.02.bam')
        split_input_bams = split_bam(sample.input_read.bam,
                                     f'consistency/{sample.input_read.name}.split.01.bam',
                                     f'consistency/{sample.input_read.name}.split.02.bam')
        split_peak_beds = [clipper_peaks(split_ip_bams[0]), clipper_peaks(split_ip_bams[1])]

        bed = f'consistency/{sample.ip_read.name}.split.01.vs.{sample.ip_read.name}.split.02.reproducible.peaks.bed'
        peak(split_ip_bams, split_input_bams, split_peak_beds, bed,  'consistency', cwd=options.outdir)
        counts.append(count_lines(bed))

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
        cmder.run(cmd, msg=f'Checking reads in {fastq} {size(fastq)} using falco ...')
        cmder.run(f'mv {tmp}/fastqc_data.txt {txt}')
    finally:
        shutil.rmtree(tmp)


@task(inputs=[], outputs=[f'{fastq}.gz' for fastq in glob.glob(f'{ECLIP}/*.trim.sort.r*.fastq')], parent=falco)
def pigz(fastq, gz):
    print(f'Compressing {fastq}, which {os.path.exists(gz)}')
    if options.debug:
        cmd = f'pigz -p {PROCESSES} {fastq}'
    else:
        cmd = f'pigz --processes {PROCESSES} {fastq}'
    cmder.run(cmd, msg=f'Compressing {fastq} ...', pmt=True)
    return gz


@task(inputs=[], outputs=f'{QC}/summary.html', parent=falco, mkdir_before_run=[QC])
def summary(inputs, outputs):
    report(ECLIP, QC, READS, SAMPLES, TYPE)


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
    pass
    main()
    # with open('/storage/vannostrand/software/eclip/test/paired/qc/summary.json') as f:
    #     counts = json.load(f)
    #     report(counts)

