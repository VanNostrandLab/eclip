#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
A pipeline for processing eCLIP data to identify genomic locations of RNA binding proteins (RBPs).
"""

import argparse
import glob
import os
import time
import subprocess
import json
import yaml
import datetime
import tempfile
import shutil
import math
import itertools
import pysam as pysam
import pandas as pd
import cmder
from seqflow import Flow, task, logger


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
sp.add_argument('--track', type=str, help="Name for the UCSC Genome Browser track, default: eCLIP", default='eCLIP')
sp.add_argument('--track_label', type=str, help="Label for the track.", default='eCLIP')
sp.add_argument('--track_genome', type=str, help="Genome name (a short name code) for the track.", default='hg19')
sp.add_argument('--l2fc', type=int, help="Only consider peaks at or above this log2 fold change cutoff.", default=3)
sp.add_argument('--l10p', type=int, help="Only consider peaks at or above this log10p value cutoff.", default=3)
sp.add_argument('--enrichment_filter', type=int, help="Pre-filter peaks that are enriched over input.", default=0)

js = parser.add_argument_group('Job Submit')
js.add_argument('--job', type=str, help="Name of your job", default='eCLIP')
js.add_argument('--email', type=str, help='Email address for notifying you the start, end, and abort of you job.')
js.add_argument('--time', type=int, help='Time (in integer hours) for running your job.', default=24)
js.add_argument('--memory', type=int, help='Amount of memory (in GB) for all CPU cores.', default=32)
js.add_argument('--cores', type=int, help='Maximum number of CPU cores can be used for your job.', default=8)
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
    
    setattr(options, 'randomer_length', options.randomer_length or manifest.get('randomer_length', '5'))
    
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
        if self.fastq1:
            if not os.path.isfile(self.fastq1):
                raise ValueError(f'{self.fastq1} may not be a file or does not exist.')
        else:
            raise KeyError(f'No fastq or fastq1 was assigned for {self.type} read\n{record}.')
        self.fastq2 = record.get('fastq2', '')
        if self.fastq2 and not os.path.isfile(self.fastq2):
            raise ValueError(f'{self.fastq2} may not be a file or does not exist.')
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
        max_size = max(sizes) / (1000 * 1000 * 1000) * 2
        n = int(options.memory / max_size)
        if n == 0:
            n = 1
        elif n > options.cores:
            n = options.cores
        return n

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

fastqs, links = [], []
for read in READS.values():
    fastqs.append(read.fastq1)
    links.append(read.link1)


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


def prepare_reads_outputs(expand=False):
    outputs = set()
    for r in READS.values():
        barcodes = r.barcodes if r.paired else ['umi', 'umi']
        outputs.add(f'{ECLIP}/{r.key}.{barcodes[0]}.r1.fastq.gz')
        if expand:
            outputs.add(f'{ECLIP}/{r.key}.{barcodes[1]}.r1.fastq.gz')
    return list(outputs)


def prepare_reads_cleanup():
    kps = prepare_reads_outputs(expand=True) + [r.link1 for r in READS.values()] + [r.link2 for r in READS.values()]
    rms = [p for p in glob.iglob(f'{ECLIP}/*.fastq.gz') if p not in kps]
    return rms


@task(inputs=soft_link, outputs=prepare_reads_outputs(),
      checkpoint=False, processes=options.cores, cleanup_after_run=prepare_reads_cleanup())
def prepare_reads(link, output):
    """Extract UMIs for single-end reads or demultiplex paired-end reads."""
    read = READS[os.path.basename(link.replace('.r1.fastq.gz', ''))]
    fastq1, fastq2 = read.link1, read.link2
    if fastq2:
        barcodes = read.barcodes
        message = f'Demultiplexing paired-end reads {fastq1} {size(fastq1)} and\n{" " * 43}{fastq2} {size(fastq2)} ...'
        cmd = ['demux',
               '--fastq_1', os.path.basename(fastq1),
               '--fastq_2', os.path.basename(fastq2),
               '--newname', read.name,
               '--expectedbarcodeida', barcodes[0],
               '--expectedbarcodeidb', barcodes[1],
               '--barcodesfile', options.barcodes_fasta,
               '--length', options.randomer_length,
               '--metrics', f'{options.dataset}.{read.name}.demux.metrics']
        cmder.run(cmd, msg=message, pmt=True, cwd=os.path.join(options.outdir, ECLIP))
    else:
        message = f'Extract UMIs for single-end read {fastq1} {size(fastq1)} ...'
        cmd = ['umi_tools', 'extract',
               '--random-seed', 1,
               '--stdin', fastq1,
               '--bc-pattern', 'NNNNNNNNNN',
               '--log', fastq1.replace('.fastq.gz', '.extract.metrics'),
               '--stdout', fastq1.replace('.r1.fastq.gz', '.umi.r1.fastq.gz')]
        cmder.run(cmd, msg=message, pmt=True)


@task(inputs=prepare_reads_outputs(expand=True), outputs=lambda i: i.replace('.r1.fastq.gz', '.trim.r1.fastq'),
      parent=prepare_reads,
      cleanup_after_run=[i.replace('.trim.', '.clean,') for i in prepare_reads_outputs(expand=True)] +
                        [i.replace('.trim.', '.clean,').replace('.r1.', '.r2.')
                         for i in prepare_reads_outputs(expand=True)])
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
            folder = tempfile.mkdtemp(dir=os.path.join(options.outdir, 'eclip'))
            cmder.run(cmd, msg='Parsing barcodes and finding adapters ...', cwd=folder)
            adapters = parse_adapters('-a', os.path.join(folder, 'a_adapters.fasta'))
            adapters += parse_adapters('-A', os.path.join(folder, 'A_adapters.fasta'))
            adapters += parse_adapters('-g', os.path.join(folder, 'g_adapters.fasta'))
            adapters1 = adapters
            adapters2 = parse_adapters('-A', os.path.join(folder, 'A_adapters.fasta'))
            overlap1 = parse_overlap(os.path.join(folder, 'trim_first_overlap_length.txt'))
            overlap2 = parse_overlap(os.path.join(folder, 'trim_again_overlap_length.txt'))
            shutil.rmtree(folder)
        return adapters1, adapters2, overlap1, overlap2

    def get_ios(input1, input2, out1, out2):
        if '.umi.' in out1:
            tmp = out1.replace('.trim.', '.clean.')
            ios1 = ['-o', tmp, input1, '>', out1.replace('.r1.trim.fastq', '.trim.first.metrics')]
            msg1 = f"Cutting adapters for single read {input1} {size(input1)} (first round) ..."
            
            ios2 = ['-o', out1, tmp, '>', out1.replace('.r1.clean.fastq', '.trim.second.metrics')]
            msg2 = f"Cutting adapters for single read {tmp} {size(tmp)} (second round) ..."
        else:
            tmp1, tmp2 = input1.replace('.trim.', '.clean.'), input2.replace('.trim.', '.clean.')
            tmp_metrics = input1.replace('.r1.fastq.gz', '.trim.first.metrics')
            ios1 = ['-o', out1, '-p', out2, tmp1, tmp2, '>', tmp_metrics]
            msg1 = (f"Cutting adapters for paired reads {input1} {size(input1)} and\n{' ' * 45}"
                    f"{input2} {size(input2)} (first round) ...")

            metrics = input1.replace('.r1.fastq.gz', '.trim.second.metrics')
            ios2 = ['-o', out1, '-p', out2, tmp1, tmp2, '>', metrics]
            msg2 = (f"Cutting adapters for paired reads {tmp1} {size(tmp1)} and\n{' ' * 45}"
                    f"{tmp2} {size(tmp2)} (second round) ...")
        return ios1, ios2, msg1, msg2

    def trim_adapters(adapters, overlap, ios, message):
        cmd = ['cutadapt', '-O', overlap, '--times', '2', '-e', '0.0', '-j', options.cores, '-m', '18',
               '--quality-cutoff', '6', '--match-read-wildcards'] + adapters + ios
        cmder.run(cmd, msg=message, pmt=True)

    key = os.path.basename(r1).rsplit('.', maxsplit=4)[0]
    read = READS[key]
    adapters1, adapters2, overlap1, overlap2 = get_adapters(read)
    ios1, ios2, msg1, msg2 = get_ios(r1, r1.replace('.r1.', '.r2.') if read.paired else '',
                                     fastq, fastq.replace('.r1.', '.r2.') if read.paired else '')
    trim_adapters(adapters1, overlap1, ios1, msg1)
    trim_adapters(adapters2, overlap2, ios2, msg2)
    return fastq


@task(inputs=cut_adapt, outputs=lambda i: i.replace('.trim.r1.fastq', '.trim.sort.r1.fastq'), processes=options.cores)
def sort_fastq(fastq, output):
    cmd = ['fastq-sort', '--id', fastq, '>', output]
    cmder.run(cmd, msg=f'Sort fastq {fastq} {size(fastq)} ...', pmt=True)

    fastq2 = fastq.replace('.r1.', '.r2.')
    if os.path.isfile(fastq2):
        output2 = output.replace('.r1.clean.sort.fastq', '.r2.clean.sort.fastq')
        cmd = ['fastq-sort', '--id', fastq2, '>', output2]
        cmder.run(cmd, msg=f'Sort fastq {fastq2} {size(fastq2)} ...', pmt=True)
    return output


@task(inputs=[i.replace('.r1.fastq.gz', '.trim.sort.r1.fastq') for i in prepare_reads_outputs(expand=True)],
      parent=sort_fastq, mkdir_before_run=[f'{ECLIP}/repeat.elements.map'],
      outputs=lambda i: (f'{ECLIP}/repeat.elements.map/{os.path.basename(i).replace(".trim.sort.r1.fastq", "")}'
                         f'/Unmapped.out.mate1'))
def map_to_repeat_elements(fastq, mate1):
    fastq1, fastq2 = fastq, fastq.replace('.trim.sort.r1.fastq', '.trim.sort.r2.fastq')
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
    cmd = ['fastq-sort', '--id', mate1, '>', output]
    cmder.run(cmd, msg=f'Sort mate1 {mate1} {size(mate1)} ...', pmt=True)

    mate2, output2 = mate1.replace('.mate1', '.mate2'), output.replace('.mate1', '.mate2')
    if os.path.isfile(mate2):
        output2 = output.replace('.out.mate1.sort.fastq', '.out.mate2.sort.fastq')
        cmd = ['fastq-sort', '--id', mate2, '>', output2]
        cmder.run(cmd, msg=f'Sort mate2 {mate2} {size(mate2)} ...', pmt=True)
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
    cmder.run(f'samtools sort -n -@ {options.cores} -m {min(4, int(options.memory / PROCESSES))}G -o {out} {bam}',
              msg=f'Name sorting {bam} {size(bam)} ...')


def position_sort_bam(bam, out):
    cmder.run(f'samtools sort -@ {options.cores} -m {min(4, int(options.memory / PROCESSES))}G -o {out} {bam}',
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


@task(inputs=prepare_bam, outputs=lambda i: i.replace('.sort.bam', '.sort.dedup.bam'))
def dedup_bam(bam, out):
    """Collapse barcodes of paired-end bam or umi_tools dedup single-end bam."""
    if TYPE == 'single':
        cmd = ['umi_tools', 'dedup', '--random-seed', 1, '--stdin', bam, '--method', 'unique',
               '--output-stats', out.replace(".dedup.bam", ".dedup.metrics"), '--stdout', out]
        message = f'Deduplicating {bam} {size(bam)} by umi_tools dedup ...'
    else:
        cmd = f'barcode_collapse.py -o {out} -m {out.replace(".bam", ".collapse.metrics")} -b {bam}'
        message = f'Deduplicating {bam} {size(bam)} by collapsing barcodes ...'
    cmder.run(cmd, msg=message, pmt=True)


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
    if options.strand_direction in ('f', 'forward'):
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
    
    
def call_peaks(bam, bed=''):
    bed = bed if bed else bam.replace('.bam', '.peak.clusters.bed')
    cmd = f'clipper --species {options.species} --processors {options.cores} --bam {bam} --outfile {bed}'
    cmder.run(cmd, msg=f'Calling peaks from {bam} {size(bam)} using clipper ...', pmt=True)
    return bed


@task(inputs=[read.bam for read in READS.values() if read.type == 'IP'],
      outputs=lambda i: i.replace('.bam', '.peak.clusters.bed'), parent=index_bam)
def clipper(bam, bed):
    return call_peaks(bam, bed)


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


@task(inputs=[], outputs=f'{".vs.".join([s.key for s in SAMPLES])}.reproducible.peaks.bed', parent=clipper)
def reproducible_peaks(inputs, outputs):
    ip_bams, input_bams, peak_beds = ['--ip_bams'], ['--input_bams'], ['--peak_beds']
    for sample in SAMPLES:
        input_bams.append(sample.ip_read.bam)
        input_bams.append(sample.input_read.bam)
        input_bams.append(sample.bed)
    cmd = ['peak.py'] + ip_bams + input_bams + peak_beds
    cmd += ['--cores', options.cores] + ['--species', options.species] + ['--outdir', options.outdir]
    cmd.run(cmd)
    
    
def split_bam(bam, bam1, bam2):
    with pysam.AlignmentFile(bam, 'rb') as sam:
        half_lines = int(sam.mapped / 2)
    cmd = f'samtools view {bam} | shuf | split -d -I {half_lines} - {bam}'
    cmder.run(cmd, msg=f'Shuffling and splitting {bam} ...')
    
    cmd = f'samtools view -H {bam} | cat - {bam}00 | samtools view -bS - > {bam1}'
    cmder.run(cmd, msg=f'Adding headers for {bam1} ...')
    cmd = f'samtools view -H {bam} | cat - {bam}00 | samtools view -bS - > {bam2}'
    cmder.run(cmd, msg=f'Adding headers for {bam2} ...')
    return bam1, bam2
    
    
@task(inputs=[], outputs=f'{ECLIP}/rescue.ratio.txt')
def rescue_ratio(inputs, outputs):
    def prepare_pseudo_bam(bam1, bam2, basename):
        pseudo_bam = f'{basename}.bam'
        cmd = f'samtools merge {pseudo_bam} {bam1} {bam2}'
        cmder.run(cmd, msg=f'Merging {bam1} and {bam2} ...')

        bam1, bam2 = split_bam(pseudo_bam, f'{basename}.pseudo.01.bam', f'{basename}.pseudo.02.bam')
        return bam1, bam2
    
    def count_lines(file):
        with open(file) as f:
            return sum(1 for line in f)
    
    pseudo_ip_bams, pseudo_input_bams, pseudo_peak_beds, keys, names = [], [], [], [], []
    for sample1, sample2 in itertools.combinations(SAMPLES, 2):
        ip_bam1, ip_bam2 = sample1.ip_read.bam, sample2.ip_read.bam
        ip_basename = f'{ECLIP}/{sample1.key}.{sample2.key}'
        keys.append(sample1.key)
        names.append(os.path.basename(ip_basename))
        pseudo_ip_bam = prepare_pseudo_bam(ip_bam1, ip_bam2, ip_basename)
        pseudo_ip_bams.extend(pseudo_ip_bam)
        
        pseudo_peak_beds.extend([call_peaks(bam) for bam in pseudo_ip_bam])
            
        input_bam1, input_bam2 = sample1.input_read.bam, sample2.input_read.bam
        input_basename = f'{ECLIP}/{sample1.key}.{sample2.key}'
        pseudo_input_bam = prepare_pseudo_bam(input_bam1, input_bam2, input_basename)
        pseudo_input_bams.extend(pseudo_input_bam)

    cmd = ['peak.py'] + pseudo_ip_bams + pseudo_input_bams + pseudo_peak_beds
    cmd += ['--cores', options.cores] + ['--species', options.species] + ['--outdir', options.outdir]
    cmd.run(cmd)
    
    count = count_lines(f'{ECLIP}/{".vs.".join(keys)}.reproducible.peaks.bed')
    pseudo_count = count_lines(f'{ECLIP}/{".vs.".join(names)}.reproducible.peaks.bed')
    ratio = max(count, pseudo_count) / min(count, pseudo_count)
    with open(outputs, 'w') as o:
        o.write(f'{ratio}')


@task(inputs=[], outputs=f'{ECLIP}/self.consistency.ratio.txt')
def consistency_ration(inputs, outputs):
    counts = []
    for (ip_read, input_read) in SAMPLES:
        names = [f'{ip_read.key}.split.{i:02d}' for i in range(2)]
        ip_bam1, ip_bam2 = [f'{ECLIP}/{name}.bam' for name in names]
        split_ip_bams = split_bam(ip_read.bam, ip_bam1, ip_bam2)
        
        split_peak_beds = [call_peaks(split_ip_bams[0]), call_peaks(split_ip_bams[1])]
        
        input_bam1, input_bam2 = [f'{ECLIP}/{name}.bam' for name in names]
        split_input_bams = split_bam(input_read.bam, input_bam1, input_bam2)

        cmd = ['peak.py'] + split_ip_bams + split_input_bams + split_peak_beds
        cmd += ['--cores', options.cores] + ['--species', options.species] + ['--outdir', options.outdir]
        cmd.run(cmd)
        
        counts.append(f'{ECLIP}/{".vs.".join(names)}.reproducible.peaks.bed')
        
    ratio = counts[0] / counts[1]
    with open(outputs, 'w') as o:
        o.write(f'{ratio}')
        

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
        flow = Flow('eCLIP', description=__doc__.strip())
        flow.run(dry=options.dry_run, processes=options.cores, verbose=options.verbose)
        logger.debug('')
        run_time = str(datetime.timedelta(seconds=int(time.perf_counter() - START_TIME)))
        logger.trace(FOOTER.format(hh_mm_ss=f'time consumed: {run_time}'.upper().center(118)))


if __name__ == '__main__':
    main()
