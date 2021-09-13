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

import datetime
import tempfile
import shutil
import itertools
import math

import pysam
import cmder
from seqflow import Flow, task, logger


def file_path(p):
    if p and not os.path.isfile(p):
        logger.error(f'File "{p}" may not be a file or does not exist.')
        sys.exit(1)
    return p


def dir_path(p):
    if p and not os.path.isdir(p):
        logger.error(f'Path "{p}" may not be a directory or does not exist.')
        sys.exit(1)
    return p


parser = argparse.ArgumentParser(description=__doc__, prog='eclip')
parser.add_argument('--input_fastqs', nargs='+', required=True, type=file_path, help='Path(s) to INPUT FASTQ file(s).')
parser.add_argument('--ip_fastqs', nargs='+', required=True, type=file_path, help='Path(s) to IP FASTQ file(s).')
parser.add_argument('--names', nargs='+', required=True, help='Shortnames for each sample, e.g., rep1, rep2.')
parser.add_argument('--adapters_fasta', help="Path to the FASTA file contains adapters and their sequences.")
parser.add_argument('--barcodes_pattern', help="Pattern of barcodes for umi-tools using to extract UMIs (for "
                                               "single-end dataset only, default: NNNNNNNNNN", default='NNNNNNNNNN')
parser.add_argument('--dataset', help='Name of the eCLIP dataset, default: eCLIP.', default='eCLIP')
parser.add_argument('--outdir', help="Path to the output directory, an pre-existing directory or current directory,"
                                     "default: current work directory.")
parser.add_argument('--species', help="Species name (short name code) the dataset associated with, e.g., hg19, mm10.",
                    default='hg19')
parser.add_argument('--genome', help="Path to STAR reference genome index directory.", type=dir_path)
parser.add_argument('--repeat', help="Path to STAR repeat elements index directory.", type=dir_path)
parser.add_argument('--blacklist_bed', help="Path to the bed file contains blacklist.")
parser.add_argument('--track', help="Name for the UCSC Genome Browser track, default: same as dataset")
parser.add_argument('--track_genome', help="Genome name (a short name code) for the track, default: same as species.")
parser.add_argument('--l2fc', type=int, help="Only consider peaks at or above this log2 fold change cutoff.", default=3)
parser.add_argument('--l10p', type=int, help="Only consider peaks at or above this log10 p value cutoff.", default=3)

parser.add_argument('--job', help="Name of your job", default='eCLIP')
parser.add_argument('--email', help='Email address for notifying you the start, end, and abort of you job.')
parser.add_argument('--time', type=int, help='Time (in integer hours) for running your job.', default=36)
parser.add_argument('--memory', type=int, help='Amount of memory (in GB) for all CPU cores.', default=32)
parser.add_argument('--cpus', type=int, help='Maximum number of CPU cores can be used for your job.', default=16)
parser.add_argument('--scheduler', choices=('pbs', 'qsub', 'slurm', 'sbatch', 'local'), default='slurm',
                    help='Name (case insensitive) of the scheduler on your cluster.')
parser.add_argument('--hold_submit', action='store_true',
                    help='Generate the submit script but hold it without submitting to the job scheduler.')

parser.add_argument('--debug', action='store_true', help='Invoke debug mode (for development use only).')
parser.add_argument('--dry_run', action='store_true',
                    help='Print out steps and files involved in each step without actually running the pipeline.')

START_TIME = time.perf_counter()


options = parser.parse_args()
setattr(options, 'outdir', options.outdir or os.getcwd())
if not os.path.isdir(options.outdir):
    try:
        os.mkdir(options.outdir)
    except OSError as e:
        logger.error(f'Create outdir failed: {e}.')
        sys.exit(1)
os.chdir(options.outdir)

adapters = '/storage/vannostrand/software/eclip/data/se.adapters.fasta'
setattr(options, 'adapters_fasta', options.adapters_fasta or adapters)
file_path(options.adapters_fasta)

setattr(options, 'repeat', options.repeat or '/storage/vannostrand/reference_data/hg19/repbase_v2_star_index')
dir_path(options.repeat)

setattr(options, 'genome', options.genome or '/storage/vannostrand/reference_data/hg19/genome_star_index')
dir_path(options.genome)

setattr(options, 'track', options.track or options.dataset)
setattr(options, 'track_genome', options.track_genome or options.species)

ips, inputs, names = options.ip_fastqs, options.input_fastqs, options.names
if len(ips) == len(names):
    if len(ips) == len(inputs):
        pass
    else:
        if len(inputs) == 1:
            inputs = inputs * len(ips)
        else:
            logger.error('Wrong number of input_fastqs were provided.')
            sys.exit(1)
else:
    logger.error('Number of items in ip_fastqs and names are not equal.')
    sys.exit(1)


class Sample:
    def __init__(self, sample_name):
        self.name = sample_name
        self.ip_fastq = f'{sample_name}.ip.fastq.gz'
        self.input_fastq = f'{sample_name}.input.fastq.gz'
        self.ip_bam = f'{sample_name}.ip.bam'
        self.input_bam = f'{sample_name}.input.bam'
        self.peak_bed = f'{sample_name}.peak.clusters.bed'
        self.cross_bed = f'{sample_name}.crosslink.sites.bed'
        
        
FASTQS, SAMPLES = {}, []
for ip_fastq, input_fastq, name in zip(ips, inputs, names):
    FASTQS[ip_fastq] = f'{name}.ip.fastq.gz'
    FASTQS[input_fastq] = f'{name}.input.fastq.gz'
    sample = Sample(name)
    SAMPLES.append(sample)

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


@task(inputs=[fastq for fastq in FASTQS], outputs=lambda i: FASTQS[i])
def soft_link(fastq, link):
    """ Create soft links for original fastq files. """
    path = os.path.abspath(fastq)
    if path == os.path.abspath(link):
        logger.warning(f"No symbolic link was made for {path}! You are directly working on the original file!")
    else:
        os.symlink(path, link)
    return link


@task(inputs=soft_link, outputs=lambda i: i.replace('.fastq.gz', '.umi.fastq.gz'), cpus=options.cpus)
def extract_umi(fastq, umi):
    message = f'Extract UMIs for {fastq} ...'
    cmd = ['umi_tools', 'extract',
           '--random-seed', 1,
           '--stdin', fastq,
           '--bc-pattern', 'NNNNNNNNNN',
           '--log', fastq.replace('.fastq.gz', '.umi.metrics'),
           '--stdout', umi]
    cmder.run(cmd, msg=message, pmt=True)


@task(inputs=extract_umi, outputs=lambda i: i.replace('.umi.fastq.gz', '.trim.fastq.gz'))
def cut_adapt(fastq, trimmed_fastq):
    """ Trimming adapters (2 rounds) and 10 NTs from 3' using cutadapt (reads < 24 NTs will be discarded!)"""
    
    trim_tmp = trimmed_fastq.replace('.trim.fastq.gz', '.trim.tmp.fastq.gz')
    trim_trim_tmp = trimmed_fastq.replace('.trim.fastq.gz', '.trim.trim.tmp.fastq.gz')
    try:
        trim_metrics = trimmed_fastq.replace('.trim.fastq.gz', '.trim.metrics')
        trim_trim_metrics = trimmed_fastq.replace('.trim.fastq.gz', '.trim.trim.metrics')
        trim_trim_trim_metrics = trimmed_fastq.replace('.trim.fastq.gz', '.trim.trim.trim.metrics')
        cmd = ['cutadapt',
               '-j', options.cpus,
               '--match-read-wildcards',
               '--times', 1,
               '-e', 0.1,
               '--quality-cutoff', 6,
               '-m', 24,
               '-a', f'file:{options.adapters_fasta}']
        
        cmder.run(cmd + ['-O', 1, '-o', trim_tmp, fastq, '>', trim_metrics],
                  msg=f"Trimming adapters (1st round) for {fastq} ...", pmt=False)
        cmder.run(cmd + ['-O', 5, '-o', trim_trim_tmp, trim_tmp, '>', trim_trim_metrics],
                  msg=f"Trimming adapters (2nd round) for {fastq} ...", pmt=False)
        
        cmd = ['cutadapt', '-j', options.cpus, '-u', -10, '-o', trimmed_fastq, trim_trim_tmp,
               '>', trim_trim_trim_metrics]
        cmder.run(cmd, msg=f"Trimming adapters (3d round, 3'-UMIs) for {fastq} ...", pmt=False)
    finally:
        if os.path.exists(trim_tmp):
            os.unlink(trim_tmp)
        if os.path.exists(trim_trim_tmp):
            os.unlink(trim_trim_tmp)


@task(inputs=cut_adapt, outputs=lambda i: i.replace('.trim.fastq.gz', '.repeat.unmap.fastq.gz'))
def map_to_repeat_elements(fastq, mate):
    prefix = fastq.replace('.trim.fastq.gz', '.repeat.map')
    try:
        if not os.path.isdir(prefix):
            os.mkdir(prefix)
        cmd = ['STAR',
               '--runMode', 'alignReads',
               '--runThreadN', options.cpus,
               '--alignEndsType', 'EndToEnd',
               '--genomeDir', options.repeat,
               '--genomeLoad', 'NoSharedMemory',
               '--outBAMcompression', 10,
               '--outFileNamePrefix', f"{prefix}/",
               '--outFilterMultimapNmax', 100,
               '--outFilterMultimapScoreRange', 1,
               '--outFilterScoreMin', 10,
               '--outFilterType', 'BySJout',
               '--outReadsUnmapped', 'Fastx',
               '--outSAMattrRGline', 'ID:foo',
               '--outSAMattributes', 'All',
               '--outSAMmode', 'Full',
               '--outSAMtype', 'BAM', 'Unsorted',
               '--outSAMunmapped', 'None',
               '--outStd', 'Log',
               '--readFilesCommand', 'zcat',
               '--readFilesIn', fastq]
        message = f'Map SE reads in {fastq} to repeat elements ...'
        cmder.run(cmd, msg=message)
        cmder.run(f'mv {prefix}/Log.final.out {mate.replace(".repeat.unmap.fastq.gz", ".repeat.map.log")}')
        cmder.run(f'mv {prefix}/Aligned.out.bam {mate.replace(".repeat.unmap.fastq.gz", ".repeat.map.bam")}')
        cmder.run(f'pigz -c -p {options.cpus} {prefix}/Unmapped.out.mate1 > {mate}')
    finally:
        shutil.rmtree(prefix)


@task(inputs=map_to_repeat_elements, outputs=lambda i: i.replace('.repeat.unmap.fastq.gz', '.genome.map.bam'))
def map_to_reference_genome(mate, bam):
    prefix = mate.replace('.repeat.unmap.fastq.gz', '.genome.map')
    try:
        if not os.path.isdir(prefix):
            os.mkdir(prefix)
        cmd = ['STAR',
               '--runMode', 'alignReads',
               '--runThreadN', options.cpus,
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
               '--outSAMunmapped', 'None',
               '--outStd', 'Log',
               '--readFilesCommand', 'zcat',
               '--readFilesIn', mate]
        message = f'Map SE repeat elements unmapped reads in {mate} to reference genome ...'
        cmder.run(cmd, msg=message)
        cmder.run(f'mv {prefix}/Log.final.out {bam.replace(".genome.map.bam", ".genome.map.log")}')
        unmap = bam.replace('.genome.map.bam', '.genome.unmap.fastq.gz')
        cmder.run(f'pigz -c -p {options.cpus} {prefix}/Unmapped.out.mate1 > {unmap}')
        cmder.run(f'samtools sort -@ {options.cpus} -m 2G -o {bam} {prefix}/Aligned.out.bam',
                  msg=f'Sorting {prefix}/Aligned.out.bam to {bam} ...')
        cmder.run(f'samtools index {bam}', msg=f'Indexing {bam} ...')
    finally:
        shutil.rmtree(prefix)
    return bam


@task(inputs=map_to_reference_genome, outputs=lambda i: i.replace('.genome.map.bam', '.bam'))
def dedup_bam(bam, out):
    """Deduplicate SE BAM using umi_tools dedup."""
    cmd = ['umi_tools', 'dedup', '--random-seed', 1, '--stdin', bam, '--method', 'unique', '--stdout', out]
    cmder.run(cmd, msg=f'Deduplicating {bam} by umi_tools dedup ...')
    cmder.run(f'samtools index {out}', msg=f'Indexing {bam} ...')


@task(inputs=dedup_bam, cpus=options.cpus, outputs=lambda i: i.replace('.bam', '.plus.bw'))
def make_bigwig_files(bam, bigwig):
    def bam_to_bigwig(bam, scale, strand, bw):
        bg = bw.replace('.bw', '.bg')
        cmd = f'genomeCoverageBed -ibam {bam} -bg -scale {scale} -strand {strand} -du -split | sort -k1,1 -k2,2n > {bg}'
        cmder.run(cmd, msg=f'Calculating genome coverage for {bam} ({strand} strand) ...')
        cmd = f'bedGraphToBigWig {bg} {options.genome}/chrNameLength.txt {bw}'
        cmder.run(cmd, msg=f'Converting {bg} to {bw} ...')
        cmder.run(f'rm {bg}')
    
    logger.info(f'Make BigWig files for {bam} ...')
    pos_bw, neg_bw = bigwig, bigwig.replace('.plus.bw', '.minus.bw')
    with pysam.AlignmentFile(bam, 'rb') as sam:
        total_reads = sam.mapped
    try:
        scale = 1000000.0 / total_reads
    except ZeroDivisionError:
        logger.warning(f'No reads was found in BAM {bam}, empty BigWig file was created.')
        with open(bigwig, 'w') as o:
            o.write('')
        return bigwig
    bam_to_bigwig(bam, scale, '+', pos_bw)
    bam_to_bigwig(bam, -1 * scale, '-', neg_bw)
    logger.info(f'Make BigWig files for {bam} complete.')
    return bigwig


@task(inputs=[], outputs=['hub.txt'], parent=make_bigwig_files)
def make_hub_file(inputs, output):
    logger.info('Make hub track file ...')
    header = f"""hub {options.track.replace(' ', '_')}
shortLabel {options.track}
longLabel {options.track}
useOneFile on
email {options.email if options.email else 'fei.yuan@bcm.edu'}

genome {options.track_genome}

track {options.track.replace(' ', '_')}
shortLabel {options.track}
longLabel {options.track}
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
        for bw in glob.iglob('*.plus.bw'):
            key = bw.replace('.plus.bw', '')
            plus, minus = f'{key}.pos.bw', f'{key}.neg.bw'
            name1, name2 = f'{key} plus', f'{key} minus'
            o.write(block.format(track=track, name1=name1, name2=name2, basename=key, plus=plus, minus=minus))
    logger.info('Make hub track file complete.')


def clipper_peaks(bam, bed=''):
    bed = bed if bed else bam.replace('.ip.bam', '.peak.clusters.bed')
    if os.path.isfile(bed):
        logger.info(f'Clipper bed {bed} already exists.')
    else:
        cmd = f'clipper --species {options.species} --processors {options.cpus} --bam {bam} --outfile {bed}'
        cmder.run(cmd, msg=f'Calling peaks from {bam} using clipper ...', pmt=True)
    return bed


@task(inputs=[], outputs=[sample.peak_bed for sample in SAMPLES], parent=dedup_bam)
def clipper(bam, bed):
    bam = bed.replace('.peak.clusters.bed', '.ip.bam')
    return clipper_peaks(bam, bed)


@task(inputs=[],
      outputs=[sample.cross_bed for sample in SAMPLES], parent=dedup_bam)
def pureclip(bam, bed):
    ip_bam, input_bam = [[sample.ip_bam, sample.input_bam] for sample in SAMPLES
                         if sample.cross_bed == bed][0]
    cmd = ['pureclip', '-i', ip_bam, '-bai', f'{ip_bam}.bai', '-g', f'{options.genome}/genome.fa',
           '-nt', options.cpus, '-ibam', input_bam, '-ibai', f'{input_bam}.bai',
           '-o', bed, '-or', bed.replace('.crosslink.sites.bed', '.binding.regions.bed'),
           '>', bed.replace('.crosslink.sites.bed', '.pureclip.log')]
    cmder.run(cmd, msg=f'Calling peaks from {bam} using pureCLIP ...')


def peak(ip_bams, input_bams, peak_beds, ids, reproducible_bed, outdir, cwd=''):
    cmd = ['peak', '--ip_bams', ' '.join(ip_bams),
           '--input_bam', ' '.join(input_bams),
           '--peak_beds', ' '.join(peak_beds),
           '--ids', ' '.join(ids),
           '--read_type', 'SE',
           '--species', 'hg19' if options.species in ('hg19', 'hg19chr19') else options.species,
           '--outdir', outdir, '--cores', options.cpus,
           '--l2fc', options.l2fc, '--l10p', options.l10p]
    cwd = cwd if cwd else os.path.dirname(reproducible_bed)
    cmder.run(cmd, cwd=cwd, stdout=sys.stdout, stderr=sys.stderr)
    return reproducible_bed


@task(inputs=[], outputs=[f'{".vs.".join([s.name for s in SAMPLES])}.reproducible.peaks.bed'], parent=clipper)
def reproducible_peaks(inputs, outputs):
    ip_bams, input_bams, peak_beds, ids = [], [], [], []
    for sample in SAMPLES:
        ip_bams.append(sample.ip_bam)
        input_bams.append(sample.input_bam)
        peak_beds.append(sample.peak_bed)
        ids.append(sample.name)
    peak(ip_bams, input_bams, peak_beds, ids, outputs, options.outdir, cwd=options.outdir)


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
        cmder.run(f'samtools sort -@ {options.cpus} -o {bam1} {tmp_bam1}')
        cmd = f'samtools view -H {bam} | cat - {bam}01 | samtools view -bS - > {tmp_bam2}'
        
        cmder.run(cmd, msg=f'Creating headers for {bam2} ...')
        cmder.run(f'samtools sort -@ {options.cpus} -o {bam2} {tmp_bam2}')
        cmder.run(f'rm {bam}00 {bam}01 {tmp_bam1} {tmp_bam2}')
    return bam1, bam2


def count_lines(file):
    lines = int(cmder.run(f'wc -l {file}').stdout.read().split()[0])
    logger.info(f'Found {lines:,} lines in {file}.')
    return lines


@task(inputs=[], outputs=['rescue.ratio.txt'], parent=clipper, mkdir=['rescue'])
def rescue_ratio(inputs, outputs):
    def prepare_pseudo_bam(bam1, bam2, basename):
        if os.path.exists(bam1) and os.path.exists(bam2):
            logger.info('Split BAM files already exist.')
        else:
            pseudo_bam = f'{basename}.bam'
            tmp_pseudo_bam = pseudo_bam.replace('.bam', '.tmp.bam')
            cmd = f'samtools merge {tmp_pseudo_bam} {bam1} {bam2}'
            cmder.run(cmd, msg=f'Merging {bam1} and {bam2} ...')
            
            cmder.run(f'samtools sort -@ {options.cpus} -m 2G -o {pseudo_bam} {tmp_pseudo_bam}')
            cmder.run(f'rm {tmp_pseudo_bam}')
            
            bam1, bam2 = split_bam(pseudo_bam, f'{basename}.pseudo.01.bam', f'{basename}.pseudo.02.bam')
        return bam1, bam2
    
    if len(SAMPLES) == 1:
        return
    pseudo_ip_bams, pseudo_input_bams, pseudo_peak_beds = [], [], []
    for i, (sample1, sample2) in enumerate(itertools.combinations(SAMPLES, 2), start=1):
        pseudo_ip_bam = prepare_pseudo_bam(sample1.ip_bam, sample2.ip_bam, f'rescue/{sample1.name}.{sample2.name}.ip')
        pseudo_ip_bams.extend(pseudo_ip_bam)
        
        pseudo_input_bam = prepare_pseudo_bam(sample1.input_bam, sample2.input_bam,
                                              f'rescue/{sample1.name}.{sample2.name}.input')
        pseudo_input_bams.extend(pseudo_input_bam)
        
        pseudo_peak_beds.extend([clipper_peaks(bam) for bam in pseudo_ip_bam])
    
    key = ".".join([sample.name for sample in SAMPLES])
    pseudo_reproducible_bed = f'rescue/{key}.ip.pseudo.01.vs.{key}.ip.pseudo.02.reproducible.peaks.bed'
    if not os.path.exists(pseudo_reproducible_bed):
        peak(pseudo_ip_bams, pseudo_input_bams, pseudo_peak_beds, pseudo_reproducible_bed, 'rescue', cwd=options.outdir)
    pseudo_count = count_lines(pseudo_reproducible_bed)
    
    count = count_lines(f'{".vs.".join([sample.name for sample in SAMPLES])}.reproducible.peaks.bed')
    try:
        ratio = max(count, pseudo_count) / min(count, pseudo_count)
    except ZeroDivisionError:
        ratio = 0
        logger.error(f'No peaks found in reproducible peaks or pseudo reproducible peaks, return ratio 0.')
    with open(outputs, 'w') as o:
        o.write(f'{ratio}\n')


@task(inputs=[], outputs=['consistency.ratio.txt'], parent=clipper, mkdir=['consistency'])
def consistency_ratio(inputs, outputs):
    if len(SAMPLES) == 1:
        return
    counts = []
    for i, sample in enumerate(SAMPLES, start=1):
        split_ip_bams = split_bam(sample.ip_bam,
                                  f'consistency/{sample.name}.ip.split.01.bam',
                                  f'consistency/{sample.name}.ip.split.02.bam')
        split_input_bams = split_bam(sample.input_bam,
                                     f'consistency/{sample.name}.input.split.01.bam',
                                     f'consistency/{sample.name}.input.split.02.bam')
        split_peak_beds = [clipper_peaks(split_ip_bams[0]), clipper_peaks(split_ip_bams[1])]
        
        bed = f'consistency/{sample.name}.ip.split.01.vs.{sample.name}.ip.split.02.reproducible.peaks.bed'
        if not os.path.exists(bed):
            peak(split_ip_bams, split_input_bams, split_peak_beds, bed,  'consistency', cwd=options.outdir)
        counts.append(count_lines(bed))
    
    try:
        ratio = counts[0] / counts[1]
    except ZeroDivisionError:
        ratio = 0
        logger.error(f'No peaks found in one of the split reproducible peaks, return ratio 0.')
    with open(outputs, 'w') as o:
        o.write(f'{ratio}\n')


def schedule():
    sbatch = """#!/usr/bin/env bash

#SBATCH -n {cpus}                       # Number of cores (-n)
#SBATCH -N 1                        # Ensure that all cores are on one Node (-N)
#SBATCH -t {runtime}                  # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH --mem={memory}G                   # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --job-name={job}            # Short name for the job
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

    if options.scheduler.upper() in ('PBS', 'QSUB'):
        runtime, directive, exe, mail = f'{options.time}:00:00', pbs, 'qsub', pbs_email
        project = '/project/vannostrand'
    elif options.scheduler.upper() in ('SLURM', 'SBATCH'):
        days, hours = divmod(options.time, 24)
        runtime, directive, exe, mail = f'{days}-{hours:02}:00', sbatch, 'sbatch', sbatch_email
        project = '/storage/vannostrand'
    else:
        raise ValueError(f'Unsupported scheduler: {options.scheduler}, see help for supported schedulers.')
    if options.debug:
        root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    else:
        root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    program = os.path.join(root, 'seclip')
    setattr(options, 'ip_fastqs', ' \\\n              '.join(ips))
    setattr(options, 'input_fastqs', ' \\\n                 '.join(inputs))
    setattr(options, 'names', ' '.join(names))
    setattr(options, 'runtime', runtime)
    setattr(options, 'project', project)
    setattr(options, 'scheduler', 'local')
    if options.debug:
        setattr(options, 'debug', '--debug ')
        setattr(options, 'program', os.path.join(root, 'rnaseq', parser.prog))
    else:
        setattr(options, 'debug', '')
        setattr(options, 'program', os.path.join(root, parser.prog))
    arguments = ('ip_fastqs', 'input_fastqs', 'names', 'adapters_fasta', 'barcodes_pattern', 'species',
                 'repeat', 'genome', 'outdir', 'dataset', 'blacklist_bed', 'track', 'track_genome',
                 'l2fc', 'l10p', 'enrichment_filter', 'cpus', 'email', 'scheduler')
    arguments = ' \\\n  '.join([f'--{k} {getattr(options, k)}' for k in arguments if getattr(options, k, '')])
    setattr(options, 'arguments', arguments)

    code = fr"""
export TMPDIR={project}/tmp
export TEMP={project}/tmp
export TMP={project}/tmp

{program} \
  {arguments}
"""
    
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
        print(f'Job {options.job} was successfully submitted with the following resources:')
        data = {'Job name:': options.job, 'Output directory:': options.outdir,
                'Number of cores:': options.cpus, 'Job memory:': options.memory,
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
    if options.scheduler and options.scheduler != 'local':
        schedule()
    else:
        keys = ('outdir', 'genome', 'repeat', 'processes')
        d = vars(options).copy()
        d['processes'] = options.cpus
        setting = '\n'.join([f'{k.title():>22}: {v}' for k, v in d.items() if k in keys])
        logger.debug(HEADER)
        logger.trace(f'\nRunning eCLIP using the following settings:\n\n{setting}\n')
        flow = Flow('eCLIP', description=__doc__.strip())
        flow.run(dry_run=options.dry_run, cpus=options.cpus)
        logger.debug('')
        run_time = str(datetime.timedelta(seconds=int(time.perf_counter() - START_TIME)))
        logger.trace(FOOTER.format(hh_mm_ss=f'time consumed: {run_time}'.upper().center(118)))


if __name__ == '__main__':
    pass
    main()
