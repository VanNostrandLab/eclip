#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
A pipeline for estimate self-consistency ratio of eCLIP dataset.
"""
import glob
import itertools
import os
import sys
import argparse
import tempfile

import cmder
from seqflow import Flow, task, logger

parser = argparse.ArgumentParser(description=__doc__, prog='eclip')
parser.add_argument('--names', nargs='+', required=True, help='Shortnames for each sample, e.g., rep1, rep2.')
parser.add_argument('--wd', required=True, help='Path to the work directory that contains eCLIP analysis results.')
parser.add_argument('--species', help="Species name (short name code) the dataset associated with, e.g., hg19, mm10.",
                    default='hg19')
parser.add_argument('--l2fc', type=int, help="Only consider peaks at or above this log2 fold change cutoff.", default=3)
parser.add_argument('--l10p', type=int, help="Only consider peaks at or above this log10 p value cutoff.", default=3)
parser.add_argument('--cpus', type=int, help='Maximum number of CPU cores can be used for your job.', default=16)
parser.add_argument('--dry_run', action='store_true',
                    help='Print out steps and files involved in each step without actually running the pipeline.')
options = parser.parse_args()
try:
    os.chdir(options.wd)
except OSError as e:
    logger.error(e)
tmp = 'consistency'
if not os.path.isdir(tmp):
    os.mkdir(tmp)
files = {name: (f'{name}.ip.bam', f'{name}.input.bam') for name in options.names}
key = '.'.join(options.names)

    
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


@task(inputs=[], outputs=[os.path.join(tmp, f'{name}.ip.split.01.bam') for name in files.keys()])
def half_bam(in_bam, out):
    basename = os.path.basename(out).replace('.ip.split.01.bam', '')
    ip_bam, input_bam = files[basename]
    input_split_bam1 = os.path.join(tmp, f'{basename}.input.split.01.bam')
    input_split_bam2 = os.path.join(tmp, f'{basename}.input.split.02.bam')
    split_bam(input_bam, input_split_bam1, input_split_bam2)

    ip_split_bam1 = os.path.join(tmp, f'{basename}.ip.split.01.bam')
    ip_split_bam2 = os.path.join(tmp, f'{basename}.ip.split.02.bam')
    split_bam(ip_bam, ip_split_bam1, ip_split_bam2)


def clipper_peaks(bam, bed=''):
    bed = bed if bed else bam.replace('.ip.bam', '.peak.clusters.bed')
    if os.path.isfile(bed):
        logger.info(f'Clipper bed {bed} already exists.')
    else:
        cmd = f'clipper --species {options.species} --processors {options.cpus} --bam {bam} --outfile {bed}'
        cmder.run(cmd, msg=f'Calling peaks from {bam} using clipper ...', pmt=True)
    return bed
    

@task(inputs=half_bam, outputs=lambda i: i.replace('.bam', '.peak.clusters.bed'))
def pseudo_clipper_peak(bam, bed):
    bam1 = bam.replace('.ip.split.01.bam', '.input.split.01.bam')
    bam2 = bam.replace('.ip.split.01.bam', '.input.split.02.bam')
    bed1 = bam.replace('.ip.split.01.bam', '.input.split.01.peak.clusters.bed')
    bed2 = bam.replace('.ip.split.01.bam', '.input.split.02.peak.clusters.bed')
    clipper_peaks(bam2, bed2)
    clipper_peaks(bam1, bed1)

    bam1, bam2 = bam, bam.replace('.ip.split.01.bam', '.ip.split.02.bam')
    bed1, bed2 = bed, bam.replace('.ip.split.01.bam', '.ip.split.02.peak.clusters.bed')
    clipper_peaks(bam2, bed2)
    clipper_peaks(bam1, bed1)


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


@task(inputs=[], parent=pseudo_clipper_peak,
      outputs=[os.path.join(tmp, f'{".vs.".join([f"{n}.ip.split.{i}" for n in options.names])}.reproducible.peaks.bed')
               for i in ('01', '02')])
def split_reproducible_peak(bed, out):
    for i in ('01', '02'):
        split_ip_bams = [os.path.join(tmp, f'{name}.ip.split.{i}.bam') for name in options.names]
        split_input_bams = [os.path.join(tmp, f'{name}.input.split.{i}.bam') for name in options.names]
        split_peak_beds = [os.path.join(tmp, f'{name}.ip.split.{i}.peak.clusters.bed') for name in options.names]
        peak(split_ip_bams, split_input_bams, split_peak_beds, options.names, out, tmp)


def count_lines(file):
    lines = int(cmder.run(f'wc -l {file}').stdout.read().split()[0])
    logger.info(f'Found {lines:,} lines in {file}.')
    return lines


@task(inputs=[], outputs=['consistency.ratio.txt'], parent=split_reproducible_peak)
def rescue_ratio(bed, txt):
    bed1, bed2 = glob.glob(os.path.join(tmp, '*.reproducible.peaks.bed'))
    count1, count2 = count_lines(bed1), count_lines(bed2)
    try:
        ratio = count1 / count2
    except ZeroDivisionError:
        ratio = 0
        logger.error(f'No peaks found in one of the split reproducible peaks, return ratio 0.')
    with open(txt, 'w') as o:
        o.write(f'{ratio}\n')


if __name__ == '__main__':
    flow = Flow('eCLIP_Consistency', description=__doc__.strip())
    flow.run(dry_run=options.dry_run, cpus=options.cpus)
