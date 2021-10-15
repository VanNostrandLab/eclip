#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
A pipeline for estimate rescue ratio of eCLIP dataset.
"""

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
tmp = 'rescue'
if not os.path.isdir(tmp):
    os.mkdir(tmp)
ip_bams = [f'{name}.ip.bam' for name in options.names]
input_bams = [f'{name}.input.bam' for name in options.names]
files = {}
for name1, name2 in itertools.combinations(options.names, 2):
    files[f'{name1}.{name2}'] = (f'{name1}.ip.bam', f'{name2}.ip.bam', f'{name1}.input.bam', f'{name2}.input.bam')
key = '.'.join(options.names)


def merge_bam(bam1, bam2, bam):
    if os.path.isfile(bam):
        logger.info(f'BAM file {bam} already exist.')
    else:
        cmd = f'samtools merge {bam} {bam1} {bam2}'
        cmder.run(cmd, msg=f'Merging {bam1} and {bam2} ...')
    return bam
    
    
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


@task(inputs=[], outputs=[os.path.join(tmp, 'ip.pseudo.01.bam'), os.path.join(tmp, 'input.pseudo.01.bam')])
def pseudo_bam(in_bam, out):
    basename = os.path.basename(out).replace('.ip.pseudo.01.bam', '')
    ip_bam1, ip_bam2, input_bam1, input_bam2 = files[basename]
    input_pseudo_bam = merge_bam(input_bam1, input_bam2, os.path.join(tmp, f'{basename}.input.pseudo.bam'))
    input_pseudo_bam1 = os.path.join(tmp, f'{basename}.input.pseudo.01.bam')
    input_pseudo_bam2 = os.path.join(tmp, f'{basename}.input.pseudo.02.bam')
    split_bam(input_pseudo_bam, input_pseudo_bam1, input_pseudo_bam2)
    
    ip_pseudo_bam = merge_bam(ip_bam1, ip_bam2, os.path.join(tmp, f'{basename}.ip.pseudo.bam'))
    ip_pseudo_bam1 = os.path.join(tmp, f'{basename}.ip.pseudo.01.bam')
    ip_pseudo_bam2 = os.path.join(tmp, f'{basename}.ip.pseudo.02.bam')
    split_bam(ip_pseudo_bam, ip_pseudo_bam1, ip_pseudo_bam2)


def clipper_peaks(bam, bed=''):
    bed = bed if bed else bam.replace('.ip.bam', '.peak.clusters.bed')
    if os.path.isfile(bed):
        logger.info(f'Clipper bed {bed} already exists.')
    else:
        cmd = f'clipper --species {options.species} --processors {options.cpus} --bam {bam} --outfile {bed}'
        cmder.run(cmd, msg=f'Calling peaks from {bam} using clipper ...', pmt=True)
    return bed
    

@task(inputs=pseudo_bam, outputs=lambda i: i.replace('.bam', '.peak.clusters.bed'))
def pseudo_clipper_peak(bam, bed):
    bam1, bam2 = bam, bam.replace('.ip.pseudo.01.bam', '.ip.pseudo.02.bam')
    bed1, bed2 = bed, bam.replace('.ip.pseudo.01.bam', '.ip.pseudo.02.peak.clusters.bed')
    clipper_peaks(bam2, bed2)
    clipper_peaks(bam1, bed1)


def peak(ip_bams, input_bams, peak_beds, reproducible_bed, outdir, cwd):
    cmd = ['peak', '--ip_bams', ' '.join(ip_bams),
           '--input_bam', ' '.join(input_bams),
           '--peak_beds', ' '.join(peak_beds),
           '--read_type', 'SE',
           '--species', 'hg19' if options.species in ('hg19', 'hg19chr19') else options.species,
           '--outdir', outdir, '--cores', options.cpus,
           '--l2fc', options.l2fc, '--l10p', options.l10p]
    if ids:
        cmd.extend(['--ids', ' '.join(ids)])
    cwd = cwd if cwd else os.path.dirname(reproducible_bed)
    cmder.run(cmd, cwd=cwd, stdout=sys.stdout, stderr=sys.stderr)
    return reproducible_bed


@task(inputs=[], parent=pseudo_clipper_peak,
      outputs=[os.path.join(tmp, f'{key}.ip.pseudo.01.vs.{key}.ip.pseudo.02.reproducible.peaks.bed')])
def pseudo_reproducible_peak(bed, out):
    pseudo_ip_bams, pseudo_input_bams, pseudo_peak_beds, ids = [], [], [], []
    for name in files.keys():
        pseudo_ip_bams.extend([os.path.join(tmp, f'{name}.ip.pseudo.01.bam'),
                               os.path.join(tmp, f'{name}.ip.pseudo.02.bam')])
        pseudo_input_bams.extend([os.path.join(tmp, f'{name}.input.pseudo.01.bam'),
                                  os.path.join(tmp, f'{name}.input.pseudo.02.bam')])
        pseudo_peak_beds.extend([os.path.join(tmp, f'{name}.ip.pseudo.01.peak.clusters.bed'),
                                 os.path.join(tmp, f'{name}.ip.pseudo.02.peak.clusters.bed')])
        ids.extend([f'{name}.pseudo.01', f'{name}.pseudo.02'])
    peak(pseudo_ip_bams, pseudo_input_bams, pseudo_peak_beds, out, tmp, options.wd)


def count_lines(file):
    lines = int(cmder.run(f'wc -l {file}').stdout.read().split()[0])
    logger.info(f'Found {lines:,} lines in {file}.')
    return lines


@task(inputs=pseudo_reproducible_peak, outputs=['rescue.ratio.txt'])
def rescue_ratio(bed, txt):
    pseudo_count = count_lines(bed)
    actual_count = count_lines(os.path.join(tmp, f'{key}.ip.pseudo.01.vs.{key}.ip.pseudo.02.reproducible.peaks.bed'))
    try:
        ratio = max(actual_count, pseudo_count) / min(actual_count, pseudo_count)
    except ZeroDivisionError:
        ratio = 0
        logger.error(f'No peaks found in reproducible peaks or pseudo reproducible peaks, return ratio 0.')
    with open(txt, 'w') as o:
        o.write(f'{ratio}\n')


if __name__ == '__main__':
    flow = Flow('eCLIP_Rescue', description=__doc__.strip())
    flow.run(dry_run=options.dry_run, cpus=options.cpus)
