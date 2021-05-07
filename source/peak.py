#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
A pipeline for identifying reproducible eCLIP peaks using IDR.
"""

import argparse
import itertools
import math
import os

import pysam
import cmder
import pandas as pd
from seqflow import Flow, task, logger

parser = argparse.ArgumentParser(description=__doc__, prog='eclip')
parser.add_argument('--ip_bams', help='Path to IP bam files (used for calling peaks).', nargs='*')
parser.add_argument('--input_bams', help='Path to INPUT bam files.', nargs='*')
parser.add_argument('--beds', help='Path to peak bed files.', nargs='*')
parser.add_argument('--species', type=str, help="Species name (short name code) the dataset associated with.",
                    default='hg19')
parser.add_argument('--cores', help='Maximum number of CPU cores can be used.', type=int, default=1)
parser.add_argument('--outdir', help='Path to the output directory.')

args = parser.parse_args()
outdir = args.outdir or os.getcwd()
os.chdir(outdir)

data = {}
for ip_bam, input_bam, peak_bed in zip(args.ip_bams, args.input_bams, args.beds):
    assert os.path.isfile(ip_bam), f'IP BAM {ip_bam} is not a file or does not exist.'
    assert os.path.isfile(input_bam), f'INPUT BAM {input_bam} is not a file or does not exist.'
    assert os.path.isfile(peak_bed), f'PEAK BED {peak_bed} is not a file or does not exist.'
    data[os.path.basename(peak_bed.replace('.peak.clusters.bed', ''))] = (ip_bam, input_bam, peak_bed)


def mapped_read_count(bam):
    with pysam.AlignmentFile(bam, 'rb') as sam:
        return sam.mapped


@task(inputs=args.beds, processes=args.cores,
      outputs=lambda i: os.path.basename(i).replace('.peak.clusters.bed', '.peak.clusters.normalized.bed'))
def overlap_peaks(peak_bed, norm_bed):
    ip_bam, input_bam, _ = data[os.path.basename(peak_bed).replace('.peak.clusters.bed', '')]
    cmd = ['overlap_peak.pl', ip_bam, input_bam, peak_bed,
           mapped_read_count(ip_bam), mapped_read_count(input_bam), norm_bed]
    cmder.run(cmd, msg=f'Normalizing peaks in {peak_bed} ...', pmt=True)
    return norm_bed


@task(inputs=overlap_peaks, processes=args.cores,
      outputs=lambda i: i.replace('.normalized.bed', '.normalized.compressed.bed'))
def compress_peaks(norm_bed, compress_bed):
    bed = norm_bed.replace('.bed', '.full.bed')
    cmd = ['compress_peak.pl', bed, compress_bed]
    cmder.run(cmd, msg=f'Compressing peaks in {norm_bed} ...', pmt=True)
    return compress_bed


@task(inputs=compress_peaks, outputs=lambda i: i.replace('.compressed.bed', '.compressed.annotated.full.bed'),
      processes=args.cores)
def annotate_peaks(bed, output):
    bed = bed.replace('.bed', '.full.bed')
    cmd = ['annotate_peak.pl', bed, output, args.species, 'full']
    cmder.run(cmd, msg=f'Annotating peaks in {bed} ...', pmt=True)
    return output


@task(inputs=annotate_peaks, processes=args.cores,
      outputs=lambda i: i.replace('.annotated.full.bed', '.annotated.entropy.bed'))
def calculate_entropy(bed, output):
    logger.info(f'Calculating entropy for {bed} ...')
    name = bed.replace('.peak.clusters.normalized.compressed.annotated.full.bed', '')
    ip_bam, input_bam, peak_bed = data[name]
    ip_mapped_read_count, input_mapped_read_count = mapped_read_count(ip_bam), mapped_read_count(input_bam)
    columns = ['chrom', 'start', 'end', 'peak', 'ip_read_number', 'input_read_number',
               'p', 'v', 'method', 'status', 'l10p', 'l2fc',
               'ensg_overlap', 'feature_type', 'feature_ensg', 'gene', 'region']
    df = pd.read_csv(bed, sep='\t', header=None, names=columns)
    df = df[df.l2fc > 0]
    df['pi'] = df['ip_read_number'] / ip_mapped_read_count
    df['qi'] = df['input_read_number'] / input_mapped_read_count
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


@task(inputs=[], parent=calculate_entropy, kind='create', processes=args.cores,
      outputs=[f'{key1}.vs.{key2}.idr.out.bed' for key1, key2 in itertools.combinations(data.keys(), 2)])
def idr_peaks(inputs, outputs):
    key1, key2 = outputs.replace('.idr.out.bed', '').split('.vs.')
    idr_out = f'{key1}.vs.{key2}.idr.out'
    idr_bed = outputs
    
    peak1 = f'{key1}.peak.clusters.normalized.compressed.annotated.entropy.bed'
    peak2 = f'{key2}.peak.clusters.normalized.compressed.annotated.entropy.bed'
    cmd = ['idr', '--sample', peak1, peak2, '--input-file-type', 'bed', '--rank', '5',
           '--peak-merge-method', 'max', '--plot', '-o', idr_out]
    cmder.run(cmd, msg=f'Running IDR to rank peaks in {peak1} and\n{" " * 40}{peak2} ...', pmt=True)
    
    cmd = ['parse_idr_peaks.pl', idr_out,
           peak1.replace('.bed', '.full.bed'), peak2.replace('.bed', '.full.bed'), idr_bed]
    cmder.run(cmd, msg=f'Parsing IDR peaks in {idr_out} ...', pmt=True)


@task(inputs=[], parent=idr_peaks, kind='create', outputs=f'{".vs.".join(data.keys())}.idr.out.bed',
      processes=args.cores)
def intersect_idr_peaks(inputs, bed):
    tmp_bed = bed.replace('.bed', '.tmp.bed')
    bed1, bed2, bed3 = [f'{key1}.vs.{key2}.idr.out.bed'
                        for key1, key2 in itertools.combinations(data.keys(), 2)]
    cmd = f"bedtools intersect -a {bed1} -b {bed2} > {tmp_bed}"
    cmder.run(cmd)
    cmd = f"bedtools intersect -a {tmp_bed} -b {bed3} > {bed}"
    cmder.run(cmd)
    cmder.run(f'rm {tmp_bed}')
    return bed
    
    
@task(inputs=[], parent=intersect_idr_peaks, outputs=[f'{name}.idr.normalized.bed' for name in data],
      processes=args.cores, kind='create')
def overlap_idr_peaks(idr_out, idr_normalized_bed):
    name = idr_normalized_bed.replace('.idr.normalized.bed', '')
    idr_out = f'{".vs.".join(data.keys())}.idr.out.bed'
    ip_bam, input_bam, peak_bed = data[name]
    cmd = ['overlap_peak.pl', ip_bam, input_bam, idr_out,
           mapped_read_count(ip_bam), mapped_read_count(input_bam), idr_normalized_bed]
    cmder.run(cmd, msg=f'Normalizing IDR peaks of sample {name} ...', pmt=True)
    return idr_normalized_bed
    

@task(inputs=[], kind='create',
      parent=overlap_idr_peaks, outputs=[f'{".vs.".join(data.keys())}.reproducible.peaks.bed'])
def reproducible_peaks(inputs, outputs):
    inputs = [f'{name}.idr.normalized.full.bed' for name in data]
    script = 'reproducible_peaks.pl' if len(data) == 2 else 'reproducible_peaks_3.pl'
    normalized_beds, entropy_full_beds, reproducible_full_beds = [], [], []
    for name in data:
        entropy_full_beds.append(f'{name}.peak.clusters.normalized.compressed.annotated.entropy.full.bed')
        reproducible_full_beds.append(f'{name}.reproducible.peaks.full.bed')
    custom_bed = outputs.replace('.peaks.bed', '.peaks.custom.bed')

    cmd = [script] + inputs + reproducible_full_beds
    cmd += [outputs, custom_bed] + entropy_full_beds
    cmd += [outputs.replace('.reproducible.peaks.bed', '.idr.out')]
    cmder.run(cmd, msg='Identifying reproducible peaks ...', pmt=True)
    

def main():
    flow = Flow('Peak', description=__doc__.strip())
    setattr(args, 'dry_run', False)
    flow.run(dry=args.dry_run, verbose=False, processes=args.cores)
    
    
if __name__ == '__main__':
    main()
