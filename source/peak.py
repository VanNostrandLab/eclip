#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Pipeline for using IDR to identify a set of reproducible peaks given eClIP dataset with two or three replicates.
"""

import os
import sys
import math
import argparse
import itertools

import cmder
import inflect
import pandas as pd
from seqflow import Flow, task, logger

parser = argparse.ArgumentParser(description=__doc__, prog='peak')
parser.add_argument('--ip_bams', nargs='+', help='Space separated IP bam files.')
parser.add_argument('--input_bams', nargs='+', help='Space separated INPUT bam files.')
parser.add_argument('--peak_beds', nargs='+', help="Space separated peak bed files.")
parser.add_argument('--outdir', type=str, help="Path to output directory.")
parser.add_argument('--species', type=str, help="Short code for species, e.g., hg19, mm10.")
parser.add_argument('--l2fc', type=float, help="Only consider peaks at or above this l2fc cutoff.", default=3)
parser.add_argument('--l10p', type=float, help="Only consider peaks at or above this l10p cutoff.", default=3)
parser.add_argument('--idr', type=float, help="Only consider peaks at or above this idr score cutoff.", default=0.01)
parser.add_argument('--dry_run', action='store_true',
                    help='Print out steps and inputs/outputs of each step without actually running the pipeline.')
parseradd_argument('--debug', action='store_true', help='Invoke debug mode (only for develop purpose).')


def validate_paths():
    def files_exist(files, tag):
        engine, paths = inflect.engine(), []
        for i, file in enumerate(files, start=1):
            if os.path.exists(file):
                if not os.path.isfile(file):
                    logger.error(f'The {engine.ordinal(i)} file in {tag} "{file}" is not a file.')
                    sys.exit(1)
                else:
                    paths.append(file)
            else:
                logger.error(f'The {engine.ordinal(i)} file in {tag} "{file}" does not exist.')
        return paths

    args = parser.parse_args()
    ip_bams = files_exist(args.ip_bams, 'IP bams')
    input_bams = files_exist(args.input_bams, 'INPUT bams')
    peak_beds = files_exist(args.peak_beds, 'Peak beds')
    outdir = args.outdir or os.getcwd()
    if os.path.exists(outdir):
        if not os.path.isdir(outdir):
            logger.error(f'Outdir "{outdir}" is a file not a directory.')
            sys.exit(1)
    else:
        logger.error(f'Outdir "{outdir}" does not exist.')
        os.mkdir(outdir)

    files = {}
    if len(ip_bams) == len(input_bams) == len(peak_beds):
        for ip_bam, input_bam, peak_bed in zip(ip_bams, input_bams, peak_beds):
            key = os.path.basename(peak_bed).replace('.peak.clusters.bed', '')
            files[key] = (ip_bam, input_bam, peak_bed,
                          f'{key}.peak.clusters.normalized.compressed.annotated.entropy.bed')
    else:
        logger.error('Unequal number of files provided!')
        sys.exit()
    return files, outdir, args


files, outdir, options = validate_paths()
os.chdir(outdir)


@task(inputs=options.ip_bams + options.input_bams, outputs=lambda i: i.replace('.bam', '.mapped.reads.count.txt'))
def count_mapped_reads(bam, txt):
    cmd = f'samtools view -c -F 0x4 {bam} > {txt}'
    cmder.run(cmd, msg=f'Count mapped reads in {bam} ...', pmt=True)


def get_mapped_reads(bam):
    with open(bam.replace('.bam', '.mapped.reads.count.txt')) as f:
        return int(f.read().strip())


@task(inputs=[v[2] for v in files.values()], outputs=lambda i: i.replace('.bed', '.normalized.bed'),
      parent=count_mapped_reads)
def normalize_peak(bed, normalized_bed):
    ip_bam, input_bam, peak_bed, _ = files[bed.replace('.peak.clusters.bed', '')]
    ip_read_count, input_read_count = get_mapped_reads(ip_bam), get_mapped_reads(input_bam)
    cmd = ['overlap_peak.pl', ip_bam, input_bam, peak_bed, ip_read_count, input_read_count, normalized_bed]
    cmder.run(cmd, msg=f'Normalizing peaks in {peak_bed} ...', pmt=True)
    return normalized_bed


@task(inputs=normalize_peak, outputs=lambda i: i.replace('.bed', '.compressed.bed'))
def compress_peak(normalized_bed, compressed_bed):
    cmd = ['compress_peak.pl', normalized_bed.replace('.bed', '.tsv'), compressed_bed]
    cmder.run(cmd, msg=f'Compressing peaks in {normalized_bed} ...', pmt=True)
    return compressed_bed


@task(inputs=compress_peak, outputs=lambda i: i.replace('.bed', '.annotated.bed'))
def annotate_peak(compressed_bed, annotated_bed):
    cmd = ['annotate_peak.pl', compressed_bed.replace('.bed', '.tsv'), annotated_bed, options.species, 'full']
    cmder.run(cmd, msg=f'Annotating peaks in {compressed_bed} ...', pmt=True)
    return annotated_bed


def calculate_entropy(bed, output, ip_read_count, input_read_count):
    logger.info(f'Calculating entropy for {bed} ...')
    columns = ['chrom', 'start', 'end', 'peak', 'ip_read_number', 'input_read_number',
               'p', 'v', 'method', 'status', 'l10p', 'l2fc',
               'ensg_overlap', 'feature_type', 'feature_ensg', 'gene', 'region']
    df = pd.read_csv(bed, sep='\t', header=None, names=columns)
    df = df[df.l2fc >= 0]
    # df = df[(df.l2fc >= options.l2fc) & (df.l10p >= options.l10p)]
    if df.empty:
        logger.error(f'No valid peaks found in {bed} (l2fc > 0 failed).')
        sys.exit(1)
    df['pi'] = df['ip_read_number'] / ip_read_count
    df['qi'] = df['input_read_number'] / input_read_count

    df['entropy'] = df.apply(lambda row: 0 if row.pi <= row.qi else row.pi * math.log2(row.pi / row.qi), axis=1)
    df['excess_reads'] = df['pi'] - df['qi']
    entropy = output.replace('.entropy.bed', '.entropy.tsv')
    df.to_csv(entropy, index=False, columns=columns + ['entropy'], sep='\t', header=False)

    excess_read = output.replace('.bed', '.excess.reads.tsv')
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


@task(inputs=annotate_peak, outputs=lambda i: i.replace('.bed', '.entropy.bed'))
def entropy_peak(annotated_bed, entropy_bed):
    basename = annotated_bed.replace('.peak.clusters.normalized.compressed.annotated.bed', '')
    ip_bam, input_bam, peak_bed, _ = files[basename]
    ip_read_count, input_read_count = get_mapped_reads(ip_bam), get_mapped_reads(input_bam)
    calculate_entropy(annotated_bed, entropy_bed, ip_read_count, input_read_count)
    return entropy_bed


@task(inputs=[], parent=entropy_peak,
      outputs=[f'{key1}.vs.{key2}.idr.out' for key1, key2 in itertools.combinations(files.keys(), 2)])
def run_idr(bed, out):
    key1, key2 = out.replace('.idr.out', '').split('.vs.')
    entropy_bed1, entropy_bed2 = files[key1][3], files[key2][3]
    cmd = ['idr', '--sample', entropy_bed1, entropy_bed2, '--input-file-type', 'bed', '--rank', '5',
           '--peak-merge-method', 'max', '--plot', '-o', out]
    cmder.run(cmd, msg=f'Running IDR to rank peaks in {entropy_bed1} and\n{" " * 40}{entropy_bed2} ...',
              pmt=True)


@task(inputs=[], parent=run_idr,
      outputs=[f'{key1}.vs.{key2}.idr.out.bed' for key1, key2 in itertools.combinations(files.keys(), 2)])
def parse_idr(out, bed):
    key1, key2 = bed.replace('.idr.out.bed', '').split('.vs.')
    idr_out, idr_bed = f'{key1}.vs.{key2}.idr.out', f'{key1}.vs.{key2}.idr.out.bed'
    if len(files) == 2:
        entropy_bed1, entropy_bed2 = files[key1][3], files[key2][3]
        cmd = ['parse_idr_peaks.pl', idr_out,
               entropy_bed1.replace('.bed', '.tsv'), entropy_bed2.replace('.bed', '.tsv'), idr_bed]
        cmder.run(cmd, msg=f'Parsing IDR peaks in {idr_out} ...', pmt=True)
    else:
        idr_cutoffs = {0.001: 1000, 0.005: 955, 0.01: 830, 0.02: 705, 0.03: 632, 0.04: 580, 0.05: 540,
                       0.06: 507, 0.07: 479, 0.08: 455, 0.09: 434,
                       0.1: 415, 0.2: 290, 0.3: 217, 0.4: 165, 0.5: 125, 1: 0}
        with open(idr_out) as f, open(idr_bed, 'w') as o:
            for line in f:
                fields = line.strip().split('\t')
                chrom, start, stop, _, idr_score, strand = fields[:6]
                if float(idr_score) >= idr_cutoffs[options.idr]:
                    o.write(f'{chrom}\t{start}\t{stop}\t.\t.\t{strand}\n')
                        

@task(inputs=[], outputs=f'{".vs.".join(files.keys())}.idr.out.bed', parent=parse_idr)
def intersect_idr(bed, intersected_bed):
    if len(files) == 2:
        idr_out, idr_bed = f'{".vs.".join(files.keys())}.idr.out', f'{".vs.".join(files.keys())}.idr.out.bed'
        idr_intersected_bed = f'{".vs.".join(files.keys())}.idr.intersected.bed'
        cmder.run(f'cp {idr_out} {idr_intersected_bed}')
    elif len(files) == 3:
        idr_intersected_bed = f'{".vs.".join(files.keys())}.idr.intersected.bed'
        idr_bed = f'{".vs.".join(files.keys())}.idr.out.bed'

        bed1, bed2, bed3 = [f'{key1}.vs.{key2}.idr.out.bed'
                            for key1, key2 in itertools.combinations(files.keys(), 2)]
        tmp_bed = idr_intersected_bed.replace('.bed', '.tmp.bed')
        cmder.run(f'bedtools intersect -a {bed1} -b {bed2} > {tmp_bed}', msg='Intersecting IDR beds ...')
        cmder.run(f'bedtools intersect -a {tmp_bed} -b {bed3} > {idr_intersected_bed}', msg='Intersecting IDR beds ...')
        cmder.run(f'rm {tmp_bed}')
        
        entropy_beds = [f'{key}.peak.clusters.normalized.compressed.annotated.entropy.tsv' for key in files]
        cmd = ['parse_idr_peaks_3.pl', idr_intersected_bed] + entropy_beds + [f'{idr_bed}']
        cmder.run(cmd, msg=f'Parsing intersected IDR peaks in {idr_bed} ...', pmt=True)


@task(inputs=[], outputs=[f'{key}.idr.normalized.bed' for key in files], parent=intersect_idr)
def normalize_idr(bed, idr_normalized_bed):
    # idr_intersected_bed = f'{".vs.".join(files.keys())}.idr.intersected.bed'
    idr_bed = f'{".vs.".join(files.keys())}.idr.out.bed'
    key = idr_normalized_bed.replace('.idr.normalized.bed', '')
    ip_bam, input_bam, peak_bed, _ = files[key]

    cmd = ['overlap_peak.pl', ip_bam, input_bam, idr_bed,
           get_mapped_reads(ip_bam), get_mapped_reads(input_bam), idr_normalized_bed]
    cmder.run(cmd, msg=f'Normalizing IDR peaks for sample {key} ...', pmt=True)
        

@task(inputs=[], outputs=f'{".vs.".join([key for key in files])}.reproducible.peaks.bed', parent=normalize_idr)
def reproducible_peak(inputs, reproducible_bed):
    script = f'reproducible_peaks_{len(files)}.pl'
    custom = reproducible_bed.replace('.peaks.bed', '.peaks.custom.tsv')
    idr_normalized_full_beds, entropy_full_beds, reproducible_txts = [], [], []
    for ip_bam, input_bam, peak_bed in zip(options.ip_bams, options.input_bams, options.peak_beds):
        name = os.path.basename(peak_bed.replace('.peak.clusters.bed', ''))
        idr_normalized_full_beds.append(f'{name}.idr.normalized.tsv')
        entropy_full_beds.append(f'{name}.peak.clusters.normalized.compressed.annotated.entropy.tsv')
        reproducible_txts.append(f'{name}.reproducible.peaks.tsv')

    cmd = [script] + idr_normalized_full_beds + reproducible_txts
    cmd += [reproducible_bed, custom] + entropy_full_beds
    cmd += [f'{".vs.".join(files.keys())}.idr.out{".bed" if len(files) == 3 else ""}']
    cmder.run(cmd, msg='Identifying reproducible peaks ...', pmt=True)

    
if __name__ == '__main__':
    flow = Flow('Peak', description=__doc__.strip())
    flow.run(dry=options.dry_run)
