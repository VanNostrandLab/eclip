#!/usr/bin/env python

import gzip
import os

import fire
import cmder
import iToolBox as it
import pysam


def b2bw(bam, basename, strand, lib_type):
    """
    Make positive and negative BigWig files from a BAM file.

    :param bam: str, path to coordinate sorted and indexed bam file.
    :param basename: str, path to the full basename of the BigWig files.
    :param strand: str, strand of read 1, either + or -.
    :param lib_type: str, library type, either single or paired.
    """

    def bam_to_bigwig(bam, scale, strand, bw, genome_length):
        bg, bg_sort = bw.replace('.bw', '.bg'), bw.replace('.bw', '.sort.bg')
        cmd = f'genomeCoverageBed -ibam {bam} -bg -scale {scale} -strand {strand} -du -split > {bg}'
        cmder.run(cmd)
        cmd = f'bedSort {bg} {bg_sort}'
        cmder.run(cmd)
        cmd = f'bedGraphToBigWig {bg_sort} {genome_length} {bw}'
        cmder.run(cmd)
        cmder.run(f'rm {bg} {bg_sort}')

    it.info(f'Making BigWig files from {bam} ...')
    pos_bw, neg_bw = f'{basename}.pos.bw', f'{basename}.neg.bw'

    with pysam.AlignmentFile(bam, 'rb') as sam:
        refs, lens = sam.references, sam.lengths
        total_reads = sam.mapped

    length = f'{basename}.reference.genome.chromosome.name.length.tsv'
    if not os.path.isfile(length):
        with open(length, 'w') as o:
            o.writelines(f'{r}\t{l}\n' for r, l in zip(refs, lens))
    total_reads = total_reads / 2 if lib_type == 'paired' else total_reads
    try:
        scale = 1000000.0 / total_reads
    except ZeroDivisionError:
        scale = 0
        it.error_and_exit(f'No reads was found in BAM {bam}, empty BigWig file was created.')
    if strand == '+':
        bam_to_bigwig(bam, scale, '+', pos_bw, length)
        bam_to_bigwig(bam, -1 * scale, '-', neg_bw, length)
    else:
        bam_to_bigwig(bam, -1 * scale, '-', pos_bw, length)
        bam_to_bigwig(bam, scale, '+', neg_bw, length)


if __name__ == '__main__':
    fire.Fire(b2bw)
