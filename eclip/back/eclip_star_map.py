#!/usr/bin/env python

import os
import shutil
import sys
import tempfile

import fire
import cmder
import iToolBox as it


def starmap(fastq1, fastq2, index, mnm, bam, cpus=1, debug=False):
    """
    Map reads to reference genome or repeat elements using STAR.

    :param fastq1: str, path to a FASTQ file (read 1).
    :param fastq2: str, path to a FASTQ file (read 2).
                 For single-end dataset, using a non-existing file path or special string
                 'none' or 'None' for fastq2 to avoid required argument error.
    :param index: str, path to a STAR genome or repeat elements index directory.
    :param mnm: int, maximum number of loci the read is allowed to map to.
    :param bam: str, path to the output BAM file, must ends with .bam extension.
    :param cpus: int, the number of CPUs can be used.
    :param debug: bool, set to True to invoke debug mode (only for development purpose).
    """

    outdir = os.path.dirname(fastq1) or os.getcwd()
    folder = tempfile.mkdtemp(prefix=os.path.basename(fastq1), suffix='.star.map', dir=outdir)
    cmd = ['STAR',
           '--runMode', 'alignReads',
           '--runThreadN', cpus,
           '--alignEndsType', 'EndToEnd',
           '--genomeDir', index,
           '--outBAMcompression', 10,
           '--outFileNamePrefix', f'{folder}/',
           '--outFilterMultimapNmax', mnm,
           '--outFilterScoreMin', 10,
           '--outReadsUnmapped', 'Fastx',
           '--outSAMattrRGline', 'ID:foo',
           '--outSAMattributes', 'All',
           '--outSAMtype', 'BAM', 'Unsorted',
           '--readFilesCommand', 'zcat' if fastq1.endswith('.gz') else '-',
           '--readFilesIn', fastq1]
    if os.path.isfile(fastq2):
        cmd.append(fastq2)
        message = f'Map paired reads {fastq1} and\n{28 * " "}{fastq2} to {label} using STAR ...'
    else:
        message = f'Map single read {fastq1} to {label} using STAR ...'
    cmder.run(cmd, msg=message, pmt=True, debug=debug)

    move = shutil.copy if debug else shutil.move
    move(os.path.join(folder, 'Aligned.out.bam'), '{name}.bam'.format(name=basename))
    move(os.path.join(folder, 'Unmapped.out.mate1'), '{name}.unmap.mate1'.format(name=basename))
    move(os.path.join(folder, 'Log.final.out'), '{name}.log.final.out'.format(name=basename))
    mate2 = os.path.join(folder, 'Unmapped.out.mate2')
    if os.path.isfile(mate2):
        move(mate2, '{name}.unmap.mate2'.format(name=basename))
    if not debug:
        shutil.rmtree(folder)


if __name__ == "__main__":
    fire.Fire(starmap)

