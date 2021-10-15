#!/usr/bin/env python
import os.path
import shutil
import sys
import tempfile

import fire
import cmder


def fastq_sort(fastq, output, debug=False):
    """
    Sort fastq file by id and save results to a gzip compressed output file.

    :param fastq: str, path to FASTQ file.
    :param output: str, path to the output file (it will always generates gzip compressed content).
    :param debug: bool, set to True for invoking debug mode.
    """

    outdir = os.path.dirname(output) or os.getcwd()
    tmp = tempfile.mkdtemp(suffix='_sort', prefix='fastq_', dir=outdir)
    if fastq.endswith('.gz'):
        cmd = f'zcat {fastq} | fastq-sort --id --temporary-directory {tmp} -S 2G {fastq} | gzip > {output}'
    else:
        cmd = f'fastq-sort --id --temporary-directory {tmp} -S 2G {fastq} | gzip > {output}'
    try:
        cmder.run(cmd, msg=f'Sorting {fastq} ...', pmt=True)
    finally:
        shutil.rmtree(tmp)


if __name__ == '__main__':
    fire.Fire(fastq_sort)
