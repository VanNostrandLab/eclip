#!/usr/bin/env python

import sys

import fire
import cmder


def umi(fastq, barcodes_pattern, output, log, debug=False):
    """
    Extract UMIs for single-end read.

    :param fastq: str, path to (gzip compressed) FASTQ file ends with .fastq.gz extension.
    :param barcodes_pattern: str, barcode pattern.
    :param output: str, path to the output file.
    :param log: str, path to the log file.
    :param debug: bool, set to True for invoking debug mode.
    """

    cmd = ['umi_tools', 'extract',
           '--random-seed', 1,
           '--stdin', fastq,
           '--bc-pattern', barcodes_pattern,
           '--log', log,
           '--stdout', output]
    message = f'Extract UMIs for single-end read {fastq} ...'
    if debug:
        cmder.run(cmd, msg=message, pmt=True, stdout=sys.stdout, stderr=sys.stderr)
    else:
        cmder.run(cmd, msg=message, pmt=True)


if __name__ == '__main__':
    fire.Fire(umi)
