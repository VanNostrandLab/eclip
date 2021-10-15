#!/usr/bin/env python

import fire
import cmder


def umi_dedup(bam, output, debug=False):
    """
    Deduplicate single-end BAM by umi-tools dedup.

    :param bam: str, path to BAM file.
    :param output: str, path to the output file.
    :param debug: bool, set to True for invoking debug mode.
    """

    cmd = ['umi_tools', 'dedup', '--random-seed', 1, '--stdin', bam, '--method', 'unique', '--stdout', output]
    cmder.run(cmd, msg=f'Deduplicating {bam} by umi_tools dedup ...', pmt=True, debug=debug)


if __name__ == '__main__':
    fire.Fire(umi_dedup)
