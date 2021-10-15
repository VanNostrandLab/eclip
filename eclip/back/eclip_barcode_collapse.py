#!/usr/bin/env python

import os.path

import fire
import cmder
import pysam


def barcode_collapse(bam, output, debug=False):
    """
    Deduplicate paired-end BAM by collapsing barcodes.

    :param bam: str, path to BAM file.
    :param output: str, path to the output file.
    :param debug: bool, set to True for invoking debug mode.
    """

    it.info(f'Deduplicating {bam} by collapsing barcodes ...')
    pysam.set_verbosity(1 if debug else 0)
    with pysam.AlignmentFile(bam, 'rb') as b1, pysam.AlignmentFile(bam, 'rb') as b2:
        results = {}
        for read1, read2 in zip(itertools.islice(b1, 0, None, 2), itertools.islice(b2, 1, None, 2)):
            if read1.query_name != read2.query_name:
                it.error_and_exit(f'Read names do not match: {read1.query_name} != {read2.query_name}.')
            if read1.is_unmapped or read2.is_unmapped or read1.reference_name != read2.reference_name:
                continue
            if not read1.is_read1:
                read1, read2 = read2, read1
            randomer = read1.query_name.split(':')[0]
            start = read1.positions[-1] if read1.is_reverse else read1.pos
            stop = read2.positions[-1] if read2.is_reverse else read2.pos
            strand = '-' if read1.is_reverse else '+'
            location = (read1.reference_name, start, stop, strand, randomer)
            if location in results:
                continue
            results[location] = (read1, read2)
        with pysam.AlignmentFile(out, 'wb', template=b1) as o:
            for (read1, read2) in results.values():
                o.write(read1)
                o.write(read2)
        it.info(f'Deduplicating {bam} by collapsing barcodes complete.')


if __name__ == '__main__':
    fire.Fire(barcode_collapse)
