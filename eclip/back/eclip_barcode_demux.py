#!/usr/bin/env python

import gzip
import os

import fire
import iToolBox as it
from Bio.SeqIO.QualityIO import FastqGeneralIterator


def demux(fastq1, fastq2, barcode1, barcode2, outdir, basename, randomer_length=10, allowed_mismatch=1):
    """
    Demultiplex eCLIP paired-end reads (2 or 4 files will be generated depends on barcode1 and barcode2).

    :param fastq1: str, path to read 1 (gzip compressed) FASTQ file ends with .fastq.gz extension.
    :param fastq2: str, path to read 2 (gzip compressed) FASTQ file ends with .fastq.gz extension.
    :param barcode1: str, short name or code for barcode 1, e.g., NIL, A01, B06, ... .
    :param barcode2: str, short name or code for barcode 2, e.g., NIL, A01, B06, ... .
    :param outdir: str, path to the output directory.
    :param basename: str, basename for the dataset.
    :param randomer_length: int, length of randomer.
    :param allowed_mismatch: int, number of allowed mismatches for barcode.
    """

    def hamming(key, barcode, seq, allow_mismatch):
        mismatch = len(barcode) - sum(x == y or x == 'N' or y == 'N' for x, y in zip(barcode, seq))
        return (key, len(barcode), mismatch) if mismatch <= allow_mismatch else None

    it.mkdir(outdir)
    it.info(f'Demultiplexing {fastq1} and {fastq2} with barcodes {barcode1} and {barcode2} ...')
    barcodes_dict = {'A01': 'AAGCAAT',
                     'A03': 'ATGACCNNNNT',
                     'A04': 'CAGCTTNNNNT',
                     'B06': 'GGCTTGT',
                     'C01': 'ACAAGTT',
                     'D8f': 'TGGTCCT',
                     'F05': 'GGATACNNNNT',
                     'G07': 'TCCTGTNNNNT',
                     'X1A': 'NNNNNCCTATAT',
                     'X1B': 'NNNNNTGCTATT',
                     'X2A': 'NNNNNTATACTT',
                     'X2B': 'NNNNNATCTTCT'}
    max_barcode_length = max(len(barcode) for barcode in barcodes_dict.values())
    writers = {}
    for barcode in {barcode1, barcode2}:
        writers[barcode] = (gzip.open(os.path.join(outdir, f'{basename}.{barcode}.r1.fastq.gz'), 'wt'),
                            gzip.open(os.path.join(outdir, f'{basename}.{barcode}.r2.fastq.gz'), 'wt'))

    with gzip.open(fastq1, 'rt') as f1, gzip.open(fastq2, 'rt') as f2:
        for i, (read1, read2) in enumerate(zip(FastqGeneralIterator(f1), FastqGeneralIterator(f2))):
            (name1, seq1, quality1), (name2, seq2, quality2) = read1, read2
            n1, n2 = name1.split()[0], name2.split()[0]
            if n1 != n2:
                it.error_and_exit(f'ValueError: Paired-End reads have mismatch names: {name1} != {name2}')

            matches = (hamming(key, barcode, seq1[:max_barcode_length], allowed_mismatch)
                       for key, barcode in barcodes_dict.items())
            matches = [match for match in matches if match]
            if matches:
                barcode, barcode_length, _ = sorted(matches, key=lambda x: x[2])[0]
                r1 = f'@{seq2[:randomer_length]}:{name1}\n{seq1[barcode_length:]}\n+\n{quality1[barcode_length:]}\n'
            else:
                barcode = 'NIL'
                r1 = f'@{seq2[:randomer_length]}:{name1}\n{seq1}\n+\n{quality1}\n'
            r2 = f'@{seq2[:randomer_length]}:{name2}\n{seq2[randomer_length:]}\n+\n{quality2[randomer_length:]}\n'

            if barcode in writers:
                writer1, writer2 = writers[barcode]
                writer1.write(r1)
                writer2.write(r2)
    _ = [[v[0].close(), v[1].close()] for v in writers.values()]
    it.info(f'Demultiplexing {fastq1} and {fastq2} with barcodes {barcode1} and {barcode2} complete.')


if __name__ == '__main__':
    fire.Fire(demux)
