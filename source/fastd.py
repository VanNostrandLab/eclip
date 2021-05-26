#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
A fast tool for eCLIP paired-end FASTQ file demuliplexing.
"""

import sys
import gzip
import argparse
import re
# from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
# from functools import partial

from Bio.SeqIO.QualityIO import FastqGeneralIterator

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--fastqs', nargs='+', help='Path to the first read (r1) and second read (r2) fastq files.')
parser.add_argument('--barcodes', nargs='+', help='Names of the barcodes.')
parser.add_argument('--basename', type=str, help='Basename of the output file.')
parser.add_argument('--cores', type=int, help='Number of CPU cores can be used.', default=8)
parser.add_argument('--randomer_length', type=int, help='Length of the randomer, default: 5.', default=5)
parser.add_argument('--allow_mismatch', type=int, help='Maximum number of base allowed for barcode mismatching.',
                    default=1)
parser.add_argument('--gz', help='Output gzipped output file.', action='store_true')
parser.add_argument('--verbose', help='Print out process details.', action='store_true')

args = parser.parse_args()

try:
    from loguru import logger
    logger.remove()
    logger.add(sys.stderr, level="DEBUG" if args.verbose else "INFO",
               format="<light-green>[{time:HH:mm:ss}]</light-green> <level>{message}</level>")
except ImportError:
    import logging
    logging.basicConfig(format='[%(asctime)s] %(message)s', datefmt='%H:%M:%S',
                        level=logging.DEBUG if args.verbose else logging.INFO)
    logger = logging.getLogger()

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
barcodes_pattern = {'A01': 'AAGCAAT',
                    'A03': 'ATGACC[ATGCN]{4}T',
                    'A04': 'CAGCTT[ATGCN]{4}T',
                    'B06': 'GGCTTGT',
                    'C01': 'ACAAGTT',
                    'D8f': 'TGGTCCT',
                    'F05': 'GGATAC[ATGCN]{4}T',
                    'G07': 'TCCTGT[ATGCN]{4}T',
                    'X1A': '[ATGCN]{5}CCTATAT',
                    'X1B': '[ATGCN]{5}TGCTATT',
                    'X2A': '[ATGCN]{5}TATACTT',
                    'X2B': '[ATGCN]{5}ATCTTCT'}
barcodes_length = {k: len(v) for k, v in barcodes_dict.items()}


def opener(filename, mode):
    if mode == 'w':
        handler, mode = (gzip.open, 'wt') if filename.endswith('.gz') else (open, 'w')
    else:
        handler, mode = (gzip.open, 'rt') if filename.endswith('.gz') else (open, 'r')
    return handler(filename, mode)


def hamming(key, barcode, seq, allow_mismatch=0):
    mismatch = len(barcode) - sum(x == y or x == 'N' or y == 'N' for x, y in zip(barcode, seq))
    return (key, len(barcode), mismatch) if mismatch <= allow_mismatch else None


def barcode_count(read, pattern):
    pass


def barcoding(read, max_barcode_length=0, randomer_length=0, allow_mismatch=0):
    (name1, seq1, quality1), (name2, seq2, quality2) = read
    n1, n2 = name1.split()[0], name2.split()[0]
    assert n1 == n2, ValueError(f'Paired-End reads have mismatch names: {n1} != {n2}')

    matches = (hamming(key, barcode, seq1[:max_barcode_length], allow_mismatch=allow_mismatch)
               for key, barcode in barcodes_dict.items())
    matches = [match for match in matches if match]
    if matches:
        barcode, barcode_length, _ = sorted(matches, key=lambda x: x[2])[0]
    else:
        barcode, barcode_length = 'NIL', randomer_length
    r1 = f'@{seq2[:randomer_length]}:{name1}\n{seq1[barcode_length:]}\n+\n{quality1[barcode_length:]}\n'
    r2 = f'@{seq2[:randomer_length]}:{name2}\n{seq2[randomer_length:]}\n+\n{quality2[randomer_length:]}\n'
    return barcode, r1, r2


def fastd():
    (fastq1, fastq2), barcodes = args.fastqs, args.barcodes
    logger.debug(f'Demultiplexing {fastq1} and {fastq2} with barcodes {" and ".join(args.barcodes)} ...')
    if barcodes[0] == 'NIL':
        barcodes = barcodes_dict
        pattern = re.compile(fr'{"|".join([f"(?P<{k}>^{v})" for k, v in barcodes_pattern.items()])}')
        writers = {'NIL': (opener(f'{args.basename}.NIL.r1.fastq{".gz" if args.gz else ""}', 'w'),
                           opener(f'{args.basename}.NIL.r2.fastq{".gz" if args.gz else ""}', 'w'))}
    else:
        barcodes = {k: v for k, v in barcodes_dict.items() if k in barcodes}
        pattern = re.compile(fr'{"|".join([f"(?P<{k}>^{v})" for k, v in barcodes_pattern.items() if k in barcodes])}')
        writers = {barcode: (opener(f'{args.basename}.{barcode}.r1.fastq{".gz" if args.gz else ""}', 'w'),
                             opener(f'{args.basename}.{barcode}.r2.fastq{".gz" if args.gz else ""}', 'w'))
                   for barcode in barcodes}

    max_barcode_length = max([len(v) for v in barcodes.values()])
    pattern = None if args.allow_mismatch else pattern

    # f1, f2 = FastqGeneralIterator(opener(fastq1, 'r')), FastqGeneralIterator(opener(fastq2, 'r'))
    # with ProcessPoolExecutor(max_workers=args.cores) as executor:
    #     reads = executor.map(partial(barcoding, max_barcode_length=max_barcode_length,
    #                          randomer_length=args.randomer_length, allow_mismatch=args.allow_mismatch),
    #                          zip(f1, f2), chunksize=400000)
    #     for barcode, read1, read2 in reads:
    #         if barcode in writers:
    #             writer1, writer2 = writers[barcode]
    #             writer1.write(read1)
    #             writer2.write(read2)
    for r1, r2 in zip(FastqGeneralIterator(opener(fastq1, 'r')), FastqGeneralIterator(opener(fastq2, 'r'))):
        barcode, read1, read2 = barcoding((r1, r2),
                                          max_barcode_length=max_barcode_length,
                                          randomer_length=args.randomer_length,
                                          allow_mismatch=args.allow_mismatch)
        if barcode in writers:
            writer1, writer2 = writers[barcode]
            writer1.write(read1)
            writer2.write(read2)
    _ = [[v[0].close(), v[1].close()] for v in writers.values()]
    logger.debug(f'Demultiplexing {fastq1} and {fastq2} with barcodes {" and ".join(args.barcodes)} complete.')


if __name__ == '__main__':
    fastd()
