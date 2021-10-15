#!/usr/bin/env python
import os
import shutil
import sys
import tempfile

import fire
import cmder
import iToolBox as it


def cut_adapt(fastq1, fastq2, output1, output2, metric1, metric2, adapters='',
              barcode1='', barcode2='', barcodes_fasta='',
              randomer_length=10, cpus=1, debug=False):
    """
    Trim off adapters (2 rounds) for eCLIP dataset.

    :param fastq1: str, path to (gzip compressed) FASTQ file for read 1 ends with .fastq.gz extension.
    :param fastq2: str, path to (gzip compressed) FASTQ file for read 2 ends with .fastq.gz extension.
        For single-end dataset, using a special string 'none' or 'None' for fastq2 to avoid required argument error.
    :param output1: str, path to output file for read 1.
    :param output2: str, path to output file for read 2. For single-end dataset, using a special string
        'none' or 'None' for output2 to avoid required argument error.
    :param metric1: str, path to the first round trimming metric file.
    :param metric2: str, path to the second round trimming metric file.
    :param adapters: str, path to the adapters file (required for single-end dataset).
    :param barcode1: str, short name or code for first barcode, e.g., NIL, A01, ... (required for paired-end dataset).
    :param barcode2: str, short name or code for second barcode, e.g., NIL, A01, ... (required for paired-end dataset).
    :param barcodes_fasta: str, path to the barcodes fasta file (required for paired-end dataset).
    :param randomer_length: int, length of randomer.
    :param cpus: int, number of CPUs can be used for trimming adapters.
    :param debug: bool, set to True for invoking debug mode.
    """

    def parse_adapters(flag, fasta):
        adapters = []
        with open(fasta) as f:
            for line in f:
                if line.strip() and not line.startswith('>'):
                    adapters.extend([flag, line.strip()])
        return adapters

    def parse_overlap(txt):
        with open(txt) as f:
            return f.read().strip()

    def get_adapters():
        if adapters:
            adapters1, adapters2 = parse_adapters('-a', adapters), parse_adapters('-a', adapters)
            overlap1, overlap2 = '1', '5'
        else:
            cmd = ['parse_barcodes.sh', randomer_length, barcodes_fasta, barcode1, barcode2]
            folder = tempfile.mkdtemp(dir=outdir)
            cmder.run(cmd, msg='Parsing barcodes and finding adapters ...', cwd=folder)
            adapters1 = parse_adapters('-g', os.path.join(folder, 'g_adapters.fasta'))
            adapters1 += parse_adapters('-A', os.path.join(folder, 'A_adapters.fasta'))
            adapters1 += parse_adapters('-a', os.path.join(folder, 'a_adapters.fasta'))
            adapters2 = parse_adapters('-A', os.path.join(folder, 'A_adapters.fasta'))
            overlap1 = parse_overlap(os.path.join(folder, 'trim_first_overlap_length.txt'))
            overlap2 = parse_overlap(os.path.join(folder, 'trim_again_overlap_length.txt'))
            shutil.rmtree(folder)
        return adapters1, adapters2, overlap1, overlap2

    def get_ios(input1, input2, out1, out2):
        if input2 and out2:
            tmp1, tmp2 = f'{out1}.tmp', f'{out2}.tmp'
            ios1 = ['-o', tmp1, '-p', tmp2, input1, input2, '>', metric1]
            msg1 = f"Cutting adapters for paired reads {input1} and\n{' ' * 45}{input2} (first round) ..."

            ios2 = ['-o', out1, '-p', out2, tmp1, tmp2, '>', metric2]
            msg2 = f"Cutting adapters for paired reads {tmp1} and\n{' ' * 45}{tmp2} (second round) ..."
            need_to_delete.extend([tmp1, tmp2, out1, out2])
        else:
            tmp = f'{out1}.tmp'
            ios1 = ['-o', tmp, input1, '>', metric1]
            msg1 = f"Cutting adapters for single read {input1} (first round) ..."

            ios2 = ['-o', out1, tmp, '>', metric2]
            msg2 = f"Cutting adapters for single read {tmp} (second round) ..."
            need_to_delete.extend([tmp, out1])
        return ios1, ios2, msg1, msg2

    def trim_adapters(adapters, overlap, ios, message):
        cmd = ['cutadapt', '-O', overlap, '-j', cpus, '--match-read-wildcards', '--times', 1,
               '-e', 0.1, '--quality-cutoff', 6, '-m', 18] + adapters + ios
        if debug:
            cmder.run(cmd, msg=message, pmt=True, stdout=sys.stdout, stderr=sys.stderr)
        else:
            cmder.run(cmd, msg=message, pmt=True)

    need_to_delete = []
    outdir = os.path.dirname(output1) or os.getcwd()
    adapters1, adapters2, overlap1, overlap2 = get_adapters()
    fastq2 = '' if fastq2 in ('none', 'None', None) else fastq2
    output2 = '' if output2 in ('none', 'None', None) else output2
    ios1, ios2, msg1, msg2 = get_ios(fastq1, fastq2, output1, output2)
    trim_adapters(adapters1, overlap1, ios1, msg1)
    trim_adapters(adapters2, overlap2, ios2, msg2)
    it.rm(need_to_delete)


if __name__ == '__main__':
    fire.Fire(cut_adapt)
