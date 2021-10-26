#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
A pipeline for processing eCLIP data to identify genomic locations of RNA binding proteins (RBPs).
"""

import argparse
import glob
import os
import gzip
import sys
import time
import subprocess
import datetime
import tempfile
import shutil
import itertools
import math

import pysam
import cmder
import pandas as pd
from seqflow import Flow, task, logger
from bokeh.models import ColumnDataSource, FactorRange, NumeralTickFormatter, Legend
from bokeh.plotting import figure
from bokeh.models.widgets import Tabs, Panel
from bokeh.embed import components


def file_path(p):
    if p and not os.path.isfile(p):
        logger.error(f'File "{p}" may not be a file or does not exist.')
        sys.exit(1)
    return p


def dir_path(p):
    if p and not os.path.isdir(p):
        logger.error(f'Path "{p}" may not be a directory or does not exist.')
        sys.exit(1)
    return p


parser = argparse.ArgumentParser(description=__doc__, prog='eclip')
parser.add_argument('--input_fastqs', nargs='+', required=True, type=file_path, help='Path(s) to INPUT FASTQ file(s).')
parser.add_argument('--ip_fastqs', nargs='+', required=True, type=file_path, help='Path(s) to IP FASTQ file(s).')
parser.add_argument('--names', nargs='+', required=True, help='Shortnames for each sample, e.g., rep1, rep2.')
parser.add_argument('--adapters_fasta', help="Path to the FASTA file contains adapters and their sequences.")
parser.add_argument('--barcodes_pattern', help="Pattern of barcodes for umi-tools using to extract UMIs (for "
                                               "single-end dataset only, default: NNNNNNNNNN", default='NNNNNNNNNN')
parser.add_argument('--dataset', help='Name of the eCLIP dataset, default: eCLIP.', default='eCLIP')
parser.add_argument('--outdir', help="Path to the output directory, an pre-existing directory or current directory,"
                                     "default: current work directory.")
parser.add_argument('--species', help="Species name (short name code) the dataset associated with, e.g., hg19, mm10.",
                    default='hg19')
parser.add_argument('--genome', help="Path to STAR reference genome index directory.", type=dir_path)
parser.add_argument('--repeat', help="Path to STAR repeat elements index directory.", type=dir_path)
parser.add_argument('--blacklist_bed', help="Path to the bed file contains blacklist.")
parser.add_argument('--track', help="Name for the UCSC Genome Browser track, default: same as dataset")
parser.add_argument('--track_genome', help="Genome name (a short name code) for the track, default: same as species.")
parser.add_argument('--l2fc', type=int, help="Only consider peaks at or above this log2 fold change cutoff.", default=3)
parser.add_argument('--l10p', type=int, help="Only consider peaks at or above this log10 p value cutoff.", default=3)

parser.add_argument('--job', help="Name of your job", default='eCLIP')
parser.add_argument('--email', help='Email address for notifying you the start, end, and abort of you job.')
parser.add_argument('--time', type=int, help='Time (in integer hours) for running your job.', default=48)
parser.add_argument('--memory', type=int, help='Amount of memory (in GB) for all CPU cores.', default=32)
parser.add_argument('--cpus', type=int, help='Maximum number of CPU cores can be used for your job.', default=32)
parser.add_argument('--scheduler', choices=('pbs', 'qsub', 'slurm', 'sbatch', 'local'), default='slurm',
                    help='Name (case insensitive) of the scheduler on your cluster.')
parser.add_argument('--hold_submit', action='store_true',
                    help='Generate the submit script but hold it without submitting to the job scheduler.')

parser.add_argument('--debug', action='store_true', help='Invoke debug mode (for development use only).')
parser.add_argument('--dry_run', action='store_true',
                    help='Print out steps and files involved in each step without actually running the pipeline.')

START_TIME = time.perf_counter()


options = parser.parse_args()
setattr(options, 'outdir', options.outdir or os.getcwd())
if not os.path.isdir(options.outdir):
    try:
        os.mkdir(options.outdir)
    except OSError as e:
        logger.error(f'Create outdir failed: {e}.')
        sys.exit(1)
os.chdir(options.outdir)

adapters = '/storage/vannostrand/software/eclip/data/se.adapters.fasta'
setattr(options, 'adapters_fasta', options.adapters_fasta or adapters)
file_path(options.adapters_fasta)

setattr(options, 'repeat', options.repeat or '/storage/vannostrand/reference_data/hg19/repbase_v2_star_index')
dir_path(options.repeat)

setattr(options, 'genome', options.genome or '/storage/vannostrand/reference_data/hg19/genome_star_index')
dir_path(options.genome)

setattr(options, 'track', options.track or options.dataset)
setattr(options, 'track_genome', options.track_genome or options.species)

ips, inputs, names = options.ip_fastqs, options.input_fastqs, options.names
if len(ips) == len(names):
    if len(ips) == len(inputs):
        pass
    else:
        if len(inputs) == 1:
            inputs = inputs * len(ips)
        else:
            logger.error('Wrong number of input_fastqs were provided.')
            sys.exit(1)
else:
    logger.error('Number of items in ip_fastqs and names are not equal.')
    sys.exit(1)


class Sample:
    def __init__(self, sample_name):
        self.name = sample_name
        self.ip_fastq = f'{sample_name}.ip.fastq.gz'
        self.input_fastq = f'{sample_name}.input.fastq.gz'
        self.ip_bam = f'{sample_name}.ip.bam'
        self.input_bam = f'{sample_name}.input.bam'
        self.peak_bed = f'{sample_name}.ip.peak.clusters.bed'
        self.cross_bed = f'{sample_name}.ip.crosslink.sites.bed'
        
        
FASTQS, SAMPLES = {}, []
for ip_fastq, input_fastq, name in zip(ips, inputs, names):
    FASTQS[ip_fastq] = f'{name}.ip.fastq.gz'
    FASTQS[input_fastq] = f'{name}.input.fastq.gz'
    sample = Sample(name)
    SAMPLES.append(sample)

HASH = 100000
VERSION = 1.0
HEADER = fr"""
                           ____   _       ___   ____
                    ___   / ___| | |     |_ _| |  _ \
                   / _ \ | |     | |      | |  | |_) |
                  |  __/ | |___  | |___   | |  |  __/
                   \___|  \____| |_____| |___| |_|

                               VERSION {VERSION}
"""


@task(inputs=[fastq for fastq in FASTQS], outputs=lambda i: FASTQS[i])
def soft_link(fastq, link):
    """ Create soft links for original fastq files. """
    path = os.path.abspath(fastq)
    if path == os.path.abspath(link):
        logger.warning(f"No symbolic link was made for {path}! You are directly working on the original file!")
    else:
        if not os.path.islink(link):
            os.symlink(path, link)
    return link


@task(inputs=soft_link, outputs=lambda i: i.replace('.fastq.gz', '.umi.fastq.gz'), cpus=options.cpus)
def extract_umi(fastq, umi):
    message = f'Extract UMIs for {fastq} ...'
    cmd = ['umi_tools', 'extract',
           '--random-seed', 1,
           '--stdin', fastq,
           '--bc-pattern', 'NNNNNNNNNN',
           '--stdout', umi, '>', '/dev/null']
    cmder.run(cmd, msg=message, pmt=True)


@task(inputs=extract_umi, outputs=lambda i: i.replace('.umi.fastq.gz', '.trim.fastq.gz'))
def cut_adapt(fastq, trimmed_fastq):
    """ Trimming adapters (2 rounds) and 10 NTs from 3' using cutadapt (reads < 24 NTs will be discarded!)"""
    
    trim_tmp = trimmed_fastq.replace('.trim.fastq.gz', '.trim.tmp.fastq.gz')
    trim_trim_tmp = trimmed_fastq.replace('.trim.fastq.gz', '.trim.trim.tmp.fastq.gz')
    try:
        trim_metrics = trimmed_fastq.replace('.trim.fastq.gz', '.trim.metrics')
        trim_trim_metrics = trimmed_fastq.replace('.trim.fastq.gz', '.trim.trim.metrics')
        trim_trim_trim_metrics = trimmed_fastq.replace('.trim.fastq.gz', '.trim.trim.trim.metrics')
        cmd = ['cutadapt',
               '-j', options.cpus,
               '--match-read-wildcards',
               '--times', 1,
               '-e', 0.1,
               '--quality-cutoff', 6,
               '-m', 24,
               '-a', f'file:{options.adapters_fasta}']
        
        cmder.run(cmd + ['-O', 1, '-o', trim_tmp, fastq, '>', trim_metrics],
                  msg=f"Trimming adapters (1st round) for {fastq} ...", pmt=False)
        cmder.run(cmd + ['-O', 5, '-o', trim_trim_tmp, trim_tmp, '>', trim_trim_metrics],
                  msg=f"Trimming adapters (2nd round) for {fastq} ...", pmt=False)
        
        cmd = ['cutadapt', '-j', options.cpus, '-u', -10, '-o', trimmed_fastq, trim_trim_tmp,
               '>', trim_trim_trim_metrics]
        cmder.run(cmd, msg=f"Trimming adapters (3d round, 3'-UMIs) for {fastq} ...", pmt=False)
    finally:
        if os.path.exists(trim_tmp):
            os.unlink(trim_tmp)
        if os.path.exists(trim_trim_tmp):
            os.unlink(trim_trim_tmp)


@task(inputs=cut_adapt, outputs=lambda i: i.replace('.trim.fastq.gz', '.repeat.unmap.fastq.gz'))
def map_to_repeat_elements(fastq, mate):
    prefix = fastq.replace('.trim.fastq.gz', '.repeat.map')
    try:
        if not os.path.isdir(prefix):
            os.mkdir(prefix)
        cmd = ['STAR',
               '--runMode', 'alignReads',
               '--runThreadN', options.cpus,
               '--alignEndsType', 'EndToEnd',
               '--genomeDir', options.repeat,
               '--genomeLoad', 'NoSharedMemory',
               '--outBAMcompression', 10,
               '--outFileNamePrefix', f"{prefix}/",
               '--outFilterMultimapNmax', 100,
               '--outFilterMultimapScoreRange', 1,
               '--outFilterScoreMin', 10,
               '--outFilterType', 'BySJout',
               '--outReadsUnmapped', 'Fastx',
               '--outSAMattrRGline', 'ID:foo',
               '--outSAMattributes', 'All',
               '--outSAMmode', 'Full',
               '--outSAMtype', 'BAM', 'Unsorted',
               '--outSAMunmapped', 'None',
               '--outStd', 'Log',
               '--readFilesCommand', 'zcat',
               '--readFilesIn', fastq]
        message = f'Map SE reads in {fastq} to repeat elements ...'
        cmder.run(cmd, msg=message)
        cmder.run(f'mv {prefix}/Log.final.out {mate.replace(".repeat.unmap.fastq.gz", ".repeat.map.log")}')
        # cmder.run(f'mv {prefix}/Aligned.out.bam {mate.replace(".repeat.unmap.fastq.gz", ".repeat.map.bam")}')
        cmder.run(f'pigz -c -p {options.cpus} {prefix}/Unmapped.out.mate1 > {mate}')
    finally:
        shutil.rmtree(prefix)


@task(inputs=map_to_repeat_elements, outputs=lambda i: i.replace('.repeat.unmap.fastq.gz', '.genome.map.bam'))
def map_to_reference_genome(mate, bam):
    prefix = mate.replace('.repeat.unmap.fastq.gz', '.genome.map')
    try:
        if not os.path.isdir(prefix):
            os.mkdir(prefix)
        cmd = ['STAR',
               '--runMode', 'alignReads',
               '--runThreadN', options.cpus,
               '--alignEndsType', 'EndToEnd',
               '--genomeDir', options.genome,
               '--genomeLoad', 'NoSharedMemory',
               '--outBAMcompression', 10,
               '--outFileNamePrefix', f"{prefix}/",
               '--outFilterMultimapNmax', 1,
               '--outFilterMultimapScoreRange', 1,
               '--outFilterScoreMin', 10,
               '--outFilterType', 'BySJout',
               '--outReadsUnmapped', 'Fastx',
               '--outSAMattrRGline', 'ID:foo',
               '--outSAMattributes', 'All',
               '--outSAMmode', 'Full',
               '--outSAMtype', 'BAM', 'Unsorted',
               '--outSAMunmapped', 'None',
               '--outStd', 'Log',
               '--readFilesCommand', 'zcat',
               '--readFilesIn', mate]
        message = f'Map SE repeat elements unmapped reads in {mate} to reference genome ...'
        cmder.run(cmd, msg=message)
        cmder.run(f'mv {prefix}/Log.final.out {bam.replace(".genome.map.bam", ".genome.map.log")}')
        unmap = bam.replace('.genome.map.bam', '.genome.unmap.fastq.gz')
        cmder.run(f'pigz -c -p {options.cpus} {prefix}/Unmapped.out.mate1 > {unmap}')
        cmder.run(f'mv {prefix}/Aligned.out.bam {bam}')
    finally:
        shutil.rmtree(prefix)
    return bam


@task(inputs=map_to_reference_genome, outputs=lambda i: i.replace('.genome.map.bam', '.genome.map.sort.bam'))
def sort_index_bam(bam, out):
    cmder.run(f'samtools sort -@ {options.cpus} -m 2G -o {out} {bam}',
              msg=f'Sorting {bam} to {out} ...')
    cmder.run(f'samtools index {out}', msg=f'Indexing {out} ...')


@task(inputs=sort_index_bam, outputs=lambda i: i.replace('.genome.map.sort.bam', '.bam'))
def dedup_bam(bam, out):
    """Deduplicate SE BAM using umi_tools dedup."""
    cmd = ['umi_tools', 'dedup', '--random-seed', 1, '--stdin', bam, '--method', 'unique', '--stdout', out]
    cmder.run(cmd, msg=f'Deduplicating {bam} by umi_tools dedup ...')
    cmder.run(f'samtools index {out}', msg=f'Indexing {bam} ...')
    
    
@task(inputs=sort_index_bam,
      outputs=lambda i: i.replace('.genome.map.sort.bam', 'repetitive.elements.combine.with.unique.map.tsv.gz'))
def repetitive_elements_map(bam, tsv):
    cmd = ['repeat-maps', '--fastq', bam.replace('.genome.map.sort.bam', '.trim.fastq.gz'),
           '--bam', bam, '--dataset', bam.replace('.genome.map.sort.bam', ''), '--scheduler', 'local',
           '--cpus', options.cpus, '--species', options.species, '--outdir', options.outdir]
    cmder.run(cmd, msg=f'Mapping {bam.replace(".genome.map.sort.bam", "")} repetitive elements ...')


@task(inputs=dedup_bam, cpus=options.cpus, outputs=lambda i: i.replace('.bam', '.plus.bw'))
def make_bigwig_files(bam, bigwig):
    def bam_to_bigwig(bam, scale, strand, bw):
        bg = bw.replace('.bw', '.bg')
        cmd = f'genomeCoverageBed -ibam {bam} -bg -scale {scale} -strand {strand} -du -split | sort -k1,1 -k2,2n > {bg}'
        cmder.run(cmd, msg=f'Calculating genome coverage for {bam} ({strand} strand) ...')
        cmd = f'bedGraphToBigWig {bg} {options.genome}/chrNameLength.txt {bw}'
        cmder.run(cmd, msg=f'Converting {bg} to {bw} ...')
        cmder.run(f'rm {bg}')
    
    logger.info(f'Make BigWig files for {bam} ...')
    pos_bw, neg_bw = bigwig, bigwig.replace('.plus.bw', '.minus.bw')
    with pysam.AlignmentFile(bam, 'rb') as sam:
        total_reads = sam.mapped
    try:
        scale = 1000000.0 / total_reads
    except ZeroDivisionError:
        logger.warning(f'No reads was found in BAM {bam}, empty BigWig file was created.')
        with open(bigwig, 'w') as o:
            o.write('')
        return bigwig
    bam_to_bigwig(bam, scale, '+', pos_bw)
    bam_to_bigwig(bam, -1 * scale, '-', neg_bw)
    logger.info(f'Make BigWig files for {bam} complete.')
    return bigwig


@task(inputs=[], outputs=[f'{options.dataset}.hub.txt'], parent=make_bigwig_files)
def make_hub_file(inputs, output):
    logger.info('Make hub track file ...')
    header = f"""hub {options.track.replace(' ', '_')}
shortLabel {options.track}
longLabel {options.track}
useOneFile on
email {options.email if options.email else 'fei.yuan@bcm.edu'}

genome {options.track_genome}

track {options.track.replace(' ', '_')}
shortLabel {options.track}
longLabel {options.track}
type bigWig
superTrack on
"""
    block = """
track {basename}
shortLabel {basename}
longLabel {basename}
type bigWig
visibility full
alwaysZero on
autoScale on
aggregate transparentOverlay
showSubtrackColorOnUi on
parent {track}
container multiWig

    track {name1}
    bigDataUrl {plus}
    shortLabel {basename} Plus strand
    longLabel {basename} Plus strand
    type bigWig
    color 0,100,0
    parent {basename}

    track {name2}
    bigDataUrl {minus}
    shortLabel {basename} Minus strand
    longLabel {basename} Minus strand
    type bigWig
    color 100,0,0
    parent {basename}
    """
    
    track = options.track.replace(' ', '_')
    with open(output, 'w') as o:
        o.write(header)
        for bw in glob.iglob('*.plus.bw'):
            key = bw.replace('.plus.bw', '')
            plus, minus = f'{key}.pos.bw', f'{key}.neg.bw'
            name1, name2 = f'{key} plus', f'{key} minus'
            o.write(block.format(track=track, name1=name1, name2=name2, basename=key, plus=plus, minus=minus))
    logger.info('Make hub track file complete.')


def clipper_peaks(bam, bed=''):
    bed = bed if bed else bam.replace('.bam', '.peak.clusters.bed')
    if os.path.isfile(bed):
        logger.info(f'Clipper bed {bed} already exists.')
    else:
        cmd = f'clipper --species {options.species} --processors {options.cpus} --bam {bam} --outfile {bed}'
        cmder.run(cmd, msg=f'Calling peaks from {bam} using clipper ...', pmt=True)
    return bed


@task(inputs=[], outputs=[sample.peak_bed for sample in SAMPLES], parent=dedup_bam)
def clipper(bam, bed):
    bam = bed.replace('.peak.clusters.bed', '.bam')
    return clipper_peaks(bam, bed)


# @task(inputs=[],
#       outputs=[sample.cross_bed for sample in SAMPLES], parent=dedup_bam)
def pureclip(bam, bed):
    ip_bam, input_bam = [[sample.ip_bam, sample.input_bam] for sample in SAMPLES
                         if sample.cross_bed == bed][0]
    header = cmder.run(f'samtools view -H {ip_bam}').stdout.read()
    refs = [line.split()[1].replace('SN:', '') for line in header.splitlines() if line.startswith('@SQ')][:3]
    refs = ';'.join(refs)
    cmd = ['pureclip', '-i', ip_bam, '-bai', f'{ip_bam}.bai', '-g', f'{options.genome}/genome.fa',
           '-nt', options.cpus, '-ibam', input_bam, '-ibai', f'{input_bam}.bai', '-iv', f'"{refs};"',
           '-o', bed, '-or', bed.replace('.crosslink.sites.bed', '.binding.regions.bed'),
           '>', bed.replace('.crosslink.sites.bed', '.pureclip.log')]
    try:
        cmder.run(cmd, msg=f'Calling peaks from {ip_bam} and {input_bam} using pureCLIP ...')
    except Exception as e:
        logger.error(f'Running pureclip failed: {e}.')


def peak(ip_bams, input_bams, peak_beds, reproducible_bed, outdir):
    cmd = ['peak', '--ip_bams', ' '.join(ip_bams),
           '--input_bam', ' '.join(input_bams),
           '--peak_beds', ' '.join(peak_beds),
           '--read_type', 'SE',
           '--species', 'hg19' if options.species in ('hg19', 'hg19chr19') else options.species,
           '--outdir', outdir, '--cores', options.cpus,
           '--l2fc', options.l2fc, '--l10p', options.l10p]
    cmder.run(cmd, cwd=options.outdir, stdout=sys.stdout, stderr=sys.stderr)
    return reproducible_bed


@task(inputs=[],
      outputs=[f'{".vs.".join([f"{s.name}.ip" for s in SAMPLES])}.annotated.reproducible.peaks.bed'
               if len(SAMPLES) >= 2 else f'{SAMPLES[0].name}.ip.peak.clusters.normalized.compressed.annotated.bed'],
      parent=clipper)
def reproducible_peaks(inputs, outputs):
    ip_bams, input_bams, peak_beds, ids = [], [], [], []
    for sample in SAMPLES:
        ip_bams.append(sample.ip_bam)
        input_bams.append(sample.input_bam)
        peak_beds.append(sample.peak_bed)
        ids.append(sample.name)
    peak(ip_bams, input_bams, peak_beds, outputs, options.outdir)


def compile_motif_html(name, output=''):
    def parse_motif_html(html):
        lines = []
        with open(html) as h:
            for line in h:
                if line.startswith('Total target sequences'):
                    lines.append(line)
                elif line.startswith('Total background sequences'):
                    lines.append(line)
                    break
            for line in h:
                if line.startswith('</TABLE>'):
                    lines.append(line)
                    break
                else:
                    if '<BR/><A target=' in line:
                        line = line.split('<BR/><A target=')[0] + '</TD></TR>\n'
                    lines.append(line)
        s = ''.join(lines).replace('<TD>Motif File</TD></TR>', '</TR>')
        s = s.replace('border="1" cellpading="0" cellspacing="0"', f'class="text-center motif"')
        s = s.replace('possible false positive</FONT><BR/>', 'possible false positive</FONT><BR/><BR/>')
        s = s.replace('<TR><TD>Rank</TD>', '<thead><TR><TD>Rank</TD>')
        s = s.replace('STD(Bg STD)', 'STD (Bg STD)')
        s = s.replace('<TD>Best Match/Details</TD></TR>', '<TD>Best Match/Details</TD></TR></thead>')
        return s
    
    ul = """<ul class="nav nav-pills" id="MotifTab" role="tablist">
    {lis}
    </ul>
    """
    li = """<li class="nav-item mx-1 my-1" role="presentation">
    <button class="nav-link border border-primary py-1{active}" id="{tid}-tab" data-bs-toggle="pill"
    data-bs-target="#{tid}"
    type="button" role="tab">{name}</button>
    </li>
    """
    div = """<div class="tab-content" id="MotifTabContent">
    {divs}
    </div>
    """
    lis, divs, table_ids = [], [], []
    folder = f'{name}.motifs.40.min.plus.5p.20.3p.5'
    for i, html in enumerate(glob.iglob(f'{folder}/*.html')):
        region = html.split('.')[-4]
        rid = 'region_' + region
        text = parse_motif_html(html)
        if region == 'all':
            lis.append(li.format(tid=rid, active=' active', name=region))
            divs.append(f'<div class="py-3 tab-pane fade show active" id="{rid}" role="tabpanel">{text}</div>')
        else:
            lis.append(li.format(tid=rid, active='', name=region))
            divs.append(f'<div class="py-5 tab-pane fade" id="{rid}" role="tabpanel">{text}</div>')
    text = ul.format(lis='\n'.join(lis)) + '\n' + div.format(divs='\n'.join(divs))
    template = '/storage/vannostrand/software/eclip/data/motif.template.html'
    if output:
        with open(template) as f, open(output, 'w') as o:
            o.write(f.read().format(title=f'Motifs found in {name}', content=text))
    return text
    # cmder.run(f'zip -r -q {folder}.zip {folder}/*')
    # cmder.run(f'rm -rf {folder}')
    
    
@task(inputs=reproducible_peaks,
      outputs=lambda i: i.split('.annotated.reproducible.')[0] + '.motifs.40.min.plus.5p.20.3p.5.html'
      if '.annotated.reproducible.' in i
      else i.split('.ip.peak.clusters.')[0] + '.motifs.40.min.plus.5p.20.3p.5.html')
def motif_analysis(bed, output):
    basename = output.split('.motifs.')[0]
    cmd = ['motif', bed, options.species, options.outdir, basename, options.l10p, options.l2fc, options.cpus]
    cmder.run(cmd, msg=f'Finding motifs in {bed} ...')
    logger.info(f'Parsing and compiling motifs for {basename} ...')
    compile_motif_html(basename, output)
    logger.info(f'Parsing and compiling motifs for {basename} complete.')
    

def merge_bam(bams, bam):
    if os.path.isfile(bam):
        logger.info(f'BAM file {bam} already exist.')
    else:
        cmd = f'samtools merge {bam} {" ".join(bams)}'
        cmder.run(cmd, msg=f'Merging {" ".join(bams)} to {bam} ...')
    return bam


def split_bam(bam, basename, n):
    def count_mapped_reads(bam):
        count = int(cmder.run(f'samtools view -c -F 0x4 {bam}', msg='').stdout.read())
        logger.info(f'Found {count:,} mapped reads in {bam}.')
        return count

    bams = [f'{basename}{i}.bam' for i in range(n)]
    if all([os.path.isfile(b) for b in bams]):
        logger.info(f'Split bams already exist.')
    else:
        lines = int(count_mapped_reads(bam) / n) + 1
        cmd = f'samtools view {bam} | shuf | split - -a 1 --additional-suffix=.bam -d -l {lines} {basename}'
        cmder.run(cmd, msg=f'Shuffling and splitting {bam} ...')
        for b in bams:
            tmp = b.replace(".bam", ".tmp.bam")
            cmder.run(f'samtools view -H {bam} | cat - {b} | samtools view -bS - > {tmp}')
            cmder.run(f'samtools sort -@ {options.cpus} -o {b} {tmp}')
            cmder.run(f'rm {tmp}')
    return bams


def count_lines(file):
    lines = int(cmder.run(f'wc -l {file}').stdout.read().split()[0])
    logger.info(f'Found {lines:,} lines in {file}.')
    return lines


# @task(inputs=[], outputs=[f'{options.dataset}.rescue.ratio.txt'], parent=reproducible_peaks, mkdir=['rescue'])
def rescue_ratio(inputs, txt):
    if len(SAMPLES) == 1:
        logger.warning('No enough samples (n = 1 < 2) to calculate rescue ratio!')
        shutil.rmtree('rescue')
        return ''
    ip_bams, input_bams = [s.ip_bam for s in SAMPLES], [s.input_bam for s in SAMPLES]
    ip_pseudo_bam = merge_bam(ip_bams, os.path.join('rescue', 'ip.pseudo.bam'))
    ip_pseudo_bams = split_bam(ip_pseudo_bam, os.path.join('rescue', 'ip.pseudo.'), len(ip_bams))
    os.unlink(ip_pseudo_bam)
    
    input_pseudo_bam = merge_bam(input_bams, os.path.join('rescue', 'input.pseudo.bam'))
    input_pseudo_bams = split_bam(input_pseudo_bam, os.path.join('rescue', 'input.pseudo.'), len(input_bams))
    os.unlink(input_pseudo_bam)
    
    pseudo_peak_beds = [clipper_peaks(bam, bam.replace('.bam', '.peak.clusters.bed')) for bam in ip_pseudo_bams]
    basename = ".vs.".join([os.path.basename(bam).replace('.bam', '') for bam in ip_pseudo_bams])
    pseudo_peak_bed = os.path.join('rescue', f'{basename}.reproducible.peaks.bed')
    peak(ip_pseudo_bams, input_pseudo_bams, pseudo_peak_beds, pseudo_peak_bed, 'rescue')

    pseudo_count = count_lines(pseudo_peak_bed)
    basename = ".vs.".join([f'{name}.ip' for name in options.names])
    actual_count = count_lines(f'{basename}.reproducible.peaks.bed')
    try:
        ratio = max(actual_count, pseudo_count) / min(actual_count, pseudo_count)
    except ZeroDivisionError:
        ratio = 0
        logger.error(f'No peaks found in reproducible peaks or pseudo reproducible peaks, return ratio 0.')
    with open(txt, 'w') as o:
        o.write(f'{ratio}\n')
    

# @task(inputs=[], outputs=[f'{options.dataset}.consistency.ratio.txt'], parent=reproducible_peaks,
# mkdir=['consistency'])
def consistency_ratio(inputs, txt):
    if len(SAMPLES) == 1:
        logger.warning('No enough samples (n = 1 < 2) to calculate self-consistency ratio!')
        shutil.rmtree('consistency')
        return ''
    ip_bam1, ip_bam2, input_bam1, input_bam2, peak_bed1, peak_bed2 = [], [], [], [], [], []
    for s in SAMPLES:
        ip_b1, ip_b2 = split_bam(s.ip_bam, os.path.join('consistency', f'{s.name}.ip.split.'), 2)
        ip_bam1.append(ip_b1), ip_bam2.append(ip_b2)
        input_b1, input_b2 = split_bam(s.input_bam, os.path.join('consistency', f'{s.name}.input.split.'), 2)
        input_bam1.append(input_b1), input_bam2.append(input_b2)
        bed1 = clipper_peaks(ip_b1, ip_b1.replace('.bam', '.peak.clusters.bed'))
        bed2 = clipper_peaks(ip_b2, ip_b2.replace('.bam', '.peak.clusters.bed'))
        peak_bed1.append(bed1), peak_bed2.append(bed2)

    basename = ".vs.".join([os.path.basename(bam).replace('.bam', '') for bam in ip_bam1])
    split_peak_bed1 = os.path.join('consistency', f'{basename}.reproducible.peaks.bed')
    peak(ip_bam1, input_bam1, peak_bed1, split_peak_bed1, 'consistency')

    basename = ".vs.".join([os.path.basename(bam).replace('.bam', '') for bam in ip_bam2])
    split_peak_bed2 = os.path.join('consistency', f'{basename}.reproducible.peaks.bed')
    peak(ip_bam2, input_bam2, peak_bed2, split_peak_bed2, 'consistency')

    count1, count2 = count_lines(split_peak_bed1), count_lines(split_peak_bed2)
    try:
        ratio = count1 / count2
    except ZeroDivisionError:
        ratio = 0
        logger.error(f'No peaks found in one of the split reproducible peaks, return ratio 0.')
    with open(txt, 'w') as o:
        o.write(f'{ratio}\n')


def parse_cutadapt_metrics(metrics):
    reads, reads_too_short = 0, 0
    with open(metrics) as f:
        for line in f:
            if ('Reads written (passing filters)' in line) or ('Pairs written (passing filters)' in line):
                reads = int(line.strip().split()[-2].replace(',', ''))
            elif ('Reads that were too short' in line) or ('Pairs that were too short' in line):
                reads_too_short = int(line.strip().split()[-2].replace(',', ''))
    return reads, reads_too_short


def parse_star_log(log):
    counts = []
    with open(log) as f:
        for line in f:
            if 'Number of input reads' in line:
                reads = int(line.strip().split('\t')[-1])
            elif 'Uniquely mapped reads number' in line:
                counts.append(int(line.strip().split('\t')[-1]))
            elif 'Number of reads mapped to multiple loci' in line:
                counts.append(int(line.strip().split('\t')[-1]))
            elif 'Number of reads mapped to too many loci' in line:
                counts.append(int(line.strip().split('\t')[-1]))
            elif '% of reads unmapped: too many mismatches' in line:
                counts.append(int(float(line.strip().split('\t')[-1].replace('%', '')) * reads / 100))
            elif '% of reads unmapped: too short' in line:
                counts.append(int(float(line.strip().split('\t')[-1].replace('%', '')) * reads / 100))
            elif '% of reads unmapped: other' in line:
                counts.append(int(float(line.strip().split('\t')[-1].replace('%', '')) * reads / 100))
    return counts


def get_usable_reads(bam1, bam2):
    x = int(cmder.run(f'samtools view -c -F 0x4 {bam1}', msg='').stdout.read())
    y = int(cmder.run(f'samtools view -c -F 0x4 {bam2}', msg='').stdout.read())
    return [y, x - y]


def count_bar_plot(data, data_type, tooltips, colors):
    samples, categories = data['Samples'], [c for c in data.columns if c != 'Samples']
    colors = colors[:len(categories)]
    fig = figure(y_range=FactorRange(factors=samples[::-1]), plot_height=250, plot_width=1250,
                 toolbar_location=None, tooltips=tooltips)
    fig.add_layout(Legend(), 'right')
    fig.hbar_stack(categories, y='Samples', height=0.8, color=colors, legend_label=categories,
                   source=ColumnDataSource(data))
    fig.x_range.start = 0
    fig.x_range.range_padding = 0.1
    formatter = NumeralTickFormatter(format="0 a") if data_type == 'counts' else NumeralTickFormatter(format="0%")
    fig.xaxis[0].formatter = formatter
    fig.xaxis.axis_label = 'Number of reads' if data_type == 'counts' else 'Percent of reads'
    fig.xaxis.axis_line_color = None
    fig.y_range.range_padding = 0.1
    fig.yaxis.axis_line_color = None
    fig.ygrid.grid_line_color = None
    fig.legend.border_line_color = None
    fig.axis.minor_tick_line_color = None
    fig.axis.major_tick_line_color = None
    fig.outline_line_color = None
    return fig


def count_percent_plot(counts, tooltips=None):
    if tooltips:
        count_tooltips, percent_tooltips = tooltips
    else:
        count_tooltips = [("Sample", "@Samples"), ("Unique Reads", "@{Unique Reads}{0.00 a}"),
                          ("Duplicate Reads", "@{Duplicate Reads}{0.00 a}")]
        percent_tooltips = [("Sample", "@Samples"), ("Unique Reads (%)", "@{Unique Reads}{0.00}"),
                            ("Duplicate Reads (%)", "@{Duplicate Reads}{0.00}")]
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b',
              '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
    count_figure = count_bar_plot(counts, 'counts', count_tooltips, colors)
    df = counts.copy()
    df = df.set_index('Samples')
    df = df.div(df.sum(axis=1), axis=0)
    df = df.reset_index()
    percent_figure = count_bar_plot(df, 'percent', percent_tooltips, colors)
    
    count_panel = Panel(child=count_figure, title='Count')
    percent_panel = Panel(child=percent_figure, title='Percentage')
    
    tabs = Tabs(tabs=[count_panel, percent_panel])
    script, div = components(tabs)
    return script, div


def peak_table_tabs():
    ul = """<ul class="nav nav-pills" id="PeakClusterTab" role="tablist">
    {lis}
    </ul>
    """
    li = """<li class="nav-item mx-1 my-1" role="presentation">
    <button class="nav-link border border-primary{active}" id="{tid}-tab" data-bs-toggle="pill" data-bs-target="#{tid}"
    type="button" role="tab">{name}</button>
    </li>
    """
    div = """<div class="tab-content" id="PeakClusterTabContent">
    {divs}
    </div>
    """
    lis, divs, table_ids = [], [], []
    for i, name in enumerate(options.names):
        tid = f'{name}.peak.clusters'.replace('.', '-')
        table_ids.append(tid)
        tsv = f'{name}.ip.peak.clusters.normalized.compressed.annotated.tsv'
        columns = ['peak', 'ip_reads', 'input_reads', 'status', 'l10p', 'l2fc', 'gene']
        df = pd.read_csv(tsv, sep='\t', usecols=columns)
        df['peak'] = df['peak'].str.rsplit(':', expand=True, n=1)[0]
        df = df[(df.l2fc >= 3) & (df.l10p >= 3)]
        s = df.to_html(index=False, escape=True, table_id=f'{tid}-table', justify='center',
                       classes='table table-bordered table-stripped table-responsive text-center')
        if i == 0:
            lis.append(li.format(tid=tid, active=' active', name=name))
            divs.append(f'<div class="py-5 tab-pane fade show active" id="{tid}" role="tabpanel">{s}</div>')
        else:
            lis.append(li.format(tid=tid, active='', name=name))
            divs.append(f'<div class="py-5 tab-pane fade" id="{tid}" role="tabpanel">{s}</div>')
    s = ul.format(lis='\n'.join(lis)) + '\n' + div.format(divs='\n'.join(divs))
    return s, table_ids


def dt_js(table_ids):
    s = """<script>
        $(document).ready( function () {{
            {tables}
        }} );
    </script>
    """.format(tables='\n'.join([f"$('#{i}-table').DataTable();" for i in table_ids]))
    return s


@task(inputs=[], outputs=[f'{options.dataset}.report.html'], parent=reproducible_peaks)
def report(txt, out):
    reads, cutadapt, repeat_map, genome_map, usable_reads = {}, [], [], [], []
    for name in options.names:
        for read in ('ip', 'input'):
            sample = f'{name}.{read}'
            n, x1 = parse_cutadapt_metrics(f'{sample}.trim.trim.trim.metrics')
            _, x2 = parse_cutadapt_metrics(f'{sample}.trim.trim.metrics')
            t, x3 = parse_cutadapt_metrics(f'{sample}.trim.metrics')
            reads[f'{name}.{read}'] = t + x3
            cutadapt.append([sample, n, x3, x2, x1])
            repeat_map.append([sample] + parse_star_log(f'{sample}.repeat.map.log'))
            genome_map.append([sample] + parse_star_log(f'{sample}.genome.map.log'))
            usable_reads.append([sample] + get_usable_reads(f'{sample}.genome.map.bam', f'{sample}.bam'))
    
    dataset = []
    for ip_fastq, input_fastq, name in zip(options.ip_fastqs, options.input_fastqs, options.names):
        dataset.append([name, 'IP', ip_fastq, f"{reads[f'{name}.ip']:,}"])
        dataset.append([name, 'INPUT', input_fastq, f"{reads[f'{name}.input']:,}"])
    dataset = pd.DataFrame(dataset, columns=['Sample', 'Read', 'FASTQ', 'Reads'])
    
    columns = ['Samples', 'Reads written (passing filters)', 'Reads that were too short (round 1)',
               'Reads that were too short (round 2)', 'Reads that were too short (round 3)']
    cutadapt = pd.DataFrame(cutadapt, columns=columns)
    count_tips = [("Sample", "@Samples"),
                  ("# reads written (passing filters)", "@{Reads written (passing filters)}{0.00 a}"),
                  ("# reads that were too short (round 1)", "@{Reads that were too short (round 1)}{0.00 a}"),
                  ("# reads that were too short (round 2)", "@{Reads that were too short (round 2)}{0.00 a}"),
                  ("# reads that were too short (round 3)", "@{Reads that were too short (round 3)}{0.00 a}")]
    percent_tips = [("Sample", "@Samples"),
                    ("% reads written (passing filters)", "@{Reads written (passing filters)}{0.00 a}"),
                    ("% reads that were too short (round 1)", "@{Reads that were too short (round 1)}{0.00}"),
                    ("% reads that were too short (round 2)", "@{Reads that were too short (round 2)}{0.00}"),
                    ("% reads that were too short (round 3)", "@{Reads that were too short (round 3)}{0.00}")]
    
    scripts = []
    script, cutadapt_div = count_percent_plot(cutadapt, tooltips=(count_tips, percent_tips))
    scripts.append(script)
    
    columns = ['Samples', 'Uniquely mapped reads', 'Reads mapped to multiple loci',
               'Reads mapped to too many loci', 'Reads unmapped: too many mismatches',
               'Reads unmapped: too short', 'Reads unmapped: other']
    count_tips = [("Sample", "@Samples"),
                  ("# reads mapped uniquely", "@{Uniquely mapped reads}{0.00 a}"),
                  ("# reads mapped to multiple loci (round 2)", "@{Reads mapped to multiple loci}{0.00 a}"),
                  ("# reads mapped to too many loci", "@{Reads mapped to too many loci}{0.00 a}"),
                  ("# reads unmapped: too many mismatches", "@{Reads unmapped: too many mismatches}{0.00 a}"),
                  ("# reads unmapped: too short", "@{Reads unmapped: too short}{0.00 a}"),
                  ("# reads unmapped: other", "@{Reads unmapped: other}{0.00 a}")]
    percent_tips = [("Sample", "@Samples"),
                    ("% reads mapped uniquely", "@{Uniquely mapped reads}{0.00}"),
                    ("% reads mapped to multiple loci", "@{Reads mapped to multiple loci}{0.00}"),
                    ("% reads mapped to too many loci", "@{Reads mapped to too many loci}{0.00}"),
                    ("% reads unmapped: too many mismatches", "@{Reads unmapped: too many mismatches}{0.00}"),
                    ("% reads unmapped: too short", "@{Reads unmapped: too short}{0.00}"),
                    ("% reads unmapped: other", "@{Reads unmapped: other}{0.00}")]
    repeat_map = pd.DataFrame(repeat_map, columns=columns)
    script, repeat_map_div = count_percent_plot(repeat_map, tooltips=(count_tips, percent_tips))
    scripts.append(script)
    
    genome_map = pd.DataFrame(genome_map, columns=columns)
    script, genome_map_div = count_percent_plot(genome_map, tooltips=(count_tips, percent_tips))
    scripts.append(script)
    
    usable_reads = pd.DataFrame(usable_reads, columns=['Samples', 'Usable Reads', 'Duplicated Reads'])
    count_tips = [("Sample", "@Samples"),
                  ("# usable reads", "@{Usable Reads}{0.00 a}"),
                  ("# duplicated reads", "@{Duplicated Reads}{0.00 a}")]
    percent_tips = [("Sample", "@Samples"),
                    ("% usable reads", "@{Usable Reads}{0.00}"),
                    ("% duplicated reads", "@{Duplicated Reads}{0.00}")]
    script, usable_reads_div = count_percent_plot(usable_reads, tooltips=(count_tips, percent_tips))
    scripts.append(script)
    
    peak_clusters, table_ids = peak_table_tabs()
    template = '/storage/vannostrand/software/eclip/data/report.html'
    data = {'dataset': dataset.to_html(index=False, escape=True, justify='center',
                                       classes='table table-bordered table-stripped text-center'),
            'cutadapt': cutadapt_div,
            'repeat_map': repeat_map_div,
            'genome_map': genome_map_div,
            'usable_reads': usable_reads_div,
            'peak_clusters': peak_clusters,
            'peak_annotation': 'peak_annotation',
            'reproducibility': 'reproducibility',
            'scripts': '\n    '.join(scripts),
            'data_table_js': dt_js(table_ids)
            }
    with open(template) as f, open(out, 'w') as o:
        o.write(f.read().format(**data))
        
      
# @task(inputs=[], outputs=['.log.metrics.zip'], parent=report)
def cleanup(inputs, output):
    cmder.run('rm *.umi.fastq.gz *.trim.fastq.gz *.unmap.fastq.gz || true')
    fastq_links = [p for p in glob.iglob("*.fastq.gz") if os.path.islink(p)]
    if fastq_links:
        cmder.run(f'rm {" ".join(fastq_links)}')
    cmder.run('rm *.repeat.map.bam *.genome.map.bam* *.genome.map.sort.bam* || true ')
    cmder.run('rm *.mapped.reads.count.txt || true')
    bam_links = [p for p in glob.iglob("*.bam") if os.path.islink(p)]
    if bam_links:
        cmder.run(f'rm {" ".join(bam_links)}')
    cmder.run('rm *.normalized.bed *.normalized.tsv')
    cmder.run('rm *.compressed.bed *.compressed.tsv')
    cmder.run(f'zip {output} *.log *.metrics')
    cmder.run('rm *.log *.metrics || true')


def schedule():
    sbatch = """#!/usr/bin/env bash

#SBATCH -n {cpus}                       # Number of cores (-n)
#SBATCH -N 1                        # Ensure that all cores are on one Node (-N)
#SBATCH -t {runtime}                  # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH --mem={memory}G                   # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --job-name={job}            # Short name for the job
#SBATCH --output=%j.{job}.log
"""
    sbatch_email = """
#SBATCH --mail-user={email}
#SBATCH --mail-type=ALL
"""
    pbs = """ #!/usr/bin/env bash

#PBS -l nodes=1:ppn={cores}
#PBS -l walltime={runtime}
#PBS -l vmem={memory}gb
#PBS -j oe
#PBS -N {jobname}
"""
    pbs_email = """
#PBS -M {email}
#PBS -m abe
"""

    if options.scheduler.upper() in ('PBS', 'QSUB'):
        runtime, directive, exe, mail = f'{options.time}:00:00', pbs, 'qsub', pbs_email
        project = '/project/vannostrand'
    elif options.scheduler.upper() in ('SLURM', 'SBATCH'):
        days, hours = divmod(options.time, 24)
        runtime, directive, exe, mail = f'{days}-{hours:02}:00', sbatch, 'sbatch', sbatch_email
        project = '/storage/vannostrand'
    else:
        raise ValueError(f'Unsupported scheduler: {options.scheduler}, see help for supported schedulers.')
    if options.debug:
        root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    else:
        root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    program = os.path.join(root, 'seclip')
    setattr(options, 'ip_fastqs', ' \\\n              '.join(ips))
    setattr(options, 'input_fastqs', ' \\\n                 '.join(inputs))
    setattr(options, 'names', ' '.join(names))
    setattr(options, 'runtime', runtime)
    setattr(options, 'project', project)
    setattr(options, 'scheduler', 'local')
    if options.debug:
        setattr(options, 'debug', '--debug ')
        setattr(options, 'program', os.path.join(root, 'rnaseq', parser.prog))
    else:
        setattr(options, 'debug', '')
        setattr(options, 'program', os.path.join(root, parser.prog))
    arguments = ('ip_fastqs', 'input_fastqs', 'names', 'adapters_fasta', 'barcodes_pattern', 'species',
                 'repeat', 'genome', 'outdir', 'dataset', 'blacklist_bed', 'track', 'track_genome',
                 'l2fc', 'l10p', 'enrichment_filter', 'cpus', 'email', 'scheduler')
    arguments = ' \\\n  '.join([f'--{k} {getattr(options, k)}' for k in arguments if getattr(options, k, '')])
    setattr(options, 'arguments', arguments)

    code = fr"""
export TMPDIR={project}/tmp
export TEMP={project}/tmp
export TMP={project}/tmp

{program} \
  {arguments}
"""
    
    text = [directive, mail, code] if options.email else [directive, code]
    text = ''.join(text).format(**vars(options))
    
    submitter = os.path.join(options.outdir, 'submit.sh')
    with open(submitter, 'w') as o:
        o.write(text)
    
    print(f'Job submit script was saved to:\n    {submitter}')
    if options.hold_submit:
        print(f'Job {options.job} was not submitted yet, submit it after carefully review the submit script using:')
        print(f'    {exe} {submitter}')
    else:
        subprocess.run([exe, submitter])
        print(f'Job {options.job} was successfully submitted with the following resources:')
        data = {'Job name:': options.job, 'Output directory:': options.outdir,
                'Number of cores:': options.cpus, 'Job memory:': options.memory,
                'Job runtime:': f'{runtime} (D-HH:MM)'}
        for k, v in data.items():
            print(f'{k:>20} {v}')


FOOTER = """
+======================================================================================================================+
|                                                                                                                      |
|                                                 MISSION ACCOMPLISHED                                                 |
|{hh_mm_ss}|
|                                                                                                                      |
+======================================================================================================================+
"""


@logger.catch()
def main():
    if options.scheduler and options.scheduler != 'local':
        schedule()
    else:
        keys = ('outdir', 'genome', 'repeat', 'processes')
        d = vars(options).copy()
        d['processes'] = options.cpus
        setting = '\n'.join([f'{k.title():>22}: {v}' for k, v in d.items() if k in keys])
        logger.debug(HEADER)
        logger.trace(f'\nRunning eCLIP using the following settings:\n\n{setting}\n')
        flow = Flow('eCLIP', description=__doc__.strip())
        flow.run(dry_run=options.dry_run, cpus=options.cpus)
        # find_motif()
        logger.debug('')
        run_time = str(datetime.timedelta(seconds=int(time.perf_counter() - START_TIME)))
        logger.trace(FOOTER.format(hh_mm_ss=f'time consumed: {run_time}'.upper().center(118)))


if __name__ == '__main__':
    main()
    # flow = Flow('eCLIP', description=__doc__.strip())
    # flow.flow_chart('seclip.jpg')
