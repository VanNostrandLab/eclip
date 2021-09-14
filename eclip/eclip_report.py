#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
A module for summarize QC of eCLIP dataset.
"""

import itertools
import os
import argparse
import tempfile

import cmder
import pandas as pd
from seqflow import Flow, task, logger
from bokeh.models import ColumnDataSource, FactorRange, NumeralTickFormatter, Legend
from bokeh.plotting import figure
from bokeh.models.widgets import Tabs, Panel
import datapane as dp

parser = argparse.ArgumentParser(description=__doc__, prog='eclip_report')
parser.add_argument('--input_fastqs', nargs='+', required=True, help='Path(s) to INPUT FASTQ file(s).')
parser.add_argument('--ip_fastqs', nargs='+', required=True, help='Path(s) to IP FASTQ file(s).')
parser.add_argument('--names', nargs='+', required=True, help='Shortnames for each sample, e.g., rep1, rep2.')
parser.add_argument('--wd', required=True, help='Path to the work directory that contains eCLIP analysis results.')
parser.add_argument('--cpus', type=int, help='Maximum number of CPU cores can be used for your job.', default=16)
parser.add_argument('--dry_run', action='store_true',
                    help='Print out steps and files involved in each step without actually running the pipeline.')
options = parser.parse_args()
try:
    os.chdir(options.wd)
except OSError as e:
    logger.error(e)


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
    
    return tabs


def report():
    dataset = [(name, ip_fastq, input_fastq)
               for ip_fastq, input_fastq, name in zip(options.ip_fastqs, options.input_fastqs, options.names)]
    dataset = pd.DataFrame(dataset, columns=['Samples', 'IP FASTQ', 'INPUT FASTQ'])
    
    cutadapt, repeat_map, genome_map, usable_reads = [], [], [], []
    for name in options.names:
        for read in ('ip', 'input'):
            sample = f'{name}.{read}'
            n, x1 = parse_cutadapt_metrics(f'{sample}.trim.trim.trim.metrics')
            _, x2 = parse_cutadapt_metrics(f'{sample}.trim.trim.metrics')
            _, x3 = parse_cutadapt_metrics(f'{sample}.trim.metrics')
            cutadapt.append([sample, n, x3, x2, x1])
            repeat_map.append([sample] + parse_star_log(f'{sample}.repeat.map.log'))
            genome_map.append([sample] + parse_star_log(f'{sample}.genome.map.log'))
            usable_reads.append([sample] + get_usable_reads(f'{sample}.genome.map.bam', f'{sample}.bam'))
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
    cutadapt = count_percent_plot(cutadapt, tooltips=(count_tips, percent_tips))
    
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
    repeat_map = count_percent_plot(repeat_map, tooltips=(count_tips, percent_tips))
    
    genome_map = pd.DataFrame(genome_map, columns=columns)
    genome_map = count_percent_plot(genome_map, tooltips=(count_tips, percent_tips))
    
    usable_reads = pd.DataFrame(usable_reads, columns=['Samples', 'Usable Reads', 'Duplicated Reads'])
    count_tips = [("Sample", "@Samples"),
                  ("# usable reads", "@{Usable Reads}{0.00 a}"),
                  ("# duplicated reads", "@{Duplicated Reads}{0.00 a}")]
    percent_tips = [("Sample", "@Samples"),
                    ("% usable reads", "@{Usable Reads}{0.00}"),
                    ("% duplicated reads", "@{Duplicated Reads}{0.00}")]
    usable_reads = count_percent_plot(usable_reads, tooltips=(count_tips, percent_tips))

    template = '/storage/vannostrand/software/eclip/data/report.md'
    data = {'dataset': dp.Table(dataset.style.hide_index()),
            'cutadapt': dp.Plot(cutadapt, responsive=False),
            'repeat_map': dp.Plot(repeat_map, responsive=False),
            'genome_map': dp.Plot(genome_map, responsive=False),
            'usable_reads': dp.Plot(usable_reads, responsive=False),
            'peak_clusters': dp.Text('peak_clusters'),
            'peak_annotation': dp.Text('peak_annotation'),
            'reproducibility': dp.Text('reproducibility'),
            }
    dp.Report(dp.Text(file=template).format(**data)).save(path='report.html', open=False)


if __name__ == '__main__':
    # flow = Flow('eCLIP_Rescue', description=__doc__.strip())
    # flow.run(dry_run=options.dry_run, cpus=options.cpus)
    report()
