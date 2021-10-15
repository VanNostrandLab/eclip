#!/usr/bin/env perl

use warnings;
use strict;

unless ($ARGV[0] && $ARGV[1]) {
	print STDERR "Usage: perl compress_peak.pl normalized_full_bed compressed_bed\n\n";
	exit;
}

my $hashing_value = 100000;
my $bed = $ARGV[0];
my $compressed_bed = $ARGV[1];
my $compressed_full_bed = $compressed_bed;
$compressed_full_bed =~ s/\.bed/.tsv/;
open(O, ">$compressed_bed") || die "Cannot open $compressed_bed for writing!";
open(FULL, ">$compressed_full_bed") || die "Cannot open $compressed_full_bed for writing!";

my %peaks2size;
my %peaks2l2fc;
my %peaks2l10p;
my %peaks2start;
my %read_hash;
my %peak_hash;
my %long_lines;
&read($bed);

for my $chrom (keys %read_hash) {
    for my $strand ("+", "-") {
		my %discard_peaks;
		my @sorted_peaks = sort {$peaks2l10p{$chrom}{$strand}{$b} <=> $peaks2l10p{$chrom}{$strand}{$a} or
								$peaks2l2fc{$chrom}{$strand}{$b} <=> $peaks2l2fc{$chrom}{$strand}{$a} or
								$peaks2size{$chrom}{$strand}{$b} <=> $peaks2size{$chrom}{$strand}{$a} or
								$peaks2start{$chrom}{$strand}{$b} <=> $peaks2start{$chrom}{$strand}{$a}}
								keys %{$peaks2l10p{$chrom}{$strand}};

		for my $peak1 (@sorted_peaks) {
			my $verbose_flag = 1;
			next if (exists $discard_peaks{$peak1});
			my ($p1chrom, $p1position, $p1strand, $p1l10p, $p1l2fc) = split(/\:/, $peak1);
			my ($p1start, $p1stop) = split(/\-/, $p1position);
			my $p1x = int($p1start / $hashing_value);
			my $p1y = int( $p1stop / $hashing_value);
			for my $p1i ($p1x..$p1y) {
				for my $peak (@{$read_hash{$chrom}{$strand}{$p1i}}) {
					next if (exists $discard_peaks{$peak});
					next if ($peak eq $peak1);
					my ($p2chrom, $p2position, $p2strand, $p2l10p, $p2l2fc) = split(/\:/, $peak);
					my ($p2start, $p2stop) = split(/\-/, $p2position);
					next if ($p2stop <= $p1start);
					next if ($p1stop <= $p2start);

					if ($p1l10p >= $p2l10p) {
						# print STDERR "Comparing $peak1 and $peak, $peak was discarded\n" if ($verbose_flag == 1);
						$discard_peaks{$peak} = 1;
					} else {
						$discard_peaks{$peak1} = 1;
						# print STDERR "Comparing $peak1 and $peak, $peak1 was discarded\n" if ($verbose_flag == 1);
					}
				}
			}
		}

		for my $peak (@sorted_peaks) {
			next if (exists $discard_peaks{$peak});
			my ($p1chrom,$p1position,$p1strand,$p1l10p,$p1l2fc) = split(/\:/, $peak);
				my ($p1start, $p1stop) = split(/\-/, $p1position);
			# if ($p1l2fc > 0) {
			# 	print O "$p1chrom\t$p1start\t$p1stop\t$p1l10p\t$p1l2fc\t$p1strand\n";
			# }
			print O "$p1chrom\t$p1start\t$p1stop\t$p1l10p\t$p1l2fc\t$p1strand\n";
			print FULL "".$long_lines{$p1chrom.":".$p1position.":".$p1strand}."\n";
		}
    }
}
close(O);

sub read {
    my $file = shift;
    open(F, $file) || die "Cannot open $file\n";
    for my $line (<F>) {
		chomp($line);
		my @fields = split(/\t/, $line);
		my ($chrom, $position, $strand, $p) = split(/\:/, $fields[3]);
		my $start = $fields[1];
		my $stop = $fields[2];
		my $l10p = $fields[10];
		my $l2fc = $fields[11];

		my $peak_id = $chrom.":".$start."-".$stop.":".$strand.":".$l10p.":".$l2fc;
		push @{$peak_hash{$chrom}{$strand}},$peak_id;
		$peaks2start{$chrom}{$strand}{$peak_id} = $start;
		$peaks2l10p{$chrom}{$strand}{$peak_id} = $l10p;
		$peaks2l2fc{$chrom}{$strand}{$peak_id} = $l2fc;
		$peaks2size{$chrom}{$strand}{$peak_id} = $stop - $start;

		my $long_line = join("\t", @fields);
		$long_lines{$chrom.":".$start."-".($stop).":".$strand} = $long_line;

		my $x = int($start / $hashing_value);
		my $y = int( $stop / $hashing_value);
		for my $i ($x..$y) {
			push @{$read_hash{$chrom}{$strand}{$i}}, $peak_id
		}
    }
    close(F);
}
