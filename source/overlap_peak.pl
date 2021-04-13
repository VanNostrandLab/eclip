#!/usr/bin/env perl

# Adopted from https://github.com/YeoLab/eclip/blob/master/bin/overlap_peak_bed_with_bam_PE.pl

use warnings;
use strict;
use Statistics::Basic qw(:all);
use Statistics::Distributions;
use Statistics::R;
my $R = Statistics::R->new();

unless ($ARGV[0] && $ARGV[1] && $ARGV[2] && $ARGV[3] && $ARGV[4]) {
    print STDERR "usage: perl overlap_peak.pl IP_bam INPUT_bam peak_bed ip_mapped_read_number input_mapped_read_number output_file\n\n";
    exit;
}

my $ip_bam = $ARGV[0];
my $input_bam = $ARGV[1];
my $peak_bed = $ARGV[2];
my %mapped_read_count;
$mapped_read_count{"ip"} = $ARGV[3];
$mapped_read_count{"input"} = $ARGV[4];
my $output = $ARGV[5];

my $verbose = 0;
my %precalculated_fisher;

unless (exists $mapped_read_count{"ip"}) {
    print STDERR "Fatal error: missing experimental read count\n";
    exit;
}
unless (exists $mapped_read_count{"input"}) {
    print STDERR "Fatal error: missing input read count\n";
    exit;
}

my %peaks;
my $hashing_value = 100000;
my %peak_read_counts;
&read_peak_bed($peak_bed);
&read_bam($ip_bam,"ip");
&read_bam($input_bam,"input");

my $output_full = $output;
$output_full =~ s/\.bed/.full.bed/;
open(OUT,">$output");
open(OUTFULL,">$output_full");

for my $peak (keys %peak_read_counts) {
    unless (exists $peak_read_counts{$peak}{"ip"}) {
	    $peak_read_counts{$peak}{"ip"} = 1;
    }
    $peak_read_counts{$peak}{"input"}++;
    my $l2fc = log(($peak_read_counts{$peak}{"ip"}/$mapped_read_count{"ip"}) /
               ($peak_read_counts{$peak}{"input"}/$mapped_read_count{"input"})) / log(2);
    my ($p,$v,$method,$status) = &fisher_or_chi($peak_read_counts{$peak}{"ip"},
                                                $mapped_read_count{"ip"}-$peak_read_counts{$peak}{"ip"},
                                                $peak_read_counts{$peak}{"input"},
                                                $mapped_read_count{"input"}-$peak_read_counts{$peak}{"input"});
    my $l10p = $p > 0 ? -1 * log($p)/log(10) : 400 ;
    my ($chrom, $position, $strand, $p_value) = split(/\:/, $peak);
    my ($start, $stop) = split(/\-/, $position);
    print OUT "$chrom\t$start\t$stop\t$l10p\t$l2fc\t$strand\n";
    print OUTFULL "$chrom\t$start\t$stop\t$peak\t".$peak_read_counts{$peak}{"ip"}."\t"
                  .$peak_read_counts{$peak}{"input"}."\t$p\t$v\t$method\t$status\t$l10p\t$l2fc\n";
}
close(OUT);

sub max {
    my $x = shift;
    my $y = shift;
    if ($x > $y) {
        return($x);
    } else {
        return($y);
    }
}

sub min {
    my $x = shift;
    my $y = shift;
    if ($x < $y) {
        return($x);
    } else {
        return($y);
    }
}

sub fisher_or_chi {
    my ($a, $b, $c, $d) = @_;
    unless ($a && $b && $c && $d) {
    }
    my $tot = $a + $b + $c + $d;
    my $exp_a = ($a+$c)*($a+$b)/$tot;
    my $exp_b = ($b+$d)*($a+$b)/$tot;
    my $exp_c = ($a+$c)*($c+$d)/$tot;
    my $exp_d = ($b+$d)*($c+$d)/$tot;

    my $direction = "enriched";
    if ($a < $exp_a) {
        $direction = "depleted";
        return(1,"DEPL","N",$direction);
    }

    if ($exp_a < 5 || $exp_b < 5 || $exp_c < 5 || $exp_d < 5 || $a < 5 || $b < 5 || $c < 5 || $d < 5) {
        if (exists $precalculated_fisher{$a."|".$c}) {
            return($precalculated_fisher{$a."|".$c}{p}, $precalculated_fisher{$a."|".$c}{v}, "F", $direction);
        } else {
            my ($p,$v) = &fisher_exact($a,$b,$c,$d);
            $precalculated_fisher{$a."|".$c}{p} = $p;
            $precalculated_fisher{$a."|".$c}{v} = $v;
            return($p, $v, "F", $direction);
        }
    } else {
        my ($p, $v) = &chi_square($a, $b, $c, $d);
        return($p, $v,"C", $direction);
    }
}

sub chi_square {
    my ($a, $b, $c, $d) = @_;
    return(0) unless ($a && $b && $c && $d);
    my $tot = $a + $b + $c + $d;
    my $exp_a = ($a+$c)*($a+$b)/$tot;
    my $exp_b = ($b+$d)*($a+$b)/$tot;
    my $exp_c = ($a+$c)*($c+$d)/$tot;
    my $exp_d = ($b+$d)*($c+$d)/$tot;

    if ($exp_a >= 5 || $exp_b >= 5 || $exp_c >= 5 || $exp_d >= 5) {
        my $v = &square(&abs($a-$exp_a) - 0.5) / $exp_a +
                &square(&abs($b - $exp_b)-0.5) / $exp_b +
                &square(&abs($c - $exp_c)-0.5) / $exp_c +
                &square(&abs($d - $exp_d)-0.5) / $exp_d;
        my $p = Statistics::Distributions::chisqrprob(1,&abs($v));
        if ($a < $exp_a) {
            $v = $v * -1;
        }
        return ($p, $v);
    } else {
        print STDERR "ERROR: shouldn't get to this - should have been shunted into fisher exact test\n";
        return(1);
    }
}

sub fisher_exact {
    my ($x1, $x2, $y1, $y2) = @_;
    $R->run("rm(list = ls())");
    $R->run("blah <- matrix(c(".$x1.",".$x2.",".$y1.",".$y2."), nrow=2)");
    $R->run("foo <- fisher.test(blah)");
    my $p = $R->get('foo$p.value');
    my $v = "F";
    return($p, $v);
}

sub abs {
    my $x = shift;
    if ($x > 0) {
        return($x);
    } else {
        return(-1 * $x);
    }
}

sub square {
    my $x = shift;
    return($x * $x);
}

sub read_bam {
    my $bam = shift;
    my $label = shift;
    print STDERR "Reading $label $bam\n" if ($verbose == 1);
    if ($bam =~ /\.bam/) {
	    open(B,"samtools view -h $bam |") || die "no $bam\n";
    } elsif ($bam =~ /\.sam/) {
	    open(B,$bam) || die "no sam $bam\n";
    } else {
	    print STDERR "ERROR: file format error, file does not end with either .sam or .bam!\n";
	    exit;
    }
    while (<B>) {
	my $r1 = $_;
	if ($r1 =~ /^\@/) {
	    next;
	}
	my @fields = split(/\t/, $r1);
	my $flag = $fields[1];
	my $chrom = $fields[2];
	my $start = $fields[3];
    my $cigar = $fields[5];

	my $strand;
    if ($flag == 147 || $flag == 16) {
       $strand = "-";
    } elsif ($flag == 163 || $flag == 0) {
       $strand = "+";
    } else {
       print STDERR "R1 strand error, unexpected flag $flag\n";
    }

	my @read_regions = &parse_cigar_string($start, $cigar, $chrom, $strand);

	my %tmp_hash;
	for my $region (@read_regions) {
	    my ($r_chrom, $r_strand, $r_position) = split(/\:/, $region);
	    my ($r_start, $r_stop) = split(/\-/, $r_position);
	    my $rx = int($r_start / $hashing_value);
	    my $ry = int($r_stop  / $hashing_value);
	    for my $i ($rx..$ry) {
            for my $peak (@{$peaks{$r_chrom}{$r_strand}{$i}}) {
                my ($p_chr, $p_pos, $p_str, $p) = split(/\:/, $peak);
                my ($p_start, $p_stop) = split(/\-/, $p_pos);
                # my ($p_start, $p_stop) = split(/\-/, split(/\:/, $peak)[1]);
                next if ($p_start >= $r_stop || $p_stop <= $r_start);
                $tmp_hash{$peak} = 1;
            }
	    }
	}

	for my $peak (keys %tmp_hash) {
	    $peak_read_counts{$peak}{$label}++;
	}
    }
    close(B);
}

sub parse_cigar_string {
    my $region_start_pos = shift;
    my $flags = shift;
    my $chrom = shift;
    my $strand = shift;

    my $current_pos = $region_start_pos;
    my @regions;

    while ($flags =~ /(\d+)([A-Z])/g) {
        if ($2 eq "N") {
            push @regions,$chrom.":".$strand.":".($region_start_pos-1)."-".($current_pos-1);
            $current_pos += $1;
            $region_start_pos = $current_pos;
        } elsif ($2 eq "M") {
            #read and genome match
            $current_pos += $1;
        } elsif ($2 eq "S") {
            #beginning of read is soft-clipped; mapped pos is actually start pos of mapping not start of read
        } elsif ($2 eq "I") {
            #read has insertion relative to genome; does not change genome position
        } elsif ($2 eq "D") {
            $current_pos += $1;
            #read has deletion relative to genome; genome position has to increase
        } else {
            print STDERR "ERROR: unexpected operator and length: $1 and $2 in cigar string $flags.\n";

        }
    }
    push @regions,$chrom.":".$strand.":".($region_start_pos-1)."-".($current_pos-1);
    return(@regions);
}

sub read_peak_bed {
    my $bed = shift;
    print STDERR "Reading peak file $bed\n" if ($verbose == 1);
    open(F, $bed) || die "Cannot open peak file $bed\n";
    for my $line (<F>) {
        chomp($line);
        my @fields = split(/\t/,$line);
        my $chrom = shift(@fields);
        my $start = shift(@fields);
        my $stop = shift(@fields);
        my $gene = shift(@fields);
        my $p = shift(@fields);
        my $strand = shift(@fields);

        print STDERR "ERROR: missing start stop in $line\n" unless ($start && $stop);
        next unless ($start && $stop);

        my $x = int($start / $hashing_value);
        my $y = int($stop  / $hashing_value);

        my $peak = $chrom.":".$start."-".$stop.":".$strand.":".$p;
        $peak_read_counts{$peak}{allpeaks} = 1;


        for my $i ($x..$y) {
            push @{$peaks{$chrom}{$strand}{$i}},$chrom.":".$start."-".$stop.":".$strand.":".$p;
        }
    }
    close(F);
}