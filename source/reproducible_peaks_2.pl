#!/usr/bin/env perl

use warnings;
use strict;

print $ARGV[4];
print $ARGV[5];

my $hashing_value = 10000;
my $l2fc_cutoff = 3;
my $l10p_cutoff = 3;

my $idr_cutoff = 0.01;
my %idr_cutoffs = ("0.001" => "1000", "0.005" => "955", "0.01" => "830", "0.02" => "705", "0.03" => "632",
    "0.04" => "580", "0.05" => "540", "0.06" => "507", "0.07" => "479", "0.08" => "455", "0.09" => "434",
    "0.1" => "415", "0.2" => "290", "0.3" => "217", "0.4" => "165", "0.5" => "125", "1" => "0");

my $rep1_idr_merged_full_bed = $ARGV[0];
my $rep2_idr_merged_full_bed = $ARGV[1];

my $rep1_full_out = $ARGV[2];
my $rep2_full_out = $ARGV[3];
open(REP1FULL,">$rep1_full_out") || die "Cannot open $rep1_full_out for writing!\n";
open(REP2FULL,">$rep2_full_out") || die "Cannot open $rep2_full_out for writing!\n";

my $bed_output = $ARGV[4];
my $custom_bed_output = $ARGV[5];
open(BEDOUT,">$bed_output") || die "Cannot open $bed_output for writing!\n";
open(CUSTOMOUT,">$custom_bed_output") || die "Cannot open $custom_bed_output for writing!\n";

my $idr_file = $ARGV[8];
my %idr_output;
&parse_idr_file($idr_file);

my $file1 = $ARGV[6];
my $file2 = $ARGV[7];

my %idr_region2peaks;
&parse_file($file1);
&parse_file($file2);
my %peak_info;

&parse_input_norm_full_file($rep1_idr_merged_full_bed);
&parse_input_norm_full_file($rep2_idr_merged_full_bed);


my $count_significant=0;
for my $idr_region (keys %idr_region2peaks) {
    my %peak_geo_mean;
    for my $peak (keys %{$idr_region2peaks{$idr_region}}) {
        my $geometric_mean = log(sqrt((2 ** $peak_info{$peak}{$rep1_idr_merged_full_bed}{l2fc}) *
                                 (2 ** $peak_info{$peak}{$rep2_idr_merged_full_bed}{l2fc}) ))/log(2);
        $peak_geo_mean{$peak} = $geometric_mean;
    }

    my @peaks_sorted = sort {$peak_geo_mean{$b} <=> $peak_geo_mean{$a}} keys %peak_geo_mean;
    my %already_used_peaks;
    for my $peak (@peaks_sorted) {
        my ($chr, $pos, $str) = split(/\:/, $peak);
        my ($start, $stop) = split(/\-/, $pos);

        my $flag = 0;
        for my $peak2 (keys %already_used_peaks) {
            my ($chr2, $pos2, $str2) = split(/\:/, $peak2);
            my ($start2,$stop2) = split(/\-/, $pos2);

            next if ($start2 >= $stop);
            next if ($start >= $stop2);
            $flag = 1;
        }
        next if ($flag == 1);

        next unless ($peak_info{$peak}{$rep1_idr_merged_full_bed}{l10p} >= $l10p_cutoff &&
                     $peak_info{$peak}{$rep2_idr_merged_full_bed}{l10p} >= $l10p_cutoff);

        $already_used_peaks{$peak} = 1;
        print CUSTOMOUT "".$idr_region."\t".$peak."\t".$peak_geo_mean{$peak}."\t".
            $peak_info{$peak}{$rep1_idr_merged_full_bed}{l2fc}."\t".
            $peak_info{$peak}{$rep2_idr_merged_full_bed}{l2fc}."\t".
            $peak_info{$peak}{$rep1_idr_merged_full_bed}{l10p}."\t".
            $peak_info{$peak}{$rep2_idr_merged_full_bed}{l10p}."\n";

	if ($peak_geo_mean{$peak} >= $l2fc_cutoff &&
        $peak_info{$peak}{$rep1_idr_merged_full_bed}{l10p} >= $l10p_cutoff &&
        $peak_info{$peak}{$rep2_idr_merged_full_bed}{l10p} >= $l10p_cutoff) {
	    print BEDOUT "".$chr."\t".$start."\t".$stop."\t".&min($peak_info{$peak}{$rep1_idr_merged_full_bed}{l10p},
            $peak_info{$peak}{$rep2_idr_merged_full_bed}{l10p})."\t".$peak_geo_mean{$peak}."\t".$str."\n";
	}

	my @rep1_full = split(/\t/,$peak_info{$peak}{$rep1_idr_merged_full_bed}{full});
	$rep1_full[3] .= ":".$peak_geo_mean{$peak};
	my $rep1_full_join = join("\t",@rep1_full);
        print REP1FULL "".$rep1_full_join."\n";

	my @rep2_full = split(/\t/,$peak_info{$peak}{$rep2_idr_merged_full_bed}{full});
	$rep2_full[3] .= ":".$peak_geo_mean{$peak};
	my $rep2_full_join = join("\t",@rep2_full);
        print REP2FULL "".$rep2_full_join."\n";

        $count_significant++ if ($peak_geo_mean{$peak} >= 3);
    }
    for my $peak (@peaks_sorted) {
        my ($chr,$pos,$str) = split(/\:/,$peak);
        my ($start,$stop) = split(/\-/,$pos);

        my $flag = 0;
        next if (exists $already_used_peaks{$peak});

        for my $peak2 (keys %already_used_peaks) {
            my ($chr2,$pos2,$str2) = split(/\:/,$peak2);
            my ($start2,$stop2) = split(/\-/,$pos2);

            next if ($start2 >=$stop);
            next if ($start >= $stop2);
            $flag = 1;
        }
        next if ($flag == 1);
        $already_used_peaks{$peak} = 1;
        print CUSTOMOUT "".$idr_region."\t".$peak."\t".$peak_geo_mean{$peak}."\t".
            $peak_info{$peak}{$rep1_idr_merged_full_bed}{l2fc}."\t".
            $peak_info{$peak}{$rep2_idr_merged_full_bed}{l2fc}."\t".
            $peak_info{$peak}{$rep1_idr_merged_full_bed}{l10p}."\t".
            $peak_info{$peak}{$rep2_idr_merged_full_bed}{l10p}."\n";

	if ($peak_geo_mean{$peak} >= $l2fc_cutoff &&
        $peak_info{$peak}{$rep1_idr_merged_full_bed}{l10p} >= $l10p_cutoff &&
        $peak_info{$peak}{$rep2_idr_merged_full_bed}{l10p} >=$l10p_cutoff) {
            print BEDOUT "".$chr."\t".$start."\t".$stop."\t".&min($peak_info{$peak}{$rep1_idr_merged_full_bed}{l10p},
                $peak_info{$peak}{$rep2_idr_merged_full_bed}{l10p})."\t".$peak_geo_mean{$peak}."\t".$str."\n";

        }

	my @rep1_full = split(/\t/,$peak_info{$peak}{$rep1_idr_merged_full_bed}{full});
	$rep1_full[3] .= ":".$peak_geo_mean{$peak};
	my $rep1_full_join = join("\t",@rep1_full);
	print REP1FULL "".$rep1_full_join."\n";

	my @rep2_full = split(/\t/,$peak_info{$peak}{$rep2_idr_merged_full_bed}{full});
	$rep2_full[3] .= ":".$peak_geo_mean{$peak};
	my $rep2_full_join = join("\t",@rep2_full);
        print REP2FULL "".$rep2_full_join."\n";
    }
}

# print STDERR "IDR and geometric mean(fc) >= 3 && p-value >= 3 in both reps: $count_significant\n";
close(REP1FULL);
close(REP2FULL);
close(CUSTOMOUT);
close(BEDOUT);


sub parse_input_norm_full_file {
    my $input_norm_file = shift;
    open(INF,$input_norm_file) || die "Cannot open $input_norm_file for reading!\n";
    for my $line (<INF>) {
	chomp($line);
	my @tmp = split(/\t/,$line);
	my $chr = $tmp[0];
	my $start = $tmp[1];
	my $stop = $tmp[2];
	my ($chr2,$pos2,$str,$del) = split(/\:/,$tmp[3]);
	$tmp[3] = $chr2.":".$pos2.":".$str;
	my $l2fc = $tmp[11];
	my $l10p = $tmp[10];

	$peak_info{$chr.":".$start."-".$stop.":".$str}{$input_norm_file}{l2fc} = $l2fc;
        $peak_info{$chr.":".$start."-".$stop.":".$str}{$input_norm_file}{l10p} = $l10p;
	$peak_info{$chr.":".$start."-".$stop.":".$str}{$input_norm_file}{full} = join("\t",@tmp);

    }
    close(INF);

}

sub parse_file {
    my $file = shift;
    open(F,$file) || die "Cannot open $file for reading!\n";
    for my $line (<F>) {
        chomp($line);

        my @tmp = split(/\t/,$line);
        my $chr = $tmp[0];
        my $start = $tmp[1];
        my $stop = $tmp[2];

        my ($chr2,$pos2,$str,$p) = split(/\:/,$tmp[3]);
        my $entropy = $tmp[12];
        my $l2fc = $tmp[11];
        my $l10p = $tmp[10];

        next unless ($l2fc >= 3 && $l10p >= 3);

        my $x = int($start / $hashing_value);
        my $y = int($stop / $hashing_value);

        my %overlapping_idrs;

        for my $i ($x..$y) {
            for my $idr_peak (@{$idr_output{$chr}{$str}{$i}}) {
                my ($ichr,$ipos,$istr,$iidr) = split(/\:/,$idr_peak);
                my ($istart,$istop) = split(/\-/,$ipos);
                next if ($istart >= $stop);
                next if ($istop <= $start);
                $overlapping_idrs{$idr_peak} = $iidr;
            }
        }
	
        if (scalar(keys %overlapping_idrs) > 0) {

            if (scalar(keys %overlapping_idrs) > 1) {
                print STDERR "This should NEVER be hit - peak overlaps with more than one IDR region $line\n".join("\t",keys %overlapping_idrs)."\n---\n";
            }

            my @sorted_idr = keys %overlapping_idrs;
            my $overlapping_idrpeak = $sorted_idr[0];
            my ($ichr,$ipos,$istr,$iidr) = split(/\:/,$overlapping_idrpeak);
            if ($iidr >= $idr_cutoffs{$idr_cutoff}) {
                $idr_region2peaks{$overlapping_idrpeak}{$chr.":".$start."-".$stop.":".$str} = 1;
            }
        } else {
            # peak not in IDR list
        }
    }
    close(F);
}


sub parse_idr_file {
    my $idr_file = shift;
    open(ID,$idr_file) || die "Cannot open $idr_file for reading!\n";
    for my $line (<ID>) {
	chomp($line);
	my @tmp = split(/\t/,$line);
	
	my $chr = $tmp[0];
	my $start = $tmp[1];
	my $stop = $tmp[2];
	my $str = $tmp[5];
	
	my $idr_score = $tmp[4];
	
	my $x = int($start / $hashing_value);
	my $y = int($stop / $hashing_value);
        
	for my $i ($x..$y) {
	    push @{$idr_output{$chr}{$str}{$i}},$chr.":".$start."-".$stop.":".$str.":".$idr_score;
	}
    }
    close(ID);
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
