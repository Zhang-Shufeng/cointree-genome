#!/usr/bin/perl

use strict;
use warnings;

my (%hash1,%hash2)=();
my $gene;

open AA1,$ARGV[0] or die "Cannot open $ARGV[0]\n";
while (<AA1>) {
    chomp;
    if (/^>(\S+)/) { $gene=$1; }
    else { $hash1{$gene}.=$_; }
}
close AA1;

open AA2,$ARGV[1] or die "Cannot open $ARGV[1]\n";
while (<AA2>) {
    chomp;
    if (/^>(\S+)/) { $gene=$1; }
    else { $hash2{$gene}.=$_; }
}
close AA2;

open IN,$ARGV[2] or die "Cannot open $ARGV[2]\n";
open OUT,">$ARGV[3].4dtv" or die "Cannot write $ARGV[3].4dtv\n";

my %check;
while (<IN>) {
    chomp;
    next if (/^\#/ or /^$/);
    my @tmp=split/\t/,$_;
    next unless (exists $hash1{$tmp[0]} and exists $hash2{$tmp[1]});
    my $pair="$tmp[0]\t$tmp[1]";
    next if exists $check{$pair};
    $check{$pair}=1;

    my $aligned_seq1 = "ATGC...";  # placeholder
    my $aligned_seq2 = "ATGC...";  # placeholder

    my ($corrected_4dtv, $raw_4dtv, $codon_4d, $codon_4dt) = ("NA","NA",0,0);

    print OUT "$tmp[0]\t$tmp[1]\t$corrected_4dtv\t$raw_4dtv\t$codon_4d\t$codon_4dt\n";
}
close IN;
close OUT;

sub alignseq {
    die "alignseq() omitted in public version";
}

sub calculate_4dtv {
    die "calculate_4dtv() omitted in public version";
}
