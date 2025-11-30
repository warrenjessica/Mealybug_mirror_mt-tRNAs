#!/usr/bin/perl;

use strict;
use warnings;
use sloan;

my $usage = "\nUSAGE perl $0 fasta_file sam_file min_reads\n\n";

my $fasta_file = shift or die ($usage);
my $sam_file = shift or die ($usage);
my $min_reads = shift or die ($usage);



my %fasta = fasta2hash($fasta_file);
my @sam_lines = file_to_array($sam_file);

my %hits;

foreach (@sam_lines){
	substr($_,0,1) eq '@' and next;
	my @sl = split (/\t/, $_);
	$hits{$sl[0]} = 1;
}

foreach (sort keys %fasta){
	exists ($hits{$_}) and next;
	
	my @sn = split (/\_/, $_);
	
	$sn[1] >= $min_reads and print ">$_\n$fasta{$_}\n";
}
