#!/usr/bin/perl env

use strict;
use warnings;
use sloan;
use Bio::SearchIO;

my $usage = "\nUSAGE: perl $0 blast_file fasta_query_file hit_names_file\n\n";

my $blast_file = shift or die ($usage);
my $fasta_file = shift or die ($usage);
my $hit_names_file = shift or die ($usage);

my @sfn = split (/\./, $blast_file);
my $lib_name = $sfn[0];

my %fasta = fasta2hash($fasta_file);

my @hit_list = file_to_array($hit_names_file);
my %hit_hash;
foreach (@hit_list){
	chomp $_;
	$hit_hash{$_} = 1;
}

my $SearchIO_obj = new Bio::SearchIO(-format => 'blast', -file   => "$blast_file");
while( my $result_obj = $SearchIO_obj->next_result ) {
	my $query_name = $result_obj->query_name;
	if ( my $hit_obj = $result_obj->next_hit ) {
		my $hit_name = $hit_obj->name;
		if (exists $hit_hash{$hit_name}){
			print ">$hit_name\.$lib_name\.$query_name\n$fasta{$query_name}\n";
		}
	}
}
