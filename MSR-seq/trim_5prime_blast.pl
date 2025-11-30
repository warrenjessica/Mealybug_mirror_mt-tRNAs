#!/usr/bin/env perl

use strict;
use warnings;
use sloan;
use Bio::SearchIO;

my $usage = "\nUSAGE: perl $0 blast_file fasta_blast_file fasta_for_screening_file > trimmed_fasta_file\n\n";

my $blast_file = shift or die ($usage);
my $fasta1_file = shift or die ($usage); #file used to generate blast output
my $fasta2_file = shift or die ($usage); #file you want to trim (can be same or different as above)

my %fasta1 = fasta2hash($fasta1_file);

my %trim_hash;
my %exclude_hash;


my $SearchIO_obj = new Bio::SearchIO(-format => 'blast', -file   => "$blast_file");
while( my $result_obj = $SearchIO_obj->next_result ) {
	my $query_name = $result_obj->query_name;
	if ( my $hit_obj = $result_obj->next_hit ) {
		my $hsp_obj = $hit_obj->next_hsp;
		my $query_start = $hsp_obj->start('query');
		my $hit_start = $hsp_obj->start('hit');
		if ($query_start > 1 and $hit_start == 1){
			$trim_hash{$fasta1{$query_name}} = $query_start - 1;
		}
	}else{
		$exclude_hash{$fasta1{$query_name}}=1;
	}
}

my ($ref1, $ref2) = get_fasta_names_and_seqs($fasta2_file);

my @headers = @{$ref1};
my @seqs = @{$ref2};

for (my $i = 0; $i < scalar (@headers); ++$i){

	exists($exclude_hash{$seqs[$i]}) and next;

	if (exists ($trim_hash{$seqs[$i]})){
		print ">$headers[$i]\n", substr($seqs[$i], $trim_hash{$seqs[$i]}), "\n";
	}else{
		print ">$headers[$i]\n$seqs[$i]\n";
	}
}