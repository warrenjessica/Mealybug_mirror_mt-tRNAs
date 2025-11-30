#!/usr/bin/env perl

use strict;
use warnings;
use sloan;
use Bio::SearchIO;

my $usage = "\nUSAGE: perl $0 blast_file\n\n";

my $blast_file = shift or die ($usage);


my $total_count;
my $total_reads;
my $overhang_count;
my $overhang_reads;


my $SearchIO_obj = new Bio::SearchIO(-format => 'blast', -file   => "$blast_file");
while( my $result_obj = $SearchIO_obj->next_result ) {
	my $query_name = $result_obj->query_name;
	my $query_length = $result_obj->query_length;
	if ( my $hit_obj = $result_obj->next_hit ) {
		++$total_count;
		my @sn = split (/\_/, $query_name);
		$total_reads += $sn[1];		
		my $hsp_obj = $hit_obj->next_hsp;
		my $query_end = $hsp_obj->end('query');
		my $hit_end = $hsp_obj->end('hit');
		my $hit_length = $hit_obj->length;
		if ($query_end < $query_length and $hit_end == $hit_length){
			++$overhang_count;
			$overhang_reads += $sn[1];			
		}
	}
}

print "\nTotal mapped unique seqs:\t$total_count\n";
print "Total mapped_reads:\t$total_reads\n";
print "Mapped seqs with 3' overhangs:\t$overhang_count\n";
print "Mapped reads with 3' overhangs:\t$overhang_reads\n\n";
