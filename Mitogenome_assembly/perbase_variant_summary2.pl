#!/usr/bin/env perl

use strict;
use warnings;
use sloan;

my $usage = "\nUSAGE: perl $0 file_list_doc reference_fasta min_freq min_cov\n\n";

my $file_list_doc = shift or die ($usage);
my $fasta_file = shift or die ($usage);
my $min_freq = shift or die ($usage);
my $min_cov = shift or die ($usage);

my %fasta = fasta2hash($fasta_file);

my @file_list = file_to_array($file_list_doc);
my %freq_HoHoH; # key structure: file > reference_seq > nuc_position; values: variant freq 
my %vartype_HoHoH; # key structure: file > reference_seq > nuc_position; values: most commonn variant type
my %recorded_data_HoH; #key structure: reference_seq > nuc_position; values: 1 if any library has enough coverage and recorded a freq

foreach (@file_list){
	chomp $_;
	my @perbase_lines = file_to_array($_);
	shift @perbase_lines;
	
	foreach my $var_line (@perbase_lines){
		chomp $var_line;

		my @sl = split(/\t/, $var_line);
		my $ref_base = uc(substr($fasta{$sl[0]}, $sl[1] - 1, 1));
		
		my $coverage = $sl[2];
		my $variants = $sl[3] + $sl[4] + $sl[5] + $sl[6] + $sl[8] + $sl[9];
		
		if ($ref_base eq "A"){
			$variants -= $sl[3];
		}elsif ($ref_base eq "C"){
			$variants -= $sl[4];
		}elsif ($ref_base eq "G"){
			$variants -= $sl[5];
		}elsif ($ref_base eq "T"){
			$variants -= $sl[6];
		}else{
			print STDERR "WARNING: non-ACTG base at $sl[0] position $sl[1]\n";
		}


		my $max_var_count = 0;
		my $max_var = "";
		
		if ($ref_base ne "A" and $sl[3] > $max_var_count){
			$max_var_count = $sl[3];
			$max_var = "$ref_base\>A";
		}

		if ($ref_base ne "C" and $sl[4] > $max_var_count){
			$max_var_count = $sl[4];
			$max_var = "$ref_base\>C";
		}

		if ($ref_base ne "G" and $sl[5] > $max_var_count){
			$max_var_count = $sl[5];
			$max_var = "$ref_base\>G";
		}

		if ($ref_base ne "T" and $sl[6] > $max_var_count){
			$max_var_count = $sl[6];
			$max_var = "$ref_base\>T";
		}

		if ($sl[8] > $max_var_count){
			$max_var_count = $sl[8];
			$max_var = "INS";
		}

		if ($sl[9] > $max_var_count){
			$max_var_count = $sl[9];
			$max_var = "DEL";
		}


		if ($coverage >= $min_cov){
			$freq_HoHoH{$_}->{$sl[0]}->{$sl[1]} = $variants / $coverage;
			$recorded_data_HoH{$sl[0]}->{$sl[1]} = 1;
			$vartype_HoHoH{$_}->{$sl[0]}->{$sl[1]} = $max_var;
		}
	}
}

print "Reference_Seq\tPosition";
foreach (@file_list){
	chomp $_;
	print "\t$_";
}
print "\n";

foreach (sort keys %recorded_data_HoH){
	my %pos_hash = %{$recorded_data_HoH{$_}};
	foreach my $pos (sort {$a <=> $b} keys %pos_hash){
		
		my $variant_detected = 0;
		foreach my $file (@file_list){
			chomp $file;
			if (exists ($freq_HoHoH{$file}->{$_}->{$pos})){
				$freq_HoHoH{$file}->{$_}->{$pos} >= $min_freq and $variant_detected = 1 and last;
			}
		}
		
		if ($variant_detected){
			print "$_\t$pos";
			foreach my $file2 (@file_list){
				print "\t";
				if (exists ($freq_HoHoH{$file2}->{$_}->{$pos})){
					print sprintf("%.3f",$freq_HoHoH{$file2}->{$_}->{$pos});
				}else{
					print "N/A";
				}
			}

			foreach my $file3 (@file_list){
				print "\t";
				if (exists ($freq_HoHoH{$file3}->{$_}->{$pos})){
					if ($freq_HoHoH{$file3}->{$_}->{$pos} >= $min_freq){
						print $vartype_HoHoH{$file3}->{$_}->{$pos};
					}else{
						print "N/A";
					}
				}else{
					print "N/A";
				}
			}
			print "\n";
		}
	}	
}
