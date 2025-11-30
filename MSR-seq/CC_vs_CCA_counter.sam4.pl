#!/usr/bin/perl;

use strict;
use warnings;

my $file = shift or die ("\nUSAGE: provide a sam filename on the command line. Optionally, add --unique flag after the filename to restrict analysis to reads that had a single top mapping hit (no ties).\n\n");
my $unique;
if ($unique = shift){
	$unique eq "--unique" or die ("\nERROR: the only option that can be provided after the same file name is --unique.\n\n")
}


my $FH = open_file ($file);


my %ref_CC_count;
my %ref_CCA_count;
my %ref_other_count;

my %ref_length_hash;

while (my $line = <$FH>){

	if ($line =~ /^\@SQ\tSN\:([\w\-]+)\tLN\:(\d+)$/){
	
		$ref_length_hash{$1} = $2;
		$ref_CC_count{$1} = 0;
		$ref_CCA_count{$1} = 0;
		$ref_other_count{$1} = 0;
	
	}

	$line =~ /^\@/ and next;
	
	my @sl = split (/\t/, $line);

	my $seq = $sl[9];
	my $name = $sl[0];
	my $ref = $sl[2];
	my $start = $sl[3];
	my $cigar = $sl[5];
	my $AS = $sl[11];
	my $XS = $sl[12];
	
	if ($unique){
		my $top_score;
		my $second_score;
		
		if ($AS =~ /AS\:i\:([\-\d]+)$/){
			$top_score = $1;
		}else{
			die ("\nERROR: could not parse AS score value: $AS\n\n");
		}

		if ($XS =~ /XS\:i\:([\-\d]+)$/){
			$second_score = $1;
			$top_score > $second_score or next;
		}
	}
	
	my $end_pos = $start - 1;
	
	while (length ($cigar) > 0){
		
		if ($cigar =~ /^(\d+[DIM])/){
			my $aln_piece = $1;
			my $aln_len = substr ($aln_piece, 0, -1);
			unless (substr ($aln_piece, -1) eq "I"){
				$end_pos += $aln_len;
			}
		
			if (length ($cigar) > length($aln_piece)){
				$cigar = substr ($cigar, length($aln_piece));
			}else{
				$cigar = "";
			}
		}else{
			die ("\nERROR: Could not parse CIGAR string $cigar\.\n\n")
		}		
	}
	
	
	if (substr($seq, -2) eq "CC" and $end_pos == $ref_length_hash{$ref} - 1){
		++$ref_CC_count{$ref};
	}elsif (substr($seq, -3) eq "CCA" and $end_pos == $ref_length_hash{$ref}){
		++$ref_CCA_count{$ref};
	}else{
		++$ref_other_count{$ref};
	}
}

print "Ref_tRNA\tCCA_count\tCC_count\tOther_count\n";

my $CCA = 0;
my $CC = 0;
my $other = 0;

foreach (sort keys %ref_length_hash){

	print "$_\t$ref_CCA_count{$_}\t$ref_CC_count{$_}\t$ref_other_count{$_}\n";
	$CCA += $ref_CCA_count{$_};
	$CC += $ref_CC_count{$_};
	$other += $ref_other_count{$_};

}

print "Total\t$CCA\t$CC\t$other\n";

sub open_file {
	use strict;
	use warnings;

    my($filename) = @_;
    my $fh;

    unless(open($fh, $filename)) {
        print "Cannot open file $filename\n";
        exit;
    }
    return $fh;
}
