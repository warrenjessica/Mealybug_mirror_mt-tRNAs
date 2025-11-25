#!/usr/bin/perl
use strict;
use warnings;
use Bio::SearchIO; 



my $usage = "\nUSAGE: $0  INPUTFILESUMMARY   tRNA_DB   E-VALUE   OUTPUT_BASE_NAME HIT_COV_CUTOFF\n\n";

#Summary file needs to be a tab delimited list with samplename	blastfilehits	readsfile


my $input = shift or die ($usage); #first item in the usage variable
my $DB = shift or die ($usage); #fasta reference file from blast search 
my $evalue_threshold = shift or die ($usage);
my $output = shift or die ($usage);
my $percent_cov = shift or die ($usage);

my %matrixHoH;  #just creating a hash up here to fill with the tRNA names 
my %matrixHoH_ties;  #just creating a hash up here to fill with the tRNA names and store counts of ties


my %hashDB= fasta2hash($DB); # turns fasta reference file into a hash (dictionary)
my @DBnames= sort(keys(%hashDB)); #sorts the hash keys (name entries in the dictionary) 

my @libraries= file_to_array($input); #an array of the input file summary


my $FH_SUMM = open_output($output.".summary.txt"); #subdirectory, ask Dan what this does again...

print $FH_SUMM "Library\tQuery_Name\tRead_Count\tSequence\tHit_name\tFilter\tCov_hit_filter\tGenome\tAA\tAnticodon\ttRNA_ID\tEvalue\tBitscore\tPercentID\tSNPs\tGaps\tQuery_Len\tQuery_Cov\tHit_Len\tHit_Cov\tStrand\tTied\n";

my $FH_Ns = open_output ($output."excluded_Ns.txt");

print $FH_Ns "Library\tExcluded YAMAT Families with Ns\tExcluded YAMAT read count\n";

foreach my $line (@libraries){

	chomp $line; #remove the newline char from each line in the input file summary

	my @library_array=split(/\t/, $line); #makes an array of each of the lines in the input file summary file, splitting on tab. This will allow us to call each element in the array separately using [1]

	foreach my $trna (@DBnames){ #for each tRNA in the DBname array (all the arabidopsis tRNA names in the db blast file used, which was turned into an array, and the keys were sorted above)
		$matrixHoH{$library_array[0]}->{$trna} = 0; # populating matrix with 0s
		$matrixHoH_ties{$library_array[0]}->{$trna} = 0; # populating matrix with 0s
	}

	my %queryhash = fasta2hash($library_array[2]); #making a hash of the second element in the library array, i.e. all the collapsed reads file i.e. the query sequences. We did this to grab the sequence of each hit
	
	my $SearchIO_obj = new Bio::SearchIO(-format => 'blast',-file => $library_array[1]); # Bioperl object, feeds in the first element of the summary file (blast file)

	my $N_fams = 0;
	my $N_count = 0;

	while( my $result_obj = $SearchIO_obj->next_result ) { #this loops through all the result objects, assigns the result obj variable to the parcer object structure.
		my $query_name = $result_obj->query_name; #creating wanted variables from parser
		
		my @split_name = split (/\-/, $query_name); #the collapser program links the sequence name and the number of reads in a single line, just splitting that here so I can sort on count.
		
		my $query_length = $result_obj->query_length; #variable

		
		if ($queryhash{$query_name} =~ /N/){
			$N_count += $split_name[-1];
			++$N_fams;
			next;
		}

		my $hit_bitscore = 0; #starting a counter to compare for bitscore ties
		my $hit_evalue; # wanted variables from parser to print
		my $hit_length;
		my $hsp_percent;
		my $hsp_length_query;
		my $hsp_length_hit;
		my $SNPs;
		my $gaps;
		my $strand;
		my @hit_names;
		
		while(my $hit_obj = $result_obj->next_hit ) { #loops through the hits
			my $new_hit_bitscore = $hit_obj->raw_score;	#creates of another hit bitscore variable for comparison to the next to account for bit score blast ties
			
			if ($hit_bitscore){ # loop through bitscores
				$hit_bitscore == $new_hit_bitscore or last; #this is the comparison command if the hit bitscore numerically (==) equals the new bit score continue, if not, or last (break loop)
			}
			
			$hit_bitscore = $new_hit_bitscore; #resets the bitscore number to the new bitscore so it can be compared again in the loop
			push (@hit_names, $hit_obj->name); # push
					
			if (scalar (@hit_names) == 1){ # loop if hit_name = 1 (converted ?? into scalar) equals a numeric 1 do following. 
				$hit_evalue = $hit_obj->significance; 
				$hit_length = $hit_obj->length;
				my $hsp_obj = $hit_obj->next_hsp;
				$hsp_percent = $hsp_obj->percent_identity;
				$hsp_length_hit = $hsp_obj->length('hit');
				$hsp_length_query = $hsp_obj->length('query');
				$gaps = $hsp_obj->gaps;
				my $hsp_length_aln =  $hsp_obj->length('total');
				my @match_array = $hsp_obj->matches('query');
				$SNPs = $hsp_length_aln - $match_array[0] - $gaps;
				$strand = $hsp_obj->strand("hit");
				
			}
		}
				
		
		print $FH_SUMM "$library_array[0]\t$query_name\t$split_name[-1]\t$queryhash{$query_name}\t";

		if (@hit_names){
			
			my @split_hit = split (/\-/, $hit_names[0]);
			
			my $genome = $split_hit[0];
			my $trna_species = $split_hit[1];
			my $trna_ID = $split_hit[2];
			
			my $anticodon = substr($trna_species, -3);
			my $AA = substr($trna_species, 0, length($trna_species) - 3);
			
			print $FH_SUMM "$hit_names[0]\t";
						
			if ($hit_evalue <= $evalue_threshold){
				print $FH_SUMM ".\t";
				
				if (scalar($hsp_length_hit/$hit_length) >= $percent_cov){
					print $FH_SUMM ".\t";
					if (scalar(@hit_names) >= 2){ # if the number of hit names (meaning a tie) is 2 or greater, do the following foreach loop
					
						foreach my $tie_hit (@hit_names){ #
						
							$matrixHoH{$library_array[0]}->{$tie_hit} += $split_name[-1] / scalar(@hit_names);
							$matrixHoH_ties{$library_array[0]}->{$tie_hit} += $split_name[-1] / scalar(@hit_names); #counting ties that will be printed out in a separate matrix
						}
					}else{
						$matrixHoH{$library_array[0]}->{$hit_names[0]} += $split_name[-1];
						}
				}else{
					print $FH_SUMM "cov_fail\t";
					}
			}else{
				print $FH_SUMM "FILTER\t\t";
				
			}
			
			unless ($trna_ID){
				print STDERR "Problem strong: $hit_names[0]\n";
			}
			print $FH_SUMM "$genome\t$AA\t$anticodon\t$trna_ID\t$hit_evalue\t$hit_bitscore\t$hsp_percent\t$SNPs\t$gaps\t$hsp_length_query\t";
			
			print $FH_SUMM $hsp_length_query/$query_length, "\t$hit_length\t", $hsp_length_hit/$hit_length, "\t$strand\t";
			
			if (scalar(@hit_names) >= 2){
				print $FH_SUMM "TIE:$hit_names[0]";
				for (my $i=1; $i < scalar (@hit_names); ++$i){
					print $FH_SUMM ";$hit_names[$i]";
				}
				print $FH_SUMM "\n";
			}else{
				print $FH_SUMM ".\n";
			}
		
		}else{
			print $FH_SUMM "NO_HIT\tFILTER\n";
		}	
	}
	print "$library_array[0]\t$N_fams\t$N_count\n";
}
	

	

my $FH_MAT = open_output ($output."countmatrix.txt");
my $FH_MAT_TIES = open_output ($output."countmatrix_ties.txt");

print $FH_MAT "Library";
print $FH_MAT_TIES "Library";

my @libs = sort keys %matrixHoH;

my @trnas = sort keys %{$matrixHoH{$libs[0]}};
	
foreach my $column_head (@trnas){
	print $FH_MAT "\t$column_head";
	print $FH_MAT_TIES "\t$column_head";
}
print $FH_MAT "\n";	
print $FH_MAT_TIES "\n";	

foreach (@libs){
	print $FH_MAT $_;
	print $FH_MAT_TIES $_;
	
	foreach my $second_key (@trnas){
		print $FH_MAT "\t", $matrixHoH{$_}->{$second_key};
		print $FH_MAT_TIES "\t", $matrixHoH_ties{$_}->{$second_key};
	}
	
	print $FH_MAT "\n";
	print $FH_MAT_TIES "\n";
	
}	









###############################################################################
#get_fasta_names_and_seqs
#returns matching arrays of fasta heders and seqs given a fastname filename

sub get_fasta_names_and_seqs {
	use strict;
	use warnings;
	my ($inputfilename) = @_;
	my @fasta_names = ();
	my @fasta_seqs= ();

		   
	unless ( open(FILEDATA, $inputfilename) ) {
		print STDERR "Cannot open file \"$inputfilename\"\n\n"; #print error message
		exit; #exit the program
	}	

	my @filedata = <FILEDATA>; #Read the lines of the file into an array
	close FILEDATA;
	
	my $seq_count = 0; #this will be used to keep track of the number of sequences
	foreach my $line (@filedata){
		if ($line =~ /^>/) { #if the line is a header line (begins with ">")...
			if ($line =~ /^>.*[\w]+/){
				my $partialLine = substr ($&, 1);
				push (@fasta_names, $partialLine); #add that line to an array of fasta names
				push (@fasta_seqs, ""); #and add a new blank element to an array of sequences
				++$seq_count; #also increment our counter which keeps track of sequence number
			}
		}else { #if the line's not blank or a header, add it to the current sequence
			$fasta_seqs[$seq_count-1] .= $line;
		}
	}
	for (my $i = 0; $i < scalar (@fasta_seqs); ++$i){
		$fasta_seqs[$i] =~s/\s//g;
	}
	
	return (\@fasta_names, \@fasta_seqs);
}

#BEGIN fasta2hash
#A subroutine that takes a fasta file and returns a hash with headers as keys and sequences as values

sub fasta2hash {
	my $fasta = shift @_ or die ("\nERROR: No fasta file name provided to fasta2hash\n\n");
	my %fastaHash = arrays2hash (get_fasta_names_and_seqs($fasta));
	return %fastaHash;
}

#BEGIN arrays2hash
#a subroutine to take two arrays of equal size and convert them to a hash
#keys taken from first array, values from the second.
#Note that arrays must be passed by reference (\@array1, \@array2)

sub arrays2hash {
	use strict;
	use warnings;

	(my $keyarray, my $valuearray) = @_;
	if (scalar(@$keyarray) != scalar(@$valuearray)) {
		die "Arrays differ in size: Mismatched number of keys and values"; 
	}
	
	my %newhash = ( );
	
	@newhash{ @$keyarray } = @$valuearray;


	return (%newhash);

}

###############################################################
# BEGIN file_to_array
#
# A subroutine to get data from a file given its filename

sub file_to_array {
	use strict;
	use warnings;

    my($filename) = @_;

    # Initialize variables
    my @filedata = (  );

    unless( open(GET_FILE_DATA, $filename) ) {
        print STDERR "Cannot open file \"$filename\"\n\n";
        exit;
    }

    @filedata = <GET_FILE_DATA>;

    close GET_FILE_DATA;

    return @filedata;
}
#END file_to_array
###############################################################

sub open_output {
	use strict;
	use warnings;

    my($filename) = @_;
    my $fh_output;

    unless(open($fh_output, ">$filename")) {
        print "Cannot open file $filename\n";
        exit;
    }
    return $fh_output;
}
