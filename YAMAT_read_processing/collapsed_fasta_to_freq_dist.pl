#!/usr/bin/perl

use warnings;
use strict;

my $usage = "\nUSAGE: perl $0 fasta_file\n\n";

# $0 is the name of the script

# my $file_name = shift or die "\nUsage: $0 FILENAME\n"; When you expect the user to provide a single filename on the command line.

# myfile = shift : shift removes and returns the first (left) element of the array; no array is provided as its parameter here, so
	#shift default works on @ARGV so the code will move the first value of @ARGV to the $file_name variable. i.e. if a file name is provided
	#then the script will run, if the @ARGV was empty then the die will be executed. Simply checks if a value was provided on the command line
	# and is then copied to $file_name. 	


#read filename from command line argument and store data in an array
my $file_name = shift or die ($usage);
my @fasta_lines = file_to_array ($file_name);
 
print ($file_name);
 
# @ array variable. This makes the file into an array, see subroutine below.
# $ scalar
 
print "Read_Count\tNumber_of_Unique_Seqs\n";

# define hash (=dictionary) to store count data
my %count_hash;

# loop over each line in the fasta file (now stored in array)
# extract the read count for that line and add 1 to the count of all sequences that share that read count [stored in a hash (=dictionary)].
 for (my $i = 0; $i < scalar (@fasta_lines); $i += 2){
 	
# 	remove newline character
 	chomp $fasta_lines[$i];
 	
# 	split the header line on hyphens
	my @split_line = split (/\-/, $fasta_lines[$i]);
	
# 	\ character escape
 	
# 	store last value from split line as our readcount;
	my $readcount = $split_line[-1];
# 	
	++$count_hash{$readcount};  
 }
 
 
 my @keys = keys %count_hash;
 
# sort keys numerically before printing and print output
 foreach my $this_key (sort {$a <=> $b} @keys){
 	print "$this_key\t$count_hash{$this_key}\n";
 }
# 
# ##subroutines

sub file_to_array {
	use strict;
	use warnings;

    my($filename) = @_;

    #Initialize variables
    my @filedata = (  );

    unless( open(GET_FILE_DATA, $filename) ) {
        print STDERR "Cannot open file \"$filename\"\n\n";
        exit;
    }

    @filedata = <GET_FILE_DATA>;

    close GET_FILE_DATA;

    return @filedata;
}