#!/usr/bin/perl;

use strict;
use warnings;

my $usage = "\nUSAGE: perl $0 input_fastq end_seq beginning_trim_length > output_fasta\n\n";

my $file = shift or die($usage);
my $end_seq = shift or die($usage);
my $beginning_trim_length = shift or die($usage);

my $FH = open_file($file);

while (my $line1 = <$FH>){

	chomp $line1;
	my $line2 = <$FH>;
	chomp $line2;
	my $line3 = <$FH>;
	my $line4 = <$FH>;
	
	substr($line2, -1 * length($end_seq)) eq $end_seq or next;
	length($line2) > length($end_seq) + $beginning_trim_length or next;
	
	print ">", substr($line1, 1), "\n", substr($line2, $beginning_trim_length, length($line2) - length($end_seq) - $beginning_trim_length), "\n";

}


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