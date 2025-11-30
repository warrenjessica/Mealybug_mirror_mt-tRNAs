#1/usr/bin/perl;

use strict;
use warnings;

my $usage = "\nUSAGE: perl $0 input_fasta > output_fasta\n\n";

my $file = shift or die($usage);
my $FH = open_file($file);

my %seq_hash;

while (my $line1 = <$FH>){

	chomp $line1;
	my $line2 = <$FH>;
	chomp $line2;
	
	++$seq_hash{$line2};
}

my $seq_num = 0;
foreach my $key (sort { $seq_hash{$b} <=> $seq_hash{$a} } keys(%seq_hash)){
	++$seq_num;
	print ">Seq$seq_num\_$seq_hash{$key}\n$key\n";
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
