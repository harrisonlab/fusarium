#!/usr/bin/perl
use strict;
use warnings;

# rename_fasta.pl
# A module that appends given strings to the beginning of fasta accessions.

my $usage = "rename_fasta.pl <infile.fa> <text_to_append> <more_text_to_append>";
my $infile = shift or die $usage;
my $addstring = join ('_', @ARGV) or die $usage;

open INFILE, $infile;

while (my $line = <INFILE>) {
	if ($line =~ s/>//) { print ">$addstring"."_$line";	}
	else { print "$line"; }
}

exit;





