#!/usr/bin/perl
use strict;

# This program will parse .tiling files outputted from MUMmer to .gff files.

my $usage = "tiling2gff.pl <MUMmer_outfile.tiling> > <outfile.gff>";

my $infile = shift;
my @ao_line;
my $col1 = "";
my $col2 = "MUMmer";
my $col3 = "ls_homolog";
my $col4 = "";
my $col5 = "";
my $col6 = "";
my $col7 = "";
my $col8 = ".";
my $col9 = "";
 
open (INFILE, "$infile") or die "\nERROR: $infile could not be opened\n"; 

while (my $line = <INFILE>) {
	@ao_line = split ('\s', $line);
	if ($ao_line[0] =~ m/^>.*/) { 
		$col1 = substr $ao_line[0], 1;
		print "\n";
	}
	else {
		$col4 = @ao_line[0];
		$col5 = @ao_line[1];
		$col6 = @ao_line[5];
		$col7 = @ao_line[6];
		$col9 = @ao_line[7];
	print "$col1\t$col2\t$col3\t$col4\t$col5\t$col6\t$col7\t$col8\t$col9\n";
	}
}

exit;

