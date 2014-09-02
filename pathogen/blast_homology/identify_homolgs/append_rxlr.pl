#!/usr/bin/perl -w
use strict;
use Cwd;

my $usage = "append_rxlr.pl <infile_sp_rxlr.fa> > appended_rxlr_sp.fa\n";
my $line;
my $organism = "F.oxysporum_fsp_cepae";

for (@ARGV) {
	my @nameparts = split ('/', $_);
	my $genome = $nameparts[-3];
	$genome =~ s/\./_/g; 
	my $isolate = $nameparts[-2];
	my $file = $_;
	open (INFILE, "$file") or die $usage;
	while ($line = <INFILE>) {
		$genome =~ s/\|/_/g;
		$line =~ s/>/>$genome--$isolate--/g;
		$line =~ s/\s//g;
		$line =~ s/gi\|(.*)\|gb\|(.*)\|$organism/$1--$2--/g; 
		print "$line\n";
	}
	close (INFILE);
}

exit
		