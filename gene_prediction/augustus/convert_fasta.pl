#!/usr/bin/perl

use strict;
use Bio::SeqIO;

my $usage = "convert_fasta.pl <fasta_file.fa> <[dna/protein]>";
my $infile = shift or die $usage;
my $seq_type = shift or die $usage;
my $itteration = 0;
my $seqio_obj;
my $seq_obj;

$seqio_obj = Bio::SeqIO->new(-file => "$infile", -format => "fasta" -alphabet => $seq_type );

while ($seq_obj = $seqio_obj->next_seq){	#	reads a line in a loop
    print ">".$seq_obj->id,"\n";
    print $seq_obj->seq,"\n";
}

exit

 