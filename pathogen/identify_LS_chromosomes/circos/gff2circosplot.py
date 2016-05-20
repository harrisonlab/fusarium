#!/usr/bin/python

'''
This program is used to convert gff files into scatterplot data for circos
'''


#-----------------------------------------------------
# Step 1
# Import variables & load input files
#-----------------------------------------------------

import sys
import argparse
import re
import numpy
from sets import Set
from Bio import SeqIO
from collections import defaultdict
from operator import itemgetter

ap = argparse.ArgumentParser()
ap.add_argument('--genome',required=True,type=str,help='Assembly file')
ap.add_argument('--gff',required=True,type=str,help='A gff file of annotations')
ap.add_argument('--feature',required=True,type=str,help='The feature you want to extract / plot e.g. gene')

conf = ap.parse_args()
feat = conf.feature

with open(conf.gff) as f:
    gff_lines = f.readlines()


#-----------------------------------------------------
# Step 2
# Identify the length and gene density of FoL chromosomes
#-----------------------------------------------------

seq_len_dict = defaultdict(list)
genome_file = open(conf.genome, 'r')
for cur_record in SeqIO.parse(genome_file,"fasta"):
    seq_id = cur_record.id
    seq_len = len(cur_record.seq)
    seq_len_dict[seq_id] = seq_len

#-----------------------------------------------------
# Step 2
# Identify the length and gene density of FoL chromosomes
#-----------------------------------------------------


# contig_lgth = "10000000"
last_contig = "first"
for line in gff_lines:
    line = line.rstrip()
    if line.startswith("#"):
        continue
    split_line = line.split("\t")
    if str(feat) == split_line[2]:
        contig = str(split_line[0])
        if last_contig == "first":
            last_contig = str(contig)
            last_stop = str(0)
        elif last_contig != contig:
            contig_lgth = str(seq_len_dict[seq_id])
            if int(last_stop) >= int(contig_lgth):
                continue
            intermediate_space = "\t".join([last_contig, last_stop, contig_lgth, "0"])
            # intermediate_space = "\t".join([last_contig, last_stop, contig_lgth, "0"])
            print(intermediate_space)
            last_stop = str(0)
        start = str(split_line[3])
        stop = str(split_line[4])
        value = str("1")

        intermediate_space = "\t".join([contig, last_stop, start, "0"])
        print(intermediate_space)

        feature = "\t".join([contig, start, stop, value])
        print(feature)

        last_contig = contig
        last_stop = stop

intermediate_space = "\t".join([last_contig, last_stop, contig_lgth, "0"])
# intermediate_space = "\t".join([last_contig, last_stop, contig_lgth, "0"])
if int(last_stop) >= int(contig_lgth):
    print(intermediate_space)
