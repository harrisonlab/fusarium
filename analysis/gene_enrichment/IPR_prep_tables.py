#!/usr/bin/python

'''
This program extracts Interproscan domains for genes and builds Fishcers
contingeny tables showing numbers of annotated genes on specified contigs
in comparison to other specified contigs
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
ap.add_argument('--interpro',required=True,type=str,help='The genome assembly')
ap.add_argument('--set1_name',required=True,type=str,help='Name for contigs in set 1')
ap.add_argument('--contig_set1',required=True,nargs='+',type=str,help='List of contigs in set 1')
ap.add_argument('--set2_name',required=True,type=str,help='Name for contigs in set 2')
ap.add_argument('--contig_set2',required=True,nargs='+',type=str,help='List of contigs in set 2')
ap.add_argument('--outdir',required=True,type=str,help='Output directory for results')

conf = ap.parse_args()

with open(conf.interpro) as f:
    interpro_lines = f.readlines()


set1_name = conf.set1_name
set1_contigs = set(conf.contig_set1)
set2_name = conf.set1_name
set2_contigs = set(conf.contig_set2)

print set1_name
print set1_contigs
print set2_name
print set2_contigs

# gff1_dict = defaultdict(list)
# specified_genes_dict = defaultdict(list)

#-----------------------------------------------------
# Step 2
# Store gene locations for organism 1 in a dictionary
#-----------------------------------------------------

for line in interpro_lines:
    line = line.rstrip()
    split_line = line.split("\t")
    gene_id = split_line[0]
    m = re.findall("IPR:......", line)
    IPR_set = set(m)
    if split_line[1] in set1_contigs:
        for IPR in IPR_set
            set1_dict[IPR] = +1
            set1_count = +1
    elif split_line[1] in set2_contigs:
        for IPR in IPR_set
            set2_dict[IPR] = +1
            set2_count = +1

print set1_count
print set2_count
