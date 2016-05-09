#!/usr/bin/python

'''
This program is used to parse the fo47 gtf annotation into a format that can be
read by the Fo_build-annot_table.py script. It is not intended as a full gtf to
gff converter.
'''


#-----------------------------------------------------
# Step 1
# Import variables & load input files
#-----------------------------------------------------

import sys
import argparse
import re
from sets import Set
from collections import defaultdict
from operator import itemgetter

ap = argparse.ArgumentParser()
ap.add_argument('--gtf',required=True,type=str,help='Input gtf file')

conf = ap.parse_args()


with open(conf.gtf) as f:
    gtf_lines = f.readlines()


#-----------------------------------------------------
# Step 2
#
#-----------------------------------------------------

store_line = ["#gff-version3"]
for line in gtf_lines:
    line = line.rstrip()
    split_line = line.split("\t")
    if "start_codon" in split_line[2]:
        print "\t".join(store_line)
        gene_start = split_line[3]
        store_line = split_line
        store_line[2] = "mRNA"
        feature_info = store_line[8].split('"')
        gene_id = feature_info[3]
        store_line[8] = "ID=" + gene_id + ";"
    elif "stop_codon" in split_line[2]:
        # print line
        gene_stop = split_line[4]
        store_line[4] = gene_stop

print "\t".join(store_line)
