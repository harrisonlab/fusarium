#!/usr/bin/python

'''
These commands were used to parse results from Blast searches of FoC pathogen
unique genes against the FoL genome as well as BLast searches of known SIX genes
against the FoL genome. These commands take information on location of blast
hits & suppliment this information with information from other files detailing
genes intersected in the target genome as well as information on the query gene
such as its location in the queries' genome.
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
ap.add_argument('--results_table',required=True,type=str,help='The information on BLAST hits of Fo path genes against the FoC genome')
ap.add_argument('--FoC_interescted_reblast',required=True,type=str,help='A bed file of FoC genes intersecting the location of reciprocal blast hits')

conf = ap.parse_args()

with open(conf.results_table) as f:
    results_table_lines = f.readlines()

with open(conf.FoC_interescted_reblast) as f:
    FoC_interescted_lines = f.readlines()


#-----------------------------------------------------
# Step 2
# Build a dictionary of FoC genes intersecting the hit
# location of reciprocal blast queries.
#-----------------------------------------------------

FoC_reblast_dict = defaultdict(list)
for line in FoC_interescted_lines:
    line = line.rstrip()
    split_line = line.split("\t")
    if "gene" in split_line[11]:
        query_id = split_line[8].strip('";')
        query_id = query_id.replace('ID=', '').replace('_BlastHit_1', '').replace('extracted_hit_of_', '')
        column_list = []
        # print split_line
        intersect_id = split_line[17]
        FoC_reblast_dict[query_id]=[query_id, intersect_id]


#-----------------------------------------------------
# Step 3
# Read Blast csv file and store in a dictionary of hits.
#-----------------------------------------------------
# The csv file may contain columns showing relationship
# between query sequences. The first line of the file
# contains column headers. The number of column headers
# containg 'Grp' can be counted and this number of
# columns skipped when reading lines.
#
# Not all columns will contain blast hits. The number
# of columns in a line must be counted to check if a
# blast hit is present.


first=True
for line in results_table_lines:
    line = line.rstrip("\n")
    split_line = line.split("\t")
    if first == True:
        split_line.extend(["FoC_reblast_hit", "match"])
        print "\t".join(split_line)
        first = False
        continue
    query_id = split_line[0]
    FoC_gene_ID = split_line[8]
    FoC_reblast_line = ["", ""]
    if FoC_reblast_dict[query_id]:
        FoC_reblast_line=FoC_reblast_dict[query_id]
    reblast_hit = FoC_reblast_line[1]
    # print reblast_hit
    # print FoC_gene_ID
    if FoC_gene_ID in reblast_hit and FoC_gene_ID != "":
        split_line.extend([reblast_hit, "match"])
    else:
        split_line.extend([reblast_hit, ""])
    print ("\t".join(split_line))
