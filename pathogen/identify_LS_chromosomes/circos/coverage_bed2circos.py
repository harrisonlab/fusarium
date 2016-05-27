#!/usr/bin/python

'''
This program is used to convert coverage .bed files into scatterplot data for circos
'''

#-----------------------------------------------------
# Step 1
# Import variables & load input files
#-----------------------------------------------------

import sys
import argparse
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument('--bed',required=True,type=str,help='A bed file output from bedtools coverage')

conf = ap.parse_args()

with open(conf.bed) as f:
    bed_lines = f.readlines()

#-----------------------------------------------------
# Step 2
# Extract feature start and stop positions, and coverage
# per Kb and print before parsing into circos format.
#-----------------------------------------------------

for line in bed_lines:
    line = line.rstrip()
    split_line = line.split("\t")
    contig = split_line[0]
    start = split_line[3]
    stop = split_line[4]
    reads = float(split_line[9])
    lgth = float(split_line[11])
    # print (str(reads) + "\t" + str(lgth))
    reads_per_bp = np.divide(reads, lgth)
    # print(reads_per_bp)
    reads_per_kb = np.multiply(reads_per_bp, 1000)
    reads_per_kb = int(np.round_(reads_per_kb, decimals=0,out=None))
    feature = "\t".join([str(contig), str(start), str(stop), str(reads_per_kb)])
    print(feature)
