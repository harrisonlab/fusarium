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

ap = argparse.ArgumentParser()
ap.add_argument('--gff',required=True,type=str,help='A gff file of annotations')
ap.add_argument('--feature',required=True,type=str,help='The feature you want to extract / plot e.g. gene')

conf = ap.parse_args()
feat = conf.feature

with open(conf.gff) as f:
    gff_lines = f.readlines()


#-----------------------------------------------------
# Step 2
# Extract feature start and stop positions and print
# them in circos format.
#-----------------------------------------------------

for line in gff_lines:
    if line.startswith("#"):
        continue
    line = line.rstrip()
    split_line = line.split("\t")
    if str(feat) == split_line[2]:
        contig = str(split_line[0])
        start = str(split_line[3])
        stop = str(split_line[4])
        value = str("1")
        feature = "\t".join([contig, start, stop, value])
        print(feature)
