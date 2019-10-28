#!/usr/bin/python

'''
This script reorders a set of contigs based on FoL contigs and outputs locations
for contig breaks to put into the circos .conf file.
'''

#-----------------------------------------------------
# Step 1
# Import variables & load input files
#-----------------------------------------------------

import sys
import argparse
import re
from collections import defaultdict
from natsort import natsorted, ns
# pip install https://github.com/SethMMorton/natsort/archive/6.0.0.tar.gz --user

ap = argparse.ArgumentParser()
ap.add_argument('--karotype',required=True,type=str,help='The kartype file provided to circos')
ap.add_argument('--links',required=True,type=str,help='A links file provided to circos')

conf = ap.parse_args()

with open(conf.karotype) as f:
    karotype_lines = f.readlines()
with open(conf.links) as f:
    links_lines = f.readlines()


#-----------------------------------------------------
# Step 2
# Create objects
#-----------------------------------------------------


class ContigLinksObj(object):
    """
    A dictionary of contigs linked to this contig by whole genome alignments,
    including information on total link length. This allows identification of
    the contig with the greatest homology
    """
    def __init__(self):
        """set the object, defining the dictionary."""
        dict = defaultdict(int)
        best_match = ''
    def set_match(self):
        contig_list = dict.keys()
        sorted_contig_list = natsorted(contig_list, key=lambda contig: self.dict[contig])
        best_match = sorted_contig_list[-1]

#-----------------------------------------------------
# Step 3
# Set FoL contig order
#-----------------------------------------------------

# FoL_order = ["CM000589.1", "CM000590.1", "CM000591.1", "CM000592.1", "CM000593.1", "CM000594.1", "CM000595.1", "CM000596.1", "CM000597.1", "CM000598.1", "CM000599.1", "CM000600.1", "CM000601.1", "CM000602.2", "CM000603.1", "DS231725.1", "DS231729.1", "DS231735.1", "DS231739.1", "DS231741.1", "DS231743.1", "DS231744.1", "DS231745.1", "DS231747.1", "DS231749.1", "DS231750.1", "DS231751.1", "DS231752.1", "DS231753.1", "DS231754.1", "DS231755.1", "DS231756.1", "DS231757.1", "DS231758.1", "DS231759.1", "DS231760.1", "DS231761.1", "DS231762.1", "DS231763.1", "DS231764.1", "DS231765.1", "DS231766.1", "DS231767.1", "DS231768.1", "DS231769.1", "DS231770.1", "DS231771.1", "DS231772.1", "DS231773.1", "DS231774.1", "DS231775.1", "DS231776.1", "DS231777.1", "DS231778.1", "DS231779.1", "DS231780.1", "DS231781.1", "DS231782.1", "DS231783.1", "DS231784.1", "DS231785.1", "DS231786.1", "DS231787.1", "DS231788.1", "DS231789.1", "DS231790.1", "DS231791.1", "DS231792.1", "DS231793.1", "DS231794.1", "DS231795.1", "DS231796.1", "DS231797.1", "DS231798.1", "DS231799.1", "DS231800.1", "DS231801.1", "DS231802.1", "DS231803.1", "DS231804.1", "DS231805.1", "DS231806.1", "DS231807.1", "DS231808.1", "DS231809.1", "DS231810.1", "DS231811.1", "DS231812.1"]
FoL_order = ["CM000589.1", "CM000590.1",  "CM000592.1", "CM000593.1",  "CM000595.1", "CM000596.1", "CM000597.1", "CM000598.1", "CM000599.1", "CM000600.1", "CM000601.1", "CM000591.1", "CM000594.1", "CM000602.2", "CM000603.1", "DS231725.1", "DS231729.1", "DS231735.1", "DS231739.1", "DS231741.1", "DS231743.1", "DS231744.1", "DS231745.1", "DS231747.1", "DS231749.1", "DS231750.1", "DS231751.1", "DS231752.1", "DS231753.1", "DS231754.1", "DS231755.1", "DS231756.1", "DS231757.1", "DS231758.1", "DS231759.1", "DS231760.1", "DS231761.1", "DS231762.1", "DS231763.1", "DS231764.1", "DS231765.1", "DS231766.1", "DS231767.1", "DS231768.1", "DS231769.1", "DS231770.1", "DS231771.1", "DS231772.1", "DS231773.1", "DS231774.1", "DS231775.1", "DS231776.1", "DS231777.1", "DS231778.1", "DS231779.1", "DS231780.1", "DS231781.1", "DS231782.1", "DS231783.1", "DS231784.1", "DS231785.1", "DS231786.1", "DS231787.1", "DS231788.1", "DS231789.1", "DS231790.1", "DS231791.1", "DS231792.1", "DS231793.1", "DS231794.1", "DS231795.1", "DS231796.1", "DS231797.1", "DS231798.1", "DS231799.1", "DS231800.1", "DS231801.1", "DS231802.1", "DS231803.1", "DS231804.1", "DS231805.1", "DS231806.1", "DS231807.1", "DS231808.1", "DS231809.1", "DS231810.1", "DS231811.1", "DS231812.1"]




#-----------------------------------------------------
# Step 4
# Order contigs
#-----------------------------------------------------


links_dict = defaultdict(list)
contigobj_dict = defaultdict(list)

contig_list = []
for line in karotype_lines:
    line = line.rstrip()
    split_line = line.split(" ")
    contig = split_line[2]
    contig_list.append(contig)

for line in links_lines:
    line = line.rstrip()
    # print(line)
    split_line = line.split("\t")
    contig = split_line[3]
    start = split_line[1]
    end = split_line[2]
    ref_contig = split_line[0]
    distance = (int(end) - int(start))
    if distance >= 15000:
        links_dict[contig].append(line)
    if contigobj_dict[contig]:
        contigobj_dict[contig][0].dict[ref_contig] += distance
    else:
        contigobj_dict[contig] = [ContigLinksObj()]
        contigobj_dict[contig][0].dict[ref_contig] += distance

# print("monkeys")
# print links_dict
# print links_dict.keys()
seen_list = []
contig_breaks = []
prev_contig = ''
orientation_dict = defaultdict(list)
for contig in FoL_order:
    # print(contig)
    contig_lines = links_dict["FoL_" + contig]
    # print contig_lines
    sorted_contig_lines = natsorted(contig_lines, key=lambda line: line.split()[4])
    first = True
    for line in sorted_contig_lines:
        # print(line)
        query_contig = line.split()[0]
        query_start = int(line.split()[1])
        if query_contig in seen_list:
            seen = True
        else:
            seen = False
            seen_list.append(query_contig)
        if first and not seen:
            contig_breaks.append(" ".join(["\t<pairwise", prev_contig, query_contig + ">"]))
            contig_breaks.append("\t\tspacing = 5000r")
            contig_breaks.append("\t</pairwise>")
            first = False
        if prev_contig == query_contig:
            if query_start >= prev_start:
                orientation = '+'
            else:
                orientation = '-'
            orientation_dict[query_contig].append(orientation)
        prev_contig = query_contig
        prev_start = query_start



unseen_contigs = []
for contig in contig_list:
    if contig in seen_list:
        continue
    elif contig.startswith("FoL"):
        continue
    else:
        unseen_contigs.append(contig)

rev_Fol = FoL_order[::-1]
rev_Fol = ["FoL_" + contig for contig in rev_Fol]

for contig in rev_Fol:
    if contig.startswith("FoL_CM"):
        contig_breaks.append(" ".join(["\t<pairwise", prev_contig, contig + ">"]))
        contig_breaks.append("\t\tspacing = 5000r")
        contig_breaks.append("\t</pairwise>")
    prev_contig = contig

# print(rev_Fol)
# print(FoL_order[::-1])
print("contig order:")
print(", ".join(seen_list + unseen_contigs + rev_Fol))

def most_common(lst):
    return max(set(lst), key=lst.count)

for_list = []
rev_list = []
for contig in (seen_list + unseen_contigs):
    if orientation_dict[contig]:
        orientation_list = orientation_dict[contig]
        orientation = most_common(orientation_list)
    else:
        orientation = '+'
    # print contig + "\t" + orientation
    if orientation == '+':
        for_list.append(contig)
    else:
        rev_list.append(contig)

print("contig breaks")
print("\n".join(contig_breaks))

print("contigs in reverse orientation:")
print ", ".join(rev_list + rev_Fol)
# sorted_links_lines = sorted(links_lines, key=lambda line: line.split()[0]):
#
# for line in
