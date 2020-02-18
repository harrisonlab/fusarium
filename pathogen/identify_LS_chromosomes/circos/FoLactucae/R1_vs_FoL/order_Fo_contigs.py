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
        self.length_dict = defaultdict(int)
        self.locations_dict = defaultdict(list)
        self.best_match = ''
        self.sorted_contig_list = []
    def set_match(self):
        contig_list = self.length_dict.keys()
        self.sorted_contig_list = natsorted(contig_list, key=lambda contig: self.length_dict[contig])
        self.best_match = self.sorted_contig_list[-1]

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

# use link information to identify contigs with greatest alignment lengths
for line in links_lines:
    line = line.rstrip()
    # print(line)
    split_line = line.split("\t")
    contig = split_line[3]
    start = split_line[1]
    end = split_line[2]
    query_contig = split_line[0]
    ref_start = split_line[4]
    # print(ref_start)
    distance = (int(end) - int(start))
    if distance >= 15000:
        links_dict[contig].append(line)
    if not contigobj_dict[query_contig]:
        contigobj_dict[query_contig] = [ContigLinksObj()]
    contigobj_dict[query_contig][0].length_dict[contig] += distance
    contigobj_dict[query_contig][0].locations_dict[contig].append(int(ref_start))

fol_contig_dict = defaultdict(list)
# identify greatest alignment lengths in contigs
# input this into a dictionary of Fol contigs, allowing sorting by tupules of
# ref-contig and location of first hit
for contig in contigobj_dict.keys():
    contigobj_dict[contig][0].set_match()
    best_match = contigobj_dict[contig][0].best_match
    best_start = contigobj_dict[contig][0].locations_dict[best_match][0]
    # print("\t".join([contig, best_match, str(best_start)]))
    fol_contig_dict[best_match].append([contig, best_start])

# iterate through FoL contigs, sorting best matching query contigs by alignment
# start positions
contig_order = []
first = ''
last = ''
contig_breaks = []
for ref_contig in FoL_order:
    ref_contig = "FoL_" + ref_contig
    query_contigs_tup = fol_contig_dict[ref_contig]
    # print fol_contig_dict[ref_contig]
    # print(query_contigs_tup)
    sorted_query_contig_tup = natsorted(query_contigs_tup, key=lambda x: x[1])
    # print(sorted_query_contig_tup)
    sorted_query_contigs = [x[0] for x in sorted_query_contig_tup]
    # print(ref_contig)
    # print(", ".join(sorted_query_contigs))
    contig_order.extend(sorted_query_contigs)
    prev_last = last
    if sorted_query_contigs:
        first = sorted_query_contigs[0]
        contig_breaks.append(" ".join(["\t<pairwise", prev_last, first + ">"]))
        contig_breaks.append("\t\tspacing = 5000r")
        contig_breaks.append("\t</pairwise>")
        last = sorted_query_contigs[-1]

unseen_contigs = []
for contig in contig_list:
    if contig not in set(contig_order) and not contig.startswith('FoL_'):
        unseen_contigs.append(contig)

print("contig_order:")
rev_Fol = FoL_order[::-1]
rev_Fol = ["FoL_" + contig for contig in rev_Fol]
print(", ".join(contig_order + unseen_contigs + rev_Fol))



# print("monkeys")
# print links_dict
# print links_dict.keys()
seen_list = []
# contig_breaks = []
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
            # contig_breaks.append(" ".join(["\t<pairwise", prev_contig, query_contig + ">"]))
            # contig_breaks.append("\t\tspacing = 5000r")
            # contig_breaks.append("\t</pairwise>")
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

# print("contig order:")
# print(", ".join(seen_list + unseen_contigs + rev_Fol))

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
#
print("contigs in reverse orientation:")
print ", ".join(rev_list + rev_Fol)
