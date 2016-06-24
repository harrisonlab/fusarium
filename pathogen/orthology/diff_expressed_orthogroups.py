#!/usr/bin/python

'''
This program accepts files detailing the results of differential expression analysis
and identifies orthogroups that these genes belongs to. It parses the output
into a summary file.
'''

#-----------------------------------------------------
# Step 1
# Import variables & load input files
#-----------------------------------------------------

import sys
import argparse
import re
import decimal
from sets import Set
from collections import defaultdict
from operator import itemgetter

ap = argparse.ArgumentParser()
ap.add_argument('--RNAseq_tab',required=True,type=str,help='A tsv file containing logfold change and p-values of DEGs')
ap.add_argument('--FoC_orthogroup',required=True,type=str,help='A file containing results of orthology analysis')
ap.add_argument('--OrthoMCL_id',required=True,type=str,help='The identifier of this strain as used in the orthology analysis')
ap.add_argument('--OrthoMCL_all',required=True,type=str,nargs='+',help='The identifiers of all strains used in the orthology analysis')
ap.add_argument('--OrthoMCL_path',required=True,type=str,nargs='+',help='The identifiers of pathogenic strains used in the orthology analysis')
ap.add_argument('--OrthoMCL_nonpath',required=True,type=str,nargs='+',help='The identifiers of non-pathogenic strains used in the orthology analysis')


conf = ap.parse_args()


with open(conf.RNAseq_tab) as f:
    RNAseq_lines = f.readlines()

with open(conf.FoC_orthogroup) as f:
    FoC_orthogroup_lines = f.readlines()


#-----------------------------------------------------
# Step 2
# Collect information on contig length from the genome
# assembly file. This can be used to determine if a
# gene has 2kb of sequence assembled up and downstream.
# This important for knock out design.
#-----------------------------------------------------

gene_id_set = Set([])
RNAseq_dict = defaultdict(list)

first = True
for line in RNAseq_lines:
    if first == True:
        first = False
    else:
        line = line.rstrip()
        split_line = line.split("\t")
        gene_id=split_line[0]
        gene_id=gene_id.replace('.', '_')
        gene_id_set.add(gene_id)
        column_list=[]
        column_list=itemgetter(0,1,2,6)(split_line)
        RNAseq_dict[gene_id].extend(column_list)
        if column_list[3] == 'NA':
            RNAseq_dict[gene_id].append("")
        elif decimal.Decimal(column_list[3]) <= 0.05:
            RNAseq_dict[gene_id].append("P<0.05")
        else:
            RNAseq_dict[gene_id].append("")

#-----------------------------------------------------
# Step 11
# Build a dictionary of orthogroups
#-----------------------------------------------------

strain_id = conf.OrthoMCL_id
all_isolates = conf.OrthoMCL_all
path_isolates = conf.OrthoMCL_path
non_path_isolates = conf.OrthoMCL_nonpath

strain_id = strain_id + "|"

FoC_orthogroup_dict = defaultdict(list)
orthogroup_content_dict = defaultdict(list)

for line in FoC_orthogroup_lines:
    line = line.rstrip("\n")
    split_line = line.split(" ")
    orthogroup_id = split_line[0].replace(":", "")
    orthogroup_contents = []
    orthogroup_content_dict.clear()
    for isolate in all_isolates:
        num_genes = line.count((isolate + "|"))
        orthogroup_contents.append(str(isolate) + "(" + str(num_genes) + ")")
        content_str = ":".join(orthogroup_contents)
        orthogroup_content_dict[isolate] = num_genes

    path_numbers = []
    for isolate in path_isolates:
        path_numbers.append(orthogroup_content_dict[isolate])
    max_path = max(path_numbers)
    min_path = min(path_numbers)
    non_path_numbers = []
    for isolate in non_path_isolates:
        non_path_numbers.append(orthogroup_content_dict[isolate])
    max_non_path = max(non_path_numbers)
    min_non_path = min(non_path_numbers)
    if min_path > max_non_path:
        expansion_status = "pathogen_expanded"
    elif min_non_path > max_path:
        expansion_status = "non-pathogen_expanded"
    else:
        expansion_status = ""

    for gene_id in split_line[1:]:
        if strain_id in gene_id:
            gene_id = gene_id.replace(strain_id, "")
            gene_id = re.sub(r'\.t\d+$','', gene_id)
            if FoC_orthogroup_dict[gene_id] and FoC_orthogroup_dict[gene_id][0] == orthogroup_id:
                continue
            elif all(x in line for x in all_isolates):
                FoC_orthogroup_dict[gene_id].extend([orthogroup_id, "all_isolates", content_str, expansion_status])
            elif all(x not in line for x in non_path_isolates) and all(x in line for x in path_isolates):
                FoC_orthogroup_dict[gene_id].extend([orthogroup_id, "path_isolates_all", content_str, expansion_status])
            elif all(x in line for x in non_path_isolates) and all(x not in line for x in path_isolates):
                FoC_orthogroup_dict[gene_id].extend([orthogroup_id, "non_path_isolates_all", content_str, expansion_status])
            elif any(x in line for x in all_isolates):
                FoC_orthogroup_dict[gene_id].extend([orthogroup_id, "some_isolates", content_str, expansion_status])


for gene_id in gene_id_set:
    outline = RNAseq_dict[gene_id]
    outline.extend(FoC_orthogroup_dict[gene_id])
    print "\t".join(outline)
