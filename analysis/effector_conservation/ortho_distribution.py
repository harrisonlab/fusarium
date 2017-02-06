#!/usr/bin/python

'''

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
# ap.add_argument('--FoC_genes_gff',required=True,type=str,help='A gff file of the genes from the FoC')
ap.add_argument('--FoC_effectorP',required=True,type=str,help='A file containing results of effectorP')
ap.add_argument('--FoC_orthogroup',required=True,type=str,help='A file containing results of orthology analysis')
ap.add_argument('--OrthoMCL_id',required=True,type=str,help='The identifier of this strain as used in the orthology analysis')
ap.add_argument('--OrthoMCL_all',required=True,type=str,nargs='+',help='The identifiers of all strains used in the orthology analysis')
# ap.add_argument('--OrthoMCL_path',required=True,type=str,nargs='+',help='The identifiers of pathogenic strains used in the orthology analysis')
# ap.add_argument('--OrthoMCL_nonpath',required=True,type=str,nargs='+',help='The identifiers of non-pathogenic strains used in the orthology analysis')
ap.add_argument('--contig_list',required=True,type=str,nargs='+',help='Contigs to run this analysis on')


conf = ap.parse_args()

# with open(conf.FoC_genes_gff) as f:
#     FoC_genes_lines = f.readlines()

with open(conf.FoC_effectorP) as f:
    FoC_effectorP_lines = f.readlines()

with open(conf.FoC_orthogroup) as f:
    FoC_orthogroup_lines = f.readlines()

column_list=[]


#-----------------------------------------------------
# Step 9
# Build a dictionary of FoC genes that have been
# predicted as effectors by effectorP
#-----------------------------------------------------

# Core_contigs = set(["contig_1_pilon", "contig_2_pilon", "contig_3_pilon", "contig_4_pilon", "contig_5_pilon", "contig_6_pilon"])
# PS_contigs = set(["contig_14_pilon", "contig_20_pilon", "contig_22_pilon"])
# Var_contigs = set(["contig_14_pilon", "contig_20_pilon", "contig_22_pilon"])
contig_set = set(conf.contig_list)
FoC_effectorP_set = set()
First = True
for line in FoC_effectorP_lines:
    if First == True:
        First = False
        continue
    line = line.rstrip("\n")
    split_line = line.split("\t")
    if 'mRNA' in split_line[2]:
        contig = split_line[0]
        if contig in contig_set:
        # if contig in Core_contigs:
            features = split_line[8]
            split_features = features.split(";")
            gene_id = split_features[0]
            # print gene_id
            gene_id = gene_id.replace("ID=", "")
            FoC_effectorP_set.add("Fus2|" + gene_id)
            # print FoC_effectorP_set
            # print len(FoC_effectorP_set)
        # m = re.findall("ID=.*;", line)
        # print line
        # print m
        # if m:
        #     gene_id = pop(m)
        #     print "hit"
        #     FoC_effectorP_set.add(gene_id)
        # else:
        #     IPR_set = set(["no_annotation"])

    # if 'Non-effector' in line:
    #     continue
    # elif 'Effector' in line:
    #     FoC_effectorP_set.add(gene_id)
    # gene_id = line


#-----------------------------------------------------
# Step 11
# Build a dictionary of orthogroups
#-----------------------------------------------------

strain_id = conf.OrthoMCL_id
all_isolates = conf.OrthoMCL_all
# path_isolates = conf.OrthoMCL_path
# non_path_isolates = conf.OrthoMCL_nonpath

strain_id = strain_id + "|"
gene_count = 0
single_Fus2_gene_orthogroup = 0
single_gene_orthogroup = 0
equal_number_count = 0

FoC_orthogroup_dict = defaultdict(list)
orthogroup_content_dict = defaultdict(list)

for line in FoC_orthogroup_lines:
    line = line.rstrip("\n")
    split_line = line.split(" ")
    orthogroup_id = split_line[0].replace(":", "")
    orthogroup_contents = []
    orthogroup_content_dict.clear()
    if any(x in line for x in FoC_effectorP_set):
        for isolate in all_isolates:
            num_genes = line.count((isolate + "|"))
            orthogroup_contents.append(str(isolate) + "(" + str(num_genes) + ")")
            content_str = ":".join(orthogroup_contents)
            orthogroup_content_dict[isolate] = num_genes

        for gene_id in split_line[1:]:
            if strain_id in gene_id:
                # gene_id = gene_id.replace(strain_id, "")
                if gene_id in FoC_effectorP_set:
                    # print "\t".join([gene_id, orthogroup_id, content_str])
                    gene_count += 1

        if orthogroup_content_dict["Fus2"] == 1:
            single_Fus2_gene_orthogroup += 1
        if all(orthogroup_content_dict[x] == 1 for x in all_isolates):
            single_gene_orthogroup += 1
        Fus2_count = orthogroup_content_dict["Fus2"]
        if all(orthogroup_content_dict[x] == Fus2_count for x in all_isolates):
            equal_number_count += 1

print "\t".join(["Total number of genes:", str(gene_count)])
print "\t".join(["Number of genes in a single copy orthogroups in Fus2:", str(single_Fus2_gene_orthogroup)])
print "\t".join(["Number of genes in a single copy orthogroups in all isolates:", str(single_gene_orthogroup)])
print "\t".join(["Number of orthogroups with the same number of genes:", str(equal_number_count)])
            #
            # if all(x in line for x in all_isolates):
            #     FoC_orthogroup_dict[gene_id].extend([orthogroup_id, "all_isolates", content_str])
            # elif all(x not in line for x in non_path_isolates) and all(x in line for x in path_isolates):
            #     FoC_orthogroup_dict[gene_id].extend([orthogroup_id, "path_isolates_all", content_str])
            # elif all(x in line for x in non_path_isolates) and all(x not in line for x in path_isolates):
            #     FoC_orthogroup_dict[gene_id].extend([orthogroup_id, "non_path_isolates_all", content_str])
            # elif any(x in line for x in all_isolates):
            #     FoC_orthogroup_dict[gene_id].extend([orthogroup_id, "some_isolates", content_str])



    # print ("\t".join(useful_columns))
