#!/usr/bin/python

'''
This program is used to extract information on NLR genes containing STAND domains
from summary tables of FoL genes. It outputs an oragsnims genes containing
NLR genes (identified by N-terminal STAND domain, a nucleotide binding domain
that is either NACHT or NB-ARC and a C-terminal domain that is an Ankyrin repeat,
a WD domain or a TPR domain)
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
ap.add_argument('--FoL_genome',required=True,type=str,help='The genome assembly of FoL 4287')
ap.add_argument('--summary_table',required=True,type=str,help='The summary table from Fo)build_annot_table.py')
ap.add_argument('--N_term_out',required=True,type=str,help='The output file detailing N terminal domains by chromosome')
ap.add_argument('--NBD_out',required=True,type=str,help='The output file detailing NBD domains by chromosome')
ap.add_argument('--C_term_out',required=True,type=str,help='The output file detailing C terminal domains by chromosome')

conf = ap.parse_args()

file_N_term_out = open(conf.N_term_out, 'w')
file_NBD_out = open(conf.NBD_out, 'w')
file_C_term_out = open(conf.C_term_out, 'w')


with open(conf.summary_table) as f:
    summary_table_lines = f.readlines()

#-----------------------------------------------------
# Step 1.5
# Build a dictionary detailing all the annotations
# that will be pulled out from interproscan annotations
#-----------------------------------------------------


pfam_domain_dict = {
# N-terminal domains:
'PF00168' : "C2",
'PF17109' : "Goodbye",
'PF14479': "Helo",
'PF17111' : "Helo-like",
'PF06985' : "HET",
'PF01734' : "Patatin",
'PF00082' : "Peptidase S8",
'PF11558' : "PFD",
'PF00069' : "PKinase",
'PF01048' : "PND UDP",
'PF04607' : "RelA Spot",
'PF17107' : 'SesA',
'PF17046' : 'sesB-like',
# Nucleotide binding domains:
'PF00931' : 'NB-ARC',
'PF05729' : 'NACHT',
'SSF52540' : 'P-loop NTPase', # P-loop NTPase superfamily (CL0023)
'IPR027417' : 'P-loop NTPase',
'PF13401' : 'AAA',
# C-terminal repeat domains:
"SSF48452" : 'TPR', # TPR superfamily (CL0020)
"PF00515" : 'TPR', # TPR1
"PF07719" : 'TPR', # TPR2
"PF07720" : 'TPR', # TPR3
"PF07721" : 'TPR', # TPR4
"PF12688" : 'TPR', # TPR5
"PF13174" : 'TPR', # TPR6
"PF13176" : 'TPR', # TPR7
"PF13181" : 'TPR', # TPR8
'PF13371' : 'TPR', # TPR9
'PF13374' : 'TPR', # TPR10
"PF13414" : 'TPR', # TPR11
"PF13424" : 'TPR', # TPR12
# "PF07719" : 'TPR', # Can't find entry for TPR_13
"PF13428" : 'TPR', # TPR14
"PF13429" : 'TPR', # TPR15
"PF13432" : 'TPR', # TPR16
'PF13431' : 'TPR', # TPR17
"PF13512" : 'TPR', # TPR18
'PF14559' : 'TPR', # TPR19
"PF14561" : 'TPR', # TPR20
"PF09976" : 'TPR', # TPR21
"PS50005" : 'TPR', # Prosite profile for TPR
'IPR013026' : 'TPR', # Interproscan annotation for TPR
'SM00386' : 'half TPR',
'SSF48403' : 'Ankyrin', # Ankyrin superfamily (CL0465)
'PF12796' : 'Ankyrin',
'SSF50978' : 'WD40-like', # WD40-like superfamily (no code)
'PF00400' : 'WD40'
}

#-----------------------------------------------------
# Step 2
# Extract interproscan information on NLR candidate
# proteins
# This section also records the number of genes on
# chromsome / region
#-----------------------------------------------------

NLR_prot_dict = defaultdict(list)
contig_dict = defaultdict(list)
gene_count_dict = defaultdict(int)
i = 0
for line in summary_table_lines:
    line = line.rstrip()
    split_line = line.split("\t")
    gene_id = split_line[0]
    if len(split_line) >= 34:
        interpro_info = split_line[33]
        orthogroup = split_line[15]
        orthogroup_content = split_line[17]
        Fol_contig = split_line[1]
        gene_start = split_line[2]
        gene_end = split_line[3]

        FoL_Chr1_LS_start = int(5963944)
        FoL_Chr2_LS_start = int(4903160)
        if str(1) == Fol_contig:
            if (int(gene_end) <= FoL_Chr1_LS_start):
                Fol_contig = "1_core"
            else:
                Fol_contig = "1_LS"
        elif str(2) == Fol_contig:
            if (int(gene_end) <= FoL_Chr2_LS_start):
                Fol_contig = "2_core"
            else:
                Fol_contig = "2_LS"
        gene_count_dict[Fol_contig] += 1

        N_terminal_list = ['PF00168', 'PF17109', 'PF14479', 'PF17111', 'PF06985', 'PF01734', 'PF00082', 'PF11558', 'PF00069', 'PF01048', 'PF04607', 'PF17046', 'PF17107']
        NBD_list = ['SSF52540', 'IPR027417', 'PF13401', 'PF00931', 'PF05729'] # The P-loop containing protein is tested first so that NACHT or NBARC annotations will replace it
        C_terminal_list = [
         'SSF48452', 'SSF48403', 'SSF50978',
         'SM00386', "PF00515", "PF07719", "PF07720",
         "PF07721", "PF12688", "PF13174", "PF13176", "PF13181", 'PF13371',
         'PF13374', "PF13414", "PF13424", "PF13428", "PF13429", "PF13432",
         'PF13431', "PF13512", 'PF14559', "PF14561", "PF09976", "PS50005",
         'IPR013026', 'PF12796',
         'PF00400'
         ]
        N_terminal = 'NA'
        NBD = 'NA'
        C_terminal = 'NA'
        for x in N_terminal_list:
            if x in interpro_info:
                domain = pfam_domain_dict[x]
                N_terminal = domain
        for x in NBD_list:
            if x in interpro_info:
                domain = pfam_domain_dict[x]
                NBD = domain
        for x in C_terminal_list:
            if x in interpro_info:
                domain = pfam_domain_dict[x]
                C_terminal = domain
        combination = N_terminal + "---" + NBD + "---" + C_terminal
        store_line = gene_id + "\t" + orthogroup + "\t" + orthogroup_content

        if any( y != 'NA' for y in [N_terminal, NBD, C_terminal] ):
            NLR_prot_dict[combination].append(store_line)
            i += 1
            contig_dict[Fol_contig].append((combination + "\t" + store_line))


#-----------------------------------------------------
# Step 2
# Identify the length and gene density of FoL chromosomes
#-----------------------------------------------------

contig_length_dict = defaultdict(list)

genome_file = open(conf.FoL_genome, 'r')
for cur_record in SeqIO.parse(genome_file,"fasta"):
    seq_id = cur_record.id
    if seq_id == '1':
        chr1_core = cur_record.seq[0:FoL_Chr1_LS_start]
        chr1_LS = cur_record.seq[FoL_Chr1_LS_start:]
        seq_len_core = len(chr1_core)
        seq_len_LS = len(chr1_LS)
        N_count_core = chr1_core.count('N')
        N_count_LS = chr1_LS.count('N')
        seq_len_no_Ns_core = int(seq_len_core) - int(N_count_core)
        seq_len_no_Ns_LS = int(seq_len_LS) - int(N_count_LS)
        gene_density_core = numpy.divide(float(gene_count_dict["1_core"]), float(seq_len_core)) * 1000000
        gene_density_LS = numpy.divide(float(gene_count_dict["1_LS"]), float(seq_len_LS)) * 1000000
        gene_density_no_Ns_core = numpy.divide(float(gene_count_dict["1_core"]), float(seq_len_no_Ns_core)) * 1000000
        gene_density_no_Ns_LS = numpy.divide(float(gene_count_dict["1_LS"]), float(seq_len_no_Ns_LS)) * 1000000

        contig_length_dict["1_core"].extend([seq_len_core, seq_len_no_Ns_core, gene_density_core, gene_density_no_Ns_core])
        contig_length_dict["1_LS"].extend([seq_len_LS, seq_len_no_Ns_LS, gene_density_LS, gene_density_no_Ns_LS])
    elif seq_id == '2':
        chr2_core = cur_record.seq[0:FoL_Chr2_LS_start]
        chr2_LS = cur_record.seq[FoL_Chr2_LS_start:]
        seq_len_core = len(chr2_core)
        seq_len_LS = len(chr2_LS)
        N_count_core = chr2_core.count('N')
        N_count_LS = chr2_LS.count('N')
        seq_len_no_Ns_core = int(seq_len_core) - int(N_count_core)
        seq_len_no_Ns_LS = int(seq_len_LS) - int(N_count_LS)
        # contig_length_dict["2_core"].extend([seq_len_core, seq_len_no_Ns_core])
        # contig_length_dict["2_LS"].extend([seq_len_LS, seq_len_no_Ns_LS])
        gene_density_core = numpy.divide(float(gene_count_dict["2_core"]), float(seq_len_core)) * 1000000
        gene_density_LS = numpy.divide(float(gene_count_dict["2_LS"]), float(seq_len_LS)) * 1000000
        gene_density_no_Ns_core = numpy.divide(float(gene_count_dict["2_core"]), float(seq_len_no_Ns_core)) * 1000000
        gene_density_no_Ns_LS = numpy.divide(float(gene_count_dict["2_LS"]), float(seq_len_no_Ns_LS)) * 1000000

        contig_length_dict["2_core"].extend([seq_len_core, seq_len_no_Ns_core, gene_density_core, gene_density_no_Ns_core])
        contig_length_dict["2_LS"].extend([seq_len_LS, seq_len_no_Ns_LS, gene_density_LS, gene_density_no_Ns_LS])
    else:
        seq_len = len(cur_record.seq)
        N_count = cur_record.seq.count('N')
        seq_len_no_Ns = int(seq_len) - int(N_count)
        # contig_length_dict[seq_id].extend([seq_len, seq_len_no_Ns])
        # print (str(gene_count_dict[seq_id]) + " - " + str(seq_len))
        gene_density = numpy.divide(float(gene_count_dict[seq_id]), float(seq_len)) * 1000000
        # print (gene_density)
        gene_density_no_Ns = numpy.divide(float(gene_count_dict[seq_id]), float(seq_len_no_Ns)) * 1000000
        # print (gene_density_no_Ns)
        # print(seq_id)
        contig_length_dict[seq_id].extend([seq_len, seq_len_no_Ns, gene_density, gene_density_no_Ns])





#-----------------------------------------------------
# Step 3
# Print genes containing NLR domains
#-----------------------------------------------------

for k in sorted(NLR_prot_dict, key=lambda k: len(NLR_prot_dict[k]), reverse=True):
    num_genes = len(NLR_prot_dict[k])
    print k + "\t" + str(num_genes)

print( "\n\nthe total number of NLR genes was:\t" + str(i) )

#-----------------------------------------------------
# Step 4
# Print genes by domain and by chromosome
#-----------------------------------------------------

print("\n\nthe following genes were found in FoL regions:")

# for each chromsome in FoL
for FoL_contig in sorted(contig_dict):
    # Get stats on the number of genes containing NLR domains
    genes_in_contig = contig_dict[FoL_contig]
    total_genes = gene_count_dict[FoL_contig]
    # print(FoL_contig)
    # print(contig_length_dict[FoL_contig])
    contig_lgth = int(contig_length_dict[FoL_contig][0])
    gene_density = int(contig_length_dict[FoL_contig][2])
    gene_density_adj = int(contig_length_dict[FoL_contig][3])
    print("\n\nFoL contig " + str(FoL_contig) + " (" + str(len(genes_in_contig)) + "/ " + str(total_genes) +" genes; with a density of " + str(gene_density) + " genes.Mb-1 or " + str(gene_density_adj) + " adjusted)")
    file_N_term_out.write("\n\nFoL contig " + str(FoL_contig) + " (" + str(len(genes_in_contig)) + "/ " + str(total_genes) +" genes; with a density of " + str(gene_density) + " genes.Mb-1 or " + str(gene_density_adj) + " adjusted)\n")
    file_NBD_out.write("\n\nFoL contig " + str(FoL_contig) + " (" + str(len(genes_in_contig)) + "/ " + str(total_genes) +" genes; with a density of " + str(gene_density) + " genes.Mb-1 or " + str(gene_density_adj) + " adjusted)\n")
    file_C_term_out.write("\n\nFoL contig " + str(FoL_contig) + " (" + str(len(genes_in_contig)) + "/ " + str(total_genes) +" genes; with a density of " + str(gene_density) + " genes.Mb-1 or " + str(gene_density_adj) + " adjusted)\n")
    # For each chromosome build a dictionary for all the N-terminal, NBD and C-terminal domains
    N_terminal_dict = defaultdict(list)
    NBD_dict = defaultdict(list)
    C_terminal_dict = defaultdict(list)
    for line in genes_in_contig:
        split_line = line.split("\t")
        combination = split_line[0].split("---")
        N_terminal = combination[0]
        NBD = combination[1]
        C_terminal = combination[2]
        if N_terminal != 'NA':
            N_terminal_dict[N_terminal].append(line)
        if NBD != 'NA':
            NBD_dict[NBD].append(line)
        if C_terminal != 'NA':
            C_terminal_dict[C_terminal].append(line)
    # extracting N-terminal domains
    for keys in ["PKinase", "HET", "Peptidase S8", "PND UDP", "Helo-like", "Helo", "C2", "Patatin", "Goodbye", "PFD", "RelA Spot", 'SesA', 'sesB-like']:
        num_genes = len(N_terminal_dict[keys])
        frac_genes = numpy.divide(float(num_genes), float(total_genes)) * 100
        genes_per_Mb = numpy.divide(float(num_genes), float(contig_lgth)) * 1000000
        frac_density_adj = numpy.divide(float(num_genes), float(gene_density_adj))
        file_N_term_out.write("\t" + "\t".join([FoL_contig, str(keys), str(num_genes), str(frac_genes), str(frac_density_adj), str(genes_per_Mb)]) + "\n")
    # extracting NBD domains
    for keys in ['P-loop NTPase', 'NACHT', 'NB-ARC', 'AAA']:
        num_genes = len(NBD_dict[keys])
        frac_genes = numpy.divide(float(num_genes), float(total_genes)) * 100
        genes_per_Mb = numpy.divide(float(num_genes), float(contig_lgth)) * 1000000
        frac_density_adj = numpy.divide(float(num_genes), float(gene_density_adj))
        file_NBD_out.write("\t" + "\t".join([FoL_contig, str(keys), str(num_genes), str(frac_genes), str(frac_density_adj), str(genes_per_Mb)]) + "\n")
    # extracting C-terminal domains
    for keys in ['Ankyrin', 'TPR', 'WD40', 'WD40-like', 'half TPR']:
        num_genes = len(C_terminal_dict[keys])
        frac_genes = numpy.divide(float(num_genes), float(total_genes)) * 100
        genes_per_Mb = numpy.divide(float(num_genes), float(contig_lgth)) * 1000000
        frac_density_adj = numpy.divide(float(num_genes), float(gene_density_adj))
        file_C_term_out.write("\t" + "\t".join([FoL_contig, str(keys), str(num_genes), str(frac_genes), str(frac_density_adj), str(genes_per_Mb)]) + "\n")

    del N_terminal_dict
    del NBD
    del C_terminal_dict
