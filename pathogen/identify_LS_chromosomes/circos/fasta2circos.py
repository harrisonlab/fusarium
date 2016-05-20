#!/usr/bin/python

'''
This program is used to convert fasta files into input format for circos
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

conf = ap.parse_args()



#-----------------------------------------------------
# Step 2
# Identify the length and gene density of FoL chromosomes
#-----------------------------------------------------

FoL_Chr1_LS_start = int(5963944)
FoL_Chr2_LS_start = int(4903160)
i = 0

contig_length_dict = defaultdict(list)

genome_file = open(conf.FoL_genome, 'r')
for cur_record in SeqIO.parse(genome_file,"fasta"):
    seq_id = cur_record.id

    # if str(1) == Fol_contig:
    #     if (int(gene_end) <= FoL_Chr1_LS_start):
    #         Fol_contig = "1_core"
    #     else:
    #         Fol_contig = "1_LS"
    # elif str(2) == Fol_contig:
    #     if (int(gene_end) <= FoL_Chr2_LS_start):
    #         Fol_contig = "2_core"
    #     else:
    #         Fol_contig = "2_LS"

    # if seq_id == '1':
    #     chr1_core = cur_record.seq[0:FoL_Chr1_LS_start]
    #     chr1_LS = cur_record.seq[FoL_Chr1_LS_start:]
    #     seq_len_core = len(chr1_core)
    #     seq_len_LS = len(chr1_LS)
    #     contig_length_dict["1_core"].extend([seq_len_core, seq_len_no_Ns_core, gene_density_core, gene_density_no_Ns_core])
    #     contig_length_dict["1_LS"].extend([seq_len_LS, seq_len_no_Ns_LS, gene_density_LS, gene_density_no_Ns_LS])
    # elif seq_id == '2':
    #     chr2_core = cur_record.seq[0:FoL_Chr2_LS_start]
    #     chr2_LS = cur_record.seq[FoL_Chr2_LS_start:]
    #     outline = " ".join(["chr", "-", str(seq_id), "FoL" + str(seq_id), "0", str(seq_len)], "chr" + str(i))
    # else:
    #     seq_len = len(cur_record.seq)
    #     outline = " ".join(["chr", "-", str(seq_id), "FoL" + str(seq_id), "0", str(seq_len)], "chr" + str(i))
    i += 1
    seq_len = len(cur_record.seq)
    outline = " ".join(["chr", "-", str(seq_id), "FoL" + str(seq_id), "0", str(seq_len), "chr" + str(i)])
    print(outline)
