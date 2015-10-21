
cd /home/groups/harrisonlab/project_files/fusarium/timecourse/tmp
samtools view -bS ../4fus.sam > 24_hr_alignment.bam
samtools sort 24_hr_alignment.bam 24_hr_alignment_sorted
samtools index 24_hr_alignment_sorted.bam
cp ../../repeat_masked/F.oxysporum_fsp_cepae/Fus2/Fus2_combined_49/Fus2_combined_49_contigs_hardmasked.fa
samtools faidx Fus2_combined_49_contigs_hardmasked.fa
samtools tview 24_hr_alignment_sorted.bam Fus2_combined_49_contigs_hardmasked.fa
cp ../../gene_pred/augustus/F.oxysporum_fsp_cepae/Fus2/Fus2_aug_out.gff
bedtools intersect -c -a Fus2_aug_out.gff -b 24_hr_alignment_sorted.bam > 24hr_intersect3.bed
cat 24hr_intersect3.bed | grep -E -v "\s0$" | grep -w "gene" > Fus2_expressed_genes.bed
cat Fus2_expressed_genes.bed | sort -r -n -k 10 | less
# This top lines of this file are:
# NODE_9700_length_34512_cov_86.442451    AUGUSTUS        gene    9014    11091   0.32    -       .       g13840  2629
# NODE_11428_length_3064_cov_2466.822754  AUGUSTUS        gene    1       2042    0.06    +       .       g15230  1759
# NODE_1097_length_595114_cov_83.826759   AUGUSTUS        gene    518067  521767  0.17    +       .       g2042   413
# NODE_11451_length_942_cov_2469.252686   AUGUSTUS        gene    248     613     0.72    -       .       g15334  372
# NODE_256_length_160675_cov_85.558762    AUGUSTUS        gene    86821   88592   0.36    +       .       g8812   309
# NODE_149_length_43877_cov_84.839462     AUGUSTUS        gene    20431   21783   0.96    +       .       g13469  295
# NODE_956_length_295543_cov_83.945419    AUGUSTUS        gene    250771  251290  0.99    -       .       g5661   241
# NODE_11459_length_659_cov_2500.394531   AUGUSTUS        gene    1       203     0.65    -       .       g15359  207
# NODE_724_length_250420_cov_83.169510    AUGUSTUS        gene    209007  210676  0.83    +       .       g6901   199
# NODE_8489_length_811399_cov_84.416985   AUGUSTUS        gene    296996  298062  0.94    +       .       g891    165
# NODE_4583_length_445765_cov_84.597122   AUGUSTUS        gene    274107  275696  1       +       .       g3043   134
# NODE_11477_length_312246_cov_83.660294  AUGUSTUS        gene    281911  283153  0.97    +       .       g5276   133
# NODE_4974_length_475916_cov_84.779022   AUGUSTUS        gene    137868  139358  0.39    -       .       g2735   114
# NODE_990_length_149944_cov_81.983368    AUGUSTUS        gene    121307  125898  0.61    +       .       g9290   107

# Of note the gene g13840 is highly expressed. We can look for its interposcan annotation
INTERPRO_DIR=/home/groups/harrisonlab/project_files/fusarium/gene_pred/interproscan/Fus2/split
cat $INTERPRO_DIR/Fus2_interproscan.gff3 | grep 'g13840.' | less

# 1. Results indicate that gene g13840 it is a transmembrane protein involved
# in peptide transport.
# g13840.t1       Pfam    protein_match   181     476     7.6E-27 +       .       Name=PF00854;signature_desc=POT family;Target=g13840.t1 181 476;status=T;ID=match$196_181_476;Ontology_term="GO:0005215","GO:0006810","GO:0016020";date=06-01-2015;Dbxref="InterPro:IPR000109","Reactome:REACT_15518","Reactome:REACT_19419"
# g13840.t1       SUPERFAMILY     protein_match   144     327     5.75E-13        +       .       Name=SSF103473;Target=g13840.t1 144 327;status=T;ID=match$201_144_327;date=06-01-2015;Dbxref="InterPro:IPR020846"
# 1. Visual inspection of g12840 shows all the alignments
# between contig positions 9331 and 9381
# This is probably part of an onion read.

# 2. g11428 appears to be a viral gene
# 2. Visual inspection shows a high level of heterogeneity
# between aligned reads and the FUS2 genome. Certainly many
# reads present covering the gene and the neighbouring region.

# 3. g2042 is a ubiquitin family protein. Which is foundubiquitously in Eukaryotes.
# - homology to onion?
# Visual inspection of the reads indicate that there is high
# heterogeneity between the aligned reads and the Fus2 genome.

# 4. g15334 is a viral capsid protein.
# Visual inspection shows that thsis gene is larger than predicted by augustus.
# There is are some notable differences between the aligned reads and Fus2 genome.
# There is apporximately 10 genomic 'errors' in the fist 983bp of the contig.

for TIMEPOINT in 0 4 8 16 24.1 24.2 36 48 72 96; do
	mkdir -p expreiment1/F.oxysporum_fsp_cepae/Fus2/$TIMEPOINT
done

cp timecourse/expreiment1/1_R1_fus_q30_l50.fastq timecourse/expreiment1/F.oxysporum_fsp_cepae/Fus2/0/.
cp timecourse/expreiment1/4hr-root_R1_fus_q30_l50.fastq timecourse/expreiment1/F.oxysporum_fsp_cepae/Fus2/4/.
cp timecourse/expreiment1/8hr-root_R1_fus_q30_l50.fastq timecourse/expreiment1/F.oxysporum_fsp_cepae/Fus2/8/.
cp timecourse/expreiment1/3_fus_R1_fus_q30_l50.fastq timecourse/expreiment1/F.oxysporum_fsp_cepae/Fus2/16/.
cp timecourse/expreiment1/4_R1_fus_q30_l50.fastq timecourse/expreiment1/F.oxysporum_fsp_cepae/Fus2/24.1/.
cp timecourse/expreiment1/24hr-root_R1_fus_q30_l50.fastq timecourse/expreiment1/F.oxysporum_fsp_cepae/Fus2/24.2/.
cp timecourse/expreiment1/36h-root_R1_fus_q30_l50.fastq timecourse/expreiment1/F.oxysporum_fsp_cepae/Fus2/36/.
cp timecourse/expreiment1/6_R1_fus_q30_l50.fastq timecourse/expreiment1/F.oxysporum_fsp_cepae/Fus2/48/.
cp timecourse/expreiment1/7_R1_fus_q30_l50.fastq timecourse/expreiment1/F.oxysporum_fsp_cepae/Fus2/72/.
cp timecourse/expreiment1/8_R1_fus_q30_l50.fastq timecourse/expreiment1/F.oxysporum_fsp_cepae/Fus2/96/.

cp timecourse/expreiment1/1fus.sam timecourse/expreiment1/F.oxysporum_fsp_cepae/Fus2/0/.
cp timecourse/expreiment1/4hfus.sam timecourse/expreiment1/F.oxysporum_fsp_cepae/Fus2/4/.
cp timecourse/expreiment1/8hfus.sam timecourse/expreiment1/F.oxysporum_fsp_cepae/Fus2/8/.
cp timecourse/expreiment1/3fus.sam timecourse/expreiment1/F.oxysporum_fsp_cepae/Fus2/16/.
cp timecourse/expreiment1/4fus.sam timecourse/expreiment1/F.oxysporum_fsp_cepae/Fus2/24.1/.
cp timecourse/expreiment1/24hfus.sam timecourse/expreiment1/F.oxysporum_fsp_cepae/Fus2/24.2/.
cp timecourse/expreiment1/36hfus.sam timecourse/expreiment1/F.oxysporum_fsp_cepae/Fus2/36/.
cp timecourse/expreiment1/6fus.sam timecourse/expreiment1/F.oxysporum_fsp_cepae/Fus2/48/.
cp timecourse/expreiment1/7fus.sam timecourse/expreiment1/F.oxysporum_fsp_cepae/Fus2/72/.
cp timecourse/expreiment1/8fus.sam timecourse/expreiment1/F.oxysporum_fsp_cepae/Fus2/96/.


# qsub /home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq/gene_expression.sh timecourse/expreiment1/F.oxysporum_fsp_cepae/Fus2/4/4hfus.sam repeat_masked/F.oxysporum_fsp_cepae/Fus2/Fus2_combined_49/Fus2_combined_49_contigs_hardmasked.fa gene_pred/augustus/F.oxysporum_fsp_cepae/Fus2/Fus2_aug_out.gff

for IN_SAM in $(ls timecourse/expreiment1/F.oxysporum_fsp_cepae/Fus2/*/*.sam); do
	SCRIPT_DIR=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
	GENOME=repeat_masked/F.oxysporum_fsp_cepae/Fus2/Fus2_combined_49/Fus2_combined_49_contigs_hardmasked.fa
	PRED_GENES=gene_pred/augustus/F.oxysporum_fsp_cepae/Fus2/Fus2_aug_out.gff
	echo $IN_SAM
	qsub $SCRIPT_DIR/gene_expression.sh $IN_SAM $GENOME $PRED_GENES
done


for FILE in $(ls timecourse/expreiment1/F.oxysporum_fsp_cepae/Fus2/*/*expressed_genes.bed); do
	printf "$FILE\n" >> tmp3/top_50_genes.csv
	cat $FILE | head -n50 >> tmp3/top_50_genes.csv
done

# identify is any of the top 100 expressed genes in Fus2 come from the same contig.
cat timecourse/expreiment1/F.oxysporum_fsp_cepae/Fus2/96/Fus2_96_expressed_genes.bed | head -n 100 | cut -f1 | sort | uniq -c | sort -r | less
# 	  4 NODE_9643_length_1511754_cov_84.114418
#       4 NODE_4583_length_445765_cov_84.597122
#       4 NODE_1105_length_605632_cov_83.266510
#       3 NODE_566_length_732593_cov_83.236145
#       3 NODE_2374_length_16459_cov_83.673859
#       3 NODE_1097_length_595114_cov_83.826759
#       2 NODE_956_length_295543_cov_83.945419
#       2 NODE_731_length_218737_cov_86.751648
#       2 NODE_7096_length_311910_cov_84.936836
#       2 NODE_6283_length_399598_cov_84.774551
#       2 NODE_608_length_519399_cov_83.445702
#       2 NODE_4974_length_475916_cov_84.779022
#       2 NODE_4873_length_67562_cov_92.001808
#       2 NODE_4575_length_487109_cov_84.320221
#       2 NODE_256_length_160675_cov_85.558762
#       2 NODE_2122_length_266153_cov_81.544640
#       2 NODE_1420_length_454186_cov_84.973450
#       2 NODE_10255_length_363655_cov_83.498451
#       1 NODE_979_length_171603_cov_83.119377

cat analysis/blast_homology/F.oxysporum_fsp_cepae/Fus2/Fus2_six-appended_parsed.fa_homologs.csv | cut -f 1,9,14,15 | sort -k2 | grep 'NODE' | cat
# Fusarium_oxysporum_f._sp._lycopersici_SIX3_gene_for_Secreted_in_xylem_3_protein	NODE_1364_length_4671_cov_108.113892	1100	1657
# Fusarium_oxysporum_f._sp._lycopersici_isolate_14844_secreted_in_xylem_3_(SIX3)_gene,_complete_cds	NODE_1364_length_4671_cov_108.113892	1136	1624
# Fusarium_oxysporum_f._sp._lycopersici_isolate_BFOL-51_secreted_in_xylem_3_(SIX3)_gene,_complete_cds	NODE_1364_length_4671_cov_108.113892	1136	1624
# Fusarium_oxysporum_f._sp._lycopersici_isolate_FOL-MM10_secreted_in_xylem_3_(SIX3)_gene,_complete_cds	NODE_1364_length_4671_cov_108.113892	1136	1624
# Fusarium_oxysporum_f._sp._lycopersici_isolate_IPO3_secreted_in_xylem_3_(SIX3)_gene,_complete_cds	NODE_1364_length_4671_cov_108.113892	1136	1624
# Fusarium_oxysporum_f._sp._lycopersici_isolate_BFOL-51_secreted_in_xylem_5_(SIX5)_gene,_partial_cds	NODE_2080_length_2535_cov_50.341618	310	836
# Fusarium_oxysporum_f._sp._lycopersici_Six5_mRNA,_complete_cds	NODE_2080_length_2535_cov_50.341618	706	843
# Fusarium_oxysporum_f._sp._lilii_isolate_NRRL_28395_secreted_in_xylem_7-like_protein_(SIX7)_gene,_complete_cds	NODE_2938_length_1682_cov_52.927467	601	1266
# Fusarium_oxysporum_f._sp._lycopersici_isolate_BFOL-51_secreted_in_xylem_7_(SIX7)_gene,_complete_cds	NODE_2938_length_1682_cov_52.927467	601	1266
# Fusarium_oxysporum_f._sp._lycopersici_secreted_in_xylem_Six7_(SIX7)_mRNA,_complete_cds	NODE_2938_length_1682_cov_52.927467	601	1266
# Fusarium_oxysporum_f._sp._lycopersici_strain_Fol007_Six10_(SIX10)_mRNA,_complete_cds	NODE_788_length_4929_cov_54.397850	4328	4697
# Fusarium_oxysporum_f._sp._lycopersici_strain_Fol007_Six12_(SIX12)_mRNA,_complete_cds	NODE_788_length_4929_cov_54.397850	750	1091
#
#In Summary:
# Fusarium_oxysporum_f._sp._lycopersici_isolate_14844_secreted_in_xylem_3_(SIX3)_gene,_complete_cds
# NODE_1364_length_4671_cov_108.113892	1136	1624
# Fusarium_oxysporum_f._sp._lycopersici_isolate_BFOL-51_secreted_in_xylem_5_(SIX5)_gene,_partial_cds
# NODE_2080_length_2535_cov_50.341618	310	836
# Fusarium_oxysporum_f._sp._lilii_isolate_NRRL_28395_secreted_in_xylem_7-like_protein_(SIX7)_gene,_complete_cds
# NODE_2938_length_1682_cov_52.927467	601	1266
# Fusarium_oxysporum_f._sp._lycopersici_strain_Fol007_Six10_(SIX10)_mRNA,_complete_cds
# NODE_788_length_4929_cov_54.397850	4328	4697
# Fusarium_oxysporum_f._sp._lycopersici_strain_Fol007_Six12_(SIX12)_mRNA,_complete_cds	NODE_788_length_4929_cov_54.397850	750	1091

#Identify expression of Fusarium SIX genes.
qsub /home/armita/git_repos/emr_repos/tools/pathogen/blast/blast_pipe.sh analysis/blast_homology/six_genes/six-appended_parsed.fa dna repeat_masked/F.oxysporum_fsp_cepae/Fus2/Fus2_combined_49/Fus2_combined_49_contigs_unmasked.fa
/home/armita/git_repos/emr_repos/tools/pathogen/blast/blast2gff.pl six_gene analysis/blast_homology/F.oxysporum_fsp_cepae/Fus2/Fus2_six-appended_parsed.fa_homologs.csv > analysis/blast_homology/F.oxysporum_fsp_cepae/Fus2/Fus2_six-appended_parsed.fa_homologs.gff
cat timecourse/expreiment1/F.oxysporum_fsp_cepae/Fus2/96/Fus2_96_expressed_genes.bed | grep 'NODE_4546' | less
bedtools intersect -c -a analysis/blast_homology/F.oxysporum_fsp_cepae/Fus2/Fus2_six-appended_parsed.fa_homologs.gff -b timecourse/expreiment1/F.oxysporum_fsp_cepae/Fus2/96/Fus2_96_sorted.bam > timecourse/expreiment1/F.oxysporum_fsp_cepae/Fus2/96/Fus2_96_expressed_six.bed
