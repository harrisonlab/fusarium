#!/bin/bash

# Commands used to perform comparative genomics in Fusarium

ls -d gene_pred/augustus/F.oxysporum_fsp_cepae/*
# gene_pred/augustus/F.oxysporum_fsp_cepae/125  gene_pred/augustus/F.oxysporum_fsp_cepae/D2
# gene_pred/augustus/F.oxysporum_fsp_cepae/55   gene_pred/augustus/F.oxysporum_fsp_cepae/Fus2
# gene_pred/augustus/F.oxysporum_fsp_cepae/A23  gene_pred/augustus/F.oxysporum_fsp_cepae/HB17
# gene_pred/augustus/F.oxysporum_fsp_cepae/A28  gene_pred/augustus/F.oxysporum_fsp_cepae/PG

set -- 125 55 A23 A28 D2 Fus2 HB17 PG
for a; do 
	shift
	for b; do
	printf "%s - %s\n" "$a" "$b"
	qsub /home/armita/git_repos/emr_repos/tools/pathogen/orthology/inparanoid/sub_inparanoid.sh gene_pred/augustus/F.oxysporum_fsp_cepae/$a/"$a"_augustus_preds.aa gene_pred/augustus/F.oxysporum_fsp_cepae/$b/"$b"_augustus_preds.aa gene_pred/augustus/F.oxysporum_fsp_cepae/$a/"$a"_augustus_preds.gtf gene_pred/augustus/F.oxysporum_fsp_cepae/$b/"$b"_augustus_preds.gtf
	done 
done

# Commands not Used yet:
# 
mkdir analysis/inparanoid/summary_tables
# cat analysis/inparanoid/*/*_seqs.txt | cut -f1 | sort | uniq | less
cat analysis/inparanoid/*/*_seqs.txt | cut -f1 | sort | uniq > analysis/inparanoid/summary_tables/all_genes.txt
for FILEZ in $(ls analysis/inparanoid/*/sqltable.*); do 
	FILE_LIST="$FILE_LIST $FILEZ"
done

/home/armita/git_repos/emr_repos/tools/pathogen/orthology/inparanoid/build_orthology_tab.pl analysis/inparanoid/summary_tables/all_genes.txt $FILE_LIST > analysis/inparanoid/summary_tables/orthology_tab.csv
mkdir -p analysis/inparanoid/summary_tables/gene_orthologs
cd analysis/inparanoid/summary_tables/gene_orthologs
/home/armita/git_repos/emr_repos/tools/pathogen/orthology/inparanoid/analyse_orthology_tab.pl ../../../../analysis/inparanoid/summary_tables/orthology_tab.csv -file ../../../../analysis/inparanoid/summary_tables/all_genes.txt -deep -print_once
cd ../../../../
/home/armita/git_repos/emr_repos/tools/pathogen/orthology/inparanoid/build_final_ortho_tab.pl analysis/inparanoid/summary_tables/gene_orthologs/ 125 55 A23 A28 D2 Fus2 HB17 PG > analysis/inparanoid/summary_tables/FoC_gene_orthologs.csv
cat analysis/inparanoid/summary_tables/FoC_gene_orthologs.csv| sed 's/\+/1/g' | sed 's/\-/0/g' > analysis/inparanoid/summary_tables/FoC_gene_orthologs2.csv

# The number of inparalogous genes was recorded for each strain
printf "" > analysis/inparanoid/summary_tables/inparalog_groups.txt
for FILEZ in $(ls analysis/inparanoid/summary_tables/gene_orthologs/*); do
THISFILE=$(basename $FILEZ)
printf "$THISFILE\n" >> analysis/inparanoid/summary_tables/inparalog_groups.txt
cat $FILEZ | cut -f1 -d"|" | sort | uniq -c >> analysis/inparanoid/summary_tables/inparalog_groups.txt
done
cat analysis/inparanoid/summary_tables/inparalog_groups.txt | grep -v '.txt' | grep -v '1 ' | grep '125' | wc -l
cat analysis/inparanoid/summary_tables/inparalog_groups.txt | grep -v '.txt' | grep -v '1 ' | grep '55' | wc -l
cat analysis/inparanoid/summary_tables/inparalog_groups.txt | grep -v '.txt' | grep -v '1 ' | grep 'A23' | wc -l
cat analysis/inparanoid/summary_tables/inparalog_groups.txt | grep -v '.txt' | grep -v '1 ' | grep 'A28' | wc -l
cat analysis/inparanoid/summary_tables/inparalog_groups.txt | grep -v '.txt' | grep -v '1 ' | grep 'D2' | wc -l
cat analysis/inparanoid/summary_tables/inparalog_groups.txt | grep -v '.txt' | grep -v '1 ' | grep 'Fus2' | wc -l
cat analysis/inparanoid/summary_tables/inparalog_groups.txt | grep -v '.txt' | grep -v '1 ' | grep 'HB17' | wc -l
cat analysis/inparanoid/summary_tables/inparalog_groups.txt | grep -v '.txt' | grep -v '1 ' | grep 'PG' | wc -l

# Venn diagrams were made to show orthology between genomes

R
library(VennDiagram, lib.loc="/home/armita/R-packages/")
data1 <- read.delim(file="analysis/inparanoid/summary_tables/FoC_gene_orthologs2.csv", header=T, sep="\t")
attach(data1)
colnames(data1) <- c("ortho_group", "gene_name", "strain1", "strain2", "strain3", "strain4", "strain5", "strain6", "strain7", "strain8")
name1<-"125"
name2<-"55"
name3<-"A23"
name4<-"A28"
name5<-"D2"
name6<-"Fus2"
name7<-"HB17"
name8<-"PG"
area1=sum(data1$strain1)
area2=sum(data1$strain2)
area3=sum(data1$strain3)
area4=sum(data1$strain4)
area5=sum(data1$strain5)
area6=sum(data1$strain6)
area7=sum(data1$strain7)
area8=sum(data1$strain8)
area1
area2
area3
area4
area5
area6
area7
area8
head(subset(data1, (strain1==1) & (strain2==0|1) & (strain3==1) & (strain4==0) & (strain5==0) & (strain6==1) & (strain7==0|1) & (strain8==0)))
head(subset(data1, (strain1==1) & (strain3==1) & (strain4==0) & (strain5==0) & (strain6==1)))

nrow(subset(data1, strain1==1 & strain2==0|1 & strain3==1 & strain4==0 & strain5==0 & strain6==1 & strain7==0|1 & strain8==0))

# Analyse a subset of two pathogenic and two non-pathogenic isolates
/home/armita/git_repos/emr_repos/tools/pathogen/orthology/inparanoid/build_final_ortho_tab.pl analysis/inparanoid/summary_tables/gene_orthologs/ A23 A28 D2 Fus2 > analysis/inparanoid/summary_tables/A23_A28_D2_Fus2_orthologs.csv
cat analysis/inparanoid/summary_tables/A23_A28_D2_Fus2_orthologs.csv| sed 's/\+/1/g' | sed 's/\-/0/g' > analysis/inparanoid/summary_tables/A23_A28_D2_Fus2_orthologs2.csv


# Making a 4-way Venn Diagram:
# Open R
R
# These commabds were executed in R

library(VennDiagram, lib.loc="/home/armita/R-packages/")
data1 <- read.delim(file="analysis/inparanoid/summary_tables/A23_A28_D2_Fus2_orthologs2.csv", header=T, sep="\t")
attach(data1)

colnames(data1) <- c("ortho_group", "gene_name", "strain1", "strain2", "strain3", "strain4")
name1<-"A23"
name2<-"A28"
name3<-"D2"
name4<-"Fus2"
area1=sum(data1$strain1)
area2=sum(data1$strain2)
area3=sum(data1$strain3)
area4=sum(data1$strain4)
label1 <- paste(name1, ' (', area1, ')', sep ="")
label2 <- paste(name2, ' (', area2, ')', sep ="")
label3 <- paste(name3, ' (', area3, ')', sep ="")
label4 <- paste(name4, ' (', area4, ')', sep ="")

n12=nrow(subset(data1, strain1==1 & strain2==1))
n13=nrow(subset(data1, strain1==1 & strain3))
n14=nrow(subset(data1, strain1==1 & strain4==1))
n23=nrow(subset(data1, strain2==1 & strain3==1))
n24=nrow(subset(data1, strain2==1 & strain4==1))
n34=nrow(subset(data1, strain3==1 & strain4==1))
n123=nrow(subset(data1, strain1==1 & strain2==1 & strain3==1))
n124=nrow(subset(data1, strain1==1 & strain2==1 & strain4==1))
n134=nrow(subset(data1, strain1==1 & strain3==1 & strain4==1))
n234=nrow(subset(data1, strain2==1 & strain3==1 & strain4==1))
n1234=nrow(subset(data1, strain1==1 & strain2==1 & strain3==1 & strain4==1))

pdf('analysis/inparanoid/summary_tables/phytoph_ortholog_venn.pdf')

draw.quad.venn(area1, area2, area3, area4, n12, n13, n14, n23, n24, n34, n123, n124, n134, n234, n1234, category = c(label1, label2, label3, label4), lwd = rep(2, 4), lty = rep("solid", 4), col = rep("black", 4), fill = NULL, alpha = rep(0.5, 4),
label.col = rep("black", 15), cex = rep(1, 15), fontface = rep("plain", 15), fontfamily = rep("serif", 15), cat.pos = c(-15, 15, 0, 0),
    cat.dist = c(0.22, 0.22, 0.11, 0.11), cat.col = rep("black", 4),
    cat.cex = rep(1, 4), cat.fontface = rep("plain", 4),
    cat.fontfamily = rep("serif", 4), cat.just = rep(list(c(0.5, 0.5)), 4),
    rotation.degree = 0, rotation.centre = c(0.5, 0.5), ind = TRUE)
dev.off()
q()

# Analyse a subset of two pathogenic and two non-pathogenic isolates
/home/armita/git_repos/emr_repos/tools/pathogen/orthology/inparanoid/build_final_ortho_tab.pl analysis/inparanoid/summary_tables/gene_orthologs/ A23 A28 PG Fus2 > analysis/inparanoid/summary_tables/A23_A28_PG_Fus2_orthologs.csv
cat analysis/inparanoid/summary_tables/A23_A28_PG_Fus2_orthologs.csv| sed 's/\+/1/g' | sed 's/\-/0/g' > analysis/inparanoid/summary_tables/A23_A28_PG_Fus2_orthologs2.csv


# Making a 4-way Venn Diagram:
# Open R
R
# These commabds were executed in R

library(VennDiagram, lib.loc="/home/armita/R-packages/")
data1 <- read.delim(file="analysis/inparanoid/summary_tables/A23_A28_PG_Fus2_orthologs2.csv", header=T, sep="\t")
attach(data1)

colnames(data1) <- c("ortho_group", "gene_name", "strain1", "strain2", "strain3", "strain4")
name1<-"A23"
name2<-"A28"
name3<-"PG"
name4<-"Fus2"
area1=sum(data1$strain1)
area2=sum(data1$strain2)
area3=sum(data1$strain3)
area4=sum(data1$strain4)
label1 <- paste(name1, ' (', area1, ')', sep ="")
label2 <- paste(name2, ' (', area2, ')', sep ="")
label3 <- paste(name3, ' (', area3, ')', sep ="")
label4 <- paste(name4, ' (', area4, ')', sep ="")

n12=nrow(subset(data1, strain1==1 & strain2==1))
n13=nrow(subset(data1, strain1==1 & strain3))
n14=nrow(subset(data1, strain1==1 & strain4==1))
n23=nrow(subset(data1, strain2==1 & strain3==1))
n24=nrow(subset(data1, strain2==1 & strain4==1))
n34=nrow(subset(data1, strain3==1 & strain4==1))
n123=nrow(subset(data1, strain1==1 & strain2==1 & strain3==1))
n124=nrow(subset(data1, strain1==1 & strain2==1 & strain4==1))
n134=nrow(subset(data1, strain1==1 & strain3==1 & strain4==1))
n234=nrow(subset(data1, strain2==1 & strain3==1 & strain4==1))
n1234=nrow(subset(data1, strain1==1 & strain2==1 & strain3==1 & strain4==1))

pdf('analysis/inparanoid/summary_tables/A23_A28_PG_Fus2_venn_diagram.pdf')

draw.quad.venn(area1, area2, area3, area4, n12, n13, n14, n23, n24, n34, n123, n124, n134, n234, n1234, category = c(label1, label2, label3, label4), lwd = rep(2, 4), lty = rep("solid", 4), col = rep("black", 4), fill = NULL, alpha = rep(0.5, 4),
label.col = rep("black", 15), cex = rep(1, 15), fontface = rep("plain", 15), fontfamily = rep("serif", 15), cat.pos = c(-15, 15, 0, 0),
    cat.dist = c(0.22, 0.22, 0.11, 0.11), cat.col = rep("black", 4),
    cat.cex = rep(1, 4), cat.fontface = rep("plain", 4),
    cat.fontfamily = rep("serif", 4), cat.just = rep(list(c(0.5, 0.5)), 4),
    rotation.degree = 0, rotation.centre = c(0.5, 0.5), ind = TRUE)
dev.off()
q()

# Neither of these Venn diagrams show a major difference in gene compliment between Fusarium genomes
# The distribution of six genes in the genomes shall be investigated further.

# Identifying which predicted genes relate to BLAST results of Six genes

bedtools intersect -wa -wb -a gene_pred/augustus/F.oxysporum_fsp_cepae/Fus2/Fus2_augustus_preds.gtf -b analysis/blast_homology/F.oxysporum_fsp_cepae/Fus2/Fus2_six-appended_parsed.fa_homologs.gff > six_intersect_Fus2_genes.bed
cat six_intersect_Fus2_genes.bed | grep -a1 'gene' | sort | uniq | less
# This shows that Fus2|g10783 relates to SIX10
# The distribution of SIX10 orthologs through isolates shall be identified.
for FILEZ in $(ls analysis/inparanoid/summary_tables/gene_orthologs/*.txt); do 
printf "$FILEZ\n" >> tmp.out; 
cat $FILEZ | grep 'Fus2|g10783' >> tmp.out; 
done
cat tmp.out | grep -B1 'Fus2|g10783' | less
# This shows that Fus2|g10783 is a member of the orthogroup in file:
# analysis/inparanoid/summary_tables/gene_orthologs/125_g10783.t1.txt
cat analysis/inparanoid/summary_tables/gene_orthologs/125_g10783.t1.txt | cut -f1 | wc -l
cat analysis/inparanoid/summary_tables/gene_orthologs/125_g10783.t1.txt | cut -f1 | cut -f1 -d '|' | sort | uniq -c
# This orthogroup contains 94 genes from all eight fusarium strains.
#      15 125
#      11 55
#      13 A23
#      14 A28
#       5 D2
#      13 Fus2
#      14 HB17
#       9 PG
# To identify how similar that the orthologous genes called by INPARANOID are
# This family was downloaded and aligned.
for PATHZ in $(ls -d gene_pred/augustus/F.oxysporum_fsp_cepae/*); do
STRAIN=$(echo $PATHZ | rev | cut -f1 -d '/' | rev)
INFILE=analysis/inparanoid/summary_tables/gene_orthologs/125_g10783.t1.txt
cat  "$INFILE" | cut -f1 | grep "$STRAIN" | cut -f2 -d'|' > analysis/inparanoid/summary_tables/Foc_SIX10/"$STRAIN"_SIX10_orthologs.txt
printf "" > analysis/inparanoid/summary_tables/Foc_SIX10/"$STRAIN"_SIX10_orthologs.fa 
while read line; do
cat gene_pred/augustus/F.oxysporum_fsp_cepae/"$STRAIN"/"$STRAIN"_augustus_preds.aa | sed -e 's/\(^>.*$\)/#\1#/' | tr -d "\r" | tr -d "\n" | sed -e 's/$/#/' | tr "#" "\n" | sed -e '/^$/d' | grep -w -A1 "$line"  >> analysis/inparanoid/summary_tables/Foc_SIX10/"$STRAIN"_SIX10_orthologs.fa
done<analysis/inparanoid/summary_tables/Foc_SIX10/"$STRAIN"_SIX10_orthologs.txt
done
# Sequences were downloaded to Geneious, Renamed and aligned.
