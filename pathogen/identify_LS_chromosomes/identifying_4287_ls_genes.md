

```bash
  Fo_Fo47_assembly=assembly/external_group/F.oxysporum/fo47/broad/fusarium_oxysporum_fo47_1_supercontigs.fasta
  Fo_Fo47_genes=assembly/external_group/F.oxysporum/fo47/broad/fusarium_oxysporum_fo47_1_proteins.fasta
  Fo_Fo47_gff=assembly/external_group/F.oxysporum/fo47/broad/fusarium_oxysporum_fo47_1_transcripts.gtf

  FoL_4287_assembly=assembly/external_group/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl/Fusarium_oxysporum.FO2.31.dna.chromosome.fa
  FoL_4287_genes=assembly/external_group/F.oxysporum_fsp_lycopersici/4287/Fusox1/Fusox1_GeneCatalog_proteins_20110522.aa.fasta
  FoL_4287_gff=assembly/external_group/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl/Fusarium_oxysporum.FO2.31.chr.gff3

  FoL_LS_gene_headers=

  OutDir=analysis/FoL_ls_genes
  mkdir -p $OutDir

for Chr in 3 6 14 15; do
echo "Extracting genes on FoL chromosome: $Chr"
cat $FoL_4287_gff | grep -P -w "^$Chr" | grep -w -P '\tgene\t' | cut -f9 | cut -f1 -d ';' | cut -f2 -d ':' > $OutDir/chr_"$Chr"_gene_headers.txt
rm -f $OutDir/chr_"$Chr"_gene_orthogroups.txt
echo "The following number of genes were on this chromosome:"
cat $OutDir/chr_"$Chr"_gene_headers.txt | wc -l
echo ""
for gene in $(cat $OutDir/chr_"$Chr"_gene_headers.txt); do
# echo "$gene"
Orthogroup=$(cat analysis/orthology/orthomcl/FoC_vs_Fo_vs_FoL/FoC_vs_Fo_vs_FoL_orthogroups.txt | grep "$gene" | cut -f1 -d ':')
echo "$Orthogroup" >> $OutDir/chr_"$Chr"_gene_orthogroups.txt
done
paste $OutDir/chr_"$Chr"_gene_headers.txt $OutDir/chr_"$Chr"_gene_orthogroups.txt > $OutDir/chr_"$Chr"_ls_genes.tab
# cat analysis/orthology/orthomcl/FoC_vs_Fo_vs_FoL/FoC_vs_Fo_vs_FoL_orthogroups.txt | grep -f $OutDir/chr_"$Chr"_gene_headers.txt > $OutDir/chr_"$Chr"_gene_orthogroups.txt
sed -r "s/FOXG/$Chr\tFOXG/g" -i $OutDir/chr_"$Chr"_ls_genes.tab
done
cat $OutDir/chr_*_ls_genes.tab | sort -n > $OutDir/chr_LS_genes.tab
cat $OutDir/chr_LS_genes.tab | cut -f1,3 | sort | uniq -c | sort -n -r > $OutDir/chr_LS_orthogroup_abundance.tab
```

```
  Extracting genes on FoL chromosome: 3
  The following number of genes were on this chromosome:
  1073

  Extracting genes on FoL chromosome: 6
  The following number of genes were on this chromosome:
  824

  Extracting genes on FoL chromosome: 14
  The following number of genes were on this chromosome:
  256

  Extracting genes on FoL chromosome: 15
  The following number of genes were on this chromosome:
  451
```

```r
mydata <- read.table("analysis/FoL_ls_genes/chr_LS_orthogroup_abundance.csv", header=TRUE, sep=',')
# Stacked Bar Plot with Colors and Legend
counts <- table(mydata$Chromosome, mydata$Orthogroup)
col_sums <- colSums(counts)
counts <- rbind(counts, col_sums)

ordered_counts <- counts[,order(-counts[nrow(counts),])]

pdf("analysis/FoL_ls_genes/chr_LS_orthogroup_abundance.pdf")
plot <- barplot( ordered_counts[1:4,1:50], main="Orthogorup distribution of LS genes",
  xlab="Orthogroup", col=c("darkblue","red", "green", "yellow"),
 	legend = rownames(ordered_counts[1:4,]))
# text(plot, labels = rownames(ordered_counts), srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=.9)
# dev.copy(pdf,"analysis/FoL_ls_genes/chr_LS_orthogroup_abundance.pdf")
dev.off()
```



```bash
mydata <- read.table("gene_pred/annotations/F.oxysporum_fsp_cepae/Fus2_edited_v2/Fus2_edited_v2_gene_annotations.tab", header=TRUE, sep='\t')
```


Laura's RNAseq analysis identified differentially expresse genes between pathogenic and non-pathogenic fusarium isolates.

These files were downloaded ontot the cluster and parsed.

Note an additional column had been added to the csv file in microsoft excel using the command: =IF(F2<0.05, 1, 0)

```bash
InCsv=analysis/expression/warwick/F.oxysporum_fsp_cepae/Fus2/04_16/Fus2_path_vs_non_path.csv
OutTab=analysis/expression/warwick/F.oxysporum_fsp_cepae/Fus2/04_16/Fus2_path_vs_non_path.tab
cat $InCsv | sed 's/^M/\n/g' | sed 's/,/\t/g' > $OutTab
# when pasting the above command you will have to edit the ^M on the command line# - it has to be manually entered using ctrl+V +M

OrthoMCL_id=Fus2
OrthoMCL_id_list="125 A23 Fus2 55 A1_2 CB3 HB6 A13 A28 D2 PG fo47 4287"
OrthoMCL_path_ids="125 A23 Fus2"
OrthoMCL_nonpath_ids="A13 A28 D2 PG fo47"
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/orthology
$ProgDir/diff_expressed_orthogroups.py \
--FoC_orthogroup $OrthogroupsTxt \
--OrthoMCL_id $OrthoMCL_id \
--OrthoMCL_all $OrthoMCL_id_list \
--OrthoMCL_path $OrthoMCL_path_ids \
--OrthoMCL_nonpath $OrthoMCL_nonpath_ids \
--RNAseq_tab $OutTab \
> analysis/expression/warwick/F.oxysporum_fsp_cepae/Fus2/04_16/Fus2_path_vs_non_path_orthogroups.tab
```
