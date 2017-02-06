# Enrichment of annotation terms of Fus2 genes

## GO enrichment using topGO

### Chromosomal comparisons

GO enrichment of terms in PS contigs:

```bash
OutDir=analysis/enrichment/F.oxysporum_fsp_cepae/Fus2_canu_new/PS_vs_core
mkdir -p $OutDir
InterProTSV=gene_pred/interproscan/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_interproscan.tsv
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/analysis/gene_enrichment
$ProgDir/GO_prep_table.py --interpro $InterProTSV > $OutDir/Fus2_gene_GO_annots.tsv

AnnotTable=gene_pred/annotations/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_gene_annotations.tab
AllGenes=$OutDir/Fus2_all_genes.txt
cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
Set1Genes=$OutDir/Fus2_Set1_genes.txt
Set2Genes=$OutDir/Fus2_Set2_genes.txt
AllGenes=$OutDir/Fus2_all_genes.txt
cat $AnnotTable | tail -n+2 | grep -e 'contig_10' -e 'contig_16' -e 'contig_19' -e 'contig_21' | cut -f1 | sed -e 's/$/\t0.001/g'> $Set1Genes
cat $AnnotTable | tail -n+2 | grep -v -e 'contig_10' -e 'contig_16' -e 'contig_19' -e 'contig_21' | cut -f1 | sed -e 's/$/\t1.00/g' > $Set2Genes
cat $Set1Genes $Set2Genes > $AllGenes

$ProgDir/GO_enrichment.r --all_genes $AllGenes --GO_annotations $OutDir/Fus2_gene_GO_annots.tsv --out_dir $OutDir > $OutDir/output.txt
```

GO enrichment of terms in other LS contigs:

```bash
OutDir=analysis/enrichment/F.oxysporum_fsp_cepae/Fus2_canu_new/other_LS_vs_core
mkdir -p $OutDir
InterProTSV=gene_pred/interproscan/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_interproscan.tsv
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/analysis/gene_enrichment
$ProgDir/GO_prep_table.py --interpro $InterProTSV > $OutDir/Fus2_gene_GO_annots.tsv

AnnotTable=gene_pred/annotations/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_gene_annotations.tab
AllGenes=$OutDir/Fus2_all_genes.txt
cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
Set1Genes=$OutDir/Fus2_Set1_genes.txt
Set2Genes=$OutDir/Fus2_Set2_genes.txt
AllGenes=$OutDir/Fus2_all_genes.txt
cat $AnnotTable | tail -n+2 | grep -e 'contig_14' -e 'contig_20' -e 'contig_22' | cut -f1 | sed -e 's/$/\t0.001/g'> $Set1Genes
cat $AnnotTable | tail -n+2 | grep -v -e 'contig_14' -e 'contig_20' -e 'contig_22' | cut -f1 | sed -e 's/$/\t1.00/g' > $Set2Genes
cat $Set1Genes $Set2Genes > $AllGenes

$ProgDir/GO_enrichment.r --all_genes $AllGenes --GO_annotations $OutDir/Fus2_gene_GO_annots.tsv --out_dir $OutDir > $OutDir/output.txt
```

GO enrichment of terms in LS contigs:

```bash
OutDir=analysis/enrichment/F.oxysporum_fsp_cepae/Fus2_canu_new/LS_vs_core
mkdir -p $OutDir
InterProTSV=gene_pred/interproscan/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_interproscan.tsv
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/analysis/gene_enrichment
$ProgDir/GO_prep_table.py --interpro $InterProTSV > $OutDir/Fus2_gene_GO_annots.tsv

AnnotTable=gene_pred/annotations/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_gene_annotations.tab
AllGenes=$OutDir/Fus2_all_genes.txt
cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
Set1Genes=$OutDir/Fus2_Set1_genes.txt
Set2Genes=$OutDir/Fus2_Set2_genes.txt
AllGenes=$OutDir/Fus2_all_genes.txt
cat $AnnotTable | tail -n+2 | grep -e 'contig_10' -e 'contig_16' -e 'contig_19' -e 'contig_21' -e 'contig_14' -e 'contig_20' -e 'contig_22' | cut -f1 | sed -e 's/$/\t0.001/g'> $Set1Genes
cat $AnnotTable | tail -n+2 | grep -v -e 'contig_10' -e 'contig_16' -e 'contig_19' -e 'contig_21' -e 'contig_14' -e 'contig_20' -e 'contig_22' | cut -f1 | sed -e 's/$/\t1.00/g' > $Set2Genes
cat $Set1Genes $Set2Genes > $AllGenes

$ProgDir/GO_enrichment.r --all_genes $AllGenes --GO_annotations $OutDir/Fus2_gene_GO_annots.tsv --out_dir $OutDir > $OutDir/output.txt
```

GO enrichment of terms in PS vs other LS contigs:

```bash
OutDir=analysis/enrichment/F.oxysporum_fsp_cepae/Fus2_canu_new/PS_vs_LS
mkdir -p $OutDir
InterProTSV=gene_pred/interproscan/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_interproscan.tsv
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/analysis/gene_enrichment
$ProgDir/GO_prep_table.py --interpro $InterProTSV > $OutDir/Fus2_gene_GO_annots.tsv

AnnotTable=gene_pred/annotations/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_gene_annotations.tab
AllGenes=$OutDir/Fus2_all_genes.txt
cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
Set1Genes=$OutDir/Fus2_Set1_genes.txt
Set2Genes=$OutDir/Fus2_Set2_genes.txt
AllGenes=$OutDir/Fus2_all_genes.txt
cat $AnnotTable | tail -n+2 | grep -e 'contig_10' -e 'contig_16' -e 'contig_19' -e 'contig_21' | cut -f1 | sed -e 's/$/\t0.001/g'> $Set1Genes
cat $AnnotTable | tail -n+2 | grep -e 'contig_14' -e 'contig_20' -e 'contig_22' | cut -f1 | sed -e 's/$/\t1.00/g' > $Set2Genes
cat $Set1Genes $Set2Genes > $AllGenes

$ProgDir/GO_enrichment.r --all_genes $AllGenes --GO_annotations $OutDir/Fus2_gene_GO_annots.tsv --out_dir $OutDir > $OutDir/output.txt
```

GO enrichment of terms in LS vs other PS contigs:

```bash
OutDir=analysis/enrichment/F.oxysporum_fsp_cepae/Fus2_canu_new/LS_vs_PS
mkdir -p $OutDir
InterProTSV=gene_pred/interproscan/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_interproscan.tsv
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/analysis/gene_enrichment
$ProgDir/GO_prep_table.py --interpro $InterProTSV > $OutDir/Fus2_gene_GO_annots.tsv

AnnotTable=gene_pred/annotations/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_gene_annotations.tab
AllGenes=$OutDir/Fus2_all_genes.txt
cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
Set1Genes=$OutDir/Fus2_Set1_genes.txt
Set2Genes=$OutDir/Fus2_Set2_genes.txt
AllGenes=$OutDir/Fus2_all_genes.txt
cat $AnnotTable | tail -n+2 | grep -e 'contig_14' -e 'contig_20' -e 'contig_22' | cut -f1 | sed -e 's/$/\t0.001/g'> $Set1Genes
cat $AnnotTable | tail -n+2 | grep -e 'contig_10' -e 'contig_16' -e 'contig_19' -e 'contig_21' | cut -f1 | sed -e 's/$/\t1.00/g' > $Set2Genes
cat $Set1Genes $Set2Genes > $AllGenes

$ProgDir/GO_enrichment.r --all_genes $AllGenes --GO_annotations $OutDir/Fus2_gene_GO_annots.tsv --out_dir $OutDir > $OutDir/output.txt
```

### Duplicated gene comparisons


GO enrichment of terms in duplicated vs single genes:

```bash
OutDir=analysis/enrichment/F.oxysporum_fsp_cepae/Fus2_canu_new/duplicated_vs_single
mkdir -p $OutDir
InterProTSV=gene_pred/interproscan/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_interproscan.tsv
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/analysis/gene_enrichment
$ProgDir/GO_prep_table.py --interpro $InterProTSV > $OutDir/Fus2_gene_GO_annots.tsv

DupGenesTab=/home/sobczm/popgen/codon/blast/dagchainer/testing/no_transposon/non-transposon_duplications_summaryf_ann
DupGenes=$OutDir/duplicated_genes.txt
cat $DupGenesTab | tail -n+2 | cut -f1,3,5 | sed 's/\t/\n/g' | sed 's/,/\n/g' | grep -v "^$" | grep 'Fus2' | sed 's/Fus2_//g' | sort | uniq > $DupGenes

# AnnotTable=gene_pred/annotations/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_gene_annotations.tab
# AllGenes=$OutDir/Fus2_all_genes.txt
# UniqGenes=$OutDir/Fus2_unique_genes.txt
# cat $AnnotTable | tail -n+2  | cut -f1 | grep -v -f $DupGenes > $UniqGenes

Set1Genes=$OutDir/Fus2_Set1_genes.txt
Set2Genes=$OutDir/Fus2_Set2_genes.txt
AllGenes=$OutDir/Fus2_all_genes.txt
cat $AnnotTable | tail -n+2 | grep -f $DupGenes | cut -f1 | sed -e 's/$/\t0.001/g'> $Set1Genes
cat $AnnotTable | tail -n+2 | grep -v -f $DupGenes | cut -f1 | sed -e 's/$/\t1.00/g' > $Set2Genes
cat $Set1Genes $Set2Genes > $AllGenes

$ProgDir/GO_enrichment.r --all_genes $AllGenes --GO_annotations $OutDir/Fus2_gene_GO_annots.tsv --out_dir $OutDir > $OutDir/output.txt
```


GO enrichment of duplicated genes from LS vs core contigs:

```bash
OutDir=analysis/enrichment/F.oxysporum_fsp_cepae/Fus2_canu_new/duplicated_LS_vs_core
mkdir -p $OutDir
InterProTSV=gene_pred/interproscan/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_interproscan.tsv
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/analysis/gene_enrichment
$ProgDir/GO_prep_table.py --interpro $InterProTSV > $OutDir/Fus2_gene_GO_annots.tsv

DupGenesTab=/home/sobczm/popgen/codon/blast/dagchainer/testing/no_transposon/non-transposon_duplications_summaryf_ann
DupGenes=$OutDir/duplicated_genes.txt
cat $DupGenesTab | tail -n+2 | cut -f1,3,5 | sed 's/\t/\n/g' | sed 's/,/\n/g' | grep -v "^$" | grep 'Fus2' | sed 's/Fus2_//g' | sort | uniq > $DupGenes

# AnnotTable=gene_pred/annotations/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_gene_annotations.tab
# AllGenes=$OutDir/Fus2_all_genes.txt
# UniqGenes=$OutDir/Fus2_unique_genes.txt
# cat $AnnotTable | tail -n+2  | cut -f1 | grep -v -f $DupGenes > $UniqGenes

Set1Genes=$OutDir/Fus2_Set1_genes.txt
Set2Genes=$OutDir/Fus2_Set2_genes.txt
AllGenes=$OutDir/Fus2_all_genes.txt
cat $AnnotTable | tail -n+2 | grep -e 'contig_10' -e 'contig_16' -e 'contig_19' -e 'contig_21' -e 'contig_14' -e 'contig_20' -e 'contig_22' | grep -f $DupGenes | cut -f1 | sed -e 's/$/\t0.001/g'> $Set1Genes
cat $AnnotTable | grep -v -e 'contig_10' -e 'contig_16' -e 'contig_19' -e 'contig_21' -e 'contig_14' -e 'contig_20' -e 'contig_22' | tail -n+2 | grep -f $DupGenes | cut -f1 | sed -e 's/$/\t1.00/g' > $Set2Genes
cat $Set1Genes $Set2Genes > $AllGenes

$ProgDir/GO_enrichment.r --all_genes $AllGenes --GO_annotations $OutDir/Fus2_gene_GO_annots.tsv --out_dir $OutDir > $OutDir/output.txt
```

GO enrichment of duplicated genes from PS vs LS contigs:

```bash
OutDir=analysis/enrichment/F.oxysporum_fsp_cepae/Fus2_canu_new/duplicated_PS_vs_LS
mkdir -p $OutDir
InterProTSV=gene_pred/interproscan/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_interproscan.tsv
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/analysis/gene_enrichment
$ProgDir/GO_prep_table.py --interpro $InterProTSV > $OutDir/Fus2_gene_GO_annots.tsv

DupGenesTab=/home/sobczm/popgen/codon/blast/dagchainer/testing/no_transposon/non-transposon_duplications_summaryf_ann
DupGenes=$OutDir/duplicated_genes.txt
cat $DupGenesTab | tail -n+2 | cut -f1,3,5 | sed 's/\t/\n/g' | sed 's/,/\n/g' | grep -v "^$" | grep 'Fus2' | sed 's/Fus2_//g' | sort | uniq > $DupGenes

# AnnotTable=gene_pred/annotations/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_gene_annotations.tab
# AllGenes=$OutDir/Fus2_all_genes.txt
# UniqGenes=$OutDir/Fus2_unique_genes.txt
# cat $AnnotTable | tail -n+2  | cut -f1 | grep -v -f $DupGenes > $UniqGenes

Set1Genes=$OutDir/Fus2_Set1_genes.txt
Set2Genes=$OutDir/Fus2_Set2_genes.txt
AllGenes=$OutDir/Fus2_all_genes.txt
AnnotTable=gene_pred/annotations/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_gene_annotations.tab
cat $AnnotTable | tail -n+2 | grep -e 'contig_10' -e 'contig_16' -e 'contig_19' -e 'contig_21' | grep -f $DupGenes | cut -f1 | sed -e 's/$/\t0.001/g'> $Set1Genes
cat $AnnotTable | tail -n+2 | grep -e 'contig_14' -e 'contig_20' -e 'contig_22' | grep -f $DupGenes | cut -f1 | sed -e 's/$/\t1.00/g' > $Set2Genes
cat $Set1Genes $Set2Genes > $AllGenes

$ProgDir/GO_enrichment.r --all_genes $AllGenes --GO_annotations $OutDir/Fus2_gene_GO_annots.tsv --out_dir $OutDir > $OutDir/output.txt
```

GO enrichment of duplicated genes from LS vs PS contigs:

```bash
OutDir=analysis/enrichment/F.oxysporum_fsp_cepae/Fus2_canu_new/duplicated_LS_vs_PS
mkdir -p $OutDir
InterProTSV=gene_pred/interproscan/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_interproscan.tsv
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/analysis/gene_enrichment
$ProgDir/GO_prep_table.py --interpro $InterProTSV > $OutDir/Fus2_gene_GO_annots.tsv

DupGenesTab=/home/sobczm/popgen/codon/blast/dagchainer/testing/no_transposon/non-transposon_duplications_summaryf_ann
DupGenes=$OutDir/duplicated_genes.txt
cat $DupGenesTab | tail -n+2 | cut -f1,3,5 | sed 's/\t/\n/g' | sed 's/,/\n/g' | grep -v "^$" | grep 'Fus2' | sed 's/Fus2_//g' | sort | uniq > $DupGenes

# AnnotTable=gene_pred/annotations/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_gene_annotations.tab
# AllGenes=$OutDir/Fus2_all_genes.txt
# UniqGenes=$OutDir/Fus2_unique_genes.txt
# cat $AnnotTable | tail -n+2  | cut -f1 | grep -v -f $DupGenes > $UniqGenes

Set1Genes=$OutDir/Fus2_Set1_genes.txt
Set2Genes=$OutDir/Fus2_Set2_genes.txt
AllGenes=$OutDir/Fus2_all_genes.txt
AnnotTable=gene_pred/annotations/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_gene_annotations.tab
cat $AnnotTable | tail -n+2 | grep -e 'contig_14' -e 'contig_20' -e 'contig_22' | grep -f $DupGenes | cut -f1 | sed -e 's/$/\t0.001/g'> $Set1Genes
cat $AnnotTable | tail -n+2 | grep -e 'contig_10' -e 'contig_16' -e 'contig_19' -e 'contig_21' | grep -f $DupGenes | cut -f1 | sed -e 's/$/\t1.00/g' > $Set2Genes
cat $Set1Genes $Set2Genes > $AllGenes

$ProgDir/GO_enrichment.r --all_genes $AllGenes --GO_annotations $OutDir/Fus2_gene_GO_annots.tsv --out_dir $OutDir > $OutDir/output.txt
```

Enrichment of genes lacking GO terms:

PS_vs_LS
```bash
OutDir=analysis/enrichment/F.oxysporum_fsp_cepae/Fus2_canu_new/fishers_GO_PS_vs_LS
mkdir -p $OutDir
InterProTSV=gene_pred/interproscan/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_interproscan.tsv
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/analysis/gene_enrichment
$ProgDir/GO_prep_table.py --interpro $InterProTSV > $OutDir/Fus2_gene_GO_annots.tsv

Annotated_list=$OutDir/Fus2_gene_GO_annots.txt
cat $OutDir/Fus2_gene_GO_annots.tsv | cut -f1 > $Annotated_list

Var_set='-e "contig_14_pilon" -e "contig_20_pilon" -e "contig_22_pilon"'
PS_set='-e "contig_10_pilon" -e "contig_16_pilon" -e "contig_19_pilon" -e "contig_21_pilon"'
Core_set='-e "contig_1_pilon" -e "contig_2_pilon" -e "contig_3_pilon" -e "contig_4_pilon" -e "contig_5_pilon" -e "contig_6_pilon"'

AnnotTable=gene_pred/annotations/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_gene_annotations.tab
Set1GO=$(cat $AnnotTable | tail -n+2 | grep -e "contig_14_pilon" -e "contig_20_pilon" -e "contig_22_pilon" | grep -f $Annotated_list | wc -l)
Set1NoGO=$(cat $AnnotTable | tail -n+2 | grep -e "contig_14_pilon" -e "contig_20_pilon" -e "contig_22_pilon" | cut -f1 | grep -v -f $Annotated_list | wc -l)
Set2GO=$(cat $AnnotTable | tail -n+2 | grep -e "contig_10_pilon" -e "contig_16_pilon" -e "contig_19_pilon" -e "contig_21_pilon" | cut -f1 | grep -f $Annotated_list | wc -l)
Set2NoGO=$(cat $AnnotTable | tail -n+2 | grep -e "contig_10_pilon" -e "contig_16_pilon" -e "contig_19_pilon" -e "contig_21_pilon" | cut -f1 | grep -v -f $Annotated_list | wc -l)
printf "Go_annotated\t$Set1GO\t$Set2GO\nUnannotated\t$Set1NoGO\t$Set2NoGO\n" > $OutDir/GO_PS_vs_LS_fischertable.txt
$ProgDir/fisherstest.r --in_dir $OutDir  --out_file $OutDir/results.tsv
```

LS_vs_Core
```bash
OutDir=analysis/enrichment/F.oxysporum_fsp_cepae/Fus2_canu_new/fishers_GO_LS_vs_core
mkdir -p $OutDir
InterProTSV=gene_pred/interproscan/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_interproscan.tsv
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/analysis/gene_enrichment
$ProgDir/GO_prep_table.py --interpro $InterProTSV > $OutDir/Fus2_gene_GO_annots.tsv

Annotated_list=$OutDir/Fus2_gene_GO_annots.txt
cat $OutDir/Fus2_gene_GO_annots.tsv | cut -f1 > $Annotated_list

Var_set='-e "contig_14_pilon" -e "contig_20_pilon" -e "contig_22_pilon"'
PS_set='-e "contig_10_pilon" -e "contig_16_pilon" -e "contig_19_pilon" -e "contig_21_pilon"'
Core_set='-e "contig_1_pilon" -e "contig_2_pilon" -e "contig_3_pilon" -e "contig_4_pilon" -e "contig_5_pilon" -e "contig_6_pilon"'

AnnotTable=gene_pred/annotations/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_gene_annotations.tab
Set1GO=$(cat $AnnotTable | tail -n+2 | grep -e "contig_10_pilon" -e "contig_16_pilon" -e "contig_19_pilon" -e "contig_21_pilon" -e "contig_14_pilon" -e "contig_20_pilon" -e "contig_22_pilon" | cut -f1 | grep -f $Annotated_list | wc -l)
Set1NoGO=$(cat $AnnotTable | tail -n+2 | grep  -e "contig_10_pilon" -e "contig_16_pilon" -e "contig_19_pilon" -e "contig_21_pilon" -e "contig_14_pilon" -e "contig_20_pilon" -e "contig_22_pilon" | cut -f1 |grep -v -f $Annotated_list | wc -l)
Set2GO=$(cat $AnnotTable | tail -n+2 | grep -e "contig_1_pilon" -e "contig_2_pilon" -e "contig_3_pilon" -e "contig_4_pilon" -e "contig_5_pilon" -e "contig_6_pilon" | cut -f1 | grep -f $Annotated_list | wc -l)
Set2NoGO=$(cat $AnnotTable | tail -n+2 | grep -e "contig_1_pilon" -e "contig_2_pilon" -e "contig_3_pilon" -e "contig_4_pilon" -e "contig_5_pilon" -e "contig_6_pilon" | cut -f1 | grep -v -f $Annotated_list | wc -l)
printf "Go_annotated\t$Set1GO\t$Set2GO\nUnannotated\t$Set1NoGO\t$Set2NoGO\n" > $OutDir/GO_LS_vs_core_fischertable.txt
$ProgDir/fisherstest.r --in_dir $OutDir  --out_file $OutDir/results.tsv
```

## Interproscan term enrichment using Fishers exact test in R

```bash
Var_set="contig_14_pilon contig_20_pilon contig_22_pilon"
PS_set="contig_10_pilon contig_16_pilon contig_19_pilon contig_21_pilon"
Core_set="contig_1_pilon contig_2_pilon contig_3_pilon contig_4_pilon contig_5_pilon contig_6_pilon contig_7_pilon"
EffRich_Core="contig_8_pilon contig_9_pilon contig_13_pilon contig_15_pilon contig_18_pilon contig_11_pilon contig_12_pilon"
AnnotTable=gene_pred/annotations/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_gene_annotations.tab

# LS vs Core
OutDir=analysis/enrichment/F.oxysporum_fsp_cepae/Fus2_canu_new/LS_vs_core
mkdir -p $OutDir/tables
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/analysis/gene_enrichment
$ProgDir/IPR_prep_tables.py --interpro $AnnotTable --set1_name LS_contigs --set2_name core_contigs --contig_set1 $PS_set $Var_set --contig_set2 $Core_set $EffRich_Core --outdir $OutDir/tables
$ProgDir/fisherstest.r --in_dir $OutDir/tables  --out_file $OutDir/results.tsv
# Run fishertest.r commands
cat $OutDir/results.tsv | cut -f1,3 | grep -e "e-" -e "0\.00" -e "0\.01" -e "0\.02" -e "0\.03" -e "0\.04"  > $OutDir/significant_terms.txt
SigIPR=$(cat $OutDir/significant_terms.txt | cut -f1)
# SigIPR="IPR012337 IPR008906 IPR004875 IPR006600 IPR007889 IPR009057 IPR013103 IPR003840 IPR025476 IPR001584 IPR031052 IPR010285 IPR005135 IPR000477 IPR025724 IPR018289 IPR002156 IPR001878 IPR029526 IPR002492 IPR005804 IPR004330 IPR022698 IPR007021 IPR004102 IPR013762 IPR022210 IPR031872 IPR003656 IPR012317 IPR012171 IPR021711 IPR006400 IPR001199 IPR002403 IPR000845 IPR006913 IPR001303 IPR008854 IPR008893 IPR013218 IPR018333 IPR032696 IPR032697 IPR011057 IPR008266 IPR002830 IPR009784 IPR012999 IPR022257 IPR000262 IPR011583 IPR000070 IPR001283 IPR008930 IPR001506 IPR001522 IPR006026 IPR006322 IPR009160 IPR015876 IPR022272 IPR023074 IPR030155 IPR002182 IPR001223 IPR013752 IPR000581 IPR001360 IPR014044 IPR020558 IPR000953 IPR023780"
for IPR in $SigIPR; do
  # echo $IPR
  cat $OutDir/tables/"$IPR"_fischertable.txt | grep "$IPR"
done

# PS vs LS
OutDir=analysis/enrichment/F.oxysporum_fsp_cepae/Fus2_canu_new/PS_vs_LS
mkdir -p $OutDir/tables
rm $OutDir/results.tsv
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/analysis/gene_enrichment
$ProgDir/IPR_prep_tables.py --interpro $AnnotTable --set1_name PS_contigs --set2_name LS_contigs --contig_set1 $PS_set --contig_set2 $Var_set --outdir $OutDir/tables
$ProgDir/fisherstest.r --in_dir $OutDir/tables  --out_file $OutDir/results.tsv
# Run fishertest.r commands
cat $OutDir/results.tsv | cut -f1,3 | grep -e "e-" -e "0\.00" -e "0\.01" -e "0\.02" -e "0\.03" -e "0\.04"  > $OutDir/PS_enriched_terms.txt
cat $OutDir/results.tsv | cut -f1,4 | grep -e "e-" -e "0\.00" -e "0\.01" -e "0\.02" -e "0\.03" -e "0\.04"  > $OutDir/LS_enriched_terms.txt
#SigIPR="IPR025476 IPR003840 IPR010285 IPR013103 IPR012337 IPR027417 IPR001584 IPR007889 IPR025724 IPR009057 IPR008906"
echo "PS enriched terms:"
SigIPR=$(cat $OutDir/PS_enriched_terms.txt| cut -f1)
for IPR in $SigIPR; do
  # echo $IPR
  cat $OutDir/tables/"$IPR"_fischertable.txt | grep "$IPR"
done
echo "LS enriched terms:"
SigIPR=$(cat $OutDir/LS_enriched_terms.txt| cut -f1)
for IPR in $SigIPR; do
  # echo $IPR
  cat $OutDir/tables/"$IPR"_fischertable.txt | grep "$IPR"
done

# # PS vs rest
# OutDir=analysis/enrichment/F.oxysporum_fsp_cepae/Fus2_canu_new/PS_vs_non-PS
# mkdir -p $OutDir/tables
# rm $OutDir/results.tsv
# ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/analysis/gene_enrichment
# $ProgDir/IPR_prep_tables.py --interpro $AnnotTable --set1_name PS_contigs --set2_name non-PS_contigs --contig_set1 $PS_set --contig_set2 $Var_set $Core_set --outdir $OutDir/tables
# $ProgDir/fisherstest.r --in_dir $OutDir/tables  --out_file $OutDir/results.tsv
# # Run fishertest.r commands
# cat $OutDir/results.tsv | cut -f1,3 | grep -e "e-" -e "0\.00" -e "0\.01" -e "0\.02" -e "0\.03" -e "0\.04"  > $OutDir/PS_enriched_terms.txt
# cat $OutDir/results.tsv | cut -f1,4 | grep -e "e-" -e "0\.00" -e "0\.01" -e "0\.02" -e "0\.03" -e "0\.04"  > $OutDir/non-PS_enriched_terms.txt
# #SigIPR="IPR025476 IPR003840 IPR010285 IPR013103 IPR012337 IPR027417 IPR001584 IPR007889 IPR025724 IPR009057 IPR008906"
# echo "PS enriched terms:"
# SigIPR=$(cat $OutDir/PS_enriched_terms.txt| cut -f1)
# for IPR in $SigIPR; do
#   # echo $IPR
#   cat $OutDir/tables/"$IPR"_fischertable.txt | grep "$IPR"
# done
# echo "non-PS enriched terms:"
# SigIPR=$(cat $OutDir/non-PS_enriched_terms.txt| cut -f1)
# for IPR in $SigIPR; do
#   # echo $IPR
#   cat $OutDir/tables/"$IPR"_fischertable.txt | grep "$IPR"
# done
#
# # other_LS vs rest
# OutDir=analysis/enrichment/F.oxysporum_fsp_cepae/Fus2_canu_new/other_LS_vs_rest
# mkdir -p $OutDir/tables
# rm $OutDir/results.tsv
# ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/analysis/gene_enrichment
# $ProgDir/IPR_prep_tables.py --interpro $AnnotTable --set1_name other_LS_contigs --set2_name rest_contigs --contig_set1 $Var_set --contig_set2 $PS_set $Core_set --outdir $OutDir/tables
# $ProgDir/fisherstest.r --in_dir $OutDir/tables  --out_file $OutDir/results.tsv
# # Run fishertest.r commands
# cat $OutDir/results.tsv | cut -f1,3 | grep -e "e-" -e "0\.00" -e "0\.01" -e "0\.02" -e "0\.03" -e "0\.04"  > $OutDir/other_LS_enriched_terms.txt
# cat $OutDir/results.tsv | cut -f1,4 | grep -e "e-" -e "0\.00" -e "0\.01" -e "0\.02" -e "0\.03" -e "0\.04"  > $OutDir/rest_enriched_terms.txt
# #SigIPR="IPR025476 IPR003840 IPR010285 IPR013103 IPR012337 IPR027417 IPR001584 IPR007889 IPR025724 IPR009057 IPR008906"
# echo "other_LS enriched terms:"
# SigIPR=$(cat $OutDir/other_LS_enriched_terms.txt| cut -f1)
# for IPR in $SigIPR; do
#   # echo $IPR
#   cat $OutDir/tables/"$IPR"_fischertable.txt | grep "$IPR"
# done
# echo "rest enriched terms:"
# SigIPR=$(cat $OutDir/rest_enriched_terms.txt| cut -f1)
# for IPR in $SigIPR; do
#   # echo $IPR
#   cat $OutDir/tables/"$IPR"_fischertable.txt | grep "$IPR"
# done

# PS vs core
OutDir=analysis/enrichment/F.oxysporum_fsp_cepae/Fus2_canu_new/PS_vs_core
mkdir -p $OutDir/tables
rm $OutDir/results.tsv
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/analysis/gene_enrichment
$ProgDir/IPR_prep_tables.py --interpro $AnnotTable --set1_name PS_contigs --set2_name core_contigs --contig_set1 $PS_set --contig_set2 $Core_set --outdir $OutDir/tables
$ProgDir/fisherstest.r --in_dir $OutDir/tables  --out_file $OutDir/results.tsv
# Run fishertest.r commands
cat $OutDir/results.tsv | cut -f1,3 | grep -e "e-" -e "0\.00" -e "0\.01" -e "0\.02" -e "0\.03" -e "0\.04"  > $OutDir/PS_enriched_terms.txt
cat $OutDir/results.tsv | cut -f1,4 | grep -e "e-" -e "0\.00" -e "0\.01" -e "0\.02" -e "0\.03" -e "0\.04"  > $OutDir/core_enriched_terms.txt
#SigIPR="IPR025476 IPR003840 IPR010285 IPR013103 IPR012337 IPR027417 IPR001584 IPR007889 IPR025724 IPR009057 IPR008906"
echo "PS enriched terms:"
SigIPR=$(cat $OutDir/PS_enriched_terms.txt| cut -f1)
for IPR in $SigIPR; do
  # echo $IPR
  cat $OutDir/tables/"$IPR"_fischertable.txt | grep "$IPR"
done
echo "core enriched terms:"
SigIPR=$(cat $OutDir/core_enriched_terms.txt| cut -f1)
for IPR in $SigIPR; do
  # echo $IPR
  cat $OutDir/tables/"$IPR"_fischertable.txt | grep "$IPR"
done

# other_LS vs core
OutDir=analysis/enrichment/F.oxysporum_fsp_cepae/Fus2_canu_new/other_LS_vs_core
mkdir -p $OutDir/tables
rm $OutDir/results.tsv
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/analysis/gene_enrichment
$ProgDir/IPR_prep_tables.py --interpro $AnnotTable --set1_name other_LS_contigs --set2_name core_contigs --contig_set1 $Var_set --contig_set2 $Core_set --outdir $OutDir/tables
$ProgDir/fisherstest.r --in_dir $OutDir/tables  --out_file $OutDir/results.tsv
# Run fishertest.r commands
cat $OutDir/results.tsv | cut -f1,3 | grep -e "e-" -e "0\.00" -e "0\.01" -e "0\.02" -e "0\.03" -e "0\.04"  > $OutDir/other_LS_enriched_terms.txt
cat $OutDir/results.tsv | cut -f1,4 | grep -e "e-" -e "0\.00" -e "0\.01" -e "0\.02" -e "0\.03" -e "0\.04"  > $OutDir/core_enriched_terms.txt
#SigIPR="IPR025476 IPR003840 IPR010285 IPR013103 IPR012337 IPR027417 IPR001584 IPR007889 IPR025724 IPR009057 IPR008906"
echo "other_LS enriched terms:"
SigIPR=$(cat $OutDir/other_LS_enriched_terms.txt| cut -f1)
for IPR in $SigIPR; do
  # echo $IPR
  cat $OutDir/tables/"$IPR"_fischertable.txt | grep "$IPR"
done
echo "rest enriched terms:"
SigIPR=$(cat $OutDir/core_enriched_terms.txt| cut -f1)
for IPR in $SigIPR; do
  # echo $IPR
  cat $OutDir/tables/"$IPR"_fischertable.txt | grep "$IPR"
done

# Effector rich core contigs (4813) vs other contigs (14265):
OutDir=analysis/enrichment/F.oxysporum_fsp_cepae/Fus2_canu_new/effector_rich_vs_all
mkdir -p $OutDir/tables
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/analysis/gene_enrichment
$ProgDir/IPR_prep_tables.py --interpro $AnnotTable --set1_name Effector_rich --set2_name core_contigs --contig_set1 $EffRich_Core --contig_set2 $Core_set $PS_set $Var_set --outdir $OutDir/tables
$ProgDir/fisherstest.r --in_dir $OutDir/tables  --out_file $OutDir/results.tsv
# Run fishertest.r commands
cat $OutDir/results.tsv | cut -f1,3 | grep -e "e-" -e "0\.00" -e "0\.01" -e "0\.02" -e "0\.03" -e "0\.04"  > $OutDir/significant_terms.txt
# SigIPR="IPR025476 IPR003840 IPR010285 IPR013103 IPR012337 IPR027417 IPR001584 IPR007889 IPR025724 IPR009057 IPR008906"
SigIPR=$(cat $OutDir/significant_terms.txt | cut -f1)
for IPR in $SigIPR; do
  # echo $IPR
  cat $OutDir/tables/"$IPR"_fischertable.txt | grep "$IPR"
done

```

Taken from: http://www.cs.tau.ac.il/~rshamir/ge/09/scribe/lec14a.pdf
'Since the signi􏰃cance test is performed for many groups, a multiple testing correction must be carried out in order to limit false positives. Both the Bonferroni and FDR methods are too stringent since there exist strong dependencies between groups (since they are often members of the same hierarchy). To get around these limitations, TANGO instead calculates the empirical p value distribution. For a given cluster Tj, TANGO samples many random sets of the same size & computes their p-values vs. each of the annotation sets Ai. Next, it also permutes gene IDs to eliminate dependency between annotation sets and target sets. This correction also applies for testing multiple clusters.'


```bash
cp gene_pred/interproscan/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_interproscan.tsv $OutDir/.
cp gene_pred/annotations/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_gene_annotations.tab $OutDir/.
```
