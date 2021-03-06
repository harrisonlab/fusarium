

## Interproscan term enrichment using Fishers exact test in R

LS vs Core

```bash
Var_set="contig_14_pilon contig_20_pilon contig_22_pilon"
PS_set="contig_10_pilon contig_16_pilon contig_19_pilon contig_21_pilon"
Core_set="contig_1_pilon contig_2_pilon contig_3_pilon contig_4_pilon contig_5_pilon contig_6_pilon contig_7_pilon"
EffRich_Core="contig_8_pilon contig_9_pilon contig_13_pilon contig_15_pilon contig_18_pilon contig_11_pilon contig_12_pilon"
AnnotTable=gene_pred/annotations/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_gene_annotations.tab

# LS vs Core
OutDir=analysis/enrichment/F.oxysporum_fsp_cepae/Fus2_canu_new/NLR/LS_vs_core
rm -r $OutDir
mkdir -p $OutDir/tables
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/analysis/gene_enrichment
$ProgDir/NLR_enrichment_prep_tables.py --interpro $AnnotTable --set1_name LS_contigs --set2_name core_contigs --contig_set1 $PS_set $Var_set --contig_set2 $Core_set $EffRich_Core --outdir $OutDir/tables
$ProgDir/fisherstest.r --in_dir $OutDir/tables  --out_file $OutDir/results.tsv
# Run fishertest.r commands
cat $OutDir/results.tsv | cut -f1,3 | grep -e "e-" -e "0\.00" -e "0\.01" -e "0\.02" -e "0\.03" -e "0\.04"  > $OutDir/significant_terms.txt
SigIPR=$(cat $OutDir/significant_terms.txt | cut -f1)
for IPR in $(ls $OutDir/tables/*_fischertable.txt | rev | cut -f1 -d '/' | rev | sed 's/_fischertable.txt//g'); do
  # echo $IPR
  cat $OutDir/tables/"$IPR"_fischertable.txt | grep "$IPR"
done
```

PS vs Core

```bash
Var_set="contig_14_pilon contig_20_pilon contig_22_pilon"
PS_set="contig_10_pilon contig_16_pilon contig_19_pilon contig_21_pilon"
Core_set="contig_1_pilon contig_2_pilon contig_3_pilon contig_4_pilon contig_5_pilon contig_6_pilon contig_7_pilon"
EffRich_Core="contig_8_pilon contig_9_pilon contig_13_pilon contig_15_pilon contig_18_pilon contig_11_pilon contig_12_pilon"
AnnotTable=gene_pred/annotations/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_gene_annotations.tab

# PS vs Core
OutDir=analysis/enrichment/F.oxysporum_fsp_cepae/Fus2_canu_new/NLR/PS_vs_core
rm -r $OutDir
mkdir -p $OutDir/tables
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/analysis/gene_enrichment
$ProgDir/NLR_enrichment_prep_tables.py --interpro $AnnotTable --set1_name LS_contigs --set2_name core_contigs --contig_set1 $PS_set --contig_set2 $Core_set $EffRich_Core --outdir $OutDir/tables
$ProgDir/fisherstest.r --in_dir $OutDir/tables  --out_file $OutDir/results.tsv
# Run fishertest.r commands
cat $OutDir/results.tsv | cut -f1,3 | grep -e "e-" -e "0\.00" -e "0\.01" -e "0\.02" -e "0\.03" -e "0\.04"  > $OutDir/significant_terms.txt
SigIPR=$(cat $OutDir/significant_terms.txt | cut -f1)
for IPR in $(ls $OutDir/tables/*_fischertable.txt | rev | cut -f1 -d '/' | rev | sed 's/_fischertable.txt//g'); do
  # echo $IPR
  cat $OutDir/tables/"$IPR"_fischertable.txt | grep "$IPR"
done
```

non-PS vs Core

```bash
Var_set="contig_14_pilon contig_20_pilon contig_22_pilon"
PS_set="contig_10_pilon contig_16_pilon contig_19_pilon contig_21_pilon"
Core_set="contig_1_pilon contig_2_pilon contig_3_pilon contig_4_pilon contig_5_pilon contig_6_pilon contig_7_pilon"
EffRich_Core="contig_8_pilon contig_9_pilon contig_13_pilon contig_15_pilon contig_18_pilon contig_11_pilon contig_12_pilon"
AnnotTable=gene_pred/annotations/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_gene_annotations.tab

# PS vs Core
OutDir=analysis/enrichment/F.oxysporum_fsp_cepae/Fus2_canu_new/NLR/non-PS_vs_core
rm -r $OutDir
mkdir -p $OutDir/tables
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/analysis/gene_enrichment
$ProgDir/NLR_enrichment_prep_tables.py --interpro $AnnotTable --set1_name LS_contigs --set2_name core_contigs --contig_set1 $Var_set --contig_set2 $Core_set $EffRich_Core --outdir $OutDir/tables
$ProgDir/fisherstest.r --in_dir $OutDir/tables  --out_file $OutDir/results.tsv
# Run fishertest.r commands
cat $OutDir/results.tsv | cut -f1,3 | grep -e "e-" -e "0\.00" -e "0\.01" -e "0\.02" -e "0\.03" -e "0\.04"  > $OutDir/significant_terms.txt
SigIPR=$(cat $OutDir/significant_terms.txt | cut -f1)
for IPR in $(ls $OutDir/tables/*_fischertable.txt | rev | cut -f1 -d '/' | rev | sed 's/_fischertable.txt//g'); do
  # echo $IPR
  cat $OutDir/tables/"$IPR"_fischertable.txt | grep "$IPR"
done
```

83 P-loop NTPase proteins were identified in LS regions (SSF52540, IPR027417). Of these 22 were
identified as helitron helicasse transposbable elements through carrying associated domains (IPR025476, 'PF14214') and further 25 of these proteins carried 'helicase' domains (PF05970, IPR003840, IPR010285, PF14214, PS51194, IPR014001).


```bash
cat gene_pred/annotations/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_gene_annotations.tab | grep 'SSF52540' | grep -e 'contig_10' -e 'contig_14' -e 'contig_16' -e 'contig_19' -e 'contig_20' -e 'contig_21' -e 'contig_22' | wc -l

cat gene_pred/annotations/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_gene_annotations.tab | grep 'SSF52540' | grep -e 'contig_10' -e 'contig_14' -e 'contig_16' -e 'contig_19' -e 'contig_20' -e 'contig_21' -e 'contig_22' | grep -e 'IPR025476' -e 'PF14214' | wc -l

cat gene_pred/annotations/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_gene_annotations.tab \
| grep 'SSF52540' | grep -e 'contig_10' -e 'contig_14' -e 'contig_16' -e 'contig_19' -e 'contig_20' -e 'contig_21' -e 'contig_22' \
| grep -v  -e 'IPR025476' -e 'PF14214' \
| grep -e 'PF05970' -e 'IPR003840' -e 'IPR010285' -e 'PS51194' -e 'IPR014001'| wc -l

cat gene_pred/annotations/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_gene_annotations.tab \
| grep 'SSF52540' | grep -e 'contig_10' -e 'contig_14' -e 'contig_16' -e 'contig_19' -e 'contig_20' -e 'contig_21' -e 'contig_22' \
| grep -v  -e 'IPR025476' -e 'PF14214' \
| grep -v -e 'PF05970' -e 'IPR003840' -e 'IPR010285' -e 'PS51194' -e 'IPR014001'| wc -l

```
