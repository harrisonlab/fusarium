

```bash
OutDir=analysis/enrichment/F.oxysporum_fsp_cepae/Fus2_canu_new
mkdir -p $OutDir
InterProTSV=gene_pred/interproscan/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_interproscan.tsv
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/analysis/gene_enrichment
$ProgDir/GO_prep_table.py --interpro $InterProTSV > $OutDir/Fus2_gene_GO_annots.tsv

AnnotTable=gene_pred/annotations/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_gene_annotations.tab
AllGenes=$OutDir/Fus2_all_genes.txt
cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
LSgenes=$OutDir/Fus2_LS_genes.txt
Coregenes=$OutDir/Fus2_core_genes.txt
AllGenes=$OutDir/Fus2_all_genes.txt
cat $AnnotTable | tail -n+2 | grep -e 'contig_10' -e 'contig16' -e 'contig_19' -e 'contig_21' | cut -f1 | sed -e 's/$/\t0.001/g'> $LSgenes
cat $AnnotTable | tail -n+2 | grep -v -e 'contig_10' -e 'contig16' -e 'contig_19' -e 'contig_21' | cut -f1 | sed -e 's/$/\t1.00/g' > $Coregenes
cat $LSgenes $Coregenes > $AllGenes

$ProgDir/GO_enrichment.r --all_genes $AllGenes --GO_annotations $OutDir/Fus2_gene_GO_annots.tsv --out_dir $OutDir

```


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
cat $OutDir/results.tsv | cut -f1,3 | grep "0\.00" > $OutDir/significant_terms.txt
SigIPR="IPR012337 IPR008906 IPR004875 IPR006600 IPR007889 IPR009057 IPR013103 IPR003840 IPR025476 IPR001584 IPR031052 IPR010285 IPR005135 IPR000477 IPR025724 IPR018289 IPR002156 IPR001878 IPR029526 IPR002492 IPR005804 IPR004330 IPR022698 IPR007021 IPR004102 IPR013762 IPR022210 IPR031872 IPR003656 IPR012317 IPR012171 IPR021711 IPR006400 IPR001199 IPR002403 IPR000845 IPR006913 IPR001303 IPR008854 IPR008893 IPR013218 IPR018333 IPR032696 IPR032697 IPR011057 IPR008266 IPR002830 IPR009784 IPR012999 IPR022257 IPR000262 IPR011583 IPR000070 IPR001283 IPR008930 IPR001506 IPR001522 IPR006026 IPR006322 IPR009160 IPR015876 IPR022272 IPR023074 IPR030155 IPR002182 IPR001223 IPR013752 IPR000581 IPR001360 IPR014044 IPR020558 IPR000953 IPR023780"
for IPR in $SigIPR; do
  # echo $IPR
  cat $OutDir/tables/"$IPR"_fischertable.txt | grep "$IPR"
done

# PS vs LS
OutDir=analysis/enrichment/F.oxysporum_fsp_cepae/Fus2_canu_new/PS_vs_LS
mkdir -p $OutDir/tables
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/analysis/gene_enrichment
$ProgDir/IPR_prep_tables.py --interpro $AnnotTable --set1_name PS_contigs --set2_name LS_contigs --contig_set1 $PS_set --contig_set2 $Var_set --outdir $OutDir/tables
$ProgDir/fisherstest.r --in_dir $OutDir/tables  --out_file $OutDir/results.tsv
# Run fishertest.r commands
cat $OutDir/results.tsv | cut -f1,3 | grep "0\.00" > $OutDir/PS_enriched_terms.txt
cat $OutDir/results.tsv | cut -f1,4 | grep "0\.00" > $OutDir/LS_enriched_terms.txt
SigIPR="IPR025476 IPR003840 IPR010285 IPR013103 IPR012337 IPR027417 IPR001584 IPR007889 IPR025724 IPR009057 IPR008906"
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
cat $OutDir/results.tsv | cut -f1,3 | grep "0\.00" > $OutDir/significant_terms.txt
SigIPR="IPR025476 IPR003840 IPR010285 IPR013103 IPR012337 IPR027417 IPR001584 IPR007889 IPR025724 IPR009057 IPR008906"
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
