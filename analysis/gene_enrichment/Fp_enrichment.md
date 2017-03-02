# Enrichment of annotation terms of Fp genes

## GO enrichment using topGO


GO enrichment of terms in Fp vs FoC:

```bash
  OutDir=analysis/enrichment/F.proliferatum/A8_ncbi/Fp_vs_FoC
  mkdir -p $OutDir
  Append_InterProTSV=$OutDir/Fp_FoC_interproscan.tsv
  Fp_InterProTSV=$(ls gene_pred/interproscan/F.proliferatum/A8_ncbi/A8_ncbi_interproscan.tsv)
  cat $Fp_InterProTSV | sed -e 's/^/Fp_/g' > $Append_InterProTSV
  FoC_InterProTSV=$(ls gene_pred/interproscan/F.proliferatum/A8_ncbi/A8_ncbi_interproscan.tsv)
  cat $FoC_InterProTSV | sed -e 's/^/FoC_/g' >> $Append_InterProTSV

  ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/analysis/gene_enrichment
  $ProgDir/GO_prep_table.py --interpro $Append_InterProTSV > $OutDir/Fp_FoC_GO_annots.tsv

  Set1Genes=$OutDir/Fp_genes.txt
  Set2Genes=$OutDir/FoC_genes.txt
  cat $OutDir/Fp_FoC_GO_annots.tsv | cut -f1 | grep 'Fp' | sed -e 's/$/\t0.001/g'> $Set1Genes
  cat $OutDir/Fp_FoC_GO_annots.tsv | cut -f1 | grep 'FoC' | sed -e 's/$/\t1.00/g' > $Set2Genes
  AllGenes=$OutDir/Fp_FoC_all_genes.txt
  cat $Set1Genes $Set2Genes > $AllGenes

  $ProgDir/GO_enrichment.r --all_genes $AllGenes --GO_annotations $OutDir/Fp_FoC_GO_annots.tsv --out_dir $OutDir > $OutDir/output.txt
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

## Interproscan term enrichment using Fishers exact test in R

```bash
  OutDir=analysis/enrichment/F.proliferatum/A8_ncbi/Fp_vs_FoC
  mkdir -p $OutDir/tables

  Fp_AnnotTable=gene_pred/annotations/F.proliferatum/A8_ncbi/A8_ncbi_gene_annotations.tab
  FoC_AnnotTable=gene_pred/annotations/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_gene_annotations.tab
  Appended_AnnotTable=$OutDir/Fp_FoC_gene_annotations.tab
  cat $Fp_AnnotTable | head -n1 > $Appended_AnnotTable
  cat $Fp_AnnotTable | tail -n+2 | sed -e "s/^/Fp_/g" | sed 's/contig/Fp_contig/g' >> $Appended_AnnotTable
  cat $FoC_AnnotTable | tail -n+2 | sed -e "s/^/FoC_/g" | sed 's/contig/FoC_contig/g' >> $Appended_AnnotTable

  Set1Contigs=$(cat $Appended_AnnotTable | tail -n+2 | cut -f2 | grep 'Fp_' | sort | uniq | tr -d '\n' | sed 's/Fp/ Fp/g')
  Set2Contigs=$(cat $Appended_AnnotTable | tail -n+2 | cut -f2 | grep 'FoC_' | sort | uniq | tr -d '\n' | sed 's/FoC/ FoC/g')

  # Fp vs FoC

  ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/analysis/gene_enrichment
  $ProgDir/IPR_prep_tables.py --interpro $Appended_AnnotTable --set1_name Fp_genes --set2_name FoC_genes --contig_set1 $Set1Contigs  --contig_set2 $Set2Contigs --outdir $OutDir/tables
  $ProgDir/fisherstest.r --in_dir $OutDir/tables  --out_file $OutDir/results.tsv
  # Run fishertest.r commands
  cat $OutDir/results.tsv | cut -f1,2 | grep -e "e-" -e "0\.00" -e "0\.01" -e "0\.02" -e "0\.03" -e "0\.04"  > $OutDir/significant_terms.txt
  SigIPR=$(cat $OutDir/significant_terms.txt | cut -f1)
  for IPR in $SigIPR; do
    # echo $IPR
    cat $OutDir/tables/"$IPR"_fischertable.txt | grep "$IPR"
  done
```

Taken from: http://www.cs.tau.ac.il/~rshamir/ge/09/scribe/lec14a.pdf
'Since the signi􏰃cance test is performed for many groups, a multiple testing correction must be carried out in order to limit false positives. Both the Bonferroni and FDR methods are too stringent since there exist strong dependencies between groups (since they are often members of the same hierarchy). To get around these limitations, TANGO instead calculates the empirical p value distribution. For a given cluster Tj, TANGO samples many random sets of the same size & computes their p-values vs. each of the annotation sets Ai. Next, it also permutes gene IDs to eliminate dependency between annotation sets and target sets. This correction also applies for testing multiple clusters.'
