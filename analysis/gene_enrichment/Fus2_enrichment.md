

```bash
OutDir=analysis/enrichment/F.oxysporum_fsp_cepae/Fus2_canu_new
mkdir -p $OutDir
/home/armita/git_repos/emr_repos/scripts/fusarium/analysis/gene_enrichment/GO_prep_table.py --interpro gene_pred/interproscan/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_interproscan.tsv > $OutDir/Fus2_gene__GO_annots.tsv

AnnotTable=gene_pred/annotations/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_gene_annotations.tab
AllGenes=$OutDir/Fus2_all_genes.txt
cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
LSgenes=$OutDir/Fus2_LS_genes.txt
cat $AnnotTable | tail -n+2 | grep -e 'contig_10' -e 'contig16' -e 'contig_19' -e 'contig_21' | cut -f1 > $LSgenes
```


```bash
Var_set="contig_14_pilon contig_20_pilon contig_22_pilon"
PS_set="contig_10_pilon contig_16_pilon contig_19_pilon contig_21_pilon"
Core_set="contig_1_pilon contig_2_pilon contig_3_pilon contig_4_pilon contig_5_pilon contig_6_pilon contig_7_pilon contig_8_pilon contig_9_pilon contig_11_pilon contig_12_pilon contig_13_pilon contig_15_pilon contig_18_pilon"
/home/armita/git_repos/emr_repos/scripts/fusarium/analysis/gene_enrichment/IPR_prep_tables.py --interpro $AnnotTable --set1_name LS_contigs --set2_name core_contigs --contig_set1 $PS_set $Var_set --contig_set2 $Core_set --outdir $OutDir

# Run fishertest.r commands
cat results.txt | cut -f1,3 | grep "0\.00" | less

```


```bash
cp gene_pred/interproscan/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_interproscan.tsv $OutDir/.
cp gene_pred/annotations/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_gene_annotations.tab $OutDir/.
```
