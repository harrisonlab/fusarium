These commands were used to identify highly expressed genes in the Fus2 assembly


# 1) QC

QC was performed as part of gene prediction. Commands are Documented in the Fusarium README.md

# 2) Align reads vs. Fus2 genome
Alignments of RNAseq reads were made against the Fus2 Genome using tophat as
part of gene prediction. Commands are Documented in the Fusarium README.md

# 3) Quantification of gene model expression


A list of transposons were made for Fus2:

```bash
mkdir -p analysis/transposons
cat gene_pred/interproscan/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_interproscan.tsv | grep -e 'IPR000477' -e 'IPR012337' -e 'IPR018289' -e 'PF03221' -e 'PF00078' -e 'IPR025476' -e 'IPR008906' -e 'transpos' -e 'integrase' | cut -f1 | sort | uniq > analysis/transposons/transposon_headers.txt
```

Expression was determined for aligned reads against predcited gene models using
cufflinks:

for Fus2:
```bash
  for Alignment in $(ls alignment/F.oxysporum_fsp_cepae/Fus2_canu_new/*/accepted_hits.bam | grep -e 'Fus2_72hrs_rep' -e 'control_72hrs_rep' -e 'FO47_72hrs_rep' | grep 'FO47'); do
    Gff=gene_pred/final_genes/F.oxysporum_fsp_cepae/Fus2_canu_new/final/final_genes_appended.gff3
    OutDir=$(dirname $Alignment)
    mkdir -p $OutDir/fpkm
    cufflinks -p 8 -o $OutDir/fpkm -G $Gff $Alignment
  done

fpkm_files=$(ls alignment/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_72hrs_rep*/fpkm/genes.fpkm_tracking | sed -r 's/\n/ /g')
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/RNAseq
Cazy=gene_pred/CAZY/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_CAZY_secreted_headers.txt
EffectorP=analysis/effectorP/F.oxysporum_fsp_cepae/Fus2_canu_new/F.oxysporum_fsp_cepae_Fus2_canu_
Cazy=gene_pred/CAZY/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_CAZY_headers.txt
EffectorP=analysis/effectorP/F.oxysporum_fsp_cepae/Fus2_canu_new/F.oxysporum_fsp_cepae_Fus2_canu_new_EffectorP_headers.txt
EffectorP=analysis/effectorP/F.oxysporum_fsp_cepae/Fus2_canu_new/F.oxysporum_fsp_cepae_Fus2_canu_new_EffectorP_secreted_headers.txt
Transposons=analysis/transposons/transposon_headers.txt
Metabolites=analysis/antismash/F.oxysporum_fsp_cepae/Fus2_canu_new/metabolite_cluster_gene_headers.txt
Annotations=gene_pred/annotations/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_gene_annotations.tab
$ProgDir/fpkm_by_gene.py --fpkm_file $fpkm_files --CAZY_headers $Cazy --effectorP_headers $EffectorP  --transposon_headers $Transposons --metabolite_headers $Metabolites --annotation_table $Annotations > analysis/expression/Fus2_expressed_genes.tsv



echo "Number of variaible core genes in top 50 expressed genes:"
cat analysis/expression/Fus2_expressed_genes.tsv | head -n 50 | grep -w -v -e 'contig_10_pilon' -e 'contig_14_pilon' -e 'contig_16_pilon' -e 'contig_19_pilon' -e 'contig_20_pilon' -e 'contig_21_pilon' -e 'contig_22_pilon' | grep -w -e 'contig_9_pilon' -e 'contig_11_pilon' -e 'contig_12_pilon' -e 'contig_13_pilon' -e 'contig_15_pilon' -e 'contig_17_pilon'  -e 'contig_18_pilon'| wc -l

# in Top 50 genes
echo "Number of effectorP genes in top 50 expressed genes:"
cat analysis/expression/Fus2_expressed_genes.tsv | head -n 50 | grep 'EffP' | wc -l
echo "Number of secreted effectorP genes in top 50 expressed genes:"
cat analysis/expression/Fus2_expressed_genes.tsv | head -n 50 | grep 'EffP' | cut -f17 | grep 'Yes' | wc -l
echo "Number of CAZY genes in top 50 expressed genes:"
cat analysis/expression/Fus2_expressed_genes.tsv | head -n 50 | grep -i 'cazy' | wc -l
echo "Number of CAZY genes in top 50 expressed genes:"
cat analysis/expression/Fus2_expressed_genes.tsv | head -n 50 | grep -i 'cazy' | cut -f17 | grep 'Yes' | wc -l
echo "Number of secreted genes in top 50 expressed genes:"
cat analysis/expression/Fus2_expressed_genes.tsv | head -n 50 | cut -f17 | grep 'Yes' | wc -l
echo "Number of ribosomal genes in top 50 expressed genes:"
cat analysis/expression/Fus2_expressed_genes.tsv | head -n 50 | grep -i 'ribosomal' | wc -l
echo "Number of transposon genes in top 50 expressed genes:"
cat analysis/expression/Fus2_expressed_genes.tsv | head -n 50 | grep -i 'Transposon' | wc -l

# in Top 50 LS genes
echo "Number of LS genes in top 50 expressed genes:"
cat analysis/expression/Fus2_expressed_genes.tsv | grep -w -e 'contig_10_pilon' -e 'contig_14_pilon' -e 'contig_16_pilon' -e 'contig_19_pilon' -e 'contig_20_pilon' -e 'contig_21_pilon' -e 'contig_22_pilon' > analysis/expression/Fus2_expressed_LS_genes.tsv

echo "Number of effectorP genes in top 50 expressed genes:"
cat analysis/expression/Fus2_expressed_LS_genes.tsv | head -n 50 | grep 'EffP' | wc -l
echo "Number of CAZY genes in top 50 expressed genes:"
cat analysis/expression/Fus2_expressed_LS_genes.tsv | head -n 50 | grep -i 'cazy' | wc -l
echo "Number of secreted genes in top 50 expressed genes:"
cat analysis/expression/Fus2_expressed_LS_genes.tsv | head -n 50 | cut -f17 | grep 'Yes' | wc -l
echo "Number of ribosomal genes in top 50 expressed genes:"
cat analysis/expression/Fus2_expressed_LS_genes.tsv | head -n 50 | grep -i 'ribosomal' | wc -l
echo "Number of transposon genes in top 50 expressed genes:"
cat analysis/expression/Fus2_expressed_LS_genes.tsv | head -n 50 | grep -i 'Transposon' | wc -l

# in Top 50 core genes
echo "Number of LS genes in top 50 expressed genes:"
cat analysis/expression/Fus2_expressed_genes.tsv | grep -v -w -e 'contig_10_pilon' -e 'contig_14_pilon' -e 'contig_16_pilon' -e 'contig_19_pilon' -e 'contig_20_pilon' -e 'contig_21_pilon' -e 'contig_22_pilon' > analysis/expression/Fus2_expressed_core_genes.tsv

echo "Number of effectorP genes in top 50 expressed genes:"
cat analysis/expression/Fus2_expressed_core_genes.tsv | head -n 50 | grep 'EffP' | wc -l
echo "Number of CAZY genes in top 50 expressed genes:"
cat analysis/expression/Fus2_expressed_core_genes.tsv | head -n 50 | grep -i 'cazy' | wc -l
echo "Number of secreted genes in top 50 expressed genes:"
cat analysis/expression/Fus2_expressed_core_genes.tsv | head -n 50 | cut -f17 | grep 'Yes' | wc -l
echo "Number of ribosomal genes in top 50 expressed genes:"
cat analysis/expression/Fus2_expressed_core_genes.tsv | head -n 50 | grep -i 'ribosomal' | wc -l
echo "Number of transposon genes in top 50 expressed genes:"
cat analysis/expression/Fus2_expressed_core_genes.tsv | head -n 50 | grep -i 'Transposon' | wc -l
```
The top 50 core genes covered a range between genes 1 and 71. As such the function of those additional 21 LS genes within this set were investigated.

```bash
echo "Number of effectorP genes in top 21 expressed genes:"
cat analysis/expression/Fus2_expressed_LS_genes.tsv | head -n 21 | grep 'EffP' | wc -l
echo "Number of CAZY genes in top 21 expressed genes:"
cat analysis/expression/Fus2_expressed_LS_genes.tsv | head -n 21 | grep -i 'cazy' | wc -l
echo "Number of secreted genes in top 21 expressed genes:"
cat analysis/expression/Fus2_expressed_LS_genes.tsv | head -n 21 | cut -f17 | grep 'Yes' | wc -l
echo "Number of ribosomal genes in top 21 expressed genes:"
cat analysis/expression/Fus2_expressed_LS_genes.tsv | head -n 21 | grep -i 'ribosomal' | wc -l
echo "Number of transposon genes in top 21 expressed genes:"
cat analysis/expression/Fus2_expressed_LS_genes.tsv | head -n 21 | grep -i 'Transposon' | wc -l
```

Checking control reads aligning to Fus2

```bash
fpkm_files=$(ls alignment/F.oxysporum_fsp_cepae/Fus2_canu_new/control_72hrs_rep*/fpkm/genes.fpkm_tracking | sed -r 's/\n/ /g')
# Each individual file was checked for contamination
# fpkm_files=$(ls alignment/F.oxysporum_fsp_cepae/Fus2_canu_new/control_72hrs_rep*/fpkm/genes.fpkm_tracking | sed -r 's/\n/ /g' | head -n3 | tail -n1)
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/RNAseq
Cazy=gene_pred/CAZY/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_CAZY_secreted_headers.txt
EffectorP=analysis/effectorP/F.oxysporum_fsp_cepae/Fus2_canu_new/F.oxysporum_fsp_cepae_Fus2_canu_new_EffectorP_secreted_headers.txt
Transposons=analysis/transposons/transposon_headers.txt
Metabolites=analysis/antismash/79c1471f-4a2b-41f7-ba36-18ba94675f59/metabolite_cluster_gene_headers.txt
Annotations=gene_pred/annotations/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_gene_annotations.tab
$ProgDir/fpkm_by_gene.py --fpkm_file $fpkm_files --CAZY_headers $Cazy --effectorP_headers $EffectorP  --transposon_headers $Transposons --metabolite_headers $Metabolites --annotation_table $Annotations > analysis/expression/Control_vs_Fus2_expressed_genes.tsv
```

Checking fo47 reads aligning to Fus2

```bash
fpkm_files=$(ls alignment/F.oxysporum_fsp_cepae/Fus2_canu_new/FO47_72hrs_rep*/fpkm/genes.fpkm_tracking | sed -r 's/\n/ /g')
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/RNAseq
Cazy=gene_pred/CAZY/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_CAZY_secreted_headers.txt
EffectorP=analysis/effectorP/F.oxysporum_fsp_cepae/Fus2_canu_new/F.oxysporum_fsp_cepae_Fus2_canu_new_EffectorP_secreted_headers.txt
Transposons=analysis/transposons/transposon_headers.txt
Metabolites=analysis/antismash/79c1471f-4a2b-41f7-ba36-18ba94675f59/metabolite_cluster_gene_headers.txt
Annotations=gene_pred/annotations/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_gene_annotations.tab
$ProgDir/fpkm_by_gene.py --fpkm_file $fpkm_files --CAZY_headers $Cazy --effectorP_headers $EffectorP  --transposon_headers $Transposons --metabolite_headers $Metabolites --annotation_table $Annotations > analysis/expression/fo47_vs_Fus2_expressed_genes.tsv

```



for Fo47:

```bash
for Alignment in $(ls alignment/F.oxysporum/fo47/*/accepted_hits.bam | grep -e 'FO47_72hrs_rep' -e 'control_72hrs_rep' | grep 'control'); do
Gff=assembly/external_group/F.oxysporum/fo47/broad/fusarium_oxysporum_fo47_1_transcripts.gtf
OutDir=$(dirname $Alignment)
mkdir -p $OutDir/fpkm
cufflinks -p 8 -o $OutDir/fpkm -G $Gff $Alignment
done

fpkm_files=$(ls alignment/F.oxysporum/fo47/FO47_72hrs_rep*/fpkm/genes.fpkm_tracking | sed -r 's/\n/ /g')
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/RNAseq
Cazy=gene_pred/CAZY/F.oxysporum/fo47/fo47_CAZY_secreted_headers.txt
EffectorP=analysis/effectorP/F.oxysporum/fo47/F.oxysporum_fo47_EffectorP_secreted_headers.txt
# Transposons=analysis/transposons/transposon_headers.txt
Metabolites=analysis/antismash/04337fe1-dce8-460b-9111-873d28f4f8e2/metabolite_cluster_gene_headers.txt
Annotations=gene_pred/annotations/F.oxysporum/fo47/fo47_gene_annotations.tab
$ProgDir/fpkm_by_gene.py --fpkm_file $fpkm_files --CAZY_headers $Cazy --effectorP_headers $EffectorP --annotation_table $Annotations --metabolite_headers $Metabolites > analysis/expression/fo47_expressed_genes.tsv


  echo "Number of effectorP genes in top 50 expressed genes:"
  cat analysis/expression/fo47_expressed_genes.tsv | head -n 50 | grep 'EffP' | wc -l
  echo "Number of CAZY genes in top 50 expressed genes:"
  cat analysis/expression/fo47_expressed_genes.tsv | head -n 50 | grep 'CAZY' | wc -l
  echo "Number of secreted genes in top 50 expressed genes:"
  cat analysis/expression/fo47_expressed_genes.tsv | head -n 50 | cut -f17 | grep 'Yes' | wc -l
  echo "Number of ribosomal genes in top 50 expressed genes:"
  cat analysis/expression/fo47_expressed_genes.tsv | head -n 50 | grep -i 'ribosomal' | wc -l
```


Annotaions for expressed secreted CAZY genes were identified:

```bash
cat analysis/expression/Fus2_expressed_genes.tsv | grep -w -f gene_pred/CAZY/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_CAZY_secreted_headers.txt > gene_pred/CAZY/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_CAZY_secreted_annot.tsv
```

<!-- ```bash
cat gene_pred/CAZY/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_CAZY_secreted_headers.txt | cut -f6 > tmp.txt
cat tmp.txt gene_pred/CAZY/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_CAZY_secreted_headers.txt | sort | uniq -u
cat gene_pred/CAZY/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_CAZY_secreted_headers.txt | grep 'g10836.t1'
``` -->

Annotaions for expressed secreted antismash genes were identified:
```bash
cat analysis/expression/Fus2_expressed_genes.tsv | grep -w -f analysis/antismash/F.oxysporum_fsp_cepae/Fus2_canu_new/metabolite_cluster_gene_headers.txt > analysis/antismash/F.oxysporum_fsp_cepae/Fus2_canu_new/metabolite_cluster_gene_annot.tsv
```

Annotaions for genes within 2kb of MIMPs were :
```bash
cat analysis/expression/Fus2_expressed_genes.tsv | grep -w -f analysis/mimps/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_genes_in_2kb_mimp.txt > analysis/mimps/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_genes_in_2kb_mimp_annot.tsv
```

Annotations of all other gene families:

```bash
cat analysis/effectorP/F.oxysporum_fsp_cepae/Fus2_canu_new/F.oxysporum_fsp_cepae_Fus2_canu_new_EffectorP_secreted_headers.txt  | cut -f1 > analysis/effectorP/F.oxysporum_fsp_cepae/Fus2_canu_new/F.oxysporum_fsp_cepae_Fus2_canu_new_EffectorP_secreted_headers_parsed.txt

cat analysis/expression/Fus2_expressed_genes.tsv \
| grep -v -w -f analysis/effectorP/F.oxysporum_fsp_cepae/Fus2_canu_new/F.oxysporum_fsp_cepae_Fus2_canu_new_EffectorP_secreted_headers_parsed.txt | grep -v -w -f gene_pred/CAZY/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_CAZY_secreted_headers.txt \
| grep -v -w -f analysis/antismash/F.oxysporum_fsp_cepae/Fus2_canu_new/metabolite_cluster_gene_headers.txt \
| grep -v -w -f analysis/mimps/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_genes_in_2kb_mimp.txt \
> analysis/expression/Fus2_expressed_genes_non-effector.tsv
```
