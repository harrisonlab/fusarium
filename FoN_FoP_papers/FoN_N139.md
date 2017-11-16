# FoN N139 publication commands


Most commands in assembly and annotation of the N139 genome are found in the
main README.md document for this repository (along with commands used in the
FoC paper). Those commands documented here are those that deviate from that
analysis.


## Removal of GPI anchor containing proteins from secreted proteins

Proteins containing GPI anchors were also removed using GPIsom


These proteins were identified through submitting the combined protein file to
the webserver at: http://gpi.unibe.ch


An output directory was made to download the file to:
"GPI anchored (C&N-term signal) (SignalP):"

```bash
  for Proteome in $(ls gene_pred/final_genes/*/*/*/final_genes_combined.pep.fasta | grep 'N139_ncbi'); do
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=gene_pred/trans_mem/$Organism/$Strain/GPIsom
    mkdir -p $OutDir
  done
```

The link to the proteins results was followed embedded in the text of
"GPI anchored (C&N-term signal) (SignalP):" thes +GPI anchor proteins
were copied and pasted into:

<!-- The link to the results was followed embedded in the text of
"Seqs with C-terminal signal (GPI-SOM):" thes +GPI anchor proteins
were copied and pasted into: -->

```bash
  nano gene_pred/trans_mem/$Organism/$Strain/GPIsom/GPI_pos.fa
```

Those proteins with GPI anchors were removed from lists of Signal
peptide containing proteins

```bash
for File in $(ls gene_pred/trans_mem/*/*/GPIsom/GPI_pos.fa | grep '10300'); do
Strain=$(echo $File | rev | cut -f3 -d '/' | rev)
Organism=$(echo $File | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
TmHeaders=$(echo "$File" | sed 's/.fa/.txt/g')
cat $File | grep '>' | cut -f1 -d ' ' | sed 's/>//g' > $TmHeaders
SigP=$(ls gene_pred/combined_sigP/$Organism/$Strain/*_sp_no_trans_mem.aa)
SigPHeaders=gene_pred/combined_sigP/$Organism/$Strain/"$Strain"_sp_no_trans_mem_headers.txt
cat $SigP | grep '>' | cut -f1 | sed 's/>//g'> $SigPHeaders
GoodHeaders=$(echo "$File" | sed 's/_pos.fa/_neg.txt/g')
cat $SigPHeaders | grep -v -f $TmHeaders > $GoodHeaders
OutDir=$(dirname $SigP)
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
# cat $SigP | grep -v -A1 -f $TmHeaders > $OutDir/"$Strain"_final_sp_no_trans_mem_no_GPI.aa
$ProgDir/extract_from_fasta.py --fasta $SigP --headers $GoodHeaders  > $OutDir/"$Strain"_final_sp_no_trans_mem_no_GPI.aa
echo "Number of SigP proteins:"
cat $SigP | grep '>' |wc -l
echo "Number with GPI anchors in entire proteome:"
cat $TmHeaders | wc -l
echo "Number without GPI anchors:"
cat $OutDir/"$Strain"_final_sp_no_trans_mem_no_GPI.aa | grep '>' | wc -l
echo "Number of gene models:"
cat $OutDir/"$Strain"_final_sp_no_trans_mem_no_GPI.aa | grep '>' | cut -f1 -d't' | sort | uniq |wc -l
done
```

```
  P.cactorum - 10300
  Number of SigP proteins:
  2244
  Number with GPI anchors in entire proteome:
  483
  Number without GPI anchors:
  1984
  Number of gene models:
  1973
```



## D) Secondary metabolites (Antismash and SMURF)

Antismash was run to identify clusters of secondary metabolite genes within
the genome. Antismash was run using the weserver at:
http://antismash.secondarymetabolites.org


Results of web-annotation of gene clusters within the assembly were downloaded to
the following directories:

```bash
  for Assembly in $(ls repeat_masked/*/*/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep -v 'HB17' | grep -e 'cepae' -e 'proliferatum' -e 'narcissi' | grep -e 'canu_new' -e 'ncbi' | grep 'N139'); do
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    OutDir=analysis/secondary_metabolites/antismash/$Organism/$Strain
    mkdir -p $OutDir
  done
```

```bash
  for Zip in $(ls analysis/secondary_metabolites/antismash/*/*/*.zip); do
    OutDir=$(dirname $Zip)
    unzip -d $OutDir $Zip
  done
```

```bash
  for AntiSmash in $(ls analysis/secondary_metabolites/antismash/*/*/*/*.final.gbk); do
    Organism=$(echo $AntiSmash | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $AntiSmash | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=analysis/secondary_metabolites/antismash/$Organism/$Strain
    Prefix=$OutDir/WT_antismash
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/secondary_metabolites
    $ProgDir/antismash2gff.py --inp_antismash $AntiSmash --out_prefix $Prefix

    # Identify secondary metabolites within predicted clusters
    printf "Number of secondary metabolite detected:\t"
    cat "$Prefix"_secmet_clusters.gff | wc -l
    GeneGff=gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed.gff3
    bedtools intersect -u -a $GeneGff -b "$Prefix"_secmet_clusters.gff > "$Prefix"_secmet_genes.gff
    cat "$Prefix"_secmet_genes.gff | grep -w 'mRNA' | cut -f9 | cut -f2 -d '=' | cut -f1 -d ';' > "$Prefix"_antismash_secmet_genes.txt
    bedtools intersect -wo -a $GeneGff -b "$Prefix"_secmet_clusters.gff | grep 'mRNA' | cut -f9,10,12,18 | sed "s/ID=//g" | perl -p -i -e "s/;Parent=g\w+//g" | perl -p -i -e "s/;Notes=.*//g" > "$Prefix"_secmet_genes.tsv
    printf "Number of predicted proteins in secondary metabolite clusters:\t"
    cat "$Prefix"_secmet_genes.txt | wc -l
    printf "Number of predicted genes in secondary metabolite clusters:\t"
    cat "$Prefix"_secmet_genes.gff | grep -w 'gene' | wc -l

      # Identify cluster finder additional non-secondary metabolite clusters
      printf "Number of cluster finder non-SecMet clusters detected:\t"
      cat "$Prefix"_clusterfinder_clusters.gff | wc -l
      GeneGff=gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed.gff3
      bedtools intersect -u -a $GeneGff -b "$Prefix"_clusterfinder_clusters.gff > "$Prefix"_clusterfinder_genes.gff
      cat "$Prefix"_clusterfinder_genes.gff | grep -w 'mRNA' | cut -f9 | cut -f2 -d '=' | cut -f1 -d ';' > "$Prefix"_clusterfinder_genes.txt

      printf "Number of predicted proteins in cluster finder non-SecMet clusters:\t"
      cat "$Prefix"_clusterfinder_genes.txt | wc -l
      printf "Number of predicted genes in cluster finder non-SecMet clusters:\t"
      cat "$Prefix"_clusterfinder_genes.gff | grep -w 'gene' | wc -l
  done
```

These clusters represented the following genes. Note that these numbers just
show the number of intersected genes with gff clusters and are not confirmed by
function

```
```

SMURF was also run to identify secondary metabolite gene clusters.

Genes needed to be parsed into a specific tsv format prior to submission on the
SMURF webserver.

```bash
for Gff in $(ls gene_pred/final_genes/F.oxysporum_fsp_narcissi/N139_ncbi/final/final_genes_appended.gff3); do
	Organism=$(echo $Gff | rev | cut -f4 -d '/' | rev)
	Strain=$(echo $Gff | rev | cut -f3 -d '/' | rev)
	echo "$Organism - $Strain"
  OutDir=analysis/secondary_metabolites/smurf/$Organism/$Strain
  mkdir -p $OutDir
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/secondary_metabolites
  $ProgDir/gff2smurf.py --gff $Gff > $OutDir/"$Strain"_genes_smurf.tsv
done
```

SMURF output was received by email and downloaded to the cluster in the output
directory above.

Output files were parsed into gff format:

```bash
  OutDir=analysis/secondary_metabolites/smurf/F.venenatum/WT
  Prefix="WT"
  GeneGff=gene_pred/final/F.venenatum/WT/final/final_genes_appended_renamed.gff3
  SmurfClusters=$OutDir/Secondary-Metabolite-Clusters.txt
  SmurfBackbone=$OutDir/Backbone-genes.txt
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/secondary_metabolites
  $ProgDir/smurf2gff.py --smurf_clusters $SmurfClusters --smurf_backbone $SmurfBackbone > $OutDir/Smurf_clusters.gff
  bedtools intersect -wo -a $GeneGff -b $OutDir/Smurf_clusters.gff | grep 'mRNA' | cut -f9,10,12,18 | sed "s/ID=//g" | perl -p -i -e "s/;Parent=g\w+//g" | perl -p -i -e "s/;Notes=.*//g" > $OutDir/"$Prefix"_smurf_secmet_genes.tsv
```

Total number of secondary metabolite clusters:

```bash
for Assembly in $(ls repeat_masked/*/*/illumina_assembly_ncbi/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep -w 'WT'); do
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
OutDir=analysis/secondary_metabolites/antismash/$Organism/$Strain
mkdir -p $OutDir
GeneGff=$(ls gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed.gff3)
AntismashClusters=$(ls analysis/secondary_metabolites/antismash/$Organism/$Strain/*_secmet_clusters.gff)
SmurfClusters=$(ls analysis/secondary_metabolites/smurf/$Organism/$Strain/Smurf_clusters.gff)
echo "Total number of Antismash clusters"
cat $AntismashClusters | wc -l
echo "Total number of SMURF clusters"
cat $SmurfClusters | wc -l
echo "number of Antismash clusters intersecting Smurf clusters"
bedtools intersect -a $AntismashClusters -b $SmurfClusters | wc -l
echo "number of Antismash clusters not intersecting Smurf clusters"
bedtools intersect -v -a $AntismashClusters -b $SmurfClusters | wc -l
echo "number of smurf clusters intersecting antismash clusters"
bedtools intersect -a $SmurfClusters -b $AntismashClusters | wc -l
echo "number of smurf clusters not intersecting antismash clusters"
bedtools intersect -v -a $SmurfClusters -b $AntismashClusters | wc -l
done
```


# F) Genes with transcription factor annotations:


A list of PFAM domains, superfamily annotations used as part of the DBD database
and a further set of interproscan annotations listed by Shelest et al 2017 were made
http://www.transcriptionfactor.org/index.cgi?Domain+domain:all
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5415576/

```bash
  for Interpro in $(ls gene_pred/interproscan/*/*/*_interproscan.tsv | grep 'N139' | grep 'ncbi'); do
    Organism=$(echo $Interpro | rev | cut -f3 -d '/' | rev)
    Strain=$(echo $Interpro | rev | cut -f2 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=analysis/transcription_factors/$Organism/$Strain
    mkdir -p $OutDir
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/transcription_factors
    $ProgDir/interpro2TFs.py --InterPro $Interpro > $OutDir/"$Strain"_TF_domains.tsv
    echo "total number of transcription factors"
    cat $OutDir/"$Strain"_TF_domains.tsv | cut -f1 | sort | uniq > $OutDir/"$Strain"_TF_gene_headers.txt
    cat $OutDir/"$Strain"_TF_gene_headers.txt | wc -l
  done
```










## Annotation Table


linked contigs were identified:


```bash


```

## G) Summarising annotation in annotation table

```bash
for GeneGff in $(ls gene_pred/final_genes/F.oxysporum_fsp_narcissi/N139_ncbi/final/final_genes_appended.gff3); do
Strain=$(echo $GeneGff | rev | cut -f3 -d '/' | rev)
Organism=$(echo $GeneGff | rev | cut -f4 -d '/' | rev)
Assembly=$(ls repeat_masked/$Organism/$Strain/*/*_contigs_unmasked.fa | grep 'ncbi')
GeneConversions=$(ls genome_submission/$Organism/$Strain/gag/edited/genome_gene_conversions.tsv)
# Antismash=$(ls analysis/secondary_metabolites/antismash/F.venenatum/WT/WT_antismash_secmet_genes.tsv)
# Smurf=$(ls analysis/secondary_metabolites/smurf/F.venenatum/WT/WT_smurf_secmet_genes.tsv)
TFs=$(ls analysis/transcription_factors/$Organism/$Strain/"$Strain"_TF_domains.tsv)
SigP=$(ls gene_pred/final_genes_signalp-4.1/F.oxysporum_fsp_narcissi/N139_ncbi/N139_ncbi_final_sp.aa)
TM_out=$(ls gene_pred/trans_mem/F.oxysporum_fsp_narcissi/N139_ncbi/N139_ncbi_TM_genes_pos.txt)
GPI_out=$(ls gene_pred/trans_mem/$Organism/$Strain/GPIsom/GPI_pos.fa)
MIMP_list=$(ls analysis/mimps/$Organism/$Strain/"$Strain"_genes_in_2kb_mimp.txt)
EffP_list=$(ls analysis/effectorP/$Organism/$Strain/"$Organism"_"$Strain"_EffectorP_headers.txt)
CAZY_list=$(ls gene_pred/CAZY/$Organism/$Strain/"$Strain"_CAZY_headers.txt)
InterPro=$(ls gene_pred/interproscan/$Organism/$Strain/*_interproscan.tsv)
SwissProt=$(ls gene_pred/swissprot/$Organism/$Strain/swissprot_vJul2016_tophit_parsed.tbl)
Orthology=$(ls analysis/orthology/orthomcl/FoN_vs_FoC_vs_FoL_vs_Fo/FoN_vs_FoC_vs_FoL_vs_Fo_orthogroups.txt)
OrthoStrainID='FoN'
OrthoStrainAll='FoN FoL FoC fo47'
OutDir=gene_pred/annotation/$Organism/$Strain
mkdir -p $OutDir
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes
# $ProgDir/build_annot_Fv.py --genome $Assembly --genes_gff $GeneGff --Antismash $Antismash --Smurf $Smurf --vitamins $Vitamins --TFs $TFs --InterPro $InterPro --Swissprot $SwissProt --orthogroups_PH1 $PH1_orthology --orthogroups_GR1 $GR1_orthology > $OutDir/"$Strain"_annotation_ncbi.tsv
$ProgDir/FoN_build_gene_annot_table.py --genome $Assembly --genes_gff $GeneGff --genes_renamed $GeneConversions --TFs $TFs --SigP $SigP --TM_list $TM_out --GPI_list $GPI_out --MIMP_list $MIMP_list --EffP_list $EffP_list --CAZY_list $CAZY_list --InterPro $InterPro --Swissprot $SwissProt --orthogroups $Orthology --strain_id $OrthoStrainID --OrthoMCL_all $OrthoStrainAll > $OutDir/"$Strain"_annotation_ncbi.tsv
done
```

```bash
for AnnotTab in $(ls gene_pred/annotation/F.venenatum/WT/WT_annotation_ncbi.tsv); do
  Strain=$(echo $AnnotTab | rev | cut -f2 -d '/' | rev)
  Organism=$(echo $AnnotTab | rev | cut -f3 -d '/' | rev)
  OutDir=$(ls -d analysis/transcription_factors/$Organism/$Strain)
  echo "Total number of clusters:"
  cat $AnnotTab | grep 'SecMet_cluster_' | cut -f7 | sort | uniq | wc -l
  cat $AnnotTab | grep 'SecMet_cluster_' | cut -f1,7 | cut -f1 > $OutDir/WT_SecMet_headers.txt
  cat $OutDir/WT_SecMet_headers.txt | cut -f1 -d '.' | sort | uniq | wc -l
  cat $AnnotTab | grep -e 'AS_' -e 'SM_' | cut -f1,14 | grep -v -e "\s$" | cut -f1 > $OutDir/WT_TF_SecMet_headers.txt
  cat $OutDir/WT_TF_SecMet_headers.txt | wc -l
done
```
