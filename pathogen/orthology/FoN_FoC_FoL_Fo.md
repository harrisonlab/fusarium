
# Methodology 4


```bash
  ProjDir=/home/groups/harrisonlab/project_files/fusarium
  cd $ProjDir
  IsolateAbrv=FoN_vs_FoC_vs_FoL_vs_Fo
  WorkDir=analysis/orthology/orthomcl/$IsolateAbrv
  mkdir -p $WorkDir
  mkdir -p $WorkDir/formatted
  mkdir -p $WorkDir/goodProteins
  mkdir -p $WorkDir/badProteins  
```

## 4.1 Format fasta files

### 3.1.a Pathogenic isolates:

### for FoN N139
```bash
  Taxon_code=FoN
  Fasta_file=$(ls gene_pred/final_genes/F.oxysporum_fsp_narcissi/N139_ncbi/*/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for FoC Fus2
```bash
  Taxon_code=FoC
  Fasta_file=$(ls gene_pred/final_genes/F.oxysporum_fsp_cepae/Fus2_canu_new/*/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for FoL 4287

```bash
  Taxon_code=FoL
  Fasta_file=$(ls assembly/external_group/F.oxysporum_fsp_lycopersici/4287_v2/fungidb/FungiDB-29_Foxysporum4287_AnnotatedProteins.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### 3.1.b Non-pathogenic isolates:

### for Fo Fo47
```bash
  Taxon_code=fo47
  Fasta_file=$(ls assembly/external_group/F.oxysporum/fo47/broad/fusarium_oxysporum_fo47_1_proteins.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```




## 3.2 Filter proteins into good and poor sets.

```bash
  Input_dir=$WorkDir/formatted
  Min_length=10
  Max_percent_stops=20
  Good_proteins_file=$WorkDir/goodProteins/goodProteins.fasta
  Poor_proteins_file=$WorkDir/badProteins/poorProteins.fasta
  orthomclFilterFasta $Input_dir $Min_length $Max_percent_stops $Good_proteins_file $Poor_proteins_file
```

## 3.3.a Perform an all-vs-all blast of the proteins

```bash
  BlastDB=$WorkDir/blastall/$IsolateAbrv.db

  makeblastdb -in $Good_proteins_file -dbtype prot -out $BlastDB
  BlastOut=$WorkDir/all-vs-all_results.tsv
  mkdir -p $WorkDir/splitfiles

  SplitDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
  $SplitDir/splitfile_500.py --inp_fasta $Good_proteins_file --out_dir $WorkDir/splitfiles --out_base goodProteins

  ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/orthology  
  for File in $(find $WorkDir/splitfiles); do
    Jobs=$(qstat | grep 'blast_500' | grep 'qw' | wc -l)
    while [ $Jobs -gt 1 ]; do
      sleep 3
      printf "."
      Jobs=$(qstat | grep 'blast_500' | grep 'qw' | wc -l)
    done
    printf "\n"
    echo $File
    BlastOut=$(echo $File | sed 's/.fa/.tab/g')
    qsub $ProgDir/blast_500.sh $BlastDB $File $BlastOut
  done
```

## 3.3.b Merge the all-vs-all blast results  
```bash  
  MergeHits="$IsolateAbrv"_blast.tab
  printf "" > $MergeHits
  for Num in $(ls $WorkDir/splitfiles/*.tab | rev | cut -f1 -d '_' | rev | sort -n); do
    File=$(ls $WorkDir/splitfiles/*_$Num)
    cat $File
  done > $MergeHits
```

## 3.4 Perform ortholog identification

```bash
  ProgDir=~/git_repos/emr_repos/tools/pathogen/orthology/orthoMCL
  MergeHits="$IsolateAbrv"_blast.tab
  GoodProts=$WorkDir/goodProteins/goodProteins.fasta
  qsub $ProgDir/qsub_orthomcl.sh $MergeHits $GoodProts 5
```
<!--
## 3.5.a Manual identification of numbers of orthologous and unique genes


```bash
  for num in 1; do
    echo "The number of ortholog groups found in pathogen but absent in non-pathogens is:"
    cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep -v -e 'A28|' -e 'PG|' -e 'A13|' -e 'fo47|' -e 'CB3|' | grep 'Fus2|' | grep '125|' | grep 'A23|' |  wc -l
    echo "The number of ortholog groups unique to pathogens are:"
    cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep -v -e 'A28|' -e 'PG|' -e 'A13|' -e 'fo47|' -e 'CB3|' | grep 'Fus2|' | grep '125|' | grep 'A23|' | wc -l
    echo "The number of ortholog groups unique to non-pathogens are:"
    cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep -v -e 'Fus2|' -e '125|' -e 'A23|' | grep 'A13|' | grep 'A28|' | grep 'PG|' | grep 'fo47|' | grep 'CB3|' | wc -l
    echo "The number of ortholog groups common to all F. oxysporum isolates are:"
    cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep 'Fus2|' | grep '125|' | grep 'A23|' | grep 'A28|' | grep 'PG|' | grep 'A13|' | grep 'fo47|' | grep 'CB3|' | grep '4287|' |wc -l
  done
```

```
  The number of ortholog groups found in pathogen but absent in non-pathogens is:
  277
  The number of ortholog groups unique to pathogens are:
  277
  The number of ortholog groups unique to non-pathogens are:
  64
  The number of ortholog groups common to all F. oxysporum isolates are:
  10396
```

The number of ortholog groups shared between FoN and FoL was identified:

```bash
  echo "The number of ortholog groups common to FoN and FoL are:"
  cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep 'N139|' | grep '4287|' | wc -l
  cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep -v -e 'A28|' -e 'PG|' -e 'A13|' -e 'fo47|' -e 'CB3' | grep 'Fus2|' | grep '125|' | grep 'A23|' | grep '4287|' |  wc -l
```

```
  The number of ortholog groups common to FoC and FoL are:
  10629
  37
``` -->

## 3.5.b Plot venn diagrams:

```bash
  ProgDir=~/git_repos/emr_repos/tools/pathogen/orthology/venn_diagrams
  $ProgDir/venn_diag_4_way.r --inp $WorkDir/"$IsolateAbrv"_orthogroups.tab --out $WorkDir/"$IsolateAbrv"_orthogroups.pdf
```

Output was a pdf file of the venn diagram.

The following additional information was also provided. The format of the
following lines is as follows:

Isolate name (total number of orthogroups)
number of unique singleton genes
number of unique groups of inparalogs


```
[1] "FoC (13583)"
[1] 602
[1] 38
[1] "FoL (14751)"
[1] 1693
[1] 392
[1] "FoN (14152)"
[1] 1081
[1] 92
[1] "fo47 (14128)"
[1] 1056
[1] 57
```


#### 6.3) Extracting fasta files for all orthogroups

```bash
# IsolateAbrv=FoC_vs_Fo_vs_FoL_publication
WorkDir=analysis/orthology/orthomcl/$IsolateAbrv
ProgDir=~/git_repos/emr_repos/tools/pathogen/orthology/orthoMCL
GoodProt=$WorkDir/goodProteins/goodProteins.fasta
OutDir=$WorkDir/orthogroups_fasta
mkdir -p $OutDir
$ProgDir/orthoMCLgroups2fasta.py --orthogroups $WorkDir/"$IsolateAbrv"_orthogroups.txt --fasta $GoodProt --out_dir $OutDir > $OutDir/extractionlog.txt
# for File in $(ls -v $OutDir/orthogroup*.fa); do
# cat $File | grep '>' | tr -d '> '
# done > $OutDir/orthogroup_genes.txt
# cat $GoodProt | grep '>' | tr -d '> ' | grep -v -f $OutDir/orthogroup_genes.txt > $OutDir/singleton_genes.txt
# ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
# $ProgDir/extract_from_fasta.py --fasta $GoodProt --headers $OutDir/singleton_genes.txt > $OutDir/singleton_genes.fa
# echo "The numbe of singleton genes extracted is:"
# cat $OutDir/singleton_genes.fa | grep '>' | wc -l

```
