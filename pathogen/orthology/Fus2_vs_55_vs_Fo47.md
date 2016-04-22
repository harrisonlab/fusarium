

# Methodology 1


```bash
  ProjDir=/home/groups/harrisonlab/project_files/fusarium
  cd $ProjDir
  IsolateAbrv=Fus2_vs_55_vs_fo47
  WorkDir=analysis/orthology/orthomcl/$IsolateAbrv
  mkdir -p $WorkDir
  mkdir -p $WorkDir/formatted
  mkdir -p $WorkDir/goodProteins
  mkdir -p $WorkDir/badProteins  
```

## Format fasta files



### for FoC 55
```bash
  Taxon_code=55
  Fasta_file=$(ls gene_pred/codingquary/F.oxysporum_fsp_cepae/55/final/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for FoC Fus2
```bash
  Taxon_code=Fus2
  Fasta_file=$(ls gene_pred/codingquary/F.oxysporum_fsp_cepae/Fus2_edited_v2/final/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```
### for Fo Fo47
```bash
Taxon_code=fo47
Fasta_file=$(ls assembly/external_group/F.oxysporum/fo47/broad/fusarium_oxysporum_fo47_1_proteins.fasta)
Id_field=1
# Fasta_file=$(ls gene_pred/codingquary/F.oxysporum/fo47/final/final_genes_combined.pep.fasta)
# Id_field=1
orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```


## Filter proteins into good and poor sets.

```bash
  Input_dir=$WorkDir/formatted
  Min_length=10
  Max_percent_stops=20
  Good_proteins_file=$WorkDir/goodProteins/goodProteins.fasta
  Poor_proteins_file=$WorkDir/badProteins/poorProteins.fasta
  orthomclFilterFasta $Input_dir $Min_length $Max_percent_stops $Good_proteins_file $Poor_proteins_file
```

## Perform an all-vs-all blast of the proteins

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

## Merge the all-vs-all blast results  
```bash  
  MergeHits="$IsolateAbrv"_blast.tab
  printf "" > $MergeHits
  for Num in $(ls $WorkDir/splitfiles/*.tab | rev | cut -f1 -d '_' | rev | sort -n); do
    File=$(ls $WorkDir/splitfiles/*_$Num)
    cat $File
  done > $MergeHits
```

## Perform ortholog identification

```bash
  ProgDir=~/git_repos/emr_repos/tools/pathogen/orthology/orthoMCL
  MergeHits="$IsolateAbrv"_blast.tab
  GoodProts=$WorkDir/goodProteins/goodProteins.fasta
  qsub $ProgDir/qsub_orthomcl.sh $MergeHits $GoodProts 5
```

## Manual identification of numbers of orthologous and unique genes


```bash
  echo "The number of ortholog groups common to all isolates is:"
  cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep 'Fus2' | grep '55' | grep 'fo47' |  wc -l
```

```
  The number of ortholog groups unique to all isoaltes is:
  12285
```


## Plot venn diagrams:

```bash
ProgDir=~/git_repos/emr_repos/tools/pathogen/orthology/venn_diagrams
$ProgDir/venn_diag_3_way.r --inp $WorkDir/"$IsolateAbrv"_orthogroups.tab --out $WorkDir/"$IsolateAbrv"_orthogroups.pdf
```

Output was a pdf file of the venn diagram.

The following additional information was also provided. The format of the
following lines is as follows:

Isolate name (total number of orthogroups)
number of unique singleton genes
number of unique groups of inparalogs

[1] "55 (14072)"
[1] 14072
[1] "Fus2 (14070)"
[1] 14070
[1] "fo47 (14238)"
[1] 14238
[1] 13580


[1] "55 (14072)"
[1] 422
[1] 10
[1] "Fus2 (14070)"
[1] 447
[1] 11
[1] "fo47 (14238)"
[1] 1701
[1] 165
NULL
