# For a comparison between pathogenic isolates and non-pathogenic isolates of FoC


```bash
  ProjDir=/home/groups/harrisonlab/project_files/fusarium
  cd $ProjDir
  IsolateAbrv=FoC_path_vs_nonpath
  WorkDir=analysis/orthology/orthomcl/$IsolateAbrv
  mkdir -p $WorkDir
  mkdir -p $WorkDir/formatted
  mkdir -p $WorkDir/goodProteins
  mkdir -p $WorkDir/badProteins  
```

## Format fasta files

### for FoC pathogenic isolates 125, A23 and Fus2
```bash
  Taxon_code=Path
  Isolate1=gene_pred/augustus/F.oxysporum_fsp_cepae/125/125_augustus_preds.aa
  Isolate2=gene_pred/augustus/F.oxysporum_fsp_cepae/A23/A23_augustus_preds.aa
  Isolate3=gene_pred/augustus/F.oxysporum_fsp_cepae/Fus2/Fus2_augustus_preds.aa
  Fasta_file=$WorkDir/"$Taxon_code"_concat.aa
  printf "" > $Fasta_file
  cat $Isolate1 | sed 's/>/>125_/g' >> $Fasta_file
  cat $Isolate2 | sed 's/>/>A23_/g' >> $Fasta_file
  cat $Isolate3 | sed 's/>/>Fus2_/g' >> $Fasta_file
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for FoC non-pathogenic isolates A28, D2, PG
```bash
  Taxon_code=NonP
  Isolate1=gene_pred/augustus/F.oxysporum_fsp_cepae/A28/A28_augustus_preds.aa
  Isolate2=gene_pred/augustus/F.oxysporum_fsp_cepae/D2/D2_augustus_preds.aa
  Isolate3=gene_pred/augustus/F.oxysporum_fsp_cepae/PG/PG_augustus_preds.aa
  Fasta_file=$WorkDir/"$Taxon_code"_concat.aa
  printf "" > $Fasta_file
  cat $Isolate1 | sed 's/>/>A28_/g' >> $Fasta_file
  cat $Isolate2 | sed 's/>/>D2_/g' >> $Fasta_file
  cat $Isolate3 | sed 's/>/>PG_/g' >> $Fasta_file
  Id_field=1
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
    Jobs=$(qstat | grep 'blast_500' | wc -l)
    while [ $Jobs -gt 32 ]; do
      sleep 10
      printf "."
      Jobs=$(qstat | grep 'blast_500' | wc -l)
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
  qsub $ProgDir/qsub_orthomcl.sh $MergeHits $GoodProts
```

## Plot venn diagrams:

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
  [1] "Aarb (9054)"
  [1] 134
  [1] 166
  [1] "Agai (8848)"
  [1] 270
  [1] 5
  [1] "AtAP (9794)"
  [1] 224
  [1] 89
  [1] "Aten (9680)"
  [1] 185
  [1] 36
```
