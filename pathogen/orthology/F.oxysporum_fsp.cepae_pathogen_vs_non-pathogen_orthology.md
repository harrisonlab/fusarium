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

# Downstream analysis

Particular orthogroups were analysed for expansion in isolates.

This section details the commands used and the results observed.

### Pathogenic Fusarium unique gene families

F. oxysporum Fus2 secreted proteins were parsed to the same format as the gene
names used in the analysis:

```bash
  SigP_Fus2=analysis/sigP_rxlr/F.oxysporum_fsp_cepae/Fus2/Fus2_sigP_RxLR.fa
  SigPDir=analysis/orthology/orthomcl/FoC_path_vs_non_path/SigP
  Orthogroups=analysis/orthology/orthomcl/FoC_path_vs_non_path/FoC_path_vs_non_path_orthogroups.txt
  SigP_ID_Fus2=$SigPDir/Fus2_aug_SigP_IDs.txt
  mkdir -p $SigPDir
  cat $SigP_Fus2 | grep '>' | cut -f1 | sed 's/>//g' | sed 's/ //g' | sed 's/g/Path|Fus2_g/g' > $SigP_ID_Fus2
```

Ortholog groups containing RxLR proteins were identified using the following
commands:
```bash
  SigP_Orthogroup_Path=$SigPDir/Path_SigP_Orthogroups.txt
  cat $Orthogroups | grep -w -f $SigP_ID_Fus2 > $SigP_Orthogroup_Path
  RxLR_Orthogroup_hits_10300=$RxLR_Dir/Pcac_RxLR_Orthogroups_hits.txt
  cat $Orthogroups | grep -o -w -f $SigP_ID_Fus2 > $RxLR_Orthogroup_hits_10300
```

All 61 predicted RxLRs were found to have orthologs in other taxa. The 61 RxLRs
were distributed through 41 orthogroups.

Orthogroup 8 was a large gene family of 298 genes. This included 46 P. cactorum
gene including 11 P. cactorum RxLRs.
```bash
  cat $Orthogroups | grep -w 'orthogroup8' | sed 's/ /\n/g' | sort | wc -l
  Orthogroup8_Pcac_ID=$RxLR_Dir/Pcac_RxLR_Orthogroups8_IDs.txt
  cat $Orthogroups | grep -w 'orthogroup8' | sed 's/ /\n/g' | sort | grep 'Pcac' > $Orthogroup8_Pcac_ID.txt
  cat $Orthogroups | grep -w 'orthogroup8' | sed 's/ /\n/g' | sort | cut -f1 -d'|' | uniq -c
```
```
  1 orthogroup8:
  46 Pcac
  100 Pinf
  73 Pram
  74 Psoj
```
