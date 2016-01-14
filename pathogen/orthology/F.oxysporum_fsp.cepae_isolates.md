# For a comparison between all isolates of Fusarium


```bash
  ProjDir=/home/groups/harrisonlab/project_files/fusarium
  cd $ProjDir
  IsolateAbrv=FoC_125_55_A23_A28_D2_Fus2_PG_FoL_4287
  WorkDir=analysis/orthology/orthomcl/$IsolateAbrv
  mkdir -p $WorkDir
  mkdir -p $WorkDir/formatted
  mkdir -p $WorkDir/goodProteins
  mkdir -p $WorkDir/badProteins  
```

## Format fasta files

### for FoC 125
```bash
  Taxon_code=125
  Fasta_file=gene_pred/augustus/F.oxysporum_fsp_cepae/125/125_augustus_preds.aa
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for FoC 55
```bash
  Taxon_code=55
  Fasta_file=gene_pred/augustus/F.oxysporum_fsp_cepae/55/55_augustus_preds.aa
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for FoC A23
```bash
  Taxon_code=A23
  Fasta_file=gene_pred/augustus/F.oxysporum_fsp_cepae/A23/A23_augustus_preds.aa
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for FoC A28
```bash
  Taxon_code=A28
  Fasta_file=gene_pred/augustus/F.oxysporum_fsp_cepae/A28/A28_augustus_preds.aa
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```
### for FoC D2
```bash
  Taxon_code=D2
  Fasta_file=gene_pred/augustus/F.oxysporum_fsp_cepae/D2/D2_augustus_preds.aa
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```
### for FoC Fus2
```bash
  Taxon_code=Fus2
  Fasta_file=gene_pred/augustus/F.oxysporum_fsp_cepae/Fus2/Fus2_augustus_preds.aa
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```
### for FoC PG
```bash
  Taxon_code=PG
  Fasta_file=gene_pred/augustus/F.oxysporum_fsp_cepae/PG/PG_augustus_preds.aa
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for FoL 4287
```bash
  Taxon_code=4287
  Fasta_file=gene_pred/external_group/F.oxysporum_fsp_lycopersici/4287/fusarium_oxysporum_f._sp._lycopersici_4287_2_proteins.fasta
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


The normal method for venn diagram plotting could not be used due to the number
of isolates used. The following commands were used instead:

```bash
cat analysis/orthology/orthomcl/FoC_path_vs_non_path/FoC_path_vs_non_path_orthogroups.txt | grep -v -e 'A28_' -e 'D2_' -e 'PG_'| grep 'Fus2_' | grep '125_' | grep 'A23_' |  wc -l
cat analysis/orthology/orthomcl/FoC_path_vs_non_path/FoC_path_vs_non_path_orthogroups.txt | grep -v -e 'Fus2_' -e '125_' -e 'A23_' | grep 'A28_' | grep 'D2_' | grep 'PG_' |  wc -l
cat analysis/orthology/orthomcl/FoC_path_vs_non_path/FoC_path_vs_non_path_orthogroups.txt | grep 'Fus2_' | grep '125_' | grep 'A23_' | grep 'A28_' | grep 'D2_' | grep 'PG_' | wc -l
```

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

```

# Downstream analysis

Particular orthogroups were analysed for expansion in isolates.

This section details the commands used and the results observed.

### SIX9

FoL 4287 gene FOXG_14223 was identified as SIX9. This gene was observed to be a
member of orthogroup 240:

```bash
  Orthogroups=analysis/orthology/orthomcl/FoC_125_55_A23_A28_D2_Fus2_PG_FoL_4287/FoC_125_55_A23_A28_D2_Fus2_PG_FoL_4287_orthogroups.txt
  cat $Orthogroups | grep 'FOXG_14223'
```  
orthogroup 240 was noted to be a large gene family. The following number of
genes were associated with the following organisms:

```bash
  cat $Orthogroups | grep -w 'orthogroup240' | cut -f2 -d':' | sed 's/ /\n/g' | sort | less
  cat $Orthogroups | grep -w 'orthogroup240' | cut -f2 -d':' | sed 's/ /\n/g' | cut -f1 -d '|' | sort | uniq -c
```
output -
```
  8 125
  12 4287
  5 55
  6 A23
  2 A28
  2 D2
  5 Fus2
  3 PG
```
Pathogens
  8 125
  6 A23
  5 Fus2
non-pathogenic
  2 A28
  2 D2
  3 PG
Intermediate
  5 55
Fol
  12 4287 # two of these represent alternatively spliced transcripts (10)


### SIX7

Fus2 gene FUS2_g10785 was identified as SIX7. This gene was observed to be a
member of orthogroup 6993:

```bash
  Orthogroups=analysis/orthology/orthomcl/FoC_125_55_A23_A28_D2_Fus2_PG_FoL_4287/FoC_125_55_A23_A28_D2_Fus2_PG_FoL_4287_orthogroups.txt
  cat $Orthogroups | grep 'Fus2|g10785'
```  
orthogroup 6993 was noted to be expanded in Fus2, containing two genes.
The following number of genes were associated with the following organisms:

```bash
  cat $Orthogroups | grep -w 'orthogroup6993' | cut -f2 -d':' | sed 's/ /\n/g' | sort | less
  cat $Orthogroups | grep -w 'orthogroup6993' | cut -f2 -d':' | sed 's/ /\n/g' | cut -f1 -d '|' | sort | uniq -c
```

### RxLR_463

Fus2 gene FUS2_g10608 was identified as RxLR_463. This gene was observed to be a
member of orthogroup 453:

```bash
  Orthogroups=analysis/orthology/orthomcl/FoC_125_55_A23_A28_D2_Fus2_PG_FoL_4287/FoC_125_55_A23_A28_D2_Fus2_PG_FoL_4287_orthogroups.txt
  cat $Orthogroups | grep 'Fus2|g10608'
```  
orthogroup 453 was noted to be a large gene family. The following number of
genes were associated with the following organisms:

```bash
  cat $Orthogroups | grep -w 'orthogroup453' | cut -f2 -d':' | sed 's/ /\n/g' | sort | less
  cat $Orthogroups | grep -w 'orthogroup453' | cut -f2 -d':' | sed 's/ /\n/g' | cut -f1 -d '|' | sort | uniq -c
```

Pathogens
  4 125
  4 A23
  4 Fus2
non-pathogenic
  4 A28
  2 D2
  2 PG
Intermediate
  4 55
Fol
  3 4287



```
