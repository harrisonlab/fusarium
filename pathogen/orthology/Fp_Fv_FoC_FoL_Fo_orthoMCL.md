
# Methodology 4


```bash
  ProjDir=/home/groups/harrisonlab/project_files/fusarium
  cd $ProjDir
  IsolateAbrv=Fp_Fv_FoC_FoL_Fo
  WorkDir=analysis/orthology/orthomcl/$IsolateAbrv
  mkdir -p $WorkDir
  mkdir -p $WorkDir/formatted
  mkdir -p $WorkDir/goodProteins
  mkdir -p $WorkDir/badProteins  
```

## 4.1 Format fasta files

### for Fp A8

```bash
  Taxon_code=Fp
  Fasta_file=$(ls gene_pred/final_genes/F.proliferatum/A8_ncbi/final/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for Fv 7600

```bash
  Taxon_code=Fv
  Fasta_file=$(ls assembly/external_group/F.verticilloides/7600/v3/Fusarium_verticillioides.ASM14955v1.pep_parsed.fa)
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

### for Fo Fo47
```bash
  Taxon_code=Fo
  Fasta_file=$(ls assembly/external_group/F.oxysporum/fo47/broad/fusarium_oxysporum_fo47_1_proteins.fasta)
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


## 4.2 Filter proteins into good and poor sets.

```bash
  Input_dir=$WorkDir/formatted
  Min_length=10
  Max_percent_stops=20
  Good_proteins_file=$WorkDir/goodProteins/goodProteins.fasta
  Poor_proteins_file=$WorkDir/badProteins/poorProteins.fasta
  orthomclFilterFasta $Input_dir $Min_length $Max_percent_stops $Good_proteins_file $Poor_proteins_file
```

## 4.3.a Perform an all-vs-all blast of the proteins

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

## 4.3.b Merge the all-vs-all blast results  
```bash  
  MergeHits=$WorkDir/"$IsolateAbrv"_blast.tab
  printf "" > $MergeHits
  for Num in $(ls $WorkDir/splitfiles/*.tab | rev | cut -f1 -d '_' | rev | sort -n); do
    File=$(ls $WorkDir/splitfiles/*_$Num)
    cat $File
  done > $MergeHits
```

## 4.4 Perform ortholog identification

```bash
  MergeHits=$WorkDir/"$IsolateAbrv"_blast.tab
  printf "" > $MergeHits
  for Num in $(ls $WorkDir/splitfiles/*.tab | rev | cut -f1 -d '_' | rev | sort -n); do
    File=$(ls $WorkDir/splitfiles/*_$Num)
    cat $File
  done > $MergeHits
  ProgDir=~/git_repos/emr_repos/tools/pathogen/orthology/orthoMCL
  MergeHits=$WorkDir/"$IsolateAbrv"_blast.tab
  GoodProts=$WorkDir/goodProteins/goodProteins.fasta
  qsub $ProgDir/qsub_orthomcl.sh $MergeHits $GoodProts 5
```
