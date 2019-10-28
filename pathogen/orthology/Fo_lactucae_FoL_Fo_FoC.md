
# Methodology 4


```bash
  ProjDir=/home/groups/harrisonlab/project_files/fusarium
  cd $ProjDir
  IsolateAbrv=Fo_lactucae_FoL_Fo_FoC
  WorkDir=analysis/orthology/orthomcl/$IsolateAbrv
  mkdir -p $WorkDir
  mkdir -p $WorkDir/formatted
  mkdir -p $WorkDir/goodProteins
  mkdir -p $WorkDir/badProteins  
```

## 4.1 Format fasta files

### 3.1.a Pathogenic isolates:

### for Fo lactucae R1 AJ520
```bash
  Taxon_code=FoLaR1
  Fasta_file=$(ls gene_pred/final/F.oxysporum_fsp_lactucae/AJ520/publication/final_genes_appended_renamed_ncbi.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for Fo lactucae R4 AJ516
```bash
  Taxon_code=FoLaR4
  Fasta_file=$(ls gene_pred/final/F.oxysporum_fsp_lactucae/AJ516/publication/final_genes_appended_renamed_ncbi.pep.fasta)
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
  Taxon_code=FoLy
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

## 3.5.a Manual identification of numbers of orthologous and unique genes

```bash
for num in 1; do
echo "The number of ortholog groups shared between all isoaltes is:"
cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep -e 'FoLaR1|' | grep -e 'FoLaR4|' | grep -e 'FoC|' | grep -e 'FoLy|' | grep -e 'fo47|' |  wc -l
echo "This represented the following number of proteins:"
cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep -e 'FoLaR1|' | grep -e 'FoLaR4|' |  grep -e 'FoC|' | grep -e 'FoLy|' | grep -e 'fo47|' |  grep -o '|' | wc -l
echo "This represented the following number of proteins from FoLaR1:"
cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep -e 'FoLaR1|' | grep -e 'FoLaR4|' | grep -e 'FoC|' | grep -e 'FoLy|' | grep -e 'fo47|' |  grep -o 'FoLaR1' | wc -l
echo "This represented the following number of proteins from FoLaR4:"
cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep -e 'FoLaR1|' | grep -e 'FoLaR4|' | grep -e 'FoC|' | grep -e 'FoLy|' | grep -e 'fo47|' |  grep -o 'FoLaR4' | wc -l
echo "This Number of orthogroups unique to FoLaR1:"
cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep -e 'FoLaR1|' | grep -v -e 'FoC|' -e 'FoLy|' -e 'fo47|' -e 'FoLaR4' | wc -l
echo "This Number of orthogroups unique to FoLaR4:"
cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep -e 'FoLaR4|' | grep -v -e 'FoC|' -e 'FoLy|' -e 'fo47|' -e 'FoLaR1' | wc -l
echo "This Number of orthogroups unique and present in both Fo fsp Lactucae:"
cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep -e 'FoLaR1|' | grep -e 'FoLaR4|' | grep -v -e 'FoC|' -e 'FoLy|' -e 'fo47|' | wc -l
echo "This represented the following number of proteins from FoLaR1"
cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep -e 'FoLaR1|' | grep -e 'FoLaR4|' | grep -v -e 'FoC|' -e 'FoLy|' -e 'fo47|' |  grep -o 'FoLaR1|' | wc -l
echo "This represented the following number of proteins from FoLaR4"
cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep -e 'FoLaR1|' | grep -e 'FoLaR4|' | grep -v -e 'FoC|' -e 'FoLy|' -e 'fo47|' |  grep -o 'FoLaR4|' | wc -l
done
```

```
The number of ortholog groups shared between all isolates is::
11082
This represented the following number of proteins:
94649
This represented the following number of proteins from FoLaR1:
17129
This represented the following number of proteins from FoLaR4:
18035
This Number of orthogroups unique to FoLaR1:
12
This Number of orthogroups unique to FoLaR4:
57
This Number of orthogroups unique and present in both Fo fsp Lactucae:
520
This represented the following number of proteins from FoLaR1
638
This represented the following number of proteins from FoLaR4
660
```

```bash
echo "The number of orthogroups present:"
cat $WorkDir/"$IsolateAbrv"_orthogroups.tab | head -n1 | sed 's/ /\n/g' | wc -l
```

```
The number of orthogroups present:
19316
```

## 3.5.b Plot venn diagrams:

Using OrthoMCL output:

```bash
  IsolateAbrv=Fo_lactucae_FoL_Fo_FoC
  WorkDir=analysis/orthology/orthomcl/$IsolateAbrv
  ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/orthology
  $ProgDir/${IsolateAbrv}_venn_diag.r --inp $WorkDir/"$IsolateAbrv"_orthogroups.tab --out $WorkDir/"$IsolateAbrv"_orthogroups.pdf
```

Output was a pdf file of the venn diagram.

The following additional information was also provided. The format of the
following lines is as follows:

Isolate name (total number of orthogroups)
number of unique singleton genes
number of unique groups of inparalogs


```
[1] "FoLaR1"
[1] 13585
[1] 413
[1] 12
[1] "FoLaR4"
[1] 13864
[1] 580
[1] 57
[1] "FoC"
[1] 13523
[1] 684
[1] 34
[1] "FoLy"
[1] 14734
[1] 1598
[1] 345
[1] "fo47"
[1] 14090
[1] 1032
[1] 49
```
<!-- 
using orthofinder output:

```bash
  IsolateAbrv=Fo_lactucae_FoL_Fo_FoC
  WorkDir=analysis/orthology/orthomcl/$IsolateAbrv
  ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/orthology
  $ProgDir/${IsolateAbrv}_venn_diag.r --inp $WorkDir/formatted/*/ --out $WorkDir/"$IsolateAbrv"_orthogroups.pdf
``` -->


## 3.6) Extracting fasta files for all orthogroups


```bash
IsolateAbrv=Fo_lactucae_FoL_Fo_FoC
WorkDir=analysis/orthology/orthomcl/$IsolateAbrv
ProgDir=~/git_repos/emr_repos/tools/pathogen/orthology/orthoMCL
GoodProt=$WorkDir/goodProteins/goodProteins.fasta
OutDir=$WorkDir/orthogroups_fasta
mkdir -p $OutDir
$ProgDir/orthoMCLgroups2fasta.py --orthogroups $WorkDir/"$IsolateAbrv"_orthogroups.txt --fasta $GoodProt --out_dir $OutDir > $OutDir/extractionlog.txt
```

A combined dataset for nucleotide data was made for all gene models:

```bash
  mkdir -p $WorkDir/FTF
  nuc_file=$(ls gene_pred/final_genes/F.oxysporum_fsp_narcissi/N139_ncbi/*/final_genes_combined.gene.fasta)
  cat $nuc_file | sed "s/>/>FoN|/g" | cut -f1 -d ' ' > $WorkDir/goodProteins/nucleotide_seq.fa
  nuc_file=$(ls gene_pred/final_genes/F.oxysporum_fsp_cepae/Fus2_canu_new/*/final_genes_combined.gene.fasta)
  cat $nuc_file | sed "s/>/>FoC|/g" | cut -f1 -d ' ' >> $WorkDir/goodProteins/nucleotide_seq.fa
  nuc_file=$(ls assembly/external_group/F.oxysporum/fo47/broad/fusarium_oxysporum_fo47_1_genes.fasta)
  cat $nuc_file | sed "s/FOZG/fo47|FOZG/g" |cut -f1 -d ' ' >> $WorkDir/goodProteins/nucleotide_seq.fa
  nuc_file=$(ls assembly/external_group/F.oxysporum_fsp_lycopersici/4287/Fusox1/Fusox1_GeneCatalog_transcripts_20110522.nt.fasta)
  cat $nuc_file | sed "s/>.*FOXG/>FoL|FOXG/g" | cut -f1 -d ' ' >> $WorkDir/goodProteins/nucleotide_seq.fa
  GoodProt=$WorkDir/goodProteins/nucleotide_seq.fa
  OutDir=$WorkDir/orthogroups_fasta_nuc
  mkdir -p $OutDir
  cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | sed -r "s/-t\S+? / /g" | sed -r "s/.t. / /g" | sed -r "s/T. / /g" > $OutDir/"$IsolateAbrv"_orthogroups_mod.txt
  ProgDir=~/git_repos/emr_repos/tools/pathogen/orthology/orthoMCL
  $ProgDir/orthoMCLgroups2fasta.py --orthogroups $OutDir/"$IsolateAbrv"_orthogroups_mod.txt --fasta $GoodProt --out_dir $OutDir > $OutDir/extractionlog.txt
  rm $OutDir/"$IsolateAbrv"_orthogroups_mod.txt
```


## Also try using orthofinder

```bash
qlogin -pe smp 16 -l virtual_free=1G -l h=blacklace08.blacklace

#16 threads used
ProjDir=/home/groups/harrisonlab/project_files/fusarium
cd $ProjDir
IsolateAbrv=Fo_lactucae_FoL_Fo_FoC
WorkDir=analysis/orthology/orthomcl/$IsolateAbrv
orthofinder -f $WorkDir/formatted -t 16 -a 16
```

orthofinder results:

```

```

output files are in:
```bash
ls $WorkDir/formatted/Results_Apr10
```

```bash
OutDir=$(ls -d analysis/orthology/orthomcl/Pcac_Pinf_publication)
Orthogroups=$(ls $OutDir/Pcac_Pinf_publication_orthogroups.txt)
ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/orthology
$ProgDir/summarise_orthogroups.py --orthogroups $Orthogroups > $OutDir/summarised_orthogroups.tsv
```


## 3.7) Extracting ortholog groups for SIX gene homologs:

```bash
IsolateAbrv=FoN_vs_FoC_vs_FoL_vs_Fo
WorkDir=analysis/orthology/orthomcl/$IsolateAbrv
cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep -w -e 'FoN|NS_09200' -e 'FoN|g17549' -e 'FoN|g16592' -e 'FoN|g16591'
```


## 3.7) Extracting ortholog groups for secreted MIMP genes:

This script still has problems with Fo47 genes not being extracted from the commands above

```bash
  IsolateAbrv=FoN_vs_FoC_vs_FoL_vs_Fo
  WorkDir=analysis/orthology/orthomcl/$IsolateAbrv
  OutDir=$WorkDir/sec_mimps
  mkdir $OutDir
  cat gene_pred/annotation/F.oxysporum_fsp_narcissi/N139_ncbi/N139_ncbi_annotation_ncbi_MIMP_secreted.tsv | cut -f22 | grep 'orthogroup' > $OutDir/N139_ncbi_MIMP_secreted_orthogroups.txt

  for Orthogroup in $(cat $OutDir/N139_ncbi_MIMP_secreted_orthogroups.txt); do
    echo $Orthogroup
    cp $WorkDir/orthogroups_fasta_nuc/$Orthogroup.fa $OutDir/.
  done
```
