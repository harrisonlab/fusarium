
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


## 4.5 Plot venn diagrams:

```bash
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/orthology
  $ProgDir/Fp_Fv_FoC_FoL_Fo_venn_diag.r --inp $WorkDir/"$IsolateAbrv"_orthogroups.tab --out $WorkDir/"$IsolateAbrv"_orthogroups.pdf
```

Output was a pdf file of the venn diagram.

The following additional information was also provided. The format of the
following lines is as follows:

Isolate name (total number of orthogroups)
number of unique singleton genes
number of unique groups of inparalogs

```
[1] "Fp"
[1] 12110
[1] 991
[1] 62
[1] "Fv"
[1] 11521
[1] 1248
[1] 15
[1] "FoC"
[1] 13516
[1] 914
[1] 52
[1] "FoL"
[1] 14738
[1] 1650
[1] 438
[1] "Fo"
[1] 14078
[1] 1087
[1] 61
```


## 4.6 Manual identification of numbers of orthologous and unique genes



```bash
for num in 1; do
echo "The number of ortholog groups common to all F. oxysporum isolates are:"
cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep 'Fp|' | grep 'Fv|' | grep 'FoC|' | grep 'FoL|' | grep 'Fo|' | wc -l
echo "This contained the following number of genes"
cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep 'Fp|' | grep 'Fv|' | grep 'FoC|' | grep 'FoL|' | grep 'Fo|' | grep -o '|' | wc -l
echo "This represented the following numebr of genes from Fp"
cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep 'Fp|' | grep 'Fv|' | grep 'FoC|' | grep 'FoL|' | grep 'Fo|' | grep -o 'Fp|' | wc -l
echo "The number of proteins in Fp-unique orthogroups are:"
cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep 'Fp|' | grep -v -e 'Fv|' -e 'FoC|' -e 'FoL|' -e 'Fo|' | grep -o 'Fp|' | wc -l
echo "The number of proteins in a 1:1 relationship between Fp and FoL are:"
AnnotTab=gene_pred/annotations/F.proliferatum/A8_ncbi/A8_ncbi_gene_annotations.tab
cat $AnnotTab | grep 'Fp(1)' | grep 'FoL(1)' | wc -l
echo "These were present in the following number of orthogroups"
cat $AnnotTab | grep 'Fp(1)' | grep 'FoL(1)' | cut -f16 | sort | uniq |  wc -l
done
```

```
The number of ortholog groups common to all F. oxysporum isolates are:
9261
This contained the following number of proteins
75055
This represented the following numebr of proteins from Fp
12301
```

The numbers of effector candidates in different regions of the venn diagram was
determined:

For Fp unique proteins:
```bash
echo "The number of proteins in Fp-unique orthogroups are:"
cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep 'Fp|' | grep -v -e 'Fv|' -e 'FoC|' -e 'FoL|' -e 'Fo|' | grep -o 'Fp|' | wc -l
echo "The number of secreted effectorP proteins in Fp-unique orthogroups are:"
cat analysis/effectorP/F.proliferatum/A8_ncbi/F.proliferatum_A8_ncbi_EffectorP_secreted_headers.txt | cut -f1 > tmp.txt
cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep 'Fp|' | grep -v -e 'Fv|' -e 'FoC|' -e 'FoL|' -e 'Fo|' | grep -f tmp.txt
cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep 'Fp|' | grep -v -e 'Fv|' -e 'FoC|' -e 'FoL|' -e 'Fo|' | grep -o -f tmp.txt
echo "The number of secreted CAZY proteins in Fp-unique orthogroups are:"
CAZY=gene_pred/CAZY/F.proliferatum/A8_ncbi/A8_ncbi_CAZY_secreted_headers.txt
cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep 'Fp|' | grep -v -e 'Fv|' -e 'FoC|' -e 'FoL|' -e 'Fo|' | grep -f $CAZY
cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep 'Fp|' | grep -v -e 'Fv|' -e 'FoC|' -e 'FoL|' -e 'Fo|' | grep -o -f $CAZY
echo "The number of Antismash proteins in Fp-unique orthogroups are:"
Antismash=analysis/antismash/F.proliferatum/A8_ncbi/metabolite_cluster_gene_headers.txt
cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep 'Fp|' | grep -v -e 'Fv|' -e 'FoC|' -e 'FoL|' -e 'Fo|' | grep -f $Antismash
cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep 'Fp|' | grep -v -e 'Fv|' -e 'FoC|' -e 'FoL|' -e 'Fo|' | grep -o -f $Antismash
echo "The number of MIMP proteins in Fp-unique orthogroups are:"
MIMPs=analysis/mimps/F.proliferatum/A8_ncbi/A8_ncbi_prots_in_2kb_mimp.txt
cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep 'Fp|' | grep -v -e 'Fv|' -e 'FoC|' -e 'FoL|' -e 'Fo|' | grep -f $MIMPs
cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep 'Fp|' | grep -v -e 'Fv|' -e 'FoC|' -e 'FoL|' -e 'Fo|' | grep -o -f $MIMPs
```


For Fp FoC unique proteins:

```bash
echo "The number of Fp proteins in Fp FoC-unique orthogroups are:"
cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep 'Fp|' | grep 'FoC' | grep -v -e 'Fv|' -e 'FoL|' -e 'Fo|' | grep -o 'Fp|' | wc -l
echo "The number of secreted effectorP proteins in Fp FoC-unique orthogroups are:"
cat analysis/effectorP/F.proliferatum/A8_ncbi/F.proliferatum_A8_ncbi_EffectorP_secreted_headers.txt | cut -f1 | sed -e 's/^/Fp|/g' > tmp.txt
cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep 'Fp|' | grep 'FoC' | grep -v -e 'Fv|' -e 'FoL|' -e 'Fo|' | grep -f tmp.txt
cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep 'Fp|' | grep 'FoC' | grep -v -e 'Fv|' -e 'FoL|' -e 'Fo|' | grep -o -f tmp.txt
echo "The number of secreted CAZY proteins in Fp FoC-unique orthogroups are:"
cat gene_pred/CAZY/F.proliferatum/A8_ncbi/A8_ncbi_CAZY_secreted_headers.txt | cut -f1 | sed -e 's/^/Fp|/g' > tmp.txt
cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep 'Fp|' | grep 'FoC' | grep -v -e 'Fv|' -e 'FoL|' -e 'Fo|' | grep -f tmp.txt
cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep 'Fp|' | grep 'FoC' | grep -v -e 'Fv|' -e 'FoL|' -e 'Fo|' | grep -o -f tmp.txt
echo "The number of Antismash proteins in Fp FoC-unique orthogroups are:"
cat analysis/antismash/F.proliferatum/A8_ncbi/metabolite_cluster_gene_headers.txt | cut -f1 | sed -e 's/^/Fp|/g' > tmp.txt
cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep 'Fp|' | grep 'FoC' | grep -v -e 'Fv|' -e 'FoL|' -e 'Fo|' | grep -f tmp.txt
cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep 'Fp|' | grep 'FoC' | grep -v -e 'Fv|' -e 'FoL|' -e 'Fo|' | grep -o -f tmp.txt
echo "The number of MIMP proteins in Fp FoC-unique orthogroups are:"
cat analysis/mimps/F.proliferatum/A8_ncbi/A8_ncbi_prots_in_2kb_mimp.txt | cut -f1 | sed -e 's/^/Fp|/g' > tmp.txt
cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep 'Fp|' | grep 'FoC' | grep -v -e 'Fv|' -e 'FoL|' -e 'Fo|' | grep -f tmp.txt
cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep 'Fp|' | grep 'FoC' | grep -v -e 'Fv|' -e 'FoL|' -e 'Fo|' | grep -o -f tmp.txt
```
