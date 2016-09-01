
# For a comparison between pathogenic isolates and non-pathogenic isolates of FoC



# Methodology 3


```bash
  ProjDir=/home/groups/harrisonlab/project_files/fusarium
  cd $ProjDir
  IsolateAbrv=FoC_vs_Fo_vs_FoL
  WorkDir=analysis/orthology/orthomcl/$IsolateAbrv
  mkdir -p $WorkDir
  mkdir -p $WorkDir/formatted
  mkdir -p $WorkDir/goodProteins
  mkdir -p $WorkDir/badProteins  
```

## 3.1 Format fasta files

### 3.1.a Pathogenic isolates:

### for FoC 125
```bash
  Taxon_code=125
  Fasta_file=$(ls gene_pred/final_genes/F.oxysporum_fsp_cepae/125/*/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for FoC A23
```bash
  Taxon_code=A23
  Fasta_file=$(ls  gene_pred/final_genes/F.oxysporum_fsp_cepae/A23/*/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for FoC Fus2
```bash
  Taxon_code=Fus2
  Fasta_file=$(ls gene_pred/final_genes/F.oxysporum_fsp_cepae/Fus2_canu_new/*/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### 3.1.b Intermediate isolates:

### for FoC 55
```bash
  Taxon_code=55
  Fasta_file=$(ls gene_pred/final_genes/F.oxysporum_fsp_cepae/55/*/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for FoC A1/2
```bash
  Taxon_code=A1_2
  Fasta_file=$(ls gene_pred/final_genes/F.oxysporum_fsp_cepae/A1-2/*/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for FoC CB3
```bash
  Taxon_code=CB3
  Fasta_file=$(ls gene_pred/final_genes/F.oxysporum_fsp_cepae/CB3/*/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```
<!--
HB17 was identified as a contaminated sequence.

### for FoC HB17
```bash
  Taxon_code=HB17
  Fasta_file=$(ls gene_pred/final_genes/F.oxysporum_fsp_cepae/HB17/*/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
``` -->

### for FoC HB6
```bash
  Taxon_code=HB6
  Fasta_file=$(ls gene_pred/final_genes/F.oxysporum_fsp_cepae/HB6/*/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### 3.1.c Non-pathogenic isolates:

### for FoC A13
```bash
  Taxon_code=A13
  Fasta_file=$(ls gene_pred/final_genes/F.oxysporum_fsp_cepae/A13/*/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for FoC A28
```bash
  Taxon_code=A28
  Fasta_file=$(ls gene_pred/final_genes/F.oxysporum_fsp_cepae/A28/*/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for FoC D2
```bash
  Taxon_code=D2
  Fasta_file=$(ls gene_pred/final_genes/F.oxysporum_fsp_cepae/D2/*/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for FoC PG
```bash
  Taxon_code=PG
  Fasta_file=$(ls gene_pred/final_genes/F.oxysporum_fsp_cepae/PG/*/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for Fo Fo47
```bash
  Taxon_code=fo47
  Fasta_file=$(ls assembly/external_group/F.oxysporum/fo47/broad/fusarium_oxysporum_fo47_1_proteins.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### 3.1.d Non-FoC isolates:
<!--
### for FoN N139
```bash
  Taxon_code=N139
  Fasta_file=$(ls gene_pred/final_genes/F.oxysporum_fsp_narcissi/N139/*/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for FoP FOP5
```bash
  Taxon_code=FOP5
  Fasta_file=$(ls gene_pred/final_genes/F.oxysporum_fsp_pisi/FOP5/*/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```
-->
<!--
### for FoP L5
```bash
  Taxon_code=FOP5
  Fasta_file=$(ls gene_pred/final_genes/F.oxysporum_fsp_pisi/L5/*/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
``` -->
<!--
### for FoP PG3
```bash
  Taxon_code=PG3
  Fasta_file=$(ls gene_pred/final_genes/F.oxysporum_fsp_pisi/PG3/*/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for FP A8
```bash
  Taxon_code=A8
  Fasta_file=$(ls gene_pred/final_genes/F.proliferatum/A8/*/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```
-->


### for FoL 4287

```bash
  Taxon_code=4287
  Fasta_file=$(ls assembly/external_group/F.oxysporum_fsp_lycopersici/4287/Fusox1/Fusox1_GeneCatalog_proteins_20110522.aa.fasta)
  Id_field=4
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
echo "The number of ortholog groups found in pathogen but absent in non-pathogens is:"
cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep -v -e 'A28|' -e 'D2|' -e 'PG|' -e 'A13|' -e 'fo47|' -e 'CB3|' | grep 'Fus2|' | grep '125|' | grep 'A23|' |  wc -l
echo "The number of ortholog groups unique to pathogens are:"
cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep -v -e 'A28|' -e 'D2|' -e 'PG|' -e 'A13|' -e 'fo47|' -e 'CB3|' | grep 'Fus2|' | grep '125|' | grep 'A23|' | grep -v -e '55|' -e 'A1_2|' -e 'HB6|' | wc -l
echo "The number of genes shared between pathogens and an intermediate isolate:"
cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep -v -e 'A28|' -e 'D2|' -e 'PG|' -e 'A13|' -e 'fo47|' -e 'CB3|' | grep 'Fus2|' | grep '125|' | grep 'A23|' | grep -e '55|' -e 'A1_2|' -e 'HB6|'| wc -l
echo "The number of genes shared between pathogens and the intermediate 55:"
cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep -v -e 'A28|' -e 'D2|' -e 'PG|' -e 'A13|' -e 'fo47|' -e 'CB3|' | grep 'Fus2|' | grep '125|' | grep 'A23|' | grep '55|' | wc -l
echo "Excluding 55, the other intermediate isolates contribute the" \
"following number of orthogroups to the shared pathogen/intermediate pool:"
cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep -v -e 'A28|' -e 'D2|' -e 'PG|' -e 'A13|' -e 'fo47|' -e 'CB3|'| grep 'Fus2|' | grep '125|' | grep 'A23|' | grep -e 'A1_2|' -e 'HB6|'| wc -l
echo "The number of ortholog groups unique to non-pathogens are:"
cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep -v -e 'Fus2|' -e '125|' -e 'A23|' | grep 'A13|' | grep 'A28|' | grep 'D2|' | grep 'PG|' | grep 'fo47|' | grep 'CB3|' | wc -l
echo "The number of ortholog groups common to all F. oxysporum isolates are:"
cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep 'Fus2|' | grep '125|' | grep 'A23|' | grep 'A28|' | grep 'D2|' | grep 'PG|' | grep 'A13|' | grep 'fo47|' | grep 'CB3|' | grep '55|' | grep 'A1_2|' | grep 'HB6|' | grep '4287|' |wc -l
done
```

```
  The number of ortholog groups found in pathogen but absent in non-pathogens is:
  142
  The number of ortholog groups unique to pathogens are:
  12
  The number of genes shared between pathogens and an intermediate isolate:
  130
  The number of genes shared between pathogens and the intermediate 55:
  130
  Excluding 55, the other intermediate isolates contribute the following number of orthogroups to the shared pathogen/intermediate pool:
  5
  The number of ortholog groups unique to non-pathogens are:
  46
  The number of ortholog groups common to all F. oxysporum isolates are:
  10115
```

The number of ortholog groups shared between FoC and FoL was identified:

```bash
  echo "The number of ortholog groups common to FoC and FoL are:"
  cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep 'Fus2|' | grep '125|' | grep 'A23|' | grep '4287|' | wc -l
  cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep -v -e 'A28|' -e 'D2|' -e 'PG|' -e 'A13|' -e 'fo47|' | grep 'Fus2|' | grep '125|' | grep 'A23|' | grep '4287|' |  wc -l
```

```
The number of ortholog groups common to FoC and FoL are:
10629
33
```

## 3.5.b Plot venn diagrams:

```bash
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/orthology
  $ProgDir/FoC_path_vs_non_path_venn_diag.r --inp $WorkDir/"$IsolateAbrv"_orthogroups.tab --out $WorkDir/"$IsolateAbrv"_orthogroups.pdf
```

Output was a pdf file of the venn diagram.

The following additional information was also provided. The format of the
following lines is as follows:

Isolate name (total number of orthogroups)
number of unique singleton genes
number of unique groups of inparalogs


```
[1] 11878 <- non-pathogen orthogroups (3 non-pathogens)
[1] 11681 <- pathogen orthogroups (3 pathogens)
[1] "A28"
[1] "The total number of orthogroups and singleton genes in this isolate:  13664"
[1] "The total number of orthogroups and singleton genes not in the venn diagram:  1786"
[1] "The total number of singleton genes not in the venn diagram:  369"
[1] "D2"
[1] "The total number of orthogroups and singleton genes in this isolate:  13206"
[1] "The total number of orthogroups and singleton genes not in the venn diagram:  1328"
[1] "The total number of singleton genes not in the venn diagram:  235"
[1] "PG"
[1] "The total number of orthogroups and singleton genes in this isolate:  13530"
[1] "The total number of orthogroups and singleton genes not in the venn diagram:  1652"
[1] "The total number of singleton genes not in the venn diagram:  436"
[1] "Fus2"
[1] "The total number of orthogroups and singleton genes in this isolate:  13735"
[1] "The total number of orthogroups and singleton genes not in the venn diagram:  2054"
[1] "The total number of singleton genes not in the venn diagram:  281"
[1] "125"
[1] "The total number of orthogroups and singleton genes in this isolate:  13930"
[1] "The total number of orthogroups and singleton genes not in the venn diagram:  2249"
[1] "The total number of singleton genes not in the venn diagram:  327"
[1] "A23"
[1] "The total number of orthogroups and singleton genes in this isolate:  13556"
[1] "The total number of orthogroups and singleton genes not in the venn diagram:  1875"
[1] "The total number of singleton genes not in the venn diagram:  240"
```

Generating a Venn diagram for Pathogens, Nonpathogens and FoL

```bash
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/orthology
  $ProgDir/FoC_path_vs_non_path_FoL_venn_diag.r --inp $WorkDir/"$IsolateAbrv"_orthogroups.tab --out $WorkDir/"$IsolateAbrv"_orthogroups_APS.pdf
```

#### 3.6.a Extracting fasta files orthogroups
```bash
  ProgDir=~/git_repos/emr_repos/tools/pathogen/orthology/orthoMCL
  OrthogroupTxt=analysis/orthology/orthomcl/FoC_vs_Fo_vs_FoL/FoC_vs_Fo_vs_FoL_orthogroups.txt
  GoodProt=analysis/orthology/orthomcl/FoC_vs_Fo_vs_FoL//goodProteins/goodProteins.fasta
  OutDir=analysis/orthology/orthomcl/FoC_vs_Fo_vs_FoL/fasta/all_orthogroups
  mkdir -p $OutDir
  $ProgDir/orthoMCLgroups2fasta.py --orthogroups $OrthogroupTxt --fasta $GoodProt --out_dir $OutDir > $OutDir/extractionlog.txt
```


#### 3.6.b Extracting Fus2 pathogen unique genes for circos



```bash
  IsolateAbrv=FoC_vs_Fo_vs_FoL
  WorkDir=analysis/orthology/orthomcl/$IsolateAbrv
  OrthogroupTxt=analysis/orthology/orthomcl/FoC_vs_Fo_vs_FoL/FoC_vs_Fo_vs_FoL_orthogroups.txt
  Fus2Gff=gene_pred/final_genes/F.oxysporum_fsp_cepae/Fus2/final/final_genes_appended.gff3
  mkdir -p $WorkDir/Fus2_genes
```

Genes in orthogroups unique to all pathogens were extracted:
```bash
  cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep -v -e 'A28|' -e 'D2|' -e 'PG|' -e 'A13|' -e 'fo47|' | grep 'Fus2|' | grep '125|' | grep 'A23|' | grep -v -e '55|' -e 'A1_2|' -e 'CB3|' -e 'HB6|' > $WorkDir/Fus2_genes/Fus2_path_shared_orthgroups.txt
  cat $WorkDir/Fus2_genes/Fus2_path_shared_orthgroups.txt | grep -o "Fus2|\w*\.t." | sed 's/Fus2|//g' | cut -f1 -d '.' > $WorkDir/Fus2_genes/Fus2_path_shared_genes.txt
  cat $Fus2Gff | grep -w 'gene' | grep -w -f $WorkDir/Fus2_genes/Fus2_path_shared_genes.txt > $WorkDir/Fus2_genes/Fus2_path_shared_genes.gff
```

Genes in orthogroups shared between all pathogens and intermediate isolate 55, and absent from non-pathogens were extracted:
```bash
  cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep -v -e 'A28|' -e 'D2|' -e 'PG|' -e 'A13|' -e 'fo47|' | grep 'Fus2|' | grep '125|' | grep 'A23|' | grep '55|' > $WorkDir/Fus2_genes/Fus2_path_55_shared_orthgroups.txt
  cat $WorkDir/Fus2_genes/Fus2_path_55_shared_orthgroups.txt | grep -o "Fus2|\w*\.t." | sed 's/Fus2|//g' | cut -f1 -d '.' > $WorkDir/Fus2_genes/Fus2_path_55_shared_genes.txt
  cat $Fus2Gff | grep -w 'gene' | grep -w -f $WorkDir/Fus2_genes/Fus2_path_55_shared_genes.txt > $WorkDir/Fus2_genes/Fus2_path_55_shared_genes.gff
```

Genes in orthogroups shared between all pathogens and at least one intermediate, excluding 55, and absent from non-pathogens were extracted:
```bash
  cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep -v -e 'A28|' -e 'D2|' -e 'PG|' -e 'A13|' -e 'fo47|' | grep 'Fus2|' | grep '125|' | grep 'A23|' | grep -e 'A1_2|' -e 'CB3|' -e 'HB6|' > $WorkDir/Fus2_genes/Fus2_path_inter_no_55_shared_orthgroups.txt
  cat $WorkDir/Fus2_genes/Fus2_path_inter_no_55_shared_orthgroups.txt | grep -o "Fus2|\w*\.t." | sed 's/Fus2|//g' | cut -f1 -d '.' > $WorkDir/Fus2_genes/Fus2_path_inter_no_55_shared_genes.txt
  cat $Fus2Gff | grep -w 'gene' | grep -w -f $WorkDir/Fus2_genes/Fus2_path_inter_no_55_shared_genes.txt > $WorkDir/Fus2_genes/Fus2_path_inter_no_55_shared_genes.gff
```

Genes in orthogroups found in pathogens but not in nonpathogens were extracted:
```bash
  cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep -v -e 'A28|' -e 'D2|' -e 'PG|' -e 'A13|' -e 'fo47|' | grep 'Fus2|' | grep '125|' | grep 'A23|' > $WorkDir/Fus2_genes/Fus2_path_not_non-path_orthgroups.txt
  cat $WorkDir/Fus2_genes/Fus2_path_not_non-path_orthgroups.txt| grep -o "Fus2|\w*\.t." | sed 's/Fus2|//g' | cut -f1 -d '.' > $WorkDir/Fus2_genes/Fus2_path_not_non-path_genes.txt
  cat $Fus2Gff | grep -w 'gene' | grep -w -f $WorkDir/Fus2_genes/Fus2_path_not_non-path_genes.txt > $WorkDir/Fus2_genes/Fus2_path_not_non-path_genes.gff
```


# Methodology 4


```bash
  ProjDir=/home/groups/harrisonlab/project_files/fusarium
  cd $ProjDir
  IsolateAbrv=FoC_vs_Fo_vs_FoL_publication
  WorkDir=analysis/orthology/orthomcl/$IsolateAbrv
  mkdir -p $WorkDir
  mkdir -p $WorkDir/formatted
  mkdir -p $WorkDir/goodProteins
  mkdir -p $WorkDir/badProteins  
```

## 4.1 Format fasta files

### 4.1.a Pathogenic isolates:

### for FoC 125
```bash
  Taxon_code=125
  Fasta_file=$(ls gene_pred/final_genes/F.oxysporum_fsp_cepae/125/*/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for FoC A23
```bash
  Taxon_code=A23
  Fasta_file=$(ls  gene_pred/final_genes/F.oxysporum_fsp_cepae/A23/*/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for FoC Fus2
```bash
  Taxon_code=Fus2
  Fasta_file=$(ls gene_pred/final_genes/F.oxysporum_fsp_cepae/Fus2_canu_new/*/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### 4.1.c Non-pathogenic isolates:

### for FoC A13
```bash
  Taxon_code=A13
  Fasta_file=$(ls gene_pred/final_genes/F.oxysporum_fsp_cepae/A13/*/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for FoC A28
```bash
  Taxon_code=A28
  Fasta_file=$(ls gene_pred/final_genes/F.oxysporum_fsp_cepae/A28/*/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for FoC CB3
```bash
  Taxon_code=CB3
  Fasta_file=$(ls gene_pred/final_genes/F.oxysporum_fsp_cepae/CB3/*/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for FoC PG
```bash
  Taxon_code=PG
  Fasta_file=$(ls gene_pred/final_genes/F.oxysporum_fsp_cepae/PG/*/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for Fo Fo47
```bash
  Taxon_code=fo47
  Fasta_file=$(ls assembly/external_group/F.oxysporum/fo47/broad/fusarium_oxysporum_fo47_1_proteins.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### 4.1.d Non-FoC isolates:

### for FoL 4287

```bash
  Taxon_code=4287
  Fasta_file=$(ls assembly/external_group/F.oxysporum_fsp_lycopersici/4287/Fusox1/Fusox1_GeneCatalog_proteins_20110522.aa.fasta)
  Id_field=4
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
  MergeHits="$IsolateAbrv"_blast.tab
  printf "" > $MergeHits
  for Num in $(ls $WorkDir/splitfiles/*.tab | rev | cut -f1 -d '_' | rev | sort -n); do
    File=$(ls $WorkDir/splitfiles/*_$Num)
    cat $File
  done > $MergeHits
```

## 4.4 Perform ortholog identification

```bash
  ProgDir=~/git_repos/emr_repos/tools/pathogen/orthology/orthoMCL
  MergeHits="$IsolateAbrv"_blast.tab
  GoodProts=$WorkDir/goodProteins/goodProteins.fasta
  qsub $ProgDir/qsub_orthomcl.sh $MergeHits $GoodProts 5
```

## 4.5.a Manual identification of numbers of orthologous and unique genes


```bash
  for num in 1; do
    echo "The number of ortholog groups found in pathogen but absent in non-pathogens is:"
    cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep -v -e 'A28|' -e 'PG|' -e 'A13|' -e 'fo47|' -e 'CB3|' | grep 'Fus2|' | grep '125|' | grep 'A23|' |  wc -l
    echo "The number of ortholog groups unique to pathogens are:"
    cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep -v -e 'A28|' -e 'PG|' -e 'A13|' -e 'fo47|' -e 'CB3|' | grep 'Fus2|' | grep '125|' | grep 'A23|' | wc -l
    echo "This represents the following number of genes:"
    cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep -v -e 'A28|' -e 'PG|' -e 'A13|' -e 'fo47|' -e 'CB3|' | grep 'Fus2|' | grep '125|' | grep 'A23|' | grep -o '|' | wc -l
    echo "This represents the following number of genes from Fus2:"
    cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep -v -e 'A28|' -e 'PG|' -e 'A13|' -e 'fo47|' -e 'CB3|' | grep 'Fus2|' | grep '125|' | grep 'A23|' | grep -o 'Fus2|' | wc -l
    echo "The number of genes in the largest orthogroup is:"
    cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep -v -e 'A28|' -e 'PG|' -e 'A13|' -e 'fo47|' -e 'CB3|' | grep 'Fus2|' | grep '125|' | grep 'A23|' | head -n1 | grep -o '|' | wc -l
    echo "The number of ortholog groups unique to non-pathogens are:"
    cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep -v -e 'Fus2|' -e '125|' -e 'A23|' | grep 'A13|' | grep 'A28|' | grep 'PG|' | grep 'fo47|' | grep 'CB3|' | wc -l
    echo "The number of ortholog groups common to all F. oxysporum isolates are:"
    cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep 'Fus2|' | grep '125|' | grep 'A23|' | grep 'A28|' | grep 'PG|' | grep 'A13|' | grep 'fo47|' | grep 'CB3|' | grep '4287|' |wc -l
  done
```

```
The number of ortholog groups found in pathogen but absent in non-pathogens is:
255
The number of ortholog groups unique to pathogens are:
255
The number of ortholog groups unique to non-pathogens are:
61
The number of ortholog groups common to all F. oxysporum isolates are:
10391
```

The number of ortholog groups shared between FoC and FoL was identified:

```bash
  echo "The number of ortholog groups common to FoC and FoL are:"
  cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep 'Fus2|' | grep '125|' | grep 'A23|' | grep '4287|' | wc -l
  cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep -v -e 'A28|' -e 'PG|' -e 'A13|' -e 'fo47|' -e 'CB3' | grep 'Fus2|' | grep '125|' | grep 'A23|' | grep '4287|' |  wc -l
```

```
  The number of ortholog groups common to FoC and FoL are:
  10916
  43
```

## 4.5.b Plot venn diagrams:

```bash
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/orthology
  $ProgDir/FoC_path_vs_non_path_venn_diag.r --inp $WorkDir/"$IsolateAbrv"_orthogroups.tab --out $WorkDir/"$IsolateAbrv"_orthogroups.pdf
```

Output was a pdf file of the venn diagram.

The following additional information was also provided. The format of the
following lines is as follows:

Isolate name (total number of orthogroups)
number of unique singleton genes
number of unique groups of inparalogs


```
  [1]  <- non-pathogen orthogroups (5 non-pathogens)
  [1]  <- pathogen orthogroups (3 pathogens)
  [1] "Fus2"
  [1] "The total number of orthogroups and singleton genes in this isolate:  13846"
  [1] "The total number of orthogroups and singleton genes not in the venn diagram:  2412"
  [1] "The total number of singleton genes not in the venn diagram:  293"
  [1] "125"
  [1] "The total number of orthogroups and singleton genes in this isolate:  13975"
  [1] "The total number of orthogroups and singleton genes not in the venn diagram:  2541"
  [1] "The total number of singleton genes not in the venn diagram:  342"
  [1] "A23"
  [1] "The total number of orthogroups and singleton genes in this isolate:  13616"
  [1] "The total number of orthogroups and singleton genes not in the venn diagram:  2182"
  [1] "The total number of singleton genes not in the venn diagram:  256"
  [1] "A13"
  [1] "The total number of orthogroups and singleton genes in this isolate:  13333"
  [1] "The total number of orthogroups and singleton genes not in the venn diagram:  1705"
  [1] "The total number of singleton genes not in the venn diagram:  647"
  [1] "fo47"
  [1] "The total number of orthogroups and singleton genes in this isolate:  14272"
  [1] "The total number of orthogroups and singleton genes not in the venn diagram:  2644"
  [1] "The total number of singleton genes not in the venn diagram:  1363"
  NULL
```


#### 4.6.a Extracting fasta files orthogroups
```bash
  ProgDir=~/git_repos/emr_repos/tools/pathogen/orthology/orthoMCL
  OrthogroupTxt=analysis/orthology/orthomcl/$IsolateAbrv/"$IsolateAbrv"_orthogroups.txt
  GoodProt=analysis/orthology/orthomcl/$IsolateAbrv/goodProteins/goodProteins.fasta
  OutDir=analysis/orthology/orthomcl/$IsolateAbrv/fasta/all_orthogroups
  mkdir -p $OutDir
  $ProgDir/orthoMCLgroups2fasta.py --orthogroups $OrthogroupTxt --fasta $GoodProt --out_dir $OutDir > $OutDir/extractionlog.txt
```

A combined dataset for nucleotide data was made for all gene models:

```bash

for nuc_file in $(ls gene_pred/final_genes/F.*/*/final/final_genes_combined.gene.fasta | grep -e 'Fus2' -e '125' -e 'A23' -e 'PG' -e 'A28' -e 'CB3' -e 'A13' -e 'PG'); do
  Strain=$(echo $nuc_file | rev | cut -f3 -d '/' | rev)
  cat analysis/orthology/orthomcl/FoC_vs_Fo_vs_FoL_publication/FoC_vs_Fo_vs_FoL_publication_orthogroups.txt | grep -e 'Fus2|g16859.t' -e 'Fus2|g10474.t' | sed 's/ /\n/g' | grep -v 'orthogroup' | sed 's/\.t*//g' | sed 's/T0//g' | grep "$Strain" | cut -f2 -d '|' > tmp.txt
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
  $ProgDir/extract_from_fasta.py --fasta $nuc_file  --headers tmp.txt > $WorkDir/FTF/"$Strain"_FTF.nuc
done

nuc_file=assembly/external_group/F.oxysporum/fo47/broad/fusarium_oxysporum_fo47_1_genes.fasta
cat $nuc_file | sed "s/FOZG/fo47|FOZG/g" >> $WorkDir/goodProteins/nucleotide_seq.fa
nuc_file=assembly/external_group/F.oxysporum_fsp_lycopersici/4287/Fusox1/Fusox1_GeneCatalog_transcripts_20110522.nt.fasta
cat $nuc_file | sed "s/FOXG/4287|FOXG/g" >> $WorkDir/goodProteins/nucleotide_seq.fa

```


The FTF ortholog group was extracted from the main table. Constituant genes were
extracted from the nucleotide file:
```bash
  mkdir -p $WorkDir/FTF
	cat analysis/orthology/orthomcl/FoC_vs_Fo_vs_FoL_publication/FoC_vs_Fo_vs_FoL_publication_orthogroups.txt | grep -e 'Fus2|g16859.t' -e 'Fus2|g10474.t' | sed 's/ /\n/g' | grep -v 'orthogroup' | sed 's/\.t*//g' | sed 's/T0//g' > $WorkDir/FTF/FTF_list.txt
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
  $ProgDir/extract_from_fasta.py --fasta $WorkDir/goodProteins/nucleotide_seq.fa --headers $WorkDir/FTF/FTF_list.txt > $WorkDir/FTF/orthogroup506_nuc.fa
```


<!-- #### 6.3) Extracting fasta files for all orthogroups

```bash
IsolateAbrv=FoC_vs_Fo_vs_FoL_publication
WorkDir=analysis/orthology/orthomcl/$IsolateAbrv
ProgDir=~/git_repos/emr_repos/tools/pathogen/orthology/orthoMCL
GoodProt=$WorkDir/goodProteins/goodProteins.fasta
OutDir=$WorkDir/orthogroups_fasta
mkdir -p $OutDir
$ProgDir/orthoMCLgroups2fasta.py --orthogroups $WorkDir/"$IsolateAbrv"_orthogroups.txt --fasta $GoodProt --out_dir $OutDir > $OutDir/extractionlog.txt
for File in $(ls -v $OutDir/orthogroup*.fa); do
cat $File | grep '>' | tr -d '> '
done > $OutDir/orthogroup_genes.txt
cat $GoodProt | grep '>' | tr -d '> ' | grep -v -f $OutDir/orthogroup_genes.txt > $OutDir/singleton_genes.txt
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $GoodProt --headers $OutDir/singleton_genes.txt > $OutDir/singleton_genes.fa
echo "The numbe of singleton genes extracted is:"
cat $OutDir/singleton_genes.fa | grep '>' | wc -l

``` -->
