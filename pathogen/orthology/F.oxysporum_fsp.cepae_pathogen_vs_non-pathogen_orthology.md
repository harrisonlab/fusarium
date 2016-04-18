
# For a comparison between pathogenic isolates and non-pathogenic isolates of FoC


# Methodology 1
<!--
```bash
  ProjDir=/home/groups/harrisonlab/project_files/fusarium
  cd $ProjDir
  IsolateAbrv=FoC_path_vs_non_path
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
  Isolate1=$(ls gene_pred/braker/F.oxysporum_fsp_cepae/125/*/125_augustus.aa)
  Isolate2=$(ls gene_pred/braker/F.oxysporum_fsp_cepae/A23/*/A23_augustus.aa)
  Isolate3=$(ls gene_pred/braker/F.oxysporum_fsp_cepae/Fus2/*/Fus2_augustus.aa)
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
  Isolate1=$(ls gene_pred/braker/F.oxysporum_fsp_cepae/A28/*/A28_augustus.aa)
  Isolate2=$(ls gene_pred/braker/F.oxysporum_fsp_cepae/D2/*/D2_augustus.aa)
  Isolate3=$(ls gene_pred/braker/F.oxysporum_fsp_cepae/PG/*/PG_augustus.aa)
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
  $ProgDir/venn_diag_2_way.r --inp $WorkDir/"$IsolateAbrv"_orthogroups.tab --out $WorkDir/"$IsolateAbrv"_orthogroups.pdf
```

Output was a pdf file of the venn diagram.

The following additional information was also provided. The format of the
following lines is as follows:

Isolate name (total number of orthogroups)
number of unique singleton genes
number of unique groups of inparalogs

```
[1] "NonP (7238)"
[1] 219
[1] 86
[1] "Path (7082)"
[1] 61
[1] 88
NULL
```

# Downstream analysis

Particular orthogroups were analysed for expansion in isolates.

This section details the commands used and the results observed.



### Pathogenic Fusarium unique gene families

Of the 7107 orthogroups identified in the analysis (not including singleton
genes - found in a single pathogen without paralogs) 88 orthogroups were only
found in pathogic isolates. Within these 88 orthogroups, 70 were represented by
genes from each of the three pathogens. This compares to 86 orthogroups only
found in non-pathogenic isolates, within which only 28 orthogroups were common
to each of the non-pathogens.

```bash
  less gene_pred/interproscan/F.oxysporum_fsp_cepae/Fus2/Fus2_interpro.gff3
  less gene_pred/interproscan/F.oxysporum_fsp_cepae/Fus2/Fus2_interproscan.tsv
  less -S $PathOrthogroupsFus2Anno
  PathgeneDir=analysis/orthology/orthomcl/FoC_path_vs_non_path/path_unique
  Orthogroups=analysis/orthology/orthomcl/FoC_path_vs_non_path/FoC_path_vs_non_path_orthogroups.txt
  mkdir -p $PathgeneDir
  PathOrthogroupsFus2=$PathgeneDir/Fus2_pathgene_orthogroups.txt
  cat $Orthogroups | grep -v 'NonP' | wc -l
  cat $Orthogroups | grep -v 'Path' | wc -l
  cat $Orthogroups | grep -v 'NonP' | grep 'Fus2' | grep 'A23' | grep '125' | wc -l
  cat $Orthogroups | grep -v 'Path' | grep 'A28' | grep 'D2' | grep 'PG' | wc -l
```

Fus2 predicted genes were used to investigate the 70 pathogen-uniqu orthogroups.
77 Fus2 genes were contained within these 70 orthogroups (orthogroups may
contain paralogs).

Gff tracks were extracted for the 77 genes. Furthermore interproscan annotations
were extracted for these genes.

```bash  
  cat $Orthogroups | grep -v 'NonP' | grep 'Fus2' | grep 'A23' | grep '125' | grep -P -o 'Fus2_g.*?\.t.' | sed 's/Fus2_//g' | sed 's/\.t.//g' > $PathOrthogroupsFus2
  PathOrthogroupsFus2Gff=$PathgeneDir/Fus2_pathgene_orthogroups.gff
  Fus2Genes=$(ls gene_pred/braker/F.oxysporum_fsp_cepae/Fus2/*/Fus2_augustus_preds.gtf)
  Fus2Annotations=gene_pred/interproscan/F.oxysporum_fsp_cepae/Fus2/Fus2_interproscan.tsv
  PathOrthogroupsFus2Gff=$PathgeneDir/Fus2_pathgene_orthogroups.gff
  PathOrthogroupsFus2Anno=$PathgeneDir/Fus2_pathgene_orthogroup_annotation.tsv
  cat $Fus2Genes | grep -w -f $PathOrthogroupsFus2 > $PathOrthogroupsFus2Gff
  cat $Fus2Annotations | grep -w -f $PathOrthogroupsFus2 > $PathOrthogroupsFus2Anno
```

The non-pathogen A28 was used to investigate non-pathogen unique genes.
Non-pathogens contain 28 orthogroups, within which the non-pathogen A28 has 34
genes.

```bash
  NonPathgeneDir=analysis/orthology/orthomcl/FoC_path_vs_non_path/common_genes
  NonPathOrthogroupsA28=$NonPathgeneDir/A28_common_orthogroup.txt
  NonPathOrthogroupsA28Anno=$NonPathgeneDir/A28_common_orthogroup_annotation.tsv
  mkdir -p $NonPathgeneDir
  cat $Orthogroups | grep -v 'Path' | grep 'A28' | grep 'D2' | grep 'PG' | grep -P -o 'A28_g.*?\.t.' | sed 's/A28_//g' | sed 's/\.t.//g' > $NonPathOrthogroupsA28
  NonPathOrthogroupsA28Gff=$PathgeneDir/A28_nonpathgene_orthogroups.gff
  A28Genes=$(ls gene_pred/braker/F.oxysporum_fsp_cepae/A28/*/A28_augustus_preds.gtf)
  A28Annotations=gene_pred/interproscan/F.oxysporum_fsp_cepae/A28/A28_interproscan.tsv
  NonPathOrthogroupsA28Gff=$NonPathgeneDir/A28_nonpathgene_orthogroups.gff
  NonPathOrthogroupsA28Anno=$NonPathgeneDir/A28_nonpathgene_orthogroup_annotation.tsv
  cat $A28Genes | grep -w -f $NonPathOrthogroupsA28 > $NonPathOrthogroupsA28Gff
  cat $A28Annotations | grep -w -f $NonPathOrthogroupsA28 > $NonPathOrthogroupsA28Anno
```


The interproscan annotations were studied.
 * Two heterokaryon incompatibility loci were observed in pathogen unique genes.
 * 130 genes were identified in Fus2 with a Heterokaron annotation
 * Two heterokaryon incompatibility loci were found in non-pathogen specific genes.

```bash
  cat $PathOrthogroupsFus2Anno | grep -i 'heterokaryon' | cut -f1 | sort | uniq | wc -l
  cat $Fus2Annotations | grep -i 'heterokaryon' | cut -f1 | sort | uniq | wc -l
  CommonOrthogroupsFus2=$PathgeneDir/Fus2_common_orthogroups.gff
  cat $Orthogroups | grep 'Fus2' | grep 'A23' | grep '125' | grep 'A28' | grep 'D2' | grep 'PG' | grep -P -o 'Fus2_g.*?\.t.' | sed 's/Fus2_//g' | sed 's/\.t.//g' > $CommonOrthogroupsFus2
  CommongeneDir=analysis/orthology/orthomcl/FoC_path_vs_non_path/common_genes
  CommonOrthogroupsFus2Anno=$CommongeneDir/Fus2_common_orthogroup_annotation.tsv
  cat $Fus2Annotations | grep -w -f $CommonOrthogroupsFus2 > $CommonOrthogroupsFus2Anno
```

```bash
  cat $CommonOrthogroupsFus2Anno | grep -i "heterokaryon" | wc -l
  HetDir=analysis/orthology/orthomcl/FoC_path_vs_non_path/het_loci
  mkdir -p $HetDir
  # Fus2HetGenes=$HetDir/Fus2HetGenes.txt
  cat $Fus2Annotations | grep -i 'heterokaryon' | cut -f1 | sort | uniq | sed 's/\.t.//g' | sed 's/ //g' > $HetDir/Fus2HetGenes.txt
  cat $HetDir/Fus2HetGenes.txt | sed 's/g/Fus2_g/g' > $HetDir/Fus2HetGenes_Ortho.txt
  cat $Orthogroups | grep -w -f $HetDir/Fus2HetGenes_Ortho.txt > $HetDir/Fus2HetOrthogroups.txt
  # cat $Fus2Annotations | grep -i 'heterokaryon' | cut -f1 | sort | uniq | sed 's/g/Fus2_g/g' > $Fus2HetGenes
  cat $Fus2Genes | grep -w -f $Fus2HetGenes > $HetDir/Fus2HetGenes.gff
```

The following was noted from studying ortholog groups:
 * Within the Fus2 genome, Het loci are dispesed on long and short contigs ($HetDir/Fus2HetGenes.gff viewed in geneious)
 * One of these Het loci is part of the lineage specific regions of the 2nd largest contig.
 * Within the Fus2 core genes, 63 have a heterokaryon incompatbility function.


```bash
  cat $Fus2Annotations | grep -i 'HMG' | cut -f1 | sort | uniq | wc -l
```

The annotations of HMG box domains were studied in the genome:
 * 20 Fus2 proteins were identified containing HMG domains.

### Secreted RxLR like proteins

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
 -->




# Methodology 2


```bash
  ProjDir=/home/groups/harrisonlab/project_files/fusarium
  cd $ProjDir
  IsolateAbrv=FoC_path_vs_non_path
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
  Fasta_file=$(ls gene_pred/codingquary/F.oxysporum_fsp_cepae/125/*/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  cat "$Taxon_code".fasta | sed -r 's/\./_/g' | sed -r "s/_t/.t/g" > $WorkDir/formatted/"$Taxon_code".fasta
```
<!--
### for FoC 55
```bash
  Taxon_code=55
  Fasta_file=$(ls gene_pred/braker/F.oxysporum_fsp_cepae/55/55_augustus.aa)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
``` -->

### for FoC A23
```bash
  Taxon_code=A23
  Fasta_file=$(ls gene_pred/codingquary/F.oxysporum_fsp_cepae/A23/*/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  cat "$Taxon_code".fasta | sed -r 's/\./_/g' | sed -r "s/_t/.t/g" > $WorkDir/formatted/"$Taxon_code".fasta
```

### for FoC A28
```bash
  Taxon_code=A28
  Fasta_file=$(ls gene_pred/codingquary/F.oxysporum_fsp_cepae/A28/*/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  cat "$Taxon_code".fasta | sed -r 's/\./_/g' | sed -r "s/_t/.t/g" > $WorkDir/formatted/"$Taxon_code".fasta
```
### for FoC D2
```bash
  Taxon_code=D2
  Fasta_file=$(ls gene_pred/codingquary/F.oxysporum_fsp_cepae/D2/*/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  cat "$Taxon_code".fasta | sed -r 's/\./_/g' | sed -r "s/_t/.t/g" > $WorkDir/formatted/"$Taxon_code".fasta
```
### for FoC Fus2
```bash
  Taxon_code=Fus2
  Fasta_file=$(ls gene_pred/codingquary/F.oxysporum_fsp_cepae/Fus2_edited_v2/*/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  cat "$Taxon_code".fasta | sed -r 's/\./_/g' | sed -r "s/_t/.t/g" > $WorkDir/formatted/"$Taxon_code".fasta
```
### for FoC PG
```bash
  Taxon_code=PG
  Fasta_file=$(ls gene_pred/codingquary/F.oxysporum_fsp_cepae/PG/*/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  cat "$Taxon_code".fasta | sed -r 's/\./_/g' | sed -r "s/_t/.t/g" > $WorkDir/formatted/"$Taxon_code".fasta
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
  echo "The number of ortholog groups unique to pathogens are:"
  cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep -v -e 'A28|' -e 'D2|' -e 'PG|'| grep 'Fus2|' | grep '125|' | grep 'A23|' |  wc -l
  echo "The number of ortholog groups unique to non-pathogens are:"
  cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep -v -e 'Fus2|' -e '125|' -e 'A23|' | grep 'A28|' | grep 'D2|' | grep 'PG|' |  wc -l
  echo "The number of ortholog groups common to all F. oxysporum isolates are:"
  cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep 'Fus2|' | grep '125|' | grep 'A23|' | grep 'A28|' | grep 'D2|' | grep 'PG|' | wc -l
```

```
  The number of ortholog groups unique to pathogens are:
  202
  The number of ortholog groups unique to non-pathogens are:
  72
  The number of ortholog groups common to all F. oxysporum isolates are:
  11275
```


## Plot venn diagrams:

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
[1] 11347 <- non-pathogen orthogroups
[1] 11477 <- pathogen orthogroups
  [[1] "A28"
  [1] "The total number of orthogroups and singleton genes in this isolate:  12573"
  [1] "The total number of orthogroups and singleton genes not in the venn diagram:  1226"
  [1] "The total number of singleton genes not in the venn diagram:  237"
  [1] "D2"
  [1] "The total number of orthogroups and singleton genes in this isolate:  12243"
  [1] "The total number of orthogroups and singleton genes not in the venn diagram:  896"
  [1] "The total number of singleton genes not in the venn diagram:  141"
  [1] "PG"
  [1] "The total number of orthogroups and singleton genes in this isolate:  12305"
  [1] "The total number of orthogroups and singleton genes not in the venn diagram:  958"
  [1] "The total number of singleton genes not in the venn diagram:  203"
  [1] "Fus2"
  [1] "The total number of orthogroups and singleton genes in this isolate:  12583"
  [1] "The total number of orthogroups and singleton genes not in the venn diagram:  1106"
  [1] "The total number of singleton genes not in the venn diagram:  75"
  [1] "125"
  [1] "The total number of orthogroups and singleton genes in this isolate:  12535"
  [1] "The total number of orthogroups and singleton genes not in the venn diagram:  1058"
  [1] "The total number of singleton genes not in the venn diagram:  65"
  [1] "A23"
  [1] "The total number of orthogroups and singleton genes in this isolate:  12555"
  [1] "The total number of orthogroups and singleton genes not in the venn diagram:  1078"
  [1] "The total number of singleton genes not in the venn diagram:  60"
```


### Pathogenic Fusarium unique gene families

Of the 13054 orthogroups identified by orthomcl (not including genes without a
BLAST hit), 297 orthogroups were only found in pathogenic isolates. Within these 297
orthogroups, 202 were represented by genes from each of the three pathogens. This compares to 278 orthogroups only
found in non-pathogenic isolates, within which only 72 orthogroups were common
to each of the non-pathogens.

```bash
  Orthogroups=$WorkDir/"$IsolateAbrv"_orthogroups.txt
  PathDir=$WorkDir/path_orthogroups
  NonPathDir=$WorkDir/non-path_orthogroups
  mkdir -p $WorkDir/path_orthogroups
  mkdir -p $WorkDir/non-path_orthogroups
  cat $Orthogroups | grep -v -w -e 'Fus2' -e 'A23' -e '125' | grep -w -e 'A28' -e 'D2' -e 'PG' | wc -l
  cat $Orthogroups | grep -v -w -e 'Fus2' -e 'A23' -e '125' | grep -w 'A28' | grep -w 'D2' | grep -w 'PG' > $NonPathDir/non-path_orthogroups.txt
  cat $NonPathDir/non-path_orthogroups.txt | wc -l
  cat $Orthogroups | grep -w -e 'Fus2' -e 'A23' -e '125' | grep -v -w -e 'A28' -e 'D2' -e 'PG' | wc -l
  cat $Orthogroups | grep -w 'Fus2' | grep -w 'A23' | grep -w '125' | grep -v -w -e 'A28' -e 'D2' -e 'PG' > $PathDir/path_orthogroups.txt
  cat $PathDir/path_orthogroups.txt | wc -l
```

Fus2 predicted genes were used to investigate the 202 pathogen-unique orthogroups.
231 Fus2 genes (233 incl. alternative transcripts) were contained within these
202 orthogroups (orthogroups may contain inparalogs).

Gff tracks were extracted for the 233 genes. Furthermore interproscan annotations
were extracted for these genes.

```bash  
  cat $Orthogroups | grep -w 'Fus2' | grep -w 'A23' | grep -w '125' | grep -v -w -e 'A28' -e 'D2' -e 'PG' | sort | uniq | wc -l
  cat $Orthogroups | grep -w 'Fus2' | grep -w 'A23' | grep -w '125' | grep -v -w -e 'A28' -e 'D2' -e 'PG' | grep -o -E 'Fus2\|g\w+\.\w+' | sed 's/Fus2|//g' | sed -E 's/.t\w+//g' | sort | uniq > $PathDir/Fus2_path_orthogroup_genes.txt
  cat $PathDir/Fus2_path_orthogroup_genes.txt | wc -l

  # Fasta accessions were extracted for pathogen genes:
  cat $Orthogroups | grep -w 'Fus2' | grep -w 'A23' | grep -w '125' | grep -v -w -e 'A28' -e 'D2' -e 'PG' | grep -o -E 'Fus2\|g\w+\.\w+' | sed 's/Fus2|//g' | sort | uniq > $PathDir/Fus2_path_orthogroup_transcripts.txt
  GeneModels=gene_pred/braker/F.oxysporum_fsp_cepae/Fus2/F.oxysporum_fsp_cepae_Fus2_braker/augustus.aa
  GeneList=$PathDir/Fus2_path_orthogroup_transcripts.txt
  OutFile=$PathDir/Fus2_path_orthogroup_genes.fa
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
  $ProgDir/extract_from_fasta.py --fasta $GeneModels --headers $GeneList > $OutFile
  rm $PathDir/Fus2_path_orthogroup_transcripts.txt


  # Exract gff annotations for pathogen genes
  Fus2Gff=gene_pred/braker/F.oxysporum_fsp_cepae/Fus2/F.oxysporum_fsp_cepae_Fus2_braker/augustus_extracted.gff
  PathOrthogroupsFus2Gff=$PathDir/Fus2_path_orthogroup_genes.gff
  echo "Gff annotations were extracted for the following number of genes/transcripts:"
  cat $Fus2Gff | grep -w -f $PathDir/Fus2_path_orthogroup_genes.txt > $PathOrthogroupsFus2Gff
  cat $PathOrthogroupsFus2Gff | grep -w 'gene' | wc -l
  cat $PathOrthogroupsFus2Gff | grep -w 'transcript' | wc -l

  # Extract interproscan predicted functions for pathogen genes
  Fus2Annotations=gene_pred/interproscan/F.oxysporum_fsp_cepae/Fus2/Fus2_interproscan.tsv
  PathOrthogroupsFus2Anno=$PathDir/Fus2_path_orthogroup_genes.tsv
  cat $Fus2Annotations | grep -w -f $PathDir/Fus2_path_orthogroup_genes.txt > $PathOrthogroupsFus2Anno
  echo "Functional annotations were extracted for the following number of transcripts:"
  cat $PathOrthogroupsFus2Anno | cut -f1 | sort | uniq | wc -l
  echo "This totaled the following number of annotations:"
  cat $PathOrthogroupsFus2Anno | wc -l

  # Identify number of pathogen genes that are secreted
  Fus2Secreted=gene_pred/augustus_signalp-4.1/F.oxysporum_fsp_cepae/Fus2/Fus2_aug_sp.tab
  PathOrthogroupsFus2Secreted=$PathDir/Fus2_path_orthogroup_secreted.txt
  echo "The following number of Fus2 pathogen unique transcripts were predicted to be secreted:"
  cat $Fus2Secreted | grep -w -f $PathDir/Fus2_path_orthogroup_genes.txt | grep -w 'YES' | cut -f1 | tr -d '>' > $PathOrthogroupsFus2Secreted
  cat $PathOrthogroupsFus2Secreted  | wc -l
  # Identify number of pathogen genes with coding signals implicating them as effectors
  Fus2EffectorP=analysis/effectorP/F.oxysporum_fsp_cepae/Fus2/F.oxysporum_fsp_cepae_Fus2_braker_augustus.aa_effectorP.txt
  PathOrthogroupsFus2EffectorP=$PathDir/Fus2_path_orthogroup_effectorP.txt
  echo "The following number of pathogen unique transcripts look like fungal effectors:"
  cat $Fus2EffectorP | grep -w -f $PathDir/Fus2_path_orthogroup_genes.txt | grep -v 'Effector probability' > $PathOrthogroupsFus2EffectorP
  cat $PathOrthogroupsFus2EffectorP | grep -w 'Effector' | wc -l
  #
  echo "The following number of genes are predicted as secreted and look like effectors:"
  PathOrthogroupsFus2SecretedEffectorP=$PathDir/Fus2_path_orthogroup_secreted_effectorP.txt
  cat $PathOrthogroupsFus2EffectorP $PathOrthogroupsFus2Secreted | grep -v 'Non-effector' | cut -f1 -d '.'| sort | uniq -d > $PathOrthogroupsFus2SecretedEffectorP
  cat $PathOrthogroupsFus2SecretedEffectorP | wc -l
  # The number of pathogen genes in 2kb of Mimps:
  echo "The number of pathogen genes in 2kb of Mimps:x"
  MimpGenesTxt=analysis/mimps/F.oxysporum_fsp_cepae/Fus2/Fus2_genes_in_2kb_mimp.txt
  Path
  PathOrthogroupsFus2Mimps=$PathDir/Fus2_path_orthogroup_mimp_genes.txt
  cat $MimpGenesTxt | grep -w -f $PathDir/Fus2_path_orthogroup_genes.txt > $PathOrthogroupsFus2Mimps
  cat $PathOrthogroupsFus2Mimps | wc -l
  echo "The number of pathogen common genes in the top 20 expressed Fus2 genes at 72hr"
  TopGenes=timecourse/2016_genes/Fus2/72hrs/cufflinks/Fus2_expressed_genes_top20.txt
  cat $TopGenes | grep -w -f $PathDir/Fus2_path_orthogroup_genes.txt | wc -l
```

```
  Gff annotations were extracted for the following number of genes/transcripts:
  231
  233
  Functional annotations were extracted for the following number of transcripts:
  125
  This totaled the following number of annotations:
  590
  The following number of Fus2 pathogen unique transcripts were predicted to be secreted:
  29
  The following number of pathogen unique transcripts look like fungal effectors:
  79
  The following number of genes are predicted as secreted and look like effectors:
  14
  The number of pathogen genes in 2kb of Mimps:
  14
```


The occurence of predicted functions in genes was investigated.

```bash
  echo "The following number of anontations were investigated:"
  cat $PathOrthogroupsFus2Anno | cut -f1,5,6 | wc -l
  echo "Taking into account redundancy in a gene, there were this many annotations:"
  cat $PathOrthogroupsFus2Anno | cut -f1,5,6 | sort | uniq | wc -l
  echo "The most common domains in pathogen genes were identified:"
  cat $PathOrthogroupsFus2Anno | cut -f1,5,6 | sort | uniq | cut -f2,3 | sort | uniq -c | sort -n -r | head -n 20
  echo "The number of annotations summarised is:"
  cat $PathOrthogroupsFus2Anno | cut -f1,5,6 | sort | uniq | cut -f2,3 | sort | uniq -c | sort -n -r | wc -l
```

```
The following number of anontations were investigated:
590
Taking into account redundancy in a gene, there were this many annotations:
410
The most common domains in pathogen genes were identified:
33 TMhelix	Region of a membrane-bound protein predicted to be embedded in the membrane.
24 Coil
 9 SSF56112	          <- Protein kinase superfamily
 6 SSF52540	          <- P-loop containing nucleoside triphosphate hydrolase
 5 SSF53474	          <- alpha/beta-Hydrolases superfamily
 5 SSF53098	          <- Ribonuclease H-like superfamily
 5 PTHR23272	        <- SUBFAMILY NOT NAMED - BED FINGER-RELATED family
 5 PF06985	Heterokaryon incompatibility protein (HET)
 5 G3DSA:3.40.50.1820
 5 G3DSA:1.10.510.10
 4 G3DSA:3.40.50.300
 3 SSF55797	          <- CAP domain (cysteine rich)
 3 SSF103473	        <- MFS general substrate transporter
 3 SM00220	Serine/Threonine protein kinases, catalytic domain
 3 SM00198	SCP / Tpx-1 / Ag5 / PR-1 / Sc7 family of extracellular domains.
 3 PTHR24148:SF16	    <- SUBFAMILY NOT NAMED (PTHR24148:SF16) - FAMILY NOT NAMED (PTHR24148)
 3 PTHR24148
 3 PTHR10334          <- CYSTEINE-RICH SECRETORY PROTEIN-RELATED (PTHR10334)
 3 PS50011	Protein kinase domain profile.
 3 PF00188	Cysteine-rich secretory protein family
The number of annotations summarised is:
248
```

#### Investigating Non-pathogen unique genes using PG:

PG predicted genes were used to investigate the 72 pathogen-unique orthogroups.
90 PG genes (90 incl. alternative transcripts) were contained within these
72 orthogroups (orthogroups may contain inparalogs).

Gff tracks were extracted for the 90 genes. Furthermore interproscan annotations
were extracted for these genes.

```bash
  ProjDir=/home/groups/harrisonlab/project_files/fusarium
  cd $ProjDir
  IsolateAbrv=FoC_path_vs_non_path
  WorkDir=analysis/orthology/orthomcl/$IsolateAbrv
  Orthogroups=$WorkDir/"$IsolateAbrv"_orthogroups.txt
  PathDir=$WorkDir/path_orthogroups
  NonPathDir=$WorkDir/non-path_orthogroups
```

```bash
  for i in 1; do
    cat $Orthogroups | grep -w 'PG' | grep -w 'A28' | grep -w 'D2' | grep -v -w -e 'Fus2' -e 'A23' -e '125' | sort | uniq | wc -l
    cat $Orthogroups | grep -w 'PG' | grep -w 'A28' | grep -w 'D2' | grep -v -w -e 'Fus2' -e 'A23' -e '125' |  grep -o -E 'PG\|g\w+\.\w+' | sed 's/PG|//g' | sed -E 's/.t\w+//g' | sort | uniq > $NonPathDir/PG_non-path_orthogroup_genes.txt
    cat $NonPathDir/PG_non-path_orthogroup_genes.txt | wc -l

    # Fasta accessions were extracted for non-pathogen genes:
    cat $Orthogroups | grep -w 'PG' | grep -w 'A28' | grep -w 'D2' | grep -v -w -e 'Fus2' -e 'A23' -e '125' | grep -o -E 'PG\|g\w+\.\w+' | sed 's/PG|//g' | sort | uniq > $NonPathDir/PG_non-path_orthogroup_transcripts.txt
    GeneModels=gene_pred/braker/F.oxysporum_fsp_cepae/PG/F.oxysporum_fsp_cepae_PG_braker/augustus.aa
    GeneList=$NonPathDir/PG_non-path_orthogroup_transcripts.txt
    OutFile=$NonPathDir/PG_non-path_orthogroup_genes.fa
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_from_fasta.py --fasta $GeneModels --headers $GeneList > $OutFile
    rm $NonPathDir/PG_non-path_orthogroup_transcripts.txt


    # Exract gff annotations for non-pathogen genes
    PGGff=gene_pred/braker/F.oxysporum_fsp_cepae/PG/F.oxysporum_fsp_cepae_PG_braker/augustus_extracted.gff
    NonPathOrthogroupsPGGff=$NonPathDir/PG_non-path_orthogroup_genes.gff
    echo "Gff annotations were extracted for the following number of genes/transcripts:"
    cat $PGGff | grep -w -f $NonPathDir/PG_non-path_orthogroup_genes.txt > $NonPathOrthogroupsPGGff
    cat $NonPathOrthogroupsPGGff | grep -w 'gene' | wc -l
    cat $NonPathOrthogroupsPGGff | grep -w 'transcript' | wc -l

    # Extract interproscan predicted functions for non-pathogen genes
    PGAnnotations=gene_pred/interproscan/F.oxysporum_fsp_cepae/Fus2/Fus2_interproscan.tsv
    NonPathOrthogroupsPGAnno=$NonPathDir/PG_non-path_orthogroup_genes.tsv
    cat $PGAnnotations | grep -w -f $NonPathDir/PG_non-path_orthogroup_genes.txt > $NonPathOrthogroupsPGAnno
    echo "Functional annotations were extracted for the following number of transcripts:"
    cat $NonPathOrthogroupsPGAnno | cut -f1 | sort | uniq | wc -l
    echo "This totaled the following number of annotations:"
    cat $NonPathOrthogroupsPGAnno | wc -l

    # Identify number of non-pathogen genes that are secreted
    PGSecreted=gene_pred/augustus_signalp-4.1/F.oxysporum_fsp_cepae/Fus2/Fus2_aug_sp.tab
    NonPathOrthogroupsPGSecreted=$NonPathDir/PG_non-path_orthogroup_secreted.txt
    echo "The following number of Fus2 pathogen unique transcripts were predicted to be secreted:"
    cat $PGSecreted | grep -w -f $NonPathDir/PG_non-path_orthogroup_genes.txt | grep -w 'YES' | cut -f1 | tr -d '>' > $NonPathOrthogroupsPGSecreted
    cat $NonPathOrthogroupsPGSecreted  | wc -l
    # Identify number of non-pathogen genes with coding signals implicating them as effectors
    PGEffectorP=analysis/effectorP/F.oxysporum_fsp_cepae/PG/F.oxysporum_fsp_cepae_PG_braker_augustus.aa_effectorP.txt
    NonPathOrthogroupsPGEffectorP=$NonPathDir/PG_non-path_orthogroup_effectorP.txt
    echo "The following number of pathogen unique transcripts look like fungal effectors:"
    cat $PGEffectorP | grep -w -f $NonPathDir/PG_non-path_orthogroup_genes.txt | grep -v 'Effector probability' > $NonPathOrthogroupsPGEffectorP
    cat $NonPathOrthogroupsPGEffectorP | grep -w 'Effector' | wc -l
    #
    echo "The following number of genes are predicted as secreted and look like effectors:"
    NonPathOrthogroupsPGSecretedEffectorP=$NonPathDir/PG_non-path_orthogroup_secreted_effectorP.txt
    cat $NonPathOrthogroupsPGEffectorP $NonPathOrthogroupsPGSecreted | grep -v 'Non-effector' | cut -f1 -d '.'| sort | uniq -d > $NonPathOrthogroupsPGSecretedEffectorP
    cat $NonPathOrthogroupsPGSecretedEffectorP | wc -l
    # The number of non-pathogen genes in 2kb of Mimps:
    echo "The number of pathogen genes in 2kb of Mimps:"
    MimpGenesTxt=analysis/mimps/F.oxysporum_fsp_cepae/PG/PG_genes_in_2kb_mimp.txt
    NonPathOrthogroupsPGMimps=$NonPathDir/PG_non-path_orthogroup_mimp_genes.txt
    cat $MimpGenesTxt | grep -w -f $NonPathDir/PG_non-path_orthogroup_genes.txt > $NonPathOrthogroupsPGMimps
    cat $NonPathOrthogroupsPGMimps | wc -l
    # echo "The number of pathogen common genes in the top 20 expressed Fus2 genes at 72hr"
    # TopGenes=timecourse/2016_genes/Fus2/72hrs/cufflinks/Fus2_expressed_genes_top20.txt
    # cat $TopGenes | grep -w -f $NonPathDir/PG_non-path_orthogroup_genes.txt | wc -l
  done
```

```
72
90
Gff annotations were extracted for the following number of genes/transcripts:
90
90
Functional annotations were extracted for the following number of transcripts:
84
This totaled the following number of annotations:
726
The following number of Fus2 pathogen unique transcripts were predicted to be secreted:
5
The following number of pathogen unique transcripts look like fungal effectors:
24
The following number of genes are predicted as secreted and look like effectors:
3
The number of pathogen genes in 2kb of Mimps:
0
```


The occurence of predicted functions in genes was investigated.

```bash
  for i in 1; do
    echo "The following number of anontations were investigated:"
    cat $NonPathOrthogroupsPGAnno | cut -f1,5,6 | wc -l
    echo "Taking into account redundancy in a gene, there were this many annotations:"
    cat $NonPathOrthogroupsPGAnno | cut -f1,5,6 | sort | uniq | wc -l
    echo "The most common domains in pathogen genes were identified:"
    cat $NonPathOrthogroupsPGAnno | cut -f1,5,6 | sort | uniq | cut -f2,3 | sort | uniq -c | sort -n -r | head -n 20
    echo "The number of annotations summarised is:"
    cat $NonPathOrthogroupsPGAnno | cut -f1,5,6 | sort | uniq | cut -f2,3 | sort | uniq -c | sort -n -r | wc -l
  done
```

```
The following number of anontations were investigated:
726
Taking into account redundancy in a gene, there were this many annotations:
494
The most common domains in pathogen genes were identified:
     20 Coil
     17 TMhelix	Region of a membrane-bound protein predicted to be embedded in the membrane.
      7 SSF52540
      7 G3DSA:3.40.50.300
      5 SSF51735
      5 G3DSA:3.40.50.720
      3 SSF57701
      3 SSF53335
      3 SSF51905
      3 SSF103473
      3 PTHR13789
      3 PTHR11711
      3 PS50048	Zn(2)-C6 fungal-type DNA-binding domain profile.
      3 PS00463	Zn(2)-C6 fungal-type DNA-binding domain signature.
      3 PR00420	Aromatic-ring hydroxylase (flavoprotein monooxygenase) signature
      3 PF01494	FAD binding domain
      3 PF00172	Fungal Zn(2)-Cys(6) binuclear cluster domain
      3 G3DSA:4.10.240.10
      3 G3DSA:3.50.50.60
      3 G3DSA:3.40.50.150
The number of annotations summarised is:
366

```


#### Orthogroups of FTF genes

The orthogroups containing FTF1 and FTF2 genes were identified. The top blast hit
of CDS from FTF1a, FTFb, FTFc and FTF2 against the Fus2 genome were found to be
Fus2 genes g8884 (FTF2, FTF1a, FTF1b) and g11787 (FTF1c). To identify FTF
inparalogs and orthologs between species the orthogroup containing these genes
was extracted. This identified a single orthogroup containing both identified
FTF genes.   

```bash
cat analysis/orthology/orthomcl/FoC_path_vs_non_path/FoC_path_vs_non_path_orthogroups.txt | grep 'Fus2|g8884'
cat analysis/orthology/orthomcl/FoC_path_vs_non_path/FoC_path_vs_non_path_orthogroups.txt | grep 'Fus2|g8884' | grep -o -E '\w+\|' | sort | uniq -c
```

```
  orthogroup1422: 125|g2735.t1 PG|g14594.t1 Fus2|g8884.t1 D2|g9595.t1 A23|g9320.t1 A28|g1094.t1 Fus2|g11787.t1 A23|g11656.t1 PG|g7400.t1 125|g13904.t1
  2 125|    <-pathogen
  2 A23|    <-pathogen
  1 A28|    <-non-pathogen
  1 D2|     <-non-pathogen
  2 Fus2|   <-pathogen
  2 PG|     <-non-pathogen
```

Orthogroup 1422 was extracted:


#### Orthogroups of top 20 expressed genes

```bash
Orthogroups=analysis/orthology/orthomcl/FoC_path_vs_non_path/FoC_path_vs_non_path_orthogroups.txt
ExpressDir=analysis/orthology/orthomcl/FoC_path_vs_non_path/top_expressed
mkdir -p $ExpressDir
Top20Expressed=$ExpressDir/orthogroups_top20_expressed.txt
cat analysis/orthology/orthomcl/FoC_path_vs_non_path/FoC_path_vs_non_path_orthogroups.txt | grep -w -e 'Fus2|g11792' -e 'Fus2|g4762' -e 'Fus2|g12630' -e 'Fus2|g10978' -e 'Fus2|g12344' -e 'Fus2|10974' -e 'Fus2|g10975' -e 'Fus2|g10973' -e 'Fus2|g4720' -e 'Fus2|g10346' -e 'Fus2|g11653' -e 'Fus2|g11849' -e 'Fus2|g1141' -e 'Fus2|g12149' -e 'Fus2|g2148' -e 'Fus2|g9529' -e 'Fus2|g10347' -e 'Fus2|g11958' -e 'Fus2|g13070' -e 'Fus2|g12418' >  $Top20Expressed

```


<!--
The non-pathogen A28 was used to investigate non-pathogen unique genes.
Non-pathogens contain 28 orthogroups, within which the non-pathogen A28 has 34
genes.

```bash
  NonPathgeneDir=analysis/orthology/orthomcl/FoC_path_vs_non_path/common_genes
  NonPathOrthogroupsA28=$NonPathgeneDir/A28_common_orthogroup.txt
  NonPathOrthogroupsA28Anno=$NonPathgeneDir/A28_common_orthogroup_annotation.tsv
  mkdir -p $NonPathgeneDir
  cat $Orthogroups | grep -v 'Path' | grep 'A28' | grep 'D2' | grep 'PG' | grep -P -o 'A28_g.*?\.t.' | sed 's/A28_//g' | sed 's/\.t.//g' > $NonPathOrthogroupsA28
  NonPathOrthogroupsA28Gff=$PathgeneDir/A28_nonpathgene_orthogroups.gff
  A28Genes=$(ls gene_pred/braker/F.oxysporum_fsp_cepae/A28/*/A28_augustus_preds.gtf)
  A28Annotations=gene_pred/interproscan/F.oxysporum_fsp_cepae/A28/A28_interproscan.tsv
  NonPathOrthogroupsA28Gff=$NonPathgeneDir/A28_nonpathgene_orthogroups.gff
  NonPathOrthogroupsA28Anno=$NonPathgeneDir/A28_nonpathgene_orthogroup_annotation.tsv
  cat $A28Genes | grep -w -f $NonPathOrthogroupsA28 > $NonPathOrthogroupsA28Gff
  cat $A28Annotations | grep -w -f $NonPathOrthogroupsA28 > $NonPathOrthogroupsA28Anno
``` -->




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

## Format fasta files

### Pathogenic isolates:

### for FoC 125
```bash
  Taxon_code=125
  Fasta_file=$(ls gene_pred/codingquary/F.oxysporum_fsp_cepae/125/*/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for FoC A13
```bash
  Taxon_code=A13
  Fasta_file=$(ls gene_pred/codingquary/F.oxysporum_fsp_cepae/A13/*/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for FoC A23
```bash
  Taxon_code=A23
  Fasta_file=$(ls  gene_pred/codingquary/F.oxysporum_fsp_cepae/A23/*/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for FoC Fus2
```bash
  Taxon_code=Fus2
  Fasta_file=$(ls gene_pred/codingquary/F.oxysporum_fsp_cepae/Fus2_edited_v2/*/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### Intermediate isolates:

### for FoC 55
```bash
  Taxon_code=55
  Fasta_file=$(ls gene_pred/codingquary/F.oxysporum_fsp_cepae/55/*/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for FoC A1/2
```bash
  Taxon_code=A1_2
  Fasta_file=$(ls gene_pred/codingquary/F.oxysporum_fsp_cepae/A1-2/*/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for FoC CB3
```bash
  Taxon_code=CB3
  Fasta_file=$(ls gene_pred/codingquary/F.oxysporum_fsp_cepae/CB3/*/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for FoC HB17
```bash
  Taxon_code=HB17
  Fasta_file=$(ls gene_pred/codingquary/F.oxysporum_fsp_cepae/HB17/*/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for FoC HB6
```bash
  Taxon_code=HB6
  Fasta_file=$(ls gene_pred/codingquary/F.oxysporum_fsp_cepae/HB6/*/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### Non-pathogenic isolates:

### for FoC A28
```bash
  Taxon_code=A28
  Fasta_file=$(ls gene_pred/codingquary/F.oxysporum_fsp_cepae/A28/*/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for FoC D2
```bash
  Taxon_code=D2
  Fasta_file=$(ls gene_pred/codingquary/F.oxysporum_fsp_cepae/D2/*/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for FoC PG
```bash
  Taxon_code=PG
  Fasta_file=$(ls gene_pred/codingquary/F.oxysporum_fsp_cepae/PG/*/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for Fo Fo47
```bash
  Taxon_code=fo47
  Fasta_file=$(ls gene_pred/codingquary/F.oxysporum/fo47/final/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### Non-FoC isolates:
<!--
### for FoN N139
```bash
  Taxon_code=N139
  Fasta_file=$(ls gene_pred/codingquary/F.oxysporum_fsp_narcissi/N139/*/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for FoP FOP5
```bash
  Taxon_code=FOP5
  Fasta_file=$(ls gene_pred/codingquary/F.oxysporum_fsp_pisi/FOP5/*/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```
-->
<!--
### for FoP L5
```bash
  Taxon_code=FOP5
  Fasta_file=$(ls gene_pred/codingquary/F.oxysporum_fsp_pisi/L5/*/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
``` -->
<!--
### for FoP PG3
```bash
  Taxon_code=PG3
  Fasta_file=$(ls gene_pred/codingquary/F.oxysporum_fsp_pisi/PG3/*/final_genes_combined.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for FP A8
```bash
  Taxon_code=A8
  Fasta_file=$(ls gene_pred/codingquary/F.proliferatum/A8/*/final_genes_combined.pep.fasta)
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
  echo "The number of ortholog groups unique to pathogens are:"
  cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep -v -e 'A28|' -e 'D2|' -e 'PG|'| grep 'Fus2|' | grep '125|' | grep 'A23|' |  wc -l
  echo "The number of ortholog groups unique to non-pathogens are:"
  cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep -v -e 'Fus2|' -e '125|' -e 'A23|' | grep 'A28|' | grep 'D2|' | grep 'PG|' |  wc -l
  echo "The number of ortholog groups common to all F. oxysporum isolates are:"
  cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep 'Fus2|' | grep '125|' | grep 'A23|' | grep 'A28|' | grep 'D2|' | grep 'PG|' | wc -l
```

```
  The number of ortholog groups unique to pathogens are:
  202
  The number of ortholog groups unique to non-pathogens are:
  72
  The number of ortholog groups common to all F. oxysporum isolates are:
  11275
```


## Plot venn diagrams:

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
[1] 11347 <- non-pathogen orthogroups
[1] 11477 <- pathogen orthogroups
  [[1] "A28"
  [1] "The total number of orthogroups and singleton genes in this isolate:  12573"
  [1] "The total number of orthogroups and singleton genes not in the venn diagram:  1226"
  [1] "The total number of singleton genes not in the venn diagram:  237"
  [1] "D2"
  [1] "The total number of orthogroups and singleton genes in this isolate:  12243"
  [1] "The total number of orthogroups and singleton genes not in the venn diagram:  896"
  [1] "The total number of singleton genes not in the venn diagram:  141"
  [1] "PG"
  [1] "The total number of orthogroups and singleton genes in this isolate:  12305"
  [1] "The total number of orthogroups and singleton genes not in the venn diagram:  958"
  [1] "The total number of singleton genes not in the venn diagram:  203"
  [1] "Fus2"
  [1] "The total number of orthogroups and singleton genes in this isolate:  12583"
  [1] "The total number of orthogroups and singleton genes not in the venn diagram:  1106"
  [1] "The total number of singleton genes not in the venn diagram:  75"
  [1] "125"
  [1] "The total number of orthogroups and singleton genes in this isolate:  12535"
  [1] "The total number of orthogroups and singleton genes not in the venn diagram:  1058"
  [1] "The total number of singleton genes not in the venn diagram:  65"
  [1] "A23"
  [1] "The total number of orthogroups and singleton genes in this isolate:  12555"
  [1] "The total number of orthogroups and singleton genes not in the venn diagram:  1078"
  [1] "The total number of singleton genes not in the venn diagram:  60"
```
