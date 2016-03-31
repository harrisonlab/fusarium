# 4287 Comparison

This document details the commands used to compare FoC genes or contigs to those
from scaffolded FoL assembly.

## 1) Contig homology by BLASTing predicted genes

### 1.1) BLASTing FoC (Fus2) ls genes vs FoL (4287) genome.

Pathogen unique genes were determined from orthology analysis - detailed in
this git repo under:
pathogen/orthology/F.oxysporum_fsp_cepae_pathogen_vs_non-pathogen_orthology.md

The fasta file of Fus2 pathogen unique genes is found in:

```bash
  Fus2PathFa=analysis/orthology/orthomcl/FoC_path_vs_non_path/path_orthogroups/Fus2_path_orthogroup_genes.fa
```

These genes were blasted against the FoL 4287 genome. This was downloaded and
located in:
```bash
  FoLGenomeFa=assembly/external_group/F.oxysporum_fsp_lycopersici/4287/Ma_et_al_2010/F.oxysporum_fsp.lycopersici_4287.fasta
```

Blast searching was performed using the commands:

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
  qsub $ProgDir/blast_pipe.sh $Fus2PathFa protein $FoLGenomeFa
```

Convert top blast hits into gff annotations

```bash
  BlastHitsCsv=analysis/blast_homology/F.oxysporum_fsp_lycopersici/4287/4287_Fus2_path_orthogroup_genes.fa_homologs.csv
  HitsGff=$(echo $BlastHitsCsv | sed  's/.csv/.gff/g')
  Column2=FoC_path_gene_homolog
  NumHits=1
  ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
  $ProgDir/blast2gff.pl $Column2 $NumHits $BlastHitsCsv > $HitsGff
```

<!-- Summary of results:

```bash
  ResultsCsv=analysis/blast_homology/F.oxysporum_fsp_lycopersici/4287/4287_Fus2_path_orthogroup_genes.fa_homologs.csv
  OutDir=analysis/blast_homology/F.oxysporum_fsp_lycopersici/4287
  SummarisedResultsCsv=$OutDir/4287_Fus2_path_orthogroup_genes.fa_homologs.csv
  mkdir -p $OutDir
  cat $ResultsCsv | cut -f1,37- > $SummarisedResultsCsv
``` -->

### 1.2) BLASTing Six genes vs FoL (4287) genome.

14 Six genes were Blasted against the FoL genome to identify the distribution of
SIX genes through LS chromosomes.

```bash
  FoLGenomeFa=assembly/external_group/F.oxysporum_fsp_lycopersici/4287/Ma_et_al_2010/F.oxysporum_fsp.lycopersici_4287.fasta
  Query=analysis/blast_homology/Fo_path_genes/Fo_path_genes_CRX.fa
  qsub $ProgDir/blast_pipe.sh $Query dna $FoLGenomeFa
```

Convert top blast hits into gff annotations

```bash
  BlastHitsCsv=analysis/blast_homology/F.oxysporum_fsp_lycopersici/4287/4287_Fo_path_genes_CRX.fa_homologs.csv
  HitsGff=$(echo $BlastHitsCsv | sed  's/.csv/.gff/g')
	Column2=SIX_homolog
  NumHits=1
  ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
  $ProgDir/blast2gff.pl $Column2 $NumHits $BlastHitsCsv > $HitsGff
```


## 2) Blast homology vs 4287 chromosomal sequences and parsing for plotting on LS chromosomes

Data was downloaded from ensembl

### 2.0) Copying across chromosomal data

```bash
  ProjDir=/home/groups/harrisonlab/project_files/fusarium
  SeqDir=$ProjDir/assembly/external_group/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl
  cd $SeqDir
  for File in  **2.gz; do
    NewName=$(echo $File | sed 's/ 2.gz/.gz/g');
    echo "$File -> $NewName";
    mv "$File" $NewName;
  done
  gunzip *.gz
  # Concatenate chromosomal fasta files:
  for num in $(seq 1 15); do
    File=$(ls Fusarium_oxysporum.FO2.31.dna.chromosome."$num".fa);
    cat $File;
  done > Fusarium_oxysporum.FO2.31.dna.chromosome.fa
```

### 2.1) BLASTing FoC (Fus2) ls genes vs FoL (4287) chromosomal sequences.

Pathogen unique genes were determined from orthology analysis - detailed in
this git repo under:
pathogen/orthology/F.oxysporum_fsp_cepae_pathogen_vs_non-pathogen_orthology.md

The fasta file of Fus2 pathogen unique genes is found in:

```bash
  Fus2PathFa=analysis/orthology/orthomcl/FoC_path_vs_non_path/path_orthogroups/Fus2_path_orthogroup_genes.fa
```

These genes were blasted against the FoL 4287 chromosomeal sequences. These were
prepared as detailed above:
```bash
  FoLGenomeFa=assembly/external_group/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl/Fusarium_oxysporum.FO2.31.dna.chromosome.fa
```

Blast searching was performed using the commands:

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
  qsub $ProgDir/blast_pipe.sh $Fus2PathFa protein $FoLGenomeFa
```

Convert top blast hits into gff annotations

```bash
  BlastHitsCsv=analysis/blast_homology/F.oxysporum_fsp_lycopersici/4287_chromosomal/4287_chromosomal_Fus2_path_orthogroup_genes.fa_homologs.csv
  HitsGff=$(echo $BlastHitsCsv | sed  's/.csv/.gff/g')
  Column2=FoC_path_gene_homolog
  NumHits=1
  ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
  $ProgDir/blast2gff.pl $Column2 $NumHits $BlastHitsCsv > $HitsGff
```

The number of blast hits to each chromosome were summarised:
```bash
  cat $HitsGff | cut -f1 | sort -n |  uniq -c
```

These should be divided by contig length to show relative enrichment for pathogen unique genes
```
  6 1
  8 2
  24 3
  11 4
  8 5
  54 6
  7 7
  7 8
  3 9
  4 10
  1 11
  5 12
  4 13
  31 14
  4 15
```

Results were sumarised:

```bash
  ResultsCsv=analysis/blast_homology/F.oxysporum_fsp_lycopersici/4287_chromosomal/4287_chromosomal_Fus2_path_orthogroup_genes.fa_homologs.csv
  OutDir=analysis/blast_homology/F.oxysporum_fsp_lycopersici/4287_chromosomal
  SummarisedResultsCsv=$OutDir/4287_chromosomal_Fus2_path_orthogroup_genes.fa_homologs_summarised.csv
  mkdir -p $OutDir
  cat $ResultsCsv | cut -f1,37- > $SummarisedResultsCsv
```

### 2.2) BLASTing Six genes vs FoL (4287) genome.

14 Six genes were Blasted against the FoL genome to identify the distribution of
SIX genes through LS chromosomes.

```bash
  FoLGenomeFa=assembly/external_group/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl/Fusarium_oxysporum.FO2.31.dna.chromosome.fa
  Query=analysis/blast_homology/Fo_path_genes/Fo_path_genes_CRX.fa
  qsub $ProgDir/blast_pipe.sh $Query dna $FoLGenomeFa
```

Convert top blast hits into gff annotations

```bash
  BlastHitsCsv=analysis/blast_homology/F.oxysporum_fsp_lycopersici/4287_chromosomal/4287_chromosomal_Fo_path_genes_CRX.fa_homologs.csv
  HitsGff=$(echo $BlastHitsCsv | sed  's/.csv/.gff/g')
	Column2=SIX_homolog
  NumHits=1
  ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
  $ProgDir/blast2gff.pl $Column2 $NumHits $BlastHitsCsv > $HitsGff
```


### 2.3) Intersecting BLAST hits against the FoL genome with FoL genes


For Fus2 Pathogen unique genes:
```bash
  FoLBlastHits=analysis/blast_homology/F.oxysporum_fsp_lycopersici/4287_chromosomal/4287_chromosomal_Fus2_path_orthogroup_genes.fa_homologs.gff
  FoLGenes=assembly/external_group/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl/Fusarium_oxysporum.FO2.31.chr.gff3
  FoLIntersect=analysis/blast_homology/F.oxysporum_fsp_lycopersici/4287_chromosomal/4287_chromosomal_Fus2_path_orthogroup_genes.fa_intersect.bed
  bedtools intersect -wo -a $FoLBlastHits -b $FoLGenes > $FoLIntersect
```

For SIX genes:

```bash
  FoLBlastHits=analysis/blast_homology/F.oxysporum_fsp_lycopersici/4287_chromosomal/4287_chromosomal_Fo_path_genes_CRX.fa_homologs.gff
  FoLGenes=assembly/external_group/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl/Fusarium_oxysporum.FO2.31.chr.gff3
  FoLIntersect=analysis/blast_homology/F.oxysporum_fsp_lycopersici/4287_chromosomal/4287_chromosomal_Fo_path_genes_CRX.fa_intersect.bed
  bedtools intersect -wo -a $FoLBlastHits -b $FoLGenes > $FoLIntersect
```


### 2.4 Creating an output table of blast results

For Fus2 Pathogen unique genes:

```bash
  ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes
  OutDir=analysis/blast_homology/F.oxysporum_fsp_lycopersici/4287_chromosomal
  Blast_csv=$OutDir/4287_chromosomal_Fus2_path_orthogroup_genes.fa_homologs.csv
  FoL_intersected=$OutDir/4287_chromosomal_Fus2_path_orthogroup_genes.fa_intersect.bed
  FoC_genes=gene_pred/braker/F.oxysporum_fsp_cepae/Fus2/F.oxysporum_fsp_cepae_Fus2_braker/augustus_extracted.gff
  Results_table=$OutDir/4287_chromosomal_Fus2_path_orthogroup_genes.tab

  $ProgDir/4287_comparison_blast_results_2_tab.py --blast_csv $Blast_csv --FoL_intersected_genes $FoL_intersected --FoC_genes_gff $FoC_genes > $Results_table
  # cat $Results_table | tail -n +2 | sort -n -k6 > $OutDir/4287_chromosomal_Fus2_path_orthogroup_genes_sorted.tab
```

For SIX genes:

```bash
  ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes
  OutDir=analysis/blast_homology/F.oxysporum_fsp_lycopersici/4287_chromosomal
  Blast_csv=$OutDir/4287_chromosomal_Fo_path_genes_CRX.fa_homologs.csv
  FoL_intersected=$OutDir/4287_chromosomal_Fo_path_genes_CRX.fa_intersect.bed
  Results_table=$OutDir/4287_chromosomal_Fo_path_genes_CRX.tab

  $ProgDir/Fo_path_genes_vs_FoL_chormosomal.py --blast_csv $Blast_csv --FoL_intersected_genes $FoL_intersected > $Results_table
  # cat $Results_table | tail -n +2 | sort -n -k6 > $OutDir/4287_chromosomal_Fus2_path_orthogroup_genes_sorted.tab
```


The final table for SIX genes also included information on the location of SIX
genes in Fus2. To get this information, the blast results of SIX genes against
Fus2 was used generated in the README of this project folder:

```bash
  FoBlastHits=analysis/blast_homology/F.oxysporum_fsp_cepae/Fus2/Fus2_Fo_path_genes_CRX.fa_homologs.gff
  Fus2Genes=gene_pred/braker/F.oxysporum_fsp_cepae/Fus2/F.oxysporum_fsp_cepae_Fus2_braker/augustus_extracted.gff
  Fus2Intersect=analysis/blast_homology/F.oxysporum_fsp_lycopersici/4287_chromosomal/Fus2_Fo_path_genes_CRX.fa_intersect.bed
  bedtools intersect -wo -a $FoBlastHits -b $Fus2Genes > $Fus2Intersect

  ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes
  OutDir=analysis/blast_homology/F.oxysporum_fsp_lycopersici/4287_chromosomal
  Blast_csv=analysis/blast_homology/F.oxysporum_fsp_cepae/Fus2/Fus2_Fo_path_genes_CRX.fa_homologs.csv
  Fus2Intersect=analysis/blast_homology/F.oxysporum_fsp_lycopersici/4287_chromosomal/Fus2_Fo_path_genes_CRX.fa_intersect.bed
  Results_table=$OutDir/Fus2_Fo_path_genes_CRX.tab
  $ProgDir/Fo_path_genes_vs_FoL_chormosomal.py --blast_csv $Blast_csv --FoL_intersected_genes $Fus2Intersect > $Results_table
```

Final table for SIX genes:

```bash
  Fus2Tab=$OutDir/Fus2_Fo_path_genes_CRX.tab
  FolTab=$OutDir/4287_chromosomal_Fo_path_genes_CRX.tab
  Results_table=$OutDir/4287_chromosomal_Fo_path_genes_CRX_final.tab
  paste $Fus2Tab $FolTab | sed 's/hit_FoL_contig/hit_FoC_contig/' | sed 's/FoL_gene_start/FoC_gene_start/' | sed 's/FoL_gene_end/FoC_gene_end/' | sed 's/FoL_strand/FoC_strand/' |  sed 's/FoL_gene_ID/FoC_gene_ID/' | cut -f-9,11- > $Results_table
```


## 2.5 Reciprocal blasting of hits

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
  BlastHits=analysis/blast_homology/F.oxysporum_fsp_lycopersici/4287_chromosomal/4287_chromosomal_Fus2_path_orthogroup_genes.fa_homologs.csv
  HeadCol=$((5+34))
  StartCol=$((11+34))
  StopCol=$((12+34))
  StrandCol=$((10+34))
  IdCol=1
  FastaFile=assembly/external_group/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl/Fusarium_oxysporum.FO2.31.dna.chromosome.fa
  ExtractedHits=analysis/blast_homology/F.oxysporum_fsp_lycopersici/4287_chromosomal/4287_chromosomal_Fus2_path_orthogroup_genes.fa_extracted_hits.fa
  $ProgDir/sequence_extractor.py --coordinates_file $BlastHits --header_column $HeadCol --start_column $StartCol --stop_column $StopCol --strand_column $StrandCol --id_column $IdCol --fasta_file $FastaFile > $ExtractedHits
```

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
  BlastHits=analysis/blast_homology/F.oxysporum_fsp_lycopersici/4287_chromosomal/4287_chromosomal_Fo_path_genes_CRX.fa_homologs.csv
  HeadCol=$((5+1))
  StartCol=$((11+1))
  StopCol=$((12+1))
  StrandCol=$((10+1))
  IdCol=1
  FastaFile=assembly/external_group/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl/Fusarium_oxysporum.FO2.31.dna.chromosome.fa
  ExtractedHits=analysis/blast_homology/F.oxysporum_fsp_lycopersici/4287_chromosomal/4287_chromosomal_Fo_path_genes_CRX.fa_extracted_hits.fa
  $ProgDir/sequence_extractor.py --coordinates_file $BlastHits --header_column $HeadCol --start_column $StartCol --stop_column $StopCol --strand_column $StrandCol --id_column $IdCol --fasta_file $FastaFile > $ExtractedHits
```

# Whole genome alignment

## Mauve alignment

Progressive mauve was used to align all Fusarium oxysporum genomes against the
FoL genome.

Note - For the commands below to work you must have logged into the cluster using
        'ssh -X'.

### Ordering contigs in reference to FoL

```bash
  ProjDir=/home/groups/harrisonlab/project_files/fusarium
  WorkDir=$ProjDir/analysis/genome_alignment/mauve
  mkdir -p $WorkDir
  Reference=$(ls $ProjDir/assembly/external_group/F.oxysporum_fsp_lycopersici/4287/ncbi/FoL_4287_chromosomes.fa)
  # Use move_contigs to order genomes based on reference for each phylogroup
  for Assembly in $(ls $ProjDir/repeat_masked/F.oxysporum*/*/filtered_contigs_repmask/*_contigs_softmasked.fa); do
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $Assembly| rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    MauveDir=~/prog/mauve/mauve_snapshot_2015-02-13
    OutDir=$WorkDir/"$Strain"_contigs
    ProgDir=~/git_repos/emr_repos/tools/seq_tools/genome_alignment/mauve
    rm -r $OutDir
    qsub $ProgDir/mauve_order_contigs.sh $MauveDir $Reference  $Assembly $OutDir
  done
```

### Running Progressive Mauve

```bash
  ProjDir=/home/groups/harrisonlab/project_files/fusarium
  WorkDir=$ProjDir/analysis/genome_alignment/mauve
  OutDir=$WorkDir/alignment
  GenomeList=$(ls $ProjDir/assembly/external_group/F.oxysporum_fsp_lycopersici/4287/ncbi/FoL_4287_chromosomes.fa)
  for Assembly in $(ls $ProjDir/repeat_masked/F.oxysporum*/*/filtered_contigs_repmask/*_contigs_softmasked.fa); do
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $Assembly| rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    NumAlignments=$(ls -d $WorkDir/"$Strain"_contigs/alignment* | wc -l)
    AlignedContigs=$(ls $WorkDir/"$Strain"_contigs/alignment"$NumAlignments"/"$Strain"_contigs_softmasked.fa.fas)
    echo $AlignedContigs
    GenomeList="$GenomeList ""$AlignedContigs "
  done
  ProgDir=~/git_repos/emr_repos/tools/seq_tools/genome_alignment/mauve
  qsub $ProgDir/run_progressive_mauve.sh $OutDir "$GenomeList"
```


## TBA alignment


Genome alignment was performed using the pipeline recomended by Eva Stukenbrock
at the Max Planck Institute for Evolutionary Biology, Kiel.

The TBA package http://www.bx.psu.edu/miller_lab/ was used to generate multiple
alignments.

* Note - It is important that fasta files and isolate names are formatted properly. Isolate names must begin with A-Z characters and can not contain '.' or '-' characters. The fasta file must match the name of the isolate exactly (no file extensions) and all headers within the fasta must be formatted in particular way, including starting with the isolate name.


### 2.0) Renaming FoL fasta headers

The fasta headers of the FoL genome were in genbank format. This required parsing
before commands below could work.

```bash
  ProjDir=/home/groups/harrisonlab/project_files/fusarium
  InFile=$ProjDir/assembly/external_group/F.oxysporum_fsp_lycopersici/4287/ncbi/FoL_4287_chromosomes.fa
  OutFile=$ProjDir/assembly/external_group/F.oxysporum_fsp_lycopersici/4287/ncbi/FoL_4287_chromosomes_renamed.fa
  cat $InFile | sed 's/>gi|213958560|gb|CM000589.1| Fusarium oxysporum f. sp. lycopersici 4287 chromosome />chromosome_/g' | sed 's/,.*//g' > $OutFile
```

### 2.1) generating a series of pair-wise alignments

generating a series of pair-wise alignments to “seed” the multiple alignment process

```bash
  ProjDir=/home/groups/harrisonlab/project_files/fusarium
  WorkDir=$ProjDir/analysis/genome_alignment/TBA
  mkdir -p $WorkDir
  cd $WorkDir
  Reference=$(ls $ProjDir/assembly/external_group/F.oxysporum_fsp_lycopersici/4287/ncbi/FoL_4287_chromosomes_renamed.fa)
  for Assembly in $(ls $ProjDir/repeat_masked/F.oxysporum_fsp_cepae/*/filtered_contigs_repmask/*_contigs_softmasked.fa $Reference); do
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $Assembly| rev | cut -f3 -d '/' | rev | sed 's/\./_/g' | sed 's/-/_/g')
    # The commands to define strain also replace any '.' with a '-' character.
    cat $Assembly | sed "s/>/>Fo_$Strain:/g" > Fo_$Strain.fa
    get_standard_headers Fo_$Strain.fa > Fo_"$Strain"_headers.txt
    ProgDir=~/git_repos/emr_repos/tools/pathogen/lineage_specific_regions
    $ProgDir/parse_tba_headers.py --inp_fasta Fo_$Strain.fa --new_headers Fo_"$Strain"_headers.txt > Fo_$Strain
  done


  all_bz - \
    "(((Fo_125 Fo_Fus2 Fo_A23 Fo_HB17 Fo_55 Fo_HB6) Fo_A1_2 Fo_PG Fo_CB3 Fo_A28 Fo_D2)(Fo_4287) (Fo_A13))" >& all_bz.log

  cat all_bz.log | grep 'blastzWrapper' > commands_part1.log
  while read Commands; do
    Nodes="blacklace02.blacklace|blacklace02.blacklace|blacklace03.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace"
    qsub -S /bin/bash -b y -pe smp 1 -l virtual_free=1G -N blastzWrapper -l h="$Nodes" -cwd "$Commands"
  done < commands_part1.log

  cat all_bz.log | grep 'single_cov2' > commands_part2.log
  Nodes="blacklace02.blacklace|blacklace02.blacklace|blacklace03.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace"
  while read Commands; do
    qsub -S /bin/bash -b y -pe smp 1 -l virtual_free=1G -N single_cov2 -l h="$Nodes" -cwd "$Commands"
  done < commands_part2.log

```

### 2.2) generating the multiple alignment

```
  tba E=A1177 "((Alt_648 Alt_1082 Alt_1164 Alt_24350 (Alt_635 Alt_743 Alt_1166 Alt_1177)) (Alt_675 Alt_97-0013 Alt_97-0016) (Alt_650))" A*.maf tba.maf >& tba.log
```
