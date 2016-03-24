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


# Mauve was used to align Fus2 to Fol genomes

Progressive mauve was used to align all Fusarium oxysporum genomes against the
FoL genome.

Note - For the commands below to work you must have logged into the cluster using
        'ssh -X'.

## Ordering contigs in reference to FoL

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

## Running Progressive Mauve

```bash

```
