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

# Whold genome alignment

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
