
# FoC genomes vs FoL Whole genome alignment


# 1. Whole Genome Alignment Programs

##  1.1. Mauve alignment

Progressive mauve was used to align all Fusarium oxysporum genomes.

Note - For the commands below to work you must have logged into the cluster using
        'ssh -X'.

### 1.1.a Ordering contigs in reference to FoL

```bash
  ProjDir=/home/groups/harrisonlab/project_files/fusarium
  WorkDir=$ProjDir/analysis/genome_alignment/mauve
  mkdir -p $WorkDir
  Reference=$(ls $ProjDir/repeat_masked/F.oxysporum_fsp_cepae/Fus2/filtered_contigs_repmask/Fus2_contigs_softmasked.fa)
  # Use move_contigs to order genomes based on reference for each phylogroup
  for Assembly in $( ls $ProjDir/repeat_masked/F.oxysporum*/*/filtered_contigs_repmask/*_contigs_softmasked.fa $ProjDir/assembly/external_group/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl/Fusarium_oxysporum_chromosome_and_additional_contigs.fa $ProjDir/assembly/external_group/F.oxysporum/fo47/broad/fusarium_oxysporum_fo47_1_supercontigs.fasta | grep -v -e 'Fus2' -e 'pisi' -e 'narcissi' -e 'HB17'); do
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

### 1.1.b Running Progressive Mauve

```bash
  ProjDir=/home/groups/harrisonlab/project_files/fusarium
  WorkDir=$ProjDir/analysis/genome_alignment/mauve
  OutDir=$WorkDir/alignment
  GenomeList=""
  # GenomeList=$(ls $ProjDir/assembly/external_group/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl/Fusarium_oxysporum_chromosome_and_additional_contigs.fa $ProjDir/assembly/external_group/F.oxysporum/fo47/broad/fusarium_oxysporum_fo47_1_supercontigs.fasta)
  for AlignDir in $(ls -d $WorkDir/*_contigs); do
    # echo "$AlignDir"
    NumAlignments=$(ls -d $AlignDir/alignment* | wc -l)
    AlignedContigs=$(ls $AlignDir/alignment"$NumAlignments"/*_contigs_softmasked.fa)
    echo $AlignedContigs
    GenomeList="$GenomeList""$AlignedContigs "
  done
  ProgDir=~/git_repos/emr_repos/tools/seq_tools/genome_alignment/mauve
  qsub $ProgDir/run_progressive_mauve.sh $OutDir "$GenomeList"
```


## 1.2 TBA alignment


Genome alignment was performed using the pipeline recomended by Eva Stukenbrock
at the Max Planck Institute for Evolutionary Biology, Kiel.

The TBA package http://www.bx.psu.edu/miller_lab/ was used to generate multiple
alignments.

* Note - It is important that fasta files and isolate names are formatted properly. Isolate names must begin with A-Z characters and can not contain '.' or '-' characters. The fasta file must match the name of the isolate exactly (no file extensions) and all headers within the fasta must be formatted in particular way, including starting with the isolate name.


### 1.2.a Parsing fasta file headers



```bash
  ProjDir=/home/groups/harrisonlab/project_files/fusarium
  WorkDir=$ProjDir/analysis/genome_alignment/TBA
  mkdir -p $WorkDir
  cd $WorkDir
  Reference=$(ls $ProjDir/repeat_masked/F.oxysporum_fsp_cepae/Fus2/filtered_contigs_repmask/Fus2_contigs_softmasked.fa)
  for Assembly in $(ls $ProjDir/repeat_masked/F.oxysporum_fsp_cepae/*/filtered_contigs_repmask/*_contigs_softmasked.fa $ProjDir/assembly/external_group/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl/Fusarium_oxysporum_chromosome_and_additional_contigs.fa $ProjDir/assembly/external_group/F.oxysporum/fo47/broad/fusarium_oxysporum_fo47_1_supercontigs.fasta | grep -v -e 'edited' -e 'pisi' -e 'narcissi' -e 'HB17'); do
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $Assembly| rev | cut -f3 -d '/' | rev | sed 's/\./_/g' | sed 's/-/_/g')
    # The commands to define strain also replace any '.' with a '-' character.
    cat $Assembly | sed "s/>/>Fo_$Strain:/g" > Fo_$Strain.fa
    get_standard_headers Fo_$Strain.fa > Fo_"$Strain"_headers.txt
    ProgDir=~/git_repos/emr_repos/tools/pathogen/lineage_specific_regions
    $ProgDir/parse_tba_headers.py --inp_fasta Fo_$Strain.fa --new_headers Fo_"$Strain"_headers.txt > Fo_$Strain
  done
```


### 1.2.b generating a series of pair-wise alignments

generating a series of pair-wise alignments to “seed” the multiple alignment process


```bash
  all_bz - \
    "(((Fo_125 Fo_Fus2 Fo_A23 Fo_55 Fo_HB6) Fo_A1_2 Fo_PG Fo_CB3 Fo_A28 Fo_D2 Fo_fo47) (Fo_4287) (Fo_A13))" >& all_bz.log

  cat all_bz.log | grep 'blastzWrapper' > commands_part1.log
  while read Commands; do
    Jobs=$(qstat | grep 'blastz' | grep 'qw' | wc -l)
    while [ $Jobs -gt 1 ]; do
      sleep 10
      printf "."
      Jobs=$(qstat | grep 'blastz' | grep 'qw' | wc -l)
    done
    Nodes="blacklace03.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace"
    qsub -S /bin/bash -b y -pe smp 1 -l virtual_free=1G -N blastzWrapper -l h="$Nodes" -cwd "$Commands"
  done < commands_part1.log
```

Once the lastz part 1 jobs had finished part 2 could be submitted:

```bash
  cat all_bz.log | grep 'single_cov2' > commands_part2.log
  Nodes="blacklace03.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace"
  while read Commands; do
    Jobs=$(qstat | grep 'single_cov2' | grep 'qw' | wc -l)
    while [ $Jobs -gt 1 ]; do
      sleep 10
      printf "."
      Jobs=$(qstat | grep 'single_cov2' | grep 'qw' | wc -l)
    done
    qsub -S /bin/bash -b y -pe smp 1 -l virtual_free=1G -N single_cov2 -l h="$Nodes" -cwd "$Commands"
  done < commands_part2.log
```

### 1.2.c Generating the multiple alignment

```
  tba E=A1177 "((Alt_648 Alt_1082 Alt_1164 Alt_24350 (Alt_635 Alt_743 Alt_1166 Alt_1177)) (Alt_675 Alt_97-0013 Alt_97-0016) (Alt_650))" A*.maf tba.maf >& tba.log
```


# 2. Alignment of raw reads vs the Fus2 genome

```bash
  Reference=$(ls repeat_masked/*/*/*/*_contigs_unmasked.fa | grep -w 'Fus2_canu_new')
  for StrainPath in $(ls -d qc_dna/paired/F.oxysporum_fsp_cepae/* | grep -v 'HB6' | grep -v 'HB17' | grep -v 'Fus2'); do
    Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
    Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
    echo "$Organism - $Strain"
    F_Read=$(ls $StrainPath/F/*.fq.gz)
    R_Read=$(ls $StrainPath/R/*.fq.gz)
    echo $F_Read
    echo $R_Read
    OutDir=analysis/genome_alignment/bowtie/$Organism/$Strain/vs_Fus2_unmasked
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment
    qsub $ProgDir/bowtie/sub_bowtie.sh $Reference $F_Read $R_Read $OutDir $Strain
  done
```

```bash
  Reference=$(lls repeat_masked/*/*/*/*_contigs_unmasked.fa | grep -w 'Fus2_canu_new')
  for StrainPath in $(ls -d qc_dna/paired/F.*/* | grep -e 'HB6' -e 'Fus2'); do
    echo $StrainPath
    Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
    Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
    echo "$Organism - $Strain"
    F1_Read=$(ls $StrainPath/F/*_trim.fq.gz | head -n1 | tail -n1);
    R1_Read=$(ls $StrainPath/R/*_trim.fq.gz | head -n1 | tail -n1);
    F2_Read=$(ls $StrainPath/F/*_trim.fq.gz | head -n2 | tail -n1);
    R2_Read=$(ls $StrainPath/R/*_trim.fq.gz | head -n2 | tail -n1);
    echo $F1_Read
    echo $R1_Read
    echo $F2_Read
    echo $R2_Read
    OutDir=analysis/genome_alignment/bowtie/$Organism/$Strain/vs_Fus2
    qsub $ProgDir/bowtie/sub_bowtie_2lib.sh $Reference $F1_Read $R1_Read $F2_Read $R2_Read $OutDir $Strain
  done
```
