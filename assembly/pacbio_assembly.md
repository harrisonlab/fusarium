
## Data extraction

```bash
  RawDatDir=/home/vicker/new_pacbio_data
  mkdir -p raw_dna/pacbio/F.oxysporum_fsp_cepae/Fus2
  cp $RawDatDir/Richard_Harrison_EMR.RH.ENQ-933.01rev01_S1R2_Fus2.tar.gz raw_dna/pacbio/F.oxysporum_fsp_cepae/Fus2/.
  cd raw_dna/pacbio/F.oxysporum_fsp_cepae/Fus2
  tar -zxvf Richard_Harrison_EMR.RH.ENQ-933.01rev01_S1R2_Fus2.tar.gz
  cd /home/groups/harrisonlab/project_files/fusarium

  OutDir=raw_dna/pacbio/F.oxysporum_fsp_cepae/Fus2/extracted
  mkdir -p $OutDir
  # for H5_File in $(ls raw_dna/pacbio/F.oxysporum_fsp_cepae/Fus2/*/Analysis_Results/*.bas.h5); do
  #   ReadSet=$(echo $H5_File | rev | cut -f3 -d '/' | rev)
  #   bash5tools.py $H5_File --outFile $OutDir/Fus2_pacbio_$ReadSet --outType fastq --readType unrolled --minLength 100
  # done
  cat raw_dna/pacbio/F.oxysporum_fsp_cepae/Fus2/*/Analysis_Results/*.subreads.fastq > $OutDir/concatenated_pacbio.fastq
```

## Assembly


### Canu assembly

```bash
  Reads=$(ls raw_dna/pacbio/F.oxysporum_fsp_cepae/Fus2/extracted/concatenated_pacbio.fastq)
  GenomeSz="50m"
  Strain=$(echo $Reads | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $Reads | rev | cut -f4 -d '/' | rev)
  Prefix="$Strain"_canu
  OutDir="assembly/canu-1.3/$Organism/$Strain"_canu
  ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/canu
  qsub $ProgDir/submit_canu.sh $Reads $GenomeSz $Prefix $OutDir
```

<!-- Default error rate is set at 0.025 for pacbio reads. However it is recomended to increase the error rate to 0.035 for coverage <20X and decrease it to 0.015 for reads with coverage >60X

```bash
  Reads=$(ls raw_dna/pacbio/F.oxysporum_fsp_cepae/Fus2/extracted/concatenated_pacbio.fastq)
  GenomeSz="50m"
  Strain=$(echo $Reads | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $Reads | rev | cut -f4 -d '/' | rev)
  Prefix="$Strain"_canu
  OutDir=assembly/canu/$Organism/"$Strain"_test
  ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/canu
  SpecFile=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/canu/canu_spec_sensitive.txt
  qsub $ProgDir/submit_canu.sh $Reads $GenomeSz $Prefix $OutDir $SpecFile

``` -->



Assembly stats were collected using quast

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly/canu/*/*/*_canu.contigs.fasta | grep -w 'Fus2'); do
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)  
    OutDir=assembly/canu/$Organism/$Strain/filtered_contigs
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```

Contigs were renamed in accordance with ncbi recomendations

```bash
  ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
  touch tmp.csv
  printf "contig_17\tsplit\t780978\t780971\tcanu:missassembly\n"
  for Assembly in $(ls assembly/canu/*/*/*_canu.contigs.fasta | grep -w 'Fus2'); do
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)  
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    OutDir=assembly/canu/$Organism/$Strain/filtered_contigs
    $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/"$Strain"_canu_contigs_renamed.fasta --coord_file tmp.csv
  done
  rm tmp.csv
```

After investigation it was found that contig_17 should be split.

```bash
  ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
  touch tmp.csv
  printf "contig_17\tmanual edit\tsplit\t780978\t780978\tcanu:missassembly\n" > tmp.csv
  for Assembly in $(ls assembly/canu/F.oxysporum_fsp_cepae/Fus2/filtered_contigs/Fus2_canu_contigs_renamed.fasta); do
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    OutDir=assembly/canu/$Organism/$Strain/edited_contigs
    mkdir -p $OutDir
    $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/"$Strain"_canu_contigs_modified.fasta --coord_file tmp.csv
  done
  rm tmp.csv
```

Assemblies were polished using Pilon

```bash
  for Assembly in $(ls assembly/canu/*/Fus2/edited_contigs/*_canu_contigs_modified.fasta); do
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    IlluminaDir=$(ls -d qc_dna/paired/$Organism/$Strain)
    TrimF1_Read=$(ls $IlluminaDir/F/s_6_1_sequence_trim.fq.gz);
    TrimR1_Read=$(ls $IlluminaDir/R/s_6_2_sequence_trim.fq.gz);
    OutDir=assembly/canu/$Organism/$Strain/edited_contigs
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/pilon
    qsub $ProgDir/sub_pilon.sh $Assembly $TrimF1_Read $TrimR1_Read $OutDir
  done
```


### Spades Assembly

```bash
  for PacBioDat in $(ls raw_dna/pacbio/F.oxysporum_fsp_cepae/Fus2/extracted/concatenated_pacbio.fastq); do
    Organism=$(echo $PacBioDat | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $PacBioDat | rev | cut -f3 -d '/' | rev)
    IlluminaDir=$(ls -d qc_dna/paired/$Organism/$Strain)
    TrimF1_Read=$(ls $IlluminaDir/F/s_6_1_sequence_trim.fq.gz);
    TrimR1_Read=$(ls $IlluminaDir/R/s_6_2_sequence_trim.fq.gz);
    TrimF2_Read=$(ls $IlluminaDir/F/FUS2_S2_L001_R1_001_trim.fq.gz);
    TrimR2_Read=$(ls $IlluminaDir/R/FUS2_S2_L001_R2_001_trim.fq.gz);
    OutDir=assembly/spades_pacbio/$Organism/"$Strain"
    echo $TrimR1_Read
    echo $TrimR1_Read
    echo $TrimF2_Read
    echo $TrimR2_Read
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades/multiple_libraries
    qsub $ProgDir/subSpades_2lib_pacbio.sh $PacBioDat $TrimF1_Read $TrimR1_Read $TrimF2_Read $TrimR2_Read $OutDir 20
  done
```

Contigs shorter thaan 500bp were removed from the assembly

```bash
  for Contigs in $(ls assembly/spades_pacbio/*/*/contigs.fasta); do
    AssemblyDir=$(dirname $Contigs)
    mkdir $AssemblyDir/filtered_contigs
    FilterDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/abyss
    $FilterDir/filter_abyss_contigs.py $Contigs 500 > $AssemblyDir/filtered_contigs/contigs_min_500bp.fasta
  done
```


## Merging pacbio and hybrid assemblies

```bash
  # for PacBioAssembly in $(ls assembly/canu/*/*/polished/pilon.fasta); do
  for PacBioAssembly in $(ls assembly/pacbio_test/*/*/edited_contigs/pilon.fasta); do
  # for PacBioAssembly in $(ls assembly/pacbio_test/F.oxysporum_fsp_cepae/Fus2_pacbio_merged/polished/pilon.fasta); do
    Organism=$(echo $PacBioAssembly | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $PacBioAssembly | rev | cut -f3 -d '/' | rev)
    HybridAssembly=$(ls assembly/spades_pacbio/$Organism/$Strain/contigs.fasta)
    QuastReport=$(ls assembly/canu/$Organism/$Strain/filtered_contigs/report.tsv)
    # N50=$(cat $QuastReport | grep 'N50' | cut -f2)
    # AnchorLength=$N50
    AnchorLength=20000
    # OutDir=assembly/merged_canu_spades/$Organism/$Strain
    OutDir=assembly/merged_canu_spades/$Organism/"$Strain"_edited2
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/quickmerge
    qsub $ProgDir/sub_quickmerge.sh $PacBioAssembly $HybridAssembly $OutDir $AnchorLength
  done
```

```bash
  touch tmp.csv
  for Assembly in $(ls assembly/merged_canu_spades/F.oxysporum_fsp_cepae/Fus2_edited2/merged.fasta); do
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)  
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    OutDir=$(dirname $Assembly)
    mkdir -p $OutDir
    ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
    $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/"$Strain"_contigs_renamed.fasta --coord_file tmp.csv
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
  rm tmp.csv
```

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
for Assembly in $(ls assembly/merged_canu_spades/F.oxysporum_fsp_cepae/Fus2_edited2/Fus2_edited2_contigs_renamed.fasta); do
echo $Assembly
Query=analysis/blast_homology/Fo_path_genes/Fo_path_genes_CRX.fa
qsub $ProgDir/blast_pipe.sh $Query dna $Assembly
done
```


This merged assembly was polished using Pilon

```bash
  # for Assembly in $(ls assembly/merged_canu_spades/F.oxysporum_fsp_cepae/Fus2/merged.fasta); do
  for Assembly in $(ls assembly/merged_canu_spades/F.oxysporum_fsp_cepae/Fus2_edited2/merged.fasta); do
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    # IlluminaDir=$(ls -d qc_dna/paired/$Organism/$Strain)
    IlluminaDir=$(ls -d qc_dna/paired/$Organism/Fus2)
    TrimF1_Read=$(ls $IlluminaDir/F/s_6_1_sequence_trim.fq.gz);
    TrimR1_Read=$(ls $IlluminaDir/R/s_6_2_sequence_trim.fq.gz);
    OutDir=assembly/merged_canu_spades/$Organism/$Strain/polished
    # OutDir=assembly/pacbio_test/$Organism/$Strain/polished
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/pilon
    qsub $ProgDir/sub_pilon.sh $Assembly $TrimF1_Read $TrimR1_Read $OutDir
  done
```


Contigs were renamed in accordance with ncbi recomendations.

```bash
  ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
  touch tmp.csv
  # for Assembly in $(ls assembly/merged_canu_spades/*/*/polished/pilon.fasta); do
  # for Assembly in $(ls assembly/pacbio_test/F.oxysporum_fsp_cepae/Fus2_pacbio_merged/polished/pilon.fasta); do
    for Assembly in $(ls assembly/merged_canu_spades/F.oxysporum_fsp_cepae/Fus2_edited2/polished/pilon.fasta); do
  # for Assembly in $(ls assembly/merged_canu_spades/F.oxysporum_fsp_cepae/Fus2_edited2/polished/pilon.fasta); do
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    # OutDir=assembly/merged_canu_spades/$Organism/$Strain/filtered_contigs
    OutDir=assembly/merged_canu_spades/$Organism/$Strain/filtered_contigs
    mkdir -p $OutDir
    $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/"$Strain"_contigs_renamed.fasta --coord_file tmp.csv
  done
  rm tmp.csv
```

Assembly stats were collected using quast

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  # for Assembly in $(ls assembly/merged_canu_spades/*/*/filtered_contigs/Fus2_contigs_renamed.fasta); do
  # for Assembly in $(ls assembly/pacbio_test/F.oxysporum_fsp_cepae/Fus2_pacbio_merged/filtered_contigs/Fus2_pacbio_merged_contigs_renamed.fasta); do
  for Assembly in $(ls assembly/merged_canu_spades/F.oxysporum_fsp_cepae/Fus2_edited2/Fus2_edited2_contigs_renamed.fasta); do
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)  
    OutDir=$(dirname $Assembly)
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```

## Renaming assemblies - temporarily
Fus2 was temporarily renamed for preliminary analysis

```bash
  cp -r assembly/canu/F.oxysporum_fsp_cepae/Fus2 assembly/canu/F.oxysporum_fsp_cepae/Fus2_pacbio_test_canu
  cp -r assembly/merged_canu_spades/F.oxysporum_fsp_cepae/Fus2 assembly/merged_canu_spades/F.oxysporum_fsp_cepae/Fus2_pacbio_test_merged
  cp -r assembly/pacbio_test/F.oxysporum_fsp_cepae/Fus2_pacbio_merged assembly/pacbio_test/F.oxysporum_fsp_cepae/Fus2_pacbio_merged_richards
```

# Repeatmasking assemblies

```bash
  # Fus2_pacbio_canu=$(ls assembly/canu/F.oxysporum_fsp_cepae/Fus2_pacbio_test_canu/filtered_contigs/Fus2_canu_contigs_renamed.fasta)
  Fus2_pacbio_merged=$(ls assembly/merged_canu_spades/*/Fus2/filtered_contigs/Fus2_contigs_renamed.fasta)
  # for Assembly in $(ls $Fus2_pacbio_merged $Fus2_pacbio_canu); do
  for Assembly in $(ls $Fus2_pacbio_merged); do
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking
    qsub $ProgDir/rep_modeling.sh $Assembly
    qsub $ProgDir/transposonPSI.sh $Assembly
  done
```

# Preliminary analysis

## Checking PacBio coverage against Fus2 contigs

The accuracy of PacBio assembly pipelines is currently unknown. To help identify
regions that may have been missassembled the pacbio reads were aligned back to
the assembled genome. Coverage was determined using bedtools genomecov and
regions with low coverage flagged using a python script flag_low_coverage.py.
These low coverage regions were visually inspected using IGV.

```bash
  Assembly=assembly/merged_canu_spades/F.oxysporum_fsp_cepae/Fus2/filtered_contigs/Fus2_contigs_renamed.fasta
  Reads=raw_dna/pacbio/F.oxysporum_fsp_cepae/Fus2/extracted/concatenated_pacbio.fastq
  OutDir=analysis/genome_alignment/bwa/F.oxysporum_fsp_cepae/Fus2/vs_Fus2
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/bwa
  qsub $ProgDir/sub_bwa_pacbio.sh $Assembly $Reads $OutDir

  Assembly=assembly/merged_canu_spades/F.oxysporum_fsp_cepae/Fus2_edited/filtered_contigs/Fus2_edited_contigs_renamed.fasta
  Reads=raw_dna/pacbio/F.oxysporum_fsp_cepae/Fus2/extracted/concatenated_pacbio.fastq
  OutDir=analysis/genome_alignment/bwa/F.oxysporum_fsp_cepae/Fus2/vs_Fus2_edited
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/bwa
  qsub $ProgDir/sub_bwa_pacbio.sh $Assembly $Reads $OutDir

  Assembly=assembly/canu/F.oxysporum_fsp_cepae/Fus2_edited3/Fus2_canu_manual_edits.fasta
  Reads=raw_dna/pacbio/F.oxysporum_fsp_cepae/Fus2/extracted/concatenated_pacbio.fastq
  OutDir=analysis/genome_alignment/bwa/F.oxysporum_fsp_cepae/Fus2/vs_Fus2_edited3
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/bwa
  qsub $ProgDir/sub_bwa_pacbio.sh $Assembly $Reads $OutDir
```
<!--
The same analysis was performed on the pacbio only assembly to see if errors
occurred at the merging step:

```bash
  Assembly=assembly/canu/F.oxysporum_fsp_cepae/Fus2/filtered_contigs/Fus2_canu_contigs_renamed.fasta
  Reads=raw_dna/pacbio/F.oxysporum_fsp_cepae/Fus2/extracted/concatenated_pacbio.fastq
  OutDir=analysis/genome_alignment/bwa/F.oxysporum_fsp_cepae/Fus2/vs_Fus2_canu_only
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/bwa
  qsub $ProgDir/sub_bwa_pacbio.sh $Assembly $Reads $OutDir

  AlignedBam=$OutDir/Fus2_canu_contigs_renamed.fasta_aligned_sorted.bam
  CoverageTxt=$OutDir/Fus2_bp_genome_cov.txt
  bedtools genomecov -max 5 -d -ibam $AlignedBam -g $Assembly > $CoverageTxt

  Threshold=5
  FlaggedRegions=$OutDir/Fus2_flagged_regions.txt
  $ProgDir/flag_low_coverage.py --genomecov $CoverageTxt --min $Threshold > $FlaggedRegions

```
-->


## Blast searches of LS region genes vs FoC

Some preliminary commands were used to analyse Pacbio assemblies Richard had generated

Blast searches were performed against these assembled contigs to identify which
contigs contained blast homologs from known FoL LS genes.

Headers of LS genes from FoL had previously been extracted. These were used to
extract the relevant proteins from fasta files.

```bash
  ProtFastaUnparsed=assembly/external_group/F.oxysporum_fsp_lycopersici/4287/Fusox1/Fusox1_GeneCatalog_proteins_20110522.aa.fasta
  ProtFasta=analysis/FoL_ls_genes/4287_proteins_renamed.fasta
  cat $ProtFastaUnparsed | sed -r "s/^>.*FOXG/>FOXG/g" > $ProtFasta
  for File in $(ls analysis/FoL_ls_genes/chr_*_gene_headers.txt); do
    Chr=$(echo $File | rev | cut -f3 -d '_' | rev);
    echo $File;
    echo "extracting proteins associated with chromosome: $Chr";
    OutFile=$(echo $File | sed 's/_headers.txt/.aa/g')
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_from_fasta.py --fasta $ProtFasta --headers $File | grep -v -P '^$' > $OutFile
  done
```

Searches were not performed for these genes as it was found that many are in
multiple copy gene families and therefore have many blast hits throughout the
genome. Single copy genes in the genome were used instead.


<!--
## Merging pacbio and hybrid assemblies

```bash

for PacBioAssembly in $(ls assembly/pacbio_test/F.oxysporum_fsp_cepae/Fus2_pacbio_merged/test/Fus2_pacbio_merged.fa); do
Organism=$(echo $PacBioAssembly | rev | cut -f4 -d '/' | rev)
Strain=$(echo $PacBioAssembly | rev | cut -f3 -d '/' | rev)
HybridAssembly=$(ls assembly/spades_pacbio/$Organism/Fus2/contigs.fasta)
OutDir=assembly/pacbio_test/$Organism/$Strain
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/quickmerge
qsub $ProgDir/sub_quickmerge.sh $PacBioAssembly $HybridAssembly $OutDir
done
``` -->

<!-- ```bash
  Fus2_pacbio_canu=$(ls assembly/canu/F.oxysporum_fsp_cepae/Fus2_pacbio_test/filtered_contigs/Fus2_canu_contigs_renamed.fasta)
  Fus2_pacbio_merged=$(ls assembly/canu_spades_hybrid/F.oxysporum_fsp_cepae/Fus2_pacbio_test/filtered_contigs/Fus2_contigs_renamed.fasta)
  for $Genome in $(ls $Fus2_pacbio_merged $Fus2_pacbio_canu); do
    for Proteome in $(ls analysis/FoL_ls_genes/chr_*_gene.aa); do
      Organism=$(echo $Genome | rev | cut -f4 -d '/' | rev)
      Strain=$(echo $Genome | rev | cut -f3 -d '/' | rev)
      Chr=$(echo $Proteome | rev | cut -f2 -d '_' | rev);
      echo "$Organism - $Strain - $Chr"
      OutDir=analysis/blast_homology/$Organism/$Strain
      ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
      qsub $ProgDir/run_blast2csv.sh $Proteome protein $Genome $OutDir
    done
  done
```

Convert top blast hits into gff annotations

```bash
  for BlastHitsCsv in $(ls analysis/blast_homology/*/Fus2_pacbio_test_*/*_chr_*_gene.aa_hits.csv); do
    Organism=$(echo $BlastHitsCsv | rev | cut -f3 -d '/' | rev)
    Strain=$(echo $BlastHitsCsv | rev | cut -f2 -d '/' | rev)
    echo "$Organism - $Strain"
    HitsGff=$(echo $BlastHitsCsv | sed  's/.csv/.gff/g')
    Chr=$(echo $BlastHitsCsv | rev | cut -f3 -d '_' | rev);
    Column2=Chr"$Chr"_gene_homolog
    NumHits=1
    ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
    $ProgDir/blast2gff.pl $Column2 $NumHits $BlastHitsCsv > $HitsGff
  done
```

-->



## Blast searches of LS region genes in single copy orthogroups vs FoC


Blast searches were performed against these assembled contigs to identify which
contigs contained blast homologs from known FoL LS genes.

Headers of LS genes in single copy orthogroups that were also present in a
single copy in FoC were extracted. These headers were used to extract the
relevant proteins from fasta files.

```bash
mkdir -p analysis/FoL_genes
FoLTable=gene_pred/annotations/F.oxysporum_fsp_lycopersici/4287/4287_gene_annotations.tab
for num in $(seq 1 15); do
# echo "$num"
TotalGenes=$(cat $FoLTable | cut -f1,2,17,18 |  grep -P "\t$num\t" | grep -v '4287(0)' | grep -v 'Fus2(0)' | wc -l)
cat $FoLTable | cut -f1,2,17,18 |  grep -P "\t$num\t" | grep '4287(1)' | grep 'Fus2(1)' | cut -f1 > analysis/FoL_genes/chr_"$num"_gene_single_copy_headers.txt
# ls analysis/FoL_genes/chr_"$num"_gene_single_copy_headers.txt
LsGenes=$(cat analysis/FoL_genes/chr_"$num"_gene_single_copy_headers.txt | wc -l)
echo "Chr$num - $LsGenes ($TotalGenes)"
done
```
This left the following number of genes for each chromosome to act as markers:
Chr1 - 1278 (2022)
Chr2 - 983 (1622)
Chr3 - 11 (748)
Chr4 - 1101 (1761)
Chr5 - 975 (1540)
Chr6 - 14 (602)
Chr7 - 902 (1429)
Chr8 - 788 (1249)
Chr9 - 663 (1083)
Chr10 - 550 (954)
Chr11 - 406 (878)
Chr12 - 410 (800)
Chr13 - 326 (644)
Chr14 - 15 (219)
Chr15 - 0 (343)

Note - There were no single copy orthologs from genes on Chromosome 15.

```bash
  ProtFastaUnparsed=assembly/external_group/F.oxysporum_fsp_lycopersici/4287/Fusox1/Fusox1_GeneCatalog_proteins_20110522.aa.fasta
  ProtFasta=analysis/FoL_ls_genes/4287_proteins_renamed.fasta
  cat $ProtFastaUnparsed | sed -r "s/^>.*FOXG/>FOXG/g" > $ProtFasta
  for File in $(ls analysis/FoL_genes/chr_*_gene_single_copy_headers.txt); do
    Chr=$(echo $File | rev | cut -f5 -d '_' | rev);
    echo $File;
    echo "extracting proteins associated with chromosome: $Chr";
    OutFile=$(echo $File | sed 's/_headers.txt/.aa/g')
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_from_fasta.py --fasta $ProtFasta --headers $File | grep -v -P '^$' > $OutFile
  done
```

```bash
  Fus2_pacbio_merged=$(ls assembly/canu/F.oxysporum_fsp_cepae/Fus2_edited3/Fus2_canu_manual_edits.fasta)
  # Fus2_pacbio_merged=$(ls assembly/merged_canu_spades/F.oxysporum_fsp_cepae/Fus2_edited2/Fus2_edited2_contigs_renamed.fasta)
  for Genome in $(ls $Fus2_pacbio_merged); do
    for Proteome in $(ls analysis/FoL_genes/chr_*_gene_single_copy.aa); do
      Organism=$(echo $Genome | rev | cut -f4 -d '/' | rev)
      Strain=$(echo $Genome | rev | cut -f3 -d '/' | rev)
      Chr=$(echo $Proteome | rev | cut -f4 -d '_' | rev);
      echo "$Organism - $Strain - $Chr"
      OutDir=analysis/blast_homology/$Organism/$Strain
      ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
      qsub $ProgDir/run_blast2csv.sh $Proteome protein $Genome $OutDir
    done
  done
```

Convert top blast hits into gff annotations

```bash
  for BlastHitsCsv in $(ls analysis/blast_homology/canu/F.oxysporum_fsp_cepae/*_chr_*_gene_single_copy.aa_hits.csv); do
  # for BlastHitsCsv in $(ls analysis/blast_homology/*/Fus2_edited3/*_chr_*_gene_single_copy.aa_hits.csv); do
    Organism=$(echo $BlastHitsCsv | rev | cut -f3 -d '/' | rev)
    Strain=$(echo $BlastHitsCsv | rev | cut -f2 -d '/' | rev)
    echo "$Organism - $Strain"
    HitsGff=$(echo $BlastHitsCsv | sed  's/.csv/.gff/g')
    Chr=$(echo $BlastHitsCsv | rev | cut -f5 -d '_' | rev);
    Column2=Chr"$Chr"_gene_homolog
    NumHits=1
    ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
    $ProgDir/blast2gff.pl $Column2 $NumHits $BlastHitsCsv > $HitsGff
  done
```

Copying data

```bash
  cp repeat_masked/F.oxysporum_fsp_cepae/Fus2/filtered_contigs_repmask/Fus2_contigs_unmasked.fa pacbio_assembly/.
  cp repeat_masked/F.oxysporum_fsp_cepae/Fus2/filtered_contigs_repmask/Fus2_contigs_hardmasked.gff pacbio_assembly/.
  cp repeat_masked/F.oxysporum_fsp_cepae/Fus2/filtered_contigs_repmask/Fus2_contigs_unmasked.fa.TPSI.allHits.chains.gff pacbio_assembly/.
  cp repeat_masked/F.oxysporum_fsp_cepae/Fus2/filtered_contigs_repmask/Fus2_contigs_unmasked.fa.TPSI.allHits.chains.gff3 pacbio_assembly/.
  cp gene_pred/codingquary/F.oxysporum_fsp_cepae/Fus2/final/final_genes_combined.* pacbio_assembly/.
  cp gene_pred/codingquary/F.oxysporum_fsp_cepae/Fus2/final/final_genes_appended.gff3 pacbio_assembly/.
  cp gene_pred/interproscan/F.oxysporum_fsp_cepae/Fus2/Fus2_interproscan.tsv pacbio_assembly/.
  cp analysis/mimps/F.oxysporum_fsp_cepae/Fus2/Fus2_mimps.gff pacbio_assembly/.
  cp analysis/mimps/F.oxysporum_fsp_cepae/Fus2/Fus2_genes_in_2kb_mimp.gff pacbio_assembly/.
  cp analysis/blast_homology/F.oxysporum_fsp_cepae/Fus2/Fus2_Fo_path_genes_CRX.fa_homologs.csv pacbio_assembly/.
  cp analysis/blast_homology/F.oxysporum_fsp_cepae/Fus2/Fus2_Fo_path_genes_CRX.fa_homologs.gff pacbio_assembly/.
  mkdir pacbio_assembly/orthology
  cp analysis/orthology/orthomcl/FoC_vs_Fo_vs_FoL/FoC_vs_Fo_vs_FoL_orthogroups.txt pacbio_assembly/orthology/.
  cp -r analysis/orthology/orthomcl/FoC_vs_Fo_vs_FoL/Fus2_genes pacbio_assembly/orthology/.
  cp -r analysis/orthology/orthomcl/FoC_vs_Fo_vs_FoL/fasta pacbio_assembly/orthology/.
```





<!--
# Analysis of preliminary assemblies

Analysis was performed on the genomes generated by Richard.

## Preparing data

projects/pacbio_test/fus2-auto/fus2.contigs.fasta
```bash
  ProjDir=/home/groups/harrisonlab/project_files/fusarium
  cd $ProjDir

  Fus2_pacbio_canu=assembly/pacbio_test/F.oxysporum_fsp_cepae/Fus2_pacbio_canu/test/Fus2_pacbio_canu.fa
  Fus2_pacbio_spades=assembly/pacbio_test/F.oxysporum_fsp_cepae/Fus2_pacbio_spades/test/Fus2_pacbio_spades.fa
  Fus2_pacbio_merged=assembly/pacbio_test/F.oxysporum_fsp_cepae/Fus2_pacbio_merged/test/Fus2_pacbio_merged.fa

  mkdir -p $(dirname $Fus2_pacbio_canu)
  mkdir -p $(dirname $Fus2_pacbio_spades)
  mkdir -p $(dirname $Fus2_pacbio_merged)
  cp /home/harrir/projects/pacbio_test/fus2/fus2-auto/fus2.contigs.fasta $Fus2_pacbio_canu
  cp /home/harrir/projects/pacbio_test/spades/FUS2/scaffolds.fasta $Fus2_pacbio_spades
  cp /home/harrir/projects/pacbio_test/hybrid_merge/fus2/merged.fasta $Fus2_pacbio_merged
```

## Repeatmasking

```bash
  Fus2_pacbio_canu=assembly/pacbio_test/F.oxysporum_fsp_cepae/Fus2_pacbio_canu/test/Fus2_pacbio_canu.fa
  Fus2_pacbio_spades=assembly/pacbio_test/F.oxysporum_fsp_cepae/Fus2_pacbio_spades/test/Fus2_pacbio_spades.fa
  Fus2_pacbio_merged=assembly/pacbio_test/F.oxysporum_fsp_cepae/Fus2_pacbio_merged/test/Fus2_pacbio_merged.fa
  for Assembly in $(ls $Fus2_pacbio_spades $Fus2_pacbio_merged $Fus2_pacbio_canu); do
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking
    qsub $ProgDir/rep_modeling.sh $Assembly
    qsub $ProgDir/transposonPSI.sh $Assembly
  done
```

## Blast searches of LS region genes vs FoC

Some preliminary commands were used to analyse Pacbio assemblies Richard had generated

Blast searches were performed against these assembled contigs to identify which
contigs contained blast homologs from known FoL LS genes.

Headers of LS genes from FoL had previously been extracted. These were used to
extract the relevant proteins from fasta files.

```bash
  ProtFastaUnparsed=assembly/external_group/F.oxysporum_fsp_lycopersici/4287/Fusox1/Fusox1_GeneCatalog_proteins_20110522.aa.fasta
  ProtFasta=analysis/FoL_ls_genes/4287_proteins_renamed.fasta
  cat $ProtFastaUnparsed | sed -r "s/^>.*FOXG/>FOXG/g" > $ProtFasta
  for File in $(ls analysis/FoL_ls_genes/chr_*_gene_headers.txt); do
    Chr=$(echo $File | rev | cut -f3 -d '_' | rev);
    echo $File;
    echo "extracting proteins associated with chromosome: $Chr";
    OutFile=$(echo $File | sed 's/_headers.txt/.aa/g')
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_from_fasta.py --fasta $ProtFasta --headers $File | grep -v -P '^$' > $OutFile
  done
```

```bash
  Fus2_pacbio_canu=assembly/pacbio_test/F.oxysporum_fsp_cepae/Fus2_pacbio_canu/test/Fus2_pacbio_canu.fa
  Fus2_pacbio_spades=assembly/pacbio_test/F.oxysporum_fsp_cepae/Fus2_pacbio_spades/test/Fus2_pacbio_spades.fa
  Fus2_pacbio_merged=assembly/pacbio_test/F.oxysporum_fsp_cepae/Fus2_pacbio_merged/test/Fus2_pacbio_merged.fa
  for Genome in $(ls $Fus2_pacbio_spades $Fus2_pacbio_merged $Fus2_pacbio_canu); do
    for Proteome in $(ls analysis/FoL_ls_genes/chr_*_gene.aa); do
      Strain=$(echo $Genome | rev | cut -f3 -d '/' | rev)
      Organism=$(echo $Genome | rev | cut -f4 -d '/' | rev)
      Chr=$(echo $Proteome | rev | cut -f2 -d '_' | rev);
      echo "$Organism - $Strain - $Chr"
      OutDir=analysis/blast_homology/$Organism/$Strain
      ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
      qsub $ProgDir/run_blast2csv.sh $Proteome protein $Genome $OutDir
    done
  done
```

Convert top blast hits into gff annotations

```bash
  for BlastHitsCsv in $(ls analysis/blast_homology/*/Fus2_pacbio_*/*_chr_*_gene.aa_hits.csv); do
    Organism=$(echo $BlastHitsCsv | rev | cut -f3 -d '/' | rev)
    Strain=$(echo $BlastHitsCsv | rev | cut -f2 -d '/' | rev)
    echo "$Organism - $Strain"
    HitsGff=$(echo $BlastHitsCsv | sed  's/.csv/.gff/g')
    Chr=$(echo $BlastHitsCsv | rev | cut -f3 -d '_' | rev);
    Column2=Chr"$Chr"_gene_homolog
    NumHits=1
    ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
    $ProgDir/blast2gff.pl $Column2 $NumHits $BlastHitsCsv > $HitsGff
  done
```

## Blast searches of LS region genes in single copy orthogroups vs FoC


Blast searches were performed against these assembled contigs to identify which
contigs contained blast homologs from known FoL LS genes.

Headers of LS genes in single copy orthogroups that were also present in a
single copy in FoC were extracted. These headers were used to extract the
relevant proteins from fasta files.

```bash
FoLTable=gene_pred/annotations/F.oxysporum_fsp_lycopersici/4287/4287_gene_annotations.tab
cat $FoLTable | cut -f1,2,17,18 |  grep -P '\t3\t' | grep '4287(1)' | grep 'Fus2(1)' | cut -f1 > analysis/FoL_ls_genes/chr_3_gene_single_copy_headers.txt
cat $FoLTable | cut -f1,2,17,18 |  grep -P '\t6\t' | grep '4287(1)' | grep 'Fus2(1)' | cut -f1 > analysis/FoL_ls_genes/chr_6_gene_single_copy_headers.txt
cat $FoLTable | cut -f1,2,17,18 |  grep -P '\t14\t' | grep '4287(1)' | grep 'Fus2(1)' | cut -f1 > analysis/FoL_ls_genes/chr_14_gene_single_copy_headers.txt
cat $FoLTable | cut -f1,2,17,18 |  grep -P '\t15\t' | grep '4287(1)' | grep 'Fus2(1)' | cut -f1 > analysis/FoL_ls_genes/chr_15_gene_single_copy_headers.txt
cat $FoLTable | cut -f1,2,17,18 | sed -r "s/\t/,/g" | grep -e ',3,' -e ',6,' -e ',14,' -e ',15,' | sed -r "s/,/\t/g" | grep '4287(1)' | grep 'Fus2(1)' | cut -f2 | sort | uniq -c | sort -n -r
cat $FoLTable | cut -f1,2,17,18 | sed -r "s/\t/,/g" | grep -v -e ',3,' -e ',6,' -e ',14,' -e ',15,' | sed -r "s/,/\t/g" | grep '4287(1)' | grep 'Fus2(1)' | cut -f1 > analysis/FoL_ls_genes/chr_nonpath_gene_single_copy_headers.txt
```
This left the following number of genes for each chromosome to act as markers:
20 14
18 6
 9 3

Note - There were no single copy orthologs from genes on Chromosome 15.

```bash
  ProtFastaUnparsed=assembly/external_group/F.oxysporum_fsp_lycopersici/4287/Fusox1/Fusox1_GeneCatalog_proteins_20110522.aa.fasta
  ProtFasta=analysis/FoL_ls_genes/4287_proteins_renamed.fasta
  cat $ProtFastaUnparsed | sed -r "s/^>.*FOXG/>FOXG/g" > $ProtFasta
  for File in $(ls analysis/FoL_ls_genes/chr_*_gene_single_copy_headers.txt); do
    Chr=$(echo $File | rev | cut -f5 -d '_' | rev);
    echo $File;
    echo "extracting proteins associated with chromosome: $Chr";
    OutFile=$(echo $File | sed 's/_headers.txt/.aa/g')
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_from_fasta.py --fasta $ProtFasta --headers $File | grep -v -P '^$' > $OutFile
  done
```

```bash
  Fus2_pacbio_canu=assembly/pacbio_test/F.oxysporum_fsp_cepae/Fus2_pacbio_canu/test/Fus2_pacbio_canu.fa
  Fus2_pacbio_spades=assembly/pacbio_test/F.oxysporum_fsp_cepae/Fus2_pacbio_spades/test/Fus2_pacbio_spades.fa
  Fus2_pacbio_merged=assembly/pacbio_test/F.oxysporum_fsp_cepae/Fus2_pacbio_merged/test/Fus2_pacbio_merged.fa
  for Genome in $(ls $Fus2_pacbio_spades $Fus2_pacbio_merged $Fus2_pacbio_canu); do
    for Proteome in $(ls analysis/FoL_ls_genes/chr_*_gene_single_copy.aa); do
      Strain=$(echo $Genome | rev | cut -f3 -d '/' | rev)
      Organism=$(echo $Genome | rev | cut -f4 -d '/' | rev)
      Chr=$(echo $Proteome | rev | cut -f4 -d '_' | rev);
      echo "$Organism - $Strain - $Chr"
      OutDir=analysis/blast_homology/$Organism/$Strain
      ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
      qsub $ProgDir/run_blast2csv.sh $Proteome protein $Genome $OutDir
    done
  done
```

Convert top blast hits into gff annotations

```bash
  for BlastHitsCsv in $(ls analysis/blast_homology/*/Fus2_pacbio_*/*_chr_*_gene_single_copy.aa_hits.csv); do
    Organism=$(echo $BlastHitsCsv | rev | cut -f3 -d '/' | rev)
    Strain=$(echo $BlastHitsCsv | rev | cut -f2 -d '/' | rev)
    echo "$Organism - $Strain"
    HitsGff=$(echo $BlastHitsCsv | sed  's/.csv/.gff/g')
    Chr=$(echo $BlastHitsCsv | rev | cut -f5 -d '_' | rev);
    Column2=Chr"$Chr"_gene_homolog
    NumHits=1
    ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
    $ProgDir/blast2gff.pl $Column2 $NumHits $BlastHitsCsv > $HitsGff
  done
``` -->
