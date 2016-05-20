
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
Reads=$(ls raw_dna/pacbio/F.oxysporum_fsp_cepae/Fus2/extracted/Fus2_pacbio.fastq)
GenomeSz="50m"
Strain=$(echo $Reads | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Reads | rev | cut -f4 -d '/' | rev)
Prefix="$Strain"_canu
OutDir="assembly/canu/$Organism/$Strain"
ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/canu
qsub $ProgDir/submit_canu.sh $Reads $GenomeSz $Prefix $OutDir
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
qsub $ProgDir/subSpades_2lib_pacbio.sh $PacBioDat $TrimF1_Read $TrimR1_Read $TrimF2_Read $TrimR2_Read $OutDir 50
done
```





# Analysis of preliminary assemblies

## Preparing data


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
```
