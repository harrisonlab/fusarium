
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
  for Assembly in $(ls assembly/canu-1.3/*/*/*_canu.contigs.fasta); do
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)  
    OutDir=assembly/canu-1.3/$Organism/$Strain/filtered_contigs
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```


Contigs were renamed in accordance with ncbi recomendations

```bash
  ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
  touch tmp.csv
  for Assembly in $(ls assembly/canu-1.3/*/*/*_canu.contigs.fasta); do
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)  
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    OutDir=assembly/canu/$Organism/$Strain/filtered_contigs
    $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/"$Strain"_canu_contigs_renamed.fasta --coord_file tmp.csv
  done
  rm tmp.csv
```


Assemblies were polished using Pilon

```bash
  for Assembly in $(ls assembly/canu-1.3/F.oxysporum_fsp_cepae/Fus2_canu/filtered_contigs/Fus2_canu_contigs_renamed.fasta); do
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    IlluminaDir=$(ls -d qc_dna/paired/$Organism/Fus2)
    TrimF1_Read=$(ls $IlluminaDir/F/s_6_1_sequence_trim.fq.gz);
    TrimR1_Read=$(ls $IlluminaDir/R/s_6_2_sequence_trim.fq.gz);
    OutDir=assembly/canu-1.3/$Organism/$Strain/pilon
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/pilon
    qsub $ProgDir/sub_pilon.sh $Assembly $TrimF1_Read $TrimR1_Read $OutDir
  done
```

```bash
  Assembly=assembly/canu-1.3/F.oxysporum_fsp_cepae/Fus2_canu/pilon/pilon.fasta
  Reads=raw_dna/pacbio/F.oxysporum_fsp_cepae/Fus2/extracted/concatenated_pacbio.fastq
  OutDir=analysis/genome_alignment/bwa/F.oxysporum_fsp_cepae/Fus2/vs_Fus2
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/bwa
  qsub $ProgDir/sub_bwa_pacbio.sh $Assembly $Reads $OutDir
```

After investigation, it was found that contigs didnt need to be split.


A Bioproject and Biosample was made with NCBI genbank for submission of genomes.
Following the creation of these submissions, the .fasta assembly was uploaded
through the submission portal. A note was provided requesting that the assembly
be run through the contamination screen to aid a more detailed resubmission in
future. The returned FCSreport.txt was downloaded from the NCBI webportal and
used to correct the assembly to NCBI standards.

NCBI reports (FCSreport.txt) were manually downloaded to the following loactions:

```bash
	for Assembly in $(ls assembly/canu-1.3/F.oxysporum_fsp_cepae/Fus2_canu/pilon/pilon.fasta); do
		Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
		Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
		NCBI_report_dir=genome_submission/$Organism/$Strain/initial_submission
		mkdir -p $NCBI_report_dir
	done
```

These downloaded files were used to correct assemblies:

```bash
	for Assembly in $(ls assembly/canu-1.3/F.oxysporum_fsp_cepae/Fus2_canu/pilon/pilon.fasta); do
		Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
		Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
		NCBI_report=$(ls genome_submission/$Organism/$Strain/initial_submission/FCSreport.txt)
		OutDir=assembly/canu-1.3/$Organism/$Strain/ncbi_edits
		mkdir -p $OutDir
		ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
		$ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/contigs_min_500bp_renamed.fasta --coord_file $NCBI_report > $OutDir/log.txt
	done
```


Assembly stats were collected using quast

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly/canu-1.3/F.oxysporum_fsp_cepae/Fus2_canu/pilon/pilon.fasta); do
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)  
    OutDir=assembly/canu-1.3/$Organism/$Strain/pilon
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```

<!--
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
 -->

 Preliminary analysis of these contigs allowed the quality of this assembly to be assessed:

 ```bash
   ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
   for Assembly in $(ls assembly/canu-1.3/F.oxysporum_fsp_cepae/Fus2_canu/pilon/pilon.fasta); do
     echo $Assembly
     Query=analysis/blast_homology/Fo_path_genes/Fo_path_genes_CRX.fa
     qsub $ProgDir/blast_pipe.sh $Query dna $Assembly
   done
 ```

 ```bash
   for Genome in $(ls assembly/canu-1.3/F.oxysporum_fsp_cepae/Fus2_canu/pilon/pilon.fasta); do
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
   for BlastHitsCsv in $(ls analysis/blast_homology/F.oxysporum_fsp_cepae/Fus2_canu/*_chr_*_gene_single_copy.aa_hits.csv); do
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
   for BlastHitsCsv in $(ls analysis/blast_homology/F.oxysporum_fsp_cepae/Fus2_canu/Fus2_canu_Fo_path_genes_CRX.fa_homologs.csv); do
     Organism=$(echo $BlastHitsCsv | rev | cut -f3 -d '/' | rev)
     Strain=$(echo $BlastHitsCsv | rev | cut -f2 -d '/' | rev)
     echo "$Organism - $Strain"
     HitsGff=$(echo $BlastHitsCsv | sed  's/.csv/.gff/g')
     Chr=$(echo $BlastHitsCsv | rev | cut -f5 -d '_' | rev);
     Column2=Fo_path_gene_homolog
     NumHits=12
     ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
     $ProgDir/blast2gff.pl $Column2 $NumHits $BlastHitsCsv > $HitsGff
   done
 ```



### Hybrid assembly:

#### Hybrid assembly: Spades Assembly

```bash
  for PacBioDat in $(ls raw_dna/pacbio/F.oxysporum_fsp_cepae/Fus2/extracted/concatenated_pacbio.fastq); do
    Organism=$(echo $PacBioDat | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $PacBioDat | rev | cut -f3 -d '/' | rev)
    IlluminaDir=$(ls -d qc_dna/paired/$Organism/$Strain)
    TrimF1_Read=$(ls $IlluminaDir/F/s_6_1_sequence_trim.fq.gz);
    TrimR1_Read=$(ls $IlluminaDir/R/s_6_2_sequence_trim.fq.gz);
    TrimF2_Read=$(ls $IlluminaDir/F/FUS2_S2_L001_R1_001_trim.fq.gz);
    TrimR2_Read=$(ls $IlluminaDir/R/FUS2_S2_L001_R2_001_trim.fq.gz);
    OutDir=assembly/spades_pacbio/$Organism/"$Strain"_3
    echo $TrimR1_Read
    echo $TrimR1_Read
    echo $TrimF2_Read
    echo $TrimR2_Read
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades/multiple_libraries
    qsub $ProgDir/subSpades_2lib_pacbio.sh $PacBioDat $TrimF1_Read $TrimR1_Read $TrimF2_Read $TrimR2_Read $OutDir 50
  done
```

Assemblies were polished using Pilon

```bash
  for Assembly in $(ls assembly/spades_pacbio/F.oxysporum_fsp_cepae/Fus2_3/scaffolds.fasta); do
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    IlluminaDir=$(ls -d qc_dna/paired/$Organism/Fus2)
    TrimF1_Read=$(ls $IlluminaDir/F/s_6_1_sequence_trim.fq.gz);
    TrimR1_Read=$(ls $IlluminaDir/R/s_6_2_sequence_trim.fq.gz);
    OutDir=assembly/spades_pacbio/$Organism/$Strain/pilon
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/pilon
    qsub $ProgDir/sub_pilon.sh $Assembly $TrimF1_Read $TrimR1_Read $OutDir
  done
```

```bash
  Assembly=assembly/spades_pacbio/F.oxysporum_fsp_cepae/Fus2_3/pilon/pilon.fasta
  # Assembly=assembly/spades_pacbio/F.oxysporum_fsp_cepae/Fus2_3/scaffolds.fasta
  Reads=raw_dna/pacbio/F.oxysporum_fsp_cepae/Fus2/extracted/concatenated_pacbio.fastq
  OutDir=analysis/genome_alignment/bwa/F.oxysporum_fsp_cepae/Fus2/vs_Fus2_spades_pacbio_3
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/bwa
  qsub $ProgDir/sub_bwa_pacbio.sh $Assembly $Reads $OutDir
```

Inspection of flagged regions didn't identify any contigs that needed to be broken.

<!-- A region of NODE_6 showed low coverage and was manually split.

```bash
  touch tmp.csv
  printf "NODE_7_length_1868836_cov_72.3536_ID_1093330_pilon\t\tsplit\t1703352\t1717240\tspades:missassembly\n" > tmp.csv
  for Assembly in $(ls assembly/spades_pacbio/*/*/pilon/pilon.fasta | grep -w 'Fus2'); do
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    OutDir=$(dirname $Assembly)
    ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
    $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/"$Strain"_canu_contigs_edited_renamed.fasta --coord_file tmp.csv
  done
  rm tmp.csv
``` -->


 Preliminary analysis of these contigs allowed the quality of this assembly to be assessed:

 ```bash
   ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
   for Assembly in $(ls assembly/spades_pacbio/F.oxysporum_fsp_cepae/Fus2_3/pilon/pilon.fasta); do
     echo $Assembly
     Query=analysis/blast_homology/Fo_path_genes/Fo_path_genes_CRX.fa
     qsub $ProgDir/blast_pipe.sh $Query dna $Assembly
   done
 ```

 ```bash
   for Genome in $(ls assembly/spades_pacbio/F.oxysporum_fsp_cepae/Fus2_3/pilon/pilon.fasta); do
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
   for BlastHitsCsv in $(ls analysis/blast_homology/F.oxysporum_fsp_cepae/Fus2_3/*_chr_*_gene_single_copy.aa_hits.csv); do
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
   for BlastHitsCsv in $(ls analysis/blast_homology/F.oxysporum_fsp_cepae/Fus2_3/Fus2_3_Fo_path_genes_CRX.fa_homologs.csv); do
     Organism=$(echo $BlastHitsCsv | rev | cut -f3 -d '/' | rev)
     Strain=$(echo $BlastHitsCsv | rev | cut -f2 -d '/' | rev)
     echo "$Organism - $Strain"
     HitsGff=$(echo $BlastHitsCsv | sed  's/.csv/.gff/g')
     Chr=$(echo $BlastHitsCsv | rev | cut -f5 -d '_' | rev);
     Column2=Fo_path_gene_homolog
     NumHits=12
     ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
     $ProgDir/blast2gff.pl $Column2 $NumHits $BlastHitsCsv > $HitsGff
   done
 ```


## Merging pacbio and hybrid assemblies

```bash
  for PacBioAssembly in $(ls assembly/canu-1.3/F.oxysporum_fsp_cepae/Fus2_canu/pilon/pilon.fasta); do
    Organism=$(echo $PacBioAssembly | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $PacBioAssembly | rev | cut -f3 -d '/' | rev)
    HybridAssembly=$(ls assembly/spades_pacbio/F.oxysporum_fsp_cepae/Fus2_3/scaffolds.fasta)
    AnchorLength=500000
    OutDir=assembly/merged_canu_spades/$Organism/"$Strain"_new
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/quickmerge
    qsub $ProgDir/sub_quickmerge.sh $PacBioAssembly $HybridAssembly $OutDir $AnchorLength
  done
```

This merged assembly was polished using Pilon

```bash
  for Assembly in $(ls assembly/merged_canu_spades/F.oxysporum_fsp_cepae/Fus2_canu_new/merged.fasta); do
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    IlluminaDir=$(ls -d qc_dna/paired/$Organism/Fus2)
    TrimF1_Read=$(ls $IlluminaDir/F/s_6_1_sequence_trim.fq.gz);
    TrimR1_Read=$(ls $IlluminaDir/R/s_6_2_sequence_trim.fq.gz);
    OutDir=assembly/merged_canu_spades/$Organism/$Strain/polished
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/pilon
    qsub $ProgDir/sub_pilon.sh $Assembly $TrimF1_Read $TrimR1_Read $OutDir
  done
```

```bash
  touch tmp.csv
  for Assembly in $(ls assembly/merged_canu_spades/F.oxysporum_fsp_cepae/Fus2_canu_new/polished/pilon.fasta); do
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
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
  Assembly=assembly/merged_canu_spades/F.oxysporum_fsp_cepae/Fus2_canu_new/polished/Fus2_canu_new_contigs_renamed.fasta
  Reads=raw_dna/pacbio/F.oxysporum_fsp_cepae/Fus2/extracted/concatenated_pacbio.fastq
  OutDir=analysis/genome_alignment/bwa/F.oxysporum_fsp_cepae/Fus2/vs_canu_new
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/bwa
  qsub $ProgDir/sub_bwa_pacbio.sh $Assembly $Reads $OutDir
```

Contigs were renamed in and split at a flagged region on contig 10

```bash
  ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
  touch tmp.csv
  printf "contig_10\t\tsplit\t1695182\t1696920\tcanu:missassembly\n" > tmp.csv
  for Assembly in $(ls assembly/merged_canu_spades/F.oxysporum_fsp_cepae/Fus2_canu_new/polished/Fus2_canu_new_contigs_renamed.fasta); do
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    OutDir=assembly/merged_canu_spades/$Organism/$Strain/edited_contigs
    mkdir -p $OutDir
    $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/"$Strain"_contigs_edited.fasta --coord_file tmp.csv
  done
  rm tmp.csv
```

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
  for Assembly in $(ls assembly/merged_canu_spades/F.oxysporum_fsp_cepae/Fus2_canu_new/edited_contigs/Fus2_canu_new_contigs_edited.fasta); do
    echo $Assembly
    Query=analysis/blast_homology/Fo_path_genes/Fo_path_genes_CRX.fa
    qsub $ProgDir/blast_pipe.sh $Query dna $Assembly
  done
```


```bash
  for Genome in $(ls assembly/merged_canu_spades/F.oxysporum_fsp_cepae/Fus2_canu_new/edited_contigs/Fus2_canu_new_contigs_edited.fasta); do
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
  for BlastHitsCsv in $(ls analysis/blast_homology/F.oxysporum_fsp_cepae/Fus2_canu_new/*_chr_*_gene_single_copy.aa_hits.csv); do
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



# Repeatmasking assemblies

```bash
  Fus2_pacbio_canu=$(ls assembly/canu-1.3/F.oxysporum_fsp_cepae/Fus2_canu/ncbi_edits/contigs_min_500bp_renamed.fasta)
  for Assembly in $(echo $Fus2_pacbio_canu); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    OutDir=repeat_masked/$Organism/"$Strain"_ncbi/ncbi_submission
    OutDir=repeat_masked/F.oxysporum_fsp_cepae/Fus2_canu_ncbi/edited_contigs_repmask
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking
    qsub $ProgDir/rep_modeling.sh $Assembly $OutDir
    qsub $ProgDir/transposonPSI.sh $Assembly $OutDir
  done
  Fus2_pacbio_merged=$(ls assembly/merged_canu_spades/F.oxysporum_fsp_cepae/Fus2_canu_new/edited_contigs/Fus2_canu_new_contigs_edited.fasta)
  for Assembly in $(echo $Fus2_pacbio_merged); do
    OutDir=repeat_masked/F.oxysporum_fsp_cepae/Fus2_merged/edited_contigs_repmask
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking
    qsub $ProgDir/rep_modeling.sh $Assembly $OutDir
    qsub $ProgDir/transposonPSI.sh $Assembly $OutDir
  done
```

```bash
  for File in $(ls repeat_masked/*/*/filtered_contigs_repmask/*_contigs_softmasked.fa); do
    OutDir=$(dirname $File)
    TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
    OutFile=$(echo $File | sed 's/_contigs_softmasked.fa/_contigs_softmasked_repeatmasker_TPSI_appended.fa/g')
    bedtools maskfasta -soft -fi $File -bed $TPSI -fo $OutFile
  done
```


The number of bases masked by transposonPSI and Repeatmasker were summarised
using the following commands:

```bash
for RepDir in $(ls -d repeat_masked/F.*/*/*); do
Strain=$(echo $RepDir | rev | cut -f2 -d '/' | rev)
Organism=$(echo $RepDir | rev | cut -f3 -d '/' | rev)  
RepMaskGff=$(ls $RepDir/*_contigs_hardmasked.gff)
TransPSIGff=$(ls $RepDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
printf "$Organism\t$Strain\n"
printf "The number of bases masked by RepeatMasker:\t"
sortBed -i $RepMaskGff | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
printf "The number of bases masked by TransposonPSI:\t"
sortBed -i $TransPSIGff | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
printf "The total number of masked bases are:\t"
cat $RepMaskGff $TransPSIGff | sortBed | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
echo
done | cut -f2
```


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
