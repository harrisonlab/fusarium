Assembly of the pacbio genome was optimised through assembly and merging under different parameters.

The commands used are documented here.

# 0. Preceeding commands

These commands were part of the original assembly:


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
  printf "contig_17\t\tsplit\t780978\t780971\tcanu:missassembly\n" > tmp.csv
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
  for Assembly in $(ls assembly/canu/*/*/filtered_contigs/*_canu_contigs_renamed.fasta | grep -w 'Fus2'); do
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)  
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
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


# 1. Assessment of the canu only assembly


```bash
  Fus2_pacbio_merged=$(ls assembly/canu/F.oxysporum_fsp_cepae/Fus2/filtered_contigs/Fus2_canu_contigs_renamed.fasta)
  for Genome in $(ls $Fus2_pacbio_merged); do
    for Proteome in $(ls analysis/FoL_genes/chr_*_gene_single_copy.aa); do
      Organism=$(echo $Genome | rev | cut -f4 -d '/' | rev)
      Strain=$(echo $Genome | rev | cut -f3 -d '/' | rev)
      Chr=$(echo $Proteome | rev | cut -f4 -d '_' | rev)
      echo "$Organism - $Strain - $Chr"
      OutDir=analysis/blast_homology/$Organism/"$Strain"_canu_only
      ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
      qsub $ProgDir/run_blast2csv.sh $Proteome protein $Genome $OutDir
    done
  done
```

Convert top blast hits into gff annotations

```bash
  for BlastHitsCsv in $(ls analysis/blast_homology/F.oxysporum_fsp_cepae/Fus2_canu_only/*_chr_*_gene_single_copy.aa_hits.csv); do
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

# 2. Assessment of the manually edited assembly


```bash
  Fus2_pacbio_merged=$(ls assembly/canu/F.oxysporum_fsp_cepae/Fus2_edited3/Fus2_canu_manual_edits.fasta)
  for Genome in $(ls $Fus2_pacbio_merged); do
    for Proteome in $(ls analysis/FoL_genes/chr_*_gene_single_copy.aa); do
      Organism=$(echo $Genome | rev | cut -f3 -d '/' | rev)
      Strain=$(echo $Genome | rev | cut -f2 -d '/' | rev)
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
  for BlastHitsCsv in $(ls analysis/blast_homology/F.oxysporum_fsp_cepae/Fus2_edited3/*_chr_*_gene_single_copy.aa_hits.csv); do
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


# 3. Assembly of Quickmerged assemblies


## 3.1 Quickmerge with default parameters

AnchorOverlapCuttoff=5.0
ExtensionOverlapCuttoff=1.5
AnchorContigLength=0
MinAlignmentLength=0

```bash
  for PacBioAssembly in $(ls assembly/canu/F.oxysporum_fsp_cepae/Fus2/filtered_contigs/Fus2_canu_contigs_renamed.fasta); do
    Organism=$(echo $PacBioAssembly | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $PacBioAssembly | rev | cut -f3 -d '/' | rev)
    HybridAssembly=$(ls assembly/spades_pacbio/$Organism/$Strain/contigs.fasta)
    QuastReport=$(ls assembly/canu/$Organism/$Strain/filtered_contigs/report.tsv)
    # N50=$(cat $QuastReport | grep 'N50' | cut -f2)
    # AnchorLength=$N50
    AnchorLength=0
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


```bash
  Fus2_pacbio_merged=$(ls assembly/merged_canu_spades/F.oxysporum_fsp_cepae/Fus2_edited2/Fus2_edited2_contigs_renamed.fasta)
  for Genome in $(ls $Fus2_pacbio_merged); do
    for Proteome in $(ls analysis/FoL_genes/chr_*_gene_single_copy.aa); do
      Organism=$(echo $Genome | rev | cut -f3 -d '/' | rev)
      Strain=$(echo $Genome | rev | cut -f2 -d '/' | rev)
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
  for BlastHitsCsv in $(ls analysis/blast_homology/F.oxysporum_fsp_cepae/Fus2_edited2/*_chr_*_gene_single_copy.aa_hits.csv); do
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


<!--

This merged assembly was polished using Pilon

```bash
  for Assembly in $(ls assembly/merged_canu_spades/F.oxysporum_fsp_cepae/Fus2_edited2/merged.fasta); do
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
 -->

## 3.2 Quickmerge with doubled parameters


AnchorOverlapCuttoff=10.0
ExtensionOverlapCuttoff=3.0
AnchorContigLength=20000
MinAlignmentLength=5000

```bash
  for PacBioAssembly in $(ls assembly/canu/F.oxysporum_fsp_cepae/Fus2/filtered_contigs/Fus2_canu_contigs_renamed.fasta); do
    Organism=$(echo $PacBioAssembly | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $PacBioAssembly | rev | cut -f3 -d '/' | rev)
    HybridAssembly=$(ls assembly/spades_pacbio/$Organism/$Strain/contigs.fasta)
    QuastReport=$(ls assembly/canu/$Organism/$Strain/filtered_contigs/report.tsv)
    AnchorLength=20000
    OutDir=assembly/merged_canu_spades/$Organism/"$Strain"_quickmerge2
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/quickmerge
    qsub $ProgDir/sub_quickmerge.sh $PacBioAssembly $HybridAssembly $OutDir $AnchorLength
  done
```

```bash
  touch tmp.csv
  for Assembly in $(ls assembly/merged_canu_spades/F.oxysporum_fsp_cepae/Fus2_quickmerge2/merged.fasta); do
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
  for Assembly in $(ls assembly/merged_canu_spades/F.oxysporum_fsp_cepae/Fus2_quickmerge2/Fus2_quickmerge2_contigs_renamed.fasta); do
    echo $Assembly
    Query=analysis/blast_homology/Fo_path_genes/Fo_path_genes_CRX.fa
    qsub $ProgDir/blast_pipe.sh $Query dna $Assembly
  done
```

```bash
  Fus2_pacbio_merged=$(ls assembly/merged_canu_spades/F.oxysporum_fsp_cepae/Fus2_quickmerge2/Fus2_quickmerge2_contigs_renamed.fasta)
  for Genome in $(ls $Fus2_pacbio_merged); do
    for Proteome in $(ls analysis/FoL_genes/chr_*_gene_single_copy.aa); do
      Organism=$(echo $Genome | rev | cut -f3 -d '/' | rev)
      Strain=$(echo $Genome | rev | cut -f2 -d '/' | rev)
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
  for BlastHitsCsv in $(ls analysis/blast_homology/F.oxysporum_fsp_cepae/Fus2_quickmerge2/*_chr_*_gene_single_copy.aa_hits.csv); do
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

## 3.3 Quickmerge with 5X parameters


AnchorOverlapCuttoff=25.0
ExtensionOverlapCuttoff=7.5
AnchorContigLength=20000
MinAlignmentLength=5000

```bash
  for PacBioAssembly in $(ls assembly/canu/F.oxysporum_fsp_cepae/Fus2/filtered_contigs/Fus2_canu_contigs_renamed.fasta); do
    Organism=$(echo $PacBioAssembly | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $PacBioAssembly | rev | cut -f3 -d '/' | rev)
    HybridAssembly=$(ls assembly/spades_pacbio/$Organism/$Strain/contigs.fasta)
    QuastReport=$(ls assembly/canu/$Organism/$Strain/filtered_contigs/report.tsv)
    AnchorLength=20000
    OutDir=assembly/merged_canu_spades/$Organism/"$Strain"_quickmerge3
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/quickmerge
    qsub $ProgDir/sub_quickmerge.sh $PacBioAssembly $HybridAssembly $OutDir $AnchorLength
  done
```

```bash
  touch tmp.csv
  for Assembly in $(ls assembly/merged_canu_spades/F.oxysporum_fsp_cepae/Fus2_quickmerge3/merged.fasta); do
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
  for Assembly in $(ls assembly/merged_canu_spades/F.oxysporum_fsp_cepae/Fus2_quickmerge3/Fus2_quickmerge3_contigs_renamed.fasta); do
    echo $Assembly
    Query=analysis/blast_homology/Fo_path_genes/Fo_path_genes_CRX.fa
    qsub $ProgDir/blast_pipe.sh $Query dna $Assembly
  done
```

```bash
  Fus2_pacbio_merged=$(ls assembly/merged_canu_spades/F.oxysporum_fsp_cepae/Fus2_quickmerge3/Fus2_quickmerge3_contigs_renamed.fasta)
  for Genome in $(ls $Fus2_pacbio_merged); do
    for Proteome in $(ls analysis/FoL_genes/chr_*_gene_single_copy.aa); do
      Organism=$(echo $Genome | rev | cut -f3 -d '/' | rev)
      Strain=$(echo $Genome | rev | cut -f2 -d '/' | rev)
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
  for BlastHitsCsv in $(ls analysis/blast_homology/F.oxysporum_fsp_cepae/Fus2_quickmerge3/*_chr_*_gene_single_copy.aa_hits.csv); do
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


## 3.4 Quickmerge with default but N50 anchor


AnchorOverlapCuttoff=5.0
ExtensionOverlapCuttoff=1.5
AnchorContigLength=N50
MinAlignmentLength=0

```bash
  for PacBioAssembly in $(ls assembly/canu/F.oxysporum_fsp_cepae/Fus2/filtered_contigs/Fus2_canu_contigs_renamed.fasta); do
    Organism=$(echo $PacBioAssembly | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $PacBioAssembly | rev | cut -f3 -d '/' | rev)
    HybridAssembly=$(ls assembly/spades_pacbio/$Organism/$Strain/contigs.fasta)
    QuastReport=$(ls assembly/canu/$Organism/$Strain/filtered_contigs/report.tsv)
    N50=$(cat $QuastReport | grep 'N50' | cut -f2)
    AnchorLength=$N50
    OutDir=assembly/merged_canu_spades/$Organism/"$Strain"_quickmerge4
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/quickmerge
    qsub $ProgDir/sub_quickmerge.sh $PacBioAssembly $HybridAssembly $OutDir $AnchorLength
  done
```

```bash
  touch tmp.csv
  for Assembly in $(ls assembly/merged_canu_spades/F.oxysporum_fsp_cepae/Fus2_quickmerge4/merged.fasta); do
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
  for Assembly in $(ls assembly/merged_canu_spades/F.oxysporum_fsp_cepae/Fus2_quickmerge4/Fus2_quickmerge4_contigs_renamed.fasta); do
    echo $Assembly
    Query=analysis/blast_homology/Fo_path_genes/Fo_path_genes_CRX.fa
    qsub $ProgDir/blast_pipe.sh $Query dna $Assembly
  done
```

```bash
  Fus2_pacbio_merged=$(ls assembly/merged_canu_spades/F.oxysporum_fsp_cepae/Fus2_quickmerge4/Fus2_quickmerge4_contigs_renamed.fasta)
  for Genome in $(ls $Fus2_pacbio_merged); do
    for Proteome in $(ls analysis/FoL_genes/chr_*_gene_single_copy.aa); do
      Organism=$(echo $Genome | rev | cut -f3 -d '/' | rev)
      Strain=$(echo $Genome | rev | cut -f2 -d '/' | rev)
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
  for BlastHitsCsv in $(ls analysis/blast_homology/F.oxysporum_fsp_cepae/Fus2_quickmerge4/*_chr_*_gene_single_copy.aa_hits.csv); do
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


## 3.5 Quickmerge with 5X parameters and large min alignment length


AnchorOverlapCuttoff=25.0
ExtensionOverlapCuttoff=7.5
AnchorContigLength=20000
MinAlignmentLength=10000

```bash
  for PacBioAssembly in $(ls assembly/canu/F.oxysporum_fsp_cepae/Fus2/filtered_contigs/Fus2_canu_contigs_renamed.fasta); do
    Organism=$(echo $PacBioAssembly | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $PacBioAssembly | rev | cut -f3 -d '/' | rev)
    HybridAssembly=$(ls assembly/spades_pacbio/$Organism/$Strain/contigs.fasta)
    QuastReport=$(ls assembly/canu/$Organism/$Strain/filtered_contigs/report.tsv)
    AnchorLength=20000
    OutDir=assembly/merged_canu_spades/$Organism/"$Strain"_quickmerge5
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/quickmerge
    qsub $ProgDir/sub_quickmerge.sh $PacBioAssembly $HybridAssembly $OutDir $AnchorLength
  done
```

```bash
  touch tmp.csv
  for Assembly in $(ls assembly/merged_canu_spades/F.oxysporum_fsp_cepae/Fus2_quickmerge5/merged.fasta); do
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
  for Assembly in $(ls assembly/merged_canu_spades/F.oxysporum_fsp_cepae/Fus2_quickmerge5/Fus2_quickmerge5_contigs_renamed.fasta); do
    echo $Assembly
    Query=analysis/blast_homology/Fo_path_genes/Fo_path_genes_CRX.fa
    qsub $ProgDir/blast_pipe.sh $Query dna $Assembly
  done
```

```bash
  Fus2_pacbio_merged=$(ls assembly/merged_canu_spades/F.oxysporum_fsp_cepae/Fus2_quickmerge5/Fus2_quickmerge5_contigs_renamed.fasta)
  for Genome in $(ls $Fus2_pacbio_merged); do
    for Proteome in $(ls analysis/FoL_genes/chr_*_gene_single_copy.aa); do
      Organism=$(echo $Genome | rev | cut -f3 -d '/' | rev)
      Strain=$(echo $Genome | rev | cut -f2 -d '/' | rev)
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
  for BlastHitsCsv in $(ls analysis/blast_homology/F.oxysporum_fsp_cepae/Fus2_quickmerge5/*_chr_*_gene_single_copy.aa_hits.csv); do
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


## 3.6 Quickmerge with recomended Anchor length


AnchorOverlapCuttoff=5.0
ExtensionOverlapCuttoff=1.5
AnchorContigLength=500000
MinAlignmentLength=5000

```bash
  for PacBioAssembly in $(ls assembly/canu/F.oxysporum_fsp_cepae/Fus2/filtered_contigs/Fus2_canu_contigs_renamed.fasta); do
    Organism=$(echo $PacBioAssembly | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $PacBioAssembly | rev | cut -f3 -d '/' | rev)
    HybridAssembly=$(ls assembly/spades_pacbio/$Organism/$Strain/edited_contigs/pilon.fasta)
    QuastReport=$(ls assembly/canu/$Organism/$Strain/filtered_contigs/report.tsv)
    AnchorLength=500000
    OutDir=assembly/merged_canu_spades/$Organism/"$Strain"_quickmerge6
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/quickmerge
    qsub $ProgDir/sub_quickmerge.sh $PacBioAssembly $HybridAssembly $OutDir $AnchorLength
  done
```

```bash
  touch tmp.csv
  for Assembly in $(ls assembly/merged_canu_spades/F.oxysporum_fsp_cepae/Fus2_quickmerge3/merged.fasta); do
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
  for Assembly in $(ls assembly/merged_canu_spades/F.oxysporum_fsp_cepae/Fus2_quickmerge6/Fus2_quickmerge6_contigs_renamed.fasta); do
    echo $Assembly
    Query=analysis/blast_homology/Fo_path_genes/Fo_path_genes_CRX.fa
    qsub $ProgDir/blast_pipe.sh $Query dna $Assembly
  done
```

```bash
  Fus2_pacbio_merged=$(ls assembly/merged_canu_spades/F.oxysporum_fsp_cepae/Fus2_quickmerge6/Fus2_quickmerge6_contigs_renamed.fasta)
  for Genome in $(ls $Fus2_pacbio_merged); do
    for Proteome in $(ls analysis/FoL_genes/chr_*_gene_single_copy.aa); do
      Organism=$(echo $Genome | rev | cut -f3 -d '/' | rev)
      Strain=$(echo $Genome | rev | cut -f2 -d '/' | rev)
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
  for BlastHitsCsv in $(ls analysis/blast_homology/F.oxysporum_fsp_cepae/Fus2_quickmerge6/*_chr_*_gene_single_copy.aa_hits.csv); do
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


# 4.0 Quickmerge of split contigs

# 4.1 Submission with low stringency

AnchorOverlapCuttoff=5.0
ExtensionOverlapCuttoff=1.5
AnchorContigLength=0
MinAlignmentLength=0

```bash
  for PacBioAssembly in $(ls assembly/canu/F.oxysporum_fsp_cepae/Fus2/edited_contigs/pilon.fasta); do
    Organism=$(echo $PacBioAssembly | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $PacBioAssembly | rev | cut -f3 -d '/' | rev)
    HybridAssembly=$(ls assembly/spades_pacbio/$Organism/$Strain/contigs.fasta)
    QuastReport=$(ls assembly/canu/$Organism/$Strain/filtered_contigs/report.tsv)
    # N50=$(cat $QuastReport | grep 'N50' | cut -f2)
    # AnchorLength=$N50
    AnchorLength=0
    OutDir=assembly/merged_canu_spades/$Organism/"$Strain"_edited_quickmerge1
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/quickmerge
    qsub $ProgDir/sub_quickmerge.sh $PacBioAssembly $HybridAssembly $OutDir $AnchorLength
  done
```

```bash
  touch tmp.csv
  for Assembly in $(ls assembly/merged_canu_spades/F.oxysporum_fsp_cepae/Fus2_edited_quickmerge1/merged.fasta); do
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
  for Assembly in $(ls assembly/merged_canu_spades/F.oxysporum_fsp_cepae/Fus2_edited_quickmerge1/Fus2_edited_quickmerge1_contigs_renamed.fasta); do
    echo $Assembly
    Query=analysis/blast_homology/Fo_path_genes/Fo_path_genes_CRX.fa
    qsub $ProgDir/blast_pipe.sh $Query dna $Assembly
  done
```


```bash
  Fus2_pacbio_merged=$(ls assembly/merged_canu_spades/F.oxysporum_fsp_cepae/Fus2_edited_quickmerge1/Fus2_edited_quickmerge1_contigs_renamed.fasta)
  for Genome in $(ls $Fus2_pacbio_merged); do
    for Proteome in $(ls analysis/FoL_genes/chr_*_gene_single_copy.aa); do
      Organism=$(echo $Genome | rev | cut -f3 -d '/' | rev)
      Strain=$(echo $Genome | rev | cut -f2 -d '/' | rev)
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
  for BlastHitsCsv in $(ls analysis/blast_homology/F.oxysporum_fsp_cepae/Fus2_edited_quickmerge1/*_chr_*_gene_single_copy.aa_hits.csv); do
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


# 4.2 Submission with 5X stringency

AnchorOverlapCuttoff=25.0
ExtensionOverlapCuttoff=7.5
AnchorContigLength=20000
MinAlignmentLength=5000

```bash
  for PacBioAssembly in $(ls assembly/canu/F.oxysporum_fsp_cepae/Fus2/edited_contigs/pilon.fasta); do
    Organism=$(echo $PacBioAssembly | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $PacBioAssembly | rev | cut -f3 -d '/' | rev)
    HybridAssembly=$(ls assembly/spades_pacbio/$Organism/$Strain/contigs.fasta)
    QuastReport=$(ls assembly/canu/$Organism/$Strain/filtered_contigs/report.tsv)
    # N50=$(cat $QuastReport | grep 'N50' | cut -f2)
    # AnchorLength=$N50
    AnchorLength=20000
    OutDir=assembly/merged_canu_spades/$Organism/"$Strain"_edited_quickmerge2
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/quickmerge
    qsub $ProgDir/sub_quickmerge.sh $PacBioAssembly $HybridAssembly $OutDir $AnchorLength
  done
```

```bash
  touch tmp.csv
  for Assembly in $(ls assembly/merged_canu_spades/F.oxysporum_fsp_cepae/Fus2_edited_quickmerge2/merged.fasta); do
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
  for Assembly in $(ls assembly/merged_canu_spades/F.oxysporum_fsp_cepae/Fus2_edited_quickmerge2/Fus2_edited_quickmerge2_contigs_renamed.fasta); do
    echo $Assembly
    Query=analysis/blast_homology/Fo_path_genes/Fo_path_genes_CRX.fa
    qsub $ProgDir/blast_pipe.sh $Query dna $Assembly
  done
```


```bash
  Fus2_pacbio_merged=$(ls assembly/merged_canu_spades/F.oxysporum_fsp_cepae/Fus2_edited_quickmerge2/Fus2_edited_quickmerge2_contigs_renamed.fasta)
  for Genome in $(ls $Fus2_pacbio_merged); do
    for Proteome in $(ls analysis/FoL_genes/chr_*_gene_single_copy.aa); do
      Organism=$(echo $Genome | rev | cut -f3 -d '/' | rev)
      Strain=$(echo $Genome | rev | cut -f2 -d '/' | rev)
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
  for BlastHitsCsv in $(ls analysis/blast_homology/F.oxysporum_fsp_cepae/Fus2_edited_quickmerge2/*_chr_*_gene_single_copy.aa_hits.csv); do
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


## 4.4 Quickmerge with recomended Anchor length, split pacbio and split, pilon edited hybrid assembly


Contigs were renamed in accordance with ncbi recomendations

```bash
  touch tmp.csv
  printf "NODE_6_length_1868837_cov_72.353_ID_1107516_pilon\t\tsplit\t1703381\t1717033\tspades:missassembly\n" > tmp.csv
  for Assembly in $(ls assembly/spades_pacbio/*/*/edited_contigs/pilon.fasta | grep -w 'Fus2'); do
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    OutDir=$(dirname $Assembly)
    ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
    $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/"$Strain"_canu_contigs_edited_renamed.fasta --coord_file tmp.csv
  done
  rm tmp.csv
```


AnchorOverlapCuttoff=5.0
ExtensionOverlapCuttoff=1.5
AnchorContigLength=500000
MinAlignmentLength=5000

```bash
  for PacBioAssembly in $(ls assembly/canu/F.oxysporum_fsp_cepae/Fus2/edited_contigs/pilon.fasta); do
    Organism=$(echo $PacBioAssembly | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $PacBioAssembly | rev | cut -f3 -d '/' | rev)
    HybridAssembly=$(ls assembly/spades_pacbio/$Organism/$Strain/edited_contigs/"$Strain"_canu_contigs_edited_renamed.fasta)
    QuastReport=$(ls assembly/canu/$Organism/$Strain/filtered_contigs/report.tsv)
    AnchorLength=500000
    OutDir=assembly/merged_canu_spades/$Organism/"$Strain"_edited_quickmerge3
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/quickmerge
    qsub $ProgDir/sub_quickmerge.sh $PacBioAssembly $HybridAssembly $OutDir $AnchorLength
  done
```

```bash
  touch tmp.csv
  for Assembly in $(ls assembly/merged_canu_spades/F.oxysporum_fsp_cepae/Fus2_edited_quickmerge3/merged.fasta); do
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
  for Assembly in $(ls assembly/merged_canu_spades/F.oxysporum_fsp_cepae/Fus2_edited_quickmerge3/Fus2_edited_quickmerge3_contigs_renamed.fasta); do
    echo $Assembly
    Query=analysis/blast_homology/Fo_path_genes/Fo_path_genes_CRX.fa
    qsub $ProgDir/blast_pipe.sh $Query dna $Assembly
  done
```

```bash
  for Genome in $(ls assembly/merged_canu_spades/F.oxysporum_fsp_cepae/Fus2_edited_quickmerge3/Fus2_edited_quickmerge3_contigs_renamed.fasta); do
    for Proteome in $(ls analysis/FoL_genes/chr_*_gene_single_copy.aa); do
      Organism=$(echo $Genome | rev | cut -f3 -d '/' | rev)
      Strain=$(echo $Genome | rev | cut -f2 -d '/' | rev)
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
  for BlastHitsCsv in $(ls analysis/blast_homology/F.oxysporum_fsp_cepae/Fus2_edited_quickmerge3/*_chr_*_gene_single_copy.aa_hits.csv); do
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

```bash
  Assembly=aassembly/merged_canu_spades/F.oxysporum_fsp_cepae/Fus2_edited_quickmerge3/Fus2_edited_quickmerge3_contigs_renamed.fasta
  Reads=raw_dna/pacbio/F.oxysporum_fsp_cepae/Fus2/extracted/concatenated_pacbio.fastq
  OutDir=analysis/genome_alignment/bwa/F.oxysporum_fsp_cepae/Fus2/vs_Fus2_edited_quickmerge3
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/bwa
  qsub $ProgDir/sub_bwa_pacbio.sh $Assembly $Reads $OutDir
```


# 5.0 Quickmerge of manually merged contigs

```bash
  for Assembly in $(ls assembly/canu/F.oxysporum_fsp_cepae/Fus2_edited3/Fus2_canu_manual_edits.fasta); do
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    IlluminaDir=$(ls -d qc_dna/paired/$Organism/Fus2)
    TrimF1_Read=$(ls $IlluminaDir/F/s_6_1_sequence_trim.fq.gz);
    TrimR1_Read=$(ls $IlluminaDir/R/s_6_2_sequence_trim.fq.gz);
    OutDir=assembly/canu/$Organism/$Strain/edited_contigs
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/pilon
    qsub $ProgDir/sub_pilon.sh $Assembly $TrimF1_Read $TrimR1_Read $OutDir
  done
```

Assembly stats were collected using quast

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly/canu/F.oxysporum_fsp_cepae/Fus2_edited3/edited_contigs/pilon.fasta); do
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)  
    OutDir=assembly/canu/$Organism/$Strain/pilon
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```

# 5.1 Submission with low stringency

AnchorOverlapCuttoff=5.0
ExtensionOverlapCuttoff=1.5
AnchorContigLength=0
MinAlignmentLength=0


```bash
  for PacBioAssembly in $(ls assembly/canu/F.oxysporum_fsp_cepae/Fus2_edited3/edited_contigs/pilon.fasta); do
    Organism=$(echo $PacBioAssembly | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $PacBioAssembly | rev | cut -f3 -d '/' | rev)
    # Strain=Fus2
    HybridAssembly=$(ls assembly/spades_pacbio/$Organism/Fus2/contigs.fasta)
    # QuastReport=$(ls assembly/canu/$Organism/$Strain/filtered_contigs/report.tsv)
    AnchorLength=0
    OutDir=assembly/merged_canu_spades/$Organism/Fus2_manual_edit_1
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/quickmerge
    qsub $ProgDir/sub_quickmerge.sh $PacBioAssembly $HybridAssembly $OutDir $AnchorLength
  done
```

```bash
  touch tmp.csv
  for Assembly in $(ls assembly/merged_canu_spades/*/Fus2_manual_edit_1/merged.fasta); do
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
  for Assembly in $(ls assembly/merged_canu_spades/F.oxysporum_fsp_cepae/Fus2_manual_edit_1/Fus2_manual_edit_1_contigs_renamed.fasta); do
    echo $Assembly
    Query=analysis/blast_homology/Fo_path_genes/Fo_path_genes_CRX.fa
    qsub $ProgDir/blast_pipe.sh $Query dna $Assembly
  done
```


```bash
  Fus2_pacbio_merged=$(ls assembly/merged_canu_spades/F.oxysporum_fsp_cepae/Fus2_manual_edit_1/Fus2_manual_edit_1_contigs_renamed.fasta)
  for Genome in $(ls $Fus2_pacbio_merged); do
    for Proteome in $(ls analysis/FoL_genes/chr_*_gene_single_copy.aa); do
      Organism=$(echo $Genome | rev | cut -f3 -d '/' | rev)
      Strain=$(echo $Genome | rev | cut -f2 -d '/' | rev)
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
  for BlastHitsCsv in $(ls analysis/blast_homology/F.oxysporum_fsp_cepae/Fus2_manual_edit_1/*_chr_*_gene_single_copy.aa_hits.csv); do
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
