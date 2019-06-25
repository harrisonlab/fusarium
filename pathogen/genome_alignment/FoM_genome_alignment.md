
# 5 Promer alignment of Assemblies

## 5.1 against the reference f.sp. lycopersici genome

```bash
Reference=$(ls assembly/external_group/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl/Fusarium_oxysporum_chromosome_and_additional_contigs.fa)
for Query in $(ls repeat_masked/*/*/filtered_ncbi/*_contigs_unmasked.fa | grep -e 'Stocks4'); do
Strain=$(echo $Query | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Query | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
Prefix="$Strain"_vs_4287
OutDir=analysis/genome_alignment/mummer/$Organism/$Strain/$Prefix
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/MUMmer
qsub $ProgDir/sub_nucmer.sh $Reference $Query $Prefix $OutDir
done
```

## Alignment of FoM raw reads vs reference genomes


Alignment of reads from a single run:

```bash
for Reference in $(ls assembly/external_group/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl/Fusarium_oxysporum_chromosome_and_additional_contigs.fa repeat_masked/F.oxysporum_fsp_cepae/Fus2_canu_new/edited_contigs_repmask/Fus2_canu_contigs_unmasked.fa); do
for StrainPath in $(ls -d qc_dna/paired/*/* | grep 'Stocks4'); do
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades
Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
F_Read=$(ls $StrainPath/F/*_trim.fq.gz)
R_Read=$(ls $StrainPath/R/*_trim.fq.gz)
echo $F_Read
echo $R_Read
Prefix="${Organism}_${Strain}"
OutDir=/home/groups/harrisonlab/project_files/fusarium/analysis/genome_alignment/bwa/$Organism/$Strain/vs_${Reference}
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment
qsub $ProgDir/bwa/sub_bwa.sh $Prefix $Reference $F_Read $R_Read $OutDir
done
done
```

Identify read coverage over each bp

```bash
  cd /data/scratch/armita/magnaporthe
  for Bam in $(ls alignment/*/*_sorted.bam | tail -n+2); do
    Ref=$(echo $Bam | cut -f2 -d '/' | sed 's/Align//g')
    Strain=$(basename ${Bam%_sorted.bam})
    Organism="M.grisea"
    echo "$Organism - $Strain"
    # OutDir=$(dirname $Bam)
    OutDir="alignment/genome_alignment/vs_${Ref}/$Organism/$Strain/"
    mkdir -p $OutDir
    samtools depth -aa $Bam > $OutDir/${Organism}_${Strain}_vs_${Ref}_depth.tsv
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/coverage_analysis
    $ProgDir/cov_by_window.py --cov $OutDir/${Organism}_${Strain}_vs_${Ref}_depth.tsv > $OutDir/${Organism}_${Strain}_vs_${Ref}_depth_10kb.tsv
    sed -i "s/$/\t$Strain/g" $OutDir/${Organism}_${Strain}_vs_${Ref}_depth_10kb.tsv
  done
for RefDir in $(ls -d alignment/genome_alignment/vs_*); do
Ref=$(echo $RefDir | rev | cut -f1 -d '/'| rev | sed 's/vs_//g')
echo $Ref
cat $RefDir/*/*/*_vs_${Ref}_depth_10kb.tsv > $RefDir/vs_${Ref}_grouped_depth.tsv
done
```


## Alignment of Reference genome reads vs the FoM genomes


Alignment of reads from a multiple runs:

```bash
for Reference in $(ls repeat_masked/F.oxysporum_fsp_mathioli/Stocks4/filtered_ncbi/Stocks4_contigs_unmasked.fa); do
for StrainPath in $(ls -d qc_dna/paired/*/* | grep 'Fus2'); do
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades
Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
F1_Read=$(ls $StrainPath/F/*_trim.fq.gz | head -n1 | tail -n1);
R1_Read=$(ls $StrainPath/R/*_trim.fq.gz | head -n1 | tail -n1);
F2_Read=$(ls $StrainPath/F/*_trim.fq.gz | head -n2 | tail -n1);
R2_Read=$(ls $StrainPath/R/*_trim.fq.gz | head -n2 | tail -n1);
echo $F1_Read
echo $R1_Read
echo $F2_Read
echo $R2_Read
Jobs=$(qstat | grep 'sub_bwa' | grep 'qw'| wc -l)
while [ $Jobs -gt 1 ]; do
sleep 1m
printf "."
Jobs=$(qstat | grep 'sub_bwa' | grep 'qw'| wc -l)
done
echo $F_Read
echo $R_Read
Prefix="${Organism}_${Strain}"
OutDir=/home/groups/harrisonlab/project_files/fusarium/analysis/genome_alignment/bwa/$Organism/$Strain/vs_${Reference}
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment
qsub $ProgDir/bowtie/sub_bowtie_2lib.sh $Reference $F1_Read $R1_Read $F2_Read $R2_Read $OutDir $Strain
done
done
```

Identify read coverage over each bp

```bash
  cd /data/scratch/armita/magnaporthe
  for Bam in $(ls alignment/*/*_sorted.bam | tail -n+2); do
    Ref=$(echo $Bam | cut -f2 -d '/' | sed 's/Align//g')
    Strain=$(basename ${Bam%_sorted.bam})
    Organism="M.grisea"
    echo "$Organism - $Strain"
    # OutDir=$(dirname $Bam)
    OutDir="alignment/genome_alignment/vs_${Ref}/$Organism/$Strain/"
    mkdir -p $OutDir
    samtools depth -aa $Bam > $OutDir/${Organism}_${Strain}_vs_${Ref}_depth.tsv
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/coverage_analysis
    $ProgDir/cov_by_window.py --cov $OutDir/${Organism}_${Strain}_vs_${Ref}_depth.tsv > $OutDir/${Organism}_${Strain}_vs_${Ref}_depth_10kb.tsv
    sed -i "s/$/\t$Strain/g" $OutDir/${Organism}_${Strain}_vs_${Ref}_depth_10kb.tsv
  done
for RefDir in $(ls -d alignment/genome_alignment/vs_*); do
Ref=$(echo $RefDir | rev | cut -f1 -d '/'| rev | sed 's/vs_//g')
echo $Ref
cat $RefDir/*/*/*_vs_${Ref}_depth_10kb.tsv > $RefDir/vs_${Ref}_grouped_depth.tsv
done
```
