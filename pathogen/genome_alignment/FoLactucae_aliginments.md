
# 5 Promer alignment of Assemblies

## 5.1 against AJ516 genome

MUMmer was run to align assemblies against the reference genome.

```bash
Isolate=AJ516
Reference=$(ls repeat_masked/*/*/filtered_contigs/${Isolate}_contigs_unmasked.fa)
for Query in $(ls repeat_masked/*/*/filtered_contigs/AJ520_contigs_unmasked.fa); do
Strain=$(echo $Query | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Query | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
Prefix="$Strain"_vs_${Isolate}
OutDir=analysis/genome_alignment/mummer/$Organism/$Strain/$Prefix
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/MUMmer
qsub $ProgDir/sub_nucmer.sh $Reference $Query $Prefix $OutDir
done
```

## 5.1 against AJ520 genome

```bash
Isolate=AJ520
Reference=$(ls repeat_masked/*/*/filtered_contigs/${Isolate}_contigs_unmasked.fa)
for Query in $(ls repeat_masked/*/*/filtered_contigs/AJ516_contigs_unmasked.fa); do
Strain=$(echo $Query | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Query | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
Prefix="$Strain"_vs_${Isolate}
OutDir=analysis/genome_alignment/mummer/$Organism/$Strain/$Prefix
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/MUMmer
qsub $ProgDir/sub_nucmer.sh $Reference $Query $Prefix $OutDir
done
```

## 5.1 against the reference f.sp. lycopersici genome

```bash
Reference=$(ls repeat_masked/F.oxysporum_fsp_lycopersici/4287_v2/fungidb_repmask/4287_v2_contigs_unmasked.fa)
for Query in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_unmasked.fa | grep -e 'AJ516' -e 'AJ520'); do
Strain=$(echo $Query | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Query | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
Prefix="$Strain"_vs_4287
OutDir=analysis/genome_alignment/mummer/$Organism/$Strain/$Prefix
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/MUMmer
qsub $ProgDir/sub_nucmer.sh $Reference $Query $Prefix $OutDir
done
```
