
# 5 Alignment of Assemblies

## 5.1 against the reference f.sp. lycopersici genome

```bash
Reference=$(ls repeat_masked/F.oxysporum_fsp_lycopersici/4287_v2/fungidb_repmask/4287_v2_contigs_unmasked.fa)
for Query in $(ls repeat_masked/*/*/*/*_contigs_unmasked.fa | grep 'FON_63'); do
Strain=$(echo $Query | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Query | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
Prefix="$Strain"_vs_4287
OutDir=analysis/genome_alignment/mummer/$Organism/$Strain/$Prefix
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/MUMmer
qsub $ProgDir/sub_nucmer.sh $Reference $Query $Prefix $OutDir
done
```
