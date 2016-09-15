
```bash
  # for Strain in 125 A23 Fus2_edited_v2 55 A1-2 CB3 HB6 A13 A28 D2 PG fo47 4287; do
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/domain_analysis
for SummaryTab in $(ls gene_pred/annotations/*/*/*_gene_annotations.tab); do
Organism=$(echo $SummaryTab | rev | cut -f3 -d '/' | rev)
Strain=$(echo $SummaryTab | rev | cut -f2 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=analysis/NLR/$Organism/$Strain
mkdir -p $OutDir
StandTxt=$OutDir/"$Strain"_stand.txt
$ProgDir/extract_stand_domains.py --summary_table $SummaryTab > $StandTxt
done
```

```bash
  for File in $(ls analysis/NLR/*/*/*_stand.txt); do
    echo $File;
    cat $File | grep 'FoL contig' | sort -n -t' ' -k'3';
    echo "";  
  done
```


```bash
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/domain_analysis
SummaryTab=gene_pred/annotations/F.oxysporum_fsp_lycopersici/4287/4287_gene_annotations.tab
Organism=$(echo $SummaryTab | rev | cut -f3 -d '/' | rev);
Strain=$(echo $SummaryTab | rev | cut -f2 -d '/' | rev);
echo "$Organism - $Strain";
OutDir=analysis/NLR/$Organism/$Strain
StandTxt=$OutDir/"$Strain"_stand.txt
FoL_genome=assembly/external_group/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl/Fusarium_oxysporum_chromosome_and_additional_contigs.fa
$ProgDir/extract_stand_domains_FoL.py --FoL_genome $FoL_genome --summary_table $SummaryTab --N_term_out "$StandTxt".Nterm --NBD_out "$StandTxt".NBD --C_term_out "$StandTxt".Cterm > $StandTxt
```

```bash
  for File in $(ls tmp4/*_stand_FoL_specific.txt); do
    echo $File;
    cat $File | grep 'FoL contig' | sort -n -t' ' -k'3';
    echo "";  
  done
```
