
```bash
  # for Strain in 125 A23 Fus2_edited_v2 55 A1-2 CB3 HB6 A13 A28 D2 PG fo47 4287; do
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/domain_analysis
for SummaryTab in $(ls gene_pred/annotations/*/*/*_gene_annotations.tab); do
Organism=$(echo $SummaryTab | rev | cut -f3 -d '/' | rev)
Strain=$(echo $SummaryTab | rev | cut -f2 -d '/' | rev)
echo "$Organism - $Strain"
StandTxt=tmp4/"$Strain"_stand.txt
$ProgDir/extract_stand_domains.py --summary_table $SummaryTab > $StandTxt
done
```

```bash
  for File in $(ls tmp4/*_stand.txt); do
    echo $File;
    cat $File | grep 'FoL contig' | sort -n -t' ' -k'3';
    echo "";  
  done
```


```bash
  for SummaryTab in $(ls gene_pred/annotations/*/4287*/*_gene_annotations.tab); do
    Organism=$(echo $SummaryTab | rev | cut -f3 -d '/' | rev);
    Strain=$(echo $SummaryTab | rev | cut -f2 -d '/' | rev); echo "$Organism - $Strain"; StandTxt=tmp4/"$Strain"_stand_FoL_specific.txt;
    $ProgDir/extract_stand_domains_FoL.py --summary_table $SummaryTab > $StandTxt;
  done
```

```bash
  for File in $(ls tmp4/*_stand_FoL_specific.txt); do
    echo $File;
    cat $File | grep 'FoL contig' | sort -n -t' ' -k'3';
    echo "";  
  done
```
