# FON Na5 genome


```bash
  OutDir=assembly/external_group/F.oxysporum_fsp_narcissi/Na5/GCA_002233775_1
  mkdir -p $OutDir
  CurDir=$PWD
  cd $OutDir
  wget ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/NJ/CV/NJCV01/NJCV01.1.fsa_nt.gz
  gunzip *
  cd $CurDir
```


Quast

Quast was run to collect assembly statistics

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly/external_group/F.oxysporum_fsp_narcissi/Na5/GCA_002233775_1/NJCV01.1.fsa_nt); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    OutDir=$(dirname $Assembly)
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```

BUSCO

```bash
  for Assembly in $(ls assembly/external_group/F.oxysporum_fsp_narcissi/Na5/GCA_002233775_1/NJCV01.1.fsa_nt); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
    BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
    OutDir=gene_pred/busco/$Organism/$Strain/assembly
    qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
  done
```



## Repeatmasking


Repeat masking was performed and used the following programs:
	Repeatmasker
	Repeatmodeler

The best assemblies were used to perform repeatmasking

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking
  for Assembly in $(ls assembly/external_group/F.oxysporum_fsp_narcissi/Na5/GCA_002233775_1/NJCV01.1.fsa_nt); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    OutDir=repeat_masked/$Organism/"$Strain"/ncbi_repmask
    qsub $ProgDir/rep_modeling.sh $BestAss $OutDir
    qsub $ProgDir/transposonPSI.sh $BestAss $OutDir
  done
```


The TransposonPSI masked bases were used to mask additional bases from the
repeatmasker / repeatmodeller softmasked and harmasked files.

```bash
for File in $(ls repeat_masked/*/*/*/*_contigs_softmasked.fa | grep 'Na5'); do
OutDir=$(dirname $File)
TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
OutFile=$(echo $File | sed 's/_contigs_softmasked.fa/_contigs_softmasked_repeatmasker_TPSI_appended.fa/g')
echo "$OutFile"
bedtools maskfasta -soft -fi $File -bed $TPSI -fo $OutFile
echo "Number of masked bases:"
cat $OutFile | grep -v '>' | tr -d '\n' | awk '{print $0, gsub("[a-z]", ".")}' | cut -f2 -d ' '
done
# The number of N's in hardmasked sequence are not counted as some may be present within the assembly and were therefore not repeatmasked.
for File in $(ls repeat_masked/*/*/*/*_contigs_hardmasked.fa | grep -w -e '4287_v2' -e 'fo47' -e '7600' | grep '7600'); do
OutDir=$(dirname $File)
TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
OutFile=$(echo $File | sed 's/_contigs_hardmasked.fa/_contigs_hardmasked_repeatmasker_TPSI_appended.fa/g')
echo "$OutFile"
bedtools maskfasta -fi $File -bed $TPSI -fo $OutFile
# echo "Number of masked bases:"
# cat $OutFile | grep -v '>' | tr -d '\n' | awk '{print $0, gsub("[N]", ".")}' | cut -f2 -d ' '
done
```

The number of bases masked by transposonPSI and Repeatmasker were summarised
using the following commands:

```bash
  for RepDir in $(ls -d repeat_masked/F.*/*/*  | grep 'Na5'); do
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
  done
```



### C) Identification of MIMPs

```bash
# for Genome in $(ls repeat_masked/F.*/*/*/*_contigs_unmasked.fa | grep 'Na5'); do
for Genome in $(ls assembly/external_group/F.oxysporum_fsp_narcissi/Na5/GCA_002233775_1/NJCV01.1.fsa_nt); do
Organism=$(echo "$Genome" | rev | cut -d '/' -f4 | rev)
Strain=$(echo "$Genome" | rev | cut -d '/' -f3 | rev)
# BrakerGff=$(ls gene_pred/final_genes/$Organism/"$Strain"/final/final_genes_CodingQuary.gff3)
# QuaryGff=$(ls gene_pred/final_genes/$Organism/"$Strain"/final/final_genes_Braker.gff3)
OutDir=analysis/mimps/$Organism/$Strain
mkdir -p "$OutDir"
echo "$Organism - $Strain"
ProgDir="/home/armita/git_repos/emr_repos/tools/pathogen/mimp_finder"
$ProgDir/mimp_finder.pl $Genome $OutDir/"$Strain"_mimps.fa $OutDir/"$Strain"_mimps.gff > $OutDir/"$Strain"_mimps.log
$ProgDir/gffexpander.pl +- 2000 $OutDir/"$Strain"_mimps.gff > $OutDir/"$Strain"_mimps_exp.gff
echo "The number of mimps identified:"
cat $OutDir/"$Strain"_mimps.fa | grep '>' | wc -l
# bedtools intersect -u -a $BrakerGff -b $OutDir/"$Strain"_mimps_exp.gff > $OutDir/"$Strain"_genes_in_2kb_mimp.gff
# bedtools intersect -u -a $QuaryGff -b $OutDir/"$Strain"_mimps_exp.gff >> $OutDir/"$Strain"_genes_in_2kb_mimp.gff
# echo "The following transcripts intersect mimps:"
# MimpProtsTxt=$OutDir/"$Strain"_prots_in_2kb_mimp.txt
# MimpGenesTxt=$OutDir/"$Strain"_genes_in_2kb_mimp.txt
# cat $OutDir/"$Strain"_genes_in_2kb_mimp.gff | grep -w 'mRNA' | cut -f9 | cut -f1 -d';' | cut -f2 -d'=' | sort | uniq > $MimpProtsTxt
# cat $OutDir/"$Strain"_genes_in_2kb_mimp.gff | grep -w 'mRNA' | cut -f9 | cut -f1 -d';' | cut -f2 -d'=' | cut -f1 -d '.'| sort | uniq > $MimpGenesTxt
# cat $MimpProtsTxt | wc -l
# cat $MimpGenesTxt | wc -l
echo ""
done
```

```
  F.oxysporum_fsp_narcissi - Na5
  The number of mimps identified:
  154
```
