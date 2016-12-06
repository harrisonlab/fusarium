This document details the commands used to analyse publically availabe genomes
downloaded onto the cluster.

This analysis was performed within the fusarium project directory:
/home/groups/harrisonlab/project_files/fusarium

<!--
but within a subdirectory named downloaded genomes. In essence, it made a second
project folder within the the first project folder, and therefore has its own set
of genes etc.

```bash
  mkdir -p /home/groups/harrisonlab/project_files/fusarium/downloaded_genomes
  cd /home/groups/harrisonlab/project_files/fusarium/downloaded_genomes
  mkdir -p assembly
  cp -r ../assembly/external_group assembly/.
```

Cegma analysis was run on the published genome assemblies to identify
orthologous core genes.

```bash
  WorkDir=/home/groups/harrisonlab/project_files/fusarium/downloaded_genomes
  for Genome in $(ls -d assembly/external_group/F.*/*); do
    cd $Genome;
    echo "$Genome";
    tar -xvf *.tar;
    cd $WorkDir;
  done
```

Each of these folders was confirmed to contain a set of contigs named supercontigs:

```bash
  ls -d assembly/external_group/F.*/*/ | wc -l
  ls assembly/external_group/F.*/*/*supercontigs* | grep -v 'linear'| wc -l
```

cegma was performed on assemblies to identify core eukaryotic genes within the
assemblies

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/cegma
  for Assembly in $(ls assembly/external_group/F.*/*/*supercontigs* | grep -v 'linear'); do
    echo $Assembly;
    qsub $ProgDir/sub_cegma.sh $Assembly dna;
  done
```

Outputs were summarised using the commands:
```bash
for File in $(ls gene_pred/cegma/F*/*/*_dna_cegma.completeness_report | grep -v 'cepae'); do
Strain=$(echo $File | rev | cut -f2 -d '/' | rev);
Species=$(echo $File | rev | cut -f3 -d '/' | rev);
printf "$Species\t$Strain\n";
cat $File | head -n18 | tail -n+4;printf "\n";
done >> gene_pred/cegma/cegma_results_dna_summary.txt

less gene_pred/cegma/cegma_results_dna_summary.txt
```

```
-->

```bash
  Fo_Fo47_assembly=assembly/external_group/F.oxysporum/fo47/broad/fusarium_oxysporum_fo47_1_supercontigs.fasta
  Fo_Fo47_genes=assembly/external_group/F.oxysporum/fo47/broad/fusarium_oxysporum_fo47_1_proteins.fasta
  Fo_Fo47_gff=assembly/external_group/F.oxysporum/fo47/broad/fusarium_oxysporum_fo47_1_transcripts.gtf

  FoL_4287_assembly=assembly/external_group/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl/Fusarium_oxysporum.FO2.31.dna.chromosome.fa
  FoL_4287_genes=assembly/external_group/F.oxysporum_fsp_lycopersici/4287/Fusox1/Fusox1_GeneCatalog_proteins_20110522.aa.fasta
  FoL_4287_gff=assembly/external_group/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl/Fusarium_oxysporum.FO2.31.gff3
```

an addtional FoL assembly was downloaded
```bash
CurDir=$PWD
OutDir=assembly/external_group/F.oxysporum_fsp_lycopersici/4287_v2/fungidb
mkdir -p $OutDir
cd $OutDir

wget http://fungidb.org/common/downloads/Current_Release/Foxysporum4287/fasta/data/FungiDB-29_Foxysporum4287_AnnotatedCDSs.fasta
wget http://fungidb.org/common/downloads/Current_Release/Foxysporum4287/fasta/data/FungiDB-29_Foxysporum4287_AnnotatedProteins.fasta
wget http://fungidb.org/common/downloads/Current_Release/Foxysporum4287/fasta/data/FungiDB-29_Foxysporum4287_AnnotatedTranscripts.fasta
wget http://fungidb.org/common/downloads/Current_Release/Foxysporum4287/fasta/data/FungiDB-29_Foxysporum4287_Genome.fasta
wget http://fungidb.org/common/downloads/Current_Release/Foxysporum4287/fasta/data/FungiDB-29_Foxysporum4287_ORFs_AA.fasta
wget http://fungidb.org/common/downloads/Current_Release/Foxysporum4287/gff/data/FungiDB-29_Foxysporum4287.gff
wget http://fungidb.org/common/downloads/Current_Release/Foxysporum4287/txt/FungiDB-29_Foxysporum4287_CodonUsage.txt
wget http://fungidb.org/common/downloads/Current_Release/Foxysporum4287/txt/FungiDB-29_Foxysporum4287_InterproDomains.txt
cd $CurDir
```

FoL assembly and protein sequences were parsed

```bash
  OutDir=assembly/external_group/F.oxysporum_fsp_lycopersici/4287_v2/fungidb
  # genome
  FoL_4287_assembly=$OutDir/FungiDB-29_Foxysporum4287_Genome.fasta
  FoL_4287_assembly_parsed=$OutDir/FungiDB-29_Foxysporum4287_Genome_parsed.fasta
  cat $FoL_4287_assembly | cut -f1 -d ' ' > $FoL_4287_assembly_parsed
  # gene models
  FoL_4287_genes=$OutDir/FungiDB-29_Foxysporum4287_AnnotatedProteins.fasta
  FoL_4287_genes_parsed=$OutDir/FungiDB-29_Foxysporum4287_AnnotatedProteins_parsed.fasta
  cat $FoL_4287_genes | sed -E 's/>*transcript=/>/g' | sed 's/ |.*//g' > $FoL_4287_genes_parsed
  # Gff of gene models
  FoL_4287_gff=$OutDir/FungiDB-29_Foxysporum4287.gff  
  FoL_4287_gff_parsed=$OutDir/FungiDB-29_Foxysporum4287_parsed.gff  
  cat $FoL_4287_gff | sed -r "s/\ttranscript\t/\tmRNA\t/g" > $FoL_4287_gff_parsed
```

<!-- Additional FoL contigs have been assembled that are not included in the 4287
chromosomal scaffolds. These were included in the assembly.

```bash
  AdditionalContigs=assembly/external_group/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl/Fusarium_oxysporum.FO2.31.dna_sm.nonchromosomal.fa
  AdditionalContigsParsed=assembly/external_group/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl/Fusarium_oxysporum.FO2.31.dna_sm.nonchromosomal_parsed.fa
  FoL_4287_assembly_parsed=assembly/external_group/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl/Fusarium_oxysporum.FO2.31.dna.chromosome_parsed.fa
  CompleteAssembly=assembly/external_group/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl/Fusarium_oxysporum_chromosome_and_additional_contigs.fa
  cat $AdditionalContigs | cut -f1 -d ' ' > $AdditionalContigsParsed
  cat $FoL_4287_assembly_parsed $AdditionalContigsParsed > $CompleteAssembly
``` -->

Fo Fo47 assembly and protein sequences were parsed

```bash
  # genome
  Fo_Fo47_assembly_parsed=assembly/external_group/F.oxysporum/fo47/broad/fusarium_oxysporum_fo47_1_supercontigs_parsed.fasta
  cat $Fo_Fo47_assembly | cut -f1 -d ' ' > $Fo_Fo47_assembly_parsed
  # gene models
  mkdir -p gene_pred/external_group/F.oxysporum/fo47/Fusox1
  Fo_Fo47_genes_parsed=gene_pred/external_group/F.oxysporum/fo47/Fusox1/fusarium_oxysporum_fo47_1_proteins_parsed.fa
  cat $Fo_Fo47_genes | cut -f1 -d ' ' > $Fo_Fo47_genes_parsed
  # Gff gene models
  Fo_Fo47_gff_parsed=assembly/external_group/F.oxysporum/fo47/broad/fusarium_oxysporum_fo47_1_transcripts_parsed.gff3
  ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/tools
  $ProgDir/fo47_gtf2_gff.py --gtf $Fo_Fo47_gff > $Fo_Fo47_gff_parsed
```


Quast

Quast was run to collect assembly statistics

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls $Fo_Fo47_assembly_parsed $FoL_4287_assembly_parsed); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    OutDir=$(dirname $Assembly)
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```


## Repeatmasking


Repeat masking was performed and used the following programs:
	Repeatmasker
	Repeatmodeler

The best assemblies were used to perform repeatmasking

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking
  for BestAss in $(ls assembly/external_group/F.oxysporum_fsp_lycopersici/4287_v2/fungidb/FungiDB-29_Foxysporum4287_Genome_parsed.fasta); do
    OutDir=repeat_masked/$Organism/"$Strain"/fungidb_repmask
    qsub $ProgDir/rep_modeling.sh $BestAss $OutDir
    qsub $ProgDir/transposonPSI.sh $BestAss $OutDir
  done
```


The TransposonPSI masked bases were used to mask additional bases from the
repeatmasker / repeatmodeller softmasked and harmasked files.

```bash
for File in $(ls repeat_masked/*/*/*/*_contigs_softmasked.fa | grep '4287_v2'); do
OutDir=$(dirname $File)
TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
OutFile=$(echo $File | sed 's/_contigs_softmasked.fa/_contigs_softmasked_repeatmasker_TPSI_appended.fa/g')
echo "$OutFile"
bedtools maskfasta -soft -fi $File -bed $TPSI -fo $OutFile
echo "Number of masked bases:"
cat $OutFile | grep -v '>' | tr -d '\n' | awk '{print $0, gsub("[a-z]", ".")}' | cut -f2 -d ' '
done
# The number of N's in hardmasked sequence are not counted as some may be present within the assembly and were therefore not repeatmasked.
for File in $(ls repeat_masked/*/*/*/*_contigs_hardmasked.fa | grep '4287_v2'); do
OutDir=$(dirname $File)
TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
OutFile=$(echo $File | sed 's/_contigs_hardmasked.fa/_contigs_hardmasked_repeatmasker_TPSI_appended.fa/g')
echo "$OutFile"
bedtools maskfasta -fi $File -bed $TPSI -fo $OutFile
echo "Number of masked bases:"
# cat $OutFile | grep -v '>' | tr -d '\n' | awk '{print $0, gsub("[N]", ".")}' | cut -f2 -d ' '
done
```

The number of bases masked by transposonPSI and Repeatmasker were summarised
using the following commands:

```bash
  for RepDir in $(ls -d repeat_masked/F.*/*/* | grep -e 'fo47' -e '4287_v2'); do
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

# Assessing gene space in assemblies:

- Quality of genome assemblies were assessed using Cegma to see how many core eukaryotic genes can be identified.


```bash
  Fo_Fo47_assembly_parsed=assembly/external_group/F.oxysporum/fo47/broad/fusarium_oxysporum_fo47_1_supercontigs_parsed.fasta
  FoL_4287_assembly_parsed=assembly/external_group/F.oxysporum_fsp_lycopersici/4287_v2/fungidb/FungiDB-29_Foxysporum4287_Genome_parsed.fasta
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/cegma
  cd /home/groups/harrisonlab/project_files/fusarium
  for Genome in $(ls $Fo_Fo47_assembly_parsed $FoL_4287_assembly_parsed); do
    echo $Genome;
    qsub $ProgDir/sub_cegma.sh $Genome dna;
  done
```

Outputs were summarised using the commands:
```bash
	for File in $(ls gene_pred/cegma/F*/FOP1/*_dna_cegma.completeness_report); do
		Strain=$(echo $File | rev | cut -f2 -d '/' | rev);
		Species=$(echo $File | rev | cut -f3 -d '/' | rev);
		printf "$Species\t$Strain\n";
		cat $File | head -n18 | tail -n+4;printf "\n";
	done > gene_pred/cegma/cegma_results_dna_summary.txt

	less gene_pred/cegma/cegma_results_dna_summary.txt
```

```bash
FoL_4287_genes_parsed=assembly/external_group/F.oxysporum_fsp_lycopersici/4287_v2/fungidb/FungiDB-29_Foxysporum4287_AnnotatedProteins_parsed.fasta
Fo_Fo47_genes_parsed=gene_pred/external_group/F.oxysporum/fo47/Fusox1/fusarium_oxysporum_fo47_1_proteins_parsed.fa
for Proteome in $(ls $FoL_4287_genes_parsed $Fo_Fo47_genes_parsed); do
  Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
  echo "$Organism - $Strain"
    cat $Proteome | grep '>' | wc -l
  if [ $Strain == 'fo47' ]; then
    cat $Proteome | grep '>' | cut -f1 -d 'T' | sort | uniq | wc -l
  elif [ $Strain == '4287_v2' ]; then
    cat $Proteome | grep '>' | cut -f1 -d '-' | sort | uniq | wc -l
  fi
done
```

### A) From Predicted gene models - Identifying secreted proteins

Required programs:
 * SignalP-4.1
 * TMHMM

Proteins that were predicted to contain signal peptides were identified using
the following commands:

```bash
SplitfileDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
CurPath=$PWD
FoL_4287_genes_parsed=assembly/external_group/F.oxysporum_fsp_lycopersici/4287_v2/fungidb/FungiDB-29_Foxysporum4287_AnnotatedProteins_parsed.fasta
# Fo_Fo47_genes_parsed=gene_pred/external_group/F.oxysporum/fo47/Fusox1/fusarium_oxysporum_fo47_1_proteins_parsed.fa
for Proteome in $(ls $FoL_4287_genes_parsed $Fo_Fo47_genes_parsed); do
Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
SplitDir=gene_pred/final_genes_split/$Organism/$Strain
mkdir -p $SplitDir
BaseName="$Organism""_$Strain"_final_preds
$SplitfileDir/splitfile_500.py --inp_fasta $Proteome --out_dir $SplitDir --out_base $BaseName
for File in $(ls $SplitDir/*_final_preds_*); do
Jobs=$(qstat | grep 'pred_sigP' | wc -l)
while [ $Jobs -gt 20 ]; do
sleep 10
printf "."
Jobs=$(qstat | grep 'pred_sigP' | wc -l)
done
printf "\n"
echo $File
qsub $ProgDir/pred_sigP.sh $File signalp-4.1
done
done
```

The batch files of predicted secreted proteins needed to be combined into a
single file for each strain. This was done with the following commands:
```bash
  for SplitDir in $(ls -d gene_pred/final_genes_split/*/* | grep -w -e '4287_v2' -e 'fo47'); do
    Strain=$(echo $SplitDir | rev |cut -d '/' -f1 | rev)
    Organism=$(echo $SplitDir | rev |cut -d '/' -f2 | rev)
    InStringAA=''
    InStringNeg=''
    InStringTab=''
    InStringTxt=''
    SigpDir=final_genes_signalp-4.1
    for GRP in $(ls -l $SplitDir/*_final_preds_*.fa | rev | cut -d '_' -f1 | rev | sort -n); do  
    InStringAA="$InStringAA gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_preds_$GRP""_sp.aa";  
    InStringNeg="$InStringNeg gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_preds_$GRP""_sp_neg.aa";  
    InStringTab="$InStringTab gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_preds_$GRP""_sp.tab";
    InStringTxt="$InStringTxt gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_preds_$GRP""_sp.txt";  
    done
    cat $InStringAA > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_final_sp.aa
    cat $InStringNeg > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_final_neg_sp.aa
    tail -n +2 -q $InStringTab > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_final_sp.tab
    cat $InStringTxt > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_final_sp.txt
  done
```

Some proteins that are incorporated into the cell membrane require secretion.
Therefore proteins with a transmembrane domain are not likely to represent
cytoplasmic or apoplastic effectors.

Proteins containing a transmembrane domain were identified:

```bash
  FoL_4287_genes_parsed=assembly/external_group/F.oxysporum_fsp_lycopersici/4287_v2/fungidb/FungiDB-29_Foxysporum4287_AnnotatedProteins_parsed.fasta
  Fo_Fo47_genes_parsed=gene_pred/external_group/F.oxysporum/fo47/Fusox1/fusarium_oxysporum_fo47_1_proteins_parsed.fa
  for Proteome in $(ls $FoL_4287_genes_parsed $Fo_Fo47_genes_parsed); do
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/transmembrane_helices
    qsub $ProgDir/submit_TMHMM.sh $Proteome
  done
```

Those proteins with transmembrane domains were removed from lists of Signal
peptide containing proteins

```bash
for File in $(ls gene_pred/trans_mem/*/*/*_TM_genes_neg.txt | grep -e 'fo47' -e '4287_v2'); do
Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
TmHeaders=$(echo "$File" | sed 's/neg.txt/neg_headers.txt/g')
cat $File | cut -f1 > $TmHeaders
SigP=$(ls gene_pred/final_genes_signalp-4.1/$Organism/$Strain/*_final_sp.aa)
cat $SigP | grep '>' | wc -l
OutDir=$(dirname $SigP)
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $SigP --headers $TmHeaders > $OutDir/"$Strain"_final_sp_no_trans_mem.aa
cat $OutDir/"$Strain"_final_sp_no_trans_mem.aa | grep '>' | wc -l
if [ $Strain == 'fo47' ]; then
  cat $OutDir/"$Strain"_final_sp_no_trans_mem.aa | grep '>' | cut -f1 -d 'T' | sort | uniq | wc -l
elif [ $Strain == '4287_v2' ]; then
  cat $OutDir/"$Strain"_final_sp_no_trans_mem.aa | grep '>' | cut -f1 -d '-' | sort | uniq | wc -l
fi
done
```
```
F.oxysporum - fo47
1933
1540
F.oxysporum_fsp_lycopersici - 4287_v2
2053
1638
```



### B) From Augustus gene models - Effector identification using EffectorP

Required programs:
 * EffectorP.py

```bash
  FoL_4287_genes_parsed=assembly/external_group/F.oxysporum_fsp_lycopersici/4287_v2/fungidb/FungiDB-29_Foxysporum4287_AnnotatedProteins_parsed.fasta
  Fo_Fo47_genes_parsed=gene_pred/external_group/F.oxysporum/fo47/Fusox1/fusarium_oxysporum_fo47_1_proteins_parsed.fa
  for Proteome in $(ls $FoL_4287_genes_parsed $Fo_Fo47_genes_parsed); do
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    BaseName="$Organism"_"$Strain"_EffectorP
    OutDir=analysis/effectorP/$Organism/$Strain
    ProgDir=~/git_repos/emr_repos/tools/seq_tools/feature_annotation/fungal_effectors
    qsub $ProgDir/pred_effectorP.sh $Proteome $BaseName $OutDir
  done
```

Those genes that were predicted as secreted and tested positive by effectorP
were identified:

```bash
for File in $(ls analysis/effectorP/*/*/*_EffectorP.txt | grep -e 'fo47' -e '4287_v2' ); do
Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
Headers=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_headers.txt/g')
cat $File | grep 'Effector' | cut -f1 > $Headers
Secretome=$(ls gene_pred/final_genes_signalp-4.1/$Organism/$Strain/*_final_sp_no_trans_mem.aa)
OutFile=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_secreted.aa/g')
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $Secretome --headers $Headers > $OutFile
OutFileHeaders=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_secreted_headers.txt/g')
cat $OutFile | grep '>' | tr -d '>' | cut -f1 > $OutFileHeaders
cat $OutFileHeaders | wc -l
if [ $Strain == 'fo47' ]; then
  cat $OutFileHeaders | cut -f1 -d 'T' | sort | uniq | wc -l
elif [ $Strain == '4287_v2' ]; then
  cat $OutFileHeaders | cut -f1 -d '-' | sort | uniq | wc -l
fi
# Gff=$(ls gene_pred/final_genes/$Organism/$Strain/*/final_genes_appended.gff3)
# EffectorP_Gff=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_secreted.gff/g')
# ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
# $ProgDir/extract_gff_for_sigP_hits.pl $OutFileHeaders $Gff effectorP ID > $EffectorP_Gff
done
```
```
F.oxysporum - fo47
319
F.oxysporum_fsp_lycopersici - 4287
384
```

### C) Identification of MIMP-flanking genes

```bash
Fo_Fo47_assembly_parsed=assembly/external_group/F.oxysporum/fo47/broad/fusarium_oxysporum_fo47_1_supercontigs_parsed.fasta
Fo_Fo47_genes=assembly/external_group/F.oxysporum/fo47/broad/fusarium_oxysporum_fo47_1_proteins.fasta
Fo_Fo47_gff=assembly/external_group/F.oxysporum/fo47/broad/fusarium_oxysporum_fo47_1_transcripts.gtf

# FoL_4287_assembly_parsed=assembly/external_group/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl/Fusarium_oxysporum.FO2.31.dna.chromosome_parsed.fa
# FoL_4287_assembly_parsed=assembly/external_group/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl/Fusarium_oxysporum_chromosome_and_additional_contigs.fa
# FoL_4287_assembly_parsed=assembly/external_group/F.oxysporum_fsp_lycopersici/4287/Ma_et_al_2010/F.oxysporum_fsp.lycopersici_4287.fasta
# FoL_4287_assembly_parsed=assembly/external_group/F.oxysporum_fsp_lycopersici/4287/Ma_et_al_2010/F.oxysporum_fsp.lycopersici_4287.fasta
  FoL_4287_assembly_parsed=assembly/external_group/F.oxysporum_fsp_lycopersici/4287_v2/fungidb/FungiDB-29_Foxysporum4287_Genome_parsed.fasta
  FoL_4287_genes_parsed=assembly/external_group/F.oxysporum_fsp_lycopersici/4287_v2/fungidb/FungiDB-29_Foxysporum4287_AnnotatedProteins_parsed.fasta
  FoL_4287_gff_parsed=assembly/external_group/F.oxysporum_fsp_lycopersici/4287_v2/fungidb/FungiDB-29_Foxysporum4287_parsed.gff  
# FoL_4287_genes=assembly/external_group/F.oxysporum_fsp_lycopersici/4287/Fusox1/Fusox1_GeneCatalog_proteins_20110522.aa.fasta
# FoL_4287_gff=assembly/external_group/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl/Fusarium_oxysporum.FO2.31.gff3
for Genome in $(ls $Fo_Fo47_assembly_parsed $FoL_4287_assembly_parsed ); do
Organism=$(echo "$Genome" | rev | cut -d '/' -f4 | rev)
Strain=$(echo "$Genome" | rev | cut -d '/' -f3 | rev)
if [ $Strain == 'fo47' ]; then
  BrakerGff=$(ls assembly/external_group/F.oxysporum/fo47/broad/fusarium_oxysporum_fo47_1_transcripts.gtf)
elif [ $Strain == '4287_v2' ]; then
  BrakerGff=$(ls $FoL_4287_gff_parsed)
fi
OutDir=analysis/mimps/$Organism/$Strain
mkdir -p "$OutDir"
echo "$Organism - $Strain"
ProgDir="/home/armita/git_repos/emr_repos/tools/pathogen/mimp_finder"
$ProgDir/mimp_finder.pl $Genome $OutDir/"$Strain"_mimps.fa $OutDir/"$Strain"_mimps.gff > $OutDir/"$Strain"_mimps.log
$ProgDir/gffexpander.pl +- 2000 $OutDir/"$Strain"_mimps.gff > $OutDir/"$Strain"_mimps_exp.gff
echo "The number of mimps identified:"
cat $OutDir/"$Strain"_mimps.fa | grep '>' | wc -l
bedtools intersect -u -a $BrakerGff -b $OutDir/"$Strain"_mimps_exp.gff > $OutDir/"$Strain"_genes_in_2kb_mimp.gff
echo "The following transcripts intersect mimps:"
MimpProtsTxt=$OutDir/"$Strain"_prots_in_2kb_mimp.txt
MimpGenesTxt=$OutDir/"$Strain"_genes_in_2kb_mimp.txt
if [ $Strain == 'fo47' ]; then
  cat $OutDir/"$Strain"_genes_in_2kb_mimp.gff | grep -w 'exon' | cut -f9 | cut -f4 -d'"' | sort | uniq > $MimpProtsTxt
  cat $OutDir/"$Strain"_genes_in_2kb_mimp.gff | grep -w 'exon' | cut -f9 | cut -f4 -d'"' | cut -f1 -d 'T' | sort | uniq > $MimpGenesTxt
elif [ $Strain == '4287_v2' ]; then
  cat $OutDir/"$Strain"_genes_in_2kb_mimp.gff | grep -w 'mRNA' | cut -f9 | cut -f1 -d';' | cut -f2 -d'=' | sort | uniq > $MimpProtsTxt
  cat $OutDir/"$Strain"_genes_in_2kb_mimp.gff | grep -w 'mRNA' | cut -f9 | cut -f1 -d';' | cut -f2 -d'=' | cut -f1 -d '-' | sort | uniq > $MimpGenesTxt
fi
cat $MimpProtsTxt | wc -l
cat $MimpGenesTxt | wc -l
echo ""
done
```

```
The number of mimps identified:
158
The following transcripts intersect mimps:
134
```

Those genes that were predicted as secreted and within 2Kb of a MIMP
were identified:

```bash
for File in $(ls analysis/mimps/*/*/*_genes_in_2kb_mimp.txt | grep -e 'fo47' -e '4287_v2' -e 'ncbi' -e 'Fus2_canu_new' | grep -e 'fo47' -e '4287_v2'); do
  ProtsFile=$(echo $File | sed 's/genes/prots/g')
Strain=$(echo $File | rev | cut -f2 -d '/' | rev | sed 's/_chromosomal//g')
Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
Secretome=$(ls gene_pred/final_genes_signalp-4.1/$Organism/$Strain/*_final_sp_no_trans_mem.aa)
OutFile=$(echo "$File" | sed 's/.gff/_secreted.gff/g')
SecretedHeaders=$(echo "$Secretome" | sed 's/.aa/_headers.txt/g')
cat $Secretome | grep '>' | tr -d '>' | sed 's/-p.//g' > $SecretedHeaders
# cat $Secretome | grep '>' | tr -d '>' | sed 's/-p.//g' > $SecretedUniqHeaders
cat $ProtsFile $SecretedHeaders | cut -f1 | sort | uniq -d | wc -l
if [ $Strain == 'fo47' ]; then
cat $SecretedHeaders | cut -f1 | cut -f1 -d 'T' | sort | uniq | grep -f $File | wc -l
elif [ $Strain == '4287_v2' ]; then
cat $SecretedHeaders | cut -f1 | cut -f1 -d '-' | sort | uniq | grep -f $File | wc -l
fi
# ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
# $ProgDir/extract_gff_for_sigP_hits.pl $SecretedHeaders $File secreted_mimp ID > $OutFile
# cat $OutFile | grep -w 'mRNA' | wc -l
done
```

F.oxysporum - fo47
3
F.oxysporum_fsp_lycopersici - 4287_v2
23

## A) Interproscan

Interproscan was used to give gene models functional annotations.
Annotation was run using the commands below:

Note: This is a long-running script. As such, these commands were run using
'screen' to allow jobs to be submitted and monitored in the background.
This allows the session to be disconnected and reconnected over time.

Screen ouput detailing the progress of submission of interporscan jobs
was redirected to a temporary output file named interproscan_submission.log .

```bash
  screen -a
  cd /home/groups/harrisonlab/project_files/fusarium
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
  FoL_4287_genes_parsed=assembly/external_group/F.oxysporum_fsp_lycopersici/4287_v2/fungidb/FungiDB-29_Foxysporum4287_AnnotatedProteins_parsed.fasta
  # Fo_Fo47_genes_parsed=gene_pred/external_group/F.oxysporum/fo47/Fusox1/fusarium_oxysporum_fo47_1_proteins_parsed.fa
  for Genes in $(ls $FoL_4287_genes_parsed $Fo_Fo47_genes_parsed); do
    echo $Genes
    $ProgDir/sub_interproscan.sh $Genes
  done 2>&1 | tee -a interproscan_submisison.log
```

Following interproscan annotation split files were combined using the following
commands:

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
  FoL_4287_genes_parsed=assembly/external_group/F.oxysporum_fsp_lycopersici/4287_v2/fungidb/FungiDB-29_Foxysporum4287_AnnotatedProteins_parsed.fasta
  Fo_Fo47_genes_parsed=gene_pred/external_group/F.oxysporum/fo47/Fusox1/fusarium_oxysporum_fo47_1_proteins_parsed.fa
  for Proteome in $(ls $FoL_4287_genes_parsed $Fo_Fo47_genes_parsed); do
    Strain=$(echo $Proteome | rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Proteome | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    echo $Strain
    InterProRaw=gene_pred/interproscan/$Organism/$Strain/raw
    $ProgDir/append_interpro.sh $Proteome $InterProRaw
  done
```


## B) SwissProt
```bash
  FoL_4287_genes_parsed=assembly/external_group/F.oxysporum_fsp_lycopersici/4287_v2/fungidb/FungiDB-29_Foxysporum4287_AnnotatedProteins_parsed.fasta
  # Fo_Fo47_genes_parsed=gene_pred/external_group/F.oxysporum/fo47/Fusox1/fusarium_oxysporum_fo47_1_proteins_parsed.fa
  for Proteome in $(ls $FoL_4287_genes_parsed $Fo_Fo47_genes_parsed); do
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    OutDir=gene_pred/swissprot/$Organism/$Strain
    SwissDbDir=../../uniprot/swissprot
    SwissDbName=uniprot_sprot
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/swissprot
    qsub $ProgDir/sub_swissprot.sh $Proteome $OutDir $SwissDbDir $SwissDbName
  done
```

```bash
  for SwissTable in $(ls gene_pred/swissprot/*/*/swissprot_v2015_10_hits.tbl | grep '4287_v2'); do
  # SwissTable=gene_pred/swissprot/Fus2/swissprot_v2015_10_hits.tbl
    Strain=$(echo $SwissTable | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $SwissTable | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    OutTable=gene_pred/swissprot/$Organism/$Strain/swissprot_vJul2016_tophit_parsed.tbl
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/swissprot
    $ProgDir/swissprot_parser.py --blast_tbl $SwissTable --blast_db_fasta ../../uniprot/swissprot/uniprot_sprot.fasta > $OutTable
  done
```



## C) CAZY proteins

Carbohydrte active enzymes were idnetified using CAZYfollowing recomendations
at http://csbl.bmb.uga.edu/dbCAN/download/readme.txt :

```bash
  FoL_4287_genes_parsed=assembly/external_group/F.oxysporum_fsp_lycopersici/4287_v2/fungidb/FungiDB-29_Foxysporum4287_AnnotatedProteins_parsed.fasta
  # Fo_Fo47_genes_parsed=gene_pred/external_group/F.oxysporum/fo47/Fusox1/fusarium_oxysporum_fo47_1_proteins_parsed.fa
  for Proteome in $(ls $FoL_4287_genes_parsed $Fo_Fo47_genes_parsed); do
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    OutDir=gene_pred/CAZY/$Organism/$Strain
    mkdir -p $OutDir
    Prefix="$Strain"_CAZY
    CazyHmm=../../dbCAN/dbCAN-fam-HMMs.txt
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/HMMER
    qsub $ProgDir/sub_hmmscan.sh $CazyHmm $Proteome $Prefix $OutDir
  done
```

The Hmm parser was used to filter hits by an E-value of E1x10-5 or E 1x10-e3 if they had a hit over a length of X %.

Those proteins with a signal peptide were extracted from the list and gff files
representing these proteins made.

```bash
for File in $(ls gene_pred/CAZY/*/*/*CAZY.out.dm | grep -e 'fo47' -e '4287_v2' | grep 'fo47'); do
Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
OutDir=$(dirname $File)
echo "$Organism - $Strain"
ProgDir=/home/groups/harrisonlab/dbCAN
$ProgDir/hmmscan-parser.sh $OutDir/"$Strain"_CAZY.out.dm > $OutDir/"$Strain"_CAZY.out.dm.ps
CazyHeaders=$(echo $File | sed 's/.out.dm/_headers.txt/g')
cat $OutDir/"$Strain"_CAZY.out.dm.ps | cut -f3 | sort | uniq > $CazyHeaders
echo "number of CAZY genes identified:"
cat $CazyHeaders | wc -l
if [ $Strain == 'fo47' ]; then
cat $CazyHeaders | cut -f1 -d 'T' | sort | uniq | wc -l
elif [ $Strain == '4287_v2' ]; then
cat $CazyHeaders | cut -f1 -d '-' | sort | uniq | wc -l
fi
# if [ $Strain == 'fo47' ]; then
#   Gff=$(ls assembly/external_group/F.oxysporum/fo47/broad/fusarium_oxysporum_fo47_1_transcripts.gtf)
# elif [ $Strain == '4287_chromosomal' ]; then
#   Gff=$(ls assembly/external_group/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl/Fusarium_oxysporum.FO2.31.gff3)
# fi
# CazyGff=$OutDir/"$Strain"_CAZY.gff
# ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
# $ProgDir/extract_gff_for_sigP_hits.pl $CazyHeaders $Gff CAZyme ID > $CazyGff

SecretedProts=$(ls gene_pred/final_genes_signalp-4.1/$Organism/$Strain/"$Strain"_final_sp_no_trans_mem.aa)
SecretedHeaders=$(echo $SecretedProts | sed 's/.aa/_headers.txt/g')
cat $SecretedProts | grep '>' | tr -d '>' > $SecretedHeaders
# CazyGffSecreted=$OutDir/"$Strain"_CAZY_secreted.gff
# $ProgDir/extract_gff_for_sigP_hits.pl $SecretedHeaders $CazyGff Secreted_CAZyme ID > $CazyGffSecreted
echo "number of Secreted CAZY genes identified:"
# cat $CazyGffSecreted | grep -w 'gene' | cut -f9 | tr -d 'ID=' | wc -l
cat $CazyHeaders $SecretedHeaders | cut -f1 | sort | uniq -d > $OutDir/"$Strain"_CAZY_secreted_headers.txt
cat $OutDir/"$Strain"_CAZY_secreted_headers.txt | wc -l
if [ $Strain == 'fo47' ]; then
cat $OutDir/"$Strain"_CAZY_secreted_headers.txt | cut -f1 -d 'T' | sort | uniq | wc -l
elif [ $Strain == '4287_v2' ]; then
cat $OutDir/"$Strain"_CAZY_secreted_headers.txt | cut -f1 -d '-' | sort | uniq | wc -l
fi
done
```
```
  F.oxysporum - fo47
  number of CAZY genes identified:
  1108 (928)
  number of Secreted CAZY genes identified:
  411 (382)
  F.oxysporum_fsp_lycopersici - 4287_v2
  number of CAZY genes identified:
  1156 (977)
  number of Secreted CAZY genes identified:
  419 (386)
```
Note - the CAZY genes identified may need further filtering based on e value and
cuttoff length - see below:

Cols in yourfile.out.dm.ps:
1. Family HMM
2. HMM length
3. Query ID
4. Query length
5. E-value (how similar to the family HMM)
6. HMM start
7. HMM end
8. Query start
9. Query end
10. Coverage

* For fungi, use E-value < 1e-17 and coverage > 0.45

* The best threshold varies for different CAZyme classes (please see http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4132414/ for details). Basically to annotate GH proteins, one should use a very relax coverage cutoff or the sensitivity will be low (Supplementary Tables S4 and S9); (ii) to annotate CE families a very stringent E-value cutoff and coverage cutoff should be used; otherwise the precision will be very low due to a very high false positive rate (Supplementary Tables S5 and S10)

## D) AntiSMASH

Antismash was run to identify clusters of secondary metabolite genes within
the genome. Antismash was run using the weserver at:
http://antismash.secondarymetabolites.org


Results of web-annotation of gene clusters within the assembly were downloaded to
the following directories:

```bash
for Assembly in $(ls repeat_masked/*/*/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep -e '4287_v2' -e 'fo47'); do
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
OutDir=analysis/antismash/$Organism/$Strain
mkdir -p $OutDir
done
```

```bash
for Zip in $(ls analysis/antismash/*/*/*.zip | grep -e '4287_v2' -e 'fo47'); do
OutDir=$(dirname $Zip)
unzip -d $OutDir $Zip
done
```


```bash
for AntiSmash in $(ls analysis/antismash/*/*/*/*.final.gbk | grep 'fo47'); do
Organism=$(echo $AntiSmash | rev | cut -f4 -d '/' | rev)
Strain=$(echo $AntiSmash | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=analysis/antismash/$Organism/$Strain
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/secondary_metabolites
$ProgDir/antismash2gff.py --inp_antismash $AntiSmash > $OutDir/"$Strain"_secondary_metabolite_regions.gff
printf "Number of clusters detected:\t"
cat $OutDir/"$Strain"_secondary_metabolite_regions.gff | grep 'antismash_cluster' | wc -l
GeneGff=$(ls assembly/external_group/$Organism/$Strain/*/*.gtf | grep -e 'fungidb' -e 'broad' | grep -v -e 'FungiDB-29_Foxysporum4287_parsed' -e 'fo47_1_transcripts_parsed')
bedtools intersect -u -a $GeneGff -b $OutDir/"$Strain"_secondary_metabolite_regions.gff > $OutDir/metabolite_cluster_genes.gff
cat $OutDir/metabolite_cluster_genes.gff | grep -w 'exon' | cut -f9 | cut -f4 -d '"' | sort | uniq > $OutDir/metabolite_cluster_gene_headers.txt
printf "Number of predicted proteins in clusters:\t"
cat $OutDir/metabolite_cluster_gene_headers.txt | wc -l
printf "Number of predicted genes in clusters:\t"
cat $OutDir/metabolite_cluster_genes.gff |  grep -w 'exon' | cut -f9 | cut -f2 -d '"' | sort | uniq | wc -l
done

for AntiSmash in $(ls analysis/antismash/*/*/*/*.final.gbk | grep -e '4287_v2'); do
Organism=$(echo $AntiSmash | rev | cut -f4 -d '/' | rev)
Strain=$(echo $AntiSmash | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=analysis/antismash/$Organism/$Strain
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/secondary_metabolites
$ProgDir/antismash2gff.py --inp_antismash $AntiSmash > $OutDir/"$Strain"_secondary_metabolite_regions.gff
printf "Number of clusters detected:\t"
cat $OutDir/"$Strain"_secondary_metabolite_regions.gff | grep 'antismash_cluster' | wc -l
GeneGff=$(ls assembly/external_group/F.oxysporum_fsp_lycopersici/4287_v2/fungidb/FungiDB-29_Foxysporum4287_parsed.gff)
bedtools intersect -u -a $GeneGff -b $OutDir/"$Strain"_secondary_metabolite_regions.gff > $OutDir/metabolite_cluster_genes.gff
cat $OutDir/metabolite_cluster_genes.gff | grep -w 'mRNA' | cut -f9 | cut -f2 -d '=' | cut -f1 -d ';' | sort | uniq > $OutDir/metabolite_cluster_gene_headers.txt
printf "Number of predicted proteins in clusters:\t"
cat $OutDir/metabolite_cluster_gene_headers.txt | wc -l
printf "Number of predicted genes in clusters:\t"
cat $OutDir/metabolite_cluster_genes.gff |  grep -w 'mRNA' | cut -f9 | cut -f3 -d '=' | sort | uniq | wc -l
done
```

These clusters represenyed the following genes. Note that these numbers just
show the number of intersected genes with gff clusters and are not confirmed by
function

```
F.oxysporum - fo47
Number of clusters detected:	47
Number of predicted genes in clusters:	705
F.oxysporum_fsp_lycopersici - 4287_v2
Number of clusters detected:	49
Number of predicted genes in clusters:	861
```

#Genomic analysis

## Comparison to FoL 4287

BLast searches were used to identify which genes had homologs on which
chromosomes of the Fusarium lycopersici genome.

```bash
  FoLGenomeFa=assembly/external_group/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl/Fusarium_oxysporum.FO2.31.dna.chromosome.fa
  for Proteome in $(ls assembly/external_group/F.oxysporum/fo47/broad/fusarium_oxysporum_fo47_1_proteins.fasta gene_pred/external_group/F.oxysporum_fsp_lycopersici/4287/Fusox1/Fusox1_GeneCatalog_proteins_20110522_parsed.fa); do
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=analysis/blast_homology/$Organism/$Strain
    ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
    qsub $ProgDir/run_blast2csv.sh $Proteome protein $FoLGenomeFa $OutDir
  done
```

Convert top blast hits into gff annotations


```bash
for BlastHitsCsv in $(ls analysis/blast_homology/F.oxysporum/fo47/4287_chromosomal_fusarium_oxysporum_fo47_1_proteins.fasta_hits.csv analysis/blast_homology/F.oxysporum_fsp_lycopersici/4287/4287_chromosomal_Fusox1_GeneCatalog_proteins_20110522_parsed.fa_hits.csv); do
Organism=$(echo $BlastHitsCsv | rev | cut -f3 -d '/' | rev)
Strain=$(echo $BlastHitsCsv | rev | cut -f2 -d '/' | rev)
echo "$Organism - $Strain"
HitsGff=$(echo $BlastHitsCsv | sed  's/.csv/.gff/g')
Column2="$Strain"_gene_homolog
NumHits=1
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
$ProgDir/blast2gff.pl $Column2 $NumHits $BlastHitsCsv > $HitsGff
done
```

#### 2.6.3) Intersecting blast hits with genes from FoL

```bash
for HitsGff in $(ls analysis/blast_homology/F.oxysporum/fo47/4287_chromosomal_fusarium_oxysporum_fo47_1_proteins.fasta_hits.gff analysis/blast_homology/F.oxysporum_fsp_lycopersici/4287/4287_chromosomal_Fusox1_GeneCatalog_proteins_20110522_parsed.fa_hits.gff ); do
Organism=$(echo $HitsGff | rev | cut -f3 -d '/' | rev)
Strain=$(echo $HitsGff| rev | cut -f2 -d '/' | rev)
echo "$Organism - $Strain"
HitsDir=$(dirname $HitsGff)
FoLGenes=assembly/external_group/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl/Fusarium_oxysporum.FO2.31.chr.gff3
FoLIntersect=$HitsDir/4287_chromosomal_final_genes_combined_intersect.bed
bedtools intersect -wo -a $HitsGff -b $FoLGenes > $FoLIntersect
done
```



## Identifying FTF genes

Previously published FTF genes from Sanchez et al 2016 were blasted against
Fusarium genomes.

```bash
FoL_4287_assembly_parsed=assembly/external_group/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl/Fusarium_oxysporum.FO2.31.dna.chromosome_parsed.fa
Fo_Fo47_assembly_parsed=assembly/external_group/F.oxysporum/fo47/broad/fusarium_oxysporum_fo47_1_supercontigs_parsed.fasta
for Assembly in $(ls $FoL_4287_assembly_parsed $Fo_Fo47_assembly_parsed); do
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
echo $Assembly
Query=analysis/blast_homology/Fo_path_genes/FTF_cds_Sanchez_et_al_2016.fasta
OutDir=analysis/FTF/$Organism/$Strain
mkdir -p $OutDir
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
qsub $ProgDir/run_blast2csv.sh $Query dna $Assembly $OutDir
done
```

BLAST hits were converted to Gff annotations and intersected with gene models:

```bash
for BlastHits in $(ls analysis/FTF/*/*/*_FTF_cds_Sanchez_et_al_2016.fasta_hits.csv | grep -e '4287' -e 'fo47'); do
Strain=$(echo $BlastHits | rev | cut -f2 -d '/' | rev)
Organism=$(echo $BlastHits | rev | cut -f3 -d '/' | rev)
OutDir=analysis/FTF/$Organism/$Strain
HitsGff=$(echo $BlastHits | sed  's/.csv/.gff/g')
Column2=FTF_homolog
NumHits=1
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
$ProgDir/blast2gff.pl $Column2 $NumHits $BlastHits > $HitsGff

GffAppended=$(ls gene_pred/codingquary/$Organism/$Strain/final/final_genes_appended.gff3)
bedtools intersect -wao -a $HitsGff -b $GffAppended > $OutDir/"$Strain"_FTF_hits_intersected.bed
done
```
