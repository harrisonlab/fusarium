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

FoL assembly and protein sequences were parsed

```bash
  # genome
  FoL_4287_assembly_parsed=assembly/external_group/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl/Fusarium_oxysporum.FO2.31.dna.chromosome_parsed.fa
  cat $FoL_4287_assembly | cut -f1 -d ' ' > $FoL_4287_assembly_parsed
  # gene models
  mkdir -p gene_pred/external_group/F.oxysporum_fsp_lycopersici/4287/Fusox1
  FoL_4287_genes_parsed=gene_pred/external_group/F.oxysporum_fsp_lycopersici/4287/Fusox1/Fusox1_GeneCatalog_proteins_20110522_parsed.fa
  cat $FoL_4287_genes | sed -E 's/>jgi.*FOXG/>FOXG/g' > $FoL_4287_genes_parsed
  # Gff of gene models
  FoL_4287_gff_parsed=assembly/external_group/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl/Fusarium_oxysporum.FO2.31_parsed.gff3
  cat $FoL_4287_gff | sed -r "s/\ttranscript\t/\tmRNA\t/g" | sed -r "s/=.\w*:/=/g" | grep 'Broad' | sed 's/T.;/;/g' > $FoL_4287_gff_parsed
```
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


### A) From Augustus gene models - Identifying secreted proteins

Required programs:
 * SignalP-4.1
 * TMHMM

Proteins that were predicted to contain signal peptides were identified using
the following commands:

```bash
  SplitfileDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
  CurPath=$PWD
  FoL_4287_genes_parsed=gene_pred/external_group/F.oxysporum_fsp_lycopersici/4287/Fusox1/Fusox1_GeneCatalog_proteins_20110522_parsed.fa
  Fo_Fo47_genes_parsed=gene_pred/external_group/F.oxysporum/fo47/Fusox1/fusarium_oxysporum_fo47_1_proteins_parsed.fa
  for Proteome in $(ls $FoL_4287_genes_parsed $Fo_Fo47_genes_parsed); do
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    SplitDir=gene_pred/final_genes_split/$Organism/$Strain
    mkdir -p $SplitDir
    BaseName="$Organism""_$Strain"_final_preds
    $SplitfileDir/splitfile_500.py --inp_fasta $Proteome --out_dir $SplitDir --out_base $BaseName
    for File in $(ls $SplitDir/*_final_preds_*); do
      Jobs=$(qstat | grep 'pred_sigP' | grep 'qw' | wc -l)
      while [ $Jobs -gt 1 ]; do
        sleep 10
        printf "."
        Jobs=$(qstat | grep 'pred_sigP' | grep 'qw' | wc -l)
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
for SplitDir in $(ls -d gene_pred/final_genes_split/*/* | grep -w -e '4287' -e 'fo47'); do
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
  FoL_4287_genes_parsed=gene_pred/external_group/F.oxysporum_fsp_lycopersici/4287/Fusox1/Fusox1_GeneCatalog_proteins_20110522_parsed.fa
  Fo_Fo47_genes_parsed=gene_pred/external_group/F.oxysporum/fo47/Fusox1/fusarium_oxysporum_fo47_1_proteins_parsed.fa
  for Proteome in $(ls $FoL_4287_genes_parsed $Fo_Fo47_genes_parsed); do
		Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/transmembrane_helices
		qsub $ProgDir/submit_TMHMM.sh $Proteome
	done
```



### B) From Augustus gene models - Effector identification using EffectorP

Required programs:
 * EffectorP.py

```bash
  FoL_4287_genes_parsed=gene_pred/external_group/F.oxysporum_fsp_lycopersici/4287/Fusox1/Fusox1_GeneCatalog_proteins_20110522_parsed.fa
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

### C) Identification of MIMP-flanking genes

```bash
  Fo_Fo47_assembly_parsed=assembly/external_group/F.oxysporum/fo47/broad/fusarium_oxysporum_fo47_1_supercontigs_parsed.fasta
  Fo_Fo47_genes=assembly/external_group/F.oxysporum/fo47/broad/fusarium_oxysporum_fo47_1_proteins.fasta
  Fo_Fo47_gff=assembly/external_group/F.oxysporum/fo47/broad/fusarium_oxysporum_fo47_1_transcripts.gtf

  FoL_4287_assembly_parsed=assembly/external_group/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl/Fusarium_oxysporum.FO2.31.dna.chromosome_parsed.fa
  FoL_4287_genes=assembly/external_group/F.oxysporum_fsp_lycopersici/4287/Fusox1/Fusox1_GeneCatalog_proteins_20110522.aa.fasta
  FoL_4287_gff=assembly/external_group/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl/Fusarium_oxysporum.FO2.31.gff3
  for Genome in $(ls $Fo_Fo47_assembly_parsed $FoL_4287_assembly_parsed); do
    Organism=$(echo "$Genome" | rev | cut -d '/' -f4 | rev)
    Strain=$(echo "$Genome" | rev | cut -d '/' -f3 | rev)
    if [ $Strain == 'fo47' ]; then
      BrakerGff=$(ls assembly/external_group/F.oxysporum/fo47/broad/fusarium_oxysporum_fo47_1_transcripts.gtf)
    elif [ $Strain == '4287_chromosomal' ]; then
      BrakerGff=$(ls assembly/external_group/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl/Fusarium_oxysporum.FO2.31.gff3)
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
    MimpGenesTxt=$OutDir/"$Strain"_genes_in_2kb_mimp.txt
    if [ $Strain == 'fo47' ]; then
      cat $OutDir/"$Strain"_genes_in_2kb_mimp.gff | grep -w 'exon' | cut -f9 | cut -f4 -d'"' > $MimpGenesTxt
    elif [ $Strain == '4287_chromosomal' ]; then
      cat $OutDir/"$Strain"_genes_in_2kb_mimp.gff | grep -w 'transcript' | cut -f9 | cut -f1 -d';' | cut -f2 -d':' > $MimpGenesTxt
    fi
    cat $MimpGenesTxt | wc -l
    echo ""
  done
```


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
  FoL_4287_genes_parsed=gene_pred/external_group/F.oxysporum_fsp_lycopersici/4287/Fusox1/Fusox1_GeneCatalog_proteins_20110522_parsed.fa
  Fo_Fo47_genes_parsed=gene_pred/external_group/F.oxysporum/fo47/Fusox1/fusarium_oxysporum_fo47_1_proteins_parsed.fa
  for Genes in $(ls $FoL_4287_genes_parsed $Fo_Fo47_genes_parsed); do
    echo $Genes
    $ProgDir/sub_interproscan.sh $Genes
  done 2>&1 | tee -a interproscan_submisison.log
```

Following interproscan annotation split files were combined using the following
commands:

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
  FoL_4287_genes_parsed=gene_pred/external_group/F.oxysporum_fsp_lycopersici/4287/Fusox1/Fusox1_GeneCatalog_proteins_20110522_parsed.fa
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
  FoL_4287_genes_parsed=gene_pred/external_group/F.oxysporum_fsp_lycopersici/4287/Fusox1/Fusox1_GeneCatalog_proteins_20110522_parsed.fa
  Fo_Fo47_genes_parsed=gene_pred/external_group/F.oxysporum/fo47/Fusox1/fusarium_oxysporum_fo47_1_proteins_parsed.fa
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
	for SwissTable in $(ls gene_pred/swissprot/*/*/swissprot_v2015_10_hits.tbl); do
		# SwissTable=gene_pred/swissprot/Fus2/swissprot_v2015_10_hits.tbl
		Strain=$(echo $SwissTable | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $SwissTable | rev | cut -f3 -d '/' | rev)
		echo "$Organism - $Strain"
		OutTable=gene_pred/swissprot/$Organism/$Strain/swissprot_v2015_tophit_parsed.tbl
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/swissprot
		$ProgDir/swissprot_parser.py --blast_tbl $SwissTable --blast_db_fasta ../../uniprot/swissprot/uniprot_sprot.fasta > $OutTable
	done
```


#Genomic analysis

## Comparison to FoL 4287

BLast searches were used to identify which genes had homologs on which
chromosomes of the Fusarium lycopersici genome.

```bash
	FoLGenomeFa=assembly/external_group/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl/Fusarium_oxysporum.FO2.31.dna.chromosome.fa
	for Proteome in $(ls gene_pred/codingquary/F.*/*/*/final_genes_combined.pep.fasta); do
	# for Proteome in $(ls assembly/external_group/F.oxysporum/fo47/broad/fusarium_oxysporum_fo47_1_proteins.fasta); do
	# for Proteome in $(ls gene_pred/external_group/F.oxysporum_fsp_lycopersici/4287/Fusox1/Fusox1_GeneCatalog_proteins_20110522_parsed.fa); do
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
