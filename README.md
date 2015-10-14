#HAPI
==========

Scripts used for the analysis of Fusarium genomes
Note - all this work was performed in the directory:
/home/groups/harrisonlab/project_files/fusarium

The following is a summary of the work presented in this Readme.

The following processes were applied to Alternaria genomes prior to analysis:
Data qc
Genome assembly
Repeatmasking
Gene prediction
Functional annotation

Analyses performed on these genomes involved BLAST searching for:


ls contigs were identified using:
Alignment of raw reads to assembled genomes
Assembly of remaining reads


#Building of directory structure
```bash
	RawDatDir=/home/groups/harrisonlab/raw_data/raw_seq/fusarium/HAPI_seq_3/
	ProjectDir=/home/groups/harrisonlab/project_files/fusarium
	mkdir -p $ProjectDir/raw_dna/paired/F.proliferatum/A8/F
	mkdir -p $ProjectDir/raw_dna/paired/F.proliferatum/A8/R
	mkdir -p $ProjectDir/raw_dna/paired/F.oxysporum_fsp_narcissi/N139/F
	mkdir -p $ProjectDir/raw_dna/paired/F.oxysporum_fsp_narcissi/N139/R
	mkdir -p $ProjectDir/raw_dna/paired/F.oxysporum_fsp_pisi/PG3/F
	mkdir -p $ProjectDir/raw_dna/paired/F.oxysporum_fsp_pisi/PG3/R
	mkdir -p $ProjectDir/raw_dna/paired/F.oxysporum_fsp_pisi/PG18/F
	mkdir -p $ProjectDir/raw_dna/paired/F.oxysporum_fsp_pisi/PG18/R
```
Sequence data was moved into the appropriate directories

```bash
	RawDatDir=/home/groups/harrisonlab/raw_data/raw_seq/fusarium/HAPI_seq_3
	ProjectDir=/home/groups/harrisonlab/project_files/fusarium
	cp $RawDatDir/Fproliferatum_S1_L001_R1_001.fastq.gz $ProjectDir/raw_dna/paired/F.proliferatum/A8/F/.
	cp $RawDatDir/Fproliferatum_S1_L001_R2_001.fastq.gz $ProjectDir/raw_dna/paired/F.proliferatum/A8/R/.
	cp $RawDatDir/FoxysporumN139_S2_L001_R1_001.fastq.gz $ProjectDir/raw_dna/paired/F.oxysporum_fsp_narcissi/N139/F/.
	cp $RawDatDir/FoxysporumN139_S2_L001_R2_001.fastq.gz $ProjectDir/raw_dna/paired/F.oxysporum_fsp_narcissi/N139/R/.
	cp $RawDatDir/FoxysporumPG3_S3_L001_R1_001.fastq.gz $ProjectDir/raw_dna/paired/F.oxysporum_fsp_pisi/PG3/F/.
	cp $RawDatDir/FoxysporumPG3_S3_L001_R2_001.fastq.gz $ProjectDir/raw_dna/paired/F.oxysporum_fsp_pisi/PG3/R/.
	cp $RawDatDir/FoxysporumPG18_S4_L001_R1_001.fastq.gz $ProjectDir/raw_dna/paired/F.oxysporum_fsp_pisi/PG18/F/.
	cp $RawDatDir/FoxysporumPG18_S4_L001_R2_001.fastq.gz $ProjectDir/raw_dna/paired/F.oxysporum_fsp_pisi/PG18/R/.
```

This process was repeated for RNAseq data:

```bash
	RawDatDir=/home/groups/harrisonlab/raw_data/raw_seq/fusarium/rna_seq
	ProjectDir=/home/groups/harrisonlab/project_files/fusarium
	mkdir -p $ProjectDir/raw_rna/paired/F.oxysporum_fsp_cepae/Fus2_PDB/F
	mkdir -p $ProjectDir/raw_rna/paired/F.oxysporum_fsp_cepae/Fus2_PDB/R
	mkdir -p $ProjectDir/raw_rna/paired/F.oxysporum_fsp_cepae/Fus2_CzapekDox/F
	mkdir -p $ProjectDir/raw_rna/paired/F.oxysporum_fsp_cepae/Fus2_CzapekDox/R
	mkdir -p $ProjectDir/raw_rna/paired/F.oxysporum_fsp_cepae/Fus2_GlucosePeptone/F
	mkdir -p $ProjectDir/raw_rna/paired/F.oxysporum_fsp_cepae/Fus2_GlucosePeptone/R
	mkdir -p $ProjectDir/raw_rna/paired/F.oxysporum_fsp_cepae/Fus2_PDA/F
	mkdir -p $ProjectDir/raw_rna/paired/F.oxysporum_fsp_cepae/Fus2_PDA/R

	cp $RawDatDir/4_S1_L001_R1_001.fastq.gz $ProjectDir/raw_rna/paired/F.oxysporum_fsp_cepae/Fus2_PDB/F/.
	cp $RawDatDir/4_S1_L001_R2_001.fastq.gz $ProjectDir/raw_rna/paired/F.oxysporum_fsp_cepae/Fus2_PDB/R/.
	cp $RawDatDir/6_S2_L001_R1_001.fastq.gz $ProjectDir/raw_rna/paired/F.oxysporum_fsp_cepae/Fus2_CzapekDox/F/.
	cp $RawDatDir/6_S2_L001_R2_001.fastq.gz $ProjectDir/raw_rna/paired/F.oxysporum_fsp_cepae/Fus2_CzapekDox/R/.
	cp $RawDatDir/7_S3_L001_R1_001.fastq.gz $ProjectDir/raw_rna/paired/F.oxysporum_fsp_cepae/Fus2_GlucosePeptone/F/.
	cp $RawDatDir/7_S3_L001_R2_001.fastq.gz $ProjectDir/raw_rna/paired/F.oxysporum_fsp_cepae/Fus2_GlucosePeptone/R/.
	cp $RawDatDir/9_S4_L001_R1_001.fastq.gz $ProjectDir/raw_rna/paired/F.oxysporum_fsp_cepae/Fus2_PDA/F/.
	cp $RawDatDir/9_S4_L001_R2_001.fastq.gz $ProjectDir/raw_rna/paired/F.oxysporum_fsp_cepae/Fus2_PDA/R/.
```



#Data qc

programs:
  fastqc
  fastq-mcf
  kmc

Data quality was visualised using fastqc:
```bash
	for RawData in $(ls raw_dna/paired/*/*/*/*.fastq.gz | grep -v 'cepae'); do
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
		echo $RawData;
		qsub $ProgDir/run_fastqc.sh $RawData
	done
```

Trimming was performed on data to trim adapters from
sequences and remove poor quality data. This was done with fastq-mcf

```bash
	for StrainPath in $(ls -d raw_dna/paired/*/* | grep -v 'cepae'); do
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
		IlluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/illumina_full_adapters.fa
		ReadsF=$(ls $StrainPath/F/*.fastq*)
		ReadsR=$(ls $StrainPath/R/*.fastq*)
		echo $ReadsF
		echo $ReadsR
		qsub $ProgDir/rna_qc_fastq-mcf.sh $ReadsF $ReadsR $IlluminaAdapters DNA
	done
```



Data quality was visualised once again following trimming:
```bash
	for RawData in qc_dna/paired/*/*/*/*.fastq*; do
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
		echo $RawData;
		qsub $ProgDir/run_fastqc.sh $RawData
	done
```

kmer counting was performed using kmc
This allowed estimation of sequencing depth and total genome size

```bash
	for TrimPath in qc_dna/paired/*/*; do
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
		TrimF=$(ls $TrimPath/F/*.fastq*)
		TrimR=$(ls $TrimPath/R/*.fastq*)
		echo $TrimF
		echo $TrimR
		qsub $ProgDir/kmc_kmer_counting.sh $TrimF $TrimR
	done
```

#Assembly

Assembly was performed using Velvet

A range of hash lengths were used and the best assembly selected for subsequent analysis

```bash
	for TrimPath in $(ls -d qc_dna/paired/*/* | grep -v 'cepae'); do
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/velvet
		Strain=$(printf $TrimPath | rev | cut -f1 -d '/' | rev)

		MinHash=41
		MaxHash=81
		HashStep=2
		TrimF=$(ls $TrimPath/F/*.fastq*)
		TrimR=$(ls $TrimPath/R/*.fastq*)
		GenomeSz=65

		echo $Strain
		if [ "$Strain" == 'PG18' ]; then
			ExpCov=24
			MinCov=8
			InsLgth=600
			echo "$Strain set"
		elif [ "$Strain" == 'PG3' ]; then
			ExpCov=27
			MinCov=9
			InsLgth=600
			echo "$Strain set"
		elif [ "$Strain" == 'N139' ]; then
			ExpCov=26
			MinCov=9
			InsLgth=600
			echo "$Strain set"
		elif [ "$Strain" == 'A8' ]; then
			ExpCov=20
			MinCov=7
			InsLgth=600
			echo "$Strain set"
		fi

		qsub $ProgDir/submit_velvet_range.sh $MinHash $MaxHash $HashStep \
		$TrimF $TrimR $GenomeSz $ExpCov $MinCov $InsLgth
	done

```


Assemblies were summarised to allow the best assembly to be determined by eye.
Although flashed reads and additional datasets were not used in this set of analyses, there were
some results from previous assemblies in the destination folders that needed to be excluded.
In the for loop that cycled through the stats file from each assembly correctly the following
commands was added: | grep -v 'flash' | grep -v 'combined' . This prevented assemblies
using flashed reads or additional datasets from being included in this summary. These
two grep expressions can be excluded if copying and pasting these commands for a different project.

```bash
	for StrainPath in $(ls -d assembly/velvet/F*/* ); do
		printf "N50\tMax_contig_size\tNumber of bases in contigs\tNumber of contigs\tNumber of contigs >=1kb\tNumber of contigs in N50\tNumber of bases in contigs >=1kb\tGC Content of contigs\n" > $StrainPath/assembly_stats.csv
		for StatsFile in $(ls $StrainPath/*/stats.txt | grep -v 'flash' | grep -v 'combined'); do
		cat $StatsFile | rev | cut -f1 -d ' ' | rev | paste -d '\t' -s >> $StrainPath/assembly_stats.csv
		done
	done
	tail -n+1 assembly/velvet/F*/*/assembly_stats.csv > assembly/velvet/Fusarium_assembly_stats.csv
```



#Repeatmasking

Repeat masking was performed and used the following programs:
	Repeatmasker
	Repeatmodeler

The best assemblies were used to perform repeatmasking

```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking
	BestAssPG18=assembly/velvet/F.*/PG18/F.*_PG18_53/sorted_contigs.fa
	BestAssPG3=assembly/velvet/F.*/PG3/F.*_PG3_79/sorted_contigs.fa
	BestAssN139=assembly/velvet/F.*/N139/F.*_N139_59/sorted_contigs.fa
	BestAssA8=assembly/velvet/F.*/A8/F.*_A8_51/sorted_contigs.fa

	qsub $ProgDir/rep_modeling.sh $BestAssPG18
	qsub $ProgDir/rep_modeling.sh $BestAssPG3
	qsub $ProgDir/rep_modeling.sh $BestAssN139
	qsub $ProgDir/rep_modeling.sh $BestAssA8

	qsub $ProgDir/transposonPSI.sh $BestAssPG18
	qsub $ProgDir/transposonPSI.sh $BestAssPG3
	qsub $ProgDir/transposonPSI.sh $BestAssN139
	qsub $ProgDir/transposonPSI.sh $BestAssA8

```


#Gene Prediction


Gene prediction followed three steps:
	Pre-gene prediction
		- Quality of genome assemblies were assessed using Cegma to see how many core eukaryotic genes can be identified.
	Gene model training
		- Gene models were trained for isolates 1166 and 650 using assembled RNAseq data
	Gene prediction
		- Gene models were used to predict genes in A. alternata genomes. This used RNAseq data as hints for gene models.

#Pre-gene prediction

Quality of genome assemblies was assessed by looking for the gene space in the assemblies.
```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/cegma
	cd /home/groups/harrisonlab/project_files/fusarium
	for Genome in $(ls repeat_masked/F.*/*/*/*_contigs_unmasked.fa | grep -v 'cepae'); do
		echo $Genome;
		qsub $ProgDir/sub_cegma.sh $Genome dna;
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

<!--
#Gene model training

Data quality was visualised using fastqc:
```bash
	for RawData in raw_rna/paired/*/*/*/*.fastq.gz; do
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
		echo $RawData;
		qsub $ProgDir/run_fastqc.sh $RawData
	done
```

Trimming was performed on data to trim adapters from
sequences and remove poor quality data. This was done with fastq-mcf

```bash
	for StrainPath in raw_rna/paired/*/*; do
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
		IlluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/illumina_full_adapters.fa
		ReadsF=$(ls $StrainPath/F/*.fastq.gz)
		ReadsR=$(ls $StrainPath/R/*.fastq.gz)
		echo $ReadsF
		echo $ReadsR
		qsub $ProgDir/rna_qc_fastq-mcf.sh $ReadsF $ReadsR $IlluminaAdapters RNA
	done
```

Data quality was visualised once again following trimming:
```bash
	for TrimData in qc_rna/paired/*/*/*/*.fastq.gz; do
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
		echo $RawData;
		qsub $ProgDir/run_fastqc.sh $TrimData
	done
```

RNAseq data was assembled into transcriptomes using Trinity
```bash
	for StrainPath in qc_rna/paired/*/*; do
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/transcriptome_assembly
		ReadsF=$(ls $StrainPath/F/*.fastq.gz)
		ReadsR=$(ls $StrainPath/R/*.fastq.gz)
		echo $ReadsF
		echo $ReadsR
		qsub $ProgDir/transcriptome_assembly_trinity.sh $ReadsF $ReadsR
	done
```
Gene training was performed using RNAseq data. The cluster can not run this script using qlogin. As such it was run on the head node (-naughty) using screen.
Training for 650 and 1166 was performed in two instances of screen and occassionally viewed to check progress over time.
(screen is detached after opening using ctrl+a then ctrl+d. - if just ctrl+d is pressed the instance of screen is deleted. - be careful)
```bash
	screen -a
	ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/augustus
	Assembly650=assembly/trinity/A.alternata_ssp._gaisen/650/650_rna_contigs/Trinity.fasta
	Genome650=repeat_masked/A.alternata_ssp._gaisen/650/A.alternata_ssp._gaisen_650_67_repmask/650_contigs_unmasked.fa
	$ProgDir/training_by_transcriptome.sh $Assembly650 $Genome650

	screen -a
	ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/augustus
	Assembly1166=assembly/trinity/A.alternata_ssp._tenuissima/1166/1166_rna_contigs/Trinity.fasta
	Genome1166=repeat_masked/A.alternata_ssp._tenuissima/1166/A.alternata_ssp._tenuissima_1166_43_repmask/1166_contigs_unmasked.fa
	$ProgDir/training_by_transcriptome.sh $Assembly1166 $Genome1166
```

Quality of Trinity assemblies were assessed using Cegma to assess gene-space within the transcriptome
```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/cegma
	for Transcriptome in $(ls assembly/trinity/A.*/*/*_rna_contigs/Trinity.fasta); do  
		echo $Transcriptome;  
		qsub $ProgDir/sub_cegma.sh $Transcriptome rna;
	done
```
Outputs were summarised using the commands:
```bash
	for File in $(ls gene_pred/cegma/A.alternata_ssp._*/*/*_rna_cegma.completeness_report); do
		Strain=$(echo $File | rev | cut -f2 -d '/' | rev);
		Species=$(echo $File | rev | cut -f3 -d '/' | rev);
		printf "$Species\t$Strain\n";
		cat $File | head -n18 | tail -n+4;printf "\n";
	done > gene_pred/cegma/cegma_results_rna_summary.txt

	less gene_pred/cegma/cegma_results_rna_summary.txt
```
 -->






#Gene prediction

Gene prediction was performed for Fusarium genomes. Two gene prediction
approaches were used:

Gene prediction using Augustus
Prediction of all putative ORFs in the genome using the ORF finder (atg.pl)
approach.

## Augustus


RNAseq reads were used as Hints for the location of CDS.

A concatenated dataset of RNAseq reads from F. oxysporum fsp. cepae isolate Fus2
were used as hints for these predictions.
A gene model trained for F.oxysporum fsp. cepae was used to describe the structure of a gene.

```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/augustus
	mkdir -p qc_rna/concatenated
	RnaFiles=$(ls qc_rna/paired/._*/*/*/*.fastq.gz | paste -s -d ' ')
	RnaF=qc_rna/paired/Fus2/aligned_appended/appended_paired.1.fastq
	RnaR=qc_rna/paired/Fus2/aligned_appended/appended_paired.2.fastq
	mkdir -p qc_rna/concatenated/F.oxysporum/Fus2
	ConcatRna=qc_rna/concatenated/F.oxysporum/Fus2/Fus2_RNA_timecourse_appended.fa.gz
	cat $RnaF $RnaR | gzip -fc > $ConcatRna
	GeneModel=F.oxysporum_fsp_cepae_Fus2
	for Genome in $(ls repeat_masked/F.*/*/*/*_contigs_unmasked.fa | grep -v 'cepae'); do
		qsub $ProgDir/augustus_pipe.sh $Genome $ConcatRna $GeneModel
	done
```

```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/augustus
	mkdir -p qc_rna/concatenated
	RnaFiles=$(ls qc_rna/paired/F.oxysporum_fsp_cepae/*/*/*_trim.fq.gz | paste -s -d ' ')
	RnaF=qc_rna/paired/Fus2/aligned_appended/appended_paired.1.fastq.gz
	RnaR=qc_rna/paired/Fus2/aligned_appended/appended_paired.2.fastq.gz
	mkdir -p qc_rna/concatenated/F.oxysporum/Fus2
	ConcatRna=qc_rna/concatenated/F.oxysporum/Fus2/Fus2_RNA_timecourse_CzapekDox_GlucosePeptone_PDA_PDB_appended.fa.gz
	cat $RnaF $RnaR $RnaFiles > $ConcatRna
	GeneModel=fusarium
	for Genome in $(ls repeat_masked/F.*/*/*/*_contigs_unmasked.fa); do
		Organism=$(echo $Genome | rev | cut -f4 -d'/' | rev)
		Strain=$(echo $Genome | rev | cut -f3 -d'/' | rev)
		OutDir=gene_pred/augustus/Model-fusarium_sp./$Organism/$Strain
		qsub $ProgDir/submit_augustus.sh $GeneModel $Genome false $OutDir
		OutDir=gene_pred/augustus/Model-fusarium_sp._Hints-Fus2/$Organism/$Strain
		qsub $ProgDir/augustus_pipe.sh $Genome $ConcatRna $GeneModel $OutDir
	done
```

```bash
for File in $(ls gene_pred/augustus/Model-fusarium_sp._Hints-Fus2/F.*/*/*_augustus_preds.aa); do
	echo $File; cat $File | grep '>' | wc -l;
done
# gene_pred/augustus/Model-fusarium_sp._Hints-Fus2/F.oxysporum_fsp_cepae/125/125_augustus_preds.aa
# 17261
# gene_pred/augustus/Model-fusarium_sp._Hints-Fus2/F.oxysporum_fsp_cepae/55/55_augustus_preds.aa
# 17013
# gene_pred/augustus/Model-fusarium_sp._Hints-Fus2/F.oxysporum_fsp_cepae/A23/A23_augustus_preds.aa
# 17101
# gene_pred/augustus/Model-fusarium_sp._Hints-Fus2/F.oxysporum_fsp_cepae/A28/A28_augustus_preds.aa
# 17292
# gene_pred/augustus/Model-fusarium_sp._Hints-Fus2/F.oxysporum_fsp_cepae/D2/D2_augustus_preds.aa
# 16727
# gene_pred/augustus/Model-fusarium_sp._Hints-Fus2/F.oxysporum_fsp_cepae/Fus2/Fus2_augustus_preds.aa
# 16988
# gene_pred/augustus/Model-fusarium_sp._Hints-Fus2/F.oxysporum_fsp_cepae/HB17/HB17_augustus_preds.aa
# 17077
# gene_pred/augustus/Model-fusarium_sp._Hints-Fus2/F.oxysporum_fsp_cepae/PG/PG_augustus_preds.aa
# 16961
# gene_pred/augustus/Model-fusarium_sp._Hints-Fus2/F.oxysporum_fsp_narcissi/N139/N139_augustus_preds.aa
# 21539
# gene_pred/augustus/Model-fusarium_sp._Hints-Fus2/F.oxysporum_fsp_pisi/PG18/PG18_augustus_preds.aa
# 21618
# gene_pred/augustus/Model-fusarium_sp._Hints-Fus2/F.oxysporum_fsp_pisi/PG3/PG3_augustus_preds.aa
# 18077
# gene_pred/augustus/Model-fusarium_sp._Hints-Fus2/F.proliferatum/A8/A8_augustus_preds.aa
# 20681
for File in $(ls gene_pred/augustus/Model-fusarium_sp./*/*/*_EMR_singlestrand_aug_out.aa); do
	echo $File;
	cat $File | grep '>' | wc -l;
done
# gene_pred/augustus/Model-fusarium_sp./F.oxysporum_fsp_cepae/125/125_EMR_singlestrand_aug_out.aa
# 18650
# gene_pred/augustus/Model-fusarium_sp./F.oxysporum_fsp_cepae/55/55_EMR_singlestrand_aug_out.aa
# 18293
# gene_pred/augustus/Model-fusarium_sp./F.oxysporum_fsp_cepae/A23/A23_EMR_singlestrand_aug_out.aa
# 18363
# gene_pred/augustus/Model-fusarium_sp./F.oxysporum_fsp_cepae/A28/A28_EMR_singlestrand_aug_out.aa
# 18730
# gene_pred/augustus/Model-fusarium_sp./F.oxysporum_fsp_cepae/D2/D2_EMR_singlestrand_aug_out.aa
# 18432
# gene_pred/augustus/Model-fusarium_sp./F.oxysporum_fsp_cepae/Fus2/Fus2_EMR_singlestrand_aug_out.aa
# 18108
# gene_pred/augustus/Model-fusarium_sp./F.oxysporum_fsp_cepae/HB17/HB17_EMR_singlestrand_aug_out.aa
# 18282
# gene_pred/augustus/Model-fusarium_sp./F.oxysporum_fsp_cepae/PG/PG_EMR_singlestrand_aug_out.aa
# 18431
# gene_pred/augustus/Model-fusarium_sp./F.oxysporum_fsp_narcissi/N139/N139_EMR_singlestrand_aug_out.aa
# 23374
# gene_pred/augustus/Model-fusarium_sp./F.oxysporum_fsp_pisi/PG18/PG18_EMR_singlestrand_aug_out.aa
# 23346
# gene_pred/augustus/Model-fusarium_sp./F.oxysporum_fsp_pisi/PG3/PG3_EMR_singlestrand_aug_out.aa
# 19906
# gene_pred/augustus/Model-fusarium_sp./F.proliferatum/A8/A8_EMR_singlestrand_aug_out.aa
# 22116
```

## ORF finder

The genome was searched in six reading frames for any start codon and following
translated identification of a start codon translating sequence until a stop
codon was found. This is based upon the atg.pl script used in paper describing
the P. infestans genome. Additional functionality was added to this script by
also printing ORFs in .gff format.


```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
	for Genome in $(ls repeat_masked/F.*/*/*/*_contigs_unmasked.fa); do
		qsub $ProgDir/run_ORF_finder.sh $Genome
	done
```

The Gff files produced by ORF finder have some formatting errors. The following
script was run to correct these:

```bash
	ProgDir=~/git_repos/emr_repos/tools/seq_tools/feature_annotation
	for ORF_Gff in $(ls gene_pred/ORF_finder/*/*/*_ORF.gff | grep -v '_F_atg_' | grep -v '_R_atg_'); do
		ORF_Gff_mod=$(echo $ORF_Gff | sed 's/_ORF.gff/_ORF_mod.gff/g')
		echo ""
		echo "Correcting the following file:"
		echo $ORF_Gff
		echo "Redirecting to:"
		echo $ORF_Gff_mod
		$ProgDir/gff_corrector.pl $ORF_Gff > $ORF_Gff_mod
	done
```


#Functional annotation

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
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan/
	for Genes in $(ls gene_pred/augustus/F.*/*/*_augustus_preds.aa | grep -v 'cepae'); do
		echo $Genes
		$ProgDir/sub_interproscan.sh $Genes
	done 2>&1 |  tee -a interproscan_submisison.log
```

Following interproscan annotation split files were combined using the following
commands:

```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
	for StrainPath in $(ls -d gene_pred/interproscan/F.*/*); do
		Strain=$(basename $StrainPath)
		Organism=$(echo $StrainPath | rev | cut -d "/" -f2 | rev)
		echo $Strain
		PredGenes=gene_pred/augustus/"$Organism"/"$Strain"/"$Strain"_augustus_preds.aa
		InterProRaw=gene_pred/interproscan/"$Organism"/"$Strain"/raw
		$ProgDir/append_interpro.sh $PredGenes $InterProRaw
	done
```

#Genomic analysis

## RxLR genes

Putative RxLR genes were identified within Augustus gene models using a number
of approaches:

 * A) From Augustus gene models - Signal peptide & RxLR motif  
 * B) From Augustus gene models - Hmm evidence of WY domains  
 * C) From Augustus gene models - Hmm evidence of RxLR effectors  
 * D) From ORF fragments - Signal peptide & RxLR motif  
 * E) From ORF fragments - Hmm evidence of WY domains  
 * F) From ORF fragments - Hmm evidence of RxLR effectors  


### A) From Augustus gene models - Signal peptide & RxLR motif

Required programs:
 * SigP
 * biopython


Proteins that were predicted to contain signal peptides were identified using
the following commands:

```bash
	SplitfileDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
	CurPath=$PWD
	for Proteome in $(ls gene_pred/augustus/F.*/*/*_augustus_preds.aa | grep -v '_old'); do
		Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		SplitDir=gene_pred/augustus_split/$Organism/$Strain
		mkdir -p $SplitDir
		BaseName="$Organism""_$Strain"_augustus_preds
		$SplitfileDir/splitfile_500.py --inp_fasta $Proteome --out_dir $SplitDir --out_base $BaseName
		for File in $(ls $SplitDir/*_augustus_preds_*); do
			Jobs=$(qstat | grep 'pred_sigP' | wc -l)
			while [ $Jobs -ge 32 ]; do
				sleep 10
				printf "."
				Jobs=$(qstat | grep 'pred_sigP' | wc -l)
			done
			printf "\n"
			echo $File
			qsub $ProgDir/pred_sigP.sh $File
		done
	done
```

The batch files of predicted secreted proteins needed to be combined into a
single file for each strain. This was done with the following commands:
```bash
	for SplitDir in $(ls -d gene_pred/sigP/F.*/*/split); do
		Strain=$(echo $SplitDir | cut -d '/' -f4)
		Organism=$(echo $SplitDir | cut -d '/' -f3)
		InStringAA=''
		InStringNeg=''
		InStringTab=''
		InStringTxt=''
		for GRP in $(ls -l $SplitDir/*_augustus_preds_*_sp.aa | rev | cut -d '_' -f2 | rev | sort -n); do  
		InStringAA="$InStringAA gene_pred/sigP/$Organism/$Strain/split/"$Organism"_"$Strain"_augustus_preds_$GRP""_sp.aa";  
		InStringNeg="$InStringNeg gene_pred/sigP/$Organism/$Strain/split/"$Organism"_"$Strain"_augustus_preds_$GRP""_sp_neg.aa";  
		InStringTab="$InStringTab gene_pred/sigP/$Organism/$Strain/split/"$Organism"_"$Strain"_augustus_preds_$GRP""_sp.tab";
		InStringTxt="$InStringTxt gene_pred/sigP/$Organism/$Strain/split/"$Organism"_"$Strain"_augustus_preds_$GRP""_sp.txt";  
		done
		cat $InStringAA > gene_pred/sigP/$Organism/$Strain/"$Strain"_aug_sp.aa
		cat $InStringNeg > gene_pred/sigP/$Organism/$Strain/"$Strain"_aug_neg_sp.aa
		tail -n +2 -q $InStringTab > gene_pred/sigP/$Organism/$Strain/"$Strain"_aug_sp.tab
		cat $InStringTxt > gene_pred/sigP/$Organism/$Strain/"$Strain"_aug_sp.txt
	done
```
The regular expression R.LR.{,40}[ED][ED][KR] has previously been used to identfy RxLR effectors. The addition of an EER motif is significant as it has been shown as required for host uptake of the protein.

The RxLR_EER_regex_finder.py script was used to search for this regular expression and annotate the EER domain where present.

```bash
	for Secretome in $(ls gene_pred/sigP/F.*/*/*_aug_sp.aa); do
		ProgDir=~/git_repos/emr_repos/tools/pathogen/RxLR_effectors;
		Strain=$(echo $Secretome | cut -d '/' -f4);
		Organism=$(echo $Secretome | cut -d '/' -f3) ;
		OutDir=analysis/RxLR_effectors/RxLR_EER_regex_finder/"$Organism"/"$Strain";
		mkdir -p $OutDir;
		printf "\nstrain: $Strain\tspecies: $Organism\n";
		printf "the number of SigP gene is:\t";
		cat $Secretome | grep '>' | wc -l;
		printf "the number of SigP-RxLR genes are:\t";
		$ProgDir/RxLR_EER_regex_finder.py $Secretome > $OutDir/"$Strain"_Aug_RxLR_EER_regex.fa;
		cat $OutDir/"$Strain"_Aug_RxLR_EER_regex.fa | grep '>' | cut -f1 | sed 's/>//g' | sed 's/ //g' > $OutDir/"$Strain"_Aug_RxLR_regex.txt
		cat $OutDir/"$Strain"_Aug_RxLR_regex.txt | wc -l
		printf "the number of SigP-RxLR-EER genes are:\t";
		cat $OutDir/"$Strain"_Aug_RxLR_EER_regex.fa | grep '>' | grep 'EER_motif_start' |  cut -f1 | sed 's/>//g' | sed 's/ //g' > $OutDir/"$Strain"_Aug_RxLR_EER_regex.txt
		cat $OutDir/"$Strain"_Aug_RxLR_EER_regex.txt | wc -l
		printf "\n"
		# ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
		# Col2=RxLR_EER_regex_finder.py
		# GeneNames=$OutDir/"$Strain"_Aug_RxLR_regex.txt
		# GeneModels=gene_pred/augustus/"$Organism"/"$Strain"/"$Strain"_augustus_preds.aa
		# $ProgDir/gene_list_to_gff.pl $GeneNames $GeneModels $Col2 Name > $OutDir/"$Strain"_Aug_RxLR_regex.gff3
	done
```

The results were as follows:

strain: 125	species: F.oxysporum_fsp_cepae  
the number of SigP gene is:	1340  
the number of SigP-RxLR genes are:	13  
the number of SigP-RxLR-EER genes are:	1  

strain: 55	species: F.oxysporum_fsp_cepae  
the number of SigP gene is:	1333  
the number of SigP-RxLR genes are:	12  
the number of SigP-RxLR-EER genes are:	1  

strain: A23	species: F.oxysporum_fsp_cepae  
the number of SigP gene is:	1335  
the number of SigP-RxLR genes are:	13  
the number of SigP-RxLR-EER genes are:	2  

strain: A28	species: F.oxysporum_fsp_cepae  
the number of SigP gene is:	1354  
the number of SigP-RxLR genes are:	16  
the number of SigP-RxLR-EER genes are:	1  

strain: D2	species: F.oxysporum_fsp_cepae  
the number of SigP gene is:	1156  
the number of SigP-RxLR genes are:	8  
the number of SigP-RxLR-EER genes are:	0  

strain: Fus2	species: F.oxysporum_fsp_cepae  
the number of SigP gene is:	1346  
the number of SigP-RxLR genes are:	13  
the number of SigP-RxLR-EER genes are:	1  

strain: HB17	species: F.oxysporum_fsp_cepae  
the number of SigP gene is:	1350  
the number of SigP-RxLR genes are:	13  
the number of SigP-RxLR-EER genes are:	1  

strain: PG	species: F.oxysporum_fsp_cepae  
the number of SigP gene is:	1334  
the number of SigP-RxLR genes are:	11  
the number of SigP-RxLR-EER genes are:	0  

strain: N139	species: F.oxysporum_fsp_narcissi  
the number of SigP gene is:	1765  
the number of SigP-RxLR genes are:	22  
the number of SigP-RxLR-EER genes are:	1  

strain: PG18	species: F.oxysporum_fsp_pisi  
the number of SigP gene is:	1807  
the number of SigP-RxLR genes are:	28  
the number of SigP-RxLR-EER genes are:	1  

strain: PG3	species: F.oxysporum_fsp_pisi  
the number of SigP gene is:	1386  
the number of SigP-RxLR genes are:	11  
the number of SigP-RxLR-EER genes are:	0  

strain: A8	species: F.proliferatum  
the number of SigP gene is:	1975  
the number of SigP-RxLR genes are:	40  
the number of SigP-RxLR-EER genes are:	3  



### B) From Augustus gene models - Hmm evidence of WY domains
Hmm models for the WY domain contained in many RxLRs were used to search gene models predicted with Augustus. These were run with the following commands:

```bash
	ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
	HmmModel=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer/WY_motif.hmm
	for Proteome in $(ls gene_pred/augustus/F.*/*/*_augustus_preds.aa | grep -v '_old'); do
		Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		OutDir=analysis/RxLR_effectors/hmmer_WY/$Organism/$Strain
		mkdir -p $OutDir
		HmmResults="$Strain"_Aug_WY_hmmer.txt
		hmmsearch -T 0 $HmmModel $Proteome > $OutDir/$HmmResults
		echo "$Organism $Strain"
		cat $OutDir/$HmmResults | grep 'Initial search space'
		cat $OutDir/$HmmResults | grep 'number of targets reported over threshold'
		HmmFasta="$Strain"_Aug_WY_hmmer.fa
		$ProgDir/hmmer2fasta.pl $OutDir/$HmmResults $Proteome > $OutDir/$HmmFasta
	done
```

Results were as follows:  

F.oxysporum_fsp_cepae 125
Initial search space (Z):              11203  [actual number of targets]  
Domain search space  (domZ):               0  [number of targets reported over threshold]  
F.oxysporum_fsp_cepae 55  
Initial search space (Z):              11047  [actual number of targets]  
Domain search space  (domZ):               0  [number of targets reported over threshold]  
F.oxysporum_fsp_cepae A23  
Initial search space (Z):              11108  [actual number of targets]  
Domain search space  (domZ):               0  [number of targets reported over threshold]  
F.oxysporum_fsp_cepae A28  
Initial search space (Z):              11131  [actual number of targets]  
Domain search space  (domZ):               0  [number of targets reported over threshold]  
F.oxysporum_fsp_cepae D2  
Initial search space (Z):              10632  [actual number of targets]  
Domain search space  (domZ):               0  [number of targets reported over threshold]  
F.oxysporum_fsp_cepae Fus2  
Initial search space (Z):              11104  [actual number of targets]  
Domain search space  (domZ):               0  [number of targets reported over threshold]  
F.oxysporum_fsp_cepae HB17  
Initial search space (Z):              11126  [actual number of targets]  
Domain search space  (domZ):               0  [number of targets reported over threshold]  
F.oxysporum_fsp_cepae PG  
Initial search space (Z):              11011  [actual number of targets]  
Domain search space  (domZ):               0  [number of targets reported over threshold]  
F.oxysporum_fsp_narcissi N139  
Initial search space (Z):              14757  [actual number of targets]  
Domain search space  (domZ):               0  [number of targets reported over threshold]  
F.oxysporum_fsp_pisi PG18  
Initial search space (Z):              14859  [actual number of targets]  
Domain search space  (domZ):               0  [number of targets reported over threshold]  
F.oxysporum_fsp_pisi PG3  
Initial search space (Z):              11517  [actual number of targets]  
Domain search space  (domZ):               0  [number of targets reported over threshold]  
F.proliferatum A8  
Initial search space (Z):              15034  [actual number of targets]  
Domain search space  (domZ):               1  [number of targets reported over threshold]  



### C) From Augustus gene models - Hmm evidence of RxLR effectors
```bash
	ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
	HmmModel=/home/armita/git_repos/emr_repos/SI_Whisson_et_al_2007/cropped.hmm
	for Proteome in $(ls gene_pred/augustus/F.*/*/*_augustus_preds.aa | grep -v '_old'); do
		Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		OutDir=analysis/RxLR_effectors/hmmer_RxLR/$Organism/$Strain
		mkdir -p $OutDir
		HmmResults="$Strain"_Aug_RxLR_hmmer.txt
		hmmsearch -T 0 $HmmModel $Proteome > $OutDir/$HmmResults
		echo "$Organism $Strain"
		cat $OutDir/$HmmResults | grep 'Initial search space'
		cat $OutDir/$HmmResults | grep 'number of targets reported over threshold'
		HmmFasta="$Strain"__Aug_RxLR_hmmer.fa
		$ProgDir/hmmer2fasta.pl $OutDir/$HmmResults $Proteome > $OutDir/$HmmFasta
	done
```

F.oxysporum_fsp_cepae 125  
Initial search space (Z):              11203  [actual number of targets]  
Domain search space  (domZ):               3  [number of targets reported over threshold]  
F.oxysporum_fsp_cepae 55  
Initial search space (Z):              11047  [actual number of targets]  
Domain search space  (domZ):               3  [number of targets reported over threshold]  
F.oxysporum_fsp_cepae A23  
Initial search space (Z):              11108  [actual number of targets]  
Domain search space  (domZ):               2  [number of targets reported over threshold]  
F.oxysporum_fsp_cepae A28  
Initial search space (Z):              11131  [actual number of targets]  
Domain search space  (domZ):               1  [number of targets reported over threshold]  
F.oxysporum_fsp_cepae D2  
Initial search space (Z):              10632  [actual number of targets]  
Domain search space  (domZ):               1  [number of targets reported over threshold]  
F.oxysporum_fsp_cepae Fus2  
Initial search space (Z):              11104  [actual number of targets]  
Domain search space  (domZ):               2  [number of targets reported over threshold]  
F.oxysporum_fsp_cepae HB17  
Initial search space (Z):              11126  [actual number of targets]  
Domain search space  (domZ):               2  [number of targets reported over threshold]  
F.oxysporum_fsp_cepae PG  
Initial search space (Z):              11011  [actual number of targets]  
Domain search space  (domZ):               2  [number of targets reported over threshold]  
F.oxysporum_fsp_narcissi N139  
Initial search space (Z):              14757  [actual number of targets]  
Domain search space  (domZ):               1  [number of targets reported over threshold]  
F.oxysporum_fsp_pisi PG18  
Initial search space (Z):              14859  [actual number of targets]  
Domain search space  (domZ):               2  [number of targets reported over threshold]  
F.oxysporum_fsp_pisi PG3  
Initial search space (Z):              11517  [actual number of targets]  
Domain search space  (domZ):               2  [number of targets reported over threshold]  
F.proliferatum A8  
Initial search space (Z):              15034  [actual number of targets]  
Domain search space  (domZ):               2  [number of targets reported over threshold]  

### D) From ORF fragments - Signal peptide & RxLR motif

Required programs:
 * SigP
 * biopython


Proteins that were predicted to contain signal peptides were identified using
the following commands:

```bash
	SplitfileDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
	CurPath=$PWD
	for Proteome in $(ls gene_pred/ORF_finder/F.*/*/*.aa_cat.fa); do
		Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		SplitDir=gene_pred/ORF_finder_split/$Organism/$Strain
		mkdir -p $SplitDir
		BaseName="$Organism""_$Strain"_ORF_finder_preds
		$SplitfileDir/splitfile_500.py --inp_fasta $Proteome --out_dir $SplitDir --out_base $BaseName
		for File in $(ls $SplitDir/$BaseName_*.fa); do
			Jobs=$(qstat | grep 'pred_sigP' | wc -l)
			while [ $Jobs -ge 32 ]; do
				sleep 10
				printf "."
				Jobs=$(qstat | grep 'pred_sigP' | wc -l)
			done
			printf "\n"
			echo $File
			qsub $ProgDir/pred_sigP.sh $File
		done
	done
```

The batch files of predicted secreted proteins needed to be combined into a
single file for each strain. This was done with the following commands:
```bash
	for SplitDir in $(ls -d gene_pred/ORF_finder_sigP/F.*/*/split); do
	Strain=$(echo $SplitDir | cut -d '/' -f4)
	Organism=$(echo $SplitDir | cut -d '/' -f3)
	InStringAA=''
	InStringNeg=''
	InStringTab=''
	InStringTxt=''
	for GRP in $(ls -l $SplitDir/*_ORF_finder_preds_*_sp.aa | rev | cut -d '_' -f2 | rev | sort -n); do  
	InStringAA="$InStringAA gene_pred/ORF_finder_sigP/$Organism/$Strain/split/"$Organism"_"$Strain"_ORF_finder_preds_$GRP""_sp.aa";  
	InStringNeg="$InStringNeg gene_pred/ORF_finder_sigP/$Organism/$Strain/split/"$Organism"_"$Strain"_ORF_finder_preds_$GRP""_sp_neg.aa";  
	InStringTab="$InStringTab gene_pred/ORF_finder_sigP/$Organism/$Strain/split/"$Organism"_"$Strain"_ORF_finder_preds_$GRP""_sp.tab";
	InStringTxt="$InStringTxt gene_pred/ORF_finder_sigP/$Organism/$Strain/split/"$Organism"_"$Strain"_ORF_finder_preds_$GRP""_sp.txt";  
	done
	cat $InStringAA > gene_pred/ORF_finder_sigP/$Organism/$Strain/"$Strain"_ORF_finder_sp.aa
	cat $InStringNeg > gene_pred/ORF_finder_sigP/$Organism/$Strain/"$Strain"_ORF_finder_neg_sp.aa
	tail -n +2 -q $InStringTab > gene_pred/ORF_finder_sigP/$Organism/$Strain/"$Strain"_ORF_finder_sp.tab
	cat $InStringTxt > gene_pred/ORF_finder_sigP/$Organism/$Strain/"$Strain"_ORF_finder_sp.txt
	done
```
The regular expression R.LR.{,40}[ED][ED][KR] has previously been used to identfy RxLR effectors. The addition of an EER motif is significant as it has been shown as required for host uptake of the protein.

The RxLR_EER_regex_finder.py script was used to search for this regular expression and annotate the EER domain where present.

```bash
	for Secretome in $(ls gene_pred/ORF_finder_sigP/F.*/*/*_ORF_finder_sp.aa); do
		ProgDir=~/git_repos/emr_repos/tools/pathogen/RxLR_effectors;
		Strain=$(echo $Secretome | cut -d '/' -f4);
		Organism=$(echo $Secretome | cut -d '/' -f3) ;
		OutDir=analysis/RxLR_effectors/RxLR_EER_regex_finder/"$Organism"/"$Strain";
		mkdir -p $OutDir;
		printf "\nstrain: $Strain\tspecies: $Organism\n";
		printf "the number of SigP gene is:\t";
		cat $Secretome | grep '>' | wc -l;
		printf "the number of SigP-RxLR genes are:\t";
		$ProgDir/RxLR_EER_regex_finder.py $Secretome > $OutDir/"$Strain"_ORF_RxLR_EER_regex.fa;
		cat $OutDir/"$Strain"_ORF_RxLR_EER_regex.fa | grep '>' | cut -f1 | sed 's/>//g' | sed 's/ //g' > $OutDir/"$Strain"_ORF_RxLR_regex.txt
		cat $OutDir/"$Strain"_ORF_RxLR_regex.txt | wc -l
		printf "the number of SigP-RxLR-EER genes are:\t";
		cat $OutDir/"$Strain"_ORF_RxLR_EER_regex.fa | grep '>' | grep 'EER_motif_start' |  cut -f1 | sed 's/>//g' | sed 's/ //g' > $OutDir/"$Strain"_ORF_RxLR_EER_regex.txt
		cat $OutDir/"$Strain"_ORF_RxLR_EER_regex.txt | wc -l
		printf "\n"
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
		Col2=RxLR_EER_regex_finder.py
		#GeneNames=$OutDir/"$Strain"_ORF_RxLR_regex.txt
		#GeneModels=gene_pred/ORF_finder/"$Organism"/"$Strain"/"$Strain"*.aa_cat.fa
		#$ProgDir/gene_list_to_gff.pl $GeneNames $GeneModels $Col2 Name > $OutDir/"$Strain"_ORF_RxLR_regex.gff3
	done
```

The results were as follows:
```
	strain: 125	species: F.oxysporum_fsp_cepae
	the number of SigP gene is:	15312
	the number of SigP-RxLR genes are:	685
	the number of SigP-RxLR-EER genes are:	17


	strain: 55	species: F.oxysporum_fsp_cepae
	the number of SigP gene is:	16347
	the number of SigP-RxLR genes are:	724
	the number of SigP-RxLR-EER genes are:	17


	strain: A23	species: F.oxysporum_fsp_cepae
	the number of SigP gene is:	16402
	the number of SigP-RxLR genes are:	731
	the number of SigP-RxLR-EER genes are:	17


	strain: A28	species: F.oxysporum_fsp_cepae
	the number of SigP gene is:	16617
	the number of SigP-RxLR genes are:	716
	the number of SigP-RxLR-EER genes are:	16


	strain: D2	species: F.oxysporum_fsp_cepae
	the number of SigP gene is:	16169
	the number of SigP-RxLR genes are:	686
	the number of SigP-RxLR-EER genes are:	15


	strain: Fus2	species: F.oxysporum_fsp_cepae
	the number of SigP gene is:	16262
	the number of SigP-RxLR genes are:	722
	the number of SigP-RxLR-EER genes are:	18


	strain: HB17	species: F.oxysporum_fsp_cepae
	the number of SigP gene is:	16326
	the number of SigP-RxLR genes are:	720
	the number of SigP-RxLR-EER genes are:	17


	strain: PG	species: F.oxysporum_fsp_cepae
	the number of SigP gene is:	16390
	the number of SigP-RxLR genes are:	720
	the number of SigP-RxLR-EER genes are:	15


	strain: N139	species: F.oxysporum_fsp_narcissi
	the number of SigP gene is:	22180
	the number of SigP-RxLR genes are:	986
	the number of SigP-RxLR-EER genes are:	22


	strain: PG18	species: F.oxysporum_fsp_pisi
	the number of SigP gene is:	22146
	the number of SigP-RxLR genes are:	1014
	the number of SigP-RxLR-EER genes are:	14


	strain: PG3	species: F.oxysporum_fsp_pisi
	the number of SigP gene is:	17424
	the number of SigP-RxLR genes are:	740
	the number of SigP-RxLR-EER genes are:	11


	strain: A8	species: F.proliferatum
	the number of SigP gene is:	24561
	the number of SigP-RxLR genes are:	1203
	the number of SigP-RxLR-EER genes are:	27
```

### E) From ORF fragments - Hmm evidence of WY domains
```bash
	ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
	HmmModel=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer/WY_motif.hmm
	for Proteome in $(ls gene_pred/ORF_finder/F.*/*/*.aa_cat.fa); do
		Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		OutDir=analysis/RxLR_effectors/hmmer_WY/$Organism/$Strain
		mkdir -p $OutDir
		HmmResults="$Strain"_ORF_WY_hmmer.txt
		hmmsearch -T 0 $HmmModel $Proteome > $OutDir/$HmmResults
		echo "$Organism $Strain"
		cat $OutDir/$HmmResults | grep 'Initial search space'
		cat $OutDir/$HmmResults | grep 'number of targets reported over threshold'
		HmmFasta="$Strain"_ORF_WY_hmmer.fa
		$ProgDir/hmmer2fasta.pl $OutDir/$HmmResults $Proteome > $OutDir/$HmmFasta
	done
```
Results were as follows:
```
	F.oxysporum_fsp_cepae 125
	Initial search space (Z):             357389  [actual number of targets]
	Domain search space  (domZ):              11  [number of targets reported over threshold]
	F.oxysporum_fsp_cepae 55
	Initial search space (Z):             356765  [actual number of targets]
	Domain search space  (domZ):              11  [number of targets reported over threshold]
	F.oxysporum_fsp_cepae A23
	Initial search space (Z):             357686  [actual number of targets]
	Domain search space  (domZ):              11  [number of targets reported over threshold]
	F.oxysporum_fsp_cepae A28
	Initial search space (Z):             365675  [actual number of targets]
	Domain search space  (domZ):              10  [number of targets reported over threshold]
	F.oxysporum_fsp_cepae D2
	Initial search space (Z):             352139  [actual number of targets]
	Domain search space  (domZ):               9  [number of targets reported over threshold]
	F.oxysporum_fsp_cepae Fus2
	Initial search space (Z):             353490  [actual number of targets]
	Domain search space  (domZ):              11  [number of targets reported over threshold]
	F.oxysporum_fsp_cepae HB17
	Initial search space (Z):             355148  [actual number of targets]
	Domain search space  (domZ):              11  [number of targets reported over threshold]
	F.oxysporum_fsp_cepae PG
	Initial search space (Z):             356817  [actual number of targets]
	Domain search space  (domZ):              15  [number of targets reported over threshold]
	F.oxysporum_fsp_narcissi N139
	Initial search space (Z):             449092  [actual number of targets]
	Domain search space  (domZ):               9  [number of targets reported over threshold]
	F.oxysporum_fsp_pisi PG18
	Initial search space (Z):             445079  [actual number of targets]
	Domain search space  (domZ):              13  [number of targets reported over threshold]
	F.oxysporum_fsp_pisi PG3
	Initial search space (Z):             384712  [actual number of targets]
	Domain search space  (domZ):              10  [number of targets reported over threshold]
	F.proliferatum A8
	Initial search space (Z):             456646  [actual number of targets]
	Domain search space  (domZ):               6  [number of targets reported over threshold]
```

### F) From ORF fragments - Hmm evidence of RxLR effectors
```bash
	ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
	HmmModel=/home/armita/git_repos/emr_repos/SI_Whisson_et_al_2007/cropped.hmm
	for Proteome in $(ls gene_pred/ORF_finder/F.*/*/*.aa_cat.fa); do
		Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		OutDir=analysis/RxLR_effectors/hmmer_RxLR/$Organism/$Strain
		mkdir -p $OutDir
		HmmResults="$Strain"_ORF_RxLR_hmmer.txt
		hmmsearch -T 0 $HmmModel $Proteome > $OutDir/$HmmResults
		echo "$Organism $Strain"
		cat $OutDir/$HmmResults | grep 'Initial search space'
		cat $OutDir/$HmmResults | grep 'number of targets reported over threshold'
		HmmFasta="$Strain"_ORF_RxLR_hmmer.fa
		$ProgDir/hmmer2fasta.pl $OutDir/$HmmResults $Proteome > $OutDir/$HmmFasta
	done
```

The results were as follows:

```
	F.oxysporum_fsp_cepae 125
	Initial search space (Z):             357389  [actual number of targets]
	Domain search space  (domZ):              46  [number of targets reported over threshold]
	F.oxysporum_fsp_cepae 55
	Initial search space (Z):             356765  [actual number of targets]
	Domain search space  (domZ):              46  [number of targets reported over threshold]
	F.oxysporum_fsp_cepae A23
	Initial search space (Z):             357686  [actual number of targets]
	Domain search space  (domZ):              46  [number of targets reported over threshold]
	F.oxysporum_fsp_cepae A28
	Initial search space (Z):             365675  [actual number of targets]
	Domain search space  (domZ):              45  [number of targets reported over threshold]
	F.oxysporum_fsp_cepae D2
	Initial search space (Z):             352139  [actual number of targets]
	Domain search space  (domZ):              41  [number of targets reported over threshold]
	F.oxysporum_fsp_cepae Fus2
	Initial search space (Z):             353490  [actual number of targets]
	Domain search space  (domZ):              46  [number of targets reported over threshold]
	F.oxysporum_fsp_cepae HB17
	Initial search space (Z):             355148  [actual number of targets]
	Domain search space  (domZ):              46  [number of targets reported over threshold]
	F.oxysporum_fsp_cepae PG
	Initial search space (Z):             356817  [actual number of targets]
	Domain search space  (domZ):              51  [number of targets reported over threshold]
	F.oxysporum_fsp_narcissi N139
	Initial search space (Z):             449092  [actual number of targets]
	Domain search space  (domZ):              58  [number of targets reported over threshold]
	F.oxysporum_fsp_pisi PG18
	Initial search space (Z):             445079  [actual number of targets]
	Domain search space  (domZ):              48  [number of targets reported over threshold]
	F.oxysporum_fsp_pisi PG3
	Initial search space (Z):             384712  [actual number of targets]
	Domain search space  (domZ):              35  [number of targets reported over threshold]
	F.proliferatum A8
	Initial search space (Z):             456646  [actual number of targets]
	Domain search space  (domZ):              45  [number of targets reported over threshold]
```

### G) Combined evidence file

Evidence of RxLR effectors from different sources was combined into a single gff
feature file.

This was done using the following commands:
<!--
```bash
		ProgDir=~/git_repos/emr_repos/scripts/phytophthora/pathogen/merge_gff
    # Infiles
    GffAug=gene_pred/augustus_unmasked/P.cactorum/414/414_augustus_preds.gtf
    ORF_RxLRs=analysis/rxlr_atg_unmasked/P.cactorum/414/414_ORF_sp_rxlr.gff3
    ORF_WYs=analysis/hmmer/WY/P.cactorum/414/414_ORF_WY_hmmer.gff3
    # Outfiles
    AugDB=414_aug.db
    WyDB=414_WY.db
    WyID=WY_id.txt
    WyDB_mod=414_WY_note.db
    RxlrDB=414_rxlr.db
    RxlrID=rxlr_id.txt
    RxlrDB_mod=414_rxlr_note.db
    Rxlr_Wy_DB=414_rxlr_WY.db

    OrfMerged=414_rxlr_WY_merged.db
    MergedDB=414_Aug_ORF_merged.db
    FinalDB=414_Aug_ORF.db
    FinalGff=414_Aug_ORF.gff
```

Make a db of aug genes and effector ORFs

```bash
$ProgDir/make_gff_database.py --inp $GffAug --db $AugDB
	$ProgDir/make_gff_database.py --inp $ORF_WYs --db $WyDB
	$ProgDir/make_gff_database.py --inp $ORF_RxLRs --db $RxlrDB
``` -->


## Mimps

The presence of Mimp promotors in Fusarium genomes were identified. This was
done in three steps:
 * Position of Mimps were identified in the genome
 * Genes within 1000bp downstream of the mimp were identified from Augustus
predictions
 * ORFs within 1000bp downstream of the mimp were identified from ORF
predictions

### A) Position of Mimps were identified
Position of Mimps
gff predictions for the position of mimps in the genome were identified
Fus2 was performed separately to ensure mimps are predicted from the correct assembly

```bash
	ProgDir="/home/armita/git_repos/emr_repos/tools/pathogen/mimp_finder"
	for Genome in $(ls repeat_masked/F.*/*/*/*_contigs_unmasked.fa); do
		Organism=$(echo "$Genome" | rev | cut -d '/' -f4 | rev)
		Strain=$(echo "$Genome" | rev | cut -d '/' -f3 | rev)
		OutDir=analysis/mimps/$Organism/"$Strain"
		mkdir -p "$OutDir"
		"$ProgDir"/mimp_finder.pl "$Genome" "$OutDir"/"$Strain"_mimps.fa "$OutDir"/"$Strain"_mimps.gff3 > "$OutDir"/"$Strain"_mimps.log
	done
```

```bash
ProgDir="/home/armita/git_repos/emr_repos/tools/pathogen/mimp_finder"
# for Genome in $(ls assembly/external_group/F.oxysporum_fsp_lycopersici/4287/Fusox1/Fusox1_AssemblyScaffolds.fasta); do
for Genome in $(ls assembly/external_group/F.oxysporum_fsp_lycopersici/4287/Ma_et_al_2010/fusarium_oxysporum_f._sp._lycopersici_mitochondrion_2_contigs.fasta); do
Organism=$(echo "$Genome" | rev | cut -d '/' -f4 | rev)
Strain=$(echo "$Genome" | rev | cut -d '/' -f3 | rev)
OutDir=analysis/mimps/$Organism/"$Strain"
mkdir -p "$OutDir"
"$ProgDir"/mimp_finder.pl "$Genome" "$OutDir"/"$Strain"_mimps.fa "$OutDir"/"$Strain"_mimps.gff3 > "$OutDir"/"$Strain"_mimps.log
done
```

### B) Augustus genes flanking Mimps


<!-- ```bash
	ProgDir=~/git_repos/emr_repos/tools/pathogen/mimp_finder
	for Mimps in $(ls -d analysis/mimps/F.*/*/*_mimps.gff3); do
		Organism=$(echo "$Mimps" | rev | cut -d '/' -f3 | rev)
		Strain=$(echo "$Mimps" | rev | cut -f2 -d '/' | rev)
		MimpDir=$(dirname $Mimps)
		echo "$Mimps"
		"$ProgDir"/gffexpander.pl + 1500 "$Mimps" > "$MimpDir"/"$Strain"_mimps_expanded.gff3
		StrainAugModels=$(ls gene_pred/augustus/$Organism/"$Strain"/*_augustus_preds.gtf)
		StrainORFModels=$(ls gene_pred/ORF_finder/$Organism/"$Strain"/*_ORF_mod.gff)
		bedtools intersect -s -a "$StrainAugModels"  -b "$MimpDir"/"$Strain"_mimps_expanded.gff3 > "$MimpDir"/"$Strain"_mimps_intersected_Aug_genes.gff
		cat "$MimpDir"/"$Strain"_mimps_intersected_Aug_genes.gff | grep 'gene' | rev | cut -f1 -d '=' | rev | sort | uniq > "$MimpDir"/"$Strain"_mimps_intersected_Aug_genes_names.txt
		bedtools intersect -s -a "$StrainORFModels"  -b "$MimpDir"/"$Strain"_mimps_expanded.gff3 > "$MimpDir"/"$Strain"_mimps_intersected_ORF_genes.gff
		cat "$MimpDir"/"$Strain"_mimps_intersected_ORF_genes.gff | grep 'gene' | rev | cut -f1 -d '=' | rev | sort | uniq > "$MimpDir"/"$Strain"_mimps_intersected_ORF_genes_names.txt
		printf "Augustus genes intersected:\t"
		cat "$MimpDir"/"$Strain"_mimps_intersected_Aug_genes_names.txt | wc -l
		printf "ORF fragments intersected:\t"
		cat "$MimpDir"/"$Strain"_mimps_intersected_ORF_genes_names.txt | wc -l
		echo ""
	done
```

Results were as follows:
```
	analysis/mimps/F.oxysporum_fsp_cepae/125/125_mimps.gff3
	Augustus genes intersected:	21
	ORF fragments intersected:	285

	analysis/mimps/F.oxysporum_fsp_cepae/55/55_mimps.gff3
	Augustus genes intersected:	15
	ORF fragments intersected:	254

	analysis/mimps/F.oxysporum_fsp_cepae/A23/A23_mimps.gff3
	Augustus genes intersected:	20
	ORF fragments intersected:	307

	analysis/mimps/F.oxysporum_fsp_cepae/A28/A28_mimps.gff3
	Augustus genes intersected:	5
	ORF fragments intersected:	99

	analysis/mimps/F.oxysporum_fsp_cepae/D2/D2_mimps.gff3
	Augustus genes intersected:	0
	ORF fragments intersected:	87

	analysis/mimps/F.oxysporum_fsp_cepae/Fus2/Fus2_mimps.gff3
	Augustus genes intersected:	20
	ORF fragments intersected:	260

	analysis/mimps/F.oxysporum_fsp_cepae/HB17/HB17_mimps.gff3
	Augustus genes intersected:	29
	ORF fragments intersected:	328

	analysis/mimps/F.oxysporum_fsp_cepae/PG/PG_mimps.gff3
	Augustus genes intersected:	5
	ORF fragments intersected:	170

	analysis/mimps/F.oxysporum_fsp_narcissi/N139/N139_mimps.gff3
	Augustus genes intersected:	15
	ORF fragments intersected:	395

	analysis/mimps/F.oxysporum_fsp_pisi/PG18/PG18_mimps.gff3
	Augustus genes intersected:	5
	ORF fragments intersected:	155

	analysis/mimps/F.oxysporum_fsp_pisi/PG3/PG3_mimps.gff3
	Augustus genes intersected:	3
	ORF fragments intersected:	146

	analysis/mimps/F.proliferatum/A8/A8_mimps.gff3
	Augustus genes intersected:	2
	ORF fragments intersected:	20
``` -->

```bash
	ProgDir=~/git_repos/emr_repos/tools/pathogen/mimp_finder
	for Mimps in $(ls -d analysis/mimps/F.*/*/*_mimps.gff3); do
		Organism=$(echo "$Mimps" | rev | cut -d '/' -f3 | rev)
		Strain=$(echo "$Mimps" | rev | cut -f2 -d '/' | rev)
		OutDir=analysis/mimps_+-2000bp/$Organism/$Strain
		mkdir -p $OutDir
		MimpDir=$(dirname $Mimps)
		echo "$Mimps"
		"$ProgDir"/gffexpander.pl +- 2000 "$Mimps" > "$OutDir"/"$Strain"_mimps_2000bp_expanded.gff3
		StrainAugModels=$(ls gene_pred/augustus/$Organism/"$Strain"/*_augustus_preds.gtf)
		StrainORFModels=$(ls gene_pred/ORF_finder/$Organism/"$Strain"/*_ORF_mod.gff)
		bedtools intersect -a "$StrainAugModels"  -b "$OutDir"/"$Strain"_mimps_2000bp_expanded.gff3 > "$OutDir"/"$Strain"_mimps_intersected_Aug_genes.gff
		#bedtools intersect -s -a "$StrainAugModels"  -b "$OutDir"/"$Strain"_mimps_2000bp_expanded.gff3 > "$OutDir"/"$Strain"_mimps_intersected_Aug_genes.gff
		cat "$OutDir"/"$Strain"_mimps_intersected_Aug_genes.gff | grep 'gene' | rev | cut -f1 -d '=' | rev | sort | uniq > "$OutDir"/"$Strain"_mimps_intersected_Aug_genes_names.txt
		bedtools intersect -s -a "$StrainORFModels"  -b "$OutDir"/"$Strain"_mimps_2000bp_expanded.gff3 > "$OutDir"/"$Strain"_mimps_intersected_ORF_genes.gff
		cat "$OutDir"/"$Strain"_mimps_intersected_ORF_genes.gff | grep 'gene' | rev | cut -f1 -d '=' | rev | sort | uniq > "$OutDir"/"$Strain"_mimps_intersected_ORF_genes_names.txt
		printf "Augustus genes intersected:\t"
		cat "$OutDir"/"$Strain"_mimps_intersected_Aug_genes_names.txt | wc -l
		printf "ORF fragments intersected:\t"
		cat "$OutDir"/"$Strain"_mimps_intersected_ORF_genes_names.txt | wc -l
		echo ""
	done
```

```
	analysis/mimps/F.oxysporum_fsp_cepae/125/125_mimps.gff3
	Augustus genes intersected:	62
	ORF fragments intersected:	767

	analysis/mimps/F.oxysporum_fsp_cepae/55/55_mimps.gff3
	Augustus genes intersected:	47
	ORF fragments intersected:	757

	analysis/mimps/F.oxysporum_fsp_cepae/A23/A23_mimps.gff3
	Augustus genes intersected:	64
	ORF fragments intersected:	755

	analysis/mimps/F.oxysporum_fsp_cepae/A28/A28_mimps.gff3
	Augustus genes intersected:	17
	ORF fragments intersected:	393

	analysis/mimps/F.oxysporum_fsp_cepae/D2/D2_mimps.gff3
	Augustus genes intersected:	11
	ORF fragments intersected:	326

	analysis/mimps/F.oxysporum_fsp_cepae/Fus2/Fus2_mimps.gff3
	Augustus genes intersected:	60
	ORF fragments intersected:	654

	analysis/mimps/F.oxysporum_fsp_cepae/HB17/HB17_mimps.gff3
	Augustus genes intersected:	68
	ORF fragments intersected:	830

	analysis/mimps/F.oxysporum_fsp_cepae/PG/PG_mimps.gff3
	Augustus genes intersected:	16
	ORF fragments intersected:	415

	analysis/mimps/F.oxysporum_fsp_narcissi/N139/N139_mimps.gff3
	Augustus genes intersected:	39
	ORF fragments intersected:	951

	analysis/mimps/F.oxysporum_fsp_pisi/PG18/PG18_mimps.gff3
	Augustus genes intersected:	20
	ORF fragments intersected:	470

	analysis/mimps/F.oxysporum_fsp_pisi/PG3/PG3_mimps.gff3
	Augustus genes intersected:	19
	ORF fragments intersected:	425

	analysis/mimps/F.proliferatum/A8/A8_mimps.gff3
	Augustus genes intersected:	5
	ORF fragments intersected:	58
```

### for FoL reference genome

```bash
ProgDir=~/git_repos/emr_repos/tools/pathogen/mimp_finder
for Mimps in $(ls -d analysis/mimps/F.oxysporum_fsp_lycopersici/4287/*_mimps.gff3); do
Organism=$(echo "$Mimps" | rev | cut -d '/' -f3 | rev)
Strain=$(echo "$Mimps" | rev | cut -f2 -d '/' | rev)
OutDir=analysis/mimps_+-2000bp/$Organism/$Strain
mkdir -p $OutDir
MimpDir=$(dirname $Mimps)
echo "$Mimps"
"$ProgDir"/gffexpander.pl +- 2000 "$Mimps" > "$OutDir"/"$Strain"_mimps_2000bp_expanded.gff3
StrainAugModels=assembly/external_group/F.oxysporum_fsp_lycopersici/4287/Fusox1/Fusox1.filtered_proteins.ExternalModels.gff3
bedtools intersect -a "$StrainAugModels"  -b "$OutDir"/"$Strain"_mimps_2000bp_expanded.gff3 > "$OutDir"/"$Strain"_mimps_intersected_Aug_genes.gff
cat "$OutDir"/"$Strain"_mimps_intersected_Aug_genes.gff | grep 'gene' | rev | cut -f1 -d '=' | rev | sort | uniq > "$OutDir"/"$Strain"_mimps_intersected_Aug_genes_names.txt
bedtools intersect -s -a "$StrainORFModels"  -b "$OutDir"/"$Strain"_mimps_2000bp_expanded.gff3 > 		printf "Augustus genes intersected:\t"
cat "$OutDir"/"$Strain"_mimps_intersected_Aug_genes_names.txt | wc -l
echo ""
done
```

## Merging effector evidence:

```bash
	for Strain in $(ls -d gene_pred/ORF_finder/F.*/* | rev | cut -f1 -d '/' | rev); do
		Queue=$(qstat | grep 'merge_fus' | sed -E 's/ +/ /g' |  cut -f5 -d ' ' | grep -c 'qw' | sed 's/ //g')
		while [[ $Queue -ge 1 ]]; do
			sleep 10
			printf "."
			Queue=$(qstat | grep 'merge_fus' | sed -E 's/ +/ /g' |  cut -f5 -d ' ' | grep -c 'qw')
		done
		# Node1Jobs=$(qstat | grep 'merge' | sed -E 's/ +/ /g' | cut -f 8 -d ' ' | sort | grep -c 'blacklace01' | sed 's/ //g')
		Node2Jobs=$(qstat | grep 'merge' | sed -E 's/ +/ /g' | cut -f 8 -d ' ' | sort | grep -c 'blacklace02' | sed 's/ //g')
		# Node3Jobs=$(qstat | grep 'merge' | sed -E 's/ +/ /g' | cut -f 8 -d ' ' | sort | grep -c 'blacklace03' | sed 's/ //g')
		Node4Jobs=$(qstat | grep 'merge' | sed -E 's/ +/ /g' | cut -f 8 -d ' ' | sort | grep -c 'blacklace04' | sed 's/ //g')
		#Node5Jobs=$(qstat | grep 'merge_fus' | sed -E 's/ +/ /g' | cut -f 8 -d ' ' | sort | grep -c 'blacklace05' | sed 's/ //g')
		Node6Jobs=$(qstat | grep 'merge' | sed -E 's/ +/ /g' | cut -f 8 -d ' ' | sort | grep -c 'blacklace06' | sed 's/ //g')
		Node11Jobs=$(qstat | grep 'merge' | sed -E 's/ +/ /g' | cut -f 8 -d ' ' | sort | grep -c 'blacklace11' | sed 's/ //g')
		# while [[ $Node1Jobs -ge 1  && $Node2Jobs -ge 1 && $Node3Jobs -ge 1 && $Node4Jobs -ge 1 ]]; do
		while [[ $Node2Jobs -ge 1 && $Node4Jobs -ge 1 && $Node6Jobs -ge 1 && $Node11Jobs -ge 1 ]]; do
			#while [ [ $Node1Jobs -ge 1 ] && [ $Node2Jobs -ge 1 ] && [ $Node3Jobs -ge 1 ] && [ $Node4Jobs -ge 1 ] && [ $Node5Jobs -ge 1 ] ]; do
			sleep 10
			printf "."
			# Node1Jobs=$(qstat | grep 'merge' | sed -E 's/ +/ /g' | cut -f 8 -d ' ' | sort | grep -c 'blacklace01' | sed 's/ //g')
			Node2Jobs=$(qstat | grep 'merge' | sed -E 's/ +/ /g' | cut -f 8 -d ' ' | sort | grep -c 'blacklace02' | sed 's/ //g')
			# Node3Jobs=$(qstat | grep 'merge' | sed -E 's/ +/ /g' | cut -f 8 -d ' ' | sort | grep -c 'blacklace03' | sed 's/ //g')
			Node4Jobs=$(qstat | grep 'merge' | sed -E 's/ +/ /g' | cut -f 8 -d ' ' | sort | grep -c 'blacklace04' | sed 's/ //g')
			#Node5Jobs=$(qstat | grep 'merge_fus' | sed -E 's/ +/ /g' | cut -f 8 -d ' ' | sort | grep -c 'blacklace05' | sed 's/ //g')
			Node6Jobs=$(qstat | grep 'merge' | sed -E 's/ +/ /g' | cut -f 8 -d ' ' | sort | grep -c 'blacklace06' | sed 's/ //g')
			Node11Jobs=$(qstat | grep 'merge' | sed -E 's/ +/ /g' | cut -f 8 -d ' ' | sort | grep -c 'blacklace11' | sed 's/ //g')
		done
		ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/merge_gff;
		echo $Strain;
		Orf_Gff=$(ls gene_pred/ORF_finder/*/$Strain/"$Strain"_ORF.gff);
		Aug_Gff=$(ls gene_pred/augustus/*/"$Strain"/"$Strain"_augustus_preds.gtf | grep -v 'old');
		echo $Orf_Gff;
		echo $Aug_Gff;
		# if [[ $Node1Jobs -lt 1 ]]; then
		# 	qsub -l h=blacklace01.blacklace $ProgDir/merge_fus_effectors.sh $Orf_Gff $Aug_Gff;
		# elif [[ $Node2Jobs -lt 1 ]]; then
		if [[ $Node2Jobs -lt 1 ]]; then
			qsub -l h=blacklace02.blacklace $ProgDir/merge_fus_effectors.sh $Orf_Gff $Aug_Gff;
		# elif [[ $Node3Jobs -lt 1 ]]; then
		# 	qsub -l h=blacklace03.blacklace $ProgDir/merge_fus_effectors.sh $Orf_Gff $Aug_Gff;
		elif [[ $Node4Jobs -lt 1 ]]; then
			qsub -l h=blacklace04.blacklace $ProgDir/merge_fus_effectors.sh $Orf_Gff $Aug_Gff;
			# elif [[ $Node5Jobs -lt 1 ]]; then
			# qsub -l h=blacklace05.blacklace $ProgDir/merge_fus_effectors.sh $Orf_Gff $Aug_Gff;
		elif [[ $Node6Jobs -lt 1 ]]; then
			qsub -l h=blacklace06.blacklace $ProgDir/merge_fus_effectors.sh $Orf_Gff $Aug_Gff;
		elif [[ $Node11Jobs -lt 1 ]]; then
			qsub -l h=blacklace11.blacklace $ProgDir/merge_fus_effectors.sh $Orf_Gff $Aug_Gff;
		else
			echo "Error something has jumped the queue"
		fi
	done
```

# 4. Genomic analysis

## 4.1 Identifcation of protospacers

To facilitate CriprCas editing of the Fusarium oxysporum genome target sites
known as protospacers were identified.

This was done using a published program OPTIMuS as well as a parser script. The
commands to do this were as follows:

```bash
	for GeneSeq in $(ls gene_pred/augustus/F.*/*/*_augustus_preds.codingseq | grep -v 'old'); do
	  Organism=$(echo $GeneSeq | rev | cut -f3 -d '/' | rev)
	  Strain=$(echo $GeneSeq | rev | cut -f2 -d '/' | rev)
		echo $Organism
		echo $Strain
	  ProgDir=~/git_repos/emr_repos/scripts/fusarium_venenatum/OPTIMus
	  OutDir=analysis/protospacers/$Organism/$Strain
	  mkdir -p $OutDir
	  $ProgDir/OPTIMuS_EMR.pl $GeneSeq "threshold" 1 > $OutDir/"$Strain"_protospacer_sites.txt
	  $ProgDir/Optimus2csv.py --inp $OutDir/"$Strain"_protospacer_sites.txt  --out $OutDir/"$Strain"_protospacer_by_gene.csv
	done
```


## 4.2 Orthology

Orthomcl was used to identify orthologous groups between Fusarium spp. genomes

Genomes were grouped by subspecies and orthology was determined within each
subspecies group. Orthology was also determined between subspecies groups.

| Pathogenic | non-pathogenic | Intermediate |
| ---------- | -------------- | -------------|
| 125        | A28            | 55           |
| A23        | D2             |              |
| Fus2       | PG             |              |




### 4.2.a) Orthology between pathogenic isolates

The Commands used to run this analysis are shown in
pathogen/orthology/F.oxysporum_fsp.cepae_pathogenic_orthology.md


### 4.2.b) Orthology between non-pathogenic isolates

The Commands used to run this analysis are shown in
pathogen/orthology/F.oxysporum_fsp.cepae_non-pathogenic_orthology.md


### 4.2.c) Orthology between pathogenic and non-pathogenic isolates

The Commands used to run this analysis are shown in
pathogen/orthology/F.oxysporum_fsp.cepae_pathogen_vs_non-pathogen_orthology.md


### 4.2.d) Orthology between all isolates

The Commands used to run this analysis are shown in
pathogen/orthology/F.oxysporum_fsp.cepae_isolates.md


## 5. BLAST Searches

## 5.1 Identifying SIX genes

Protein sequence of previously characterised SIX genes used to BLAST against
assemlies.

```bash
	‘grep -A1 '>' analysis/blast_homology/six_genes/six-appended.fa | cut -d '|' -f5 > analysis/blast_homology/six_genes/six-appended_parsed.fa
	vi analysis/blast_homology/six_genes/six-appended_parsed_tmp.fa # Removed initial ' ' from line and added '>'
	sed 's/\s/_/g' analysis/blast_homology/six_genes/six-appended_parsed_tmp.fa > analysis/blast_homology/six_genes/six-appended_parsed.fa

	qsub /home/armita/git_repos/emr_repos/tools/pathogen/blast/blast_pipe.sh analysis/blast_homology/six_genes/six-appended_parsed.fa dna assembly/velvet/F.oxysporum_fsp_cepae/125/125_assembly.41/sorted_contigs.fa
	qsub /home/armita/git_repos/emr_repos/tools/pathogen/blast/blast_pipe.sh analysis/blast_homology/six_genes/six-appended_parsed.fa dna assembly/velvet/F.oxysporum_fsp_cepae/55/55_assembly.41/sorted_contigs.fa
	qsub /home/armita/git_repos/emr_repos/tools/pathogen/blast/blast_pipe.sh analysis/blast_homology/six_genes/six-appended_parsed.fa dna assembly/velvet/F.oxysporum_fsp_cepae/A23/A23_assembly.41/sorted_contigs.fa
	qsub /home/armita/git_repos/emr_repos/tools/pathogen/blast/blast_pipe.sh analysis/blast_homology/six_genes/six-appended_parsed.fa dna assembly/velvet/F.oxysporum_fsp_cepae/A28/A28_assembly.41/sorted_contigs.fa
	qsub /home/armita/git_repos/emr_repos/tools/pathogen/blast/blast_pipe.sh analysis/blast_homology/six_genes/six-appended_parsed.fa dna assembly/velvet/F.oxysporum_fsp_cepae/D2/D2_assembly.41/sorted_contigs.fa
	qsub /home/armita/git_repos/emr_repos/tools/pathogen/blast/blast_pipe.sh analysis/blast_homology/six_genes/six-appended_parsed.fa dna assembly/velvet/F.oxysporum_fsp_cepae/Fus2/Fus2_assembly.41/sorted_contigs.fa
	qsub /home/armita/git_repos/emr_repos/tools/pathogen/blast/blast_pipe.sh analysis/blast_homology/six_genes/six-appended_parsed.fa dna assembly/velvet/F.oxysporum_fsp_cepae/HB17/HB17_assembly.41/sorted_contigs.fa
	qsub /home/armita/git_repos/emr_repos/tools/pathogen/blast/blast_pipe.sh analysis/blast_homology/six_genes/six-appended_parsed.fa dna assembly/velvet/F.oxysporum_fsp_cepae/PG/PG_assembly.41/sorted_contigs.fa

```

## 5.2 Identifying PHIbase homologs

The PHIbase database was searched agasinst the assembled genomes using tBLASTx.

```bash
	for Assembly in $(ls repeat_masked/F.*/*/*/*_contigs_unmasked.fa); do
		qsub /home/armita/git_repos/emr_repos/tools/pathogen/blast/blast_pipe.sh analysis/blast_homology/PHIbase/PHI_36_accessions.fa protein $Assembly
	done
```

following blasting PHIbase to the genome, the hits were filtered by effect on
virulence.

First the a tab seperated file was made in the clusters core directory containing
PHIbase. These commands were run as part of previous projects but have been
included here for completeness.

```bash
	PhibaseDir=/home/groups/harrisonlab/phibase/v3.8
	printf "header\n" > $PhibaseDir/PHI_headers.csv
	cat $PhibaseDir/PHI_accessions.fa | grep '>' | cut -f1 | sed 's/>//g' | sed 's/\r//g' >> $PhibaseDir/PHI_headers.csv
	printf "effect\n" > .$PhibaseDir/PHI_virulence.csv
	cat $PhibaseDir/PHI_accessions.fa | grep '>' | cut -f1 | sed 's/>//g' | rev | cut -f1 -d '|' | rev  >> $PhibaseDir/PHI_virulence.csv
```


```bash
	PhibaseDir=/home/groups/harrisonlab/phibase/v3.8
	PhibaseHeaders=$PhibaseDir/PHI_headers.csv
	PhibaseVirulence=$PhibaseDir/PHI_virulence.csv
	for BlastCSV in $(ls analysis/blast_homology/F*/*/*_PHI_36_accessions.fa_homologs.csv); do
		Strain=$(echo $BlastCSV | rev | cut -f2 -d'/' | rev)
		echo "$Strain"
		OutDir=$(dirname $BlastCSV)
		paste -d '\t' $PhibaseHeaders $PhibaseVirulence $BlastCSV | cut -f-3,1185- > $OutDir/"$Strain"_PHIbase_virulence.csv
		cat $OutDir/"$Strain"_PHIbase_virulence.csv | grep 'NODE_' | cut -f2 | sort | uniq -c | tee $OutDir/"$Strain"_PHIbase_virulence.txt
	done
```
results were:

```
	125
	      3 Chemistry target
	      1 Effector (plant avirulence determinant)
	      1 Increased virulence (Hypervirulence)
	     15 Lethal
	     66 Loss of pathogenicity
	     14 Mixed outcome
	      2 reduced virulence
	    131 Reduced virulence
	     87 Unaffected pathogenicity
	55
	      3 Chemistry target
	      1 Effector (plant avirulence determinant)
	      1 Increased virulence (Hypervirulence)
	     14 Lethal
	     58 Loss of pathogenicity
	     12 Mixed outcome
	      2 reduced virulence
	    121 Reduced virulence
	     83 Unaffected pathogenicity
	A23
	      3 Chemistry target
	      1 Effector (plant avirulence determinant)
	      1 Increased virulence (Hypervirulence)
	     14 Lethal
	     60 Loss of pathogenicity
	     12 Mixed outcome
	      2 reduced virulence
	    123 Reduced virulence
	     79 Unaffected pathogenicity
	A28
	      3 Chemistry target
	      1 Effector (plant avirulence determinant)
	      2 Increased virulence (Hypervirulence)
	     15 Lethal
	     58 Loss of pathogenicity
	     13 Mixed outcome
	      2 reduced virulence
	    125 Reduced virulence
	     79 Unaffected pathogenicity
	D2
	      3 Chemistry target
	      1 Effector (plant avirulence determinant)
	      2 Increased virulence (Hypervirulence)
	     16 Lethal
	     62 Loss of pathogenicity
	     12 Mixed outcome
	      2 reduced virulence
	    123 Reduced virulence
	     81 Unaffected pathogenicity
	Fus2
	      3 Chemistry target
	      1 Effector (plant avirulence determinant)
	      1 Increased virulence (Hypervirulence)
	     15 Lethal
	     59 Loss of pathogenicity
	     12 Mixed outcome
	      2 reduced virulence
	    120 Reduced virulence
	     79 Unaffected pathogenicity
	HB17
	      3 Chemistry target
	      1 Effector (plant avirulence determinant)
	      1 Increased virulence (Hypervirulence)
	     14 Lethal
	     60 Loss of pathogenicity
	     13 Mixed outcome
	      2 reduced virulence
	    123 Reduced virulence
	     84 Unaffected pathogenicity
	PG
	      3 Chemistry target
	      1 Effector (plant avirulence determinant)
	      2 Increased virulence (Hypervirulence)
	     15 Lethal
	     61 Loss of pathogenicity
	     13 Mixed outcome
	      2 reduced virulence
	    126 Reduced virulence
	     78 Unaffected pathogenicity
	N139
	      4 Chemistry target
	      1 Effector (plant avirulence determinant)
	      2 Increased virulence (Hypervirulence)
	     16 Lethal
	     72 Loss of pathogenicity
	     13 Mixed outcome
	      2 reduced virulence
	    149 Reduced virulence
	     93 Unaffected pathogenicity
	PG18
	      4 Chemistry target
	      1 Effector (plant avirulence determinant)
	      2 Increased virulence (Hypervirulence)
	     17 Lethal
	     72 Loss of pathogenicity
	     14 Mixed outcome
	      2 reduced virulence
	    153 Reduced virulence
	     97 Unaffected pathogenicity
	PG3
	      4 Chemistry target
	      1 Effector (plant avirulence determinant)
	      2 Increased virulence (Hypervirulence)
	     17 Lethal
	     65 Loss of pathogenicity
	     13 Mixed outcome
	      2 reduced virulence
	    139 Reduced virulence
	     93 Unaffected pathogenicity
	A8
	      4 Chemistry target
	      1 Effector (plant avirulence determinant)
	      4 Increased virulence (Hypervirulence)
	     15 Lethal
	     71 Loss of pathogenicity
	     14 Mixed outcome
	      2 reduced virulence
	    152 Reduced virulence
	     82 Unaffected pathogenicity

```

The analysis was also performed by blasting the predicted proteins against the
PHIbase database:

The PHIbase database was searched agasinst the assembled genomes using tBLASTx.

```bash
	for Proteins in $(ls gene_pred/augustus/F.oxysporum_fsp_cepae/*/*_augustus_preds.aa); do
		qsub /home/armita/git_repos/emr_repos/tools/pathogen/blast/blast_pipe.sh $Proteins protein ../../phibase/v3.8/PHI_accessions.fa
	done
```
