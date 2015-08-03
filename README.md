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

### B) Augustus genes flanking Mimps


```bash
	ProgDir=~/git_repos/emr_repos/tools/pathogen/mimp_finder
	for Mimps in $(ls -d analysis/mimps/F.*/*/*_mimps.gff3); do
		Organism=$(echo "$Mimps" | rev | cut -d '/' -f3 | rev)
		Strain=$(echo "$Mimps" | rev | cut -f2 -d '/' | rev)
		MimpDir=$(dirname $Mimps)
		echo "$Mimps"
		"$ProgDir"/gffexpander.pl + 1000 "$Mimps" > "$MimpDir"/"$Strain"_mimps_expanded.gff3
		StrainModels=$(ls gene_pred/augustus/$Organism/"$Strain"/*_augustus_preds.gtf)
		bedtools intersect -s -a "$StrainModels"  -b "$MimpDir"/"$Strain"_mimps_expanded.gff3 > "$MimpDir"/"$Strain"_mimps_intersected_genes.gff
		cat "$MimpDir"/*_mimps_intersected_genes.gff | grep 'gene' | rev | cut -f1 -d '=' | rev | sort | uniq | wc -l
	done
```

Results were as follows:
```
	analysis/mimps/F.oxysporum_fsp_cepae/125/125_mimps.gff3  
	18  
	analysis/mimps/F.oxysporum_fsp_cepae/55/55_mimps.gff3  
	12  
	analysis/mimps/F.oxysporum_fsp_cepae/A23/A23_mimps.gff3  
	17  
	analysis/mimps/F.oxysporum_fsp_cepae/A28/A28_mimps.gff3  
	3  
	analysis/mimps/F.oxysporum_fsp_cepae/D2/D2_mimps.gff3  
	0  
	analysis/mimps/F.oxysporum_fsp_cepae/Fus2/Fus2_mimps.gff3  
	16  
	analysis/mimps/F.oxysporum_fsp_cepae/HB17/HB17_mimps.gff3  
	23  
	analysis/mimps/F.oxysporum_fsp_cepae/PG/PG_mimps.gff3  
	3  
	analysis/mimps/F.oxysporum_fsp_narcissi/N139/N139_mimps.gff3  
	9  
	analysis/mimps/F.oxysporum_fsp_pisi/PG18/PG18_mimps.gff3  
	3  
	analysis/mimps/F.oxysporum_fsp_pisi/PG3/PG3_mimps.gff3  
	3  
	analysis/mimps/F.proliferatum/A8/A8_mimps.gff3  
	0  
```
<!--

#Genomic analysis


The first analysis was based upon BLAST searches for genes known to be involved in toxin production

#BLAST Searches
```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast/
	Query=analysis/blast_homology/CDC_genes/A.alternata_CDC_genes.fa
	BestAss675=assembly/velvet/A.alternata_ssp._arborescens/675/A.alternata_ssp._arborescens_675_69/sorted_contigs.fa
	BestAss970013=assembly/velvet/A.alternata_ssp._arborescens/97.0013/A.alternata_ssp._arborescens_97.0013_59/sorted_contigs.fa
	BestAss970016=assembly/velvet/A.alternata_ssp._arborescens/97.0016/A.alternata_ssp._arborescens_97.0016_77/sorted_contigs.fa
	BestAss650=assembly/velvet/A.alternata_ssp._gaisen/650/A.alternata_ssp._gaisen_650_67/sorted_contigs.fa
	BestAss1082=assembly/velvet/A.alternata_ssp._tenuissima/1082/A.alternata_ssp._tenuissima_1082_49/sorted_contigs.fa
	BestAss1164=assembly/velvet/A.alternata_ssp._tenuissima/1164/A.alternata_ssp._tenuissima_1164_67/sorted_contigs.fa
	BestAss1166=assembly/velvet/A.alternata_ssp._tenuissima/1166/A.alternata_ssp._tenuissima_1166_43/sorted_contigs.fa
	BestAss1177=assembly/velvet/A.alternata_ssp._tenuissima/1177/A.alternata_ssp._tenuissima_1177_63/sorted_contigs.fa
	BestAss24350=assembly/velvet/A.alternata_ssp._tenuissima/24350/A.alternata_ssp._tenuissima_24350_63/sorted_contigs.fa
	BestAss635=assembly/velvet/A.alternata_ssp._tenuissima/635/A.alternata_ssp._tenuissima_635_59/sorted_contigs.fa
	BestAss648=assembly/velvet/A.alternata_ssp._tenuissima/648/A.alternata_ssp._tenuissima_648_45/sorted_contigs.fa
	BestAss743=assembly/velvet/A.alternata_ssp._tenuissima/743/A.alternata_ssp._tenuissima_743_69/sorted_contigs.fa
	qsub $ProgDir/blast_pipe.sh $Query dna $BestAss675
	qsub $ProgDir/blast_pipe.sh $Query dna $BestAss970013
	qsub $ProgDir/blast_pipe.sh $Query dna $BestAss970016
	qsub $ProgDir/blast_pipe.sh $Query dna $BestAss650
	qsub $ProgDir/blast_pipe.sh $Query dna $BestAss1082
	qsub $ProgDir/blast_pipe.sh $Query dna $BestAss1164
	qsub $ProgDir/blast_pipe.sh $Query dna $BestAss1166
	qsub $ProgDir/blast_pipe.sh $Query dna $BestAss1177
	qsub $ProgDir/blast_pipe.sh $Query dna $BestAss24350
	qsub $ProgDir/blast_pipe.sh $Query dna $BestAss635
	qsub $ProgDir/blast_pipe.sh $Query dna $BestAss648
	qsub $ProgDir/blast_pipe.sh $Query dna $BestAss743
```
BLAST search results were summarised into a presence/absence table

The presence /absence table determines presence if a hit is present and the alignment represents >50% of the query sequence.
This thresholding means that some hits have not been summarised including AMT11, AMT15 and ALT1.
```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast/
	InFiles=$(ls analysis/blast_homology/A.alternata_ssp._*/*/*_A.alternata_CDC_genes.fa_homologs.csv | paste -s -d ' ')
	echo $InFiles
	$ProgDir/blast_differentials.pl $InFiles
	mv *.csv analysis/blast_homology/CDC_genes/.
```






#CDC Assembly

Raw reads were aligned against assembled genomes to identify contigs that were unique to a isolate or clade
```bash
	for Pathz in $(ls -d qc_dna/paired/A.alternata_ssp._*/*); do  
		Strain=$(echo $Pathz | cut -d '/' -f4)
		echo "using reads for $Strain"
		ProgPath=/home/armita/git_repos/emr_repos/tools/pathogen/lineage_specific_regions
		F_IN=$(ls $Pathz/F/*.fastq.gz)
		R_IN=$(ls $Pathz/R/*.fastq.gz)
		for Assemblyz in $(ls repeat_masked/A.alternata_ssp._*/*/*/*_contigs_unmasked.fa); do
			basename $Assemblyz
			qsub "$ProgPath"/bowtie2_alignment_pipe.sh $F_IN $R_IN $Assemblyz
		done
	done
```
A summary file was made from the alignment logs.
The percentage of reads aligning to each set of assembled contigs was determined.

```bash
	SummaryFile=analysis/ls_contigs/alignment_summaries.txt
	printf "" > "$SummaryFile"
	for OUTPUT in $(ls bowtie2_alignment_pipe.sh.e*); do
		ID=$(echo $OUTPUT | rev | cut -d 'e' -f1 | rev | less);
		cat bowtie2_alignment_pipe.sh.o"$ID" | grep -E "Trimmed .* reads .*/F/|Output files: " | sed -e 's/.*\/F\///g' | cut -f1 -d ')' | cut -f2 -d '"' >> "$SummaryFile";
		cat $OUTPUT >> "$SummaryFile";
		printf "\n" >> "$SummaryFile";
	done
	cat $SummaryFile | grep ' overall alignment rate' | cut -f1 -d ' ' | less
```
These are the reads percentages reported by bowtie but do not actually reflect the percentage unaligned reads where neither pair matched.

These were identified using SAM FLAGS to extract unaligned pairs.
```bash
	for Pathz in $(ls assembly/ls_contigs/A.alternata_ssp._*/*/vs_*/*_sorted.bam); do  
		OutFileF=$(echo $Pathz | sed 's/.bam/_unaligned_F.txt/g')
		OutFileR=$(echo $Pathz | sed 's/.bam/_unaligned_R.txt/g')
		OutFileSum=$(echo $Pathz | sed 's/.bam/_sum.txt/g')
		samtools view -f 77 "$Pathz" | cut -f1 > $OutFileF
		samtools view -f 141 "$Pathz" | cut -f1 > $OutFileR
		NoReads=$(samtools view -f 1 "$Pathz" | wc -l)
		NoReadsF=$(cat "$OutFileF" | wc -l)
		NoReadsR=$(cat "$OutFileR" | wc -l)
		printf "File\tNo. paired reads\tF reads unaligned\tR reads unaligned\n" > "$OutFileSum"
		printf "$Pathz\t$NoReads\t$NoReadsF\t$NoReadsR\n" >> "$OutFileSum"
	done

	for File in $(ls assembly/ls_contigs/A.alternata_ssp._*/*/vs_*/*_sum.txt); do
		cat $File | tail -n+2;
	done > analysis/ls_contigs/assembly_summaries2.txt
```

The number of reads aligning per bp of assembly was determined. Typical alignment values were 0.20 reads per bp. Contigs were detemined as unique to that alignment if they contained an average of 0 reads per bp.
The number of bp unique to each assembly were identified.

```bash
	mkdir -p analysis/ls_contigs
	Outfile=analysis/ls_contigs/ls_contig_size.csv
	printf "Reads" > $Outfile
	for Genome in $(ls -d assembly/ls_contigs/A.alternata_ssp._arborescens/675/vs_A.alternata_ssp._*); do
		NameG=$(basename "$Genome" | sed 's/vs_//g' | sed 's/_repmask_contigs//g')
		printf "\t$NameG" >> $Outfile
	done
	printf "\n" >> $Outfile
	for Reads in $(ls -d assembly/ls_contigs/A.alternata_ssp._*/*); do
		NameR=$(basename "$Reads")
		printf "$NameR" >> $Outfile
		for Genome in $(ls -d $Reads/*); do
			printf "\t" >> $Outfile
			cat $Genome/*_sorted_indexstats_coverage.csv | head -n-1 | grep -E -w "0.00$" | cut -f2 | awk '{s+=$1} END {printf s}' >> $Outfile
		done
		printf "\n" >> $Outfile
	done
```
This did not give clear results.

Contigs that had no reads align to them were identified.
 -->
