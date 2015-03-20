#!/bin/bash

# qc the raw reads before performing an alignment
for Strainz in $(ls -d raw_dna/paired/F.oxysporum_fsp_cepae/*); do 
	ProgPath=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
	F_IN=$(ls $Strainz/F/*.fastq.gz)
	R_IN=$(ls $Strainz/R/*.fastq.gz)
	IlluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/illumina_full_adapters.fa
	SeqType=dna
	qsub "$ProgPath"/rna_qc_fastq-mcf.sh "$F_IN" "$R_IN" "$IlluminaAdapters" "$SeqType"
done

mkdir -p repeat_masked_old/F.oxysporum_fsp_cepae/Fus2/Fus2_assembly.81_repmask
mv repeat_masked/F.oxysporum_fsp_cepae/Fus2/Fus2_assembly.81_repmask repeat_masked_old/F.oxysporum_fsp_cepae/Fus2/Fus2_assembly.81_repmask


# For each strain perform an alignment of raw reads against the genomes of all strains.
for Pathz in $(ls -d raw_dna/paired/F.oxysporum_fsp_cepae/*); do  
Strain=$(echo Strainz | cut -d '/' -f4)
echo "using reads for $Strain"
ProgPath=/home/armita/git_repos/emr_repos/tools/pathogen/lineage_specific_regions
F_IN=$(ls $Pathz/F/*.fastq.gz)
R_IN=$(ls $Pathz/R/*.fastq.gz)
for Assemblyz in $(ls repeat_masked/F.oxysporum_fsp_cepae/*/*/*_contigs_unmasked.fa); do
basename $Assemblyz
qsub "$ProgPath"/bowtie2_alignment_pipe.sh	$F_IN $R_IN $Assemblyz
done
done

SummaryFile=assembly/ls_contigs/alignment_summaries.txt
printf "" > "$SummaryFile"
for OUTPUT in $(ls bowtie2_alignment_pipe.sh.e*); do 
ID=$(echo $OUTPUT | rev | cut -d 'e' -f1 | rev | less); 
cat bowtie2_alignment_pipe.sh.o"$ID" | grep -E "Trimmed .* reads .*/F/|Output files: " | sed -e 's/.*\/F\///g' | cut -f1 -d ')' | cut -f2 -d '"' >> "$SummaryFile"; 
cat $OUTPUT >> "$SummaryFile"; 
printf "\n" >> "$SummaryFile"; 
done

Fo4287_Ec=assembly/external_group/F.oxysporum_fsp_lycopersici/4287/Ma_et_al_2010/F.oxysporum_fsp.lycopersici_4287.fasta 
Fo4287_Mt=assembly/external_group/F.oxysporum_fsp_lycopersici/4287/Ma_et_al_2010/fusarium_oxysporum_f._sp._lycopersici_mitochondrion_2_contigs.fasta
Fo4287_Out=assembly/external_group/F.oxysporum_fsp_lycopersici/4287/Ma_et_al_2010/F.oxysporum_fsp.lycopersici_4287_appended.fasta
cat $Fo4287_Ec $Fo4287_Mt > $Fo4287_Out

Fo47_Ec=assembly/external_group/F.oxysporum/fo47/broad/fusarium_oxysporum_fo47_1_supercontigs.fasta
Fo47_Linear=assembly/external_group/F.oxysporum/fo47/broad/fusarium_oxysporum_fo47_1_supercontigs_linear.fasta
cat $Fo47_Ec | sed -e 's/\(^>.*$\)/#\1#/' | tr -d "\r" | tr -d "\n" | sed -e 's/$/#/' | tr "#" "\n" | sed -e '/^$/d' > $Fo47_Linear

for Pathz in $(ls -d raw_dna/paired/F.oxysporum_fsp_cepae/*); do  
Strain=$(echo $Pathz | cut -d '/' -f4)
echo "using reads for $Strain"
ProgPath=/home/armita/git_repos/emr_repos/tools/pathogen/lineage_specific_regions
F_IN=$(ls $Pathz/F/*.fastq.gz)
R_IN=$(ls $Pathz/R/*.fastq.gz)
for Assemblyz in "$Fo4287_Out" "$Fo47_Linear"; do
basename $Assemblyz
qsub "$ProgPath"/bowtie2_alignment_pipe.sh	$F_IN $R_IN $Assemblyz
done
done
