These commands were used to identify highly expressed genes in the Fus2 assembly


# 1) QC

Perform qc of RNAseq timecourse data
```bash
  for FilePath in $(ls -d timecourse/experiment1/F.oxysporum_fsp_cepae/Fus2/*); do
    echo $FilePath
    FileF=$(ls $FilePath/F/*.gz)
    FileR=$(ls $FilePath/R/*.gz)
    IlluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa
    qsub /home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc/rna_qc_fastq-mcf.sh $FileF $FileR $IlluminaAdapters RNA
  done
```

# 2) Align reads vs. Fus2 genome
Alignments of RNAseq reads were made against the Fus2 Genome using tophat:

## 2.1) Alignment

```bash
  for FilePath in $(ls -d qc_rna/F.oxysporum_fsp_cepae/Fus2/*); do
    Genome=repeat_masked/F.oxysporum_fsp_cepae/Fus2/Fus2_combined_49/Fus2_combined_49_contigs_unmasked.fa
    FileF=$(ls $FilePath/F/*_trim.fq.gz)
    FileR=$(ls $FilePath/R/*_trim.fq.gz)
    qsub /home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq/tophat_alignment.sh $Genome $FileF $FileR
  done
  for FilePath in $(ls -d qc_rna/paired/F.oxysporum_fsp_cepae/*); do
    Genome=repeat_masked/F.oxysporum_fsp_cepae/Fus2/Fus2_combined_49/Fus2_combined_49_contigs_unmasked.fa
    FileF=$(ls $FilePath/F/*_trim.fq.gz)
    FileR=$(ls $FilePath/R/*_trim.fq.gz)
    qsub /home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq/tophat_alignment.sh $Genome $FileF $FileR
  done
```

## 2.2) Summarising alignments
Results were summarised using the following commands:

```bash
  for File in $(ls bowtie_alignment.sh.e*); do
    echo $File;
    Logfile=$(echo $File | sed "s/.sh.e/.sh.o/g");
    head -n 2 $Logfile;
    cat $File;
    echo "";  
  done
```

The text output of these alignments is stored in summised_output.txt

# 3) Move to timecourse directory

Data was copied from the alignment directory to a working directory for this
timecourse experiment.

```bash
  WorkDir=timecourse/v2_genes
  mkdir -p $WorkDir
  cp -r alignment/Fus2 $WorkDir/.
  # cd $WorkDir
```

# 4) Assemble transcripts

Cufflinks was used to assemble transcripts from reads aligned to the genome.

```bash
  for TimePoint in $(ls timecourse/v2_genes/* | grep -v 'merged' | grep -v 'quantified'); do
  echo $TimePoint;
  Alignment=$(ls timecourse/v2_genes/*/$TimePoint/accepted_hits.bam)
  echo $Alignment
  cufflinks -o timecourse/v2_genes/Fus2/$TimePoint/cufflinks -p 16 --max-intron-length 4000 $Alignment
  done
  # for Media in $(ls timecourse/v2_genes/F.oxysporum_fsp_cepae | grep -v 'merged' | grep -v 'quantified'); do
  # echo $Media;
  # Alignment=$(ls "timecourse/v2_genes/F.oxysporum_fsp_cepae/$Media/accepted_hits.bam")
  # cufflinks -o timecourse/v2_genes/Fus2/$Media/cufflinks -p 16 --max-intron-length 4000 $Alignment
  # done
```

# 5) Merge assembled transcripts

```bash
  ls timecourse/v2_genes/Fus2/*/cufflinks/transcripts.gtf | sort -g -k4 -t '/' > transcript_list.txt
  cuffmerge -o timecourse/v2_genes/Fus2/merged --num-threads 16 transcript_list.txt
  rm transcript_list.txt
```

# 6) quantify expression

```bash
  for TimePoint in $(ls timecourse/v2_genes/Fus2 | grep -v 'merged' | grep -v 'quantified' | sort -g -k4 -t '/'); do
    echo $TimePoint;
    Alignment=$(ls timecourse/v2_genes/*/$TimePoint/accepted_hits.bam)
    echo $Alignment
    cuffquant timecourse/v2_genes/Fus2/merged/merged.gtf -o timecourse/v2_genes/Fus2/quantified/$TimePoint -p 16 $Alignment
  done
```

```bash
  Alignments=$(ls timecourse/v2_genes/Fus2/quantified/*/abundances.cxb | sort -g -k5 -t '/')
  Labels=$(ls -d timecourse/v2_genes/Fus2/quantified/*/ | cut -f5 -d '/' | sort -g | tr '\n' ',')
  cuffdiff --time-series -o timecourse/v2_genes/Fus2/quantified -p 16 -L $Labels timecourse/v2_genes/Fus2/merged/merged.gtf $Alignments
```

```R
  source("http://bioconductor.org/biocLite.R")
  biocLite()
  biocLite("cummeRbund")

```

# 7) Use merged transcripts to train Augustus gene models
Instuctions were followed from:
http://transdecoder.github.io
http://jamg.sourceforge.net/tutorial.html (Step2c)

Note - this require the transdecoder utils scripts to be installed.

```bash
cufflinks_gtf_genome_to_cdna_fasta.pl $GffMerged $Assembly > tmp.fa
cufflinks_gtf_to_alignment_gff3.pl $GffMerged > tmp.gff
TransDecoder.LongOrfs -t tmp.fa
cdna_alignment_orf_to_genome_orf.pl tmp.fa.transdecoder_dir/longest_orfs.gff3 tmp.gff tmp.fa > tmp2.gff
gff2gbSmallDNA.pl tmp2.gff $Assembly 4000 tmp.gb


OutDir=gene_pred/training_augustus/Fusarium_oxysporum_fsp_cepae/Fus2
mkdir -p $OutDir
AUG_DIR=
PASA_DIR=$(which pasa | sed s%/pasa%%)
GffMerged=timecourse/v2_genes/Fus2/merged/merged.gtf
Assembly=repeat_masked/F.oxysporum_fsp_cepae/Fus2/Fus2_combined_49/Fus2_combined_49_contigs_unmasked.fa
Organism=F.oxysporum_fsp_cepae
Strain=Fus2
#------------------------------------------------------
# 		Step A		Convert .gff output to .gb format
#------------------------------------------------------
/home/armita/git_repos/emr_repos/tools/gene_prediction/augustus/gff3_2_auggff.pl $GffMerged > $OutDir/merged_mod.gff
MAX_FLANK_DNA=4000

gff2gbSmallDNA.pl $OutDir/merged_mod.gff $Assembly $MAX_FLANK_DNA $OutDir/"$ORGANISM"_"$STRAIN"_evidence.gb
# --good=$OutDir/complete_genes.txt

#------------------------------------------------------
# 		Step 5.		Extract a subset of the aligned reads to use as a test set
#------------------------------------------------------
# The genes.gb must be in the directory this command is being run from
randomSplit.pl $OutDir/"$ORGANISM"_"$STRAIN"_evidence.gb 100


#------------------------------------------------------
# 		Step 6.		Create a metafile for the new species
#------------------------------------------------------

new_species.pl --species="$ORGANISM"_"$STRAIN"

perl -pi -e "s%codingseq           off%codingseq           on%g" ~/prog/augustus-3.1/config/species/"$ORGANISM"_"$STRAIN"/"$ORGANISM"_"$STRAIN"_parameters.cfg
perl -pi -e "s%stopCodonExcludedFromCDS false%stopCodonExcludedFromCDS true%g" ~/prog/augustus-3.1/config/species/"$ORGANISM"_"$STRAIN"/"$ORGANISM"_"$STRAIN"_parameters.cfg
perl -pi -e "s%alternatives-from-evidence  false%alternatives-from-evidence  true%g" ~/prog/augustus-3.1/config/species/"$ORGANISM"_"$STRAIN"/"$ORGANISM"_"$STRAIN"_parameters.cfg


#------------------------------------------------------
# 		Step 7.		Train Augustus using aligned reads
#------------------------------------------------------

etraining --species="$ORGANISM"_"$STRAIN" "$ORGANISM"_"$STRAIN"_evidence.gb.train

augustus --species="$ORGANISM"_"$STRAIN" "$ORGANISM"_"$STRAIN"_evidence.gb.train | tee "$ORGANISM"_"$STRAIN"_sum.txt

```

<!--
```bash
  samtools view -bS ../4fus.sam > 24_hr_alignment.bam
  samtools sort 24_hr_alignment.bam 24_hr_alignment_sorted
  samtools index 24_hr_alignment_sorted.bam
  cp ../../repeat_masked/F.oxysporum_fsp_cepae/Fus2/Fus2_combined_49/Fus2_combined_49_contigs_hardmasked.fa
  samtools faidx Fus2_combined_49_contigs_hardmasked.fa
  samtools tview 24_hr_alignment_sorted.bam Fus2_combined_49_contigs_hardmasked.fa
  cp ../../gene_pred/augustus/F.oxysporum_fsp_cepae/Fus2/Fus2_aug_out.gff
  bedtools intersect -c -a Fus2_aug_out.gff -b 24_hr_alignment_sorted.bam > 24hFileRtersect3.bed
  cat 24hFileRtersect3.bed | grep -E -v "\s0$" | grep -w "gene" > Fus2_expressed_genes.bed
  cat Fus2_expressed_genes.bed | sort -r -n -k 10 | less
```
# This top lines of this file are:
# NODE_9700_length_34512_cov_86.442451    AUGUSTUS        gene    9014    11091   0.32    -       .       g13840  2629
# NODE_11428_length_3064_cov_2466.822754  AUGUSTUS        gene    1       2042    0.06    +       .       g15230  1759
# NODE_1097_length_595114_cov_83.826759   AUGUSTUS        gene    518067  521767  0.17    +       .       g2042   413
# NODE_11451_length_942_cov_2469.252686   AUGUSTUS        gene    248     613     0.72    -       .       g15334  372
# NODE_256_length_160675_cov_85.558762    AUGUSTUS        gene    86821   88592   0.36    +       .       g8812   309
# NODE_149_length_43877_cov_84.839462     AUGUSTUS        gene    20431   21783   0.96    +       .       g13469  295
# NODE_956_length_295543_cov_83.945419    AUGUSTUS        gene    250771  251290  0.99    -       .       g5661   241
# NODE_11459_length_659_cov_2500.394531   AUGUSTUS        gene    1       203     0.65    -       .       g15359  207
# NODE_724_length_250420_cov_83.169510    AUGUSTUS        gene    209007  210676  0.83    +       .       g6901   199
# NODE_8489_length_811399_cov_84.416985   AUGUSTUS        gene    296996  298062  0.94    +       .       g891    165
# NODE_4583_length_445765_cov_84.597122   AUGUSTUS        gene    274107  275696  1       +       .       g3043   134
# NODE_11477_length_312246_cov_83.660294  AUGUSTUS        gene    281911  283153  0.97    +       .       g5276   133
# NODE_4974_length_475916_cov_84.779022   AUGUSTUS        gene    137868  139358  0.39    -       .       g2735   114
# NODE_990_length_149944_cov_81.983368    AUGUSTUS        gene    121307  125898  0.61    +       .       g9290   107

# Of note the gene g13840 is highly expressed. We can look for its interposcan annotation
INTERPRO_DIR=/home/groups/harrisonlab/project_files/fusarium/gene_pred/interproscan/Fus2/split
cat $INTERPRO_DIR/Fus2_interproscan.gff3 | grep 'g13840.' | less

# 1. Results indicate that gene g13840 it is a transmembrane protein involved
# in peptide transport.
# g13840.t1       Pfam    protein_match   181     476     7.6E-27 +       .       Name=PF00854;signature_desc=POT family;Target=g13840.t1 181 476;status=T;ID=match$196_181_476;Ontology_term="GO:0005215","GO:0006810","GO:0016020";date=06-01-2015;Dbxref="InterPro:IPR000109","Reactome:REACT_15518","Reactome:REACT_19419"
# g13840.t1       SUPERFAMILY     protein_match   144     327     5.75E-13        +       .       Name=SSF103473;Target=g13840.t1 144 327;status=T;ID=match$201_144_327;date=06-01-2015;Dbxref="InterPro:IPR020846"
# 1. Visual inspection of g12840 shows all the alignments
# between contig positions 9331 and 9381
# This is probably part of an onion read.

# 2. g11428 appears to be a viral gene
# 2. Visual inspection shows a high level of heterogeneity
# between aligned reads and the FUS2 Genome. Certainly many
# reads present covering the gene and the neighbouring region.

# 3. g2042 is a ubiquitin family protein. Which is foundubiquitously in Eukaryotes.
# - homology to onion?
# Visual inspection of the reads indicate that there is high
# heterogeneity between the aligned reads and the Fus2 Genome.

# 4. g15334 is a viral capsid protein.
# Visual inspection shows that thsis gene is larger than predicted by augustus.
# There is are some notable differences between the aligned reads and Fus2 Genome.
# There is apporximately 10 genomic 'errors' in the fist 983bp of the contig.

for TIMEPOINT in 0 4 8 16 24.1 24.2 36 48 72 96; do
	mkdir -p expreiment1/F.oxysporum_fsp_cepae/Fus2/$TIMEPOINT
done

cp timecourse/expreiment1/1_R1_fus_q30_l50.fastq timecourse/expreiment1/F.oxysporum_fsp_cepae/Fus2/0/.
cp timecourse/expreiment1/4hr-root_R1_fus_q30_l50.fastq timecourse/expreiment1/F.oxysporum_fsp_cepae/Fus2/4/.
cp timecourse/expreiment1/8hr-root_R1_fus_q30_l50.fastq timecourse/expreiment1/F.oxysporum_fsp_cepae/Fus2/8/.
cp timecourse/expreiment1/3_fus_R1_fus_q30_l50.fastq timecourse/expreiment1/F.oxysporum_fsp_cepae/Fus2/16/.
cp timecourse/expreiment1/4_R1_fus_q30_l50.fastq timecourse/expreiment1/F.oxysporum_fsp_cepae/Fus2/24.1/.
cp timecourse/expreiment1/24hr-root_R1_fus_q30_l50.fastq timecourse/expreiment1/F.oxysporum_fsp_cepae/Fus2/24.2/.
cp timecourse/expreiment1/36h-root_R1_fus_q30_l50.fastq timecourse/expreiment1/F.oxysporum_fsp_cepae/Fus2/36/.
cp timecourse/expreiment1/6_R1_fus_q30_l50.fastq timecourse/expreiment1/F.oxysporum_fsp_cepae/Fus2/48/.
cp timecourse/expreiment1/7_R1_fus_q30_l50.fastq timecourse/expreiment1/F.oxysporum_fsp_cepae/Fus2/72/.
cp timecourse/expreiment1/8_R1_fus_q30_l50.fastq timecourse/expreiment1/F.oxysporum_fsp_cepae/Fus2/96/.

cp timecourse/expreiment1/1fus.sam timecourse/expreiment1/F.oxysporum_fsp_cepae/Fus2/0/.
cp timecourse/expreiment1/4hfus.sam timecourse/expreiment1/F.oxysporum_fsp_cepae/Fus2/4/.
cp timecourse/expreiment1/8hfus.sam timecourse/expreiment1/F.oxysporum_fsp_cepae/Fus2/8/.
cp timecourse/expreiment1/3fus.sam timecourse/expreiment1/F.oxysporum_fsp_cepae/Fus2/16/.
cp timecourse/expreiment1/4fus.sam timecourse/expreiment1/F.oxysporum_fsp_cepae/Fus2/24.1/.
cp timecourse/expreiment1/24hfus.sam timecourse/expreiment1/F.oxysporum_fsp_cepae/Fus2/24.2/.
cp timecourse/expreiment1/36hfus.sam timecourse/expreiment1/F.oxysporum_fsp_cepae/Fus2/36/.
cp timecourse/expreiment1/6fus.sam timecourse/expreiment1/F.oxysporum_fsp_cepae/Fus2/48/.
cp timecourse/expreiment1/7fus.sam timecourse/expreiment1/F.oxysporum_fsp_cepae/Fus2/72/.
cp timecourse/expreiment1/8fus.sam timecourse/expreiment1/F.oxysporum_fsp_cepae/Fus2/96/.


# qsub /home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq/gene_expression.sh timecourse/expreiment1/F.oxysporum_fsp_cepae/Fus2/4/4hfus.sam repeat_masked/F.oxysporum_fsp_cepae/Fus2/Fus2_combined_49/Fus2_combined_49_contigs_hardmasked.fa gene_pred/augustus/F.oxysporum_fsp_cepae/Fus2/Fus2_aug_out.gff -->

```bash
  for SamFile in $(ls $WorkDir/Fus2/*/accepted_hits.bam); do
    ScriptDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
    Genome=repeat_masked/F.oxysporum_fsp_cepae/Fus2/Fus2_combined_49/Fus2_combined_49_contigs_hardmasked.fa
    GeneGff=gene_pred/augustus/F.oxysporum_fsp_cepae/Fus2/Fus2_augustus_preds.gtf
    echo $SamFile
    qsub $ScriptDir/gene_expression.sh $SamFile $Genome $GeneGff
  done
```

```bash
  for ExpressedGenes in $(ls timecourse/v2_genes/Fus2/*/*_expressed_genes.bed); do
    OutFile=$(echo $ExpressedGenes | sed 's/_expressed_genes.bed/_expressed_genes_sorted_by_genelength2.bed/g')
    cat $ExpressedGenes | while read Line; do
    Start=$(echo $Line | cut -f4 -d ' ');
    Stop=$(echo $Line | cut -f5 -d ' ');
    Diff=$(( $Stop - $Start ));  
    Reads=$(echo $Line | cut -f10 -d ' ');
    Cov=$((( $Reads * 1000 ) / $Diff ))
    echo -e "$Line\t$Diff\t$Cov";
    done | sort -n -r -k12 > $OutFile
  done
```

```bash
  for FILE in $(ls timecourse/expreiment1/F.oxysporum_fsp_cepae/Fus2/*/*expressed_genes.bed); do
  	printf "$FILE\n" >> tmp3/top_50_genes.csv
  	cat $FILE | head -n50 >> tmp3/top_50_genes.csv
  done
```

  ```bash
  # identify is any of the top 100 expressed genes in Fus2 come from the same contig.
  cat timecourse/expreiment1/F.oxysporum_fsp_cepae/Fus2/96/Fus2_96_expressed_genes.bed | head -n 100 | cut -f1 | sort | uniq -c | sort -r | less
  # 	  4 NODE_9643_length_1511754_cov_84.114418
  #       4 NODE_4583_length_445765_cov_84.597122
  #       4 NODE_1105_length_605632_cov_83.266510
  #       3 NODE_566_length_732593_cov_83.236145
  #       3 NODE_2374_length_16459_cov_83.673859
  #       3 NODE_1097_length_595114_cov_83.826759
  #       2 NODE_956_length_295543_cov_83.945419
  #       2 NODE_731_length_218737_cov_86.751648
  #       2 NODE_7096_length_311910_cov_84.936836
  #       2 NODE_6283_length_399598_cov_84.774551
  #       2 NODE_608_length_519399_cov_83.445702
  #       2 NODE_4974_length_475916_cov_84.779022
  #       2 NODE_4873_length_67562_cov_92.001808
  #       2 NODE_4575_length_487109_cov_84.320221
  #       2 NODE_256_length_160675_cov_85.558762
  #       2 NODE_2122_length_266153_cov_81.544640
  #       2 NODE_1420_length_454186_cov_84.973450
  #       2 NODE_10255_length_363655_cov_83.498451
  #       1 NODE_979_length_171603_cov_83.119377

  cat analysis/blast_homology/F.oxysporum_fsp_cepae/Fus2/Fus2_six-appended_parsed.fa_homologs.csv | cut -f 1,9,14,15 | sort -k2 | grep 'NODE' | cat
  # Fusarium_oxysporum_f._sp._lycopersici_SIX3_gene_for_Secreted_in_xylem_3_protein	NODE_1364_length_4671_cov_108.113892	1100	1657
  # Fusarium_oxysporum_f._sp._lycopersici_isolate_14844_secreted_in_xylem_3_(SIX3)_gene,_complete_cds	NODE_1364_length_4671_cov_108.113892	1136	1624
  # Fusarium_oxysporum_f._sp._lycopersici_isolate_BFOL-51_secreted_in_xylem_3_(SIX3)_gene,_complete_cds	NODE_1364_length_4671_cov_108.113892	1136	1624
  # Fusarium_oxysporum_f._sp._lycopersici_isolate_FOL-MM10_secreted_in_xylem_3_(SIX3)_gene,_complete_cds	NODE_1364_length_4671_cov_108.113892	1136	1624
  # Fusarium_oxysporum_f._sp._lycopersici_isolate_IPO3_secreted_in_xylem_3_(SIX3)_gene,_complete_cds	NODE_1364_length_4671_cov_108.113892	1136	1624
  # Fusarium_oxysporum_f._sp._lycopersici_isolate_BFOL-51_secreted_in_xylem_5_(SIX5)_gene,_partial_cds	NODE_2080_length_2535_cov_50.341618	310	836
  # Fusarium_oxysporum_f._sp._lycopersici_Six5_mRNA,_complete_cds	NODE_2080_length_2535_cov_50.341618	706	843
  # Fusarium_oxysporum_f._sp._lilii_isolate_NRRL_28395_secreted_in_xylem_7-like_protein_(SIX7)_gene,_complete_cds	NODE_2938_length_1682_cov_52.927467	601	1266
  # Fusarium_oxysporum_f._sp._lycopersici_isolate_BFOL-51_secreted_in_xylem_7_(SIX7)_gene,_complete_cds	NODE_2938_length_1682_cov_52.927467	601	1266
  # Fusarium_oxysporum_f._sp._lycopersici_secreted_in_xylem_Six7_(SIX7)_mRNA,_complete_cds	NODE_2938_length_1682_cov_52.927467	601	1266
  # Fusarium_oxysporum_f._sp._lycopersici_strain_Fol007_Six10_(SIX10)_mRNA,_complete_cds	NODE_788_length_4929_cov_54.397850	4328	4697
  # Fusarium_oxysporum_f._sp._lycopersici_strain_Fol007_Six12_(SIX12)_mRNA,_complete_cds	NODE_788_length_4929_cov_54.397850	750	1091
  #
  #In Summary:
  # Fusarium_oxysporum_f._sp._lycopersici_isolate_14844_secreted_in_xylem_3_(SIX3)_gene,_complete_cds
  # NODE_1364_length_4671_cov_108.113892	1136	1624
  # Fusarium_oxysporum_f._sp._lycopersici_isolate_BFOL-51_secreted_in_xylem_5_(SIX5)_gene,_partial_cds
  # NODE_2080_length_2535_cov_50.341618	310	836
  # Fusarium_oxysporum_f._sp._lilii_isolate_NRRL_28395_secreted_in_xylem_7-like_protein_(SIX7)_gene,_complete_cds
  # NODE_2938_length_1682_cov_52.927467	601	1266
  # Fusarium_oxysporum_f._sp._lycopersici_strain_Fol007_Six10_(SIX10)_mRNA,_complete_cds
  # NODE_788_length_4929_cov_54.397850	4328	4697
  # Fusarium_oxysporum_f._sp._lycopersici_strain_Fol007_Six12_(SIX12)_mRNA,_complete_cds	NODE_788_length_4929_cov_54.397850	750	1091

  #Identify expression of Fusarium SIX genes.
  qsub /home/armita/git_repos/emr_repos/tools/pathogen/blast/blast_pipe.sh analysis/blast_homology/six_genes/six-appended_parsed.fa dna repeat_masked/F.oxysporum_fsp_cepae/Fus2/Fus2_combined_49/Fus2_combined_49_contigs_unmasked.fa
  /home/armita/git_repos/emr_repos/tools/pathogen/blast/blast2gff.pl six_gene analysis/blast_homology/F.oxysporum_fsp_cepae/Fus2/Fus2_six-appended_parsed.fa_homologs.csv > analysis/blast_homology/F.oxysporum_fsp_cepae/Fus2/Fus2_six-appended_parsed.fa_homologs.gff
  cat timecourse/expreiment1/F.oxysporum_fsp_cepae/Fus2/96/Fus2_96_expressed_genes.bed | grep 'NODE_4546' | less
  bedtools intersect -c -a analysis/blast_homology/F.oxysporum_fsp_cepae/Fus2/Fus2_six-appended_parsed.fa_homologs.gff -b timecourse/expreiment1/F.oxysporum_fsp_cepae/Fus2/96/Fus2_96_sorted.bam > timecourse/expreiment1/F.oxysporum_fsp_cepae/Fus2/96/Fus2_96_expressed_six.bed
```
