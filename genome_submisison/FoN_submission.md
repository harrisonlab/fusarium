# Submission Commands

Submisison of annotations with an assembly appears to be a complex process.
If a genome is to be submitted without annotation then all that is needed is the
fasta file containing the assembled contigs. If an annotated genome is to be
submitted then a number of processing steps are required before submission. The
fasta file of contigs and the gff file of annotations must be combined to form a
.asn file. The program that does this conversion (tbl2asn) requires the fasta
files and gff files to be formatted correctly. In the case of the gff file, this
means parsing it to a .tbl file.

The commands used to parse these files and prepare the F. oxysporum f. sp.
narcissi genome for submisson are shown below.



# Preliminary submission

A Bioproject and biosample number was prepared for the genome submission at:
https://submit.ncbi.nlm.nih.gov

A preliminary submission was made for the .fasta assembly to check if
any contigs needed to be split. This step was performed early in the annotation
process (prior to gene prediction) to ensure that annotation did not have to
be repeated at the end of the project.


The following note was provided in the WGS submission page on NCBI in the box
labeled "Private comments to NCBI staff":

```
I have been advised to submit my assemblies to NCBI early in my submission process to ensure that my contigs pass the contamination screen. This assembly will be revised as appropriate, including renaming of contigs where needed. Please allow me to modify this submission at a later date, including upload of the final gene models.

'For future submissions, you could send us the fasta files early
in the submission process so we can run them through our foreign
contamination screen. We will let you know if we find any
sequences to exclude or trim before you generate your final
WGS submission.'...'*IMPORTANT* Include a comment that you are submitting
the fasta files to be screened by the contamination screen
prior to creating your final annotated submission.'
```


# Submission of sequence data to SRA

Reads were submitted to the SRA at https://submit.ncbi.nlm.nih.gov/subs/sra/ .
To do this, a metadata file was provided detailing each of the files in the
bioproject. The file was downloaded in excel format and edited manually. A copy
of the edited file and the final .tsv file is present at:

```bash
  ls genome_submission/SRA_metadata_acc.txt genome_submission/SRA_metadata_acc.xlsx
```

As these files included a file > 500 Mb, a presubmission folder was requested.
This aids submission of large data files. This file was created on the ftp server
at ftp-private.ncbi.nlm.nih.gov, with a private folder named
uploads/andrew.armitage@emr.ac.uk_6L2oakBI. Ncbi provided a username a password.
Files were uploaded into a folder created within my preload folder using ftp.

```bash
ReadsF=$(ls raw_dna/paired/F.oxysporum_fsp_narcissi/N139/F/FoxysporumN139_S2_L001_R1_001.fastq.gz)
ReadsR=$(ls raw_dna/paired/F.oxysporum_fsp_narcissi/N139/R/FoxysporumN139_S2_L001_R2_001.fastq.gz)
# Bioproject="PRJNA338236"
SubFolder="FoN_PRJNA338236"
mkdir $SubFolder
cp $ReadsF $SubFolder/.
cp $ReadsR $SubFolder/.
ftp ftp-private.ncbi.nlm.nih.gov
cd uploads/andrew.armitage@emr.ac.uk_6L2oakBI
mkdir FoN_PRJNA338236
cd FoN_PRJNA338236
# put FoN_PRJNA338236
prompt
mput FoxysporumN139_S2_L001_R1_001.fastq.gz FoxysporumN139_S2_L001_R2_001.fastq.gz
put raw_dna/paired/F.oxysporum_fsp_narcissi/N139/F/FoxysporumN139_S2_L001_R1_001.fastq.gz
put raw_dna/paired/F.oxysporum_fsp_narcissi/N139/R/FoxysporumN139_S2_L001_R2_001.fastq.gz

```
## SbtFile
The genbank submission template tool was used at:
http://www.ncbi.nlm.nih.gov/WebSub/template.cgi
This produce a template file detailing the submission.

# Final Submission

These commands were used in the final submission of the FoN genome:


```bash
for Assembly in $(ls repeat_masked/*/*/*/*_contigs_unmasked.fa | grep -v 'Fus2' | grep 'ncbi' | grep -e 'N139_ncbi' | grep -v 'old'); do
# tbl2asn options:
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
OrganismOfficial=$(echo $Organism | sed 's/F./Fusarium /g' | sed 's/_fsp_/ f.sp. /g')
StrainOfficial=$(echo $Strain | sed 's/_ncbi//g')
#
ProjDir=/home/groups/harrisonlab/project_files/fusarium
cd $ProjDir
OutDir="genome_submission/$Organism/$Strain"
mkdir -p $OutDir

# Program locations:
AnnieDir="/home/armita/prog/annie/genomeannotation-annie-c1e848b"
ProgDir="/home/armita/git_repos/emr_repos/tools/genbank_submission"
# File locations:
# Assembly=$(ls repeat_masked/$Organism/$Strain/edited_contigs_repmask/Fus2_canu_contigs_unmasked.fa)
InterProTab=$(ls gene_pred/interproscan/$Organism/$Strain/"$Strain"_interproscan.tsv)
SwissProtBlast=$(ls gene_pred/swissprot/$Organism/$Strain/swissprot_vJul2016_tophit_parsed.tbl)
SwissProtFasta=$(ls /home/groups/harrisonlab/uniprot/swissprot/uniprot_sprot.fasta)
GffFile=$(ls gene_pred/final_genes/$Organism/$Strain/final/final_genes_appended.gff3)
SbtFile=genome_submission/F.oxysporum_fsp_narcissi/N139/template.sbt

SRA_metadata=$(ls genome_submission/FoC_PRJNA338256_SRA_metadata_acc.txt)
BioProject=$(cat $SRA_metadata | sed 's/PRJNA/\nPRJNA/g' | grep "$StrainOfficial" | cut -f1 | head -n1)
BioSample=$(cat $SRA_metadata | sed 's/PRJNA/\nPRJNA/g' | grep "$StrainOfficial" | cut -f2 | head -n1)


# ncbi_tbl_corrector script options:
SubmissionID="BFJ63"
LabID="ArmitageEMR"
# Final submisison file name:
FinalName="$Organism"_"$Strain"_Armitage_2016

python3 $AnnieDir/annie.py -ipr $InterProTab -g $GffFile -b $SwissProtBlast -db $SwissProtFasta -o $OutDir/annie_output.csv --fix_bad_products
$ProgDir/edit_tbl_file/annie_corrector.py --inp_csv $OutDir/annie_output.csv --out_csv $OutDir/annie_corrected_output.csv

mkdir -p $OutDir/gag/round1
gag.py -f $Assembly -g $GffFile -a $OutDir/annie_corrected_output.csv --fix_start_stop -o $OutDir/gag/round1 2>&1 | tee $OutDir/gag_log1.txt
sed -i 's/Dbxref/db_xref/g' $OutDir/gag/round1/genome.tbl

cp $Assembly $OutDir/gag/round1/genome.fsa  
cp $SbtFile $OutDir/gag/round1/genome.sbt
mkdir -p $OutDir/tbl2asn/round1
tbl2asn -p $OutDir/gag/round1/. -t $OutDir/gag/round1/genome.sbt -r $OutDir/tbl2asn/round1 -M n -X E -Z $OutDir/gag/round1/discrep.txt -j "[organism=$OrganismOfficial] [strain=$StrainOfficial]"

mkdir -p $OutDir/gag/edited
$ProgDir/edit_tbl_file/ncbi_tbl_corrector.py --inp_tbl $OutDir/gag/round1/genome.tbl --inp_val $OutDir/tbl2asn/round1/genome.val --locus_tag $SubmissionID --lab_id $LabID --gene_id "remove" --edits stop pseudo unknown_UTR correct_partial --rename_genes "vAg" --remove_product_locus_tags "True" --out_tbl $OutDir/gag/edited/genome.tbl
printf "StructuredCommentPrefix\t##Genome-Annotation-Data-START##
Annotation Provider\tHarrison Lab NIAB-EMR
Annotation Date\tSEP-2016
Annotation Version\tRelease 1.01
Annotation Method\tAb initio gene prediction: Braker 1.9 and CodingQuary 2.0; Functional annotation: Swissprot (July 2016 release) and Interproscan 5.18-57.0" \
> $OutDir/gag/edited/annotation_methods.strcmt.txt

cp $Assembly $OutDir/gag/edited/genome.fsa
cp $SbtFile $OutDir/gag/edited/genome.sbt
mkdir $OutDir/tbl2asn/final
tbl2asn -p $OutDir/gag/edited/. -t $OutDir/gag/edited/genome.sbt -r $OutDir/tbl2asn/final -M n -X E -Z $OutDir/tbl2asn/final/discrep.txt -j "[organism=$OrganismOfficial] [strain=$StrainOfficial]" -l paired-ends -a r10k -w $OutDir/gag/edited/annotation_methods.strcmt.txt
cat $OutDir/tbl2asn/final/genome.sqn | sed 's/_pilon//g' | sed 's/\. subunit/kDa subunit/g' | sed 's/, mitochondrial//g' > $OutDir/tbl2asn/final/$FinalName.sqn
done
```

```bash
for File in $(ls genome_submission/F.*/*_ncbi/tbl2asn/final/errorsummary.val | grep 'N139'); do
Organism=$(echo $File | rev | cut -f5 -d '/' | rev);
Strain=$(echo $File | rev | cut -f4 -d '/' | rev);
echo "$Organism - $Strain";
cat $File;
echo "Duplicated genes:"
cat genome_submission/$Organism/$Strain/tbl2asn/round1/genome.val | grep 'DuplicateFeat' | cut -f4 -d ':' | cut -f2 -d' '
echo "";
done > genome_submission/FoN_isolate_errors.txt
```


The final error report contained the following warnings. These were judged to be
legitimate concerns but biologically explainable.

67 WARNING: SEQ_FEAT.PartialProblem
 5 WARNING: SEQ_FEAT.ProteinNameEndsInBracket
211 WARNING: SEQ_FEAT.ShortExon
18 WARNING: SEQ_FEAT.SuspiciousFrame
 5 INFO:    SEQ_FEAT.PartialProblem

 Note -
 *SEQ_FEAT.partial problem. In this case, upon investigation these genes were hannging
 off the end of a contig but did not have an mRNA feature that went off of the
 end of the contig. This was occuring due to an intron being predicted hanging
 off the contig. An example on the ncbi guidelines here shows this to be
 acceptable:
 http://www.ncbi.nlm.nih.gov/genbank/eukaryotic_genome_submission_annotation#Partialcodingregionsinincompletegenomes
 *SEQ_FEAT.ProteinNameEndsInBracket. These gene names include brackets for good
 reason
