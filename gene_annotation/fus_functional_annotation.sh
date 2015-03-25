#!/bin/bash
# Functional annotation of Fusarium genomes

# Interproscan

SCRIPT_DIR=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan

for strain_path in $(ls -d gene_pred/augustus/F.oxysporum_fsp_cepae/*); do 
	echo $strain_path
	STRAIN=$(basename $strain_path)
	echo $STRAIN
	$SCRIPT_DIR/sub_interproscan.sh $strain_path/"$STRAIN"_augustus_preds.aa
done

# SigP
SIGP_SCRIPTS=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
for FILEZ in $(ls -d gene_pred/interproscan/F.oxysporum_fsp_cepae/*/*_augustus_preds.aa_split_*); do   
echo $FILEZ;  
qsub $SIGP_SCRIPTS/pred_sigP.sh $FILEZ; 
done

for PATHZ in $(ls -d gene_pred/interproscan/F.oxysporum_fsp_cepae/*); do
STRAIN=$(echo $PATHZ | cut -d '/' -f4)
SPECIES=$(echo $PATHZ | cut -d '/' -f3)
IN_AA_STRING=''
IN_NEG_STRING=''
IN_TAB_STRING=''
IN_TXT_STRING=''
for GRP in $(ls -l gene_pred/interproscan/$SPECIES/$STRAIN/"$STRAIN"_augustus_preds* | cut -d '_' -f8 | sort -n); do  
IN_AA_STRING="$IN_AA_STRING gene_pred/sigP/$SPECIES/$STRAIN/split/$STRAIN""_augustus_preds_split_$GRP""_sp.aa";  
IN_NEG_STRING="$IN_NEG_STRING gene_pred/sigP/$SPECIES/$STRAIN/split/$STRAIN""_augustus_preds_split_$GRP""_sp_neg.aa";  
IN_TAB_STRING="$IN_TAB_STRING gene_pred/sigP/$SPECIES/$STRAIN/split/$STRAIN""_augustus_preds_split_$GRP""_sp.tab"; 
IN_TXT_STRING="$IN_TXT_STRING gene_pred/sigP/$SPECIES/$STRAIN/split/$STRAIN""_augustus_preds_split_$GRP""_sp.txt";  
done
cat $IN_AA_STRING > gene_pred/sigP/$SPECIES/$STRAIN/"$STRAIN"_sp.aa
cat $IN_NEG_STRING > gene_pred/sigP/$SPECIES/$STRAIN/"$STRAIN"_neg_sp.aa
tail -n +2 -q $IN_TAB_STRING > gene_pred/sigP/$SPECIES/$STRAIN/"$STRAIN"_sp.tab
cat $IN_TXT_STRING > gene_pred/sigP/$SPECIES/$STRAIN/"$STRAIN"_sp.txt
done

# RxLR

# Identify sequences containing the RxLR motif in proteins containing signal peptides

for PATHZ in $(ls -d gene_pred/sigP/F.oxysporum_fsp_cepae/*); do 
echo $PATHZ; 
STRAIN=$(printf "$PATHZ" | cut -d '/' -f4); 
SPECIES=$(printf "$PATHZ" | cut -d '/' -f3); 
echo "strain: $STRAIN    species: $SPECIES";
mkdir -p analysis/sigP_rxlr/"$SPECIES"/"$STRAIN"
echo "the number of SigP gene is:"
cat gene_pred/sigP/F.oxysporum_fsp_cepae/"$STRAIN"/"$STRAIN"_sp.aa | grep '>' | wc -l
echo "the number of SigP-RxLR genes are:"
cat gene_pred/sigP/F.oxysporum_fsp_cepae/"$STRAIN"/"$STRAIN"_sp.aa | grep 'R.LR' | wc -l
cat gene_pred/sigP/F.oxysporum_fsp_cepae/"$STRAIN"/"$STRAIN"_sp.aa | grep -B1 'R.LR' | grep -vw '\-\-' > analysis/sigP_rxlr/"$SPECIES"/"$STRAIN"/"$STRAIN"_sigP_RxLR.fa
done

# Extract gene features from gff files for genes in the "$STRAIN"_sigP_RxLR.fa file.
for THIS_DIR in $(ls -d analysis/sigP_rxlr/F.oxysporum_fsp_cepae/*); do 
STRAIN=$(echo $THIS_DIR | rev | cut -f1 -d '/' | rev);
echo $STRAIN; 
GENE_MODELS=$(ls gene_pred/augustus/F.oxysporum_fsp_cepae/$STRAIN/"$STRAIN"_augustus_preds.gtf)
cat $THIS_DIR/"$STRAIN"_sigP_RxLR.fa  | cut -f1 | grep '>' | sed 's/>//g' | sed 's/\.t.*//g' > "$STRAIN"_id_tmp.txt
printf "" > $THIS_DIR/"$STRAIN"_sigP_RxLR.gtf
cat $GENE_MODELS | grep 'gene' | grep -w -f "$STRAIN"_id_tmp.txt > $THIS_DIR/"$STRAIN"_sigP_RxLR.gtf
rm "$STRAIN"_id_tmp.txt
done

cp -r gene_pred/sigP gene_pred/sigP_hints



# The analysis was repeated on gene predictions not using RNAseq data as hints
# This provides an idicaiton of the use of infection RNAseq data to predict 
# pathogenicity gene models. 

# SigP
SIGP_SCRIPTS=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
for FILEZ in $(ls -d gene_pred/augustus/F.oxysporum_fsp_cepae_not_as_old/*/*_aug_out.aa_split_*); do   
echo $FILEZ;  
qsub $SIGP_SCRIPTS/pred_sigP.sh $FILEZ; 
done

for PATHZ in $(ls -d gene_pred/augustus/F.oxysporum_fsp_cepae/*); do
STRAIN=$(echo $PATHZ | cut -d '/' -f4)
SPECIES=$(echo $PATHZ | cut -d '/' -f3)
IN_AA_STRING=''
IN_NEG_STRING=''
IN_TAB_STRING=''
IN_TXT_STRING=''
for GRP in $(ls -l gene_pred/sigP/F.oxysporum_fsp_cepae/$STRAIN/split/*_aug_out_split_* | cut -d '_' -f8 | cut -d '.' -f1 | sort -n | uniq); do  
IN_AA_STRING="$IN_AA_STRING gene_pred/sigP/$SPECIES/$STRAIN/split/$STRAIN""_aug_out_split_$GRP"".fa_sp.aa";  
IN_NEG_STRING="$IN_NEG_STRING gene_pred/sigP/$SPECIES/$STRAIN/split/$STRAIN""_aug_out_split_$GRP"".fa_sp_neg.aa";  
IN_TAB_STRING="$IN_TAB_STRING gene_pred/sigP/$SPECIES/$STRAIN/split/$STRAIN""_aug_out_split_$GRP"".fa_sp.tab"; 
IN_TXT_STRING="$IN_TXT_STRING gene_pred/sigP/$SPECIES/$STRAIN/split/$STRAIN""_aug_out_split_$GRP"".fa_sp.txt";  
done
cat $IN_AA_STRING > gene_pred/sigP/$SPECIES/$STRAIN/"$STRAIN"_sp.aa
cat $IN_NEG_STRING > gene_pred/sigP/$SPECIES/$STRAIN/"$STRAIN"_neg_sp.aa
tail -n +2 -q $IN_TAB_STRING > gene_pred/sigP/$SPECIES/$STRAIN/"$STRAIN"_sp.tab
cat $IN_TXT_STRING > gene_pred/sigP/$SPECIES/$STRAIN/"$STRAIN"_sp.txt
done

# RxLR

# Identify sequences containing the RxLR motif in proteins containing signal peptides

for PATHZ in $(ls -d gene_pred/sigP/F.oxysporum_fsp_cepae/*); do 
echo $PATHZ; 
STRAIN=$(printf "$PATHZ" | cut -d '/' -f4); 
SPECIES=$(printf "$PATHZ" | cut -d '/' -f3); 
echo "strain: $STRAIN    species: $SPECIES";
mkdir -p analysis/sigP_rxlr/"$SPECIES"/"$STRAIN"
echo "the number of SigP gene is:"
cat gene_pred/sigP/F.oxysporum_fsp_cepae/"$STRAIN"/"$STRAIN"_sp.aa | grep '>' | wc -l
echo "the number of SigP-RxLR genes are:"
cat gene_pred/sigP/F.oxysporum_fsp_cepae/"$STRAIN"/"$STRAIN"_sp.aa | grep 'R.LR' | wc -l
cat gene_pred/sigP/F.oxysporum_fsp_cepae/"$STRAIN"/"$STRAIN"_sp.aa | grep -B1 'R.LR' | grep -vw '\-\-' > analysis/sigP_rxlr/"$SPECIES"/"$STRAIN"/"$STRAIN"_sigP_RxLR.fa
done

mv gene_pred/sigP gene_pred/sigP_no_hints


# Position of Mimps
# gff predictions for the position of mimps in the genome were identified
# Fus2 was performed separately to ensure mimps are predicted from the correct assembly
# Unhash to run!
# for PATHZ in $(ls repeat_masked/F.oxysporum_fsp_cepae/Fus2/Fus2_combined_49/Fus2_combined_49_contigs_unmasked.fa); do
for PATHZ in $(ls repeat_masked/F.oxysporum_fsp_cepae/*/*/*_contigs_unmasked.fa); do
	STRAIN=$(echo "$PATHZ" | rev | cut -d '/' -f3 | rev)
	STRAIN=$(echo "$PATHZ" | rev | cut -d '/' -f3 | rev)
	OUTDIR=analysis/mimps/F.oxysporum_fsp_cepae/"$STRAIN"
	mkdir -p "$OUTDIR"
	SCRIPT_DIR="/home/armita/git_repos/emr_repos/tools/pathogen/mimp_finder"
	"$SCRIPT_DIR"/mimp_finder.pl "$PATHZ" "$OUTDIR"/"$STRAIN"_mimps.fa "$OUTDIR"/"$STRAIN"_mimps.gff3 > "$OUTDIR"/"$STRAIN"_mimps.log
done


# Presence of SIX genes
for INFILE in $(ls repeat_masked/F.oxysporum_fsp_cepae/*/*/*_contigs_unmasked.fa); do
SCRIPT_DIR="/home/armita/git_repos/emr_repos/tools/pathogen/blast"
QUERY_FILE="analysis/blast_homology/Fo_path_genes/Fo_path_genes.fa"
echo "$INFILE"
qsub "$SCRIPT_DIR"/blast_pipe.sh $QUERY_FILE dna $INFILE
done
# repeat_masked/F.oxysporum_fsp_cepae/Fus2/Fus2_assembly.81_repmask/Fus2_contigs_unmasked.fa
# Your job 6265521 ("blast_pipe.sh") has been submitted
qdel 6265521

SCRIPT_DIR="/home/armita/git_repos/emr_repos/tools/pathogen/blast"
for INFILE in $(ls analysis/blast_homology/F.oxysporum_fsp_cepae/*/*_Fo_path_genes.fa_homologs.csv); do
OUTFILE=$(echo $INFILE | sed 's/.csv/.gff/g')
echo "$OUTFILE" 
$SCRIPT_DIR/blast2gff.pl Fo_Path_gene 1 "$INFILE" > $OUTFILE
done




# Genes flanking Mimps
ProgPath=/home/armita/git_repos/emr_repos/tools/pathogen/mimp_finder
for Pathz in $(ls -d analysis/mimps/F.oxysporum_fsp_cepae/*); do
Strain=$(echo "$Pathz" | cut -f4 -d '/')
echo "$Strain"
"$ProgPath"/gffexpander.pl + 1000 "$Pathz"/"$Strain"_mimps.gff3 > "$Pathz"/"$Strain"_mimps_exp.gff3
StrainModels=$(ls gene_pred/augustus/F.oxysporum_fsp_cepae/"$Strain"/*_augustus_preds.gtf) 
bedtools intersect -s -a "$StrainModels"  -b "$Pathz"/"$Strain"_mimps_exp.gff3> "$Pathz"/"$Strain"_mimps_ass_genes.bed
cat "$Pathz"/*_mimps_ass_genes.bed | grep 'gene' | wc -l
done
