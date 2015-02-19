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


# for strain_path in $(ls -d gene_pred/augustus/F.oxysporum_fsp_cepae/*); do  
# 	echo $strain_path;
# 	STRAIN=$(basename $strain_path); 
# 	echo $STRAIN
# 	qsub /home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/crinkler/sub_crinkler.sh $strain_path/"$STRAIN"_augustus_preds.codingseq RxLR
# done