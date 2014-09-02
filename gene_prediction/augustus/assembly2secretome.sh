SCRIPT_DIR=/home/armita/git_repos/emr_repos/scripts/fusarium/gene_prediction/augustus

qsub $SCRIPT_DIR/submit_augustus.sh fusarium_graminearum assembly/velvet/F.oxysporum_fsp_cepae/125/125_assembly.71/sorted_contigs.fa 
qsub $SCRIPT_DIR/submit_augustus.sh fusarium_graminearum assembly/velvet/F.oxysporum_fsp_cepae/55/55_assembly.61/sorted_contigs.fa 
qsub $SCRIPT_DIR/submit_augustus.sh fusarium_graminearum assembly/velvet/F.oxysporum_fsp_cepae/A23/A23_assembly.61/sorted_contigs.fa 
qsub $SCRIPT_DIR/submit_augustus.sh fusarium_graminearum assembly/velvet/F.oxysporum_fsp_cepae/A28/A28_assembly.61/sorted_contigs.fa 
qsub $SCRIPT_DIR/submit_augustus.sh fusarium_graminearum assembly/velvet/F.oxysporum_fsp_cepae/D2/D2_assembly.41/sorted_contigs.fa 
qsub $SCRIPT_DIR/submit_augustus.sh fusarium_graminearum assembly/velvet/F.oxysporum_fsp_cepae/Fus2/Fus2_assembly.81/sorted_contigs.fa 
qsub $SCRIPT_DIR/submit_augustus.sh fusarium_graminearum assembly/velvet/F.oxysporum_fsp_cepae/HB17/HB17_assembly.51/sorted_contigs.fa 
qsub $SCRIPT_DIR/submit_augustus.sh fusarium_graminearum assembly/velvet/F.oxysporum_fsp_cepae/PG/PG_assembly.81/sorted_contigs.fa 


/home/armita/prog/augustus/scripts/getAnnoFasta.pl gene_pred/augustus/F.oxysporum_fsp_cepae/125/125_aug_out.txt
/home/armita/prog/augustus/scripts/getAnnoFasta.pl gene_pred/augustus/F.oxysporum_fsp_cepae/55/55_aug_out.txt
/home/armita/prog/augustus/scripts/getAnnoFasta.pl gene_pred/augustus/F.oxysporum_fsp_cepae/A23/A23_aug_out.txt
/home/armita/prog/augustus/scripts/getAnnoFasta.pl gene_pred/augustus/F.oxysporum_fsp_cepae/A28/A28_aug_out.txt
/home/armita/prog/augustus/scripts/getAnnoFasta.pl gene_pred/augustus/F.oxysporum_fsp_cepae/D2/D2_aug_out.txt
/home/armita/prog/augustus/scripts/getAnnoFasta.pl gene_pred/augustus/F.oxysporum_fsp_cepae/Fus2/Fus2_aug_out.txt
/home/armita/prog/augustus/scripts/getAnnoFasta.pl gene_pred/augustus/F.oxysporum_fsp_cepae/HB17/HB17_aug_out.txt
/home/armita/prog/augustus/scripts/getAnnoFasta.pl gene_pred/augustus/F.oxysporum_fsp_cepae/PG/PG_aug_out.txt

# signalp-2.0 -t euk -f summary -trunc 70 gene_pred/tmp/gene_preds_1-500.aa > gene_pred/tmp/gene_preds_1-500_sp.txt
# tail -n +5 gene_pred/tmp/gene_preds_1-500_sp.txt > gene_pred/tmp/gene_preds_1-500_sp.txt2
# echo '----------------------------------------------------------------------' >> gene_pred/tmp/gene_preds_1-500_sp.txt2
# /home/armita/git_repos/emr_repos/tools/pathogen/rxlr/annotate_signalP2hmm3_v3.pl gene_pred/tmp/gene_preds_1-500_sp.txt2 gene_pred/tmp/gene_preds_1-500_sp.tab gene_pred/tmp/gene_preds_1-500_sp.aa gene_pred/tmp/gene_preds_1-500_sp_neg.fa gene_pred/tmp/gene_preds_1-500.aa

SIGP_SCRIPTS=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
	
for STRAIN in gene_pred/augustus/F.oxysporum_fsp_cepae/*; do
	STRAIN=$(echo $STRAIN | cut -d '/' -f4)
	$SIGP_SCRIPTS/splitfile_500.pl gene_pred/augustus/F.oxysporum_fsp_cepae/$STRAIN/"$STRAIN"_aug_out.aa
	for FILE in gene_pred/augustus/F.oxysporum_fsp_cepae/$STRAIN/"$STRAIN"_aug_out.aa_split_*.fa; do  
		qsub $SIGP_SCRIPTS/pred_sigP.sh $FILE; 
	done
done
	
for STRAIN in gene_pred/augustus/F.oxysporum_fsp_cepae/*; do
	STRAIN=$(echo $STRAIN | cut -d '/' -f4)
	IN_AA_STRING=''
	IN_NEG_STRING=''
	IN_TAB_STRING=''
	IN_TXT_STRING=''
	for GRP in $(ls -l gene_pred/augustus/F.oxysporum_fsp_cepae/$STRAIN/"$STRAIN"_aug_out.aa_split_*.fa | cut -d '_' -f8 | sort -n); do  
		IN_AA_STRING="$IN_AA_STRING gene_pred/sigP/F.oxysporum_fsp_cepae/$STRAIN/split/$STRAIN""_aug_out_split_$GRP""_sp.aa";  
		IN_NEG_STRING="$IN_NEG_STRING gene_pred/sigP/F.oxysporum_fsp_cepae/$STRAIN/split/$STRAIN""_aug_out_split_$GRP""_sp_neg.aa";  
		IN_TAB_STRING="$IN_TAB_STRING gene_pred/sigP/F.oxysporum_fsp_cepae/$STRAIN/split/$STRAIN""_aug_out_split_$GRP""_sp.tab"; 
		IN_TXT_STRING="$IN_TXT_STRING gene_pred/sigP/F.oxysporum_fsp_cepae/$STRAIN/split/$STRAIN""_aug_out_split_$GRP""_sp.txt";  
	done
cat $IN_AA_STRING > gene_pred/sigP/F.oxysporum_fsp_cepae/$STRAIN/"$STRAIN"_sp.aa
cat $IN_NEG_STRING > gene_pred/sigP/F.oxysporum_fsp_cepae/$STRAIN/"$STRAIN"_neg_sp.aa
tail -n +2 -q $IN_TAB_STRING > gene_pred/sigP/F.oxysporum_fsp_cepae/$STRAIN/"$STRAIN"_sp.tab
cat $IN_TXT_STRING > gene_pred/sigP/F.oxysporum_fsp_cepae/$STRAIN/"$STRAIN"_sp.txt
done

# SIGP_SCRIPTS=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
# $SIGP_SCRIPTS/splitfile_500.pl gene_pred/augustus/F.oxysporum_fsp_cepae/125/125_aug_out.aa
# 
# for FILE in gene_pred/augustus/F.oxysporum_fsp_cepae/125/125_aug_out.aa_split_*.fa; do  
# 	qsub $SIGP_SCRIPTS/pred_sigP.sh $FILE; 
# done
# IN_AA_STRING=''
# IN_NEG_STRING=''
# IN_TAB_STRING=''
# IN_TXT_STRING=''
# 
# for GRP in $(ls -l gene_pred/augustus/F.oxysporum_fsp_cepae/125/125_aug_out.aa_split_*.fa | cut -d '_' -f8 | sort -n); do  
# 	IN_AA_STRING="$IN_AA_STRING gene_pred/sigP/F.oxysporum_fsp_cepae/125/split/125_aug_out_split_$GRP""_sp.aa";  
# 	IN_NEG_STRING="$IN_NEG_STRING gene_pred/sigP/F.oxysporum_fsp_cepae/125/split/125_aug_out_split_$GRP""_sp_neg.aa";  
# 	IN_TAB_STRING="$IN_TAB_STRING gene_pred/sigP/F.oxysporum_fsp_cepae/125/split/125_aug_out_split_$GRP""_sp.tab"; 
# 	IN_TXT_STRING="$IN_TXT_STRING gene_pred/sigP/F.oxysporum_fsp_cepae/125/split/125_aug_out_split_$GRP""_sp.txt";  
# done
# cat $IN_AA_STRING > gene_pred/sigP/F.oxysporum_fsp_cepae/125/125_sp.aa
# cat $IN_NEG_STRING > gene_pred/sigP/F.oxysporum_fsp_cepae/125/125_neg_sp.aa
# tail -n +2 -q $IN_TAB_STRING > gene_pred/sigP/F.oxysporum_fsp_cepae/125/125_sp.tab
# cat $IN_TXT_STRING > gene_pred/sigP/F.oxysporum_fsp_cepae/125/125_sp.txt
# 
# mkdir -p gene_pred/rxlr/F.oxysporum_fsp_cepae/125/
# /home/armita/git_repos/emr_repos/tools/pathogen/rxlr/find_rxlr_v2.pl gene_pred/sigP/F.oxysporum_fsp_cepae/125/125_sp.aa gene_pred/rxlr/F.oxysporum_fsp_cepae/125/125_sp_rxlr.aa gene_pred/rxlr/F.oxysporum_fsp_cepae/125/125_sp_rxlr.txt

for STRAIN in gene_pred/sigP/F.oxysporum_fsp_cepae/*; do
	STRAIN=$(echo $STRAIN | cut -d '/' -f4)
	mkdir -p gene_pred/rxlr/F.oxysporum_fsp_cepae/$STRAIN/
	/home/armita/git_repos/emr_repos/tools/pathogen/rxlr/find_rxlr_v2.pl gene_pred/sigP/F.oxysporum_fsp_cepae/$STRAIN/"$STRAIN"_sp.aa gene_pred/rxlr/F.oxysporum_fsp_cepae/$STRAIN/"$STRAIN"_sp_rxlr.out gene_pred/rxlr/F.oxysporum_fsp_cepae/$STRAIN/"$STRAIN"_sp_rxlr.txt
	tail -n +2 gene_pred/rxlr/F.oxysporum_fsp_cepae/$STRAIN/"$STRAIN"_sp_rxlr.out | head -n -1 > gene_pred/rxlr/F.oxysporum_fsp_cepae/$STRAIN/"$STRAIN"_sp_rxlr.tmp
	/home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/blast_homology/identify_homolgs/append_rxlr.pl gene_pred/rxlr/F.oxysporum_fsp_cepae/$STRAIN/"$STRAIN"_sp_rxlr.tmp > gene_pred/rxlr/F.oxysporum_fsp_cepae/$STRAIN/"$STRAIN"_sp_rxlr.fa
	rm gene_pred/rxlr/F.oxysporum_fsp_cepae/$STRAIN/"$STRAIN"_sp_rxlr.tmp
done

cat gene_pred/rxlr/F.oxysporum_fsp_cepae/*/*_sp_rxlr.fa > analysis/blast_homology/sp_rxlr/F.oxysporum_fsp_cepae_sp_rxlr.fa

qsub /home/armita/git_repos/emr_repos/tools/pathogen/blast/blast_pipe.sh analysis/blast_homology/sp_rxlr/F.oxysporum_fsp_cepae_sp_rxlr.fa protein assembly/velvet/F.oxysporum_fsp_cepae/125/125_assembly.71/sorted_contigs.fa
qsub /home/armita/git_repos/emr_repos/tools/pathogen/blast/blast_pipe.sh analysis/blast_homology/sp_rxlr/F.oxysporum_fsp_cepae_sp_rxlr.fa protein assembly/velvet/F.oxysporum_fsp_cepae/55/55_assembly.61/sorted_contigs.fa
qsub /home/armita/git_repos/emr_repos/tools/pathogen/blast/blast_pipe.sh analysis/blast_homology/sp_rxlr/F.oxysporum_fsp_cepae_sp_rxlr.fa protein assembly/velvet/F.oxysporum_fsp_cepae/A23/A23_assembly.61/sorted_contigs.fa
qsub /home/armita/git_repos/emr_repos/tools/pathogen/blast/blast_pipe.sh analysis/blast_homology/sp_rxlr/F.oxysporum_fsp_cepae_sp_rxlr.fa protein assembly/velvet/F.oxysporum_fsp_cepae/A28/A28_assembly.61/sorted_contigs.fa
qsub /home/armita/git_repos/emr_repos/tools/pathogen/blast/blast_pipe.sh analysis/blast_homology/sp_rxlr/F.oxysporum_fsp_cepae_sp_rxlr.fa protein assembly/velvet/F.oxysporum_fsp_cepae/D2/D2_assembly.41/sorted_contigs.fa
qsub /home/armita/git_repos/emr_repos/tools/pathogen/blast/blast_pipe.sh analysis/blast_homology/sp_rxlr/F.oxysporum_fsp_cepae_sp_rxlr.fa protein assembly/velvet/F.oxysporum_fsp_cepae/Fus2/Fus2_assembly.81/sorted_contigs.fa
qsub /home/armita/git_repos/emr_repos/tools/pathogen/blast/blast_pipe.sh analysis/blast_homology/sp_rxlr/F.oxysporum_fsp_cepae_sp_rxlr.fa protein assembly/velvet/F.oxysporum_fsp_cepae/HB17/HB17_assembly.51/sorted_contigs.fa
qsub /home/armita/git_repos/emr_repos/tools/pathogen/blast/blast_pipe.sh analysis/blast_homology/sp_rxlr/F.oxysporum_fsp_cepae_sp_rxlr.fa protein assembly/velvet/F.oxysporum_fsp_cepae/PG/PG_assembly.81/sorted_contigs.fa

/home/armita/git_repos/emr_repos/tools/pathogen/blast/blast_differentials.pl analysis/blast_homology/F.oxysporum_fsp_cepae/*/*_F.oxysporum_fsp_cepae_sp_rxlr.fa_homologs.csv
mv *.csv analysis/blast_homology/sp_rxlr/.
/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/blast_homology/identify_differentials/homolog_grp_srt.pl analysis/blast_homology/F.oxysporum_fsp_cepae/125/125_F.oxysporum_fsp_cepae_sp_rxlr.fa_homologs.csv analysis/blast_homology/F.oxysporum_fsp_cepae/55/55_F.oxysporum_fsp_cepae_sp_rxlr.fa_homologs.csv analysis/blast_homology/F.oxysporum_fsp_cepae/A23/A23_F.oxysporum_fsp_cepae_sp_rxlr.fa_homologs.csv analysis/blast_homology/F.oxysporum_fsp_cepae/A28/A28_F.oxysporum_fsp_cepae_sp_rxlr.fa_homologs.csv analysis/blast_homology/F.oxysporum_fsp_cepae/D2/D2_F.oxysporum_fsp_cepae_sp_rxlr.fa_homologs.csv analysis/blast_homology/F.oxysporum_fsp_cepae/Fus2/Fus2_F.oxysporum_fsp_cepae_sp_rxlr.fa_homologs.csv analysis/blast_homology/F.oxysporum_fsp_cepae/HB17/HB17_F.oxysporum_fsp_cepae_sp_rxlr.fa_homologs.csv analysis/blast_homology/F.oxysporum_fsp_cepae/PG/PG_F.oxysporum_fsp_cepae_sp_rxlr.fa_homologs.csv
mv *.csv analysis/blast_homology/sp_rxlr/.

grep '>' analysis/rxlr/P.cactorum/404/404_sp_rxlr.fa  | cut -f1 > tmp.fa
	while read line; do 
	grep -A1 -w "$line" analysis/rxlr/P.cactorum/404/404_nuc.fa >> analysis/rxlr/P.cactorum/404/404_sp_rxlr_nuc.fa; 
done<tmp.fa
grep -A1 '>' analysis/rxlr/P.cactorum/404/404_sp_rxlr_nuc.fa | cut -d '|' -f5 > tmp.fa
sed 's/\r/_/g' tmp.fa > tmp2.fa
sed 's/>NODE/>P.cactorum_404_NODE/g' tmp2.fa > analysis/rxlr/P.cactorum/404/404_sp_rxlr_nuc_parsed.fa