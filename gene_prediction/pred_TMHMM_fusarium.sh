#!/usr/bin/sh


for INFILE in $(ls -d gene_pred/augustus/F.oxysporum_fsp_cepae/*/*_aug_out.aa); do 
	echo $INFILE
	qsub /home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/transmembrane_helices/submit_TMHMM.sh $INFILE
done

