#!/bin/bash

# Commands used to perform comparative genomics in Fusarium

ls -d gene_pred/augustus/F.oxysporum_fsp_cepae/*
# gene_pred/augustus/F.oxysporum_fsp_cepae/125  gene_pred/augustus/F.oxysporum_fsp_cepae/D2
# gene_pred/augustus/F.oxysporum_fsp_cepae/55   gene_pred/augustus/F.oxysporum_fsp_cepae/Fus2
# gene_pred/augustus/F.oxysporum_fsp_cepae/A23  gene_pred/augustus/F.oxysporum_fsp_cepae/HB17
# gene_pred/augustus/F.oxysporum_fsp_cepae/A28  gene_pred/augustus/F.oxysporum_fsp_cepae/PG

set -- 125 55 A23 A28 D2 Fus2 HB17 PG
for a; do 
	shift
	for b; do
	printf "%s - %s\n" "$a" "$b"
	qsub /home/armita/git_repos/emr_repos/tools/pathogen/orthology/inparanoid/sub_inparanoid.sh gene_pred/augustus/*/$a/"$a"_augustus_preds.aa gene_pred/augustus/*/$b/"$b"_augustus_preds.aa gene_pred/augustus/*/$a/"$a"_augustus_preds.gtf gene_pred/augustus/*/$b/"$b"_augustus_preds.gtf
	done 
done

# Commands not Used yet:
# 
# mkdir analysis/inparanoid/summary_tables
# cat analysis/inparanoid/*/*_seqs.txt | cut -f1 | sort | uniq | less
# cat analysis/inparanoid/*/*_seqs.txt | cut -f1 | sort | uniq > analysis/inparanoid/summary_tables/all_genes.txt
# 
# /home/armita/git_repos/emr_repos/tools/pathogen/orthology/inparanoid/build_orthology_tab.pl analysis/inparanoid/summary_tables/all_genes.txt analysis/inparanoid/*/sqltable.* > analysis/inparanoid/summary_tables/rxlr_orthology_tab.csv
# 
# cat analysis/inparanoid/10300-404/*_seqs.txt analysis/inparanoid/10300-414/*_seqs.txt analysis/inparanoid/404-414/*_seqs.txt | cut -f1 | sort | uniq > analysis/inparanoid/summary_tables/all_P.cac_RxLR.txt
# /home/armita/git_repos/emr_repos/tools/pathogen/orthology/inparanoid/build_orthology_tab.pl analysis/inparanoid/summary_tables/all_P.cac_RxLR.txt analysis/inparanoid/10300-404/sqltable.10300-404 analysis/inparanoid/10300-414/sqltable.10300-414 analysis/inparanoid/404-414/sqltable.404-414 > analysis/inparanoid/summary_tables/rxlr_orthology_P.cact.csv
# 
# 
# for FILE in $(ls analysis/inparanoid/summary_tables/orthogroups4/*.txt); do 
# 	cat $FILE | cut -f1 -d '|' | uniq -dc; printf "$FILE\n"
# done | grep -B1 '414' | less
# cat analysis/inparanoid/summary_tables/new_final_tab.csv analysis/rxlr_unmasked/P.cactorum/10300/10300_sp_rxlr.fa | grep -o '10300.*' | cut -f1 | sed s'/\.txt//' | sed s'/_singleton_insl.*bp//' | sed 's/_Pcac10300//' | sed 's/_contig//' | sed 's/ //' | tail -n+2 | sort | uniq -d | wc -l
# cat analysis/inparanoid/summary_tables/new_final_tab.csv analysis/rxlr_unmasked/P.cactorum/10300/10300_sp_rxlr.fa | grep -o '10300.*' | cut -f1 | sed s'/\.txt//' | sed s'/_singleton_insl.*bp//' | sed 's/_Pcac10300//' | sed 's/_contig//' | sed 's/ //' | tail -n+2 | sort | uniq -u | wc -l