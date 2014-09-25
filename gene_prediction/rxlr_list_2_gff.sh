#!/bin/bash

# Commands used to create feature tracks for genes containing SigP/RxLRs in fusarium genomes.

SOURCE=RxLR_pipe
FEATURE=SigP_RxLR_gene
SCRIPT_DIR=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation

$SCRIPT_DIR/gene_list_to_gff.pl gene_pred/rxlr/F.oxysporum_fsp_cepae/125/125_sp_rxlr_id.txt gene_pred/augustus/F.oxysporum_fsp_cepae/125/125_aug.gff $SOURCE $FEATURE > gene_pred/rxlr/F.oxysporum_fsp_cepae/125/125_sp_rxlr.gff
$SCRIPT_DIR/gene_list_to_gff.pl gene_pred/rxlr/F.oxysporum_fsp_cepae/55/55_sp_rxlr_id.txt gene_pred/augustus/F.oxysporum_fsp_cepae/55/55_aug.gff $SOURCE $FEATURE > gene_pred/rxlr/F.oxysporum_fsp_cepae/55/55_sp_rxlr.gff
$SCRIPT_DIR/gene_list_to_gff.pl gene_pred/rxlr/F.oxysporum_fsp_cepae/A23/A23_sp_rxlr_id.txt gene_pred/augustus/F.oxysporum_fsp_cepae/A23/A23_aug.gff $SOURCE $FEATURE > gene_pred/rxlr/F.oxysporum_fsp_cepae/A23/A23_sp_rxlr.gff
$SCRIPT_DIR/gene_list_to_gff.pl gene_pred/rxlr/F.oxysporum_fsp_cepae/A28/A28_sp_rxlr_id.txt gene_pred/augustus/F.oxysporum_fsp_cepae/A28/A28_aug.gff $SOURCE $FEATURE > gene_pred/rxlr/F.oxysporum_fsp_cepae/A28/A28_sp_rxlr.gff
$SCRIPT_DIR/gene_list_to_gff.pl gene_pred/rxlr/F.oxysporum_fsp_cepae/D2/D2_sp_rxlr_id.txt gene_pred/augustus/F.oxysporum_fsp_cepae/D2/D2_aug.gff $SOURCE $FEATURE > gene_pred/rxlr/F.oxysporum_fsp_cepae/D2/D2_sp_rxlr.gff
$SCRIPT_DIR/gene_list_to_gff.pl gene_pred/rxlr/F.oxysporum_fsp_cepae/Fus2/Fus2_sp_rxlr_id.txt gene_pred/augustus/F.oxysporum_fsp_cepae/Fus2/Fus2_aug.gff $SOURCE $FEATURE > gene_pred/rxlr/F.oxysporum_fsp_cepae/Fus2/Fus2_sp_rxlr.gff
$SCRIPT_DIR/gene_list_to_gff.pl gene_pred/rxlr/F.oxysporum_fsp_cepae/HB17/HB17_sp_rxlr_id.txt gene_pred/augustus/F.oxysporum_fsp_cepae/HB17/HB17_aug.gff $SOURCE $FEATURE > gene_pred/rxlr/F.oxysporum_fsp_cepae/HB17/HB17_sp_rxlr.gff
$SCRIPT_DIR/gene_list_to_gff.pl gene_pred/rxlr/F.oxysporum_fsp_cepae/PG/PG_sp_rxlr_id.txt gene_pred/augustus/F.oxysporum_fsp_cepae/PG/PG_aug.gff $SOURCE $FEATURE > gene_pred/rxlr/F.oxysporum_fsp_cepae/PG/PG_sp_rxlr.gff
