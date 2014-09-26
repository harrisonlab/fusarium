#!/usr/bin/bash

qsub /home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking/rep_modeling.sh assembly/velvet/F.oxysporum_fsp_cepae/125/125_assembly.71/sorted_contigs.fa 
qsub /home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking/rep_modeling.sh assembly/velvet/F.oxysporum_fsp_cepae/55/55_assembly.61/sorted_contigs.fa 
qsub /home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking/rep_modeling.sh assembly/velvet/F.oxysporum_fsp_cepae/A23/A23_assembly.61/sorted_contigs.fa 
qsub /home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking/rep_modeling.sh assembly/velvet/F.oxysporum_fsp_cepae/A28/A28_assembly.61/sorted_contigs.fa 
qsub /home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking/rep_modeling.sh assembly/velvet/F.oxysporum_fsp_cepae/D2/D2_assembly.41/sorted_contigs.fa 
qsub /home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking/rep_modeling.sh assembly/velvet/F.oxysporum_fsp_cepae/Fus2/Fus2_assembly.81/sorted_contigs.fa 
qsub /home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking/rep_modeling.sh assembly/velvet/F.oxysporum_fsp_cepae/HB17/HB17_assembly.51/sorted_contigs.fa 
qsub /home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking/rep_modeling.sh assembly/velvet/F.oxysporum_fsp_cepae/PG/PG_assembly.81/sorted_contigs.fa 

# waited until job finished
# then ran

grep -v '#' repeat_masked/F.oxysporum_fsp_cepae/125/125_assembly.71_repmask/125_contigs_hardmasked.fa.out.gff > repeat_masked/F.oxysporum_fsp_cepae/125/125_assembly.71_repmask/125_contigs_hardmasked.gff
/home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking/fix_rm_gff3.py -i repeat_masked/F.oxysporum_fsp_cepae/125/125_assembly.71_repmask/125_contigs_hardmasked.gff -o repeat_masked/F.oxysporum_fsp_cepae/125/125_assembly.71_repmask/125_contigs_hardmasked_fix.gff

qsub /home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking/transposonPSI.sh assembly/velvet/F.oxysporum_fsp_cepae/125/125_assembly.71/sorted_contigs.fa 
qsub /home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking/transposonPSI.sh assembly/velvet/F.oxysporum_fsp_cepae/55/55_assembly.61/sorted_contigs.fa 
qsub /home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking/transposonPSI.sh assembly/velvet/F.oxysporum_fsp_cepae/A23/A23_assembly.61/sorted_contigs.fa 
qsub /home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking/transposonPSI.sh assembly/velvet/F.oxysporum_fsp_cepae/A28/A28_assembly.61/sorted_contigs.fa 
qsub /home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking/transposonPSI.sh assembly/velvet/F.oxysporum_fsp_cepae/D2/D2_assembly.41/sorted_contigs.fa 
qsub /home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking/transposonPSI.sh assembly/velvet/F.oxysporum_fsp_cepae/Fus2/Fus2_assembly.81/sorted_contigs.fa 
qsub /home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking/transposonPSI.sh assembly/velvet/F.oxysporum_fsp_cepae/HB17/HB17_assembly.51/sorted_contigs.fa 
qsub /home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking/transposonPSI.sh assembly/velvet/F.oxysporum_fsp_cepae/PG/PG_assembly.81/sorted_contigs.fa 
