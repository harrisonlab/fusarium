qsub /home/armita/git_repos/emr_repos/scripts/fusarium/assembly/Fus2_assembly.sh qc_dna/paired/F.oxysporum_fsp_cepae/Fus2/Fus2_70bp_F.fastq qc_dna/paired/F.oxysporum_fsp_cepae/Fus2/Fus2_70bp_R.fastq qc_dna/paired/F.oxysporum_fsp_cepae/Fus2/Fus2_250_F.fastq qc_dna/paired/F.oxysporum_fsp_cepae/Fus2/Fus2_250_R.fastq


for directory in $(ls -d assembly/velvet/F.oxysporum_fsp_cepae/Fus2/Fus2_combined_*); do process_contigs.pl -i $directory/contigs.fa -o $directory/process_contigs; done

for directory in $(ls -d assembly/velvet/F.oxysporum_fsp_cepae/Fus2/Fus2_combined_*); do cat $directory/process_contigs/stats.txt | rev | cut -d ' ' -f1 | rev; done

for directory in $(ls -d assembly/velvet/F.oxysporum_fsp_cepae/Fus2/Fus2_combined_*); do /home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_stats.pl $directory/process_contigs/sorted_contigs.fa | tail -n 3 | cut -d ':' -f2; done

grep 'Final graph' Fus2_assembly.sh.o6236768 | cut -d ' ' -f 15 | cut -d '/' -f1 

qsub /home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking/rep_modeling.sh assembly/velvet/F.oxysporum_fsp_cepae/Fus2/Fus2_combined_49/process_contigs/sorted_contigs.fa

/home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking/fix_rm_gff3.py -i repeat_masked/F.oxysporum_fsp_cepae/Fus2/Fus2_combined_49/Fus2_contigs_hardmasked.fa.out.gff -o repeat_masked/F.oxysporum_fsp_cepae/Fus2/Fus2_combined_49/Fus2_contigs_hardmasked.gff

qsub /home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking/transposonPSI.sh assembly/velvet/F.oxysporum_fsp_cepae/Fus2/Fus2_combined_49/process_contigs/sorted_contigs.fa

SCRIPT_DIR=/home/armita/git_repos/emr_repos/tools/gene_prediction/augustus
qsub $SCRIPT_DIR/submit_augustus.sh fusarium_graminearum repeat_masked/F.oxysporum_fsp_cepae/Fus2/Fus2_combined_49/Fus2_combined_49_contigs_hardmasked.fa
