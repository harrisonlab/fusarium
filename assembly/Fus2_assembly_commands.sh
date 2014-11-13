qsub /home/armita/git_repos/emr_repos/scripts/fusarium/assembly/Fus2_assembly.sh qc_dna/paired/F.oxysporum_fsp_cepae/Fus2/Fus2_70bp_F.fastq qc_dna/paired/F.oxysporum_fsp_cepae/Fus2/Fus2_70bp_R.fastq qc_dna/paired/F.oxysporum_fsp_cepae/Fus2/Fus2_250_F.fastq qc_dna/paired/F.oxysporum_fsp_cepae/Fus2/Fus2_250_R.fastq

for folder in $(ls -d Fus2_combined_*); do /home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_stats.pl $folder/process_contigs/sorted_contigs.fa | tail -n 3 | cut -d ':' -f2; done

grep 'Final graph' Fus2_assembly.sh.o6236768 | cut -d ' ' -f 15 | cut -d '/' -f1 