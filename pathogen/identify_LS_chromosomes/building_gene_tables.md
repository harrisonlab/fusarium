

```bash
  for Strain in 125 A23 Fus2_edited_v2 55 A1-2 CB3 HB6 A13 A28 D2 PG; do
    for GeneGff in $(ls gene_pred/codingquary/*/$Strain/final/final_genes_appended.gff3); do
      Organism=$(echo $GeneGff | rev | cut -f4 -d '/' | rev)
      Strain=$(echo $GeneGff | rev | cut -f3 -d '/' | rev)
      echo "$Organism - $Strain"
      Genome=$(ls repeat_masked/$Organism/$Strain/filtered_contigs_repmask/*_contigs_unmasked.fa)
      BlastCsv=$(ls analysis/blast_homology/$Organism/$Strain/4287_chromosomal_final_genes_combined.pep.fasta_hits.csv)
      FolIntersect=$(ls analysis/blast_homology/$Organism/$Strain/4287_chromosomal_final_genes_combined_intersect.bed)
      GeneGff=$(ls gene_pred/codingquary/$Organism/$Strain/final/final_genes_appended.gff3)
      SigpTab=$(ls gene_pred/final_genes_signalp-4.1/$Organism/$Strain/*_final_sp.tab)
      TmhmmTxt=$(ls gene_pred/trans_mem/$Organism/$Strain/*_tmhmm_out.txt)
      MimpTxt=$(ls analysis/mimps/$Organism/$Strain/*_genes_in_2kb_mimp.txt)
      EffectorpTxt=$(ls analysis/effectorP/$Organism/$Strain/*_EffectorP.txt)
      # OrthogroupsTxt=$(ls analysis/orthology/orthomcl/FoC_path_vs_non_path/FoC_path_vs_non_path_orthogroups.txt)
      OrthogroupsTxt=$(ls analysis/orthology/orthomcl/FoC_vs_Fo_vs_FoL/FoC_vs_Fo_vs_FoL_orthogroups.txt)
      InterProTsv=$(ls gene_pred/interproscan/$Organism/$Strain/*_interproscan.tsv)
      SwissprotTab=$(ls gene_pred/swissprot/$Organism/$Strain/swissprot_v2015_tophit_parsed.tbl)
      DEG_Orthogroups=$(ls analysis/expression/warwick/F.oxysporum_fsp_cepae/Fus2/04_16/Fus2_path_vs_non_path_orthogroups.tab)

      OrthoMCL_id="$Strain"
      OrthoMCL_id_list="125 A23 Fus2 55 A1_2 CB3 HB6 A13 A28 D2 PG fo47 4287"
      OrthoMCL_path_ids="125 A23 Fus2"
      OrthoMCL_nonpath_ids="A13 A28 D2 PG fo47"

      if [ "$Strain" == 'Fus2_edited_v2' ]; then OrthoMCL_id="Fus2"; fi
      if [ "$Strain" == 'A1-2' ]; then OrthoMCL_id="A1_2"; fi

      OutDir=gene_pred/annotations/$Organism/$Strain
      OutTable=$OutDir/"$Strain"_gene_annotations.tab

      mkdir -p $OutDir

      ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes
      $ProgDir/Fo_build_gene_annot_table.py \
      --blast_csv $BlastCsv \
      --FoL_intersected_genes $FolIntersect \
      --genome $Genome \
      --FoC_genes_gff $GeneGff \
      --FoC_SigP $SigpTab \
      --FoC_TM_list $TmhmmTxt \
      --FoC_MIMP_list $MimpTxt \
      --FoC_effectorP $EffectorpTxt \
      --FoC_orthogroup $OrthogroupsTxt \
      --OrthoMCL_id $OrthoMCL_id \
      --OrthoMCL_all $OrthoMCL_id_list \
      --OrthoMCL_path $OrthoMCL_path_ids \
      --OrthoMCL_nonpath $OrthoMCL_nonpath_ids \
      --InterPro $InterProTsv \
      --Swissprot $SwissprotTab \
      --DEG_Orthogroups $DEG_Orthogroups \
      > $OutTable
    done
  done
```

Gene tables were made for Fo fo47 and FoL 4287

```bash
  for Strain in fo47 4287; do
    if [ $Strain == '4287' ]; then
      GeneGff=$(ls assembly/external_group/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl/Fusarium_oxysporum.FO2.31_parsed.gff3)
      Genome=$(ls assembly/external_group/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl/Fusarium_oxysporum.FO2.31.dna.chromosome_parsed.fa)
      BlastCsv=$(ls analysis/blast_homology/F.oxysporum_fsp_lycopersici/4287/4287_chromosomal_Fusox1_GeneCatalog_proteins_20110522_parsed.fa_hits.csv)
      FolIntersect=$(ls analysis/blast_homology/F.oxysporum_fsp_lycopersici/4287/4287_chromosomal_final_genes_combined_intersect.bed)
    elif [ $Strain == 'fo47' ]; then
      GeneGff=$(ls assembly/external_group/F.oxysporum/fo47/broad/fusarium_oxysporum_fo47_1_transcripts_parsed.gff3)
      Genome=$(ls assembly/external_group/F.oxysporum/fo47/broad/fusarium_oxysporum_fo47_1_supercontigs_parsed.fasta)
      BlastCsv=$(ls analysis/blast_homology/F.oxysporum/fo47/4287_chromosomal_final_genes_combined.pep.fasta_hits.csv)
      FolIntersect=$(ls analysis/blast_homology/F.oxysporum/fo47/4287_chromosomal_final_genes_combined_intersect.bed)
    fi
    Organism=$(echo $GeneGff | rev | cut -f4 -d '/' | rev)
    # Strain=$(echo $GeneGff | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    # Genome=$(ls repeat_masked/$Organism/$Strain/filtered_contigs_repmask/*_contigs_unmasked.fa)
    # BlastCsv=$(ls analysis/blast_homology/$Organism/$Strain/4287_chromosomal_final_genes_combined.pep.fasta_hits.csv)
    # FolIntersect=$(ls analysis/blast_homology/$Organism/$Strain/4287_chromosomal_final_genes_combined_intersect.bed)
    # GeneGff=$(ls gene_pred/codingquary/$Organism/$Strain/final/final_genes_appended.gff3)
    SigpTab=$(ls gene_pred/final_genes_signalp-4.1/$Organism/$Strain/*_final_sp.tab)
    TmhmmTxt=$(ls gene_pred/trans_mem/$Organism/$Strain*/*_tmhmm_out.txt)
    MimpTxt=$(ls analysis/mimps/$Organism/$Strain*/*_genes_in_2kb_mimp.txt)
    EffectorpTxt=$(ls analysis/effectorP/$Organism/$Strain/*_EffectorP.txt)
    # OrthogroupsTxt=$(ls analysis/orthology/orthomcl/FoC_path_vs_non_path/FoC_path_vs_non_path_orthogroups.txt)
    OrthogroupsTxt=$(ls analysis/orthology/orthomcl/FoC_vs_Fo_vs_FoL/FoC_vs_Fo_vs_FoL_orthogroups.txt)
    InterProTsv=$(ls gene_pred/interproscan/$Organism/$Strain/*_interproscan.tsv)
    SwissprotTab=$(ls gene_pred/swissprot/$Organism/$Strain/swissprot_v2015_tophit_parsed.tbl)
    DEG_Orthogroups=$(ls analysis/expression/warwick/F.oxysporum_fsp_cepae/Fus2/04_16/Fus2_path_vs_non_path_orthogroups.tab)

    OrthoMCL_id="$Strain"
    OrthoMCL_id_list="125 A23 Fus2 55 A1_2 CB3 HB6 A13 A28 D2 PG fo47 4287"
    OrthoMCL_path_ids="125 A23 Fus2"
    OrthoMCL_nonpath_ids="A13 A28 D2 PG fo47"

    if [ "$Strain" == 'Fus2_edited_v2' ]; then OrthoMCL_id="Fus2"; fi
    if [ "$Strain" == 'A1-2' ]; then OrthoMCL_id="A1_2"; fi

    OutDir=gene_pred/annotations/$Organism/$Strain
    OutTable=$OutDir/"$Strain"_gene_annotations.tab

    mkdir -p $OutDir

    ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes
    $ProgDir/Fo_build_gene_annot_table.py \
    --blast_csv $BlastCsv \
    --FoL_intersected_genes $FolIntersect \
    --genome $Genome \
    --FoC_genes_gff $GeneGff \
    --FoC_SigP $SigpTab \
    --FoC_TM_list $TmhmmTxt \
    --FoC_MIMP_list $MimpTxt \
    --FoC_effectorP $EffectorpTxt \
    --FoC_orthogroup $OrthogroupsTxt \
    --OrthoMCL_id $OrthoMCL_id \
    --OrthoMCL_all $OrthoMCL_id_list \
    --OrthoMCL_path $OrthoMCL_path_ids \
    --OrthoMCL_nonpath $OrthoMCL_nonpath_ids \
    --InterPro $InterProTsv \
    --Swissprot $SwissprotTab \
    --DEG_Orthogroups $DEG_Orthogroups \
    > $OutTable
  done
```
