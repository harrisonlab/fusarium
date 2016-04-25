

```bash
  for Strain in 125 A23 Fus2_edited_v2 55 A1-2 CB3 HB6 A13 A28 D2 PG fo47 4287; do
    for GeneGff in $(ls gene_pred/codingquary/*/$Strain/final/final_genes_appended.gff3 | grep -v -e "fo47" -e "4287" | grep 'Fus2_edited_v2'); do
      Organism=$(echo $GeneGff | rev | cut -f4 -d '/' | rev)
      Strain=$(echo $GeneGff | rev | cut -f3 -d '/' | rev)
      Genome=$(ls repeat_masked/$Organism/$Strain/filtered_contigs_repmask/*_contigs_unmasked.fa)
      GeneGff=$(ls gene_pred/codingquary/$Organism/$Strain/final/final_genes_appended.gff3)
      SigpTab=$(ls gene_pred/final_genes_signalp-4.1/$Organism/$Strain/*_final_sp.tab)
      TmhmmTxt=$(ls gene_pred/trans_mem/$Organism/$Strain/*_tmhmm_out.txt)
      MimpTxt=$(ls analysis/mimps/$Organism/$Strain/*_genes_in_2kb_mimp.txt)
      EffectorpTxt=$(ls analysis/effectorP/$Organism/$Strain/*_EffectorP.txt)
      # OrthogroupsTxt=$(ls analysis/orthology/orthomcl/FoC_path_vs_non_path/FoC_path_vs_non_path_orthogroups.txt)
      OrthogroupsTxt=$(ls analysis/orthology/orthomcl/FoC_vs_Fo_vs_FoL/FoC_vs_Fo_vs_FoL_orthogroups.txt)

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
      > $OutTable
    done
  done


```
