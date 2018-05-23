Editors requested that figure 6 from the FoC paper was edited to show isolate
names for published genomes. This was also taken as an oppertunity to reformat
FoC gene IDs, using those gene IDs held at Ncbi rather internally.


The newick file was edited using sed commands:

These commands were performed on my local machine.

```bash
WorkDir=/Users/armita/Google\ Drive/Fusarium\ comparative\ genomics/old/Fusarium\ comparative\ genomics\ paper/Archive/Figures/Development/Fig6
cd $WorkDir
CurrentNewick=$(ls FTF_NJ_tree3.newick)
NewNewick="FTF_NJ_parsed.newick"
cat $CurrentNewick \
| sed 's/FoM_/FoM_26406|/g' \
| sed 's/FoV_/FoV_25433|/g' \
| sed 's/FoL_/FoL_4287|/g' \
| sed 's/FoR_/FoR_54005|/g' \
| sed 's/FoCub_/FoCub_54006|/g' \
| sed 's/FoCu_/FoCub_54006|/g' \
| sed 's/FoC_/FoCub_54006|/g' \
| sed 's/FTF1b_FoCub_54006|JH658347_FOIG_16630/FoCub_54006|FTF1b_FOIG_16630_JH658347/g' \
| sed 's/Fo_/Fo_FOSC_3-a|/g' \
| sed 's/FoRL_/FoRL_26381|/g' \
| sed 's/FoPh_/FoPh_FOP-SP4|/g' \
| sed 's/fo47|/Fo_Fo47|/g' \
| sed 's/Fo47_/Fo_Fo47|/g' \
| sed 's/Fus2|g16859/FoC_Fus2|BFJ65_g18330/g' \
| sed 's/Fus2|NS_02884/FoC_Fus2|BFJ65_g18276/g' \
| sed 's/Fus2|g10474/FoC_Fus2|BFJ65_g11825/g' \
| sed 's/125|g16201/FoC_125|BFJ66_g17660/g' \
| sed 's/125|g2599/FoC_125|BFJ66_g2789/g' \
| sed 's/A23|g16099/FoC_A23|BFJ67_g17550/g' \
| sed 's/A23|g3017/FoC_A23|BFJ67_g3310/g' \
| sed 's/A13|g15237/Fo_A13|BFJ69_g15942/g' \
| sed 's/A13|g9338/Fo_A13|BFJ69_g9724/g' \
| sed 's/PG|g14993/Fo_PG|BFJ70_g16057/g' \
| sed 's/PG|g9482/Fo_PG|BFJ70_g10160/g' \
| sed 's/CB3|g6501/Fo_CB3|BFJ71_g6993/g' \
| sed 's/A28|g1088/Fo_A28|BFJ68_g839/g' \
> $NewNewick
```
