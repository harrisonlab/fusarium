Editors requested that figure 1 from the FoC paper was edited to show isolate
names for published genomes.


The nexus file was edited using sed commands:

These commands were performed on my local machine.

```bash
WorkDir=/Users/armita/Google\ Drive/Fusarium\ comparative\ genomics/old/Fusarium\ comparative\ genomics\ paper/Archive/Figures/Development/Fig1
cd $WorkDir
CurrentNexus=$(ls 30loci_1471365593665_summary_pruned_bs_rooted_work.tree)
NewNexus="30loci_1471365593665_summary_pruned_bs_rooted_parsed.tree"
cat $CurrentNexus \
| sed 's/F. fujikuroi/F. fujikuroi B14/g' \
| sed 's/F. oxysporum f. sp. conglutinans race 2/F. oxysporum f. sp. conglutinans race 2 54008/g' \
| sed 's/F. oxysporum f. sp. cubense race 1/F. oxysporum f. sp. cubense race 1 Foc1/g' \
| sed 's/F. oxysporum f. sp. lycopersici/F. oxysporum f. sp. lycopersici MN25/g' \
| sed 's/F. oxysporum f. sp. melonis/F. oxysporum f. sp. melonis 26406/g' \
| sed 's/F. oxysporum f. sp. pisi/F. oxysporum f. sp. pisi HDV247/g' \
| sed 's/F. oxysporum f. sp. radices-lycopersici/F. oxysporum f. sp. radicis-lycopersici 26381/g' \
| sed 's/F. oxysporum f. sp. raphani/F. oxysporum f. sp. raphani 54005/g' \
| sed 's/F. oxysporum f. sp. vasinifectum/F. oxysporum f. sp. vasinifectum 25433/g' \
| sed 's/F. oxysporum fo47/F. oxysporum Fo47/g' \
| sed 's/F. verticillioides/F. verticillioides 7600/g' \
> $NewNexus


```
