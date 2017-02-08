

```bash
cat analysis/expression/fo47_expressed_genes.tsv analysis/expression/Fus2_expressed_genes.tsv | cut -f 1-5,7,17,18,19,24 > expressed2.tsv
```

```
expressed2 <- read.delim("~/cluster_mount/groups/harrisonlab/project_files/fusarium/expressed2.tsv", header=FALSE)
View(expressed2)
expressed2 = within(expressed2, {organism = ifelse(grepl("Supercontig*", V6), "Fo", "FoC")})
expressed2 = within(expressed2, {LS = ifelse(grepl("contig_10_pilon|contig_14_pilon|contig_16_pilon|contig_19_pilon|contig_20_pilon|contig_21_pilon|contig_22_pilon", V6), "LS", "Core")})
expressed2$log_expression = log(expressed2$V3)

fpkm_100_genes <- subset(expressed2, V3 >= 100)
qqnorm(fpkm_100_genes$log_expression)
boxplot(log_expression~organism*LS,data=fpkm_100_genes, xlab="Genomic region", ylab="fpkm")
# boxplot(log_expression~V5,data=fpkm_100_genes, xlab="Genomic region", ylab="fpkm")
fpkm_100_genes = within(fpkm_100_genes, {MIMP = ifelse(grepl("Yes", V8), "MIMP", "No mimp")})
boxplot(log_expression~organism*MIMP*LS,data=fpkm_100_genes, xlab="Genomic region", ylab="fpkm")
#boxplot(log_expression~organism*MIMP*LS,data=fpkm_100_genes, xaxt = "n", xlab="", ylab="fpkm")
#axis(1, labels = FALSE)
#text(x =  seq_along(labels), y = par("usr")[3] - 1, srt = 45, adj = 1,
#     labels = labels, xpd = TRUE)
boxplot1$FoC_core <- subset(fpkm_100_genes$V3, fpkm_100_genes$contig == 'core'*fpkm_100_genes$organism == 'FoC')
FoC_LS <- subset(fpkm_100_genes, grepl('FoC', organism) & grepl('LS', LS))
FoC_Core <- subset(fpkm_100_genes, grepl('FoC', organism) & grepl('Core', LS))
Fo_Core <- subset(fpkm_100_genes, grepl('Fo', organism) & grepl('Core', LS))
Expression$FoC_LS <- FoC_LS$
boxplot(FoC_LS$log_expression, xlab="Genomic region", ylab="fpkm")
boxplot(Fo_Core$log_expression FoC_Core$log_expression FoC_LS$log_expression, xlab="Genomic region", ylab="fpkm")


fpkm_100_genes$treatment <- paste(fpkm_100_genes$contig, fpkm_100_genes$LS)
#fpkm_100_genes <- factor(fpkm_100_genes$log_expression, levels=c("Fo Core", "FoC Core", "FoC LS"), fpkm_100_genes$treatment)
boxplot(log_expression~treatment,data=fpkm_100_genes, xlab="Genomic region", ylab="fpkm")

fpkm_100_genes$treatment_mimp <- paste(fpkm_100_genes$treatment, fpkm_100_genes$MIMP)
#fpkm_100_genes <- factor(fpkm_100_genes$log_expression, levels=c("Fo Core No mimp", "Fo Core MIMP", "FoC Core No mimp", "FoC Core MIMP", "FoC LS No mimp", "FoC LS MIMP"), fpkm_100_genes$treatment_mimp)
boxplot(log_expression~treatment_mimp,data=fpkm_100_genes, xlab="Genomic region", ylab="fpkm")
boxplot(log_expression~treatment_mimp,data=fpkm_100_genes, xaxt="n", xlab="", ylab="fpkm")
axis(1, labels = FALSE)
labs <- paste(names(table(fpkm_100_genes$treatment_mimp)))
text(x =  seq_along(labs), y = par("usr")[3] - 0.4, srt = 45, adj = 1,
     labels = labs, xpd = TRUE, cex=0.7)
text(cex=0.75, x=labels-.25, y=-1.25, labs, xpd=TRUE, srt=45, pos=2)

FoC_All <- subset(fpkm_100_genes, grepl('FoC', organism))
FoC_All$V6 <- factor(FoC_All$V6)
boxplot(log_expression~V6,data=FoC_All, xlab="contig", ylab="fpkm")

boxplot(FoC_All~V6,data=FoC_All, xaxt="n", xlab="", ylab="fpkm")
labs <- paste(names(table(FoC_All$V6)))
text(x =  seq_along(labs), y = par("usr")[3] - 0.4, srt = 45, adj = 1,
     labels = labs, xpd = TRUE, cex=0.7)
text(cex=0.75, x=labels-.25, y=-1.25, labs, xpd=TRUE, srt=45, pos=2)
```
