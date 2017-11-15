#!/usr/bin/Rscript

# Plot a 5-way Venn diagram from a tab delimited file containing a matrix showing
 # presence /absence of orthogroups within 5 isolates.

 # This is intended to be used on the output of the orthoMCL pipeline following
 # building of the matrix using:
 # ~/git_repos/emr_repos/tools/pathogen/orthology/orthoMCL/orthoMCLgroups2tab.py

# This script requires the optparse R package. This can be downloaded by opening
# R and running the following command:
# install.packages("optparse",repos="http://cran.uk.r-project.org")
# When given the option, install this package to a local library.

#get config options
library(optparse)
library(colorspace)
library(VennDiagram, lib.loc="/home/armita/R-packages/")
opt_list = list(
    make_option("--inp", type="character", help="tab seperated file containing matrix of presence of orthogroups"),
    make_option("--out", type="character", help="output venn diagram in pdf format")
#    make_option("--maxrf", type="double", default=0.2, help="max rf to consider as linked"),
#    make_option("--minlod", type="double", default=20.0, help="min LOD to consider as linked")
)
opt = parse_args(OptionParser(option_list=opt_list))
f = opt$inp
o = opt$out

orthotabs <-data.frame()
orthotabs <- read.table(f)
df1 <- t(orthotabs)
summary(df1)


area1=sum(df1[,"FoN"])
area2=sum(df1[,"fo47"])
area3=sum(df1[,"FoC"])
area4=sum(df1[,"FoL"])

# print(area1, area2, area3, area4)

colname1 <- paste("FoN")
colname2 <- paste("fo47")
colname3 <- paste("FoC")
colname4 <- paste("FoL")

#label1 <- paste(colname1, ' (', area1, ')', sep="" )
#label2 <- paste(colname2, ' (', area2, ')', sep="" )
#label3 <- paste(colname3, ' (', area3, ')', sep="" )
#label4 <- paste(colname4, ' (', area4, ')', sep="" )
label1 <- paste("FON", ' (', area1, ')', sep="" )
label2 <- paste("FO", ' (', area2, ')', sep="" )
label3 <- paste("FOC", ' (', area3, ')', sep="" )
label4 <- paste("FOL", ' (', area4, ')', sep="" )

n12=nrow(subset(df1, df1[,"FoN"] == 1 & df1[,"fo47"] == 1))
n13=nrow(subset(df1, df1[,"FoN"] == 1 & df1[,"FoC"] == 1))
n14=nrow(subset(df1, df1[,"FoN"] == 1 & df1[,"FoL"] == 1))
n23=nrow(subset(df1, df1[,"fo47"] == 1 & df1[,"FoC"] == 1))
n24=nrow(subset(df1, df1[,"fo47"] == 1 & df1[,"FoL"] == 1))
n34=nrow(subset(df1, df1[,"FoC"] == 1 & df1[,"FoL"] == 1))
n123=nrow(subset(df1, df1[,"FoN"] == 1 & df1[,"fo47"] == 1 & df1[,"FoC"] == 1))
n124=nrow(subset(df1, df1[,"FoN"] == 1 & df1[,"fo47"] == 1 & df1[,"FoL"] == 1))
n134=nrow(subset(df1, df1[,"FoN"] == 1 & df1[,"FoC"] == 1 & df1[,"FoL"] == 1))
n234=nrow(subset(df1, df1[,"fo47"] == 1 & df1[,"FoC"] == 1 & df1[,"FoL"] == 1))
n1234=nrow(subset(df1, df1[,"FoN"] == 1 & df1[,"fo47"] == 1 & df1[,"FoC"] == 1 & df1[,"FoL"] == 1))
summary(n12)
summary(n123)
summary(n1234)
pdf(o)
draw.quad.venn(area1, area2, area3, area4,
    n12, n13, n14, n23, n24, n34,
    n123, n124, n134, n234,
    n1234,
    category = c(label1, label2, label3, label4),
    lwd = rep(2, 4),
    lty = rep("solid", 4),
    col = rep("black", 4),
    fill = c(rainbow_hcl(4)),
    alpha = rep(0.5, 4),
    label.col = rep("black", 15),
    cex = rep(1, 15),
    fontface = rep("plain", 15),
    fontfamily = rep("serif", 15),
    cat.pos = c(-15, 15, 0, 0),
    cat.dist = c(0.22, 0.22, 0.11, 0.11),
    cat.col = rep("black", 4),
    cat.cex = rep(1, 4),
    cat.fontface = rep("plain", 4),
    cat.fontfamily = rep("serif", 4),
    cat.just = rep(list(c(0.5, 0.5)), 4),
    rotation.degree = 0,
    rotation.centre = c(0.5, 0.5)
    )

dev.off()
singles = df1[grepl("single*", rownames(df1)), ]
uniq_1=sum(singles[, colname1])
uniq_2=sum(singles[, colname2])
uniq_3=sum(singles[, colname3])
uniq_4=sum(singles[, colname4])
orthogroups = df1[grepl("orthogroup*", rownames(df1)), ]
inpara_1 = sum(orthogroups[,colname1] == 1 & orthogroups[,colname2] == 0 & orthogroups[,colname3] == 0 & orthogroups[,colname4] == 0)
inpara_2 = sum(orthogroups[,colname1] == 0 & orthogroups[,colname2] == 1 & orthogroups[,colname3] == 0 & orthogroups[,colname4] == 0)
inpara_3 = sum(orthogroups[,colname1] == 0 & orthogroups[,colname2] == 0 & orthogroups[,colname3] == 1 & orthogroups[,colname4] == 0)
inpara_4 = sum(orthogroups[,colname1] == 0 & orthogroups[,colname2] == 0 & orthogroups[,colname3] == 0 & orthogroups[,colname4] == 1)
label1
uniq_1
inpara_1
label2
uniq_2
inpara_2
label3
uniq_3
inpara_3
label4
uniq_4
inpara_4

warnings()
q()
