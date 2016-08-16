#!/usr/bin/Rscript

# Plot a 3-way Venn diagram from a tab delimited file containing a matrix showing
 # presence /absence of orthogroups between 3 sets of genomes.

 # This is intended to be used on the output of the orthoMCL pipeline following
 # building of the matrix using:
 # ~/git_repos/emr_repos/tools/pathogen/orthology/orthoMCL/orthoMCLgroups2tab.py

 # The script also requires the colorspace package. This can be downloaded by
 # opening R and running the following command:
 # options(download.file.method = "wget")
 # install.packages("colorspace")

 #get config options
library(optparse)
library(colorspace)
library(VennDiagram)
opt_list = list(
    make_option("--inp", type="character", help="tab seperated file containing matrix of presence of orthogroups"),
    make_option("--out", type="character", help="output venn diagram in pdf format")
)
opt = parse_args(OptionParser(option_list=opt_list))
f = opt$inp
o = opt$out

orthotabs <-data.frame()
orthotabs <- read.table(f)
df1 <- t(orthotabs)
summary(df1)


nonpath=subset(df1, df1[,"A28"] == 1 & df1[,"CB3"] == 1 & df1[,"PG"] == 1 & df1[,"fo47"] == 1 & df1[,"A1_2"] == 1 & df1[,"Fus2"] == 0 & df1[,"125"] == 0 & df1[,"A23"] == 0 & df1[,"4287"] == 0)
path=subset(df1, df1[,"A28"] == 0 & df1[,"CB3"] == 0 & df1[,"PG"] == 0 & df1[,"fo47"] == 0 & df1[,"A1_2"] == 0 & df1[,"Fus2"] == 1 & df1[,"125"] == 1 & df1[,"A23"] == 1 & df1[,"4287"] == 0)
tomato=subset(df1, df1[,"A28"] == 0 & df1[,"CB3"] == 0 & df1[,"PG"] == 0 & df1[,"fo47"] == 0 & df1[,"A1_2"] == 0 & df1[,"Fus2"] == 0 & df1[,"125"] == 0 & df1[,"A23"] == 0 & df1[,"4287"] == 1)
# orthologs=subset(df1, df1[,"A28"] == 1 & df1[,"CB3"] == 1 & df1[,"PG"] == 1 & df1[,"fo47"] == 1 & df1[,"A1_2"] == 1 & df1[,"Fus2"] == 1 & df1[,"125"] == 1 & df1[,"A23"] == 1 & df1[,"4287"] == 1)

# area1=(nrow(nonpath) + nrow(orthologs))
# area2=(nrow(path) + nrow(orthologs))
# area3=(nrow(tomato) + nrow(orthologs))
# area3
# area2
# area1



# Print labels
# label1 <- paste('Path', ' (', area2, ')', sep="" )
# label2 <- paste('NonPath', ' (', area1, ')', sep="" )
# label3 <- paste('FoL', ' (', area3, ')', sep="" )

# No labels
label1 <- paste("", sep="" )
label2 <- paste("", sep="" )
label3 <- paste("", sep="" )

n123=nrow(subset(df1, df1[,"A28"] == 1 & df1[,"CB3"] == 1 & df1[,"PG"] == 1 & df1[,"fo47"] == 1 & df1[,"A1_2"] == 1 & df1[,"Fus2"] == 1 & df1[,"125"] == 1 & df1[,"A23"] == 1 & df1[,"4287"] == 1))
n12=n123 + nrow(subset(df1, df1[,"A28"] == 1 & df1[,"CB3"] == 1 & df1[,"PG"] == 1 & df1[,"fo47"] == 1 & df1[,"A1_2"] == 1 & df1[,"Fus2"] == 1 & df1[,"125"] == 1 & df1[,"A23"] == 1 & df1[,"4287"] == 0))
n13=n123 + nrow(subset(df1, df1[,"A28"] == 1 & df1[,"CB3"] == 1 & df1[,"PG"] == 1 & df1[,"fo47"] == 1 & df1[,"A1_2"] == 1 & df1[,"Fus2"] == 0 & df1[,"125"] == 0 & df1[,"A23"] == 0 & df1[,"4287"] == 1))
n23=n123 + nrow(subset(df1, df1[,"A28"] == 0 & df1[,"CB3"] == 0 & df1[,"PG"] == 0 & df1[,"fo47"] == 0 & df1[,"A1_2"] == 0 & df1[,"Fus2"] == 1 & df1[,"125"] == 1 & df1[,"A23"] == 1 & df1[,"4287"] == 1))
summary(n12)
summary(n13)
summary(n23)
summary(n123)

area1=(nrow(nonpath) + (n12 - n123) + (n13 - n123) + n123)
area2=(nrow(path) + (n12 - n123) + (n23 - n123) + n123)
area3=(nrow(tomato) + (n13 - n123) + (n23 - n123) + n123)
#nrow(nonpath)
nrow(nonpath)
nrow(path)
nrow(tomato)
n12
n13
n23
n123
area1
#area1 - n12 - n13 + n123
area2
area3

pdf(o)
draw.triple.venn(area1, area2, area3,
    n12, n23, n13,
    n123,
    category = c(label1, label2, label3),
#    rep("", 4),
    rotation = 1,
    reverse = FALSE,
    lwd = rep(2, 3),
    lty = rep("solid", 3),
    col = rep("black", 3),
    fill = c(rainbow_hcl(3)),
    alpha = rep(0.5, 3),
    label.col = rep("black", 7),
    cex = rep(1, 7),
    fontface = rep("plain", 7),
    fontfamily = rep("serif", 7),
    cat.pos = c(-40, 40, 180),
    cat.dist = c(0.05, 0.05, 0.025),
    cat.col = rep("black", 3),
    cat.cex = rep(1, 3),
    cat.fontface = rep("plain", 3),
    cat.fontfamily = rep("serif", 3),
    cat.just = list(c(0.5, 1), c(0.5, 1), c(0.5, 0)),
    cat.default.pos = "outer",
    cat.prompts = FALSE,
    rotation.degree = 0,
    rotation.centre = c(0.5, 0.5),
    ind = TRUE, sep.dist = 0.05, offset = 0,
    )


dev.off()
singles = df1[grepl("single*", rownames(df1)), ]
print("A28")
total_1 = nrow(subset (df1, df1[,"A28"] == 1))
missing_1 = (total_1 - area2)
uniq_1=sum(singles[, "A28"])
paste('The total number of orthogroups and singleton genes in this isolate: ', total_1)
paste('The total number of orthogroups and singleton genes not in the venn diagram: ', missing_1)
paste('The total number of singleton genes not in the venn diagram: ', uniq_1)
print("CB3")
total_2 = nrow(subset (df1, df1[,"CB3"] == 1))
missing_2 = (total_2 - area2)
uniq_2=sum(singles[, "CB3"])
paste('The total number of orthogroups and singleton genes in this isolate: ', total_2)
paste('The total number of orthogroups and singleton genes not in the venn diagram: ', missing_2)
paste('The total number of singleton genes not in the venn diagram: ', uniq_2)
print("PG")
total_3 = nrow(subset (df1, df1[,"PG"] == 1))
missing_3 = (total_3 - area2)
uniq_3=sum(singles[, "PG"])
paste('The total number of orthogroups and singleton genes in this isolate: ', total_3)
paste('The total number of orthogroups and singleton genes not in the venn diagram: ', missing_3)
paste('The total number of singleton genes not in the venn diagram: ', uniq_3)
print("Fus2")
total_4 = nrow(subset (df1, df1[,"Fus2"] == 1))
missing_4 = (total_4 - area1)
uniq_4=sum(singles[, "Fus2"])
paste('The total number of orthogroups and singleton genes in this isolate: ', total_4)
paste('The total number of orthogroups and singleton genes not in the venn diagram: ', missing_4)
paste('The total number of singleton genes not in the venn diagram: ', uniq_4)
print("125")
total_5 = nrow(subset (df1, df1[,"125"] == 1))
missing_5 = (total_5 - area1)
uniq_5=sum(singles[, "125"])
paste('The total number of orthogroups and singleton genes in this isolate: ', total_5)
paste('The total number of orthogroups and singleton genes not in the venn diagram: ', missing_5)
paste('The total number of singleton genes not in the venn diagram: ', uniq_5)
print("A23")
total_6 = nrow(subset (df1, df1[,"A23"] == 1))
missing_6 = (total_6 - area1)
uniq_6=sum(singles[, "A23"])
paste('The total number of orthogroups and singleton genes in this isolate: ', total_6)
paste('The total number of orthogroups and singleton genes not in the venn diagram: ', missing_6)
paste('The total number of singleton genes not in the venn diagram: ', uniq_6)
print("4287")
total_7 = nrow(subset (df1, df1[,"4287"] == 1))
missing_7 = (total_6 - area3)
uniq_7=sum(singles[, "4287"])
paste('The total number of orthogroups and singleton genes in this isolate: ', total_6)
paste('The total number of orthogroups and singleton genes not in the venn diagram: ', missing_6)
paste('The total number of singleton genes not in the venn diagram: ', uniq_6)

#inpara_2 = sum(orthogroups[,"A28"] == 0 & orthogroups[,"CB3"] == 1)
#label1
#uniq_1
#inpara_1
#label2
#uniq_2
#inpara_2

warnings()
q()
