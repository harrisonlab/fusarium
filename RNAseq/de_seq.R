
#===============================================================================
#       Load libraries
#===============================================================================

library(ggplot2)
library(Rsubread)
library(limma)
library(DESeq2)
library(Biostrings)
library(GenomicRanges)
library(GenomicAlignments)


#==========================================================================================
#       Make featureCount compatible gene list form gff file
#				featureCounts takes input in the format:
#				gene_id	start	end	strand
#				It can also use the gtf format - but call to featureCounts would need to be changed
##=========================================================================================

bf <- BamFile(gene_pred/final_genes/F.oxysporum_fsp_cepae/Fus2_canu_new/final/final_genes_appended.gff3, asMates=TRUE)

genes <- read.table("gene_pred/final_genes/F.oxysporum_fsp_cepae/Fus2_canu_new/final/final_genes_appended.gff3",header=F,sep="\t") ###augustus gene models
names(genes) <- c("chr","source","feature","start","end","score","strand","frame","attrib")
genes.pos <- as.data.frame(cbind(GeneID=sub(";","",sub("ID=","",genes[genes$feature=="gene",9])),chr=as.character(genes[genes$feature=="gene",1]),start=genes[genes$feature=="gene",4],end=genes[genes$feature=="gene",5],Strand=genes[genes$feature=="gene",7]))

#===============================================================================
#       Load data and count features
#===============================================================================
targets <- readTargets("alignment/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_72hrs_rep1/accepted_hits.bam")
targets <- readTargets()
fc <- featureCounts(files=targets$filename,annot.ext=genes.pos,useMetaFeatures=TRUE,allowMultiOverlap=TRUE,isPairedEnd=TRUE,nthreads=8,readExtension5=4,readExtension3=4,reportReads=TRUE)

#===============================================================================
#       Data filtering
# 			not necessary (according to the DESeq paper anyway)
#===============================================================================


#===============================================================================
#       DESeq2 analysis
#				Set alpha to the required significance level. This also effects how
#				DESeq calculated FDR - setting to 0.05 and then extracting results with a
#				significance below 0.01 will give slightly different results form setting
#				alpha to 0.01
#================================================================================

countData <- fc$counts
colData <- targets
rownames(colData) <- colData[,1]
colData <- colData[,-1]
colData$id <- factor(paste(colData$sample,colData$condition,sep=""))
myobj <- list(countData=countData,colData=colData)

dds <- 	DESeqDataSetFromMatrix(myobj$countData,myobj$colData,~condition)
sizeFactors(dds) <- sizeFactors(estimateSizeFactors(dds))
dds <- collapseReplicates(dds,groupby=dds$id)
dds <- DESeq(dds, fitType="local")

alpha <- 0.01
res.2 = results(dds, alpha=alpha,contrast=c("condition","RH1","RH2"))
res.3 = results(dds, alpha=alpha,contrast=c("condition","RH1","RH3"))
res.4 = results(dds, alpha=alpha,contrast=c("condition","RH1","RH4"))
res.5 = results(dds, alpha=alpha,contrast=c("condition","RH1","RH5"))
res.6 = results(dds, alpha=alpha,contrast=c("condition","RH1","RH6"))
res.7 = results(dds, alpha=alpha,contrast=c("condition","RH1","RH7"))
res.8 = results(dds, alpha=alpha,contrast=c("condition","RH1","RH8"))

sig.res.2 <- subset(res.2,padj<=alpha)
sig.res.3 <- subset(res.3,padj<=alpha)
sig.res.4 <- subset(res.4,padj<=alpha)
sig.res.5 <- subset(res.5,padj<=alpha)
sig.res.6 <- subset(res.6,padj<=alpha)
sig.res.7 <- subset(res.7,padj<=alpha)
sig.res.8 <- subset(res.8,padj<=alpha)

sig.res.2 <- sig.res.2[order(sig.res.2$padj),]
sig.res.3 <- sig.res.3[order(sig.res.3$padj),]
sig.res.4 <- sig.res.4[order(sig.res.4$padj),]
sig.res.5 <- sig.res.5[order(sig.res.5$padj),]
sig.res.6 <- sig.res.6[order(sig.res.6$padj),]
sig.res.7 <- sig.res.7[order(sig.res.7$padj),]
sig.res.8 <- sig.res.8[order(sig.res.8$padj),]

write.csv(sig.res.2,"02793.csv",quote=F)
write.csv(sig.res.3,"F55.csv",quote=F)
write.csv(sig.res.4,"10170.csv",quote=F)
write.csv(sig.res.5,"MWT.csv",quote=F)
write.csv(sig.res.6,"MOL.csv",quote=F)
write.csv(sig.res.7,"MKO.csv",quote=F)
write.csv(sig.res.8,"TJ.csv",quote=F)

#===============================================================================
#       Get FASTA for significant genes
#===============================================================================
####get genes
mygenes <- readDNAStringSet("/home/groups/harrisonlab/project_files/fusarium_venenatum/gene_pred/final/F.venenatum/strain1/final/final_genes_Braker.gene.fasta")

s2.genes <- mygenes[rownames(sig.res.2)]
s3.genes <- mygenes[rownames(sig.res.3)]
s4.genes <- mygenes[rownames(sig.res.4)]
s5.genes <- mygenes[rownames(sig.res.5)]
s6.genes <- mygenes[rownames(sig.res.6)]
s7.genes <- mygenes[rownames(sig.res.7)]
s8.genes <- mygenes[rownames(sig.res.8)]

writeXStringSet(s2.genes,"02793.fa")
writeXStringSet(s3.genes,"F55.fa")
writeXStringSet(s4.genes,"10170.fa")
writeXStringSet(s5.genes,"MWT.fa")
writeXStringSet(s6.genes,"MOL.fa")
writeXStringSet(s7.genes,"MKO.fa")
writeXStringSet(s8.genes,"TJ.fa")


#===============================================================================
#       FPKM
#===============================================================================
library(naturalsort)
t1 <- counts(dds)
t1 <- mygenes[rownames(t1)]
rowRanges(dds) <- GRanges(t1@ranges@NAMES,t1@ranges)
myfpkm <- fpkm(dds,robust=T)
myfpkm <- cbind(myfpkm, t1@ranges@width)
rm(t1)
myfpkm <- myfpkm[naturalsort(rownames(myfpkm)),]

#===============================================================================
#       Heirachical clustering
#===============================================================================

clus <- function(X,clusters=10,m=1,name="hclust.pdf") {
	if (m==1) {d <- dist(X, method = "manhattan")}
	else if (m==2) {d <- dist(X, method = "euclidean")}
	else if (m==3) {d <- dist(X, method = "maximum")}
	else if (m==4) {d <- dist(X, method = "canberra")}
	else if (m==5) {d <- dist(X, method = "binary")}
	else d <- {dist(X, method = "minkowski")}
	hc <- hclust(d, method="ward")
	groups <- cutree(hc, k=clusters) # cut tree into n clusters
	pdf(name,height=8,width=8)
	plot(hc)
	rect.hclust(hc,k=clusters)
	dev.off()
	return(list(hc,groups,d))
}

#===============================================================================
#       Graphs
#===============================================================================

plotPCAWithLabels <- function (object, intgroup = "condition", ntop = 500,pcx = 1,pcy = 2, returnData = FALSE)
{
    suppressPackageStartupMessages(require(genefilter))
    suppressPackageStartupMessages(require(ggplot2))
    rv <- rowVars(assay(object))
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,
        length(rv)))]
    pca <- prcomp(t(assay(object)[select, ]))
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    if (!all(intgroup %in% names(object@colData))) {
        stop("the argument 'intgroup' should specify columns of colData")
    }
    intgroup.df <- as.data.frame(object@colData[, intgroup,
        drop = FALSE])
    group <- if (length(intgroup) > 1) {
        factor(apply(intgroup.df, 1, paste, collapse = " : "))
    }
    else {
        colData(object)[[intgroup]]
    }
    d <- data.frame(PC1 = pca$x[, pcx], PC2 = pca$x[, pcy], media = group,
        intgroup.df, name = object$label)
    if (returnData) {
        attr(d, "percentVar") <- percentVar[1:2]
        return(d)
    }

    ggplot() +
    geom_point(data=d, mapping=aes(x=PC1, y=PC2, colour=media),size=3) +
    geom_text(data=d, mapping=aes(x=PC1, y=PC2, label=name,colour=media), size=3, vjust=2, hjust=0.5) +
    xlab(paste0("PC",pcx,": ", round(percentVar[pcx] * 100), "% variance")) +
    ylab(paste0("PC",pcy,": ", round(percentVar[pcy] * 100), "% variance"))
}


rld <- varianceStabilizingTransformation(dds,blind=F,fitType="local")
rld$label <- dds$sample
rld$condition <- c("02780","02793","F55","10170","MWT","MOL","MKO","TJ","02780","02793","F55","10170","MWT","MOL","MKO","TJ","02780","02793","F55","10170","MWT","MOL","MKO","TJ")
pdf("quorn.pca_2.pdf",height=10,width=10)
plotPCAWithLabels(rld)
dev.off()
