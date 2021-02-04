#!/usr/bin/env Rscript

library(tximport)
library(readr)
library("DESeq2")
library("dplyr")
library("ggplot2")
library("pheatmap")
library("RColorBrewer")
library("apeglm")
library("genefilter")

# Currently set up for old analysis, revise filepaths

##################################################
# READ IN SALMON DATA
##################################################

coldata <- read.delim("R05_samples.tsv", row.names=1, stringsAsFactors=TRUE)
coldata$files <- file.path("../../../03_quantification/VsTranscriptome_Preliminary", coldata$Sample_ID, "quant.sf")
file.exists(coldata$files)
rownames(coldata) <- coldata$Sample_ID
coldata$Dormancy <- as.factor(coldata$Dormancy)
coldata$Dormancy <- relevel(coldata$Dormancy, "0")

##################################################
# SAVE A SPECIFIC GROUP
##################################################

timepoint <- "49"	# Add timepoint of interest here!
veg_rest <- coldata[which(coldata$Dormancy == "0" | coldata$Dormancy == timepoint),]
veg_rest$Dormancy <- factor(veg_rest$Dormancy, c("0",timepoint))

##################################################
# IMPORT TRANSCRIPT ABUNDANCE DATA
##################################################

txi <- tximport(coldata$files, type = "salmon", txOut = TRUE)

txi_vr <- tximport(veg_rest$files, type = "salmon", txOut = TRUE)

##################################################
# GENERATE DESeqDataSet OBJECT
##################################################

dds <- DESeqDataSetFromTximport(txi, coldata, ~Batch + Dormancy)

dds_vr <- DESeqDataSetFromTximport(txi_vr, veg_rest, ~Batch + Dormancy)

##################################################
# PRE-FILTER THE DATASET
##################################################

nrow(dds)		# 21,237
nrow(dds_vr)		# 21,237

## Slightly stricter filter - at least one count in at least three samples
## Three samples per condition, so should retain condition-specific genes

keep <- rowSums(counts(dds) >= 1 ) >= 3
dds <- dds[keep,]

nrow(dds)		# 20,177 (0 vs 49)

keep_vr <- rowSums(counts(dds_vr) >= 1 ) >= 3
dds_vr <- dds_vr[keep_vr,]

nrow(dds_vr)		# 19,164 (0 vs 49)

##################################################
# DATA TRANSFORMATION
# (FOR APPLICATIONS OTHER THAN DIFF. TESTING)
##################################################

## Variance stabilising transformation - recommended for n > 30
vsd <- vst(dds, blind = FALSE)
vsd_blind <- vst(dds, blind = TRUE)

vsd_vr <- vst(dds_vr, blind = FALSE)	# This notes that the dispersion trend was not well-captured
vsd_vr_blind <- vst(dds_vr, blind = TRUE)

## rlog transformation - recommended for n < 30
rld <- rlog(dds, blind = FALSE)
rld_blind <- rlog(dds, blind = TRUE)

rld_vr <- rlog(dds_vr, blind = FALSE)	# This notes that the dispersion trend was not well-captured
rld_vr_blind <- rlog(dds_vr, blind = TRUE)

##################################################
# ASSESS THE SIMILARITY BETWEEN SAMPLES
##################################################

## Heat map

annotation_col <- coldata[,c(4,5)]

all_dist_rld <- dist(t(assay(rld)))
all_dist_rld_blind <- dist(t(assay(rld_blind)))
all_dist_vsd <- dist(t(assay(vsd)))
all_dist_vsd_blind <- dist(t(assay(vsd_blind)))

all_dist_rld_vr <- dist(t(assay(rld_vr)))
all_dist_rld_vr_blind <- dist(t(assay(rld_vr_blind)))
all_dist_vsd_vr <- dist(t(assay(vsd_vr)))
all_dist_vsd_vr_blind <- dist(t(assay(vsd_vr_blind)))

distMatrix_rld <- as.matrix(all_dist_rld)
distMatrix_rld_blind <- as.matrix(all_dist_rld_blind)
distMatrix_vsd <- as.matrix(all_dist_vsd)
distMatrix_vsd_blind <- as.matrix(all_dist_vsd_blind)

distMatrix_rld_vr <- as.matrix(all_dist_rld_vr)
distMatrix_rld_vr_blind <- as.matrix(all_dist_rld_vr_blind)
distMatrix_vsd_vr <- as.matrix(all_dist_vsd_vr)
distMatrix_vsd_vr_blind <- as.matrix(all_dist_vsd_vr_blind)

pheatmap(distMatrix_rld, annotation_col = annotation_col)
pheatmap(distMatrix_rld_blind, annotation_col = annotation_col)
pheatmap(distMatrix_vsd, annotation_col = annotation_col)
pheatmap(distMatrix_vsd_blind, annotation_col = annotation_col)

pheatmap(distMatrix_rld_vr, annotation_col = annotation_col)
pheatmap(distMatrix_rld_vr_blind, annotation_col = annotation_col)
pheatmap(distMatrix_vsd_vr, annotation_col = annotation_col)
pheatmap(distMatrix_vsd_vr_blind, annotation_col = annotation_col)

## PCA

plotPCA(rld, intgroup = c("Dormancy","Batch"))
plotPCA(rld_blind, intgroup = c("Dormancy","Batch"))
plotPCA(vsd, intgroup = c("Dormancy","Batch"))
plotPCA(vsd_blind, intgroup = c("Dormancy","Batch"))

plotPCA(rld_vr, intgroup = c("Dormancy","Batch"))
plotPCA(rld_vr_blind, intgroup = c("Dormancy","Batch"))
plotPCA(vsd_vr, intgroup = c("Dormancy","Batch"))
plotPCA(vsd_vr_blind, intgroup = c("Dormancy","Batch"))

##################################################
# DIFFERENTIAL EXPRESSION ANALYSIS
##################################################

dds <- DESeq(dds)
dds_vr <- DESeq(dds_vr)

res <- results(dds, alpha = 0.01)
res_vr <- results(dds_vr, alpha = 0.01)

hist(res$padj, breaks=100)
hist(res_vr$padj, breaks=100)

plotMA(res, alpha = 0.01)
plotMA(res_vr, alpha = 0.01)

sum(is.na(res$padj))		# 311
sum(is.na(res_vr$padj))		# 1,115

res <- res[ !is.na(res$padj), ]			# 19,866
res_vr <- res_vr[ !is.na(res_vr$padj), ]	# 18,049

sig <- res[which(res$padj < 0.01), ]		# 11,369
sig_vr <- res_vr[which(res_vr$padj < 0.01), ]	# 11,917

sig_sorted <- sig[order(sig$log2FoldChange),]
sig_vr_sorted <- sig_vr[order(sig_vr$log2FoldChange),]

#####

# Annotation

transcriptome <- read.table("R05_Transcriptome.emapper.NoComments.annotations", sep="\t", row.names = 1, quote="", header=TRUE)

#####

# Write to table

if (!(file.exists("ExpSW_SigResults.tsv")))
	write.table(ExpSW_sig_sorted, file="ExpSW_SigResults.tsv", quote=FALSE, sep="\t")
if (!(file.exists("StatSW_SigResults.tsv")))
	write.table(StatSW_sig_sorted, file="StatSW_SigResults.tsv", quote=FALSE, sep="\t")
if (!(file.exists("SummerES_SigResults.tsv")))
	write.table(SummerES_sig_sorted, file="SummerES_SigResults.tsv", quote=FALSE, sep="\t")
if (!(file.exists("WinterES_SigResults.tsv")))
	write.table(WinterES_sig_sorted, file="WinterES_SigResults.tsv", quote=FALSE, sep="\t")

if (!(file.exists("SE_WS_SigResults.tsv")))
	write.table(SE_WS_sig_sorted, file="SE_WS_SigResults.tsv", quote=FALSE, sep="\t")
if (!(file.exists("WS_SE_SigResults.tsv")))
	write.table(WS_SE_sig_sorted, file="WS_SE_SigResults.tsv", quote=FALSE, sep="\t")


# ExpSW_sig[order(ExpSW_sig$log2FoldChange),]

ExpSW_SigResults.UpInWinter <- ExpSW_sig_sorted[which(ExpSW_sig_sorted$log2FoldChange > 0),]
ExpSW_SigResults.UpInWinter <- ExpSW_SigResults.UpInWinter[order(ExpSW_SigResults.UpInWinter$log2FoldChange, decreasing=TRUE),]
ExpSW_SigResults.UpInWinter.annot <- merge(data.frame(ExpSW_SigResults.UpInWinter), transcriptome, by=0, sort=FALSE)			# 3,597 annotations
if (!(file.exists("ExpSW_SigResults.UpInWinter.tsv")))
        write.table(ExpSW_SigResults.UpInWinter.annot, file="ExpSW_SigResults.UpInWinter.tsv", quote=FALSE, row.names=FALSE, sep="\t")

StatSW_SigResults.UpInWinter <- StatSW_sig_sorted[which(StatSW_sig_sorted$log2FoldChange > 0),]
StatSW_SigResults.UpInWinter <- StatSW_SigResults.UpInWinter[order(StatSW_SigResults.UpInWinter$log2FoldChange, decreasing=TRUE),]
StatSW_SigResults.UpInWinter.annot <- merge(data.frame(StatSW_SigResults.UpInWinter), transcriptome, by=0, sort=FALSE)			# 4,199 annotations
if (!(file.exists("StatSW_SigResults.UpInWinter.tsv")))
        write.table(StatSW_SigResults.UpInWinter.annot, file="StatSW_SigResults.UpInWinter.tsv", quote=FALSE, row.names=FALSE, sep="\t")

SummerES_SigResults.UpInStat <- SummerES_sig_sorted[which(SummerES_sig_sorted$log2FoldChange > 0),]
SummerES_SigResults.UpInStat <- SummerES_SigResults.UpInStat[order(SummerES_SigResults.UpInStat$log2FoldChange, decreasing=TRUE),]
SummerES_SigResults.UpInStat.annot <- merge(data.frame(SummerES_SigResults.UpInStat), transcriptome, by=0, sort=FALSE)			# 4,206 annotations
if (!(file.exists("SummerES_SigResults.UpInStat.tsv")))
        write.table(SummerES_SigResults.UpInStat.annot, file="SummerES_SigResults.UpInStat.tsv", quote=FALSE, row.names=FALSE, sep="\t")

WinterES_SigResults.UpInStat <- WinterES_sig_sorted[which(WinterES_sig_sorted$log2FoldChange > 0),]
WinterES_SigResults.UpInStat <- WinterES_SigResults.UpInStat[order(WinterES_SigResults.UpInStat$log2FoldChange, decreasing=TRUE),]
WinterES_SigResults.UpInStat.annot <- merge(data.frame(WinterES_SigResults.UpInStat), transcriptome, by=0, sort=FALSE)			# 1,929 annotations
if (!(file.exists("WinterES_SigResults.UpInStat.tsv")))
        write.table(WinterES_SigResults.UpInStat.annot, file="WinterES_SigResults.UpInStat.tsv", quote=FALSE, row.names=FALSE, sep="\t")

ExpSW_SigResults.UpInSummer <- ExpSW_sig_sorted[which(ExpSW_sig_sorted$log2FoldChange < 0),]
ExpSW_SigResults.UpInSummer.annot <- merge(data.frame(ExpSW_SigResults.UpInSummer), transcriptome, by=0, sort=FALSE)			# 3,109 annotations
if (!(file.exists("ExpSW_SigResults.UpInSummer.tsv")))
        write.table(ExpSW_SigResults.UpInSummer.annot, file="ExpSW_SigResults.UpInSummer.tsv", quote=FALSE, row.names=FALSE, sep="\t")

StatSW_SigResults.UpInSummer <- StatSW_sig_sorted[which(StatSW_sig_sorted$log2FoldChange < 0),]
StatSW_SigResults.UpInSummer.annot <- merge(data.frame(StatSW_SigResults.UpInSummer), transcriptome, by=0, sort=FALSE)			# 3,946 annotations
if (!(file.exists("StatSW_SigResults.UpInSummer.tsv")))
        write.table(StatSW_SigResults.UpInSummer.annot, file="StatSW_SigResults.UpInSummer.tsv", quote=FALSE, row.names=FALSE, sep="\t")

SummerES_SigResults.UpInExp <- SummerES_sig_sorted[which(SummerES_sig_sorted$log2FoldChange < 0),]
SummerES_SigResults.UpInExp.annot <- merge(data.frame(SummerES_SigResults.UpInExp), transcriptome, by=0, sort=FALSE)			# 3,512 annotations
if (!(file.exists("SummerES_SigResults.UpInExp.tsv")))
        write.table(SummerES_SigResults.UpInExp.annot, file="SummerES_SigResults.UpInExp.tsv", quote=FALSE, row.names=FALSE, sep="\t")

WinterES_SigResults.UpInExp <- WinterES_sig_sorted[which(WinterES_sig_sorted$log2FoldChange < 0),]
WinterES_SigResults.UpInExp.annot <- merge(data.frame(WinterES_SigResults.UpInExp), transcriptome, by=0, sort=FALSE)			# 2,092 annotations
if (!(file.exists("WinterES_SigResults.UpInExp.tsv")))
        write.table(WinterES_SigResults.UpInExp.annot, file="WinterES_SigResults.UpInExp.tsv", quote=FALSE, row.names=FALSE, sep="\t")

# Negative = down in WS
SE_WS_SigResults.UpInSE <- SE_WS_sig_sorted[which(SE_WS_sig_sorted$log2FoldChange < 0),]
SE_WS_SigResults.UpInSE.annot <- merge(data.frame(SE_WS_SigResults.UpInSE), transcriptome, by=0, sort=FALSE)
if (!(file.exists("SE_WS_SigResults.UpInSE.tsv")))
	write.table(SE_WS_SigResults.UpInSE.annot, file="SE_WS_SigResults.UpInSE.tsv", quote=FALSE, row.names=FALSE, sep="\t")

SE_WS_SigResults.UpInWS <- SE_WS_sig_sorted[which(SE_WS_sig_sorted$log2FoldChange > 0),]
SE_WS_SigResults.UpInWS <- SE_WS_SigResults.UpInWS[order(SE_WS_SigResults.UpInWS$log2FoldChange, decreasing=TRUE),]
SE_WS_SigResults.UpInWS.annot <- merge(data.frame(SE_WS_SigResults.UpInWS), transcriptome, by=0, sort=FALSE)
if (!(file.exists("SE_WS_SigResults.UpInWS.tsv")))
	write.table(SE_WS_SigResults.UpInWS.annot, file="SE_WS_SigResults.UpInWS.tsv", quote=FALSE, row.names=FALSE, sep="\t")

# Negative = down in SE
WS_SE_SigResults.UpInWS <- WS_SE_sig_sorted[which(WS_SE_sig_sorted$log2FoldChange < 0),]
WS_SE_SigResults.UpInWS.annot <- merge(data.frame(WS_SE_SigResults.UpInWS), transcriptome, by=0, sort=FALSE)
if (!(file.exists("WS_SE_SigResults.UpInWS.tsv")))
	write.table(WS_SE_SigResults.UpInWS.annot, file="WS_SE_SigResults.UpInWS.tsv", quote=FALSE, row.names=FALSE, sep="\t")

WS_SE_SigResults.UpInSE <- WS_SE_sig_sorted[which(WS_SE_sig_sorted$log2FoldChange > 0),]
WS_SE_SigResults.UpInSE <- WS_SE_SigResults.UpInSE[order(WS_SE_SigResults.UpInSE$log2FoldChange, decreasing=TRUE),]
WS_SE_SigResults.UpInSE.annot <- merge(data.frame(WS_SE_SigResults.UpInSE), transcriptome, by=0, sort=FALSE)
if (!(file.exists("WS_SE_SigResults.UpInSE.tsv")))
	write.table(WS_SE_SigResults.UpInSE.annot, file="WS_SE_SigResults.UpInSE.tsv", quote=FALSE, row.names=FALSE, sep="\t")

#####

# Note: 'Upregulation' is relative to the other group

length(which(ExpSW_sig_sorted$log2FoldChange > 0))	# 4,194 upregulated in winter
length(which(ExpSW_sig_sorted$log2FoldChange < 0))	# 4,058 upregulated in summer

length(which(StatSW_sig_sorted$log2FoldChange > 0))	# 5,054 upregulated in winter
length(which(StatSW_sig_sorted$log2FoldChange < 0))	# 5,005 upregulated in summer

length(which(SummerES_sig_sorted$log2FoldChange > 0))	# 4,988 upregulated in stationary phase
length(which(SummerES_sig_sorted$log2FoldChange < 0))	# 4,483 upregulated in exponential phase

length(which(WinterES_sig_sorted$log2FoldChange > 0))	# 2,335 upregulated in stationary phase
length(which(WinterES_sig_sorted$log2FoldChange < 0))	# 2,472 upregulated in exponential phase

length(which(SE_WS_sig_sorted$log2FoldChange > 0))	# 4,351 upregulated in winter
length(which(SE_WS_sig_sorted$log2FoldChange < 0))	# 4,062 upregulated in summer

length(which(WS_SE_sig_sorted$log2FoldChange > 0))	# 4,062 upregulated in summer
length(which(WS_SE_sig_sorted$log2FoldChange < 0))	# 4,351 upregulated in winter

#####

ExpSW_sig2heatmap <- assay(ExpSWrld)[rownames(assay(ExpSWrld))%in%rownames(ExpSW_sig),]
StatSW_sig2heatmap <- assay(StatSWrld)[rownames(assay(StatSWrld))%in%rownames(StatSW_sig),]
SummerES_sig2heatmap <- assay(SummerESrld)[rownames(assay(SummerESrld))%in%rownames(SummerES_sig),]
WinterES_sig2heatmap <- assay(WinterESrld)[rownames(assay(WinterESrld))%in%rownames(WinterES_sig),]

ExpSW_sig2heatmap_Mean <- ExpSW_sig2heatmap - rowMeans(ExpSW_sig2heatmap)
StatSW_sig2heatmap_Mean <- StatSW_sig2heatmap - rowMeans(StatSW_sig2heatmap)
SummerES_sig2heatmap_Mean <- SummerES_sig2heatmap - rowMeans(SummerES_sig2heatmap)
WinterES_sig2heatmap_Mean <- WinterES_sig2heatmap - rowMeans(WinterES_sig2heatmap)

SE_WS_sig2heatmap <- assay(SE_WSrld)[rownames(assay(SE_WSrld))%in%rownames(SE_WS_sig),]
SE_WS_sig2heatmap_Mean <- SE_WS_sig2heatmap - rowMeans(SE_WS_sig2heatmap)

WS_SE_sig2heatmap <- assay(WS_SErld)[rownames(assay(WS_SErld))%in%rownames(WS_SE_sig),]
WS_SE_sig2heatmap_Mean <- WS_SE_sig2heatmap - rowMeans(WS_SE_sig2heatmap)

if (!(file.exists("03_ExpSW_HeatMap.png"))) {
	png("03_ExpSW_HeatMap.png")
	pheatmap(ExpSW_sig2heatmap_Mean, main="Exponential - Summer vs Winter", annotation_col = ExpSWannotation_col, annotation_colors = annotation_colors, show_rownames = FALSE)
	dev.off()
}

if (!(file.exists("03_StatSW_HeatMap.png"))) {
	png("03_StatSW_HeatMap.png")
	pheatmap(StatSW_sig2heatmap_Mean, main="Stationary - Summer vs Winter", annotation_col = StatSWannotation_col, annotation_colors = annotation_colors, show_rownames = FALSE)
	dev.off()
}

if (!(file.exists("03_SummerES_HeatMap.png"))) {
	png("03_SummerES_HeatMap.png")
	pheatmap(SummerES_sig2heatmap_Mean, main="Summer - Exponential vs Stationary", annotation_col = SummerESannotation_col, annotation_colors = annotation_colors, show_rownames = FALSE)
	dev.off()
}

if (!(file.exists("03_WinterES_HeatMap.png"))) {
	png("03_WinterES_HeatMap.png")
	pheatmap(WinterES_sig2heatmap_Mean, main="Winter - Exponential vs Stationary", annotation_col = WinterESannotation_col, annotation_colors = annotation_colors, show_rownames = FALSE)
	dev.off()
}

if (!(file.exists("03_SE_WS_HeatMap.png"))) {
	png("03_SE_WS_HeatMap.png")
	pheatmap(SE_WS_sig2heatmap_Mean, main="Summer Exponential vs Winter Stationary", annotation_col = SE_WSannotation_col, annotation_colors = annotation_colors, show_rownames = FALSE)
	dev.off()
}

if (!(file.exists("03_WS_SE_HeatMap.png"))) {
	png("03_WS_SE_HeatMap.png")
	pheatmap(WS_SE_sig2heatmap_Mean, main="Winter Stationary vs Summer Exponential", annotation_col = WS_SEannotation_col, annotation_colors = annotation_colors, show_rownames = FALSE)
	dev.off()
}

#####

# plotCounts(???dds, gene="", intgroup = "season")
# plotCounts(???dds, gene="", intgroup = "phase")

# GOs
# EC
# KEGG_Pathway
