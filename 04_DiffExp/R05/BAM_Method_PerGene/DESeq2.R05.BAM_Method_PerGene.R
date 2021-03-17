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

##################################################
# READ IN SALMON DATA
##################################################

coldata <- read.delim("../R05_samples.tsv", row.names=1, stringsAsFactors=TRUE)
coldata$files <- file.path("../../../03_quantification", coldata$Sample_ID, "quant.sf")
file.exists(coldata$files)
rownames(coldata) <- coldata$Sample_ID
coldata$Dormancy <- as.factor(coldata$Dormancy)
coldata$Dormancy <- relevel(coldata$Dormancy, "0")

##################################################
# IMPORT TRANSCRIPT ABUNDANCE DATA
##################################################

tx2gene <- read.csv("../R05_tx2gene.lst", sep="\t", header=FALSE)
txi <- tximport(coldata$files, type = "salmon", tx2gene = tx2gene)

##################################################
# GENERATE DESeqDataSet OBJECT
##################################################

dds <- DESeqDataSetFromTximport(txi, coldata, ~Batch + Dormancy)

##################################################
# PRE-FILTER THE DATASET
##################################################

print(paste0("Number of genes: ", nrow(dds)))		# 38,865

## Slightly stricter filter - at least one count in at least three samples
## Three samples per condition, so should retain condition-specific genes

keep <- rowSums(counts(dds) >= 1 ) >= 3
dds <- dds[keep,]

print(paste0("Number of genes with one count in at least three samples: ", nrow(dds)))	# 28,702

##################################################
# DATA TRANSFORMATION
# (FOR APPLICATIONS OTHER THAN DIFF. TESTING)
##################################################

## rlog transformation - recommended for n < 30
rld <- rlog(dds, blind = FALSE)

##################################################
# ASSESS THE SIMILARITY BETWEEN SAMPLES
##################################################

## Heat map

annotation_col <- coldata[,c(4,5)]

all_dist_rld <- dist(t(assay(rld)))

distMatrix_rld <- as.matrix(all_dist_rld)

#png("R05.BAM_Method_PerGene.Heatmap.png")
#pheatmap(distMatrix_rld, annotation_col = annotation_col)
#dev.off()

## PCA

#png("R05.BAM_Method_PerGene.PCA.png")
#plotPCA(rld, intgroup = c("Dormancy","Batch"))
#dev.off()

##################################################
# DIFFERENTIAL EXPRESSION ANALYSIS
##################################################

dds <- DESeq(dds)

res49 <- results(dds, alpha = 0.05, contrast=c("Dormancy","0","49"))
res56 <- results(dds, alpha = 0.05, contrast=c("Dormancy","0","56"))
res72 <- results(dds, alpha = 0.05, contrast=c("Dormancy","0","72"))
res91 <- results(dds, alpha = 0.05, contrast=c("Dormancy","0","91"))
res126 <- results(dds, alpha = 0.05, contrast=c("Dormancy","0","126"))
res189 <- results(dds, alpha = 0.05, contrast=c("Dormancy","0","189"))

res49_56 <- results(dds, alpha = 0.05, contrast=c("Dormancy","49","56"))
res49_72 <- results(dds, alpha = 0.05, contrast=c("Dormancy","49","72"))
res49_91 <- results(dds, alpha = 0.05, contrast=c("Dormancy","49","91"))
res49_126 <- results(dds, alpha = 0.05, contrast=c("Dormancy","49","126"))
res49_189 <- results(dds, alpha = 0.05, contrast=c("Dormancy","49","189"))

res56_72 <- results(dds, alpha = 0.05, contrast=c("Dormancy","56","72"))
res56_91 <- results(dds, alpha = 0.05, contrast=c("Dormancy","56","91"))
res56_126 <- results(dds, alpha = 0.05, contrast=c("Dormancy","56","126"))
res56_189 <- results(dds, alpha = 0.05, contrast=c("Dormancy","56","189"))

res72_91 <- results(dds, alpha = 0.05, contrast=c("Dormancy","72","91"))
res72_126 <- results(dds, alpha = 0.05, contrast=c("Dormancy","72","126"))
res72_189 <- results(dds, alpha = 0.05, contrast=c("Dormancy","72","189"))

res91_126 <- results(dds, alpha = 0.05, contrast=c("Dormancy","91","126"))
res91_189 <- results(dds, alpha = 0.05, contrast=c("Dormancy","91","189"))

res126_189 <- results(dds, alpha = 0.05, contrast=c("Dormancy","126","189"))

###

res49 <- res49[ !is.na(res49$padj), ]
print(paste0("0d vs 49d - genes with adjusted p-value: ", nrow(res49)))         # 21,990
res56 <- res56[ !is.na(res56$padj), ]
print(paste0("0d vs 56d - genes with adjusted p-value: ", nrow(res56)))         # 23,103
res72 <- res72[ !is.na(res72$padj), ]
print(paste0("0d vs 72d - genes with adjusted p-value: ", nrow(res72)))         # 23,659
res91 <- res91[ !is.na(res91$padj), ]
print(paste0("0d vs 91d - genes with adjusted p-value: ", nrow(res91)))         # 23,659
res126 <- res126[ !is.na(res126$padj), ]
print(paste0("0d vs 126d - genes with adjusted p-value: ", nrow(res126)))       # 24,772
res189 <- res189[ !is.na(res189$padj), ]
print(paste0("0d vs 189d - genes with adjusted p-value: ", nrow(res189)))       # 22,547

res49_56 <- res49_56[ !is.na(res49_56$padj), ]
print(paste0("49d vs 56d - genes with adjusted p-value: ", nrow(res49_56)))	## 19,211
res49_72 <- res49_72[ !is.na(res49_72$padj), ]
print(paste0("49d vs 72d - genes with adjusted p-value: ", nrow(res49_72)))	## 19,211
res49_91 <- res49_91[ !is.na(res49_91$padj), ]
print(paste0("49d vs 91d - genes with adjusted p-value: ", nrow(res49_91)))	## 17,557
res49_126 <- res49_126[ !is.na(res49_126$padj), ]
print(paste0("49d vs 126d - genes with adjusted p-value: ", nrow(res49_126)))	## 18,657
res49_189 <- res49_189[ !is.na(res49_189$padj), ]
print(paste0("49d vs 189d - genes with adjusted p-value: ", nrow(res49_189)))	## 21,990

res56_72 <- res56_72[ !is.na(res56_72$padj), ]
print(paste0("56d vs 72d - genes with adjusted p-value: ", nrow(res56_72)))	## 16,452
res56_91 <- res56_91[ !is.na(res56_91$padj), ]
print(paste0("56d vs 91d - genes with adjusted p-value: ", nrow(res56_91)))	## 15,897
res56_126 <- res56_126[ !is.na(res56_126$padj), ]
print(paste0("56d vs 126d - genes with adjusted p-value: ", nrow(res56_126)))	## 19,764
res56_189 <- res56_189[ !is.na(res56_189$padj), ]
print(paste0("56d vs 189d - genes with adjusted p-value: ", nrow(res56_189)))	## 22,547

res72_91 <- res72_91[ !is.na(res72_91$padj), ]
print(paste0("72d vs 91d - genes with adjusted p-value: ", nrow(res72_91)))	## 14,786
res72_126 <- res72_126[ !is.na(res72_126$padj), ]
print(paste0("72d vs 126d - genes with adjusted p-value: ", nrow(res72_126)))	## 18,110
res72_189 <- res72_189[ !is.na(res72_189$padj), ]
print(paste0("72d vs 189d - genes with adjusted p-value: ", nrow(res72_189)))	## 21,434

res91_126 <- res91_126[ !is.na(res91_126$padj), ]
print(paste0("91d vs 126d - genes with adjusted p-value: ", nrow(res91_126)))	## 17,004
res91_189 <- res91_189[ !is.na(res91_189$padj), ]
print(paste0("91d vs 189d - genes with adjusted p-value: ", nrow(res91_189)))	## 19,764

res126_189 <- res126_189[ !is.na(res126_189$padj), ]
print(paste0("126d vs 189d - genes with adjusted p-value: ", nrow(res126_189)))	## 18,657

###

sig49 <- res49[which(res49$padj < 0.05), ]
print(paste0("0d vs 49d - genes with adjusted p-value < 0.05: ", nrow(sig49)))          # 14,507
sig49_sorted <- sig49[order(sig49$log2FoldChange),]
sig56 <- res56[which(res56$padj < 0.05), ]
print(paste0("0d vs 56d - genes with adjusted p-value < 0.05: ", nrow(sig56)))          # 14,316
sig56_sorted <- sig56[order(sig56$log2FoldChange),]
sig72 <- res72[which(res72$padj < 0.05), ]
print(paste0("0d vs 72d - genes with adjusted p-value < 0.05: ", nrow(sig72)))          # 14,893
sig72_sorted <- sig72[order(sig72$log2FoldChange),]
sig91 <- res91[which(res91$padj < 0.05), ]
print(paste0("0d vs 91d - genes with adjusted p-value < 0.05: ", nrow(sig91)))          # 15,254
sig91_sorted <- sig91[order(sig91$log2FoldChange),]
sig126 <- res126[which(res126$padj < 0.05), ]
print(paste0("0d vs 126d - genes with adjusted p-value < 0.05: ", nrow(sig126)))        # 15,797
sig126_sorted <- sig126[order(sig126$log2FoldChange),]
sig189 <- res189[which(res189$padj < 0.05), ]
print(paste0("0d vs 189d - genes with adjusted p-value < 0.05: ", nrow(sig189)))        # 13,859
sig189_sorted <- sig189[order(sig189$log2FoldChange),]

sig49_56 <- res49_56[which(res49_56$padj < 0.05), ]
print(paste0("49d vs 56d - genes with adjusted p-value < 0.05: ", nrow(sig49_56)))	## 6,770
sig49_56_sorted <- sig49_56[order(sig49_56$log2FoldChange),]
sig49_72 <- res49_72[which(res49_72$padj < 0.05), ]
print(paste0("49d vs 72d - genes with adjusted p-value < 0.05: ", nrow(sig49_72)))	## 6,682
sig49_72_sorted <- sig49_72[order(sig49_72$log2FoldChange),]
sig49_91 <- res49_91[which(res49_91$padj < 0.05), ]
print(paste0("49d vs 91d - genes with adjusted p-value < 0.05: ", nrow(sig49_91)))	## 4,191
sig49_91_sorted <- sig49_91[order(sig49_91$log2FoldChange),]
sig49_126 <- res49_126[which(res49_126$padj < 0.05), ]
print(paste0("49d vs 126d - genes with adjusted p-value < 0.05: ", nrow(sig49_126)))	## 5,927
sig49_126_sorted <- sig49_126[order(sig49_126$log2FoldChange),]
sig49_189 <- res49_189[which(res49_189$padj < 0.05), ]
print(paste0("49d vs 189d - genes with adjusted p-value < 0.05: ", nrow(sig49_189)))	## 9.936
sig49_189_sorted <- sig49_189[order(sig49_189$log2FoldChange),]

sig56_72 <- res56_72[which(res56_72$padj < 0.05), ]
print(paste0("56d vs 72d - genes with adjusted p-value < 0.05: ", nrow(sig56_72)))	## 3,264
sig56_72_sorted <- sig56_72[order(sig56_72$log2FoldChange),]
sig56_91 <- res56_91[which(res56_91$padj < 0.05), ]
print(paste0("56d vs 91d - genes with adjusted p-value < 0.05: ", nrow(sig56_91)))	## 3,412
sig56_91_sorted <- sig56_91[order(sig56_91$log2FoldChange),]
sig56_126 <- res56_126[which(res56_126$padj < 0.05), ]
print(paste0("56d vs 126d - genes with adjusted p-value < 0.05: ", nrow(sig56_126)))	## 7,760
sig56_126_sorted <- sig56_126[order(sig56_126$log2FoldChange),]
sig56_189 <- res56_189[which(res56_189$padj < 0.05), ]
print(paste0("56d vs 189d - genes with adjusted p-value < 0.05: ", nrow(sig56_189)))	## 9,828
sig56_189_sorted <- sig56_189[order(sig56_189$log2FoldChange),]

sig72_91 <- res72_91[which(res72_91$padj < 0.05), ]
print(paste0("72d vs 91d - genes with adjusted p-value < 0.05: ", nrow(sig72_91)))	## 334
sig72_91_sorted <- sig72_91[order(sig72_91$log2FoldChange),]
sig72_126 <- res72_126[which(res72_126$padj < 0.05), ]
print(paste0("72d vs 126d - genes with adjusted p-value < 0.05: ", nrow(sig72_126)))	## 5,123
sig72_126_sorted <- sig72_126[order(sig72_126$log2FoldChange),]
sig72_189 <- res72_189[which(res72_189$padj < 0.05), ]
print(paste0("72d vs 189d - genes with adjusted p-value < 0.05: ", nrow(sig72_189)))	## 8,711
sig72_189_sorted <- sig72_189[order(sig72_189$log2FoldChange),]

sig91_126 <- res91_126[which(res91_126$padj < 0.05), ]
print(paste0("91d vs 126d - genes with adjusted p-value < 0.05: ", nrow(sig91_126)))	## 3,068
sig91_126_sorted <- sig91_126[order(sig91_126$log2FoldChange),]
sig91_189 <- res91_189[which(res91_189$padj < 0.05), ]
print(paste0("91d vs 189d - genes with adjusted p-value < 0.05: ", nrow(sig91_189)))	## 7,158
sig91_189_sorted <- sig91_189[order(sig91_189$log2FoldChange),]

sig126_189 <- res126_189[which(res126_189$padj < 0.05), ]
print(paste0("126d vs 189d - genes with adjusted p-value < 0.05: ", nrow(sig126_189)))	## 4,773
sig126_189_sorted <- sig126_189[order(sig126_189$log2FoldChange),]

#####

# Annotation
# Need to figure out gene-level annotation; currently at transcript-level

# transcriptome <- read.table("../R05.emapper.annotations", sep="\t", row.names = 1, quote="", header=TRUE)

#####

# Write to table

#write.table(sig49_sorted, file="0_49d_SigResults.tsv", quote=FALSE, sep="\t")
#write.table(sig56_sorted, file="0_56d_SigResults.tsv", quote=FALSE, sep="\t")
#write.table(sig72_sorted, file="0_72d_SigResults.tsv", quote=FALSE, sep="\t")
#write.table(sig91_sorted, file="0_91d_SigResults.tsv", quote=FALSE, sep="\t")
#write.table(sig126_sorted, file="0_126d_SigResults.tsv", quote=FALSE, sep="\t")
#write.table(sig189_sorted, file="0_189d_SigResults.tsv", quote=FALSE, sep="\t")

#write.table(sig49_56_sorted, file="49_56d_SigResults.tsv", quote=FALSE, sep="\t")
#write.table(sig49_72_sorted, file="49_72d_SigResults.tsv", quote=FALSE, sep="\t")
#write.table(sig49_91_sorted, file="49_91d_SigResults.tsv", quote=FALSE, sep="\t")
#write.table(sig49_126_sorted, file="49_126d_SigResults.tsv", quote=FALSE, sep="\t")
#write.table(sig49_189_sorted, file="49_189d_SigResults.tsv", quote=FALSE, sep="\t")

#write.table(sig56_72_sorted, file="56_72d_SigResults.tsv", quote=FALSE, sep="\t")
#write.table(sig56_91_sorted, file="56_91d_SigResults.tsv", quote=FALSE, sep="\t")
#write.table(sig56_126_sorted, file="56_126d_SigResults.tsv", quote=FALSE, sep="\t")
#write.table(sig56_189_sorted, file="56_189d_SigResults.tsv", quote=FALSE, sep="\t")

#write.table(sig72_91_sorted, file="72_91d_SigResults.tsv", quote=FALSE, sep="\t")
#write.table(sig72_126_sorted, file="72_126d_SigResults.tsv", quote=FALSE, sep="\t")
#write.table(sig72_189_sorted, file="72_189d_SigResults.tsv", quote=FALSE, sep="\t")

#write.table(sig91_126_sorted, file="91_126d_SigResults.tsv", quote=FALSE, sep="\t")
#write.table(sig91_189_sorted, file="91_189d_SigResults.tsv", quote=FALSE, sep="\t")

#write.table(sig126_189_sorted, file="126_189d_SigResults.tsv", quote=FALSE, sep="\t")

###

# Negative LFC = upregulated in dormant

sig49_DormantUp <- sig49_sorted[which(sig49_sorted$log2FoldChange < 0),]
print(paste0("0d vs 49d - genes upregulated in dormancy: ", nrow(sig49_DormantUp)))		# 6,997
sig49_DormantUp <- sig49_DormantUp[order(sig49_DormantUp$log2FoldChange, decreasing=FALSE),]
##sig49_DormantUp.annot <- merge(data.frame(sig49_DormantUp), transcriptome, by=0, sort=FALSE, all.x = TRUE)
##print(paste0("0d vs 49d - annotated genes upregulated in dormancy: ", nrow(sig49_DormantUp.annot[ !is.na(sig49_DormantUp.annot$seed_eggNOG_ortholog), ])))		# ??? annotated
##write.table(sig49_DormantUp.annot, file="0_49d_SigResults.DormantUp.tsv", quote=FALSE, row.names=TRUE, sep="\t")
#write.table(sig49_DormantUp, file="0_49d_SigResults.DormantUp.tsv", quote=FALSE, row.names=TRUE, sep="\t")

sig49_DormantDown <- sig49_sorted[which(sig49_sorted$log2FoldChange > 0),]
print(paste0("0d vs 49d - genes downregulated in dormancy: ", nrow(sig49_DormantDown)))		# 7,510
sig49_DormantDown <- sig49_DormantDown[order(sig49_DormantDown$log2FoldChange, decreasing=TRUE),]
##sig49_DormantDown.annot <- merge(data.frame(sig49_DormantDown), transcriptome, by=0, sort=FALSE, all.x = TRUE)
##print(paste0("0d vs 49d - annotated genes downregulated in dormancy: ", nrow(sig49_DormantDown.annot[ !is.na(sig49_DormantDown.annot$seed_eggNOG_ortholog), ])))	# ??? annotated
##write.table(sig49_DormantDown.annot, file="0_49d_SigResults.DormantDown.tsv", quote=FALSE, row.names=TRUE, sep="\t")
#write.table(sig49_DormantDown, file="0_49d_SigResults.DormantDown.tsv", quote=FALSE, row.names=TRUE, sep="\t")

sig56_DormantUp <- sig56_sorted[which(sig56_sorted$log2FoldChange < 0),]
print(paste0("0d vs 56d - genes upregulated in dormancy: ", nrow(sig56_DormantUp)))		# 6,760
sig56_DormantUp <- sig56_DormantUp[order(sig56_DormantUp$log2FoldChange, decreasing=FALSE),]
##sig56_DormantUp.annot <- merge(data.frame(sig56_DormantUp), transcriptome, by=0, sort=FALSE, all.x = TRUE)
##print(paste0("0d vs 56d - annotated genes upregulated in dormancy: ", nrow(sig56_DormantUp.annot[ !is.na(sig56_DormantUp.annot$seed_eggNOG_ortholog), ])))		# ??? annotated
##write.table(sig56_DormantUp.annot, file="0_56d_SigResults.DormantUp.tsv", quote=FALSE, row.names=TRUE, sep="\t")
#write.table(sig56_DormantUp, file="0_56d_SigResults.DormantUp.tsv", quote=FALSE, row.names=TRUE, sep="\t")

sig56_DormantDown <- sig56_sorted[which(sig56_sorted$log2FoldChange > 0),]
print(paste0("0d vs 56d - genes downregulated in dormancy: ", nrow(sig56_DormantDown)))		# 7,556
sig56_DormantDown <- sig56_DormantDown[order(sig56_DormantDown$log2FoldChange, decreasing=TRUE),]
##sig56_DormantDown.annot <- merge(data.frame(sig56_DormantDown), transcriptome, by=0, sort=FALSE, all.x = TRUE)
##print(paste0("0d vs 56d - annotated genes downregulated in dormancy: ", nrow(sig56_DormantDown.annot[ !is.na(sig56_DormantDown.annot$seed_eggNOG_ortholog), ])))	# ??? annotated
##write.table(sig56_DormantDown.annot, file="0_56d_SigResults.DormantDown.tsv", quote=FALSE, row.names=TRUE, sep="\t")
#write.table(sig56_DormantDown, file="0_56d_SigResults.DormantDown.tsv", quote=FALSE, row.names=TRUE, sep="\t")

sig72_DormantUp <- sig72_sorted[which(sig72_sorted$log2FoldChange < 0),]
print(paste0("0d vs 72d - genes upregulated in dormancy: ", nrow(sig72_DormantUp)))		# 7,143
sig72_DormantUp <- sig72_DormantUp[order(sig72_DormantUp$log2FoldChange, decreasing=FALSE),]
##sig72_DormantUp.annot <- merge(data.frame(sig72_DormantUp), transcriptome, by=0, sort=FALSE, all.x = TRUE)
##print(paste0("0d vs 72d - annotated genes upregulated in dormancy: ", nrow(sig72_DormantUp.annot[ !is.na(sig72_DormantUp.annot$seed_eggNOG_ortholog), ])))		# ??? annotated
##write.table(sig72_DormantUp.annot, file="0_72d_SigResults.DormantUp.tsv", quote=FALSE, row.names=TRUE, sep="\t")
#write.table(sig72_DormantUp, file="0_72d_SigResults.DormantUp.tsv", quote=FALSE, row.names=TRUE, sep="\t")

sig72_DormantDown <- sig72_sorted[which(sig72_sorted$log2FoldChange > 0),]
print(paste0("0d vs 72d - genes downregulated in dormancy: ", nrow(sig72_DormantDown)))		# 7,750
sig72_DormantDown <- sig72_DormantDown[order(sig72_DormantDown$log2FoldChange, decreasing=TRUE),]
##sig72_DormantDown.annot <- merge(data.frame(sig72_DormantDown), transcriptome, by=0, sort=FALSE, all.x = TRUE)
##print(paste0("0d vs 72d - annotated genes downregulated in dormancy: ", nrow(sig72_DormantDown.annot[ !is.na(sig72_DormantDown.annot$seed_eggNOG_ortholog), ])))	# ??? annotated
##write.table(sig72_DormantDown.annot, file="0_72d_SigResults.DormantDown.tsv", quote=FALSE, row.names=TRUE, sep="\t")
#write.table(sig72_DormantDown, file="0_72d_SigResults.DormantDown.tsv", quote=FALSE, row.names=TRUE, sep="\t")

sig91_DormantUp <- sig91_sorted[which(sig91_sorted$log2FoldChange < 0),]
print(paste0("0d vs 91d - genes upregulated in dormancy: ", nrow(sig91_DormantUp)))		# 7,251
sig91_DormantUp <- sig91_DormantUp[order(sig91_DormantUp$log2FoldChange, decreasing=FALSE),]
##sig91_DormantUp.annot <- merge(data.frame(sig91_DormantUp), transcriptome, by=0, sort=FALSE, all.x = TRUE)
##print(paste0("0d vs 91d - annotated genes upregulated in dormancy: ", nrow(sig91_DormantUp.annot[ !is.na(sig91_DormantUp.annot$seed_eggNOG_ortholog), ])))		# ??? annotated
##write.table(sig91_DormantUp.annot, file="0_91d_SigResults.DormantUp.tsv", quote=FALSE, row.names=TRUE, sep="\t")
#write.table(sig91_DormantUp, file="0_91d_SigResults.DormantUp.tsv", quote=FALSE, row.names=TRUE, sep="\t")

sig91_DormantDown <- sig91_sorted[which(sig91_sorted$log2FoldChange > 0),]
print(paste0("0d vs 91d - genes downregulated in dormancy: ", nrow(sig91_DormantDown)))		# 8,003
sig91_DormantDown <- sig91_DormantDown[order(sig91_DormantDown$log2FoldChange, decreasing=TRUE),]
##sig91_DormantDown.annot <- merge(data.frame(sig91_DormantDown), transcriptome, by=0, sort=FALSE, all.x = TRUE)	# ??? annotated
##print(paste0("0d vs 91d - annotated genes downregulated in dormancy: ", nrow(sig91_DormantDown.annot[ !is.na(sig91_DormantDown.annot$seed_eggNOG_ortholog), ])))	# ??? annotated
##write.table(sig91_DormantDown.annot, file="0_91d_SigResults.DormantDown.tsv", quote=FALSE, row.names=TRUE, sep="\t")
#write.table(sig91_DormantDown, file="0_91d_SigResults.DormantDown.tsv", quote=FALSE, row.names=TRUE, sep="\t")

sig126_DormantUp <- sig126_sorted[which(sig126_sorted$log2FoldChange < 0),]
print(paste0("0d vs 126d - genes upregulated in dormancy: ", nrow(sig126_DormantUp)))		# 7,605
sig126_DormantUp <- sig126_DormantUp[order(sig126_DormantUp$log2FoldChange, decreasing=FALSE),]
##sig126_DormantUp.annot <- merge(data.frame(sig126_DormantUp), transcriptome, by=0, sort=FALSE, all.x = TRUE)
##print(paste0("0d vs 126d - annotated genes upregulated in dormancy: ", nrow(sig126_DormantUp.annot[ !is.na(sig126_DormantUp.annot$seed_eggNOG_ortholog), ])))		# ??? annotated
##write.table(sig126_DormantUp.annot, file="0_126d_SigResults.DormantUp.tsv", quote=FALSE, row.names=TRUE, sep="\t")
#write.table(sig126_DormantUp, file="0_126d_SigResults.DormantUp.tsv", quote=FALSE, row.names=TRUE, sep="\t")

sig126_DormantDown <- sig126_sorted[which(sig126_sorted$log2FoldChange > 0),]
print(paste0("0d vs 126d - genes downregulated in dormancy: ", nrow(sig126_DormantDown)))	# 8,192
sig126_DormantDown <- sig126_DormantDown[order(sig126_DormantDown$log2FoldChange, decreasing=TRUE),]
##sig126_DormantDown.annot <- merge(data.frame(sig126_DormantDown), transcriptome, by=0, sort=FALSE, all.x = TRUE)
##print(paste0("0d vs 126d - annotated genes downregulated in dormancy: ", nrow(sig126_DormantDown.annot[ !is.na(sig126_DormantDown.annot$seed_eggNOG_ortholog), ])))	# ??? annotated
##write.table(sig126_DormantDown.annot, file="0_126d_SigResults.DormantDown.tsv", quote=FALSE, row.names=TRUE, sep="\t")
#write.table(sig126_DormantDown, file="0_126d_SigResults.DormantDown.tsv", quote=FALSE, row.names=TRUE, sep="\t")

sig189_DormantUp <- sig189_sorted[which(sig189_sorted$log2FoldChange < 0),]
print(paste0("0d vs 189d - genes upregulated in dormancy: ", nrow(sig189_DormantUp)))		# 7,112
sig189_DormantUp <- sig189_DormantUp[order(sig189_DormantUp$log2FoldChange, decreasing=FALSE),]
##sig189_DormantUp.annot <- merge(data.frame(sig189_DormantUp), transcriptome, by=0, sort=FALSE, all.x = TRUE)
##print(paste0("0d vs 189d - annotated genes upregulated in dormancy: ", nrow(sig189_DormantUp.annot[ !is.na(sig189_DormantUp.annot$seed_eggNOG_ortholog), ])))		# ??? annotated
##write.table(sig189_DormantUp.annot, file="0_189d_SigResults.DormantUp.tsv", quote=FALSE, row.names=TRUE, sep="\t")
#write.table(sig189_DormantUp, file="0_189d_SigResults.DormantUp.tsv", quote=FALSE, row.names=TRUE, sep="\t")

sig189_DormantDown <- sig189_sorted[which(sig189_sorted$log2FoldChange > 0),]
print(paste0("0d vs 189d - genes downregulated in dormancy: ", nrow(sig189_DormantDown)))	# 6,747
sig189_DormantDown <- sig189_DormantDown[order(sig189_DormantDown$log2FoldChange, decreasing=TRUE),]
##sig189_DormantDown.annot <- merge(data.frame(sig189_DormantDown), transcriptome, by=0, sort=FALSE, all.x = TRUE)
##print(paste0("0d vs 189d - annotated genes downregulated in dormancy: ", nrow(sig189_DormantDown.annot[ !is.na(sig189_DormantDown.annot$seed_eggNOG_ortholog), ])))	# ??? annotated
##write.table(sig189_DormantDown.annot, file="0_189d_SigResults.DormantDown.tsv", quote=FALSE, row.names=TRUE, sep="\t")
#write.table(sig189_DormantDown, file="0_189d_SigResults.DormantDown.tsv", quote=FALSE, row.names=TRUE, sep="\t")

#####

# plotCounts(dds, gene="", intgroup = "Dormancy")

# GOs
# EC
# KEGG_Pathway
