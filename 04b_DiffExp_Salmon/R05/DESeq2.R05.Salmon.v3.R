#!/usr/bin/env Rscript

library("Rsubread")
library("tximport")
library("readr")
library("DESeq2")
library("dplyr")
library("ggplot2")
library("pheatmap")
library("RColorBrewer")
library("apeglm")
library("genefilter")
library("topGO")

##################################################
# READ IN SAMPLE DATA
##################################################

coldata <- read.delim("R05_samples.tsv", row.names=1, stringsAsFactors=TRUE)
coldata$files <- file.path("../../03_quantification/VsGeneModels_v1.1.2", coldata$Sample_ID, "quant.sf")
file.exists(coldata$files)
rownames(coldata) <- coldata$Sample_ID
coldata$Dormancy <- as.factor(coldata$Dormancy)
coldata$Dormancy <- relevel(coldata$Dormancy, "0")

## For plotting PCA/heatmap without 0d samples
## Will also require regenerating the feature counts (saved as allcounts.ControlExcluded.rds)
# coldata <- coldata[coldata$Dormancy != 0,]
# coldata$Dormancy <- relevel(coldata$Dormancy, "49")

##################################################
# IMPORT TRANSCRIPT ABUNDANCE DATA
##################################################

txi <- tximport(coldata$files, type = "salmon", txOut = TRUE)

##################################################
# GENERATE DESeqDataSet OBJECT
##################################################

dds <- DESeqDataSetFromTximport(txi, coldata, ~Batch + Dormancy)

##################################################
# PRE-FILTER THE DATASET
##################################################

print(paste0("Number of genes: ", nrow(dds)))		# 17,383

## Slightly stricter filter - at least one count in at least two samples
## 189d has only two samples, so should retain more condition-specific genes
## Risk of more noise, however...

keep <- rowSums(counts(dds) >= 1 ) >= 2
dds <- dds[keep,]

print(paste0("Number of genes with one count in at least two samples: ", nrow(dds)))	# 14,377

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

if (!file.exists("R05.featureCounts.Heatmap.png")){
	png("R05.featureCounts.Heatmap.png")
	p <- pheatmap(distMatrix_rld, annotation_col = annotation_col)
	print(p)
	dev.off()
}

## PCA

if (!file.exists("R05.featureCounts.PCA.png")){
	png("R05.featureCounts.PCA.png")
	p <- plotPCA(rld, intgroup = c("Dormancy","Batch"))
	print(p)
	dev.off()
}

##################################################
# DIFFERENTIAL EXPRESSION ANALYSIS
#  ITERATE THROUGH COMBINATIONS OF TIMEPOINTS
##################################################

dds <- DESeq(dds)

transcriptome <- read.table("../v1.1.2.emapper.annotations", sep="\t", row.names = 1, quote="", header=TRUE)

geneID2GO <- readMappings("../R05_v1.1.2.GOterms.AllGenes.tsv")
geneNames <- names(geneID2GO)

timepoints <- as.list(levels(coldata$Dormancy))
combos <- combn(timepoints, 2)

for (pairnum in 1:ncol(combos)){
	numberone <- toString(combos[1,pairnum])
	numbertwo <- toString(combos[2,pairnum])

	resdir <- paste0(numberone, "_", numbertwo, "d")
	dir.create(resdir)
	dir.create(paste0(resdir,"/GO"))

#	res <- results(dds, alpha = 0.05, contrast=c("Dormancy", numberone, numbertwo))
	res <- results(dds, alpha = 0.05, contrast=c("Dormancy", numberone, numbertwo), lfcThreshold = 1)

	###

	res <- res[ !is.na(res$padj), ]
	print(paste0(numberone, "d vs ", numbertwo, "d - genes with adjusted p-value: ", nrow(res)))

	###

	sig <- res[which(res$padj < 0.05), ]
	print(paste0(numberone, "d vs ", numbertwo, "d - genes with adjusted p-value < 0.05: ", nrow(sig)))
	sig_sorted <- sig[order(sig$log2FoldChange),]

	#####

	sig_annot <- merge(data.frame(sig_sorted), transcriptome, by=0, sort=FALSE, all.x = TRUE)

	# Write to table

	if (nrow(sig) == 0){
		print(paste0("No significantly differentially expressed genes between ", numberone, "d and ", numbertwo, "d"))
	} else {
		outtable <- paste0(numberone, "_", numbertwo, "d_SigResults.tsv")
		if (!file.exists(outtable)){
			write.table(sig_annot, file=paste0(resdir,"/",outtable), quote=FALSE, row.names=FALSE, sep="\t")
		}
	}

	###

	# Negative LFC = upregulated in dormant

	sig_DormantUp <- sig_annot[which(sig_annot$log2FoldChange < 0),]
	if (nrow(sig_DormantUp) == 0){
		print(paste0("No significantly upregulated genes in ", numbertwo, "d vs. ", numberone, "d"))
	} else {
		print(paste0(numberone, "d vs ", numbertwo, "d - genes upregulated in ", numbertwo, "d: ",
			nrow(sig_DormantUp)))
		sig_DormantUp <- sig_DormantUp[order(sig_DormantUp$log2FoldChange, decreasing=FALSE),]
		print(paste0(numberone, "d vs ", numbertwo, "d - annotated genes upregulated in ", numbertwo, "d: ",
			nrow(sig_DormantUp[ !is.na(sig_DormantUp$seed_eggNOG_ortholog), ])))
		outtable <- paste0(numberone, "_", numbertwo, "d_SigResults.", numbertwo, "d_Up.tsv")
		if (!file.exists(outtable)){
			write.table(sig_DormantUp, file=paste0(resdir,"/",outtable), quote=FALSE, row.names=FALSE, sep="\t")
		}
	}

	sig_DormantDown <- sig_annot[which(sig_annot$log2FoldChange > 0),]
	if (nrow(sig_DormantDown) == 0){
		print(paste0("No significantly downregulated genes in ", numbertwo, "d vs. ", numberone, "d"))
	} else {
		print(paste0(numberone, "d vs ", numbertwo, "d - genes downregulated in ", numbertwo, "d: ",
			nrow(sig_DormantDown)))
		sig_DormantDown <- sig_DormantDown[order(sig_DormantDown$log2FoldChange, decreasing=TRUE),]
		print(paste0(numberone, "d vs ", numbertwo, "d - annotated genes downregulated in ", numbertwo, "d: ",
			nrow(sig_DormantDown[ !is.na(sig_DormantDown$seed_eggNOG_ortholog), ])))
		outtable <- paste0(numberone, "_", numbertwo, "d_SigResults.", numbertwo, "d_Down.tsv")
		if (!file.exists(outtable)){
			write.table(sig_DormantDown, file=paste0(resdir,"/",outtable), quote=FALSE, row.names=FALSE, sep="\t")
		}
	}

	### GO enrichment analysis

	myInterestingGenes_Up <- sig_DormantUp$Row.names
	geneList_Up <- factor(as.integer(geneNames %in% myInterestingGenes_Up))
	names(geneList_Up) <- geneNames
	if (!"1" %in% levels(geneList_Up)){
		print(paste0("No upregulated genes in ", numbertwo, "d vs. ", numberone))
	} else {
		for (GOgroup in c("BP", "MF", "CC")) {
			GOdata <- new("topGOdata", ontology = GOgroup, allGenes = geneList_Up, annot = annFUN.gene2GO, gene2GO = geneID2GO)
			resultFis <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
			allRes <- GenTable(GOdata, classic = resultFis, orderBy = "weight", ranksOf = "classic", topNodes = 1000)
			write.table(allRes, paste0(resdir, "/GO/", numberone, "_", numbertwo, "d_", GOgroup, ".UpIn", numbertwo, "d.Results.txt"),
				quote = FALSE, sep = "\t", row.names = FALSE)
			#showSigOfNodes(GOdata, score(resultFis), firstSigNodes = 5, useInfo = 'all')
			printGraph(GOdata, resultFis, firstSigNodes = 5, fn.prefix = paste0(resdir, "/GO/", numberone, "_", numbertwo, "d_", GOgroup, ".UpIn", numbertwo, "d"),
				useInfo = "all", pdfSW = TRUE)
		}
	}

	myInterestingGenes_Down <- sig_DormantDown$Row.names
	geneList_Down <- factor(as.integer(geneNames %in% myInterestingGenes_Down))
	names(geneList_Down) <- geneNames
	if (!"1" %in% levels(geneList_Down)){
		print(paste0("No downregulated genes in ", numbertwo, "d vs. ", numberone))
	} else {
		for (GOgroup in c("BP", "MF", "CC")) {
			GOdata <- new("topGOdata", ontology = GOgroup, allGenes = geneList_Down, annot = annFUN.gene2GO, gene2GO = geneID2GO)
			resultFis <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
			allRes <- GenTable(GOdata, classic = resultFis, orderBy = "weight", ranksOf = "classic", topNodes = 1000)
			write.table(allRes, paste0(resdir, "/GO/", numberone, "_", numbertwo, "d_", GOgroup, ".DownIn", numbertwo, "d.Results.txt"),
				quote = FALSE, sep = "\t", row.names = FALSE)
			#showSigOfNodes(GOdata, score(resultFis), firstSigNodes = 5, useInfo = 'all')
			printGraph(GOdata, resultFis, firstSigNodes = 5, fn.prefix = paste0(resdir, "/GO/", numberone, "_", numbertwo, "d_", GOgroup, ".DownIn", numbertwo, "d"),
				useInfo = "all", pdfSW = TRUE)
		}
	}
}

#####

# plotCounts(dds, gene="", intgroup = "Dormancy")

# GOs
# EC
# KEGG_Pathway
