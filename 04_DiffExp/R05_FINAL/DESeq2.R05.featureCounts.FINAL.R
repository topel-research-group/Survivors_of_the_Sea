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
coldata$files <- file.path("../../01_mapping/VsReferenceGenome/StricterMapping",
			coldata$Sample_ID,
			paste(coldata$Sample_ID,"_vs_RefGenome.MaxInt10K.NoMixedDiscordant_sorted.bam", sep=""))
file.exists(coldata$files)
rownames(coldata) <- coldata$Sample_ID
coldata$Dormancy <- as.factor(coldata$Dormancy)
coldata$Dormancy <- relevel(coldata$Dormancy, "0")

## For plotting PCA/heatmap without 0d samples
## Will also require regenerating the feature counts (saved as allcounts.ControlExcluded.rds)
# coldata <- coldata[coldata$Dormancy != 0,]
# coldata$Dormancy <- relevel(coldata$Dormancy, "49")

##################################################
# QUANTIFY TRANSCRIPT ABUNDANCE
##################################################

if (file.exists("allcounts.rds")){
	allcounts <- readRDS("allcounts.rds")
} else {
	inannot <- "Sm_ManualCuration.v1.1.9.gff"
	infasta <- "../../01_mapping/References/Skeletonema_marinoi_Ref_v1.1.2.fst"

	allcounts <- featureCounts(coldata$files,
			annot.ext=inannot,
			isGTFAnnotationFile=TRUE,
			GTF.featureType="mRNA",
			GTF.attrType="ID",
			useMetaFeatures=TRUE,
			allowMultiOverlap=FALSE,
			largestOverlap=TRUE,
			countMultiMappingReads=TRUE,
			genome=infasta,
			isPairedEnd=TRUE,
			requireBothEndsMapped=TRUE,
			nthreads=10,
			verbose=TRUE)

	saveRDS(allcounts, file="allcounts.rds")
}

countdata <- allcounts$counts

##################################################
# DIFFERENTIAL EXPRESSION ANALYSIS
##################################################

if (file.exists("dds.rds")){
	dds <- readRDS("dds.rds")
} else {

	##################################################
	# GENERATE DESeqDataSet OBJECT
	##################################################

	dds <- DESeqDataSetFromMatrix(countdata, coldata, ~Batch + Dormancy)

	##################################################
	# PRE-FILTER THE DATASET
	##################################################

	print(paste0("Number of genes: ", nrow(dds)))		# 

	## Slightly stricter filter - at least one count in at least two samples
	## 189d has only two samples, so should retain more condition-specific genes
	## Risk of more noise, however...

	keep <- rowSums(counts(dds) >= 1 ) >= 2
	dds <- dds[keep,]

	print(paste0("Number of genes with one count in at least two samples: ", nrow(dds)))	# 

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

	if (!file.exists("R05.featureCounts.Heatmap.svg")){
		svg("R05.featureCounts.Heatmap.svg")
		p <- pheatmap(distMatrix_rld, annotation_col = annotation_col)
		print(p)
		dev.off()
	}

	## PCA

	if (!file.exists("R05.featureCounts.PCA.svg")){
		svg("R05.featureCounts.PCA.svg")
		p <- plotPCA(rld, intgroup = c("Dormancy","Batch"))
		print(p)
		dev.off()
	}

	## PCA zoomed in on the resting stages
	if (!file.exists("R05.featureCounts.PCA.RestingZoom.svg")){
		svg("R05.featureCounts.PCA.RestingZoom.svg")
		p <- plotPCA(rld, intgroup = c("Dormancy","Batch"))
		p <- p + coord_fixed(xlim=c(0,NA))
		print(p)
		dev.off()
	}

	##################################################
	# DIFFERENTIAL EXPRESSION ANALYSIS
	#  ITERATE THROUGH COMBINATIONS OF TIMEPOINTS
	##################################################

	dds <- DESeq(dds)

	if (!file.exists("dds.rds")){
		saveRDS(dds, file="dds.rds")
	}
}

geneID2GO <- readMappings("../R05_v1.1.9.GOterms.AllGenes.tsv")
geneNames <- names(geneID2GO)

timepoints <- as.list(levels(coldata$Dormancy))
combos <- combn(timepoints, 2)

for (threshold in c(0, 0.585, 1)){

	# Initialise subdirectory and results summary file for this LFC threshold
	resdir <- paste0("lfc_", threshold)
	dir.create(resdir)

	summaryfile <- paste0(resdir, "/ResultSummary_LFC_", threshold, ".tsv")
	# Addition of headers could likely be cleaner
	headers <- "Comparison\tGenes w/ adj. p-value\tGenes w/ adj. p-value < 0.05\tGenes upregulated in older sample\tGenes downregulated in older sample"
	write(headers, summaryfile)

	# Iterate through pairs of timepoints
	for (pairnum in 1:ncol(combos)){
		numberone <- toString(combos[1,pairnum])
		numbertwo <- toString(combos[2,pairnum])

		compdir <- paste0(resdir, "/", numberone, "_", numbertwo, "d")
		godir <- paste0(resdir, "/GO")
		dir.create(compdir)
		dir.create(godir)

		# Initiate empty list for summary results (print to summary file later)
		sumres <- list()
		sumres <- append(sumres, paste0(numberone, "d vs ", numbertwo, "d"))

		res <- results(dds, alpha = 0.05, contrast=c("Dormancy", numbertwo, numberone), lfcThreshold = threshold)

		###

		res <- res[ !is.na(res$padj), ]
		print(paste0(numberone, "d vs ", numbertwo, "d - genes with adjusted p-value: ", nrow(res)))
		sumres <- append(sumres, nrow(res))		

		###

		sig <- res[which(res$padj < 0.05), ]
		print(paste0(numberone, "d vs ", numbertwo, "d - genes with adjusted p-value < 0.05: ", nrow(sig)))
		sumres <- append(sumres, nrow(sig))
		sig_sorted <- sig[order(sig$log2FoldChange),]

		#####

		# Write to table

		if (nrow(sig) == 0){
			print(paste0("No significantly differentially expressed genes between ", numberone, "d and ", numbertwo, "d"))
		} else {
			outtable <- paste0(compdir, "/", numberone, "_", numbertwo, "d_SigResults.tsv")
			if (!file.exists(outtable)){
				write.table(sig_sorted, file=outtable, quote=FALSE, row.names=TRUE, sep="\t")
			}
		}

		###

		# Positive LFC = upregulated in dormant

		sig_DormantUp <- sig_sorted[which(sig_sorted$log2FoldChange > 0),]
		if (nrow(sig_DormantUp) == 0){
			print(paste0("No significantly upregulated genes in ", numbertwo, "d vs. ", numberone, "d"))
			sumres <- append(sumres, nrow(sig_DormantUp))
		} else {
			print(paste0(numberone, "d vs ", numbertwo, "d - genes upregulated in ", numbertwo, "d: ",
				nrow(sig_DormantUp)))
			sumres <- append(sumres, nrow(sig_DormantUp))
			sig_DormantUp <- sig_DormantUp[order(sig_DormantUp$log2FoldChange, decreasing=TRUE),]
			outtable <- paste0(compdir, "/", numberone, "_", numbertwo, "d_SigResults.", numbertwo, "d_Up.tsv")
			if (!file.exists(outtable)){
				write.table(sig_DormantUp, file=outtable, quote=FALSE, row.names=TRUE, sep="\t")
			}
		}

		sig_DormantDown <- sig_sorted[which(sig_sorted$log2FoldChange < 0),]
		if (nrow(sig_DormantDown) == 0){
			print(paste0("No significantly downregulated genes in ", numbertwo, "d vs. ", numberone, "d"))
			sumres <- append(sumres, nrow(sig_DormantDown))
		} else {
			print(paste0(numberone, "d vs ", numbertwo, "d - genes downregulated in ", numbertwo, "d: ",
				nrow(sig_DormantDown)))
			sumres <- append(sumres, nrow(sig_DormantDown))
			sig_DormantDown <- sig_DormantDown[order(sig_DormantDown$log2FoldChange, decreasing=FALSE),]
			outtable <- paste0(compdir, "/", numberone, "_", numbertwo, "d_SigResults.", numbertwo, "d_Down.tsv")
			if (!file.exists(outtable)){
				write.table(sig_DormantDown, file=outtable, quote=FALSE, row.names=TRUE, sep="\t")
			}
		}

		# Add summary statistics to summary file
		print(paste0("Summary line: ", sumres))
		write.table(sumres, summaryfile, quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t", append=TRUE)

		### GO enrichment analysis

		myInterestingGenes_Up <- rownames(sig_DormantUp)
		geneList_Up <- factor(as.integer(geneNames %in% myInterestingGenes_Up))
		names(geneList_Up) <- geneNames
		if (!"1" %in% levels(geneList_Up)){
			print(paste0("No upregulated genes in ", numbertwo, "d vs. ", numberone, "d."))
		} else {
			for (GOgroup in c("BP", "MF", "CC")) {
				GOdata <- new("topGOdata", ontology = GOgroup, allGenes = geneList_Up, annot = annFUN.gene2GO, gene2GO = geneID2GO)
				resultFis <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
				allRes <- GenTable(GOdata, classic = resultFis, orderBy = "weight", ranksOf = "classic", topNodes = 1000)
				write.table(allRes, paste0(godir, "/", numberone, "_", numbertwo, "d_", GOgroup, ".UpIn", numbertwo, "d.Results.txt"),
					quote = FALSE, sep = "\t", row.names = FALSE)
				#showSigOfNodes(GOdata, score(resultFis), firstSigNodes = 5, useInfo = 'all')
				printGraph(GOdata, resultFis, firstSigNodes = 5, fn.prefix = paste0(godir, "/", numberone, "_", numbertwo, "d_", GOgroup, ".UpIn", numbertwo, "d"),
					useInfo = "all", pdfSW = TRUE)
			}
		}

		myInterestingGenes_Down <- rownames(sig_DormantDown)
		geneList_Down <- factor(as.integer(geneNames %in% myInterestingGenes_Down))
		names(geneList_Down) <- geneNames
		if (!"1" %in% levels(geneList_Down)){
			print(paste0("No downregulated genes in ", numbertwo, "d vs. ", numberone, "d."))
		} else {
			for (GOgroup in c("BP", "MF", "CC")) {
				GOdata <- new("topGOdata", ontology = GOgroup, allGenes = geneList_Down, annot = annFUN.gene2GO, gene2GO = geneID2GO)
				resultFis <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
				allRes <- GenTable(GOdata, classic = resultFis, orderBy = "weight", ranksOf = "classic", topNodes = 1000)
				write.table(allRes, paste0(godir, "/", numberone, "_", numbertwo, "d_", GOgroup, ".DownIn", numbertwo, "d.Results.txt"),
					quote = FALSE, sep = "\t", row.names = FALSE)
				#showSigOfNodes(GOdata, score(resultFis), firstSigNodes = 5, useInfo = 'all')
				printGraph(GOdata, resultFis, firstSigNodes = 5, fn.prefix = paste0(godir, "/", numberone, "_", numbertwo, "d_", GOgroup, ".DownIn", numbertwo, "d"),
					useInfo = "all", pdfSW = TRUE)
			}
		}
	}
}

#####

# plotCounts(dds, gene="", intgroup = "Dormancy")

# GOs
# EC
# KEGG_Pathway
