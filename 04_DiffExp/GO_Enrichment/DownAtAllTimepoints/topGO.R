#!/usr/bin/env Rscript

library("topGO")

# Tab-separated file of gene ID (Sm_...) and associated GO terms
geneID2GO <- readMappings("../R05_v1.1.8.GOterms.AllGenes.IncOrganelles.tsv")
geneNames <- names(geneID2GO)

interesting_genes <- as.character(read.csv("DownAtAllTimepoints.lst", header=FALSE, sep="\t")$V1)
geneList <- factor(as.integer(geneNames %in% interesting_genes))
names(geneList) <- geneNames

for (GOgroup in c("BP", "MF", "CC")) {
	GOdata <- new("topGOdata", ontology = GOgroup, allGenes = geneList,
		annot = annFUN.gene2GO, gene2GO = geneID2GO)
	resultFis <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
	allRes <- GenTable(GOdata, classic = resultFis, orderBy = "weight",
		ranksOf = "classic", topNodes = 1000)
	write.table(allRes, paste0("DownInAll.", GOgroup, ".txt"), quote = FALSE, sep = "\t", row.names = FALSE)
	#showSigOfNodes(GOdata, score(resultFis), firstSigNodes = 5, useInfo = 'all')
	printGraph(GOdata, resultFis, firstSigNodes = 5, fn.prefix = paste0("DownInAll.", GOgroup),
		useInfo = "all", pdfSW = TRUE)
}

