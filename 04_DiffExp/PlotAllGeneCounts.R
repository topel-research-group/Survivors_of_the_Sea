#!/usr/bin/env Rscript

library("DESeq2")
library("ggplot2")

# Read in both R05 and GF04 results

r05dds <- readRDS("R05_FINAL/dds.rds")
gf04dds <- readRDS("GF04_FINAL/dds.rds")

# Get list of genes present in these analyses

r05list <- row.names(r05dds)
gf04list <- row.names(gf04dds)

# Create master dataframe for all results to be read into

masterdf <- data.frame(count=double(),
			Dormancy=factor(),
			sample=character(),
			gene=character(),
			strain=character())

# Go through all results, outputting the mean of normalised counts
# and appending them to the master dataframe

for (mygene in r05list) {
	print(mygene)

	thisgene <- plotCounts(r05dds, gene=mygene, intgroup = "Dormancy", returnData=T)
	thisgene$sample <- row.names(thisgene)
	thisgene$gene <- mygene
	thisgene$strain <- "R05AC"
	row.names(thisgene) <- NULL

	masterdf <- rbind(masterdf, thisgene)
}
for (mygene in gf04list) {
	print(mygene)

	thatgene <- plotCounts(gf04dds, gene=mygene, intgroup = "Dormancy", returnData=T)
        thatgene$sample <- row.names(thatgene)
        thatgene$gene <- mygene
        thatgene$strain <- "GF0410J"
	row.names(thatgene) <- NULL

	masterdf <- rbind(masterdf, thatgene)
}

write.table(masterdf, file="AllMeanCounts.lst", quote=FALSE, row.names=FALSE, sep="\t")

#####

# test on Sm_t00012191-RA

masterdf$gene <- as.factor(masterdf$gene)

for (mygene in levels(masterdf$gene)) {
	print(mygene)

	thisgene <- masterdf[masterdf$gene == mygene,]
	outgraph <- paste0("NormalisedCountGraphs_FINAL/", mygene, ".png")

	png(outgraph, height=500, width=500)
	p <- ggplot(thisgene, aes(x=Dormancy, y=count, colour=strain, shape=strain)) + 
		labs(x="Time (days)", y="Mean of normalised counts", title=mygene) +
		theme(plot.title = element_text(hjust=0.5)) +
		geom_point(size=2)
	print(p)
	dev.off()
}
