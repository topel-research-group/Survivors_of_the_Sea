#!/usr/bin/env Rscript

# Usage:
# Rscript list2heatmap.altered.R --input infile.lst --lfc num --output outfile.svg
# Fixed colour scale (+/- 11.5 LFC), no numbers on cells

# Note: this only looks at vegetative-vs-resting comparisons,
# so significant resting-vs-resting differences won't be reflected


library(argparser)
library("DESeq2")
library("pheatmap")

parser <- arg_parser("Arguments")
parser <- add_argument(parser, "--input", type="character", help="Input file")
parser <- add_argument(parser, "--lfc", type="double", help="LFC filter threshold")
parser <- add_argument(parser, "--output", type="character", help="Output file")
opt <- parse_args(parser)

# Define input file (make this a command line argument)

infile <- opt$input

# Import the relevant data

dds <- readRDS("dds.rds")

# Read in and format the list of genes

ingenes <- read.csv(infile, header=F, row.names = 1, sep="\t", stringsAsFactors=FALSE)
colnames(ingenes) <- "Name"
ingenes$Name <- paste0(ingenes$Name, " (", rownames(ingenes), ")")

# Add columns to the dataframe for the timepoint LFCs
ingenes$Day49 <- 0
ingenes$Day56 <- 0
ingenes$Day72 <- 0
ingenes$Day91 <- 0
ingenes$Day126 <- 0
ingenes$Day189 <- 0

# For each resting stage comparison:

timepoints <- c(49, 56, 72, 91, 126, 189)

for (resting in timepoints) {
	restingcol <- paste0("Day", resting)

	# Run the res comparison and filter results
	res <- results(dds, alpha = 0.05, contrast=c("Dormancy", resting, 0), lfcThreshold = opt$lfc)
	sig <- res[which(res$padj < 0.05), ]

	# For each gene of interest:
	for (mygene in rownames(ingenes)) {

		# Extract the log fold change and send it to the relevant part of the data frame
		ingenes[mygene,restingcol] <- sig[mygene,"log2FoldChange"]
	}
}

# Output a heatmap of these LFCs

## First, ensure the height of the plot will be correct
## Code from https://stackoverflow.com/questions/61874876/get-size-of-plot-in-pixels-in-r

get_plot_dims <- function(heat_map)
{
  plot_height <- sum(sapply(heat_map$gtable$heights, grid::convertHeight, "in"))
  plot_width  <- sum(sapply(heat_map$gtable$widths, grid::convertWidth, "in"))
  return(list(height = plot_height, width = plot_width))
}

### Code for adjusting the breaks required to ensure the colouration centres on white (prevents confusing 'negative warm' values)
### Code from https://stackoverflow.com/questions/31677923/set-0-point-for-pheatmap-in-r

#myrange <- 10
myrange <- 11.5
#lfcmap <- pheatmap(ingenes[-1], breaks=seq(-myrange, myrange, length.out = 100), cluster_rows=F, cluster_cols=F, labels_row=ingenes$Name, cellheight=20, display_numbers = TRUE)
#lfcmap <- pheatmap(ingenes[-1], breaks=seq(-myrange, myrange, length.out = 100), cluster_rows=F, cluster_cols=F, labels_row=ingenes$Name, cellheight=20, display_numbers = FALSE)
lfcmap <- pheatmap(ingenes[-1], breaks=seq(-myrange, myrange, length.out = 100), cluster_rows=F, cluster_cols=F, labels_row=ingenes$Name, cellheight=20, cellwidth=20, display_numbers = FALSE)
plot_dims <- get_plot_dims(lfcmap)

#png(opt$output, height=plot_dims$height, width=plot_dims$width, units="in", res=72)
svg(opt$output, height=plot_dims$height, width=plot_dims$width)
lfcmap
dev.off()
