* Should strict or non-strict mapping be used?
  * ~6% difference in overall mapping, but eyeballing IGV it seems to make little overall difference
  * Go with stricter mapping (--no-mixed --no-discordant), hopefully this will remove misleading mappings?
* Should allowMultiOverlap / fraction parameters be used?
  * ~2.5% more reads mapping when MultiOverlap is allowed
  * 7/9 genes from the top of the output are shared between the test outputs
  * fraction can't be used as DESeq2 only accepts integer values
  * Don't allow multi-overlap; as above, this should avoid issues with ambiguity

```
allcounts <- featureCounts(coldata$files,
		annot.ext=inannot,
		isGTFAnnotationFile=TRUE,	# Works even with GFF file
		GTF.featureType="mRNA",
		GTF.attrType="ID",
		useMetaFeatures=TRUE,		# Seems to make no difference whether TRUE or FALSE
		allowMultiOverlap=FALSE,	# Avoid ambiguous reads
		largestOverlap=TRUE,		# Attempting to avoid ambiguity
		countMultiMappingReads=TRUE,	# Avoid bias against homologs?
		genome=infasta,
		isPairedEnd=TRUE,
		requireBothEndsMapped=FALSE,	# Irrelevant, as both ends must be mapped given the new mapping options
		nthreads=10,
		verbose=TRUE)
```
