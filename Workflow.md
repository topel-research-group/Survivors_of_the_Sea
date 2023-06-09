# Differential expression analysis of vegetative vs. resting stage _S. marinoi_

## Workflow:

1. QC and trimming
  a. Cutadapt - Illumina adaptor removal, quality trimming threshold 30, minimum read length 50bp
  b. SortMeRNA - removal of rRNA sequences based on the S. marinoi genome (vs. SSU/LSU clusters from nuclear, plastid and mito-genomes)
  c. PRINSEQ - trimming of long (10+) polyA/T tails from reads, minimum read length 50bp

2. Transcriptome assembly (separately for R05 and GF04)
  a. Trinity - de novo assembly, inclusion of jaccard_clip parameter to minimise fusion transcripts
  b. Transrate - verified versus the trimmed reads
	* Run vs. Maker models?
  c. CD-HIT - removed a small number of transcripts with 100% identity (reduced number of multi-mappers)

3. Mapping reads to transcriptome
  a. Bowtie2 - `--no-mixed` and `--no-discordant` flags used, as these are used in Trinity's `align_and_estimate_abundance.pl` utility script,
     and should stop reads mapping to different transcripts

4. Quantification of transcripts
  a. `featureCounts` in R, with the following parameter choices
    i. `allowMultiOverlap=FALSE` - where two genes overlap and the reads map entirely to the overlap, this avoids ambiguity
    ii. `largestOverlap=TRUE` - in cases where the overlap isn't complete, this should resolve the ambiguity
    iii. `countMultiMappingReads=TRUE` - after Googling this parameter, `TRUE` seems like the best choice to avoid bias against homologues

5. Differential expression analysis
  a. DESeq2 pipeline
    i. 20 samples per strain, so visualising with rlog transformation to assess variance, rather than variance stabilising transformation
       The DESeq2 manual also states that "the rlog transformation inherently accounts for differences in sequencing depth",
       and in our case, there are between 9.4M and 48M read pairs per sample

### Previous steps

4. ~Quantification of transcripts~  
  a. ~`salmon quant`, attempted using two different methods:~  
    i. ~Mapping-based quantification (using the reads and an index of the transcriptome)~  
    ii. ~Alignment-based quantification (using the (unsorted) BAM file outputted by Bowtie2 above)~

