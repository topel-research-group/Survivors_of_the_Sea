# Mapping the RNAseq reads to the reference genome to get an idea of how well everything maps

Compare with mapping to the transcriptome to see how much additional data is lost

Problem - GF04 maps around 10% worse to the (R05) reference than R05
* Also need to consider heterozygous sites... Relax mapping parameters?
* Discussion on use of non-native reference genome (in bacteria): https://doi.org/10.1371/journal.pone.0180904
  * Long reads are preferred
  * Relaxed parameters could lead to false positives...

At the very least, within-strain comparisons can be performed
