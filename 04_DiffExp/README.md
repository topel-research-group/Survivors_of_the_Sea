# Differential expression analysis of vegetative vs. resting cells

## Groups

| Strain | Incubation (days) | n |
|--------|-------------------|---|
|  R05   |    0 (control)    | 3 |
|  R05   |      49 (T1)      | 3 |
|  R05   |      56 (T2)      | 3 |
|  R05   |      72 (T3)      | 3 |
|  R05   |      91 (T4)      | 3 |
|  R05   |     126 (T5)      | 3 |
|  R05   |     189 (T6)      | 2 |
|--------|-------------------|---|
|  GF04  |    0 (control)    | 3 |
|  GF04  |      49 (T1)      | 3 |
|  GF04  |      56 (T2)      | 3 |
|  GF04  |      72 (T3)      | 3 |
|  GF04  |      91 (T4)      | 3 |
|  GF04  |     126 (T5)      | 3 |
|  GF04  |     189 (T6)      | 2 |

## Which comparisons to make?

* R05 vegetative vs. R05 resting (each timepoint separately?)
* GF04 vegetative vs. GF04 resting (each timepoint separately?)
* All vegetative vs. All resting?

## Which transformation to use for heatmaps/PCA plots?
* vst or rlog?
  * Notes from [DESeq2 paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8):
    * "The rlog transformation accounts for variation in sequencing depth across samples"  
    * "... while the VST is also effective at stabilizing variance, it does not directly take into account  
      differences in size factors; and in datasets with large variation in sequencing depth (dynamic range of  
      size factors ≳4) we observed undesirable artifacts in the performance of the VST."  
    * "A disadvantage of the rlog transformation with respect to the VST is, however, that the ordering of genes  
      within a sample will change if neighboring genes undergo shrinkage of different strength."

* Blind or not?
  * Likely NOT
  * Note from [DESeq2 vignette](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html):
    * "... blind dispersion estimation is not the appropriate choice if one expects that many or the majority of  
      genes (rows) will have large differences in counts which are explainable by the experimental design..."

## Samples to check
* R05 control _103
  * Very separated from other controls in rlog, less so in vst
* GF04 49d _112 (sorts itself with 72d samples...)
  * Very separated from other samples in rlog, but not in vst...
  * Batch effect?
  * Would need to remove a GF04 control _106 too, as this would be the only one in the other batch if so...

Would I need to remove the B samples from R05 too if I removed them from GF04...?
