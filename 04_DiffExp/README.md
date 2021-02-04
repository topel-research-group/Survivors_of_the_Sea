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
* Blind or not?

## Samples to check
* R05 control _103
* GF04 49d _112 (sorts itself with 72d samples...)
  * Batch effect?
  * Would need to remove a GF04 control _106 too, as this would be the only one in the other batch if so...

Would I need to remove the B samples from R05 too if I removed them from GF04...?
