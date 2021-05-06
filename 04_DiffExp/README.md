# Differential expression analysis of vegetative vs. resting cells

* R05/GF04 - require one read in at least three samples; no LFC filtering
* R05_v2/GF04_v2 - require one read in at least three samples; LFC threshold 1/-1
* R05_v3/GF04_v3 - require one read in at least TWO samples; LFC filtering 1/-1

v2 and v3 are somewhat similar in terms of gene numbers

## Observations by category

1. What is turned off when resting (only present in vegetative)?
2. What is turned on when resting (totally absent in vegetative)?
3. Which genes increase in expression throughout resting?
4. Which genes decrease in expression throughout resting?
5. Which genes fall into subcategories of the above (e.g. only on during the final timepoint)?


## Initial observations

* Based on the heatmap and PCA, the vegetative cells have VERY different expression patterns versus the resting stages
* R05 189d relatively distant from other resting timepoints, AND from one another...
  * Worth redoing the analysis without 189d?
* GF04 resting stages seem to cluster together rather more cleanly...

Some concerns; looking at downregulated-in-resting genes...
* Sm_t00017228-RA appears to be the top downregulated-in-56d gene in R05, but looking at the plotCounts plot, 2/3 of the vegetative
  samples show no expression...

## Interesting results (based on expression of v1.1.2 models)

* Sm_t00012191-RA appears to be a fatty acid desaturase, and is upregulated in all resting stage timepoints in both strains
  * This would seem to be consistent with existing findings on resting stage storage compounds- from Ellegaard & Ribeiro, 2017:
    "In C. curvisetus (and other diatoms) changes in the composition of fatty acids were recorded,
     mainly consisting of increased amounts of neutral lipids and UNSATURATED FATTY ACIDS (e.g. Kuwata et al., 1993)."

* Other genes seemingly upregulated in all R05 resting stages vs. vegetative, and in the top ten in terms of upregulation
  * Sm_t00003059-RA
  * Sm_t00009227-RA
  * Sm_t00018719-RA - PA14
  * Sm_t00021541-RA
* Sm_t00009316-RA (heat shock factor) and Sm_t00019204-RA (chitin synthase activity) are both in the top ten of 5/6 comparisons

* Other genes seemingly upregulated in all GF04 resting stages vs. vegetative, and in the top ten in terms of upregulation
  * Sm_t00002056-RA - YjgF/chorismate_mutase-like, putative endoribonuclease
  * Sm_t00006630-RA
  * Sm_t00009980-RA
  * Sm_t00012155-RA - Rit1 N-terminal domain
  * Sm_t00021963-RA
* Sm_t00009520-RA (Bacterial protein of unknown function (DUF839)) is in the top ten of 5/6 comparisons

* Sm_t00009759-RA is the only gene in R05 which is significantly downregulated in 189d vs. 126d
  * Down in 49d, up again in 56d, then gradually down
  * Same general pattern in GF04, albeit with more spread of normalised counts so less significant
* Sm_t00002815-RA is the only gene in GF04 which is significantly downregulated in 189d vs. 126d
  * Gradual downregulation over time
  * Same general pattern in R05, albeit with more spread of normalised counts so less significant


* Regarding GO terms, ribosome-related GO terms appear to be upregulated in ALL resting stages vs. vegetative
  * Consistent with Pelosi's findings in Chaetoceros socialis (doctoral thesis, 2018)
