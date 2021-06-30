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
  * Sm_t00003059-RA - best BLAST match to T. oceanica hypothetical protein THAOC-30017 (52% ID, 84% cov)
    * No conserved domains
  * Sm_t00009227-RA - best BLAST match to T. pseudonana predicted protein (41.3% ID, 41% cov)
    * Some hits (~36% ID) to bacterial DsrE family, methyltransferase domains, UbiE; no conserved domains though
    * InterProScan notes 'Prokaryotic membrane lipoprotein lipid attachment site'
  * Sm_t00018719-RA - PA14
    * All hits are bacterial, conserved domain hit to part of FTP ('fucolectin tachylectin-4 pentaxrin-1 domain')
    * eggnog mapper hits to Aquimarina pnlA (pectin lyase?), and labels as PA14
    * EC numbers correspond to lecithinase C (phospholipase) and pectin lyase
      * Related to changes in cell membrane, if it's a phospholipase?
  * Sm_t00021541-RA
    * Best BLAST match to F. cylindrus hypothetical protein FRACYDRAFT_208194 (30.66% ID, 20% cov)
* Sm_t00009316-RA (heat shock factor) and Sm_t00019204-RA (chitin synthase activity) are both in the top ten of 5/6 comparisons
* Checking top 30 hits for all R05 vegetative vs. resting upreg. comparisons
  * Sm_t00009162-RA (6/6) - NO BLAST hits, but I have it labelled as 'LOW QUALITY PROTEIN: ecto-ADP-ribosyltransferase 5-like'
  * Sm_t00019202-RA (5/6) - Trypsin-like serine protease; T9SS C-terminal target domain-containing protein?
  * Sm_t00008968-RA (5/6) - Peptidase S8
  * Sm_t00006726-RA (5/6) - Belongs to the class I-like SAM-binding methyltransferase superfamily. RNA M5U methyltransferase family
    * Significance?
  * Sm_t00006457-RA and Sm_t00006458-RA (5/6) - Best BLAST hit to T. oceanica hypothetical protein THAOC_23441
  * Sm_t00004483-RA (5/6) - best BLAST match to T. oceanica hypothetical protein THAOC_25469 (49.3% ID, 55% cov)
  * Sm_t00000917-RA (5/6) - best BLAST match to T. pseudonana predicted protein (46.03% ID, 53% cov)

* Other genes seemingly upregulated in all GF04 resting stages vs. vegetative, and in the top ten in terms of upregulation
  * Sm_t00002056-RA - YjgF/chorismate_mutase-like, putative endoribonuclease
  * Sm_t00006630-RA - best BLAST match to T. pseudonana predicted protein (51.14% ID, 34% cov)
    * Conserved PHD finger domain; transcription factor? DNA binding?
  * Sm_t00009980-RA - best BLAST match to T. pseudonana predicted protein (70.19% ID, 13% cov)
    * Putative chitinase in F. cylindrus, albeit 7% cov...
    * The original model is much shorter than its true length, noted during manual curation
  * Sm_t00012155-RA - Rit1 N-terminal domain
    * Best BLAST hit to T. oceanica hypothetical protein THAOC_21179 (47.59% ID, 98% cov)
    * Other BLAST hits with good coverage (though ~30% ID) code for tRNA phosophoribosyl transferases (e.g. Rit1)
    * In yeast, rit1 null mutation allows initiator tRNA to be used in elongation as well
      * Consequence of overexpression/upregulation?
  * Sm_t00021963-RA - best BLAST matches to leucine-rich repeat proteins
* Sm_t00009520-RA (Bacterial protein of unknown function (DUF839)) is in the top ten of 5/6 comparisons

* Other possibly interesting genes (checking top 30 upregulated genes in all vegetative vs. resting comparisons)
  * Sm_t00017257-RA (10/12 comparisons) - best BLAST match to T. oceanica hypothetical protein THAOC_04206 (47.21% ID, 73% cov)


* Sm_t00009759-RA is the only gene in R05 which is significantly downregulated in 189d vs. 126d
  * Down in 49d, up again in 56d, then gradually down
  * Same general pattern in GF04, albeit with more spread of normalised counts so less significant
* Sm_t00002815-RA is the only gene in GF04 which is significantly downregulated in 189d vs. 126d
  * Gradual downregulation over time
  * Same general pattern in R05, albeit with more spread of normalised counts so less significant


* Regarding GO terms, ribosome-related GO terms appear to be upregulated in ALL resting stages vs. vegetative
  * Consistent with Pelosi's findings in Chaetoceros socialis (doctoral thesis, 2018)
