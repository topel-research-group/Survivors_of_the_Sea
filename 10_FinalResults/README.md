## Results summary
In each control vs. resting comparison, only around 2000-4000 genes are up/downregulated with |LFC| > 1,  
and around 3000-5000 are up/downregulated with |LFC| > 0.585, meaning that the majority do not appear to  
be differentially expressed (see [`04_DiffExp`](../04_DiffExp)),  despite resting stages being 'dormant'  
* DOUBLE-CHECK RESULTS WITHOUT LFC THRESHOLD

## Pathways

### Citric Acid Cycle
a.k.a. TCA cycle, Krebs cycle  
**Relevance:** Energy release from fats, carbohydrates and proteins, and production of NADH and AA precursors

The cycle overall seems to have a 56d uptick, slightly delayed after addition of nutrients, before dropping

| R05 | GF04 |
|-----|------|
| ![R05_TCA_Cycle_partial](../08_GenesToLookFor/GeneChecks/CitricAcidCycle/TCAcycle.R05.lfc0.png) | ![GF04_TCA_Cycle_partial](../08_GenesToLookFor/GeneChecks/CitricAcidCycle/TCAcycle.GF04.lfc0.png) |

### Urea Cycle
a.k.a. ornithine-urea cycle  
**Relevance:** Production of urea, potential nitrogen storage? _Check specifics in literature_

Parts of the cycle are upregulated at all times, some experience an uptick, others are downregulated
* What is the consequence of this regarding accumulation of intermediates?

| R05 | GF04 |
|-----|------|
| ![R05_Urea_Cycle](../08_GenesToLookFor/GeneChecks/UreaCycle/UreaCycle.R05.lfc0.png) | ![GF04_Urea_Cycle](../08_GenesToLookFor/GeneChecks/UreaCycle/UreaCycle.GF04.lfc0.png) |

### Glycolysis
**Relevance**: Ultimately leads to the production of pyruvate and energy

Need to check that all annotations are correct, but very little seems to be regulated; SOME upregulation at 56d?

_Results below are preliminary; some genes may be unaccounted for_

| R05 | GF04 |
|-----|------|
| ![R05_Glycolysis_prelim](../08_GenesToLookFor/GeneChecks/Glycolysis/Glycolysis.R05.lfc0.png) | ![GF04_Glycolysis_prelim](../08_GenesToLookFor/GeneChecks/Glycolysis/Glycolysis.GF04.lfc0.png) |

### Gluconeogenesis
**Relevance**: Production of hexoses from pyruvate as fuel for storage carbohydrate synthesis

Seems to be upregulated overall, to varying degrees

_Some genes still need to be checked_

| R05 | GF04 |
|-----|------|
| ![R05_Gluconeogenesis](../08_GenesToLookFor/GeneChecks/Gluconeogenesis/Gluconeogenesis.R05.lfc0.png) | ![GF04_Gluconeogenesis](../08_GenesToLookFor/GeneChecks/Gluconeogenesis/Gluconeogenesis.GF04.lfc0.png) |

### GS-GOGAT cycle
**Relevance**: Ammonium assimmilation (i.e. downstream from DNRA?)

Two glutamate synthases seem to be upregulated  throughout

_Are there any more glutamine synthetases I'm missing?_

| R05 | GF04 |
|-----|------|
| ![R05_GS-GOGAT_prelim](../08_GenesToLookFor/GeneChecks/GlutamateCycle/GS-GOGAT.R05.lfc0.png) | ![GF04_GS-GOGAT_prelim](../08_GenesToLookFor/GeneChecks/GlutamateCycle/GS-GOGAT.GF04.lfc0.png) |

### Ribosomes
**Relevance**: Potential nitrogen storage, or protein synthesis is in overdrive to ensure survival

Ribosomal proteins are upregulated at all timepoints during the experiment
* R05 has a peak at 49d
* GF04 has a peak at 56d

_Likely missing some ribosomal proteins at the moment_

| R05 | GF04 |
|-----|------|
| ![R05_RibosomalProteins](../08_GenesToLookFor/GeneChecks/RibosomalProteins/RibosomalProteins.R05.lfc0.png) | ![GF04_RibosomalProteins](../08_GenesToLookFor/GeneChecks/RibosomalProteins/RibosomalProteins.GF04.lfc0.png) |

### DNRA
**Relevance**: Phenomenon seen in several diatoms as a means of anaerobic energy generation

Both nitrate reductase (Sm_t00006755-RA) and nitrite reductase (Sm_t00001748-RA) show increased expression at 49d before dropping
* Nitrate reductase shows LFC 3.93 increase at 49d in R05, LFC 8.15 increase at 49d in GF04 (upregulated throughout)
* Nitrite reductase shows LFC 2.66 increase at 49d in R05, LFC 5.35 increase at 49d in GF04
  * Downregulated in all subsequent R05 timepoints, upregulated at 56/72d in GF04 (down at 189d, rest not significant)
* **This is consistent with expectations regarding resting stages being capable of DNRA**

* Nitric oxide reductase (Sm_t00000914-RA and Sm_t00013435-RA?) slightly upregulated; consequences?
* No nitrous oxide reductase found

_Any further conclusions to be drawn?_


| R05 | GF04 |
|-----|------|
| ![R05_Nitrate](../08_GenesToLookFor/GeneChecks/NitrateMetabolism/NitrateMetabolism.R05.lfc0.png) | ![GF04_Nitrate](../08_GenesToLookFor/GeneChecks/NitrateMetabolism/NitrateMetabolism.GF04.lfc0.png) |

### Calvin cycle
**Relevance**: Related to photosynthesis, so would expect downregulation

Results need checking, as many of the genes involved seem to be upregulated, some quite significantly
* This could be due to these genes being involved in other pathways, however

_Check the upregulated genes_


| R05 | GF04 |
|-----|------|
| ![R05_CalvinCycle](../08_GenesToLookFor/GeneChecks/CalvinCycle/CalvinCycle.R05.lfc0.png) | ![GF04_CalvinCycle](../08_GenesToLookFor/GeneChecks/CalvinCycle/CalvinCycle.GF04.lfc0.png) |


### Pentose phosphate pathway
**Relevance**: Parallel process to glycolysis

(Cytosolic in diatoms? [cf. plastidic in plants])


| R05 | GF04 |
|-----|------|
| ![R05_PentosePhosphatePathway](../08_GenesToLookFor/GeneChecks/PentosePhosphatePathway/PentosePhosphatePathway.R05.lfc0.png) | ![GF04_PentosePhosphatePathway](../08_GenesToLookFor/GeneChecks/PentosePhosphatePathway/PentosePhosphatePathway.GF04.lfc0.png) |
