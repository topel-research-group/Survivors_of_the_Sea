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
**Relevance:** Production of urea, potential nitrogen storage? __check specifics in literature__

Parts of the cycle are upregulated at all times, some experience an uptick, others are downregulated
* What is the consequence of this regarding accumulation of intermediates?

| R05 | GF04 |
|-----|------|
| ![R05_Urea_Cycle](../08_GenesToLookFor/GeneChecks/UreaCycle/UreaCycle.R05.lfc0.png) | ![GF04_Urea_Cycle](../08_GenesToLookFor/GeneChecks/UreaCycle/UreaCycle.GF04.lfc0.png) |
