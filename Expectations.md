# Genes/pathways we expect to change

## Fatty acid desaturase

This is one of the genes most upregulated in all resting stages, and is consistent with previous findings
* See e.g. [Kuwata et al., 1993](https://doi.org/10.3354/meps102245)


## Mechanisms for uptake of nitrate and ammonium
* DiAMT1/DiNRT2?


## Mechanisms for anaerobic fermentation?


## Upregulation of DNRA?
* Nitrate reductase - 1.7.5.1/1.9.6.1
* Nitrite reductase - 1.7.2.1/1.7.2.2
* Nitrite oxide reductase - 1.7.5.2/1.7.2.5
* Nitrous oxide reductase - 1.7.2.4

None of the above EC numbers are present in the v1.1.2 eggnog mapper file

Sm_t00011895-RA labelled as 'nitrate reductase activity'
* This is downregulated in the vegetative cells??
Sm_t00006755-RA seems to be a nitrate reductase but may be assimilatory
* Seems to be upregulated in resting stages

Nitrate/nitrite reductases - Nap/Nrf and Nar/Nir?

## Any significant downregulation signal?

##########

## Abscisic acid (ABA)?
Important in resting stages of dinoflagellates and higher plants
* Looking at NCBI, some diatom sequences are present for genes in this pathway
 * Coscinodiscus granii? (Kentzer and Mazur, 1991)

Pathway genes (noted for dinoflagellate, as starting point; see Deng et al., 2017)
* Zeaxanthin epoxidase (ZEP); present in T. pseudonana
  * Sm_g00013199? (Sm_000038F) - Downregulated in later resting stages?
  * Sm_g00012687? (Sm_000029F) - Downregulated in all resting stages?
  * Sm_g00005941? (Sm_000011F) - Downregulated in all resting stages?
  * In S. trochoidea, ZEP is UPREGULATED in cysts...

* 9-cis-epoxycarotenoid dioxygenases (NCED)
  * Sm_g00012791? (Sm_000027F) - Upregulated quite significantly in all resting stages
    * Best BLAST hit to the NCED sequence of P. tricornutum
  * Possibly Sm_g00003413? (Sm_000003F) - Downregulated in later resting stages
  * In S. trochoidea, NCED is highest in immature cysts

* Abscisic-aldehyde oxidase (AAO)
  * Sm_g00005740 (Sm_000007F) - Weird signal: downregulated in all resting stages except 189d, and upregulated between resting stages as time goes on
    * Best BLAST hit to the AAO sequence of (dinoflagellate) Scrippsiella trochoidea
    * Strong domain hits to xanthine dehydrogenase/aldehyde oxidase
  * In S. trochoidea, no significant difference between stages

* ABA-8'-hydroxylase (ABAH) [CYP707A family genes]
  * Sm_g00014680? (Sm_000034F) - Weird signal: upregulated in 49d, but downregulated between 49d and later resting stages
    * BLASTs to cytochrome P450
      * ABA-8'-hydroxylase IS a cytochrome P450 protein
  * Sm_g00000943? (Sm_000000F) - downregulated quite significantly among resting stages
  * Many other decent-length hits
  * In S. trochoidea, highest in vegetative cells and immature cysts

## Other phytohormones and their pathways?
