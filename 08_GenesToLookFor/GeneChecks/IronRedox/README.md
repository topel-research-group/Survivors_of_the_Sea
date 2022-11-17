# Oxidation of Fe2+ and reduction of Fe3+

|           Annotation            |       Gene      |                   Regulation                    |       Notes       |
|---------------------------------|-----------------|-------------------------------------------------|-------------------|
| Ferric reductase family protein | Sm_t00009056-RA*| Down throughout in GF04 (some late down in R05) | V. low expression |
|                                 | Sm_t00010479-RA | Up throughout in GF04 (some late up in R05)     | (>100 reads each) |
|---------------------------------|-----------------|-------------------------------------------------|-------------------|
| Ferric reductase-like           | Sm_t00018630-RA | 56d down in GF04; 49+189 up in R05              |                   |
|  transmembrane component        |                 |	                                              |                   |

* Sm_t00009056-RA is a very messy model, on a long region of poor assembly

## FeGenie

Results of FeGenie may need further inspection; current results seem to hit odd proteins

## Review - 'Iron Uptake Mechanisms in Marine Phytoplankton' (Sutak et al., 2020)

"The results of these studies suggest that P. tricornutum makes use of at least three different inducible iron uptake pathways:
* a phytotransferrin-mediated non-reductive uptake pathway (Fe’ model), 
* a reductive uptake pathway (reductive model) and 
* a siderophore-mediated uptake pathway."

Most of the genes mentioned don't appear to be present in T. pseudonana!

|   P. tricornutum gene    |                  |                                      |                                   |
|--------------------------|------------------|--------------------------------------|-----------------------------------|
| Ferric reductase    FRE1 | Sm_t00010479-RA? | All up in GF04; 189d up in R05       | Both of these have okay hits      |
|                     FRE2 | Sm_t00018630-RA? | 56d down in GF04; 49+189 up in R05   |  to both FRE1 (tail end) and FRE2 |
|                     FRE3 | Sm_t00013368-RA? | Some wobble in GF04; all up in R05   |                                   |
|                          | Sm_t00015390-RA? | All down in both strains (esp 49d)   |                                   |
|                          | Sm_t00002711-RA? | All down in GF04; some wobble in R05 |                                   |
|                     FRE4 | Sm_t00000378-RA  | All down (except 56d in R05)         | 56d peak in GF04, too             |
|--------------------------|------------------|--------------------------------------|-----------------------------------|
| Putative ferrichrome-    | No good hits     |                                      |                                   |
|  binding protein FBP1    |                  |                                      |                                   |
|--------------------------|------------------|--------------------------------------|-----------------------------------|
| Irt-like protein of      | Sm_t00018447-RA  | All up in R05 only                   | Annotated as zinc transporter,    |
|  the ZIP family          |                  |                                      | but could be Fe/Zn?               |
|--------------------------|------------------|--------------------------------------|-----------------------------------|
| Iron starvation-  ISIP1  | No good hits     |                                      |                                   |
|  induced proteins ISIP2a | Sm_t00018628-RA  | Late upreg in GF04 (126d, 189d)      | Assembly problems here, and       |
|                          |                  | Up at 49, 72, 91, 126d in R05        | annotated as predicted protein    |

## Other

Based on the above result regarding a possible iron/zinc transporer labelled purely as zinc transporter, *check other zinc transporters!*

|                    |                 |                                                      |                              |
|--------------------|-----------------|------------------------------------------------------|------------------------------|
| Iron transporter   | Sm_t00005368-RA | Dramatically upregulated throughout in both strains! | Definitely iron transporter? |
|--------------------|-----------------|------------------------------------------------------|------------------------------|
| 'Zinc transporter' | Sm_t00000096-RA | Dramatically upregulated throughout, esp. in R05     |
|                    | Sm_t00001624-RA | Up at all except 49d in GF04, up at 189d in R05      |
|                    | Sm_t00001766-RA | Slightly downregulation at 91d and 189d in GF04      |
|                    | Sm_t00002559-RA | No significant change                                |
|                    | Sm_t00002700-RA | Up throughout in GF04, esp. 49d (59->189 up in R05)  |
|                    | Sm_t00004040-RA | No significant change                                |
|                    | Sm_t00010759-RA | No significant change                                |
|                    | Sm_t00012102-RA | Down from 72d on in R05                              |
|                    | Sm_t00016646-RA | Down at 49/91/126d in GF04                           |
|                    | Sm_t00017696-RA | Up throughout in both strains                        |
|--------------------|-----------------|------------------------------------------------------|------------------------------|
| Flavodoxin         | Sm_t00007335-RA | Down throughout in R05; down at 49d in GF04          |
|                    | Sm_t00020531-RA | Up at 91d in R05                                     | V. low expression
| (-family protein)  | Sm_t00003819-RA | Up at all except 56d in R05; up at 49/126d in GF04   | Half cov. region

######

## Fe and Metal Uptake and Homeostasis (plus a few other highly regulated genes)

| Annotation                                                             | ThaPse | PhaTr2 |
|------------------------------------------------------------------------|--------|--------|
| FRE1 (Ferric reductase 1)                                              |    -   | 54486  |
| FRE2 (Ferric reductase 2)                                              |    -   | 46928  |
| FRE3 (Cytochrome b561 / ferric reductase transmembrane domain protein) |    -   | 54940  |
| FRE4 (Cytochrome b561 / ferric reductase transmembrane domain protein) |  5335  | 43682  |
| FBP1 (Ferrichrome binding protein)                                     |    -   | 46929  |
| zupT (Fe/Zn divalent cation permease)                                  | 268192 | 38445  |
| HMA1 (P1B-type heavy metal translocating P-type ATPase)                | 261657 | 52367  |
| CDF1 (Cation diffusion facilitator)                                    | 260789 | 45119  |
| CDF2 (Cation diffusion facilitator)                                    |  4507  | 54552  |
| Heme oxygenase (decyclizing)                                           | 17865  | 12588  |
| Flavodoxin                                                             | 19141  | 23658  |
| Flavodoxin-like encoding gene                                          | 24306  |  7755  |
| IscA1 (Iron-sulfur cluster assembly protein)                           | 39032  | 14867  |
| 2Fe-2S Ferredoxin domain protein                                       |  1665  | 46444  |
| Cytochrome c biogenesis protein                                        | 260772 | 54791  |
| Rieske-type ferredoxin subunit                                         | 29842  |  9046  |
| Cytochrome b5 electron transport protein                               | 30887  | 35395  |
| Haem peroxidase                                                        | 262753 | 47395  |
| Fe/Mn SOD (Superoxide dismutase)                                       | 40713  | 42832  |
| Fe/Mn SOD (Superoxide dismutase)                                       | 32874  | 12583  |
| ChrA (Chromate efflux pump)                                            | 262849 | 48511  |
| ISIP1 (Fe starvation induced protein 1)                                |    -   | 55031  |
| ISIP2A                                                                 |    -   | 54465  |
| ISIP2B                                                                 |    -   | 54987  |
| ISIP3 (bacterial DUF305 domain)                                        |    -   | 47674  |
