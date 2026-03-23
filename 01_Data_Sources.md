# Data Sources for Mendelian Randomization Analysis
## 1. Exposure: 91 Circulating Inflammatory Proteins 
- Dataset name: Genetics of circulating inflammatory proteins (91 plasma proteins pQTL GWAS)
- Publication PMID: 37673607
- Publication reference: Zhao JH, Stacey D, Eriksson N, et al. Genetics of circulating inflammatory proteins identifies drivers of immune-mediated disease risk and therapeutic targets. Nat Immunol. 2023 Sep;24(9):1540-1551.
- Access link: https://www.ebi.ac.uk/gwas/ (EBI GWAS Catalog, accession numbers: GCST90274758 to GCST90274848)
- Sample size: 14,824 individuals of European ancestry
- Population: European ancestry
- Platform: Olink Target platform (measuring 91 plasma inflammatory proteins)
- Notes: A total of 180 pQTLs (59 cis-pQTL and 121 trans-pQTL) were identified in the original study.

## 2. Outcome: Superficial Thrombophlebitis 
- Dataset name: FinnGen Consortium R11 release (I9_PHLETHROM)
- Publication PMID: 36702820 (FinnGen main publication)
- Publication reference: Kurki MI, Karjalainen J, Palta P, et al; FinnGen. FinnGen provides genetic insights from a well-phenotyped isolated population. Nature 2023;613(7944):508–518.
- Access link: https://www.finngen.fi/en/access_results (FinnGen public data portal, direct summary stats link: https://storage.googleapis.com/finngen-public-data-r11/summary_stats/finngen_R11_I9_PHLETHROM.gz)
- Sample size: 453,733 participants (8,602 STP cases, 392,860 controls)
- Population: Finnish (European ancestry)
- Phenotype definition: Based on ICD-8 (451), ICD-9 (4519X|4518X), and ICD-10 (I80) codes
- Notes: Excluded cases with specific ICD codes as per FinnGen R11 quality control criteria.

## 3. Data Access Instructions 
- All GWAS summary data used in this study are publicly available.
- For the 91 inflammatory proteins dataset: 
  1. Visit the EBI GWAS Catalog (https://www.ebi.ac.uk/gwas/).
  2. Search for accession numbers GCST90274758 to GCST90274848.
  3. Download the summary statistics files for each protein.
- For the STP dataset: 
  1. Visit the FinnGen public data portal (https://www.finngen.fi/en/access_results).
  2. Navigate to the R11 release.
  3. Download the summary statistics file for "I9_PHLETHROM" directly from the link above.