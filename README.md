R Version
R ≥ 4.4.0 (consistent with the original study's analysis environment)
R Packages
All required packages, version requirements, and installation instructions are detailed in requirements.txt. The pipeline automatically detects and installs missing packages via a built-in function.
Data Sources
All GWAS summary statistics used in this study are from publicly available, peer-reviewed datasets.
1. Exposure: 91 Circulating Inflammatory Proteins
Source: Zhao JH, et al. (2023) Nat Immunol (PMID: 37673607)
Access: EBI GWAS Catalog (accession numbers: GCST90274758 to GCST90274848)
Sample Size: 14,824 individuals of European ancestry
2. Outcome: Superficial Thrombophlebitis (STP)
Source: FinnGen Consortium R11 Release (phenotype code: I9_PHLETHROM)
Access: FinnGen Public Data Portal
Sample Size: 453,733 individuals of European ancestry (8,602 STP cases, 392,860 controls)
Full data access instructions are provided in 01_Data_Sources.md.
Usage Instructions
Step 1: Prepare Input Data
Download exposure and outcome GWAS data following the instructions in 01_Data_Sources.md
Prepare your id.txt file: a tab-delimited text file with a header id containing the GWAS accession numbers for your exposure datasets
Ensure the 1000 Genomes European (EUR) PLINK reference panel (for LD clumping) is downloaded to your local system
Step 2: Configure File Paths
Open MR_analysis.R and Supplementary_code.R, and update the global path parameters at the top of each script to match your local file directory structure (base directory, raw data paths, PLINK reference panel path, etc.).
Step 3: Run the Analysis Pipeline
1. Core MR Analysis
Run the full MR analysis pipeline in R:
source("MR_analysis.R")
2. Result Aggregation & Visualization
Run the supplementary code to aggregate results and generate publication-ready figures:
source("Supplementary_code.R")
Output Description
Core MR Analysis Output
exp_dat/: Filtered IV data for each exposure (after P-value filtering and LD clumping)
result_filter/: Per-exposure MR results (multi-sheet Excel files containing: filtered exposure data, harmonized data, full MR results, OR values, IVW+MR-Egger results, heterogeneity tests, pleiotropy tests, single-SNP analysis, and leave-one-out analysis)
MR_summary.xlsx: Aggregated odds ratio (OR) and 95% confidence interval (CI) results for all exposures
Visualization Output
Aggregated MR results table for heatmap generation
Circular heatmap of normalized MR results across 5 analytical methods
Forest plot of OR and 95% CI for the primary IVW analysis (labeled protective/risk factors)
Key Notes
Population Consistency: All GWAS data are restricted to individuals of European ancestry to avoid population stratification bias.
Instrument Strength: Only SNPs with an F-statistic > 10 are retained in the final analysis to eliminate weak instrument bias.
Path Configuration: All hard-coded file paths must be updated to match your local directory structure before running the scripts.
Error Handling: The pipeline uses tryCatch to avoid workflow interruption from single-file errors, with all error messages printed to the R console.
Parameter Adjustment: All analytical thresholds (P-value for IV selection, LD clumping parameters, F-statistic cutoff) are clearly annotated in the scripts and can be modified as needed.