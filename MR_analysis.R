# ==============================================================================
# Mendelian Randomization (MR) Analysis Pipeline
# Function: Exposure data filtering → Core MR analysis → Sensitivity analysis → Result summary
# Description: A complete pipeline for two-sample MR analysis, including LD clumping, 
#              harmonization, IV strength filtering, core MR tests, and sensitivity analyses.
# Notes: 
#   - All file paths must be modified according to your local storage!
#   - Ensure PLINK binary and 1000 Genomes reference panel are correctly configured.
# ==============================================================================

# -------------------------- 1. Load Required Packages --------------------------
# Install packages first if not installed:
# install.packages(c("data.table","TwoSampleMR","ggplot2","writexl","foreach","ieugwasr","tibble","plinkbinr","readxl"))
library(data.table)
library(TwoSampleMR)
library(ggplot2)
library(writexl)
library(foreach)
library(ieugwasr)
library(tibble)
library(plinkbinr)
library(readxl)

# -------------------------- 2. Configure Core Path Parameters (EDIT HERE!) --------------------------
# Base directory for all files
base_path <- "D:/BaiduSyncdisk"
# Path to raw exposure data (GWAS summary stats)
exp_raw_dir <- paste0(base_path, "/Data storage/Inflammatory factor/Inflammatory_unfiltered/")
# Output directory for filtered exposure data
exp_filtered_dir <- paste0(base_path, "/Dandelion/exp_dat/")
# Path to outcome GWAS summary stats
outcome_file <- paste0(base_path, "/Dandelion/finngen_R11_I9_PHLETHROMBDVTLOW.gz")
# Output directory for MR results
mr_result_dir <- paste0(base_path, "/Dandelion/result_filter/")
# Output path for final summary Excel
summary_excel_path <- paste0(base_path, "/Dandelion/MR_summary.xlsx")
# Path to ID file (txt with exposure GWAS IDs)
id_file <- paste0(base_path, "/Dandelion/id.txt")
# PLINK reference panel (1000 Genomes EUR) for LD clumping
plink_bfile <- paste0(base_path, "/Data storage/plink/1kg/EUR/EUR")

# -------------------------- 3. Step 1: Exposure Data Filtering (P-value + LD Clumping) --------------------------
# Create directory for filtered exposure data (ignore if exists)
dir.create(exp_filtered_dir, showWarnings = FALSE)

# Read exposure GWAS ID list
id_if <- read.table(id_file, header = TRUE, sep = "\t")
exp_ids <- as.vector(id_if$id)

# Loop through each exposure GWAS for filtering
foreach(i = exp_ids, .errorhandling = "pass") %do% {
  # Read raw exposure data
  exp_dat <- read_exposure_data(
    filename = paste0(exp_raw_dir, i, ".tsv.gz"),
    sep = "\t",
    chr_col = "chromosome",
    pos_col = "base_pair_location",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    beta_col = "beta",
    se_col = "standard_error",
    eaf_col = "effect_allele_frequency",
    pval_col = "p_value",  # Fix: remove extra space in original code
    snp_col = "rsid",
    samplesize_col = "n"
  )
  
  # 1. Filter SNPs by P-value (p < 1e-5)
  exp_dat <- exp_dat[exp_dat$pval.exposure < 1e-5, ]
  
  # 2. LD clumping to remove correlated SNPs
  exp_dat_clump <- ieugwasr::ld_clump(
    dplyr::tibble(
      rsid = exp_dat$SNP,
      pval = exp_dat$pval.exposure,
      id = exp_dat$id.exposure
    ),
    clump_kb = 10000,    # LD clumping window (10kb)
    clump_r2 = 0.001,    # R2 threshold for LD clumping
    clump_p = 1,
    plink_bin = plinkbinr::get_plink_exe(),
    bfile = plink_bfile,
    pop = "EUR"          # Population: European
  )
  
  # Retain clumped SNPs
  exp_dat <- subset(exp_dat, SNP %in% exp_dat_clump$rsid) 
  
  # Export filtered exposure data
  write.table(
    exp_dat,
    paste0(exp_filtered_dir, i, ".txt"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  cat(paste0("Exposure data ", i, " filtered successfully!\n"))
}

# -------------------------- 4. Step 2: Core MR Analysis + Sensitivity Analyses --------------------------
# Create directory for MR results
dir.create(mr_result_dir, showWarnings = FALSE)

# Re-read exposure ID list (ensure consistency)
id_if <- read.table(id_file, header = TRUE, sep = "\t")
exp_ids <- as.vector(id_if$id)

# Empty dataframe for final result collection
result <- data.frame()

# Loop through each exposure for MR analysis
foreach(i = exp_ids, .errorhandling = "pass") %do% {
  # Read filtered exposure data
  exp_dat <- read_exposure_data(
    filename = paste0(exp_filtered_dir, i, ".txt"),
    sep = "\t",
    chr_col = "chromosome",
    pos_col = "base_pair_location",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    beta_col = "beta",
    se_col = "standard_error",
    eaf_col = "effect_allele_frequency",
    pval_col = "p_value",  # Fix: remove extra space in original code
    snp_col = "rsid",
    samplesize_col = "n"
  )
  
  # Secondary filter: stricter P-value threshold (p < 5e-6)
  exp_dat <- exp_dat[exp_dat$pval.exposure < 5e-6, ]
  
  # Read outcome data
  out_dat <- read_outcome_data(
    snps = exp_dat$SNP,
    filename = outcome_file,
    sep = "\t",
    snp_col = "rsids",
    phenotype_col = "phenotype",
    beta_col = "beta",
    se_col = "sebeta",
    eaf_col = "af_alt",
    effect_allele_col = "alt",
    other_allele_col = "ref",
    pval_col = "pval"
  )
  
  # Harmonize exposure and outcome data (allele matching)
  dat <- harmonise_data(exposure_dat = exp_dat, outcome_dat = out_dat)
  
  # Calculate IV strength (F-statistic)
  # Method 1: F-statistic from R2
  dat$R2 <- (2 * (dat$beta.exposure^2) * dat$eaf.exposure * (1 - dat$eaf.exposure)) /
    (2 * (dat$beta.exposure^2) * dat$eaf.exposure * (1 - dat$eaf.exposure) +
       2 * dat$samplesize.exposure * dat$eaf.exposure * (1 - dat$eaf.exposure) * dat$se.exposure^2)
  dat$f_r2 <- dat$R2 * (dat$samplesize.exposure - 2) / (1 - dat$R2)
  
  # Method 2: F-statistic from beta/SE (used in original code)
  dat$f_beta <- (dat$beta.exposure/dat$se.exposure)^2
  
  # Filter SNPs with F-statistic > 10 (avoid weak instrument bias)
  dat <- dat[dat$f_beta > 10, ]
  
  # Core MR analysis
  mr_res <- mr(dat)                     # All MR methods
  mr_res_or <- generate_odds_ratios(mr_res)  # Convert to odds ratios (OR)
  
  # MR analysis with restricted methods (IVW + MR-Egger)
  mr_res_2method <- mr(dat, method_list = c("mr_ivw", "mr_egger_regression"))
  
  # Sensitivity analyses
  hetero_res <- mr_heterogeneity(dat)   # Heterogeneity test
  pleio_res <- mr_pleiotropy_test(dat)  # Horizontal pleiotropy test
  single_snp_res <- mr_singlesnp(dat)   # Single SNP MR analysis
  loo_res <- mr_leaveoneout(dat)        # Leave-one-out analysis
  
  # Export results to Excel (multiple sheets)
  write_xlsx(
    list(
      Sheet1 = exp_dat,          # Filtered exposure data
      Sheet2 = dat,              # Harmonized data
      Sheet3 = mr_res,           # All MR method results
      Sheet4 = mr_res_or,        # MR results (OR values)
      Sheet5 = mr_res_2method,   # IVW + MR-Egger results
      Sheet6 = hetero_res,       # Heterogeneity test results
      Sheet7 = pleio_res,        # Horizontal pleiotropy test results
      Sheet8 = single_snp_res,   # Single SNP MR results
      Sheet9 = loo_res           # Leave-one-out analysis results (fix case consistency)
    ),
    path = paste0(mr_result_dir, i, ".outcome.xlsx")
  )
  cat(paste0("MR analysis for ", i, " completed!\n"))
}

# Export empty result file (retained from original code)
write.table(result, paste0(mr_result_dir, "result.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

# -------------------------- 5. Step 3: Summarize MR Results --------------------------
# Define folder path for result Excel files
folder_path <- paste0(base_path, "/Dandelion/result_excel/")
# List all Excel files in the folder
files <- list.files(path = folder_path, pattern = "\\.xlsx$", full.names = FALSE)
# Empty dataframe for summary
all_summary_rows <- NULL

# Extract row 3 from Sheet4 (OR results) for each file
for (file in files) {
  file_path <- file.path(folder_path, file)
  tryCatch({
    # Read Sheet4 (MR results with OR values)
    excel_data <- read_excel(file_path, sheet = "Sheet4")
    # Extract 3rd row (target row from original code)
    target_row <- excel_data[3, ]
    # Add filename for traceability
    target_row$FileName <- file
    
    # Merge to summary dataframe
    if (is.null(all_summary_rows)) {
      all_summary_rows <- target_row
    } else {
      all_summary_rows <- rbind(all_summary_rows, target_row)
    }
  }, error = function(e) {
    # Capture and print errors (avoid loop interruption)
    cat(sprintf("Error reading %s: %s\n", file_path, conditionMessage(e)))
  })
}

# Export summary results to Excel
write_xlsx(all_summary_rows, summary_excel_path)
cat(paste0("All results summarized! Summary file: ", summary_excel_path, "\n"))