# Function to install/load packages (CRAN + Bioconductor)
install_load_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% c("ComplexHeatmap")) { # Bioconductor packages
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(pkg, update = FALSE)
    } else { # CRAN packages
      install.packages(pkg, repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
    }
  }
  library(pkg, character.only = TRUE)
}

# Load required packages
core_pkgs <- c("readxl", "writexl", "dplyr", "circlize", "ComplexHeatmap", 
               "RColorBrewer", "dendextend", "dendsort", "gridBase", 
               "grid", "forestploter")
invisible(lapply(core_pkgs, install_load_pkg))

# ===================== 2. Configure Global Paths (Unified Management) =====================
# Path for raw Excel files (row data extraction)
data_extract_path <- "D:\\BaiduSyncdisk\\Dandelion\\data_test"
# Path for plotting data (heatmap + forest plot)
plot_data_path <- "D:\\BaiduNetdiskDownload\\作图测试"
# Output Excel filename
output_excel_name <- "Heatmap_Data_Untidy.xlsx"

# ===================== 3. Batch Extract Rows 1-5 from Sheet3 of Excel Files =====================
# Get all xlsx files in target directory
files <- list.files(path = data_extract_path, pattern = "\\.xlsx$", full.names = FALSE)

# Initialize storage list (map row numbers to methods)
row_data_list <- list(
  ME = NULL,    # Row 1: ME method
  WM = NULL,    # Row 2: WM method
  IVW = NULL,   # Row 3: IVW method
  SM = NULL,    # Row 4: SM method
  WMODE = NULL  # Row 5: WMODE method
)

# Single loop to extract rows 1-5 (replace 5 redundant loops)
for (file in files) {
  file_fullpath <- file.path(data_extract_path, file)
  tryCatch({
    # Read Sheet3 data from Excel file
    excel_data <- read_excel(file_fullpath, sheet = "Sheet3")
    
    # Extract rows 1-5 and append filename column
    for (row_idx in 1:5) {
      current_row <- excel_data[row_idx, ]
      current_row$FileName <- file
      
      # Assign to corresponding list (ME=row1/WM=row2/IVW=row3/SM=row4/WMODE=row5)
      list_name <- names(row_data_list)[row_idx]
      if (is.null(row_data_list[[list_name]])) {
        row_data_list[[list_name]] <- current_row
      } else {
        row_data_list[[list_name]] <- rbind(row_data_list[[list_name]], current_row)
      }
    }
  }, error = function(e) {
    # Error message for failed file reading (corrected Sheet name)
    cat(sprintf("Failed to read file %s (Sheet3): %s\n", file_fullpath, conditionMessage(e)))
  })
}

# Export aggregated data to Excel (separate sheets for different methods)
write_xlsx(
  list(Sheet1 = row_data_list$ME, Sheet2 = row_data_list$WM, 
       Sheet3 = row_data_list$IVW, Sheet4 = row_data_list$SM, 
       Sheet5 = row_data_list$WMODE),
  file.path(plot_data_path, output_excel_name)
)

# ===================== 4. Plot Circular Heatmap =====================
# Read cleaned heatmap data
heatmap_data <- read.table(
  file.path(plot_data_path, "Heatmap_Data_Cleaned.csv"), # Rename to English if possible
  header = TRUE, row.names = 1, sep = ","
)

# Convert to matrix and normalize (row-wise scaling)
heatmap_mat <- as.matrix(heatmap_data)
heatmap_mat_norm <- t(scale(t(heatmap_mat))) # Row-wise normalization for machine learning
cat("Normalized data range: ", range(heatmap_mat_norm), "\n")

# Define heatmap color gradient (customizable)
color_gradient <- colorRamp2(c(-1.8, 0, 1.8), c("blue", "white", "red"))
# Alternative color scheme: color_gradient <- colorRamp2(c(-2, 0, 2), c("#003399", "white", "#cccc00"))

# Plot circular heatmap
circos.par(gap.after = 30) # Adjust gap between start/end of the circle
circos.heatmap(
  heatmap_mat_norm,
  col = color_gradient,
  dend.side = "inside",        # Dendrogram on inner side of the circle
  rownames.side = "outside",   # Row names on outer side (opposite to dendrogram)
  track.height = 0.2,          # Thickness of the circular track
  rownames.col = "black",      # Color of row names
  rownames.cex = 0.6,          # Font size of row names
  cluster = TRUE,              # Enable row clustering
  dend.track.height = 0.08,    # Height of the dendrogram track
  dend.callback = function(dend, m, si) {
    # Color dendrogram branches (k=15 = 15 clusters, colors 1:15)
    color_branches(dend, k = 15, col = 1:15)
  }
)
circos.clear() # Mandatory: avoid overlay in subsequent plots

# Add legend
legend <- Legend(
  title = "Exp", col_fun = color_gradient,
  direction = "vertical"
)
draw(legend, x = unit(0.75, "npc"), y = unit(0.7, "npc"), just = c("right", "center"))

# ===================== 5. Plot Forest Plot =====================
# Read forest plot data
forest_data <- read_excel(file.path(plot_data_path, "abbb.xlsx"))

# Data preprocessing (unified logic, reduce redundancy)
forest_data_processed <- forest_data %>%
  # Insert blank column (separator column)
  mutate(` ` = paste(rep(" ", 20), collapse = " ")) %>%
  relocate(` `, .after = 5) %>%
  # Retain 4 decimal places (specified columns)
  mutate(across(c(5, 10:15), ~ sprintf("%.4f", as.numeric(.x)))) %>%
  # Process Egger_intercept format (5 significant digits, no scientific notation)
  mutate(Egger_intercept = format(as.numeric(Egger_intercept), digits = 5, nsmall = 4, scientific = FALSE)) %>%
  # Replace NA with blank spaces
  mutate(Exposure = ifelse(is.na(Exposure), "  ", as.character(Exposure))) %>%
  mutate(across(c(Q, Q_pval, Egger_intercept, Plelotropy_pval, Test.RSSobs, Test.Pvalue),
                ~ ifelse(is.na(as.numeric(.x)), "  ", as.numeric(.x)))) %>%
  # Add OR(95%CI) column
  mutate(`OR(95%CI)` = ifelse(is.na(or), "",
                              sprintf('%.3f(%.3f to %.3f)', or, or_lci95, or_uci95))) %>%
  relocate(`OR(95%CI)`, .after = 9)

# Define forest plot theme
forest_theme <- forest_theme(
  base_size = 10,
  ci_pch = 20,                  # Symbol for confidence interval point
  ci_col = "#4575b4",           # Color of confidence interval
  ci_lty = 1,                   # Line type of confidence interval
  ci_lwd = 2.3,                 # Line width of confidence interval
  ci_Theight = 0.2,             # Height of confidence interval tick
  refline_lwd = 1.5,            # Line width of reference line (OR=1)
  refline_lty = "dashed",       # Line type of reference line
  refline_col = "red",          # Color of reference line
  summary_fill = "#4575b4",     # Fill color of summary point
  summary_col = "#4575b4",      # Border color of summary point
  footnote_cex = 1.1,           # Font size of footnote
  footnote_fontface = "italic", # Font style of footnote
  footnote_col = "blue"         # Color of footnote
)

# Plot forest plot
forest_plot <- forest(
  forest_data_processed[, c(1:6, 10:16)],
  est = forest_data_processed$or,        # Effect size (OR)
  lower = forest_data_processed$or_lci95, # 95%CI lower bound
  upper = forest_data_processed$or_uci95, # 95%CI upper bound
  sizes = 0.6,                           # Size of CI points
  ci_column = 6,                         # Column position for confidence interval
  ref_line = 1,                          # Reference line (OR=1)
  xlim = c(0.5, 1.8),                    # X-axis range
  ticks_at = c(0.5, 0.75, 1, 1.2, 1.4, 1.6, 1.8), # X-axis ticks
  arrow_lab = c('protective factor', 'risk factor'), # Arrow labels
  footnote = 'P<0.05 was considered statistically significant', # Footnote
  theme = forest_theme
)

# Display forest plot
print(forest_plot)