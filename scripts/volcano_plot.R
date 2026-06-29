# ----------------------------------------------------------
# Volcano Plot
# ----------------------------------------------------------

# ---- Load libraries ----
library(EnhancedVolcano)

# ---- User settings ----
project_dir <- getwd()

results_dir <- file.path(project_dir, "results")
plots_dir   <- file.path(project_dir, "plots")

comparison_name <- "PBS_vs_LNP"

results_file <- file.path(
  results_dir,
  paste0("DESeq2_results_", comparison_name, ".csv")
)

dir.create(plots_dir, showWarnings = FALSE)

# ---- Load DESeq2 results ----
res_df <- read.csv(
  results_file,
  row.names = 1,
  check.names = FALSE
)

# Remove genes with missing values
res_df <- res_df[!is.na(res_df$pvalue), ]

# ---- Generate volcano plot ----
volcano_plot <- EnhancedVolcano(
  res_df,
  lab = rownames(res_df),
  x = "log2FoldChange",
  y = "pvalue",
  
  pCutoff = 0.05,
  FCcutoff = 1,
  
  drawConnectors = TRUE,
  max.overlaps = 20,
  
  title = paste("Volcano Plot -", gsub("_", " ", comparison_name)),
  subtitle = NULL
)

# ---- Save figure ----
ggplot2::ggsave(
  filename = file.path(
    plots_dir,
    paste0("Volcano_", comparison_name, ".png")
  ),
  plot = volcano_plot,
  width = 8,
  height = 6,
  dpi = 300
)