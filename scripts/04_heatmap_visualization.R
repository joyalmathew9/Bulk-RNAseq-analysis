# Heatmap visualization of differentially expressed genes from DESeq2 output

# ---- Load libraries ----
library(DESeq2)
library(genefilter)
library(ComplexHeatmap)
library(grid)

# ---- User settings ----
project_dir <- getwd()

results_dir <- file.path(project_dir, "results")
plots_dir <- file.path(project_dir, "plots")

rlog_file <- file.path(results_dir, "rlog_object.rds")
filtered_results_file <- file.path(results_dir, "filtered_results_PBS_vs_LNP.rds")

conditions_to_plot <- c("PBS", "LNP")
top_n_genes <- 100

dir.create(plots_dir, showWarnings = FALSE)

# ---- Load saved objects ----
rld <- readRDS(rlog_file)
filtered_results <- readRDS(filtered_results_file)

# ---- Subset samples by condition ----
rld_subset <- rld[, rld$condition %in% conditions_to_plot]

normalized_counts <- assay(rld_subset)

# ---- Select genes present in both expression matrix and DE results ----
common_genes <- intersect(rownames(normalized_counts), filtered_results$gene)

filtered_counts <- normalized_counts[common_genes, ]

# ---- Select top variable genes ----
top_var_genes <- head(order(rowVars(filtered_counts), decreasing = TRUE), top_n_genes)

mat <- filtered_counts[top_var_genes, ]

# Center expression values by gene
mat <- mat - rowMeans(mat)

# ---- Sample annotation ----
sample_annotation <- HeatmapAnnotation(
  condition = colData(rld_subset)$condition
)

# ---- Generate heatmap ----
heatmap_object <- Heatmap(
  mat,
  name = "Expression",
  col = colorRampPalette(c("blue", "white", "red"))(100),
  show_row_names = TRUE,
  show_column_names = TRUE,
  top_annotation = sample_annotation,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  row_names_gp = grid::gpar(fontsize = 7)
)

png(
  filename = file.path(plots_dir, "heatmap_top_variable_DE_genes.png"),
  width = 1800,
  height = 1400,
  res = 200
)

draw(
  heatmap_object,
  heatmap_legend_side = "right",
  column_title = paste(
    paste(conditions_to_plot, collapse = " vs "),
    "- top",
    top_n_genes,
    "differentially expressed genes"
  )
)

dev.off()
