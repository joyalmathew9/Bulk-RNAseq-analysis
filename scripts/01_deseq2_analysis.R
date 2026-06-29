# Bulk RNA-seq differential expression analysis using Salmon quantification files and DESeq2

# ---- Load libraries ----
library(readr)
library(tximport)
library(tximeta)
library(DESeq2)
library(ggplot2)

# ---- User settings ----
project_dir <- getwd()

samples_file <- file.path(project_dir, "samples.txt")
salmon_dir <- file.path(project_dir, "salmon")
tx2gene_file <- file.path(project_dir, "salmon_tx2gene.csv")

reference_condition <- "PBS"
condition_levels <- c("PBS", "LNP", "SN", "N", "S")

contrast_condition <- c("condition", "PBS", "LNP")

min_gene_count <- 10
pvalue_cutoff <- 0.05
log2fc_cutoff <- 1

results_dir <- file.path(project_dir, "results")
plots_dir <- file.path(project_dir, "plots")

dir.create(results_dir, showWarnings = FALSE)
dir.create(plots_dir, showWarnings = FALSE)

# ---- Load sample metadata and Salmon quantification files ----
samples <- read.table(samples_file, header = TRUE)

files <- file.path(salmon_dir, samples$sample_id, "quant.sf")
names(files) <- samples$sample_id
rownames(samples) <- samples$sample_id

# ---- Import transcript-level quantifications ----
coldata <- samples
coldata$files <- files
coldata$names <- coldata$sample_id

tx2gene <- read_delim(tx2gene_file, delim = "\t", show_col_types = FALSE)

se <- tximeta(
  coldata,
  skipMeta = TRUE,
  txOut = FALSE,
  tx2gene = tx2gene
)

# ---- Create DESeq2 dataset ----
dds <- DESeqDataSet(se, design = ~ condition)

dds$condition <- factor(dds$condition, levels = condition_levels)

# ---- Filter low-count genes ----
keep <- rowSums(counts(dds)) >= min_gene_count
dds <- dds[keep, ]

# ---- Transformation and PCA ----
rld <- rlog(dds, blind = FALSE)

pca_all <- plotPCA(rld, intgroup = "condition") +
  geom_text(aes(label = name), vjust = 2) +
  ggtitle("PCA - all samples") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(
  filename = file.path(plots_dir, "PCA_all_samples.png"),
  plot = pca_all,
  width = 7,
  height = 5,
  dpi = 300
)

# ---- Differential expression analysis ----
dds <- DESeq(dds)

res <- results(dds, contrast = contrast_condition, independentFiltering = FALSE)
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

res_df <- res_df[!is.na(res_df$pvalue), ]

comparison_name <- paste(contrast_condition[2], "vs", contrast_condition[3], sep = "_")

write.csv(
  res_df,
  file.path(results_dir, paste0("DESeq2_results_", comparison_name, ".csv")),
  row.names = TRUE
)

# ---- Filter significant genes ----
filtered <- res_df[res_df$pvalue < pvalue_cutoff, ]

upregulated_genes <- filtered[
  filtered$log2FoldChange > log2fc_cutoff,
]

downregulated_genes <- filtered[
  filtered$log2FoldChange < -log2fc_cutoff,
]

write.csv(
  filtered,
  file.path(results_dir, paste0("significant_genes_", comparison_name, ".csv")),
  row.names = FALSE
)

write.csv(
  upregulated_genes,
  file.path(results_dir, paste0("upregulated_genes_", comparison_name, ".csv")),
  row.names = FALSE
)

write.csv(
  downregulated_genes,
  file.path(results_dir, paste0("downregulated_genes_", comparison_name, ".csv")),
  row.names = FALSE
)

# ---- Save normalized expression matrix and objects ----
normalized_counts <- assay(rld)

write.csv(
  normalized_counts,
  file.path(results_dir, "rlog_normalized_counts.csv"),
  row.names = TRUE
)

saveRDS(dds, file.path(results_dir, "dds_object.rds"))
saveRDS(rld, file.path(results_dir, "rlog_object.rds"))
saveRDS(filtered, file.path(results_dir, paste0("filtered_results_", comparison_name, ".rds")))
