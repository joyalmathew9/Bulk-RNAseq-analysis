Bulk RNA-seq Analysis
------------------------

This repository contains R scripts for bulk RNA-seq analysis using Salmon quantification files and DESeq2.

Workflow
--------

1.  Build a Salmon transcriptome index

2.  Quantify transcript abundance using Salmon

3.  Import Salmon quantification files using `tximeta`

4.  Create a DESeq2 dataset

5.  Filter low-count genes

6.  Perform rlog transformation

7.  Generate PCA plot

8.  Run differential expression analysis

9.  Export significant, upregulated and downregulated genes

10.  Generate volcano plots

11.  Generate heatmap of selected differentially expressed genes

Input files
-----------

Expected files:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ text
samples.txt
salmon/
  sample_1/quant.sf
  sample_2/quant.sf
salmon_tx2gene.csv
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The `samples.txt` file should contain at least:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ text
sample_id    condition
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Scripts
-------

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ text
scripts/
├── 01_salmon_index.sh
├── 02_salmon_quant.sh
├── 03_deseq2_analysis.R
├── 04_volcano_plot.R
└── 05_heatmap_visualization.R
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Outputs
-------

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ text

plots/
├── PCA.png
├── Volcano_PBS_vs_LNP.png
├── heatmap1
└── heatmap2.png
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

