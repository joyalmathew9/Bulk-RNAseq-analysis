Bulk RNA-seq Analysis
------------------------

This repository contains R scripts for bulk RNA-seq analysis using Salmon quantification files and DESeq2.

Workflow
--------

1.  Import Salmon quantification files using `tximeta`

2.  Create a DESeq2 dataset

3.  Filter low-count genes

4.  Perform rlog transformation

5.  Generate PCA plot

6.  Run differential expression analysis

7.  Export significant, upregulated and downregulated genes

8.  Generate volcano plots

9.  Generate heatmap of selected differentially expressed genes

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
├── 01_deseq2_analysis.R
├── 02_volcano_plot.R
└── 03_heatmap_visualization.R
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

