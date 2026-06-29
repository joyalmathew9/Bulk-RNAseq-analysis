#!/bin/bash

# ----------------------------------------------------------
# Quantify transcript abundance using Salmon
# ----------------------------------------------------------

# User settings
INDEX_DIR="salmon_index"

READ1="sample_R1.fastq.gz"     # sequencing reads (FASTQ - Forward reads)
READ2="sample_R2.fastq.gz"     # sequencing reads (FASTQ - Reverse reads)

OUTPUT_DIR="sample_quant"

# Quantify transcript abundance
salmon quant \
    -i "$INDEX_DIR" \
    -l A \
    -1 "$READ1" \
    -2 "$READ2" \
    -p 8 \
    -o "$OUTPUT_DIR"

echo "Quantification completed."
echo "Output file: $OUTPUT_DIR/quant.sf"

# Inside the output directory you will find:
# quant.sf - transcript-level abundance estimates