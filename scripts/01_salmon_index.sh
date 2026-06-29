#!/bin/bash

# ----------------------------------------------------------
# Build Salmon transcriptome index
# ----------------------------------------------------------

# User settings
TRANSCRIPTOME="transcripts.fa"
INDEX_DIR="salmon_index"

# Build Salmon index (run once per reference transcriptome)
salmon index \
    -t "$TRANSCRIPTOME" \
    -i "$INDEX_DIR"

echo "Salmon index created successfully."
