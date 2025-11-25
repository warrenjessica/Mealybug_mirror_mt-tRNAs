#!/bin/bash

# Adapter sequences
ADAPTER_5P="AAGCAGTGGTATCAACGCAGAGTAC"
ADAPTER_3P="GTACTCTGCGTTGATACCACTGCTT"

# Loop over all FASTQ files in the directory
for file in *.fastq; do
    base="${file%.fastq}"
    echo "Processing $file ..."

    trimmed="${base}_trimmed.fastq"
    summary="${base}_cutadapt_summary.txt"
    final="${base}_trimmed_polyT.fastq"   

    # First: adapter trimming with cutadapt
    cutadapt \
      -g "$ADAPTER_5P" \
      -a "$ADAPTER_3P" \
      --poly-a \
      --overlap 10 \
      -e 0.2 \
      --match-read-wildcards \
      -m 100 \
      -j 8 \
      -o "$trimmed" \
      "$file" \
      > "$summary" 2>&1


    python trim_polyT_sliding.py "$trimmed" > "$final"

done



