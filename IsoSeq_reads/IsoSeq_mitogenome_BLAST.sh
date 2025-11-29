#!/bin/bash


DB="assembly_I5-EDIT_Final_cox1_start.fas"

for file in *_trimmed_polyT.fastq; do
    

    base="${file%.fastq}"
    fasta="${base}.fas"
    blast_out="${base}_blast.out"

    # Convert fastq to fasta using paste
    paste - - - - < "$file" \
        | cut -f 1,2 \
        | sed 's/^@/>/' \
        | tr "\t" "\n" \
        > "$fasta"

    blastn -dust no -evalue 0.000001 \
        -qcov_hsp_perc 80 \
        -perc_identity 95 \
        -outfmt "6 qseqid sseqid pident length evalue bitscore qcovs qcovhsp" \
        -query "$fasta" \
        -db "$DB" \
        -out "$blast_out"

done

