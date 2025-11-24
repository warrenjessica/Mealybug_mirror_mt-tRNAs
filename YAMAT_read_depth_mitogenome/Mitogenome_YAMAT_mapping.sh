#!/usr/bin/sh

# fastq files for the 3 mitochondrial isolation reps ( ) were converted to fasta files and mapped to the mitogenome using BLAST

# only YAMAT reads with 2 or more reads were used for this analysis 

for f in *.fasta; do
    python filter_fasta.py "$f" > "${f%.fasta}_above1.CCA.fasta"
done

# filtered reads then had the CCA tail off of the reads removed prior to BLASTing using trim_fasta.py on each of the files 

for f in *_above1.CCA.fasta; do
    out="${f%.fasta}_trimmed.fasta"
    echo "Trimming CCA tail: $f -> $out"
    python trim_fasta.py "$f" > "$out"
done


# use BLAST to find the coordinates of each hit on the mitogenome 

for file in *_above1.CCA.fasta; do
    base="${file%.fasta}"
    output="${base}_blast_above1_CCA_fmt6.txt"
    blastn -task blastn -dust no -evalue 0.000001 \
        -outfmt "6 qseqid sstrand qlen pident length qstart qend sstart send evalue bitscore qcovhsp" \
        -query "$file" \
        -db assembly-I5_edit_final.fa \
        -out "$output"
    
    echo "Finished BLAST for $file -> $output"
done

#The BLAST output file was then loaded into R and reads were filtered so that only reads with a >95% percent identity and >95% hit coverage were mapped to the mitogenome to remove any nuclear (including NUMT) sequences. 
# see Mitogenome_YAMAT_mapping.R script for filtering work. 