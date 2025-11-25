#!/bin/bash

for fastq in *_1.fastq; 
do 
    # Define the reverse read based on the forward read
    reverse_fastq=${fastq%_1.fastq}_2.fastq
    
    # Run cutadapt only for the pair of files
    cutadapt -q 10 --discard-untrimmed --nextseq-trim=20 -m 25 -M 95 \
    -a ACTGGATACTGGN...GTATCCAGTTGGAATTCTCGGGTGCCAAGG \
    -A CTGGATAC...NCCAGTATCCAGTGATCGTCGGACTGTAGAACTCTGAAC \
    -o ${fastq%_1.fastq}_1.trimmed.fq \
    -p ${fastq%_1.fastq}_2.trimmed.fq \
    $fastq $reverse_fastq > ${fastq%_1.fastq}_trimstats.txt
done;


for fq in *_1.trimmed.fq; 
do 
    # Define the reverse read based on the forward read
    reverse_fq=${fq%_1.trimmed.fq}_2.trimmed.fq
    
    # Define the output merged file and the other output files
    output_prefix=${fq%_1.trimmed.fq}
    
    # Run bbmerge with the paired reads
    bbmerge.sh in1=$fq in2=$reverse_fq out=${output_prefix}_merged.fq \
    ihist=${output_prefix}_ihist.txt ordered=t minoverlap=20 mismatches=0 \
    2> ${output_prefix}_bbmerge_log.txt
done;

#Count CCA % for library, collapse reads with read count, rename 

for merged in *merged.fq; 

		 do perl count_CCA.pl $merged > ${merged%_merged.fq}_count_CCA.text
		 
		 fastx_collapser -v -i $merged -o ${merged%merged.fq}collapsed.fasta > ${merged%merged.fq}_collapsestats.txt
		  
		 awk '/>/{sub(">","&"FILENAME"_");sub(/\.fasta/,x)}1' ${merged%merged.fq}collapsed.fasta > ${merged%merged.fq}collapsed.ID.fasta
		 
		 perl collapsed_fasta_to_freq_dist.pl ${merged%merged.fq}collapsed.ID.fasta > ${merged%merged.fq}readcount_dist.txt
		 
done; 



