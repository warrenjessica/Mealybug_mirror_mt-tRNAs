#!/usr/bin/sh

# Porechop (v. 0.2.4) to trim Nanopore reads
porechop -i TotalMito_nanopore.fastq -o TotalMito_nanopore_trimmed.fastq --discard_middle

# The TotalMito_nanopore_trimmed.fastq file was converted to TotalMito_nanopore_trimmed.fas
paste - - - - < TotalMito_nanopore_trimmed.fastq | cut -f 1,2 | sed 's/^@/>/' | tr "\t" "\n" > TotalMito_nanopore_trimmed.fas

# Using BLAST to extract reads that map to the previously assembled Wellcome Sanger Tree of Life Programme (GenBank accession OX465514)
blastn -dust no  -evalue 0.0000001 -outfmt "6 qseqid" -query TotalMito_nanopore_trimmed.fas -db OX465514.fasta -out TotalMito_nanopore_trimmed_OX465514_fmt6_BLAST.txt

sort TotalMito_nanopore_trimmed_OX465514_fmt6_BLAST.txt | uniq > unique_TotalMito_nanopore_trimmed_OX465514_fmt6_BLAST.txt

python extract_fastq.py TotalMito_nanopore_trimmed.fas unique_TotalMito_nanopore_trimmed_OX465514_fmt6_BLAST.txt extracted_unique_TotalMito_nanopore_trimmed_OX465514_fmt6_BLAST.fastq

#flye assembly
flye --nano-raw extracted_unique_TotalMito_nanopore_trimmed_OX465514_fmt6_BLAST.fastq --out-dir flye_assembly_I5 --iterations 5

# renamed assembly-I5.fasta to renamed to mito.fa

# polishing nanopore assembly with Illumina reads. Illumina reads were trimmed with cutadapt.

bowtie2-build mito.fa mito.fa

cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -q 20 --minimum-length 50 -e 0.15 -j 24 -o Total_Mito.R1.trim.fq -p Total_Mito.R2.trim.fq Total_Mito_S231_R1_001.fastq.gz Total_Mito_S231_R2_001.fastq.gz  >> cutadapt.log.txt

bowtie2 --no-unal -p 24 -x mito.fa -1 Total_Mito.R1.trim.fq -2 Total_Mito.R2.trim.fq -S Total_Mito.sam >> bowtie2.log.txt 2>> bowtie2.err.txt

samtools sort Total_Mito.sam > Total_Mito.bam

samtools index Total_Mito.bam

conda activate perbase_env

perbase base-depth --max-depth 1000000 --threads 24 Total_Mito.bam > Total_Mito.depth.txt
echo Total_Mito.depth.txt > file_list.txt

#perbase_variant_summary2.pl used to summarize the perbase changes. 

perl perbase_variant_summary2.pl file_list.txt mito.fa 0.1 20 > perbase_summary.txt





