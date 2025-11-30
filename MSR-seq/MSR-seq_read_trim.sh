#!/bin/sh

#merge paired-end reads
for file in *1.fq.gz; do bbmerge.sh in1=$file in2=${file%1.fq.gz}2.fq.gz out=${file%_CKDL250008135-1A_22TG3WLT4_L3_1.fq.gz}.bbmerge.fq outu=${file%_CKDL250008135-1A_22TG3WLT4_L3_1.fq.gz}.unmerged1.fq outu2=${file%_CKDL250008135-1A_22TG3WLT4_L3_1.fq.gz}.unmerged2.fq ihist=${file%_CKDL250008135-1A_22TG3WLT4_L3_1.fq.gz}.bbmerge_hist.txt ordered=t qtrim=r minoverlap=30 mismatches=4 2> ${file%_CKDL250008135-1A_22TG3WLT4_L3_1.fq.gz}.bbmerge_log.txt; done

#trim MSR-seq adapter sequence
for file in *bbmerge.fq; do perl MSR-seq_trim.pl $file ACTGGAA 6 > ${file%fq}trim.fas; done

#generate a file that collapses all identical read sequences into a single sequence
for file in *trim.fas; do perl collapse_identical_seqs.pl $file > ${file%trim.fas}collapsed.fas; done



