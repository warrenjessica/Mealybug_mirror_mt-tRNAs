#!/bin/bash
for fasta in *collapsed.fasta; 
do 
    awk 'BEGIN{RS=">"}{print $1"\t"$2;}' $fasta | tail -n+2 > ${fasta%.collapsed.fasta}_collapsed.tab
done;
