#!/bin/bash



for file in *collapsed.ID.fasta;
	do blastn -task blastn -dust no -query ${file} -db pcitri_db_update3.fas -out ${file%_collapsed.ID.fasta}.collapsed.ID.BLAST.txt
	
	
done;