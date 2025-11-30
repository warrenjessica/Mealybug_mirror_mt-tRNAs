#!/bin/sh

#SBATCH --time=100000:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dan.sloan@colostate.edu
#SBATCH --error=slurm.stderr
#SBATCH --output=slurm.stdout


for file in *collapsed.fas; do blastn -task blastn -db pcitri_db4.fas -query $file -evalue 1e-12 -num_threads 12 -dust no -out ${file%collapsed.fas}blast.txt; perl trim_5prime_blast.pl ${file%collapsed.fas}blast.txt $file ${file%collapsed.fas}trim.fas > ${file%collapsed.fas}blast_processed.fas; bowtie2 --no-unal -p 12 -L 10 -i C,1 --mp 5,2 --score-min L,-0.7,-0.7 -f -x pcitri_db4.fas -U ${file%collapsed.fas}blast_processed.fas -S ${file%collapsed.fas}sam; perl CC_vs_CCA_counter.sam4.pl ${file%collapsed.fas}sam > ${file%collapsed.fas}CC_vs_CCA.txt; perl CC_vs_CCA_counter.sam4.pl ${file%collapsed.fas}sam --unique > ${file%collapsed.fas}CC_vs_CCA.unique.txt; done
