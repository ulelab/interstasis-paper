#!/bin/bash

#SBATCH -J minimap2
#SBATCH -c 8
#SBATCH --mem 64G
#SBATCH -t 48:00:00
#SBATCH --partition=ncpu
#SBATCH -o ./minimap2_%A.out

~/home/General/minimap2/minimap2 -ax splice all_reporter_sequences_with_introns.fa /camp/lab/ulej/home/users/farawar/GASR/export_reporter/nanopore/20220518/fastq/demulitplexed/rupert_all.fastq.gz > alignment.sam