#!/bin/bash

#SBATCH -J minimap2
#SBATCH -c 8
#SBATCH --mem 64G
#SBATCH -t 48:00:00
#SBATCH --partition=ncpu
#SBATCH -o ./plasmid/minimap2_%A.out


# cat /camp/lab/ulej/home/users/farawar/GASR/export_reporter/nanopore/20220429/fastq/cat_fastq/both_runs.fastq.gz /camp/lab/ulej/home/users/farawar/GASR/export_reporter/nanopore/20230426_oscarneve_plasmid/both_basecalls_merged.fastq.gz > /camp/home/farawar/home/GASR/export_reporter/nanopore/minimap2/plasmid/merged_fastq_bothruns.fq.gz
~/home/General/minimap2/minimap2 -ax splice all_reporter_sequences_with_introns.fa /camp/home/farawar/home/GASR/export_reporter/nanopore/minimap2/plasmid/merged_fastq_bothruns.fq.gz | samtools sort -@ 8 - > ./plasmid/plasmid_alignment.bam