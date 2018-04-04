#!/bin/bash
#SBATCH -o count.out
#SBATCH -e count.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16000

echo '
anno.files <- file.path("../../genomes/annotation/", c("hg38.gtf", "ERCC92.gtf"))
bam.files <- list.files("bam", full=TRUE, pattern="bam$")
stat.file <- "all_qual.tsv"
ispet <- TRUE
# NOT STRAND-SPECIFIC.
' | cat - ../../tools/counter.R > count_me.R

R CMD BATCH --no-save count_me.R 
