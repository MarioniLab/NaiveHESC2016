# Naive and primed pluripotency in human ESCs

## Overview

This repository provides the code for the paper **Transcriptional heterogeneity in naive and primed human pluripotent stem cells at single-cell resolution**,
by [Messmer _et al._ (2019)](https://doi.org/10.1016/j.celrep.2018.12.099).

## Generating the results

To reproduce the analyses from count matrices:

1. Enter `data/` and run `download.sh` to download all relevant count tables.
2. Enter `analysis/` and run `run_main.sh` to compile all of the analysis scripts.
3. Enter `figures/` and run `make_figures.sh` to create all of the figures used in the manuscript.

## Realigning the sequence data

To reproduce the analyses from the FASTQ files:

1. Follow the instructions in `genomes/builds/README.md` to build the genome indices.
Similarly, follow the instructions in `genomes/annotation/README.md` to obtain the annotation.
2. Download the FASTQ files from [ArrayExpress](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-6819/).
Files corresponding to each batch of data should be placed in `data/real-2383_20161024/fastq`, etc.
3. Run the various `mapme.sh` scripts to execute the master scripts for alignment, and `count_me.sh` for read counting.
Paths should refer to the top-level `tools/` directory obtained using `download.sh`.
