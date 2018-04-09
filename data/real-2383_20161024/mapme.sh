ispet=1
fastq=$(ls fastq/* | grep ".fastq.gz")
genome=../../genomes/builds/hg38_ERCC
source ../../tools/multi_align.sh
