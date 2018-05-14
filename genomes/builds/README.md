# Sequence sources

- ERCC sequences were obtained from https://www.thermofisher.com/order/catalog/product/4456739
- hg38 sequence was obtained from UCSC (via `/scratchb/bioinformatics/reference_data/reference_genomes/homo_sapiens/hg38/fasta/hsa.hg38.fa`)

# Genome builds

This combines the hg38 build of the human genome with the ERCC sequences:

```sh
subread-buildindex -o hg38_ERCC /scratchb/bioinformatics/reference_data/reference_genomes/homo_sapiens/hg38/fasta/hsa.hg38.fa ../sequences/spikes/ERCC92.fa
```
