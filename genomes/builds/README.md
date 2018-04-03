# Sequence sources

- ERCC sequences were obtained from https://www.thermofisher.com/order/catalog/product/4456739
- mm10 sequence was obtained from UCSC (via `/scratchb/bioinformatics/reference_data/reference_genomes/mus_musculus/mm10/fasta/mmu.mm10.fa`)

# Genome builds

This combines the mm10 build of the mouse genome with the ERCC sequences:

```sh
subread-buildindex -o mm10_ERCC /scratchb/bioinformatics/reference_data/reference_genomes/mus_musculus/mm10/fasta/mmu.mm10.fa ERCC92.fa
```
