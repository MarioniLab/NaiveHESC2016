# Annotation sources

- ERCC92.gtf was taken from http://www.thermofisher.com/order/catalog/product/4456739
- hg38.gtf was constructed from http://ftp.ensembl.org/pub/release-91/gtf/homo_sapiens/

```sh
zcat Homo_sapiens.GRCh38.91.gtf.gz | \
    sed -r "s/^([0-9MXY])/chr\1/" | \
    sed "s/^chrMT/chrM/g" | \
    awk '$3 == "exon"' > hg38.gtf
```
