---
title: "hESC methods"
author: Tobias Messmer and Aaron Lun
bibliography: ref.bib
---

## Alignment and read counting

Read pairs were aligned to a reference consisting of the hg38 build of the human genome as well sequences for the ERCC spike-in transcripts. 
This was performed using the `subread` aligner v1.5.0-p3 [@liao2013subread] in paired-end mode with unique alignment.
Each read pair was then assigned to a gene in the Ensembl GRCh38 v83 annotation or to the spike-in transcripts.
This was done using the `featureCounts` function in the _Rsubread_ package v1.24.0 [@liao2014featurecounts].
Only reads with mapping quality scores above 10 were used for counting.
Read counts from technical (sequencing) replicates of the same cell were added together prior to further analysis.

## Quality control on cells and genes 

A range of quality metrics were computed for each cell using the `calculateQCMetrics` function in the _scater_ package v1.2.0 [@mccarthy2016scater].
For each metric, outlier values were identified as those that were more than three median absolute deviations from the median.
Low quality cells were identified in each batch, as those with small outlier values for the log-transformed total count;
    small outliers for the log-transformed number of expressed genes;
    large outliers for the proportion of read pairs assigned to mitochondrial genes;
    or large outliers for the proportion of read pairs assigned to spike-in transcripts.
These cells were removed from the data set prior to further analysis.

The cell cycle phase for each cell was identified using the `cyclone` classifier [@scialdone2015computational] implemented in the _scran_ package v1.2.0. 
This was performed with a set of human marker genes, identified by training the classifier on a pre-existing hESC data set [@leng2015oscope].

Finally, low-abundance genes with mean counts less than 1 were filtered out [@lun2016step].
This removes genes with low counts that contain little information for downstream methods such as normalization and HVG detection,
    thus reducing computational work and the severity of any multiple testing corrections.

## Normalization of cell-specific biases 

For the endogenous genes, cell-specific size factors were computed using the deconvolution method [@lun2016pooling] with pre-clustering.
For each gene, the count for each cell was divided by the appropriate size factor.
A pseudo-count of 1 was added, and the value was log~2~-transformed to obtain log-normalized expression values.
This was repeated using the spike-in transcripts, where the size factor for each cell was proportional to the sum of counts for all spike-in transcripts.
Note that different sets of size factors are necessary for different features, as spike-in transcripts are not subject to biases due to total RNA content [@lun2016step].

## Detecting correlated HVGs

To represent the technical variance, a mean-dependent trend was fitted to the variances of the spike-in transcripts using the `trendVar` function in the _scran_ package.
HVGs were identified as genes with variances that were significantly greater than the trend, after using the `decomposeVar` function to decompose the variance components.
Specifically, genes detected at a FDR of 5% and with biological components above 0.5 were considered to be HVGs [@lun2016step].
This was done while blocking on the batch in which each cell was sequenced, to ensure that large variances were not driven by uninteresting batch effects.

Correlations in the log-expression values between pairs of HVGs were also identified using the `correlatePairs` function in _scran_.
Significant correlations were identified at a FDR of 5%, and the genes in the significantly correlated pairs were used to define a set of correlated HVGs.
Again, blocking on the sequencing batch was performed during calculation of the correlations.

Identification of correlated HVGs was first performed using all (high-quality) cells in the data set.
The top 1000 HVGs with the largest biological components were identified, and the `removeBatchEffect` function from the _limma_ package v3.30.4 [@ritchie2015limma] was applied to eliminate the batch effect from their log-expression values.
A PCA plot was constructed using the corrected values to visualize the structure in the cell population.
Correlated HVGs were also identifed within the naive and primed conditions separately (see below).

## Testing for differential expression

Counts for the naive and primed cells within each batch were pooled to obtain four sets of pseudo-bulk counts [@lun2016overcoming].
Genes were tested for differential expression (DE) between naive and primed conditions, using the quasi-likelihood framework in the _edgeR_ package v3.16.3 [@chen2016reads]. 
The experimental design was parameterized using an additive design containing a condition term and the batch blocking factor.
DE genes were defined as those with significant differences between conditions at a FDR of 5%.

## Detecting the transition subpopulation 

Correlated HVGs were detected as previously described, but using only the cells in the naive condition.
The `removeBatchEffect` function was applied to remove the batch effect in the log-expression values of the correlated HVGs.
The corrected values were used for hierarchical clustering of the cells with the `hclust` function in _R_, using Ward linkage on the Euclidean distances.
Clusters of cells were identified using a tree cut, where the optimal number of clusters was determined by maximizing the average silhouette width.
The cluster of cells located between the bulk of cells from the naive and primed conditions in the PCA plot was denoted as the "transition" population.

The transition population was characterized by testing for differential expression relative to the other naive cells or to the primed cells.
This was performed by treating the log-expression values as microarray intensities and performing linear modelling with methods in _limma_.
(This corresponds to the "_limma_-trend" method described by @law2014voom, and is more stringent than `voom` or _edgeR_ for marker gene identification.)
For each contrast, several candidates were chosen from the top set of DE genes for further validation by staining and FACS. 

## Comparing transcriptional heterogeneity 

Correlated HVGs were detected as previously described, but using only the cells in the primed condition.
Clustering was not attempted as no clear separation between clusters was observed in the silhouette plots.
Instead, the variability of expression was compared between naive and primed cells.
We compared the total number of HVGs detected at a FDR of 5% in each condition;
    the sizes of the biological components for HVGs detected in both conditions;
    and the distribution of biological components for HVGs unique to each condition.
The sets of shared and primed-only HVGs were also tested for enrichment of GO or KEGG terms, using the `goana` and `kegga` functions, respectively, from _limma_.

## External temporal trajectories on naïve/primed map 

Naive marker genes were defined as those that were DE relative to primed cells (using the pseudo-bulk statistics, above) at a FDR of 5% and with a log~2~-fold change of 10;
    were present in at least 25% of naive cells; and were present in no more than 5% of primed cells.
Similarly, primed marker genes were defined as those that were DE relative to naive cells at a FDR of 5% and with a log~2~-fold change of -10;
    were present in at least 25% of primed cells; and were present in no more than 5% of naive cells.

A marker gene was considered to be expressed in a cell if its (normalized) count or CPM was greater than 10.
For each cell, we calculated the proportion of naive markers that were expressed.
This was repeated for the primed markers. 
Cells were mapped onto the "naive/primed axis" based on these proportions.
Large naive proportions and small primed proportions indicate that the cell is naive, and vice versa for primed cells.

Mapping onto the naive/primed axis was performed for cells collected from human pre-implantation embryos [@petropoulos2016single],
mouse embryos [@mohammed2016transcriptional], and cynomolgus monkey embryos [@nakamura2016developmental].
Mouse homologs for the marker genes were identified using the `getLDS` function from the _biomaRt_ package [@durinck2005biomart], using the homology relationships predicted by Ensembl.
Monkey homologs for marker genes were identified as those with the same gene symbol.
As a control, we also performed remapping using the naive, primed and transition cells in our own data set.

***

