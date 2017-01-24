# Figure captions

**Figure 1. Naive and primed human ESCs exhibit strong differences in gene expression.**
(A) _Experimental setup (Ferdinand)_
(B) PCA plot of hESC expression profiles, constructed from batch-corrected and normalized log-expression values of highly variable genes detected across the entire data set.
Cells are coloured by their condition and the percentage of variance explained by the first two principal components is shown.
(C) Smear plot of log~2~-fold changes in expression between the naive and primed conditions, where DE genes were detected using _edgeR_ at a FDR of 5%.

**Figure 2. The transition population is transcriptionally distinct from the other naive and primed cells.**
(A) Heat map of the top 50 genes with the strongest differential expression between the naive and transition cells (top), or between the transition and primed cells (bottom).
The box for each cell (column) and gene (row) is coloured according to the log~2~-fold change from the average expression for each gene.
(B) Log~2~-expression profiles of selected marker genes across cells in the naive, transition and primed populations.
Each point represents a cell in the corresponding population.
(C) _FaCS (Ferdinand)_

**Figure 3. The primed population has greater transcriptional heterogeneity than the naive population.**
(A) Estimates of the biological component of the variance for highly variable genes detected in each condition.
The number below each boxplot represents the total number of detected HVGs.
Points represent HVGs with biological components that are more than 1.5 interquartile ranges from the third quartile.
(B) Biological components for cell cycle HVGs in each condition, ranked by the magnitude of the component.
Association with the cell cycle is determined by GO annotation (GO:0007049), and some key genes are highlighted.
Biological components of ERCC spike-in transcripts are shown as a control.
(C) Distribution of biological components for splicing HVGs detected in either condition.
The number of genes is shown in brackets, and some key genes are highlighted.

**Figure 4. Cells shift from a naive-like to a primed-like expression pattern during early embryonic development.**
Naive and primed markers were identified from the DE analysis of the hESC data.
Upon profiling its transcriptome, a cell can be mapped onto the naive/primed axis based on the proportions of naive and primed markers that it expresses.
This was performed for cells derived from human pre-implantation embryos (A), mouse embryos (B) and cynomolgus monkey embryos (C).
For each plot, the density of cells is represented by the colour of the pixels.
Cells on the red line have equal proportions of expressed primed and naive markers.

# Supplements

**Figure 1.** No Figure 1 yet.

**Figure 2. Cells in the transition state lie between the naive and primed populations.**
Cells in the transition populations were identified after clustering on HVGs identified in the naive condition.
These cells are highlighted in green, using the same PCA plot in Figure 1B.

**Figure 3. Increased variability in primed cells is not driven by technical effects.**
(A) The number of expressed genes, number of cells and size factors of cells in the naive and primed conditions.
(B) Left: the distribution of biological components for all HVGs in each condition, separated into quartiles for greater resolution.
Right: the distribution of the ratio of the total variance to the technical component for all HVGs in each condition, separated into quartiles.
(C) Distribution of cells across cell cycle phases in each condition, as determined using the `cyclone` method. 

**Figure 4. Proof of concept of naive/primed mapping with naive, transition and primed hESCs.** 
For each plot, the density of cells is represented by the colour of the pixels.
Cells on the red line have equal proportions of expressed primed and naive markers.

