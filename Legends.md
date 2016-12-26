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
(A) Number of highly variable genes detected in each condition at a FDR of 5% and with a biological component of the variance greater than or equal to 0.5. 
(B) Estimates of the biological component of the variance in the primed condition, plotted against the corresponding estimates in the naive condition.
Each point is a gene, coloured based on whether it is detected as a HVG in both conditions (black) and if it is involved in the cell cycle (GO:0000278, pink).
The red line represents equality between biological components.
(C) Distribution of biological components for genes detected as HVGs only in the primed or naive conditions.
Inset plots show the distribution of biological components for HVGs in each condition that are involved in cilia assembly (GO:0042384) or brain development (GO:0007420).
The number of condition-specific HVGs in each plot is also shown.

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
(A) Diagnostic plots showing the number of expressed genes, number of cells and size factors of cells in the naive and primed conditions.
(B) Changes in the biological components between naive and primed conditions, plotted against the log~2~-fold change in expression.
Density plots are shown for all genes, genes detected as HVGs in both conditions, and condition-specific HVGs.
(C) Mean-variance relationships in the naive and primed conditions.
The technical variance is represented by the red line fitted to spike-in variances.
HVGs are detected in each condition as those genes with variances significantly greater than the trend, and are coloured.
(D) Distribution of cells across cell cycle phases in each condition, as determined using the `cyclone` method. 

**Figure 4. Proof of concept of naive/primed mapping with naive, transition and primed hESCs.** 
For each plot, the density of cells is represented by the colour of the pixels.
Cells on the red line have equal proportions of expressed primed and naive markers.

