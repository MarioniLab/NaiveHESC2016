# Figure Captures

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

**Figure 2. PCA plot of log-expression values of human ESCs coloured according to their respective subpopulation indicates a transition state from naive to primed cells.**
The putative transition subpopulation was identified by clustering naive cells into three subsets based on euclidian distances and selecting the most protruded one.
Initial naive cells in that transition cluster were coloured in green, suggesting the additional subpopulation between naive and primed.  

**Figure 3. Augmented variability of gene expression in primed than in naive subpopulation is not driven by cell cycle or technical effects.**  
(A) Overall numbers of primed and naive cells in cell cycle phases G1, G2M and S. 
Cells were categorised by means of human cell cycle genes and the cyclone function in scran.  
(B) Total number of expressed genes in primed and naive ESCs and total number of expressed genes present in both subsets.
Low-abundance genes with mean counts less than 1 were neglected.  
(C) Histogram of logged means of all genes with a mean count equal to or greater than 1.   
(D) Distribution of size factors of endogenous genes with a mean count equal to or greater than 1.   

**Figure 4. Proportions of significant expression of naive and primed markers for naive, transition and primed cells demonstrate proof of concept.** 
Expression was considered significant if counts for respective marker genes were greater than 10. 
The red line indicates an equilibrium of significantly expressed primed and naive markers.
Colouring shows the number of cells with the exact same proportions.

