# Methods

## Alignment and quality control
We specified the alignment to paired end RNA-sequencing data setting a phred offset to +33 and selected only uniquely mapped reads. 
The aligned output BAM files were summarised with featureCounts in Rsubread ([PMID: 24227677](https://www.ncbi.nlm.nih.gov/pubmed/24227677)) using 4 threads. 
We merged resulting counts with SumTechReps in edgeR if they were technical replicates corresponding to the same cell.  

First, we converted the count matrix to a SCESet object, calculated relevant QC metrics and defined spike-in transcripts. 
Then, we filtered for high quality cells in both batches separately by removing cells with a total count number less than three median-absolute deviations of the total count median. 
Cells were controlled for cell cycle phases with scran's cyclone function based on human cell cycle genes ([scran](https://bioconductor.org/packages/release/bioc/html/scran.html)). 
As cells were not predominantly found in one cell cycle phase, all high quality cells were kept. 
Next, we removed genes with a mean count less than 1 in order to only select genes of reasonable abundance.
Read counts were normalised after computing size factors for endogenous genes and spike factors for spike-ins with the normalise function in scater.
With a technical component based on the variance of the spike-ins and the respective mean expression of the genes, biological components were computed by subtracting the technical factor from the total variance of log~2~ expression values, indicating true biological variability.
HVGs were considered having a biological fold change greater than or equal to 2 with a maximum FDR of 5%. We ranked genes by decreasing biological component to sort for most distinct HVGs.
Against the null hypothesis that there were no gene correlations, we computed p values for each correlated pair and selected only the correlated genes with a FDR less than 0.05.
Using the glmFit function in edgeR under the null hypothesis that there was no difference in gene expression between naïve and primed cells and setting a p value of 0.001, significantly up- or down-regulated genes were computed (Figure 1C).

Computing Euclidian distances of the log2-expression values between naive cells for the correlated HVGs allowed hierarchical clustering within the subpopulation. Here, we used the hclust function of the package stats and the agglomeration method ward.D2 ([stats citation -> standard R package](https://cran.r-project.org/doc/FAQ/R-FAQ.html#Citing-R)).
We used silhouette plots to determine the amount of clusters that resulted in most pronounced clusters.  

We characterised resulting clusters by analysing DE genes between each cluster and all remaining clusters. 
Clusters that showed only marginal differences in gene expression or whose distinction depended primarily on cell cycle effects were neglected. 
For clusters that were more prominent, DE was performed against the residual cells of the naïve subpopulation and against the primed cells. 
With the resultant DE genes, we assessed the grade of distinction of a particular cluster by visualising the expression values of cells within each subpopulation for those genes in a heat map (Figure 2A).
Additionally, obtained clusters were examined for their expression of specific marker genes, for instance previously discovered naïve markers (Figure 2B, [DOI](http://dx.doi.org/10.1016/j.stemcr.2016.02.005)) or eminently DE markers between naïve and primed subpopulation.


Additionally, we compared the genetic heterogeneity of the naïve and the primed ESCs.
Therefore, we examined the sheer number of HVGs as well as the variance of log expression for genes that were found highly variable in both phenotypes and uniquely for each phenotype. (Figure 3)
To exclude potential technical sources for gene expression variability, we counted the number of naïve and primed cells in each cell cycle and the total number of expressed genes in each population.
Besides, we compared the mean expression and the distribution of size factors of the HVGs between both conditions (Supp. Figure 3).


Here, we only focused on genes that were differentially expressed between the naïve and primed cells with a logFC greater than 10 (primed markers) and less than -10 (naïve markers).
Similarity of cells of interest to the naive and the primed subpopulation was indicated by count values greater than an arbitrary threshold of 10 for those marker genes (or orthologous genes for other species after [HomoloGene](https://www.ncbi.nlm.nih.gov/homologene?itool=toolbar) database). We compared the similarities to each subpopulation by calculating the proportions of the expressed markers.

  

***  
***
