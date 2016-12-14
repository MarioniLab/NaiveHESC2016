# Methods

## Alignment and quality control### Alignment and read countingThe raw FASTQ files were aligned to the human reference genome HG38 including ERCC sequences with the subread-align function in the subread package ([PMID: 23558742](https://www.ncbi.nlm.nih.gov/pubmed/23558742)).
We specified the alignment to paired end RNA-sequencing data setting a phred offset to +33 and selected only uniquely mapped reads. 
The aligned output BAM files were summarised with featureCounts in Rsubread ([PMID: 24227677](https://www.ncbi.nlm.nih.gov/pubmed/24227677)) using 4 threads. 
We merged resulting counts with SumTechReps in edgeR if they were technical replicates corresponding to the same cell.  
### QCQuality control and data analysis were conducted applying the single cell RNA sequencing analysis package scater following the developers suggested workflow ([doi: scater](http://dx.doi.org/10.1101/069633), [PMCID: PMC5112579](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5112579/)). 
First, we converted the count matrix to a SCESet object, calculated relevant QC metrics and defined spike-in transcripts. 
Then, we filtered for high quality cells in both batches separately by removing cells with a total count number less than three median-absolute deviations of the total count median. 
Cells were controlled for cell cycle phases with scran's cyclone function based on human cell cycle genes ([scran](https://bioconductor.org/packages/release/bioc/html/scran.html)). 
As cells were not predominantly found in one cell cycle phase, all high quality cells were kept. 
Next, we removed genes with a mean count less than 1 in order to only select genes of reasonable abundance.
Read counts were normalised after computing size factors for endogenous genes and spike factors for spike-ins with the normalise function in scater.## Data Analysis### Phenotype validation and characterisationWe identified highly variable genes (HVGs) using the biological component of logged expression variance for each gene with scran’s decomposeVar function ([scranworkflow](http://bioconductor.org/packages/devel/bioc/vignettes/scran/inst/doc/scran.html)). 
With a technical component based on the variance of the spike-ins and the respective mean expression of the genes, biological components were computed by subtracting the technical factor from the total variance of log~2~ expression values, indicating true biological variability.
HVGs were considered having a biological fold change greater than or equal to 2 with a maximum FDR of 5%. We ranked genes by decreasing biological component to sort for most distinct HVGs.The identified and sorted HVGs were used to compute correlations based on Spearman’s rho with the correlatePairs function in scran.
Against the null hypothesis that there were no gene correlations, we computed p values for each correlated pair and selected only the correlated genes with a FDR less than 0.05.This chosen set of correlated genes was then used to perform principal component analysis (PCA) in order to separate cells according to their biological variability (Figure 1B). We pooled naive and primed cells respectively, treating both subpopulations as pseudo bulk data sets, and observed differential expression (DE) between both phenotypes to identify marker genes for each condition ([PMID: 27008025](https://www.ncbi.nlm.nih.gov/pubmed/27008025)). 
Using the glmFit function in edgeR under the null hypothesis that there was no difference in gene expression between naïve and primed cells and setting a p value of 0.001, significantly up- or down-regulated genes were computed (Figure 1C).
### Naive subpopulation characterisation and cluster identificationSimilar to the identification of highly variable and significantly correlated genes for all cells, we determined HVGs and correlations for the naïve subpopulation with the only exception that we corrected for batch effects between the steps with the removeBatchEffect function in the Bioconductor package limma ([limma](https://bioconductor.org/packages/release/bioc/html/limma.html)).
Computing Euclidian distances of the log2-expression values between naive cells for the correlated HVGs allowed hierarchical clustering within the subpopulation. Here, we used the hclust function of the package stats and the agglomeration method ward.D2 ([stats citation -> standard R package](https://cran.r-project.org/doc/FAQ/R-FAQ.html#Citing-R)).
We used silhouette plots to determine the amount of clusters that resulted in most pronounced clusters.  

We characterised resulting clusters by analysing DE genes between each cluster and all remaining clusters. 
Clusters that showed only marginal differences in gene expression or whose distinction depended primarily on cell cycle effects were neglected. 
For clusters that were more prominent, DE was performed against the residual cells of the naïve subpopulation and against the primed cells. 
With the resultant DE genes, we assessed the grade of distinction of a particular cluster by visualising the expression values of cells within each subpopulation for those genes in a heat map (Figure 2A).
Additionally, obtained clusters were examined for their expression of specific marker genes, for instance previously discovered naïve markers (Figure 2B, [DOI](http://dx.doi.org/10.1016/j.stemcr.2016.02.005)) or eminently DE markers between naïve and primed subpopulation.
Finally, identified cells in the most prominent cluster were visualised along with the original subpopulations in a PCA plot that was coloured accordingly and functional enrichment analysis was performed to assess the GO functions of DE genes between the cluster and the residual naïve subpopulation.
### Primed subpopulation characterisation and variability analysisWe repeated all the steps that were conducted for the naïve subpopulation for the primed cells.  
Additionally, we compared the genetic heterogeneity of the naïve and the primed ESCs.
Therefore, we examined the overall number of HVGs as well as the variance of log expression for genes that were found highly variable in both phenotypes and uniquely for each phenotype. (Figure 3)
To exclude potential technical sources for gene expression variability, we counted the number of naïve and primed cells in each cell cycle and the total number of expressed genes in each population.
Besides, we compared the mean expression and the distribution of size factors of the HVGs between both conditions (Supp. Figure 3).
## External temporal trajectories on naïve/primed map 
We created a transcriptional map that places cells of interest relative to their similarity to the expression profiles of naïve and primed cells.
Here, we only focused on genes that were differentially expressed between the naïve and primed cells with a logFC greater than 10 (primed markers) and less than -10 (naïve markers).
Similarity of cells of interest to the naive and the primed subpopulation was indicated by count values greater than an arbitrary threshold of 10 for those marker genes (or orthologous genes for other species after [HomoloGene](https://www.ncbi.nlm.nih.gov/homologene?itool=toolbar) database). We compared the similarities to each subpopulation by calculating the proportions of the expressed markers.
After testing this model for our dataset (Supp. Figure 4), we downloaded raw count matrices of human ([PMID: 27062923](https://www.ncbi.nlm.nih.gov/pubmed/27062923)), mouse (Hisham, not yet published) and monkey ESCs ([PMID: 27556940 ](http://www.nature.com/nature/journal/v537/n7618/full/nature19096.html)) for various developmental stages and graphically contrasted naive against primed marker proportions (Figure 4).  
  

***  
***
