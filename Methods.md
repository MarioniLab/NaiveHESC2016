# Methods
## Alignment and Quality Control### Alignment and CountingThe raw FASTQ files were aligned to the human reference genome HG38 and the ERCC sequences ??? (WHICH ONES?) with the subread-align function of subread package (The Subread aligner: fast, accurate and scalable read mapping by seed-and-vote.).Liao Y, Smyth GK and Shi W (2013).e.
We specified the alignment to paired end RNA-sequencing data setting a phred offset to +33 and selected only uniquely mapped reads. 
The aligned output BAM files were summarized with featureCounts of Rsubread (Liao et al., 2013) using 4 threads. 
Finally, we merged sequences in R if they were technical replicates corresponding to the same cell.  
### QC
Quality control and data analysis were conducted applying the single cell RNA sequencing analysis Bioconductor package Scater following the developers suggested workflow (PAPER SCATER , Aaron WORKFLOW  → Paying tribute). 
First, we converted the count matrix to a SCESet object, calculated QC metrics and defined spike-in transcripts. 
Then, we filtered for high quality cells in both batches separately by removing cells with a total count value greater than three median-absolute deviations less than the total count median. 
Cells were controlled for cell cycle phases with the cyclone function based on cell cycle genes by Bioconductor package Scran. 
As cells were not predominantly found in one cell cycle phase, all high quality cells were kept. 
Next, we removed genes with a mean count less than 1 in order to select only reasonably expressed genes. 
Read counts were normalised after computing size factors for endogenous genes and spike factors for spike-ins with the normalise function (default settings).## Data Analysis### Whole dataset (OTHER TITLE: Phenotype validation?)We identified highly variable genes (HVGs) using the biological component of logged expression variance for each gene with Scran’s deceomposeVar function. 
With the technical component based on the variance of the spike-ins and the respective mean expression of a gene, biological components were computed by subtracting the technical factor from the total variance of log2 expression values, indicating true biological variability.(NECESSARY?)
HVGs were considered having a biological fold change greater than or equal to 2 with a maximum FDR of 5%. We ranked genes by decreasing biological component to sort for the most distinct HVGs.The identified HVGs were used to compute correlations based on Spearman’s rho with the correlatePairs function of Scran.
Against the null hypothesis that there were no gene correlations, we computed p-values for each correlated pair and selected only the correlated genes with a FDR less than 0.05.This chosen set of correlated genes was then used to perform principal component analysis (PCA) in order to separate cells according to their biological variability (Figure 1b). Considering/treating the naïve and primed subpopulation as bulkdata(SYNONYM), we observed differential expression (DE) between both phenotypes to identify marker genes for each condition. 
Using the glmFit function of edgeR package under the null hypothesis that there is no difference in gene expression between naïve and primed cells and with a P-value of 0.001, significantly up- or down-regulated genes were computed (Figure 1C).
### Naive Subpopulation
Similar to the identification of HVGs and significantly correlated genes for all cells, we determined HVGs and correlations for the naïve subpopulation with the only exception that we corrected for batch effects between the steps with the removeBatchEffect function of the Bioconductor package limma (SOURCE PAPER MISSING).
Computing Euclidian distances of the log2-expression values between naive cells for the correlated HVGs allowed hierarchical clustering within the subpopulation with the hclust function of the package stats and the agglomeration method ward.D2.
We used silhouette plots to determine the amount of clusters that resulted in most pronounced clusters. 
Then, we evaluated resulting clusters by analysing differentially expressed genes between each cluster and all respectively remaining clusters. 
Clusters that showed only marginal differences in gene expression or whose distinction depended primarily on cell cycle effects were neglected. 
For clusters that were more prominent, differential expression was performed against the residual cells of the naïve subpopulation and against the primed cells while each subpopulation was treated as BULK DATA. 
With the resultant differentially expressed genes, we could assess the grade of distinction of a particular cluster by visualising the expression values of cells within each subpopulation for those genes in a heat map (Figure 2a).
Additionally, obtained clusters were examined for their expression of specific marker genes, for instance naïve markers (Figure 2b, GUO PAPER) or the top differentially expressed markers between naïve and primed subpopulation (Suppl Figure 2b).  Eventually, identified cells in the most prominent cluster were visualized along with the original subpopulations in a PCA plot that was coloured accordingly (FIGURE SUPP PCA) and Gene Set enrichment analysis was performed to assess the functions of differentially expressed genes between the cluster and the residual naïve subpopulation.
### Primed Subpopulation
We repeated all the steps that were conducted for the naïve subpopulation for the primed cells.  
Additionally, we compared the genetic heterogeneity of the naïve and the primed ESCs.
Therefore, we examined the sheer number of HVGs (Figure 3a) and the variance of log expression for both genes that were found highly variable for both phenotypes (Figure 3b) and highly variably uniquely for each phenotype (Figure 3c).
NECESSARY? :To exclude potential technical sources for gene expression variability, we counted the number of naïve and primed cells in each cell cycle (Supp. Figure 3a) as well as the total number of expressed genes in each population (Supp. Figure 3b).
Besides, we compared the mean gene expression (Supp. Figure 3c) and the distribution of size factors (Supp. Figure 3d) between both conditions.))
## External temporal trajectories on naïve/primed map 
To conjecture the possible temporal position/location/placement/etc. of both phenotypes, we created a transcriptional map that places cells of interest (COI) relative to their similarity to the expression profiles of naïve and primed cells.
Here, we only focused on genes that were differentially expressed between naïve and primed with a logFC greater than 10 (primed markers) and less than -10 (naïve markers).
Subsequently, count values of COI equal to or greater than an arbitrary threshold of 10 for those marker genes (or orthologous genes for other species) indicated their similarity to naïve and primed cells, respectively.  
After testing the model for our?? cells (Supp. Figure 4), we downloaded raw count matrices of human (Petropoulos ), mouse (Hisham ) and monkey ESCs (Monkeys ) for various developmental stages and computed proportions of marker expression (Figure 4).  
  
  

***  
***

## References

* The Subread aligner: fast, accurate and scalable read mapping by seed-and-vote, Liao Y. et al.* Liao Y, Smyth GK and Shi W. featureCounts: an efficient general-purpose program for assigning sequence reads to genomic features. Bioinformatics, 2013 Nov 30.)
* McCarthy D, Wills Q and Campbell K (2016). scater: Single-cell analysis toolkit for gene expression data in R. R package version 1.2.0, https://github.com/davismcc/scater. 
* A step-by-step workflow for low-level analysis of single-cell RNA-seq data with Bioconductor - Aaron T. L. Lun1, Davis J. McCarthy2 and John C. Marioni3
* Using scran to perform basic analyses of single-cell RNA-seq data – Aaron Lun
* Using scran to perform basic analyses of single-cell RNA-seq data – Aaron Lun
* Single-Cell RNA-Seq Reveals Lineage and X Chromosome Dynamics in Human Preimplantation Embryos. – Petropoulos S
* Hisham: Not published yet
* A developmental coordinate of pluripotency among mice, monkeys and humans