# Create figure directory and specifiy path of necessary R-objects
dir.create(file.path("./figures"), showWarnings = F)
objectpath <- file.path("./objects/")

# Load Colors 
primed.col <- readRDS(file = paste0(objectpath, "primed.col"))
naive.col <- readRDS(file = paste0(objectpath, "naive.col"))

# Figure 3

# A
hvg.matrix <- readRDS(file = paste0(objectpath, "hvg.matrix"))
colnames(hvg.matrix) <- c("Primed", "Naive")

pdf(file="./figures/figure3a.pdf", width = 4)
par(mar = c(5.1, 5.5, 4.1, 2.1), las = 1)
plot <- barplot(hvg.matrix, ylab="Highly Variable Genes",
              ylim = c(0,7500), col = c(primed.col, naive.col), beside = TRUE,
              xlim = c(0,2), width = 0.7, space = c(0,0.5), cex.axis = 0.75, cex.lab = 1.3, cex.names = 1.3)
text(x = plot, y = hvg.matrix, label = hvg.matrix, cex = 0.75, pos = 3)

dev.off()


# B
library("org.Hs.eg.db")
library("limma")
hvg.out_naive <- readRDS(file = paste0(objectpath, "hvg.naive"))
hvg.out_primed <- readRDS(file = paste0(objectpath, "hvg.primed"))
shared_genes <- intersect(rownames(hvg.out_naive), rownames(hvg.out_primed))

ensembl <- readRDS(file = paste0(objectpath, "ensemblGenes"))

entrez <- ensemblGenes[match(shared_genes, ensemblGenes$external_gene_name),]
shared_genes <- shared_genes[-which(is.na(entrez$entrezgene))]
entrez <- entrez[-which(is.na(entrez$entrezgene)),]
entrez_genes <- entrez$entrezgene

GoTerms <- goana(entrez_genes, species = "Hs")

# [567] "regulation of cell proliferation"   GO:0042127         
# [629] "regulation of DNA methylation"  GO:0044030
# [688] "embryo implantation"          GO:0007566                                                                              
    
cell.prol <- AnnotationDbi::select(org.Hs.eg.db, keys="GO:0042127", keytype="GOALL", column="ENSEMBL")
DNA.meth <- AnnotationDbi::select(org.Hs.eg.db, keys="GO:0044030", keytype="GOALL", column="ENSEMBL")
embryo.implant <- AnnotationDbi::select(org.Hs.eg.db, keys="GO:0007566", keytype="GOALL", column="ENSEMBL")

ensembl_ids <- ensemblGenes[match(shared_genes, ensemblGenes$external_gene_name),]
ensembl_ids <- ensembl_ids$ensembl_gene_id
cell.prol <- intersect(ensembl_ids, cell.prol$ENSEMBL)
DNA.meth <- intersect(ensembl_ids, DNA.meth$ENSEMBL)
embryo.implant <- intersect(ensembl_ids, embryo.implant$ENSEMBL)

gene.function <- data.frame("Color"=rep("grey20", length(ensembl_ids)), row.names = ensembl_ids)
gene.function$Color <- as.character(gene.function$Color)
gene.function[cell.prol,1] <- "coral"
gene.function[DNA.meth,1] <- "darkgoldenrod1"
gene.function[embryo.implant,1] <- "deepskyblue"


pdf(file="./figures/figure3b.pdf")
par(bty='l', las = 1, mar = c(5.1, 5.5, 4.1, 2.1))
plot(hvg.out_naive[shared_genes,]$bio, hvg.out_primed[shared_genes,]$bio, xlim=c(0.5,9), ylim=c(0.5,9), 
     pch = 16, cex = 1,  ylab = "Biological Variance [Primed]", col=gene.function[,1], 
     xlab = "Biological Variance [Naive]", cex.lab=1.3)
abline(0,1, col = "indianred2", lty = 1)
legend("topright", legend = c("Cell proliferation", "DNA methylation", "Embryo implantation"), 
       col = c("coral", "darkgoldenrod1", "deepskyblue"), bty="n", pch=16, cex)
#text(hvg.out_naive[top_genes,]$bio+0.5, hvg.out_primed[top_genes,]$bio, labels = top_genes)
dev.off()


# C
primed_genes <- unique(rownames(hvg.out_naive), rownames(hvg.out_primed))
naive_genes <- unique(rownames(hvg.out_primed), rownames(hvg.out_naive))

data <- list(primed = data.frame(hvg.out_primed[primed_genes,]$bio), naive = data.frame(hvg.out_naive[naive_genes,]$bio))
data <- lapply(data, function(x) { x[1:nrow(hvg.out_primed),] })
data <- cbind(data$primed, data$naive)
data <- as.data.frame(data)
colnames(data) <- c("Primed", "Naive")

pdf(file="./figures/figure3c.pdf", width = 4)
par(bty='n', las = 1, mar = c(3.1, 5.5, 4.1, 2.1))
boxplot(data, col = c(primed.col, naive.col), boxwex = 0.5, ylim = c(-1,10), cex.lab = 1.3,
        outline=FALSE, ylab = "Biological Variance", xaxt = "n")
text(x = c(1,2), y=c(-1,-1), labels = c("Primed", "Naive"), cex=1.3)
dev.off()


# Supplements

sce_primed <- readRDS(file = paste0(objectpath, "new_sce_primed_object"))
sce_naive <- readRDS(file = paste0(objectpath, "new_sce_naive_object"))

# A

par(mfrow=c(1,1))
primed.phases <- table(sce_primed$phase)
naive.phases <- table(sce_naive$phase)

pdf(file="./figures/sup.figure3a.pdf")
barplot(rbind(primed.phases, naive.phases), beside = TRUE, col = c(primed.col, naive.col), ylab="Number of cells")
legend("top", legend = c("Primed", "Naive"), fill = c(primed.col, naive.col), bty="n", cex=1.5)
dev.off()

# B
naive.means <- log(rowMeans(counts(sce_naive)))
primed.means <- log(rowMeans(counts(sce_primed)))

naive.genenum <- table(naive.means > 1)["TRUE"]
primed.genenum <- table(primed.means > 1)["TRUE"]
shared <- length(intersect(rownames(sce_naive[naive.means >1,]), rownames(sce_primed[primed.means >1,])))

num.matrix <- cbind(primed.genenum, naive.genenum, shared)
colnames(num.matrix) <- c("Primed", "Naive", "Shared")

pdf(file="./figures/sup.figure3b.pdf")
par(mfrow=c(1,1), mar = c(5.1, 5.5, 4.1, 2.1), las = 1)
plot <- barplot(num.matrix, col = c(primed.col, naive.col, "#CC79A7"), beside = TRUE, names.arg = c("Primed", "Naive", "Shared"), 
                ylab="Number of expressed genes", xlim = c(0,3), width = 0.3, space = c(0,.5),
                cex.axis = .75, cex.lab = 1.3, cex.names = 1., ylim=c(0,14700))
text(x = plot, y = num.matrix, label = num.matrix, cex = 0.75, pos = 3)
dev.off()

# C

pdf(file="./figures/sup.figure3c.pdf")
par(mfrow=c(1,2), mar=c(5.1, 4.1, 4.1, 0.1))
hist(primed.means, xlim = c(0,11), breaks = 100, col=primed.col, xlab="log Mean", main="Primed")
hist(naive.means, xlim = c(0,11), breaks = 100, col=naive.col, xlab="log Mean", main="Naive")
dev.off()

# D

naive.size <- sce_naive$size_factor
primed.size <- sce_primed$size_factor

pdf(file="./figures/sup.figure3d.pdf")
par(mfrow=c(1,2), mar=c(5.1, 4.1, 4.1, 1), las=1)
hist(primed.size, xlim = c(0,8), breaks = 50, col=primed.col, xlab="Size Factor", main="Primed")
hist(naive.size, xlim = c(0,8), breaks = 50, col=naive.col, xlab="Size Factor", main="Naive", ylab = "")
dev.off()
