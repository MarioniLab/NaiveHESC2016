source("Figures/central.R")

# Figure 1B
sce <- readRDS("analysis/results-preprocess/sce_all.rds")
pcs <- readRDS(file.path("analysis/results-overall", "pcs.rds"))

pdf(file=file.path(figdir, "1b.pdf"), width=9, height=8, useDingbats=FALSE)
layout(cbind(1,2), width=c(5, 1))
par(mar=c(5.1, 4.2, 4.1, 1.1))
var1 <- round(attr(pcs, "percentVar")[1]*100)
var2 <- round(attr(pcs, "percentVar")[2]*100)
plot(pcs[,1], pcs[,2], col = ifelse(sce$phenotype=="naive", naive.col, primed.col),
     pch=16, xlab = paste0('Component 1: ', var1, '%'), ylab = paste0('Component 2: ', var2, '%') ,
     cex.lab = 1.5 , cex.axis=1.2)
par(mar=c(5.1, 0.1, 4.1, 0.1))
plot.new()
legend("topleft", legend=c("Naive", "Primed"), col=c(naive.col, primed.col), pch=16, cex=1.5, bty='n')
dev.off()

# Figure 1C
de.res <- read.table(file.path("analysis/results-overall", "de.tsv"), header=TRUE)
is.up <- de.res$logFC > 0 & de.res$FDR <= 0.05
is.down <- de.res$logFC < 0 & de.res$FDR <= 0.05
x <- de.res$logCPM
y <- de.res$logFC
main.col <- "grey50"
naive.genes <- c("KLF4", "KLF17", "DPPA3", "DNMT3L", "DPPA5", "GATA6", 'TBX3', "IL6ST", "KLF5")
primed.genes <- c("DUSP6", "THY1", "CD24", "ZIC2", "SFRP2")
not.de <- c("OTX2", "T", 'SOX2', "POU5F1", "NANOG")

pdf(file=file.path(figdir, "1c.pdf"), useDingbats=FALSE)
par(mar = c(5.1, 5.1, 4.1, 2.1), las=1)
plot(x, y, cex=0.3, ylim=c(-22.5,22.5), pch = 16, ylab=expression(Log[2]~"fold change"), 
     xlab=expression("Average"~Log[2]~"CPM"), cex.lab=1.5, col=main.col, cex.axis=1.2)
points(x[is.up], y[is.up], cex=0.3, col=naive.col, pch = 16)
points(x[is.down], y[is.down], cex=0.3, col=primed.col, pch = 16)
coords <- par()$usr
legend(coords[2]-0.5, coords[3]+1, legend=c(sprintf("Upregulated in naive (%i)", sum(is.up)), 
                                            sprintf("Upregulated in primed (%i)", sum(is.down)),
                                            sprintf("Not DE (%i)", sum(!is.up & !is.down))), 
       col=c(naive.col, primed.col, main.col), pch=16, yjust=0, xjust=1)

m <- match(naive.genes, rownames(de.res))
points(x[m], y[m], col='#D55E00', pch=16, cex=1.2)
pos <- c(1,3,4,4,4,2,3,3,4)
text(x[m], y[m], labels = naive.genes, col='#D55E00', pos=pos, offset=0.5)

m <- match(primed.genes, rownames(de.res))
points(x[m], y[m], col='royalblue4', pch=16, cex=1.2)
pos <- c(4,2,4,4,4)
text(x[m], y[m], labels = primed.genes, col='royalblue4', pos=pos, offset=0.5)

m <- match(not.de, rownames(de.res))
points(x[m], y[m], col='grey20', pch=16, cex=1.2)
pos <- c(4,4,4,2,4)
text(x[m], y[m], labels = not.de, col='black', pos=pos, offset=0.5)
dev.off()

# Figure S1a
library(scater)
ugly1 <- plotColData(sce, x="sample", y="log10_total_counts")
ugly2 <- plotColData(sce, x="sample", y="log10_total_features_by_counts")
ugly3 <- plotColData(sce, x="sample", y="pct_counts_ERCC")
ugly4 <- plotColData(sce, x="sample", y="pct_counts_Mt")

pdf(file=file.path(figdir, "s1a.pdf"),width = 10, height = 10, useDingbats=FALSE)
multiplot(ugly1 + geom_boxplot(width=0.1) + theme_classic() + ylab('Total Counts [log10]') + xlab(''),
ugly2 + geom_boxplot(width=0.1) + theme_classic() + ylab('Total Features [log10]') + xlab('Batch'),
ugly3 + geom_boxplot(width=0.1) + theme_classic() + ylab('PCT Counts ERCC') + xlab(''),
ugly4 + geom_boxplot(width=0.1) + theme_classic() + ylab('PCT Counts Mt') + xlab('Batch'),
cols=2)
dev.off()

# Figure S1b
pdf(file=file.path(figdir, "s1b.pdf"),width = 20, height = 12, useDingbats=FALSE)
par(las=1, mar = c(2.1, 4.5, 0.1, 0.5), mfrow=c(2,1))
markers <- c("KLF4", "HORMAD1", "KHDC3L", "ALPP", "ALPPL2", "ZNF729", "TRIM60", "SOX11", "CYTL1", "HMX2", "THY1", "DUSP6", "PTPRZ1")

for (ptype in c("naive", "primed")) {
  object <- sce[,sce$phenotype==ptype]
  ugly.plot <- plotExpression(object, features = markers, col = "phenotype")
  plot_data <- ggplot_build(ugly.plot)
  
  pcol <- switch(ptype, naive=naive.col, primed=primed.col, transition=trans.col)      
  plot(plot_data$data[[2]]$x, plot_data$data[[2]]$y, col = pcol, 
       ylab=bquote(.(stuff)~log[2]*"-expression", list(stuff=paste0(toupper(substring(ptype, 1, 1)), substring(ptype, 2)))),
       pch=16, xlab = "" , cex.lab = 1.5, xaxt = "n", cex=2, ylim=c(0,13.3), cex.axis=1.2)
  if (ptype == "primed"){
    axis(1, at=seq_along(markers), labels = levels(plot_data$plot$data$Feature), cex.axis = 1.5)
  } else {
    axis(1, at=seq_along(markers), labels=character(length(markers)))
  }
}
dev.off()

#Figure S1C: Proteomics

# Figure S1d
DE_pastor <- readRDS(file = file.path('analysis/results-bulk', 'de_pastor.rds'))
DE_theunissen <- readRDS(file = file.path('analysis/results-bulk', 'de_theunissen.rds'))

naive.markers <- rownames(de.res[de.res$logFC > 10 & de.res$FDR <= 0.05,])
primed.markers <- rownames(de.res[de.res$logFC < -10 & de.res$FDR <= 0.05,])

naive.markers <- naive.markers[naive.markers %in% rownames(DE_pastor) & naive.markers %in% rownames(DE_theunissen)]
primed.markers <- primed.markers[primed.markers %in% rownames(DE_pastor) & primed.markers %in% rownames(DE_theunissen)]

pdf(file=file.path(figdir, "s1d.pdf"),width = 12, height = 6, useDingbats=FALSE)

layout(cbind(1,2,3), width=c(5,5,1))
par(mar=c(5.1, 4.2, 4.1, 1.1))
plot(DE_pastor$logFC, -log10(DE_pastor$PValue), pch = 16, xlab = 'Log FC', ylab = '- Log10(Pval)', bty = 'n', main = 'Pastor, 2016', xlim = c(-11, 11))
points(DE_pastor[naive.markers,]$logFC, -log10(DE_pastor[naive.markers,]$PValue), pch = 16, cex = 1.5, col = naive.col)
points(DE_pastor[primed.markers,]$logFC, -log10(DE_pastor[primed.markers,]$PValue), pch = 16, cex = 1.5, col = primed.col)
max.pval <- min(-log10(DE_pastor[DE_pastor$FDR < 0.05,]$PValue))
lines(x = c(-11, 11), y = c(max.pval, max.pval), col = 'red', lty = 2, lwd = 3)
text(x = -9, y = max.pval, labels = 'FDR < 0.05', pos = 1, col = 'red')

plot(DE_theunissen$logFC, -log10(DE_theunissen$PValue), pch = 16, xlab = 'Log FC', ylab = '', bty = 'n', main = 'Theunissen, 2016', xlim = c(-15, 15))
points(DE_theunissen[naive.markers,]$logFC, -log10(DE_theunissen[naive.markers,]$PValue), pch = 16, cex = 1.5, col = naive.col)
points(DE_theunissen[primed.markers,]$logFC, -log10(DE_theunissen[primed.markers,]$PValue), pch = 16, cex = 1.5, col = primed.col)
max.pval <- min(-log10(DE_theunissen[DE_theunissen$FDR < 0.05,]$PValue))
lines(x = c(-15, 15), y = c(max.pval, max.pval), col = 'red', lty = 2, lwd = 3)
text(x = -12, y = max.pval, labels = 'FDR < 0.05', pos = 1, col = 'red')

par(mar=c(5.1, 0.1, 4.1, 0.1))
plot.new()
legend("topleft", legend=c("Naive", "Primed"), col=c(naive.col, primed.col), pch=16, cex=1.5, bty='n')

dev.off()

