source("central.R")

# Figure 1B
sce <- readRDS("sce_all.rds")
chosen <- readRDS(file.path("results-overall", "cor_hvg.rds"))
library(scater)
ugly.plot <- plotPCA(sce, exprs_values="norm_exprs", feature_set = chosen)
plot_data <- ggplot_build(ugly.plot)

pdf(file=file.path(figdir, "1b.pdf"), width=9, height=8)
layout(cbind(1,2), width=c(5, 1))
par(mar=c(5.1, 4.2, 4.1, 1.1))
plot(plot_data$data[[2]]$x, plot_data$data[[2]]$y, col = ifelse(sce$phenotype=="naive", naive.col, primed.col),
     pch=16, xlab = plot_data$plot$labels$x , ylab = plot_data$plot$labels$y, cex.lab = 1.5 , cex.axis=1.2)
par(mar=c(5.1, 0.1, 4.1, 0.1))
plot.new()
legend("topleft", legend=c("Naive", "Primed"), col=c(naive.col, primed.col), pch=16, cex=1.5)
dev.off()

# Figure 1C
de.res <- read.table(file.path("results-overall", "de.tsv"), header=TRUE)
is.up <- de.res$logFC > 0 & de.res$FDR <= 0.05
is.down <- de.res$logFC < 0 & de.res$FDR <= 0.05
x <- de.res$logCPM
y <- de.res$logFC
main.col <- "grey50"
naive.genes <- c("KLF4", "KLF17", "DPPA3", "DNMT3L", "DPPA5")
primed.genes <- c("DUSP6", "THY1")

pdf(file=file.path(figdir, "1c.pdf"))
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
points(x[m], y[m], col=naive.col, pch=16, cex=1.2)
pos <- c(4,3,4,4,4)
text(x[m], y[m], labels = naive.genes, col=naive.col, pos=pos, offset=0.5)

m <- match(primed.genes, rownames(de.res))
points(x[m], y[m], col=primed.col, pch=16, cex=1.2)
pos <- c(4,3)
text(x[m], y[m], labels = primed.genes, col=primed.col, pos=pos, offset=0.5)
dev.off()


