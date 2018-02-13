source("central.R")

# Figure 1B
sce <- readRDS("sce_all.rds")
pcs <- readRDS(file.path("results-overall", "pcs.rds"))

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
points(x[m], y[m], col=naive.col, pch=16, cex=1.2)
pos <- c(4,3,4,4,4)
text(x[m], y[m], labels = naive.genes, col=naive.col, pos=pos, offset=0.5)

m <- match(primed.genes, rownames(de.res))
points(x[m], y[m], col=primed.col, pch=16, cex=1.2)
pos <- c(4,3)
text(x[m], y[m], labels = primed.genes, col=primed.col, pos=pos, offset=0.5)
dev.off()

# Figure S1b
pdf(file=file.path(figdir, "s1b.pdf"),width = 10, height = 10, useDingbats=FALSE)
markers <- c("KLF17", "DNMT3L", "DPPA5", "DPPA3", "KLF4", "THY1", "DUSP6")
ugly.plot <- plotExpression(sce, features = markers, col = "phenotype")
plot_data <- ggplot_build(ugly.plot)
plot_data$data[[1]]$colour[plot_data$data[[1]]$colour == '#729ECE'] <- naive.col
plot_data$data[[1]]$colour[plot_data$data[[1]]$colour == '#FF9E4A'] <- primed.col

plot(plot_data$data[[1]]$x, plot_data$data[[1]]$y, col = plot_data$data[[1]]$colour, 
       ylab=bquote(~log[2]*"-expression"), pch=16, xlab = "" , cex.lab = 1.5, xaxt = "n", cex=1.2, ylim=c(0,13), cex.axis=1.2)
axis(1, at=seq_along(markers), labels = levels(plot_data$plot$data$Feature), cex.axis = 1.5)
legend("topright", legend=c("naive", "primed"), col=c(naive.col, primed.col, main.col), pch=16, cex=1.5, yjust=0, xjust=1)
dev.off()

