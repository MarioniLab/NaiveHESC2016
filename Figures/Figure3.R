source("Figures/central.R")

# Figure 3A
library(scater)
library(ggplot2)
library(grid)

tsne1 <- readRDS("analysis/results-diff/naive_tsne1.Rds")
tsne2 <- readRDS("analysis/results-diff/naive_tsne2.Rds")
tsne3 <- readRDS("analysis/results-diff/naive_tsne3.Rds")

pdf(file=file.path(figdir, "3a.pdf"), height=5, width=10, useDingbats=FALSE)

par(mfrow=c(1,3), mar = c(5.1, 5.1, 4.1, 6.1), las = 1)
plot_data <- ggplot_build(tsne1)
plot(plot_data$data[[1]]$x, plot_data$data[[1]]$y, col = plot_data$data[[1]]$fill, 
     pch=16, xlab = plot_data$plot$labels$x , ylab = plot_data$plot$labels$y, cex.lab = 1.5 , cex=1.5)
tmp <- ggplot_gtable(plot_data) 
leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
legend <- tmp$grobs[[leg]] 
legend$vp$x <- unit(.295, "npc")
legend$vp$y <- unit(.5, "npc")
grid.draw(legend)

par(mar = c(5.1, 4.1, 4.1, 6.1), las = 1)
plot_data <- ggplot_build(tsne2)
plot(plot_data$data[[1]]$x, plot_data$data[[1]]$y, main = 'naive', col = plot_data$data[[1]]$fill, 
     pch=16, xlab = plot_data$plot$labels$x , ylab ='', cex.lab = 1.5 , cex=1.5)
tmp <- ggplot_gtable(plot_data) 
leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
legend <- tmp$grobs[[leg]] 
legend$vp$x <- unit(.62, "npc")
legend$vp$y <- unit(.5, "npc")
grid.draw(legend)

par(mar = c(5.1, 4.1, 4.1, 6.1), las = 1)
plot_data <- ggplot_build(tsne3)
plot(plot_data$data[[1]]$x, plot_data$data[[1]]$y, col = plot_data$data[[1]]$fill, 
     pch=16, xlab = plot_data$plot$labels$x , ylab = '', cex.lab = 1.5 , cex=1.5)
tmp <- ggplot_gtable(plot_data) 
leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
legend <- tmp$grobs[[leg]] 
legend$vp$x <- unit(.96, "npc")
legend$vp$y <- unit(.5, "npc")
grid.draw(legend)

dev.off()

# Figure 3B
tsne1 <- readRDS("analysis/results-diff/primed_tsne1.Rds")
tsne2 <- readRDS("analysis/results-diff/primed_tsne2.Rds")
tsne3 <- readRDS("analysis/results-diff/primed_tsne3.Rds")

pdf(file=file.path(figdir, "3b.pdf"), height=5, width=10, useDingbats=FALSE)

par(mfrow=c(1,3), mar = c(5.1, 5.1, 4.1, 6.1), las = 1)
plot_data <- ggplot_build(tsne1)
plot(plot_data$data[[1]]$x, plot_data$data[[1]]$y, col = plot_data$data[[1]]$fill, 
     pch=16, xlab = plot_data$plot$labels$x , ylab = plot_data$plot$labels$y, cex.lab = 1.5 , cex=1.5)
tmp <- ggplot_gtable(plot_data) 
leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
legend <- tmp$grobs[[leg]] 
legend$vp$x <- unit(.295, "npc")
legend$vp$y <- unit(.5, "npc")
grid.draw(legend)

par(mar = c(5.1, 4.1, 4.1, 6.1), las = 1)
plot_data <- ggplot_build(tsne2)
plot(plot_data$data[[1]]$x, plot_data$data[[1]]$y, main = 'primed', col = plot_data$data[[1]]$fill, 
     pch=16, xlab = plot_data$plot$labels$x , ylab ='', cex.lab = 1.5 , cex=1.5)
tmp <- ggplot_gtable(plot_data) 
leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
legend <- tmp$grobs[[leg]] 
legend$vp$x <- unit(.62, "npc")
legend$vp$y <- unit(.5, "npc")
grid.draw(legend)

par(mar = c(5.1, 4.1, 4.1, 6.1), las = 1)
plot_data <- ggplot_build(tsne3)
plot(plot_data$data[[1]]$x, plot_data$data[[1]]$y, col = plot_data$data[[1]]$fill, 
     pch=16, xlab = plot_data$plot$labels$x , ylab ='', cex.lab = 1.5 , cex=1.5)
tmp <- ggplot_gtable(plot_data) 
leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
legend <- tmp$grobs[[leg]] 
legend$vp$x <- unit(.96, "npc")
legend$vp$y <- unit(.5, "npc")
grid.draw(legend)

dev.off()

# 3C

library(pheatmap)
library(gridExtra)

naive.epi.heat <- readRDS('analysis/results-correlations/naive_mat.Rds')

heat.dists <- dist(naive.epi.heat)
heat.tree <- hclust(heat.dists)

heat.dist2 <- dist(t(naive.epi.heat))
heat.tree2 <- hclust(heat.dist2)

naive.epi.heat[which(naive.epi.heat>0.5, arr.ind = TRUE)] <- 0.5
naive.epi.heat[which(naive.epi.heat< (-0.5), arr.ind = TRUE)] <- -0.5

pdf(file.path(figdir, "3c.pdf"), width=5, height = 5, onefile=FALSE)
pheatmap(naive.epi.heat, breaks=seq(-0.5, 0.5, length.out=101), main = 'naive', color = colorRampPalette(c("navy", "white", "orangered"))(101),
         treeheight_row = 0, treeheight_col = 0, cluster_rows = heat.tree, cluster_cols = heat.tree2)
dev.off()


# 3D

primed.epi.heat <- readRDS('analysis/results-correlations/primed_mat.Rds')

primed.epi.heat[which(primed.epi.heat>0.5, arr.ind = TRUE)] <- 0.5
primed.epi.heat[which(primed.epi.heat< (-0.5), arr.ind = TRUE)] <- -0.5   ##to save each plot into a list. note the [[4]]

pdf(file.path(figdir, "3d.pdf"), width=5, height = 5, onefile=FALSE)
pheatmap(primed.epi.heat, breaks=seq(-0.5, 0.5, length.out=101), main = 'primed', color = colorRampPalette(c("navy", "white", "orangered"))(101),
         treeheight_row = 0, treeheight_col = 0, cluster_rows = heat.tree, cluster_cols = heat.tree2)
dev.off()

### Supplements

# Figure S3A
sce <- readRDS("analysis/results-preprocess/sce_all.rds")
pops <- read.table(file.path("analysis/results-transition", "groups.tsv"), stringsAsFactors=FALSE, header=TRUE)
sce$phenotype <- pops$Type
sce_naive <- sce[,sce$phenotype == 'naive']
sce_primed <- sce[,sce$phenotype == 'primed']

fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16), title = element_text(size=12))
naive.plot <- plotPCA(sce_naive, colour_by = "CDK1") + fontsize
primed.plot <- plotPCA(sce_primed, colour_by = "CDK1") + fontsize


pdf(file=file.path(figdir, "s3a.pdf"), width=10, useDingbats=FALSE)
par(mfrow=c(1,2), mar = c(5.1, 5.1, 4.1, 4.1), las = 1)
plot_data <- ggplot_build(naive.plot)
plot(plot_data$data[[1]]$x, plot_data$data[[1]]$y, main = 'naive', col = plot_data$data[[1]]$fill, 
           pch=16, xlab = plot_data$plot$labels$x , ylab = plot_data$plot$labels$y, cex.lab = 1.5 , cex=1.5)
par(mar = c(5.1, 4.1, 4.1, 5.1), las = 1)
plot_data <- ggplot_build(primed.plot)
plot(plot_data$data[[1]]$x, plot_data$data[[1]]$y, main = 'primed', col = plot_data$data[[1]]$fill, 
     pch=16, xlab = plot_data$plot$labels$x , ylab = plot_data$plot$labels$y, cex.lab = 1.5 , cex=1.5)
tmp <- ggplot_gtable(plot_data) 
leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
legend <- tmp$grobs[[leg]] 
legend$vp$x <- unit(.95, "npc")
legend$vp$y <- unit(.5, "npc")
grid.draw(legend)
dev.off()


# Figure S3B
colour.code <- readRDS('analysis/results-correlations/fdr.corr.Rds')

colour.scheme <- c('grey', 'black', naive.col, primed.col)
colour.code <- colour.scheme[colour.code]

pdf(file.path(figdir, "s3b.pdf"), width=10, height = 10, onefile=FALSE)
plot(naive.epi.heat, primed.epi.heat, xlim = c(-0.5, 0.5) , ylim = c(-0.5, 0.5), xlab = 'Naive correlations', ylab='Primed correlations', pch=16, col = colour.code)
legend('topleft', legend = c('Total: 625', 'Sig in both: 15', 'Sig in naive: 282', 'Sig in primed: 4'), 
       pch = 16, col = c('grey', 'black', naive.col, primed.col), bty='n', cex=1.5)
dev.off()

# Figure S3C
tsne4 <- readRDS("analysis/results-diff/naive_tsne4.Rds")
pdf(file=file.path(figdir, "s3c.pdf"), width=10, useDingbats=FALSE)
par(mfrow=c(1,2), mar = c(5.1, 5.1, 4.1, 4.1), las = 1)
plot_data <- ggplot_build(tsne4)
plot(plot_data$data[[1]]$x, plot_data$data[[1]]$y, main = 'naive', col = plot_data$data[[1]]$fill, 
     pch=16, xlab = plot_data$plot$labels$x , ylab = plot_data$plot$labels$y, cex.lab = 1.5 , cex=1.5)
par(mar = c(5.1, 4.1, 4.1, 5.1), las = 1)
tsne4 <- readRDS("analysis/results-diff/primed_tsne4.Rds")
plot_data <- ggplot_build(tsne4)
plot(plot_data$data[[1]]$x, plot_data$data[[1]]$y, main = 'primed', col = plot_data$data[[1]]$fill, 
     pch=16, xlab = plot_data$plot$labels$x , ylab = plot_data$plot$labels$y, cex.lab = 1.5 , cex=1.5)
tmp <- ggplot_gtable(plot_data) 
leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
legend <- tmp$grobs[[leg]] 
legend$vp$x <- unit(.95, "npc")
legend$vp$y <- unit(.5, "npc")
grid.draw(legend)
dev.off()
